/*
  Copyright 2020 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/ChowPatelIlu.hpp>

namespace bda
{

    using Opm::OpmLog;

// if PARALLEL is 0:
//    Each row gets 1 workgroup, 1 workgroup can do multiple rows sequentially.
//    Each block in a row gets 1 workitem, all blocks are expected to be processed simultaneously,
//    except when the number of blocks in that row exceeds the number of workitems per workgroup.
//    In that case some workitems will process multiple blocks sequentially.
// else:
//    Each row gets 1 workgroup, 1 workgroup can do multiple rows sequentially
//    Each block in a row gets a warp of 32 workitems, of which 9 are always active.
//    Multiple blocks can be processed in parallel if a workgroup contains multiple warps.
//    If the number of blocks exceeds the number of warps, some warps will process multiple blocks sequentially.

// Notes:
// PARALLEL 0 should be able to run with any number of workitems per workgroup, but 8 and 16 tend to be quicker than 32.
// PARALLEL 1 should be run with at least 32 workitems per workgroup.
// The recommended number of workgroups for both options is Nb, which gives every row their own workgroup.
// PARALLEL 0 is generally faster, despite not having parallelization.
// only 3x3 blocks are supported

#define PARALLEL 0

#if PARALLEL
inline const char* chow_patel_ilu_sweep_s  = R"(

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// subtract blocks: a = a - b * c
// the output block has 9 entries, each entry is calculated by 1 thread
void blockMultSub(
    __local double * restrict a,
    __global const double * restrict b,
    __global const double * restrict c)
{
    const unsigned int block_size = 3;
    const unsigned int warp_size = 32;
    const unsigned int idx_t = get_local_id(0);                   // thread id in work group
    const unsigned int thread_id_in_warp = idx_t % warp_size;     // thread id in warp (32 threads)
    if(thread_id_in_warp < block_size * block_size){
        const unsigned int row = thread_id_in_warp / block_size;
        const unsigned int col = thread_id_in_warp % block_size;
        double temp = 0.0;
        for (unsigned int k = 0; k < block_size; k++) {
            temp += b[block_size * row + k] * c[block_size * k + col];
        }
        a[block_size * row + col] -= temp;
    }
}

// multiply blocks: resMat = mat1 * mat2
// the output block has 9 entries, each entry is calculated by 1 thread
void blockMult(
    __local const double * restrict mat1,
    __local const double * restrict mat2,
    __global double * restrict resMat)
{
    const unsigned int block_size = 3;
    const unsigned int warp_size = 32;
    const unsigned int idx_t = get_local_id(0);                   // thread id in work group
    const unsigned int thread_id_in_warp = idx_t % warp_size;     // thread id in warp (32 threads)
    if(thread_id_in_warp < block_size * block_size){
        const unsigned int row = thread_id_in_warp / block_size;
        const unsigned int col = thread_id_in_warp % block_size;
        double temp = 0.0;
        for (unsigned int k = 0; k < block_size; k++) {
            temp += mat1[block_size * row + k] * mat2[block_size * k + col];
        }
        resMat[block_size * row + col] = temp;
    }
}

// invert block: inverse = matrix^{-1}
// the output block has 9 entries, each entry is calculated by 1 thread
void invert(
    __global const double * restrict matrix,
    __local double * restrict inverse)
{
    const unsigned int block_size = 3;
    const unsigned int bs = block_size;                           // rename to shorter name
    const unsigned int warp_size = 32;
    const unsigned int idx_t = get_local_id(0);                   // thread id in work group
    const unsigned int thread_id_in_warp = idx_t % warp_size;     // thread id in warp (32 threads)
    if(thread_id_in_warp < block_size * block_size){
        // code generated by maple, copied from Dune::DenseMatrix
        double t4  = matrix[0] * matrix[4];
        double t6  = matrix[0] * matrix[5];
        double t8  = matrix[1] * matrix[3];
        double t10 = matrix[2] * matrix[3];
        double t12 = matrix[1] * matrix[6];
        double t14 = matrix[2] * matrix[6];

        double det = (t4 * matrix[8] - t6 * matrix[7] - t8 * matrix[8] +
                      t10 * matrix[7] + t12 * matrix[5] - t14 * matrix[4]);
        double t17 = 1.0 / det;

        const unsigned int r = thread_id_in_warp / block_size;
        const unsigned int c = thread_id_in_warp % block_size;
        const unsigned int r1 = (r+1) % bs;
        const unsigned int c1 = (c+1) % bs;
        const unsigned int r2 = (r+bs-1) % bs;
        const unsigned int c2 = (c+bs-1) % bs;
        inverse[c*bs+r] = ((matrix[r1*bs+c1] * matrix[r2*bs+c2]) - (matrix[r1*bs+c2] * matrix[r2*bs+c1])) * t17;
    }
}

// perform the fixed-point iteration
// all entries in L and U are updated once
// output is written to [LU]tmp
// aij and ujj are local arrays whose size is specified before kernel launch
__kernel void chow_patel_ilu_sweep(
    __global const double * restrict Ut_vals,
    __global const double * restrict L_vals,
    __global const double * restrict LU_vals,
    __global const int * restrict Ut_rows,
    __global const int * restrict L_rows,
    __global const int * restrict LU_rows,
    __global const int * restrict Ut_cols,
    __global const int * restrict L_cols,
    __global const int * restrict LU_cols,
    __global double * restrict Ltmp,
    __global double * restrict Utmp,
    const int Nb,
    __local double *aij,
    __local double *ujj)
{
    const int bs = 3;
    const unsigned int warp_size = 32;
    const unsigned int work_group_size = get_local_size(0);
    const unsigned int idx_b = get_global_id(0) / work_group_size;
    const unsigned int num_groups = get_num_groups(0);
    const unsigned int warps_per_group = work_group_size / warp_size;
    const unsigned int idx_t = get_local_id(0);                   // thread id in work group
    const unsigned int thread_id_in_warp = idx_t % warp_size;     // thread id in warp (32 threads)
    const unsigned int warp_id_in_group = idx_t / warp_size;
    const unsigned int lmem_offset = warp_id_in_group * bs * bs;  // each workgroup gets some lmem, but the workitems have to share it
                                                                  // every workitem in a warp has the same lmem_offset

    // for every row of L or every col of U
    for (int row = idx_b; row < Nb; row+=num_groups) {
        // Uij = (Aij - sum k=1 to i-1 {Lik*Ukj})
        int jColStart = Ut_rows[row];    // actually colPointers to U
        int jColEnd = Ut_rows[row + 1];
        // for every block on this column
        for (int ij = jColStart + warp_id_in_group; ij < jColEnd; ij+=warps_per_group) {
            int rowU1 = Ut_cols[ij]; // actually rowIndices for U
            // refine Uij element (or diagonal)
            int i1 = LU_rows[rowU1];
            int i2 = LU_rows[rowU1+1];
            int kk = 0;
            // LUmat->nnzValues[kk] is block Aij
            for(kk = i1; kk < i2; ++kk) {
                int c = LU_cols[kk];
                if (c >= row) {
                    break;
                }
            }

            // copy block Aij so operations can be done on it without affecting LUmat
            if(thread_id_in_warp < bs*bs){
                aij[lmem_offset+thread_id_in_warp] = LU_vals[kk*bs*bs + thread_id_in_warp];
            }

            int jk = L_rows[rowU1]; // points to row rowU1 in L
            // if row rowU1 is empty: skip row. The whole warp looks at the same row, so no divergence
            if (jk < L_rows[rowU1+1]) {
                int colL = L_cols[jk];
                // only check until block U(i,j) is reached
                for (int k = jColStart; k < ij; ++k) {
                    int rowU2 = Ut_cols[k]; // actually rowIndices for U
                    while (colL < rowU2) {
                        ++jk; // check next block on row rowU1 of L
                        colL = L_cols[jk];
                    }
                    if (colL == rowU2) {
                        // Aij -= (Lik * Ukj)
                        blockMultSub(aij+lmem_offset, L_vals + jk * bs * bs, Ut_vals + k * bs * bs);
                    }
                }
            }

            // Uij_new = Aij - sum
            // write result of this sweep
            if(thread_id_in_warp < bs*bs){
                Utmp[ij*bs*bs + thread_id_in_warp] = aij[lmem_offset + thread_id_in_warp];
            }
        }

        // update L
        // Lij = (Aij - sum k=1 to j-1 {Lik*Ukj}) / Ujj
        int iRowStart = L_rows[row];
        int iRowEnd = L_rows[row + 1];

        for (int ij = iRowStart + warp_id_in_group; ij < iRowEnd; ij+=warps_per_group) {
            int j = L_cols[ij];
            // // refine Lij element
            int i1 = LU_rows[row];
            int i2 = LU_rows[row+1];
            int kk = 0;
            // LUmat->nnzValues[kk] is block Aij
            for(kk = i1; kk < i2; ++kk) {
                int c = LU_cols[kk];
                if (c >= j) {
                    break;
                }
            }

            // copy block Aij so operations can be done on it without affecting LUmat
            if(thread_id_in_warp < bs*bs){
                aij[lmem_offset+thread_id_in_warp] = LU_vals[kk*bs*bs + thread_id_in_warp];
            }

            int jk = Ut_rows[j];    // actually colPointers, jk points to col j in U
            int rowU = Ut_cols[jk]; // actually rowIndices, rowU is the row of block jk
            // only check until block L(i,j) is reached
            for (int k = iRowStart; k < ij; ++k) {
                int colL = L_cols[k];
                while(rowU < colL) {
                    ++jk; // check next block on col j of U
                    rowU = Ut_cols[jk];
                }

                if(rowU == colL) {
                    // Aij -= (Lik * Ukj)
                    blockMultSub(aij+lmem_offset, L_vals + k * bs * bs , Ut_vals + jk * bs * bs);
                }
            }

            // calculate 1 / Ujj
            invert(Ut_vals + (Ut_rows[j+1] - 1) * bs * bs, ujj+lmem_offset);

            // Lij_new = (Aij - sum) / Ujj
            // write result of this sweep
            blockMult(aij+lmem_offset, ujj+lmem_offset, Ltmp + ij * bs * bs);
        }
    }
}
)";

#else

inline const char* chow_patel_ilu_sweep_s  = R"(

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// subtract blocks: a = a - b * c
// only one workitem performs this action
void blockMultSub(
    __local double * restrict a,
    __global const double * restrict b,
    __global const double * restrict c)
{
    const unsigned int block_size = 3;
    for (unsigned int row = 0; row < block_size; row++) {
        for (unsigned int col = 0; col < block_size; col++) {
            double temp = 0.0;
            for (unsigned int k = 0; k < block_size; k++) {
                temp += b[block_size * row + k] * c[block_size * k + col];
            }
            a[block_size * row + col] -= temp;
        }
    }
}

// multiply blocks: resMat = mat1 * mat2
// only one workitem performs this action
void blockMult(
    __local const double * restrict mat1,
    __local const double * restrict mat2,
    __global double * restrict resMat)
{
    const unsigned int block_size = 3;
    for (unsigned int row = 0; row < block_size; row++) {
        for (unsigned int col = 0; col < block_size; col++) {
            double temp = 0.0;
            for (unsigned int k = 0; k < block_size; k++) {
                temp += mat1[block_size * row + k] * mat2[block_size * k + col];
            }
            resMat[block_size * row + col] = temp;
        }
    }
}

// invert block: inverse = matrix^{-1}
// only one workitem performs this action
__kernel void inverter(
    __global const double * restrict matrix,
    __local double * restrict inverse)
{
    // code generated by maple, copied from Dune::DenseMatrix
    double t4  = matrix[0] * matrix[4];
    double t6  = matrix[0] * matrix[5];
    double t8  = matrix[1] * matrix[3];
    double t10 = matrix[2] * matrix[3];
    double t12 = matrix[1] * matrix[6];
    double t14 = matrix[2] * matrix[6];

    double det = (t4 * matrix[8] - t6 * matrix[7] - t8 * matrix[8] +
                  t10 * matrix[7] + t12 * matrix[5] - t14 * matrix[4]);
    double t17 = 1.0 / det;

    inverse[0] =  (matrix[4] * matrix[8] - matrix[5] * matrix[7]) * t17;
    inverse[1] = -(matrix[1] * matrix[8] - matrix[2] * matrix[7]) * t17;
    inverse[2] =  (matrix[1] * matrix[5] - matrix[2] * matrix[4]) * t17;
    inverse[3] = -(matrix[3] * matrix[8] - matrix[5] * matrix[6]) * t17;
    inverse[4] =  (matrix[0] * matrix[8] - t14) * t17;
    inverse[5] = -(t6 - t10) * t17;
    inverse[6] =  (matrix[3] * matrix[7] - matrix[4] * matrix[6]) * t17;
    inverse[7] = -(matrix[0] * matrix[7] - t12) * t17;
    inverse[8] =  (t4 - t8) * t17;
}

// perform the fixed-point iteration
// all entries in L and U are updated once
// output is written to [LU]tmp
// aij and ujj are local arrays whose size is specified before kernel launch
__kernel void chow_patel_ilu_sweep(
    __global const double * restrict Ut_vals,
    __global const double * restrict L_vals,
    __global const double * restrict LU_vals,
    __global const int * restrict Ut_rows,
    __global const int * restrict L_rows,
    __global const int * restrict LU_rows,
    __global const int * restrict Ut_cols,
    __global const int * restrict L_cols,
    __global const int * restrict LU_cols,
    __global double * restrict Ltmp,
    __global double * restrict Utmp,
    const int Nb,
    __local double *aij,
    __local double *ujj)
{
    const int bs = 3;

    const unsigned int warp_size = 32;
    const unsigned int work_group_size = get_local_size(0);
    const unsigned int idx_b = get_global_id(0) / work_group_size;
    const unsigned int num_groups = get_num_groups(0);
    const unsigned int warps_per_group = work_group_size / warp_size;
    const unsigned int idx_t = get_local_id(0);                   // thread id in work group
    const unsigned int thread_id_in_warp = idx_t % warp_size;     // thread id in warp (32 threads)
    const unsigned int warp_id_in_group = idx_t / warp_size;

    // for every row of L or every col of U
    for (int row = idx_b; row < Nb; row+=num_groups) {
        // Uij = (Aij - sum k=1 to i-1 {Lik*Ukj})
        int jColStart = Ut_rows[row];    // actually colPointers to U
        int jColEnd = Ut_rows[row + 1];
        // for every block on this column
        for (int ij = jColStart + idx_t; ij < jColEnd; ij+=work_group_size) {
            int rowU1 = Ut_cols[ij]; // actually rowIndices for U
            // refine Uij element (or diagonal)
            int i1 = LU_rows[rowU1];
            int i2 = LU_rows[rowU1+1];
            int kk = 0;
            // LUmat->nnzValues[kk] is block Aij
            for(kk = i1; kk < i2; ++kk) {
                int c = LU_cols[kk];
                if (c >= row) {
                    break;
                }
            }

            // copy block Aij so operations can be done on it without affecting LUmat
            for(int z = 0; z < bs*bs; ++z){
                aij[idx_t*bs*bs+z] = LU_vals[kk*bs*bs + z];
            }

            int jk = L_rows[rowU1];
            // if row rowU1 is empty: do not sum. The workitems have different rowU1 values, divergence is possible
            int colL = (jk < L_rows[rowU1+1]) ? L_cols[jk] : Nb;

            // only check until block U(i,j) is reached
            for (int k = jColStart; k < ij; ++k) {
                int rowU2 = Ut_cols[k]; // actually rowIndices for U
                while (colL < rowU2) {
                    ++jk; // check next block on row rowU1 of L
                    colL = L_cols[jk];
                }
                if (colL == rowU2) {
                    // Aij -= (Lik * Ukj)
                    blockMultSub(aij+idx_t*bs*bs, L_vals + jk * bs * bs, Ut_vals + k * bs * bs);
                }
            }

            // Uij_new = Aij - sum
            // write result of this sweep
            for(int z = 0; z < bs*bs; ++z){
                Utmp[ij*bs*bs + z] = aij[idx_t*bs*bs+z];
            }
        }

        // update L
        // Lij = (Aij - sum k=1 to j-1 {Lik*Ukj}) / Ujj
        int iRowStart = L_rows[row];
        int iRowEnd = L_rows[row + 1];

        for (int ij = iRowStart + idx_t; ij < iRowEnd; ij+=work_group_size) {
            int j = L_cols[ij];
            // // refine Lij element
            int i1 = LU_rows[row];
            int i2 = LU_rows[row+1];
            int kk = 0;
            // LUmat->nnzValues[kk] is block Aij
            for(kk = i1; kk < i2; ++kk) {
                int c = LU_cols[kk];
                if (c >= j) {
                    break;
                }
            }

            // copy block Aij so operations can be done on it without affecting LUmat
            for(int z = 0; z < bs*bs; ++z){
                aij[idx_t*bs*bs+z] = LU_vals[kk*bs*bs + z];
            }

            int jk = Ut_rows[j];    // actually colPointers, jk points to col j in U
            int rowU = Ut_cols[jk]; // actually rowIndices, rowU is the row of block jk
            // only check until block L(i,j) is reached
            for (int k = iRowStart; k < ij; ++k) {
                int colL = L_cols[k];
                while(rowU < colL) {
                    ++jk; // check next block on col j of U
                    rowU = Ut_cols[jk];
                }

                if(rowU == colL) {
                    // Aij -= (Lik * Ukj)
                    blockMultSub(aij+idx_t*bs*bs, L_vals + k * bs * bs , Ut_vals + jk * bs * bs);
                }
            }
            // calculate 1 / ujj
            inverter(Ut_vals + (Ut_rows[j+1] - 1) * bs * bs, ujj+idx_t*bs*bs);

            // Lij_new = (Aij - sum) / Ujj
            // write result of this sweep
            blockMult(aij+idx_t*bs*bs, ujj+idx_t*bs*bs, Ltmp + ij * bs * bs);
        }
    }
}
)";

#endif







void ChowPatelIlu::decomposition(
    cl::CommandQueue *queue, cl::Context *context,
    int *Ut_ptrs, int *Ut_idxs, double *Ut_vals, int Ut_nnzbs,
    int *L_rows, int *L_cols, double *L_vals, int L_nnzbs,
    int *LU_rows, int *LU_cols, double *LU_vals, int LU_nnzbs,
    int Nb, int num_sweeps, int verbosity)
{
    const int block_size = 3;

    try {
        // just put everything in the capture list
        std::call_once(initialize_flag, [&](){
            cl::Program::Sources source(1, std::make_pair(chow_patel_ilu_sweep_s, strlen(chow_patel_ilu_sweep_s)));  // what does this '1' mean? cl::Program::Sources is of type 'std::vector<std::pair<const char*, long unsigned int> >'
            cl::Program program = cl::Program(*context, source, &err);
            if (err != CL_SUCCESS) {
                OPM_THROW(std::logic_error, "ChowPatelIlu OpenCL could not create Program");
            }

            std::vector<cl::Device> devices = context->getInfo<CL_CONTEXT_DEVICES>();
            program.build(devices);

            chow_patel_ilu_sweep_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                     cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                     cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                     cl::Buffer&, cl::Buffer&,
                                                     const int, cl::LocalSpaceArg, cl::LocalSpaceArg>(cl::Kernel(program, "chow_patel_ilu_sweep", &err)));
            if (err != CL_SUCCESS) {
                OPM_THROW(std::logic_error, "ChowPatelIlu OpenCL could not create Kernel");
            }

            // allocate GPU memory
            d_Ut_vals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * Ut_nnzbs * block_size * block_size);
            d_L_vals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * L_nnzbs * block_size * block_size);
            d_LU_vals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * LU_nnzbs * block_size * block_size);
            d_Ut_ptrs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Nb+1));
            d_L_rows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Nb+1));
            d_LU_rows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Nb+1));
            d_Ut_idxs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * Ut_nnzbs);
            d_L_cols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * L_nnzbs);
            d_LU_cols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * LU_nnzbs);
            d_Ltmp = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * L_nnzbs * block_size * block_size);
            d_Utmp = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * Ut_nnzbs * block_size * block_size);

            Dune::Timer t_copy_pattern;
            events.resize(6);
            err |= queue->enqueueWriteBuffer(d_Ut_ptrs, CL_FALSE, 0, sizeof(int) * (Nb+1), Ut_ptrs, nullptr, &events[0]);
            err |= queue->enqueueWriteBuffer(d_L_rows, CL_FALSE, 0, sizeof(int) * (Nb+1), L_rows, nullptr, &events[1]);
            err |= queue->enqueueWriteBuffer(d_LU_rows, CL_FALSE, 0, sizeof(int) * (Nb+1), LU_rows, nullptr, &events[2]);
            err |= queue->enqueueWriteBuffer(d_Ut_idxs, CL_FALSE, 0, sizeof(int) * Ut_nnzbs, Ut_idxs, nullptr, &events[3]);
            err |= queue->enqueueWriteBuffer(d_L_cols, CL_FALSE, 0, sizeof(int) * L_nnzbs, L_cols, nullptr, &events[4]);
            err |= queue->enqueueWriteBuffer(d_LU_cols, CL_FALSE, 0, sizeof(int) * LU_nnzbs, LU_cols, nullptr, &events[5]);
            cl::WaitForEvents(events);
            events.clear();
            if (verbosity >= 4){
                std::ostringstream out;
                out << "ChowPatelIlu copy sparsity pattern time: " << t_copy_pattern.stop() << " s";
                OpmLog::info(out.str());
            }
            if (verbosity >= 2){
                std::ostringstream out;
                out << "ChowPatelIlu PARALLEL: " << PARALLEL;
                OpmLog::info(out.str());
            }
        });


        // copy to GPU
        Dune::Timer t_copy1;
        events.resize(3);
        err = queue->enqueueWriteBuffer(d_Ut_vals, CL_FALSE, 0, sizeof(double) * Ut_nnzbs * block_size * block_size, Ut_vals, nullptr, &events[0]);
        err |= queue->enqueueWriteBuffer(d_L_vals, CL_FALSE, 0, sizeof(double) * L_nnzbs * block_size * block_size, L_vals, nullptr, &events[1]);
        err |= queue->enqueueWriteBuffer(d_LU_vals, CL_FALSE, 0, sizeof(double) * LU_nnzbs * block_size * block_size, LU_vals, nullptr, &events[2]);
        cl::WaitForEvents(events);
        events.clear();
        if (verbosity >= 4){
            std::ostringstream out;
            out << "ChowPatelIlu copy1 time: " << t_copy1.stop() << " s";
            OpmLog::info(out.str());
        }
        if (err != CL_SUCCESS) {
            // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "ChowPatelIlu OpenCL enqueueWriteBuffer error");
        }

        // call kernel
        for (int sweep = 0; sweep < num_sweeps; ++sweep) {
            // normally, L_vals and Ltmp are swapped after the sweep is done
            // these conditionals implement that without actually swapping pointers
            // 1st sweep reads X_vals, writes to Xtmp
            // 2nd sweep reads Xtmp, writes to X_vals
            auto *Larg1 = (sweep % 2 == 0) ? &d_L_vals : &d_Ltmp;
            auto *Larg2 = (sweep % 2 == 0) ? &d_Ltmp : &d_L_vals;
            auto *Uarg1 = (sweep % 2 == 0) ? &d_Ut_vals : &d_Utmp;
            auto *Uarg2 = (sweep % 2 == 0) ? &d_Utmp : &d_Ut_vals;
            int num_work_groups = Nb;
#if PARALLEL
            int work_group_size = 32;
#else
            int work_group_size = 16;
#endif
            int total_work_items = num_work_groups * work_group_size;
            int lmem_per_work_group = work_group_size * block_size * block_size * sizeof(double);
            Dune::Timer t_kernel;
            event = (*chow_patel_ilu_sweep_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                *Uarg1, *Larg1, d_LU_vals,
                d_Ut_ptrs, d_L_rows, d_LU_rows,
                d_Ut_idxs, d_L_cols, d_LU_cols,
                *Larg2, *Uarg2, Nb, cl::Local(lmem_per_work_group), cl::Local(lmem_per_work_group));
            event.wait();
            if (verbosity >= 4){
                std::ostringstream out;
                out << "ChowPatelIlu sweep kernel time: " << t_kernel.stop() << " s";
                OpmLog::info(out.str());
            }
        }

        // copy back
        Dune::Timer t_copy2;
        events.resize(2);
        if (num_sweeps % 2 == 0) {
            err = queue->enqueueReadBuffer(d_Ut_vals, CL_FALSE, 0, sizeof(double) * Ut_nnzbs * block_size * block_size, Ut_vals, nullptr, &events[0]);
            err |= queue->enqueueReadBuffer(d_L_vals, CL_FALSE, 0, sizeof(double) * L_nnzbs * block_size * block_size, L_vals, nullptr, &events[1]);
        } else {
            err = queue->enqueueReadBuffer(d_Utmp, CL_FALSE, 0, sizeof(double) * Ut_nnzbs * block_size * block_size, Ut_vals, nullptr, &events[0]);
            err |= queue->enqueueReadBuffer(d_Ltmp, CL_FALSE, 0, sizeof(double) * L_nnzbs * block_size * block_size, L_vals, nullptr, &events[1]);
        }
        cl::WaitForEvents(events);
        events.clear();
        if (verbosity >= 4){
            std::ostringstream out;
            out << "ChowPatelIlu copy2 time: " << t_copy2.stop() << " s";
            OpmLog::info(out.str());
        }
        if (err != CL_SUCCESS) {
            // enqueueReadBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "ChowPatelIlu OpenCL enqueueReadBuffer error");
        }

    } catch (const cl::Error& error) {
        std::ostringstream oss;
        oss << "OpenCL Error: " << error.what() << "(" << error.err() << ")\n";
        oss << getErrorString(error.err()) << std::endl;
        // rethrow exception
        OPM_THROW(std::logic_error, oss.str());
    } catch (const std::logic_error& error) {
        // rethrow exception by OPM_THROW in the try{}
        throw error;
    }
}


} // end namespace bda

