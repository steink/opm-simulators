/*
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

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

#ifndef OPM_GRID_DATA_OUTPUT_HPP
#define OPM_GRID_DATA_OUTPUT_HPP

#include <opm/common/ErrorMacros.hpp>

#include <dune/grid/common/partitionset.hh>

#include <iosfwd>
#include <string>

/** @file
    @brief Allows model geometry data to be passed to external code - via a copy
   direct to input pointers.

    This data extractor provides the full set of vertices (corresponding to
   Dune::Partition::all) and then allows a user to specify Dune sub-partitions
   to get the references into the vertex array and element (aka cell) types for
   the sub-partition. This allows the full set of vertices to be reused for
    visualisation of the various sub-partitions, at the expense of copying all
   the vertices. Typically a user is interested in the interiorBorder elements
   which make use of the bulk (~80%) of the vertices. This saves having to
   renumber the indexes to the vertices for the sub-partitions. The vertex data
   can be retrieved as seperate x, y and z arrays, or as a single array of
   array of structures, or as single structure of arrays based array.


    Example:
        // From the opm-simulators repository
        #include <opm/simulators/utils/GridDataOutput.hpp>

        // N.B. does not seem to be able to be allocated with new operator.
        Opm::GridDataOutput::SimMeshDataAccessor geomData(gridView,
                                               Dune::Partition::interior );

        geomData.printGridDetails();

        int nvert = geomData.getNVertices();
        // example using seperate x, y and z arrays
        int nvert = geomData.getNVertices();
        double * x_vert = new double[nvert];
        double * y_vert = new double[nvert];
        double * z_vert = new double[nvert];
        geomData.writeGridPoints(x_vert,y_vert,z_vert, nvert);

        ... do something with vertex data x_vert, y_vert and z_vert ....

        delete [] x_vert;
        delete [] y_vert;
        delete [] z_vert;

        // example using AOS
        double * xyz_vert_aos = new double[nvert*3];
        geomData.writeGridPoints_AOS(xyz_vert_aos, nvert);

        ... do something with vertex data xyz_vert_aos....

        delete [] xyz_vert_aos;


        // example using SOA with std::vector<double>
        std::vector<double> xyz_vert_soa(nvert*3);
        geomData.writeGridPoints_SOA(xyz_vert_soa);

        ... do something with vertex data xyz_vert_soa....


*/

namespace Opm::GridDataOutput
{
/**
 * Allows selection of order of vertices in writeConnectivity()
 */
enum ConnectivityVertexOrder { DUNE = 0, VTK = 1 };

template <class GridView, unsigned int partitions>
class SimMeshDataAccessor
{
public:
    /**
     * @brief Construct a SimMeshDataAccessor working on a specific GridView and
     * specialize to a Dune::PartitionSet<>.
     *
     * @param gridView The gridView
     * @param PartitionSet<> the set of cells from which to extract geometric data
     *
     *  The PartitionSet of the data can be specified from one of:
     *   Dune::Partitions::all
     *   Dune::Partitions::interior
     *   Dune::Partitions::border
     *   Dune::Partitions::overlap
     *   Dune::Partitions::front
     *   Dune::Partitions::ghost
     *   Dune::Partitions::interiorBorder
     *   Dune::Partitions::interiorBorderOverlap
     *   Dune::Partitions::interiorBorderOverlapFront
     *   Dune::Partitions::all
     *
     * N.B. To visualise 'field' data on the extracted grid mesh then the field
     * variable should contain at least as many vlaues as the mesh has cells
     * (ncells_) or vertices (nvertices_) depending on if data is cell centred or
     * vertex centred, respectively.
     *
     *  As we are templated on the Dune::PartitionSet<partitions>, values for
     * ncorners_, nvertices_ and ncells_ cannot change
     *
     *  This class does not work with grids containing polyhedral cells (well, it
     * has not been tested with this kind of grid data). The user should call
     * polyhedralCellPresent() to test if polyhedral cells are present and decide
     * what they want to do before copying data using the data accessor methods.
     */
    explicit SimMeshDataAccessor(const GridView& gridView,
                                 Dune::PartitionSet<partitions> dunePartition);

    /**
      Checks for cells that have polyhedral type within the current partition of
      cells

      Returns true if a polyhedral sell is found. If this is the case then this
      partition is not going to be available for visualisation as this class does
      not yet handle polyhedral cells.
   */
    bool polyhedralCellPresent() const;

    /**
        Count the vertices, cells and corners.

        Count all the vertices ( the Dune::Partitions::all partition ) as then we
        do not need to renumber the vertices as all the subsets use references to
        the full set.
    */
    void countEntities();

    /**
      Write the positions of vertices - directly to the pointers given in
      parameters

      @param x_inout to be filled with x coordinate verticies
      @param y_inout to be filled with y coordinate verticies
      @param y_inout to be filled with z coordinate verticies
      @param max_size the maximum number of elements of type T that can be
             written to the input pointer memory regions.

      Returns the number of vertices written
   */
    template <typename T>
    long writeGridPoints(T* x_inout, T* y_inout, T* z_inout, long max_size = 0) const;

    /**
      Write the positions of vertices - directly to the pointers given in
      parameters

      @param x_inout to be filled with x coordinate verticies
      @param y_inout to be filled with y coordinate verticies
      @param y_inout to be filled with z coordinate verticies

      All parameters must have a size() and data() method (e.g. a std::vector<T>)
      and the current size() must be big enough

      Returns the number of vertices written
   */
    template <typename VectType>
    long writeGridPoints(VectType& x_inout, VectType& y_inout, VectType& z_inout) const;

    /**
      Write the positions of vertices - directly to the pointers given in
      parameters as Array of Structures x,y,z,x,y,z,x,y,z,...

      @param xyz_inout is the array to be filled with x,y,z coordinate verticies.
      @param max_size  is the maximum number x,y,z structures with elements of type T
                    that can be written to the input pointer memory regions.

      Returns the number of vertices written
   */
    template <typename T>
    long writeGridPoints_AOS(T* xyz_inout, long max_size = 0) const;

    /**
      Write the positions of vertices - directly to the pointers given in
      parameters as Array of Structures x,y,z,x,y,z,x,y,z,...

      @param xyz_inout is the array to be filled with x,y,z coordinate verticies.
              The object VectType must have a size() and data() method (e.g. a std::vector<T>)

      Returns the number of vertices written
   */
    template <typename VectType>
    long writeGridPoints_AOS(VectType& xyz_inout) const;

    /**
        Write the positions of vertices - directly to the pointers given in
        parameters as Structure of Arrays:  x,x,x,...,y,y,y,...,z,z,z,...

        @param xyz_inout is the array to be filled with x,y,z coordinate verticies.
        @param max_size  number of verticies (x,...y,...z,... structures) with elements of type T
                      that can be written to the input pointer memory regions.

        Returns the number of vertices written
     */
    template <typename T>
    long writeGridPoints_SOA(T* xyz_inout, long max_size = 0) const;

    /**
      Write the positions of vertices - directly to the pointers given in
      parameters as Structure of Arrays:  x,x,x,...,y,y,y,...,z,z,z,...

      @param xyz_inout is the array to be filled with x,y,z coordinate verticies.
              The object VectType must have a size() and data() method (e.g. a std::vector<T>)

      Returns the number of vertices written
   */
    template <typename VectType>
    long writeGridPoints_SOA(VectType& xyz_inout) const;

    /**
    * Write the connectivity array - directly to the pointer given in parameter 1
      Reorders the indices as selected either in DUNE order or  VTK order.

      @param connectivity_inout is the array to be filled with connectivity indexes
             (i.e. the index into the vertex array)
      @param whichOrder, is the order that verticies are traversed to create a cell (VTK or DUNE)
      @param max_size is used to check that the space available in the input pointer
             parameter will fit the number of corner values written.
      Returns the number of corner indices written.
    */
    template <typename Integer>
    long writeConnectivity(Integer* connectivity_inout,
                           ConnectivityVertexOrder whichOrder, long max_size = 0) const;

    /**
    * Write the connectivity array - directly to a VectType object given in parameter 1
      Reorders the indices as selected either in DUNE order or  VTK order.

      @param connectivity_inout is the array to be filled with connectivity indexes
             (i.e. the index into the vertex array)
             The object VectType must have a size() and data() method (e.g. a std::vector<T>)
      @param whichOrder, is the order that verticies are traversed to create a cell (VTK or DUNE)
      @param max_size is used to check that the space available in the input pointer
             parameter will fit the number of corner values written.
      Returns the number of corner indices written.
    */
    template <typename VectType>
    long writeConnectivity(VectType& connectivity_inout,
                           ConnectivityVertexOrder whichOrder) const;

    /**
     * Write the offsets values  - directly to the pointer given in parameter 1

       @param offsets_inout is the array to be filled with offsets into the connectivity array
              (i.e. the index into the connectivity array to determine the vertices used for
              the particular cell)
       @param max_size is used to check that the space available in the input pointer
              parameter will fit the number of cell offset values written.

       Returns number of offset values written + 1
   */
    template <typename Integer>
    long writeOffsetsCells(Integer* offsets_inout, long max_size = 0) const;

    /**
      * Write the offsets values  -  directly to a VectType object given in parameter 1

        @param offsets_inout is the array to be filled with offsets into the connectivity array
               (i.e. the index into the connectivity array to determine the vertices used for
               the particular cell).
               The object VectType must have a size() and data() method (e.g. a std::vector<T>)

        Returns number of offset values written + 1
    */
    template <typename VectType>
    long writeOffsetsCells(VectType& offsets_inout) const;

    /**
     * Write the cell types values  - directly to the pointer given in parameter 1

       @param types_inout is the array to be filled with the cell types (VTK defined values)

       @param max_size is used to check that the space available in the input pointer
               parameter will fit the number of cell offset values written.

        Returns number of cells type values written
    */
    template <typename Integer>
    long writeCellTypes(Integer* types_inout, long max_size = 0) const;

    /**
     * Write the cell types values  - directly to the VectType object given in parameter 1

       @param types_inout is the array to be filled with the cell types (VTK defined values)
              The object VectType must have a size() and data() method (e.g. a std::vector<T>)

        Returns number of cells type values written
    */
    template <typename VectType>
    long writeCellTypes(VectType& types_inout) const;

    std::string getPartitionTypeString() const;

    Dune::PartitionSet<partitions> getPartition(void)
    {
        return this->dunePartition_;
    }

    void printGridDetails(std::ostream& outstr) const;

    int getNCells() const
    {
        return ncells_;
    }

    int getNVertices() const
    {
        return nvertices_;
    }

    int getNCorners() const
    {
        return ncorners_;
    }

    std::string getError() const
    {
        return error_;
    }

    void clearError()
    {
        error_.clear();
    }

    bool hasError() const
    {
        return !error_.empty();
    }

protected:
    GridView gridView_; // the grid

    Dune::PartitionSet<partitions> dunePartition_;
    unsigned int partition_value_;

    /**
    Current partition grid information
    */
    int ncells_;
    /**
    Current partition grid information
    */
    int nvertices_;
    /**
    Current partition grid information
    */
    int ncorners_;

    int dimw_; // dimensions of the input grid

private:
    std::string error_;
};

} // namespace Opm::GridDataOutput

#endif
