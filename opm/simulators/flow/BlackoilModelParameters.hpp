/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED
#define OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/SubDomain.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace Opm::Properties {

namespace TTag {
struct FlowModelParameters {};
}

template<class TypeTag, class MyTypeTag>
struct EclDeckFileName {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DbhpMaxRel {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DwellFractionMax {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxResidualAllowed {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct RelaxedMaxPvFraction {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ToleranceMb {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ToleranceMbRelaxed {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ToleranceCnv {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ToleranceCnvRelaxed {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ToleranceWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ToleranceWellControl {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxWelleqIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct UseMultisegmentWell {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxSinglePrecisionDays {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MinStrictCnvIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MinStrictMbIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolveWelleqInitially {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct UpdateEquationsScaling {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct UseUpdateStabilization {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MatrixAddWellContributions {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableWellOperabilityCheck {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableWellOperabilityCheckIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DebugEmitCellPartition {
    using type = UndefinedProperty;
};
// parameters for multisegment wells
template<class TypeTag, class MyTypeTag>
struct TolerancePressureMsWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxPressureChangeMsWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxInnerIterMsWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct StrictInnerIterWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct RelaxedWellFlowTol {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct StrictOuterIterWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct RelaxedPressureTolMsw {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct RegularizationFactorWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxNewtonIterationsWithInnerWellIterations  {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ShutUnsolvableWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxInnerIterWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct AlternativeWellRateInit {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaximumNumberOfWellSwitches {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct UseAverageDensityMsWells {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LocalWellSolveControlSwitching {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct UseImplicitIpr {
    using type = UndefinedProperty;
};
// Network solver parameters
template<class TypeTag, class MyTypeTag>
struct NetworkMaxStrictIterations {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct NetworkMaxIterations {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct NonlinearSolver {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LocalSolveApproach {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MaxLocalSolveIterations {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LocalToleranceScalingMb {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LocalToleranceScalingCnv {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct NlddNumInitialNewtonIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct NumLocalDomains {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LocalDomainsPartitioningImbalance {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LocalDomainsPartitioningMethod {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LocalDomainsOrderingMeasure {
    using type = UndefinedProperty;
};
template<class TypeTag>
struct DbhpMaxRel<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct DwellFractionMax<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.2;
};
template<class TypeTag>
struct MaxResidualAllowed<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e7;
};
template<class TypeTag>
struct RelaxedMaxPvFraction<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.03;
};
template<class TypeTag>
struct ToleranceMb<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-6;
};
template<class TypeTag>
struct ToleranceMbRelaxed<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-6;
};
template<class TypeTag>
struct ToleranceCnv<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};
template<class TypeTag>
struct ToleranceCnvRelaxed<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1;
};
template<class TypeTag>
struct ToleranceWells<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-4;
};
template<class TypeTag>
struct ToleranceWellControl<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-7;
};
template<class TypeTag>
struct MaxWelleqIter<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 30;
};
template<class TypeTag>
struct UseMultisegmentWell<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct MaxSinglePrecisionDays<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 20.0;
};
template<class TypeTag>
struct MinStrictCnvIter<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 0;
};
template<class TypeTag>
struct MinStrictMbIter<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = -1;
};
template<class TypeTag>
struct SolveWelleqInitially<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct UpdateEquationsScaling<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct UseUpdateStabilization<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct MatrixAddWellContributions<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct TolerancePressureMsWells<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.01*1e5;
};
template<class TypeTag>
struct MaxPressureChangeMsWells<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10*1e5;
};
template<class TypeTag>
struct MaxNewtonIterationsWithInnerWellIterations<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 8;
};
template<class TypeTag>
struct MaxInnerIterMsWells<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 100;
};
template<class TypeTag>
struct MaxInnerIterWells<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 50;
};
template<class TypeTag>
struct ShutUnsolvableWells<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct AlternativeWellRateInit<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct StrictOuterIterWells<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 6;
};
template<class TypeTag>
struct StrictInnerIterWells<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 40;
};
template<class TypeTag>
struct RegularizationFactorWells<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100;
};
template<class TypeTag>
struct EnableWellOperabilityCheck<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableWellOperabilityCheckIter<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct DebugEmitCellPartition<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct RelaxedWellFlowTol<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};
template<class TypeTag>
struct RelaxedPressureTolMsw<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0e4;
};
template<class TypeTag>
struct MaximumNumberOfWellSwitches<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 3;
};
template<class TypeTag>
struct UseAverageDensityMsWells<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct LocalWellSolveControlSwitching<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct UseImplicitIpr<TypeTag, TTag::FlowModelParameters> {
    static constexpr bool value = true;
};

// Network solver parameters
template<class TypeTag>
struct NetworkMaxStrictIterations<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 100;
};
template<class TypeTag>
struct NetworkMaxIterations<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 200;
};
template<class TypeTag>
struct NonlinearSolver<TypeTag, TTag::FlowModelParameters> {
    static constexpr auto value = "newton";
};
template<class TypeTag>
struct LocalSolveApproach<TypeTag, TTag::FlowModelParameters> {
    static constexpr auto value = "gauss-seidel";
};
template<class TypeTag>
struct MaxLocalSolveIterations<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = 20;
};
template<class TypeTag>
struct LocalToleranceScalingMb<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct LocalToleranceScalingCnv<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.1;
};
template<class TypeTag>
struct NlddNumInitialNewtonIter<TypeTag, TTag::FlowModelParameters> {
    using type = int;
    static constexpr auto value = type{1};
};
template<class TypeTag>
struct NumLocalDomains<TypeTag, TTag::FlowModelParameters> {
    using type = int;
    static constexpr auto value = 0;
};
template<class TypeTag>
struct LocalDomainsPartitioningImbalance<TypeTag, TTag::FlowModelParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr auto value = type{1.03};
};
template<class TypeTag>
struct LocalDomainsPartitioningMethod<TypeTag, TTag::FlowModelParameters> {
    static constexpr auto value = "zoltan";
};
template<class TypeTag>
struct LocalDomainsOrderingMeasure<TypeTag, TTag::FlowModelParameters> {
    static constexpr auto value = "maxpressure";
};
// if openMP is available, determine the number threads per process automatically.
#if _OPENMP
template<class TypeTag>
struct ThreadsPerProcess<TypeTag, TTag::FlowModelParameters> {
    static constexpr int value = -1;
};
#endif

} // namespace Opm::Properties

namespace Opm
{
    /// Solver parameters for the BlackoilModel.
    template <class TypeTag>
    struct BlackoilModelParameters
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    public:
        /// Max relative change in bhp in single iteration.
        Scalar dbhp_max_rel_;
        /// Max absolute change in well volume fraction in single iteration.
        Scalar dwell_fraction_max_;
        /// Absolute max limit for residuals.
        Scalar max_residual_allowed_;
        //// Max allowed pore volume faction where CNV is violated. Below the
        //// relaxed tolerance tolerance_cnv_relaxed_ is used.
        Scalar relaxed_max_pv_fraction_;
        /// Relative mass balance tolerance (total mass balance error).
        Scalar tolerance_mb_;
        /// Relaxed mass balance tolerance (can be used when iter >= min_strict_mb_iter_).
        Scalar tolerance_mb_relaxed_;
        /// Local convergence tolerance (max of local saturation errors).
        Scalar tolerance_cnv_;
        /// Relaxed local convergence tolerance (can be used when iter >= min_strict_cnv_iter_ && cnvViolatedPV < relaxed_max_pv_fraction_).
        Scalar tolerance_cnv_relaxed_;
        /// Well convergence tolerance.
        Scalar tolerance_wells_;
        /// Tolerance for the well control equations
        //  TODO: it might need to distinguish between rate control and pressure control later
        Scalar tolerance_well_control_;
        /// Tolerance for the pressure equations for multisegment wells
        Scalar tolerance_pressure_ms_wells_;
        /// Relaxed tolerance for for the well flow residual
        Scalar relaxed_tolerance_flow_well_;

        /// Relaxed tolerance for the MSW pressure solution
        Scalar relaxed_tolerance_pressure_ms_well_;

        /// Maximum pressure change over an iteratio for ms wells
        Scalar max_pressure_change_ms_wells_;

        /// Maximum inner iteration number for ms wells
        int max_inner_iter_ms_wells_;

        /// Strict inner iteration number for wells
        int strict_inner_iter_wells_;

        /// Newton iteration where wells are stricly convergent
        int strict_outer_iter_wells_;

        /// Regularization factor for wells
        Scalar regularization_factor_wells_;

        /// Maximum newton iterations with inner well iterations
        int max_niter_inner_well_iter_;

        /// Whether to shut unsolvable well
        bool shut_unsolvable_wells_;

        /// Maximum inner iteration number for standard wells
        int max_inner_iter_wells_;

        /// Maximum iteration number of the well equation solution
        int max_welleq_iter_;

        /// Tolerance for time step in seconds where single precision can be used
        /// for solving for the Jacobian
        Scalar maxSinglePrecisionTimeStep_;

        /// Minimum number of Newton iterations before we can use relaxed CNV convergence criterion
        int min_strict_cnv_iter_;

        /// Minimum number of Newton iterations before we can use relaxed MB convergence criterion
        int min_strict_mb_iter_;

        /// Solve well equation initially
        bool solve_welleq_initially_;

        /// Update scaling factors for mass balance equations
        bool update_equations_scaling_;

        /// Try to detect oscillation or stagnation.
        bool use_update_stabilization_;

        /// Whether to use MultisegmentWell to handle multisegment wells
        /// it is something temporary before the multisegment well model is considered to be
        /// well developed and tested.
        /// if it is false, we will handle multisegment wells as standard wells, which will be
        /// the default behavoir for the moment. Later, we might set it to be true by default if necessary
        bool use_multisegment_well_;

        /// The file name of the deck
        std::string deck_file_name_;

        /// Whether to add influences of wells between cells to the matrix and preconditioner matrix
        bool matrix_add_well_contributions_;

        /// Whether to check well operability
        bool check_well_operability_;
        /// Whether to check well operability during iterations
        bool check_well_operability_iter_;

        /// Maximum number of times a well can switch to the same controt
        int max_number_of_well_switches_;

        /// Whether to approximate segment densities by averaging over segment and its outlet 
        bool use_average_density_ms_wells_;

        /// Whether to allow control switching during local well solutions 
        bool local_well_solver_control_switching_;

        /// Whether to use implicit IPR for thp stability checks and solution search
        bool use_implicit_ipr_;

        /// Maximum number of iterations in the network solver before relaxing tolerance
        int network_max_strict_iterations_;
        
        /// Maximum number of iterations in the network solver before giving up
        int network_max_iterations_;

        /// Nonlinear solver type: newton or nldd.
        std::string nonlinear_solver_;
        /// 'jacobi' and 'gauss-seidel' supported.
        DomainSolveApproach local_solve_approach_{DomainSolveApproach::Jacobi};

        int max_local_solve_iterations_;

        Scalar local_tolerance_scaling_mb_;
        Scalar local_tolerance_scaling_cnv_;

        int nldd_num_initial_newton_iter_{1};
        int num_local_domains_{0};
        Scalar local_domain_partition_imbalance_{1.03};
        std::string local_domain_partition_method_;
        DomainOrderingMeasure local_domain_ordering_{DomainOrderingMeasure::MaxPressure};

        bool write_partitions_{false};

        /// Construct from user parameters or defaults.
        BlackoilModelParameters()
        {
            dbhp_max_rel_=  Parameters::get<TypeTag, Properties::DbhpMaxRel>();
            dwell_fraction_max_ = Parameters::get<TypeTag, Properties::DwellFractionMax>();
            max_residual_allowed_ = Parameters::get<TypeTag, Properties::MaxResidualAllowed>();
            relaxed_max_pv_fraction_ = Parameters::get<TypeTag, Properties::RelaxedMaxPvFraction>();
            tolerance_mb_ = Parameters::get<TypeTag, Properties::ToleranceMb>();
            tolerance_mb_relaxed_ = std::max(tolerance_mb_, Parameters::get<TypeTag, Properties::ToleranceMbRelaxed>());
            tolerance_cnv_ = Parameters::get<TypeTag, Properties::ToleranceCnv>();
            tolerance_cnv_relaxed_ = std::max(tolerance_cnv_, Parameters::get<TypeTag, Properties::ToleranceCnvRelaxed>());
            tolerance_wells_ = Parameters::get<TypeTag, Properties::ToleranceWells>();
            tolerance_well_control_ = Parameters::get<TypeTag, Properties::ToleranceWellControl>();
            max_welleq_iter_ = Parameters::get<TypeTag, Properties::MaxWelleqIter>();
            use_multisegment_well_ = Parameters::get<TypeTag, Properties::UseMultisegmentWell>();
            tolerance_pressure_ms_wells_ = Parameters::get<TypeTag, Properties::TolerancePressureMsWells>();
            relaxed_tolerance_flow_well_ = Parameters::get<TypeTag, Properties::RelaxedWellFlowTol>();
            relaxed_tolerance_pressure_ms_well_ = Parameters::get<TypeTag, Properties::RelaxedPressureTolMsw>();
            max_pressure_change_ms_wells_ = Parameters::get<TypeTag, Properties::MaxPressureChangeMsWells>();
            max_inner_iter_ms_wells_ = Parameters::get<TypeTag, Properties::MaxInnerIterMsWells>();
            strict_inner_iter_wells_ = Parameters::get<TypeTag, Properties::StrictInnerIterWells>();
            strict_outer_iter_wells_ = Parameters::get<TypeTag, Properties::StrictOuterIterWells>();
            regularization_factor_wells_ = Parameters::get<TypeTag, Properties::RegularizationFactorWells>();
            max_niter_inner_well_iter_ = Parameters::get<TypeTag, Properties::MaxNewtonIterationsWithInnerWellIterations>();
            shut_unsolvable_wells_ = Parameters::get<TypeTag, Properties::ShutUnsolvableWells>();
            max_inner_iter_wells_ = Parameters::get<TypeTag, Properties::MaxInnerIterWells>();
            maxSinglePrecisionTimeStep_ = Parameters::get<TypeTag, Properties::MaxSinglePrecisionDays>() * 24 * 60 * 60;
            min_strict_cnv_iter_ = Parameters::get<TypeTag, Properties::MinStrictCnvIter>();
            min_strict_mb_iter_ = Parameters::get<TypeTag, Properties::MinStrictMbIter>();
            solve_welleq_initially_ = Parameters::get<TypeTag, Properties::SolveWelleqInitially>();
            update_equations_scaling_ = Parameters::get<TypeTag, Properties::UpdateEquationsScaling>();
            use_update_stabilization_ = Parameters::get<TypeTag, Properties::UseUpdateStabilization>();
            matrix_add_well_contributions_ = Parameters::get<TypeTag, Properties::MatrixAddWellContributions>();
            check_well_operability_ = Parameters::get<TypeTag, Properties::EnableWellOperabilityCheck>();
            check_well_operability_iter_ = Parameters::get<TypeTag, Properties::EnableWellOperabilityCheckIter>();
            max_number_of_well_switches_ = Parameters::get<TypeTag, Properties::MaximumNumberOfWellSwitches>();
            use_average_density_ms_wells_ = Parameters::get<TypeTag, Properties::UseAverageDensityMsWells>();
            local_well_solver_control_switching_ = Parameters::get<TypeTag, Properties::LocalWellSolveControlSwitching>();
            use_implicit_ipr_ = Parameters::get<TypeTag, Properties::UseImplicitIpr>();
            nonlinear_solver_ = Parameters::get<TypeTag, Properties::NonlinearSolver>();
            const auto approach = Parameters::get<TypeTag, Properties::LocalSolveApproach>();
            if (approach == "jacobi") {
                local_solve_approach_ = DomainSolveApproach::Jacobi;
            } else if (approach == "gauss-seidel") {
                local_solve_approach_ = DomainSolveApproach::GaussSeidel;
            } else {
                throw std::runtime_error("Invalid domain solver approach '" + approach + "' specified.");
            }

            max_local_solve_iterations_ = Parameters::get<TypeTag, Properties::MaxLocalSolveIterations>();
            local_tolerance_scaling_mb_ = Parameters::get<TypeTag, Properties::LocalToleranceScalingMb>();
            local_tolerance_scaling_cnv_ = Parameters::get<TypeTag, Properties::LocalToleranceScalingCnv>();
            nldd_num_initial_newton_iter_ = Parameters::get<TypeTag, Properties::NlddNumInitialNewtonIter>();
            num_local_domains_ = Parameters::get<TypeTag, Properties::NumLocalDomains>();
            local_domain_partition_imbalance_ = std::max(Scalar{1.0}, Parameters::get<TypeTag, Properties::LocalDomainsPartitioningImbalance>());
            local_domain_partition_method_ = Parameters::get<TypeTag, Properties::LocalDomainsPartitioningMethod>();
            deck_file_name_ = Parameters::get<TypeTag, Properties::EclDeckFileName>();
            network_max_strict_iterations_ = Parameters::get<TypeTag, Properties::NetworkMaxStrictIterations>();
            network_max_iterations_ = Parameters::get<TypeTag, Properties::NetworkMaxIterations>();
            local_domain_ordering_ = domainOrderingMeasureFromString(Parameters::get<TypeTag, Properties::LocalDomainsOrderingMeasure>());
            write_partitions_ = Parameters::get<TypeTag, Properties::DebugEmitCellPartition>();
        }

        static void registerParameters()
        {
            Parameters::registerParam<TypeTag, Properties::DbhpMaxRel>
                ("Maximum relative change of the bottom-hole pressure in a single iteration");
            Parameters::registerParam<TypeTag, Properties::DwellFractionMax>
                ("Maximum absolute change of a well's volume fraction in a single iteration");
            Parameters::registerParam<TypeTag, Properties::MaxResidualAllowed>
                ("Absolute maximum tolerated for residuals without cutting the time step size");
            Parameters::registerParam<TypeTag, Properties::RelaxedMaxPvFraction>
                ("The fraction of the pore volume of the reservoir "
                 "where the volumetric error (CNV) may be voilated "
                 "during strict Newton iterations.");
            Parameters::registerParam<TypeTag, Properties::ToleranceMb>
                ("Tolerated mass balance error relative to total mass present");
            Parameters::registerParam<TypeTag, Properties::ToleranceMbRelaxed>
                ("Relaxed tolerated mass balance error that applies for iterations "
                 "after the iterations with the strict tolerance");
            Parameters::registerParam<TypeTag, Properties::ToleranceCnv>
                ("Local convergence tolerance (Maximum of local saturation errors)");
            Parameters::registerParam<TypeTag, Properties::ToleranceCnvRelaxed>
                ("Relaxed local convergence tolerance that applies for iterations "
                 "after the iterations with the strict tolerance");
            Parameters::registerParam<TypeTag, Properties::ToleranceWells>
                ("Well convergence tolerance");
            Parameters::registerParam<TypeTag, Properties::ToleranceWellControl>
                ("Tolerance for the well control equations");
            Parameters::registerParam<TypeTag, Properties::MaxWelleqIter>
                ("Maximum number of iterations to determine solution the well equations");
            Parameters::registerParam<TypeTag, Properties::UseMultisegmentWell>
                ("Use the well model for multi-segment wells instead of the "
                 "one for single-segment wells");
            Parameters::registerParam<TypeTag, Properties::TolerancePressureMsWells>
                ("Tolerance for the pressure equations for multi-segment wells");
            Parameters::registerParam<TypeTag, Properties::RelaxedWellFlowTol>
                ("Relaxed tolerance for the well flow residual");
            Parameters::registerParam<TypeTag, Properties::RelaxedPressureTolMsw>
                ("Relaxed tolerance for the MSW pressure solution");
            Parameters::registerParam<TypeTag, Properties::MaxPressureChangeMsWells>
                ("Maximum relative pressure change for a single iteration "
                 "of the multi-segment well model");
            Parameters::registerParam<TypeTag, Properties::MaxInnerIterMsWells>
                ("Maximum number of inner iterations for multi-segment wells");
            Parameters::registerParam<TypeTag, Properties::StrictInnerIterWells>
                ("Number of inner well iterations with strict tolerance");
            Parameters::registerParam<TypeTag, Properties::StrictOuterIterWells>
                ("Number of newton iterations for which wells are checked with strict tolerance");
            Parameters::registerParam<TypeTag, Properties::MaxNewtonIterationsWithInnerWellIterations>
                ("Maximum newton iterations with inner well iterations");
            Parameters::registerParam<TypeTag, Properties::ShutUnsolvableWells>
                ("Shut unsolvable wells");
            Parameters::registerParam<TypeTag, Properties::MaxInnerIterWells>
                ("Maximum number of inner iterations for standard wells");
            Parameters::registerParam<TypeTag, Properties::AlternativeWellRateInit>
                ("Use alternative well rate initialization procedure");
            Parameters::registerParam<TypeTag, Properties::RegularizationFactorWells>
                ("Regularization factor for wells");
            Parameters::registerParam<TypeTag, Properties::MaxSinglePrecisionDays>
                ("Maximum time step size where single precision floating point "
                 "arithmetic can be used solving for the linear systems of equations");
            Parameters::registerParam<TypeTag, Properties::MinStrictCnvIter>
                ("Minimum number of Newton iterations before relaxed tolerances "
                 "can be used for the CNV convergence criterion");
            Parameters::registerParam<TypeTag, Properties::MinStrictMbIter>
                ("Minimum number of Newton iterations before relaxed tolerances "
                 "can be used for the MB convergence criterion. "
                 "Default -1 means that the relaxed tolerance is used when maximum "
                 "number of Newton iterations are reached.");
            Parameters::registerParam<TypeTag, Properties::SolveWelleqInitially>
                ("Fully solve the well equations before each iteration of the reservoir model");
            Parameters::registerParam<TypeTag, Properties::UpdateEquationsScaling>
                ("Update scaling factors for mass balance equations during the run");
            Parameters::registerParam<TypeTag, Properties::UseUpdateStabilization>
                ("Try to detect and correct oscillations or stagnation during the Newton method");
            Parameters::registerParam<TypeTag, Properties::MatrixAddWellContributions>
                ("Explicitly specify the influences of wells between cells in "
                 "the Jacobian and preconditioner matrices");
            Parameters::registerParam<TypeTag, Properties::EnableWellOperabilityCheck>
                ("Enable the well operability checking");
            Parameters::registerParam<TypeTag, Properties::EnableWellOperabilityCheckIter>
                ("Enable the well operability checking during iterations");
            Parameters::registerParam<TypeTag, Properties::MaximumNumberOfWellSwitches>
                ("Maximum number of times a well can switch to the same control");
            Parameters::registerParam<TypeTag, Properties::UseAverageDensityMsWells>
                ("Approximate segment densitities by averaging over segment and its outlet");
            Parameters::registerParam<TypeTag, Properties::LocalWellSolveControlSwitching>
                ("Allow control switching during local well solutions");
            Parameters::registerParam<TypeTag, Properties::UseImplicitIpr>
                ("Compute implict IPR for stability checks and stable solution search");
            Parameters::registerParam<TypeTag, Properties::NetworkMaxStrictIterations>
                ("Maximum iterations in network solver before relaxing tolerance");
            Parameters::registerParam<TypeTag, Properties::NetworkMaxIterations>
                ("Maximum number of iterations in the network solver before giving up");
            Parameters::registerParam<TypeTag, Properties::NonlinearSolver>
                ("Choose nonlinear solver. Valid choices are newton or nldd.");
            Parameters::registerParam<TypeTag, Properties::LocalSolveApproach>
                ("Choose local solve approach. Valid choices are jacobi and gauss-seidel");
            Parameters::registerParam<TypeTag, Properties::MaxLocalSolveIterations>
                ("Max iterations for local solves with NLDD nonlinear solver.");
            Parameters::registerParam<TypeTag, Properties::LocalToleranceScalingMb>
                ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
            Parameters::registerParam<TypeTag, Properties::LocalToleranceScalingCnv>
                ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
            Parameters::registerParam<TypeTag, Properties::NlddNumInitialNewtonIter>
                ("Number of initial global Newton iterations when running the NLDD nonlinear solver.");
            Parameters::registerParam<TypeTag, Properties::NumLocalDomains>
                ("Number of local domains for NLDD nonlinear solver.");
            Parameters::registerParam<TypeTag, Properties::LocalDomainsPartitioningImbalance>
                ("Subdomain partitioning imbalance tolerance. 1.03 is 3 percent imbalance.");
            Parameters::registerParam<TypeTag, Properties::LocalDomainsPartitioningMethod>
                ("Subdomain partitioning method. Allowed values are "
                 "'zoltan', "
                 "'simple', "
                 "and the name of a partition file ending with '.partition'.");
            Parameters::registerParam<TypeTag, Properties::LocalDomainsOrderingMeasure>
                ("Subdomain ordering measure. Allowed values are "
                 "'maxpressure', "
                 "'averagepressure' "
                 "and  'residual'.");
            Parameters::registerParam<TypeTag, Properties::DebugEmitCellPartition>
                ("Whether or not to emit cell partitions as a debugging aid.");

            Parameters::hideParam<TypeTag, Properties::DebugEmitCellPartition>();
        }
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED
