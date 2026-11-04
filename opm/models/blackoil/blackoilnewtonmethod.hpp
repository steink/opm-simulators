// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::BlackOilNewtonMethod
 */
#ifndef OPM_BLACK_OIL_NEWTON_METHOD_HPP
#define OPM_BLACK_OIL_NEWTON_METHOD_HPP

#include <opm/common/Exceptions.hpp>

#include <opm/models/blackoil/blackoilmodules.hpp>
#include <opm/models/blackoil/blackoilnewtonmethodparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/nonlinear/newtonmethod.hh>

#include <opm/models/parallel/gridcommhandles.hh>

#include <opm/models/utils/signum.hh>

#include <opm/material/common/Valgrind.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <vector>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief A newton solver which is specific to the black oil model.
 */
template <class TypeTag>
class BlackOilNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    static constexpr bool enableBioeffects = getPropValue<TypeTag, Properties::EnableBioeffects>();
    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag, enableBioeffects>;

    static const unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();

public:
    explicit BlackOilNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {
        bparams_.read();
    }

    /*!
     * \copydoc NewtonMethod::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        wasSwitched_.resize(this->model().numTotalDof(), false);
        keepRsEvent_.resize(this->model().numTotalDof(), KeepRsEvent::None);
        appearingRs_.resize(this->model().numTotalDof(), noRs);
    }

    void capturePrimaryVariableSwitches()
    {
        // Preserve the switching history across retries of this timestep.
        wasSwitchedAtStepStart_ = wasSwitched_;
    }

    void restorePrimaryVariableSwitches()
    {
        wasSwitched_ = wasSwitchedAtStepStart_;
    }

    /*!
     * \brief Register all run-time parameters for the blackoil newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        BlackoilNewtonParams<Scalar>::registerParameters();
    }

    /*!
     * \brief Returns the number of degrees of freedom for which the
     *        interpretation has changed for the most recent iteration.
     */
    unsigned numPriVarsSwitched() const
    { return numPriVarsSwitched_; }

    /// DEBUG: what keepRsOnSwitch_() last did for a cell
    enum class KeepRsEvent : std::uint8_t {
        None,      // nothing yet
        Inflow,    // oil appeared: Rs from the inflowing oil only
        OwnOil,    // oil appeared: Rs from the start-of-step oil only
        Mixed,     // oil appeared: Rs from both
        RsSat,     // oil appeared: no oil in place or inflow, RsSat kept
    };

    KeepRsEvent keepRsEvent(const unsigned globalDofIdx) const
    { return keepRsEvent_.empty() ? KeepRsEvent::None : keepRsEvent_[globalDofIdx]; }

protected:
    friend NewtonMethod<TypeTag>;
    friend ParentType;

    /*!
     * \copydoc FvBaseNewtonMethod::beginIteration_
     */
    void beginIteration_()
    {
        numPriVarsSwitched_ = 0;
        ParentType::beginIteration_();
    }

    /*!
     * \copydoc FvBaseNewtonMethod::endIteration_
     *
     * \param uCurrentIter Current solution iterator
     * \param uLastIter Last solution iterator
     */
    void endIteration_(SolutionVector& uCurrentIter,
                       const SolutionVector& uLastIter)
    {
#if HAVE_MPI
        // in the MPI enabled case we need to add up the number of DOF
        // for which the interpretation changed over all processes.
        const int localSwitched = numPriVarsSwitched_;
        MPI_Allreduce(&localSwitched,
                      &numPriVarsSwitched_,
                      /*num=*/1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);
#endif // HAVE_MPI

        this->simulator_.model().newtonMethod().endIterMsg()
            << ", num switched=" << numPriVarsSwitched_;

        ParentType::endIteration_(uCurrentIter, uLastIter);
    }

public:
    /// Whether the primary variables of a degree of freedom were switched in the
    /// current time step, which makes switching them back harder.
    bool wasSwitched(const unsigned globalDofIdx) const
    { return wasSwitched_[globalDofIdx]; }

    void setWasSwitched(const unsigned globalDofIdx, const bool switched)
    { wasSwitched_[globalDofIdx] = switched; }

    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        const auto& comm = this->simulator_.gridView().comm();
        std::ranges::fill(appearingRs_, noRs);

        int succeeded;
        try {
            ParentType::update_(nextSolution,
                                currentSolution,
                                solutionUpdate,
                                currentResidual);
            succeeded = 1;
        }
        catch (...) {
            succeeded = 0;
        }
        succeeded = comm.min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("A process did not succeed in adapting the primary variables");
        }

        applyAppearingRs_(nextSolution);

        numPriVarsSwitched_ = comm.sum(numPriVarsSwitched_);
    }

    template <class DofIndices>
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual,
                 const DofIndices& dofIndices)
    {
        const auto zero = 0.0 * solutionUpdate[0];
        for (auto dofIdx : dofIndices) {
            if (solutionUpdate[dofIdx] == zero) {
                continue;
            }
            appearingRs_[dofIdx] = noRs;
            updatePrimaryVariables_(dofIdx,
                                    nextSolution[dofIdx],
                                    currentSolution[dofIdx],
                                    solutionUpdate[dofIdx],
                                    currentResidual[dofIdx]);
            if constexpr (Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max()) {
                // NLDD subdomain solve: no ghost cells to synchronise
                if (appearingRs_[dofIdx] != noRs) {
                    auto& rs = nextSolution[dofIdx][Indices::compositionSwitchIdx];
                    rs = std::min(rs, appearingRs_[dofIdx]);
                }
            }
        }
    }

protected:
    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        static constexpr bool enableSolvent =
            Indices::solventSaturationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableExtbo =
            Indices::zFractionIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enablePolymer =
            Indices::polymerConcentrationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enablePolymerWeight =
            Indices::polymerMoleWeightIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableFullyImplicitThermal =
            Indices::temperatureIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableFoam =
            Indices::foamConcentrationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableBrine =
            Indices::saltConcentrationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableMICP = Indices::enableMICP;

        currentValue.checkDefined();
        Valgrind::CheckDefined(update);
        Valgrind::CheckDefined(currentResidual);

        // DEBUG (OPM_DEBUG_KEEP_RS, OPM_DEBUG_RS_CHOP): values before the update.
        // (nextValue and currentValue may be the same object, so copy.)
        const PrimaryVariables beforeUpdate = currentValue;

        // saturation delta for each phase
        Scalar deltaSw = 0.0;
        Scalar deltaSo = 0.0;
        Scalar deltaSg = 0.0;
        Scalar deltaSs = 0.0;

        if (currentValue.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw)
        {
            if constexpr (Indices::waterSwitchIdx != std::numeric_limits<unsigned>::max()) {
                deltaSw = update[Indices::waterSwitchIdx];
                deltaSo -= deltaSw;
            }
        }
        if (currentValue.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg)
        {
            if constexpr (Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max()) {
                deltaSg = update[Indices::compositionSwitchIdx];
                deltaSo -= deltaSg;
            }
        }
        if (currentValue.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
            if constexpr (Indices::solventSaturationIdx != std::numeric_limits<unsigned>::max()) {
                deltaSs = update[Indices::solventSaturationIdx];
                deltaSo -= deltaSs;
            }
        }

        // maximum saturation delta
        Scalar maxSatDelta = std::max(std::abs(deltaSg), std::abs(deltaSo));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSw));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSs));

        // scaling factor for saturation deltas to make sure that none of them exceeds
        // the specified threshold value.
        Scalar satAlpha = 1.0;
        if (maxSatDelta > bparams_.dsMax_) {
            satAlpha = bparams_.dsMax_ / maxSatDelta;
        }

        for (unsigned pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            // calculate the update of the current primary variable. For the black-oil
            // model we limit the pressure delta relative to the pressure's current
            // absolute value (Default: 30%) and saturation deltas to an absolute change
            // (Default: 20%). Further, we ensure that the R factors, solvent
            // "saturation" and polymer concentration do not become negative after the
            // update.
            Scalar delta = update[pvIdx];

            // limit pressure delta
            if (pvIdx == Indices::pressureSwitchIdx) {
                if (std::abs(delta) > bparams_.dpMaxRel_ * currentValue[pvIdx]) {
                    delta = signum(delta) * bparams_.dpMaxRel_ * currentValue[pvIdx];
                }
            }
            // water saturation delta
            else if (pvIdx == Indices::waterSwitchIdx)
                if (currentValue.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                    delta *= satAlpha;
                }
                else {
                    //Ensure Rvw and Rsw factor does not become negative
                    if (delta > currentValue[ Indices::waterSwitchIdx]) {
                        delta = currentValue[ Indices::waterSwitchIdx];
                    }
                }
            else if (pvIdx == Indices::compositionSwitchIdx) {
                // the switching primary variable for composition is tricky because the
                // "reasonable" value ranges it exhibits vary widely depending on its
                // interpretation since it can represent Sg, Rs or Rv. For now, we only
                // limit saturation deltas and ensure that the R factors do not become
                // negative.
                if (currentValue.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg) {
                    delta *= satAlpha;
                }
                else {
                    // Ensure Rv and Rs factor does not become negative
                    if (delta > currentValue[Indices::compositionSwitchIdx]) {
                        delta = currentValue[Indices::compositionSwitchIdx];
                    }
                }
            }
            else if (enableSolvent && pvIdx == Indices::solventSaturationIdx) {
                // solvent saturation updates are also subject to the Appleyard chop
                if (currentValue.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
                    delta *= satAlpha;
                }
                else {
                    // Ensure Rssolw factor does not become negative
                    if (delta > currentValue[Indices::solventSaturationIdx]) {
                        delta = currentValue[Indices::solventSaturationIdx];
                    }
                }
            }
            else if (enableExtbo && pvIdx == Indices::zFractionIdx) {
                // z fraction updates are also subject to the Appleyard chop
                const auto& curr = currentValue[Indices::zFractionIdx]; // or currentValue[pvIdx] given the block condition
                delta = std::clamp(delta, curr - Scalar{1.0}, curr);
            }
            else if (enablePolymerWeight && pvIdx == Indices::polymerMoleWeightIdx) {
                const double sign = delta >= 0. ? 1. : -1.;
                // maximum change of polymer molecular weight, the unit is MDa.
                // applying this limit to stabilize the simulation. The value itself is still experimental.
                const Scalar maxMolarWeightChange = 100.0;
                delta = sign * std::min(std::abs(delta), maxMolarWeightChange);
                delta *= satAlpha;
            }
            else if (enableFullyImplicitThermal && pvIdx == Indices::temperatureIdx) {
                const double sign = delta >= 0. ? 1. : -1.;
                delta = sign * std::min(std::abs(delta), bparams_.maxTempChange_);
            }
            else if (enableBrine && pvIdx == Indices::saltConcentrationIdx &&
                     enableSaltPrecipitation &&
                     currentValue.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp)
            {
                const Scalar maxSaltSaturationChange = 0.1;
                const Scalar sign = delta >= 0. ? 1. : -1.;
                delta = sign * std::min(std::abs(delta), maxSaltSaturationChange);
            }

            // do the actual update
            nextValue[pvIdx] = currentValue[pvIdx] - delta;

            // keep the solvent saturation between 0 and 1
            if (enableSolvent && pvIdx == Indices::solventSaturationIdx) {
                if (currentValue.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
                    nextValue[pvIdx] = std::min(std::max(nextValue[pvIdx], Scalar{0.0}), Scalar{1.0});
                }
            }

            // keep the z fraction between 0 and 1
            if (enableExtbo && pvIdx == Indices::zFractionIdx) {
                nextValue[pvIdx] = std::min(std::max(nextValue[pvIdx], Scalar{0.0}), Scalar{1.0});
            }

            // keep the polymer concentration above 0
            if (enablePolymer && pvIdx == Indices::polymerConcentrationIdx) {
                nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
            }

            if (enablePolymerWeight && pvIdx == Indices::polymerMoleWeightIdx) {
                nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                const double polymerConcentration = nextValue[Indices::polymerConcentrationIdx];
                if (polymerConcentration < 1.e-10) {
                    nextValue[pvIdx] = 0.0;
                }
            }

            // keep the foam concentration above 0
            if (enableFoam && pvIdx == Indices::foamConcentrationIdx) {
                nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
            }

            if (enableBrine && pvIdx == Indices::saltConcentrationIdx) {
               // keep the salt concentration above 0
                if (!enableSaltPrecipitation ||
                    currentValue.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Cs)
               {
                   nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                }
               // keep the salt saturation below upperlimit
                if (enableSaltPrecipitation &&
                    currentValue.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp)
                {
                   nextValue[pvIdx] = std::min(nextValue[pvIdx], Scalar{1.0-1.e-8});
                }
            }

            // keep the temperature within given values
            if (enableFullyImplicitThermal && pvIdx == Indices::temperatureIdx) {
                nextValue[pvIdx] = std::clamp(nextValue[pvIdx], bparams_.tempMin_, bparams_.tempMax_);
            }

            if (pvIdx == Indices::pressureSwitchIdx) {
                nextValue[pvIdx] = std::clamp(nextValue[pvIdx], bparams_.pressMin_, bparams_.pressMax_);
            }

            // keep the values above 0
            // for the biofilm and calcite, we set an upper limit equal to the initial porosity
            // minus 1e-8. This prevents singularities (e.g., one of the calcite source term is
            // evaluated at 1/(iniPoro - calcite)). The value 1e-8 is taken from the salt precipitation
            // clapping above.
            if constexpr (enableBioeffects) {
                if (pvIdx == Indices::microbialConcentrationIdx) {
                    nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                }
                if (pvIdx == Indices::biofilmVolumeFractionIdx) {
                    nextValue[pvIdx] = std::clamp(nextValue[pvIdx],
                                                  Scalar{0.0},
                                                  this->problem().referencePorosity(globalDofIdx, 0) - 1e-8);
                }
                if constexpr (enableMICP) {
                    if (pvIdx == Indices::oxygenConcentrationIdx) {
                        nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                    }
                    if (pvIdx == Indices::ureaConcentrationIdx) {
                        nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                    }
                    if (pvIdx == Indices::calciteVolumeFractionIdx) {
                        nextValue[pvIdx] = std::clamp(nextValue[pvIdx], Scalar{0.0},
                                                                        this->problem().referencePorosity(globalDofIdx, 0) - 1e-8);
                    }
                }
            }
        }

        // DEBUG (OPM_DEBUG_KEEP_RS): see keepRsOnSwitch_().
        constexpr bool keepRsPossible =
            Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max()
            && Indices::waterEnabled;
        if constexpr (keepRsPossible) {
            chopRsUpdate_(nextValue, beforeUpdate);
        }

        // switch the new primary variables to something which is physically meaningful.
        // use a threshold value after a switch to make it harder to switch back
        // immediately.
        if (wasSwitched_[globalDofIdx]) {
            wasSwitched_[globalDofIdx] = nextValue.adaptPrimaryVariables(this->problem(),
                                                                         globalDofIdx,
                                                                         bparams_.waterSaturationMax_,
                                                                         bparams_.waterOnlyThreshold_,
                                                                         bparams_.priVarOscilationThreshold_);
        }
        else {
            wasSwitched_[globalDofIdx] = nextValue.adaptPrimaryVariables(this->problem(),
                                                                         globalDofIdx,
                                                                         bparams_.waterSaturationMax_,
                                                                         bparams_.waterOnlyThreshold_);
        }

        if constexpr (keepRsPossible) {
            keepRsOnSwitch_(globalDofIdx, nextValue, beforeUpdate);
        }

        if (wasSwitched_[globalDofIdx]) {
            ++numPriVarsSwitched_;
        }
        if (bparams_.projectSaturations_) {
            nextValue.chopAndNormalizeSaturations();
        }

        nextValue.checkDefined();
    }

private:
    /// DEBUG: give oil appearing in a water-only cell the Rs of that oil.
    ///
    /// A cell with Sw >= WaterOnlyThreshold is in the Sg state with Sg = 0, which has
    /// no Rs. When oil appears there (Sg -> Rs) adaptPrimaryVariables() sets Rs to
    /// RsSat, which is only right if free gas was present. With OPM_DEBUG_KEEP_RS
    /// set, Rs is instead the gas/oil ratio of the oil the cell should hold,
    /// see appearingOilRs_(), capped at RsSat. Gas disappearing from a cell that
    /// has oil keeps RsSat. The value is only recorded here. It is applied by
    /// applyAppearingRs_(), as ghost cells have to get it from their owner.
    void keepRsOnSwitch_(const unsigned globalDofIdx,
                         const PrimaryVariables& nextValue,
                         const PrimaryVariables& beforeUpdate)
    {
        static const bool enabled = std::getenv("OPM_DEBUG_KEEP_RS") != nullptr;
        if (!enabled) {
            return;
        }

        using GasMeaning = typename PrimaryVariables::GasMeaning;
        if (beforeUpdate.primaryVarsMeaningGas() != GasMeaning::Sg
            || nextValue.primaryVarsMeaningGas() != GasMeaning::Rs
            || !hasNoOil_(beforeUpdate)) {
            return;
        }

        KeepRsEvent event = KeepRsEvent::RsSat;
        const Scalar rs = appearingOilRs_(globalDofIdx, event);
        if (event != KeepRsEvent::RsSat) {
            appearingRs_[globalDofIdx] = rs;
        }
        keepRsEvent_[globalDofIdx] = event;
    }

    /// DEBUG: apply the Rs recorded by keepRsOnSwitch_(). A ghost cell has an
    /// incomplete stencil, so its own inflow is wrong: it takes its owner's value.
    void applyAppearingRs_(SolutionVector& nextSolution)
    {
        if constexpr (Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max()
                      && Indices::waterEnabled) {
            static const bool enabled = std::getenv("OPM_DEBUG_KEEP_RS") != nullptr;
            if (!enabled) {
                return;
            }
            const auto& gridView = this->simulator_.gridView();
            if (gridView.comm().size() > 1) {
                const auto& mapper = this->model().elementMapper();
                using Handle = GridCommHandleGhostSync<Scalar, std::vector<Scalar>,
                                                       std::decay_t<decltype(mapper)>, /*commCodim=*/0>;
                Handle handle(appearingRs_, mapper);
                gridView.communicate(handle, Dune::InteriorBorder_All_Interface,
                                     Dune::ForwardCommunication);
            }
            const auto numGridDof = this->model().numGridDof();
            for (std::size_t dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                if (appearingRs_[dofIdx] != noRs) {
                    // the value here is min(RsMax, RsSat)
                    auto& rs = nextSolution[dofIdx][Indices::compositionSwitchIdx];
                    rs = std::min(rs, appearingRs_[dofIdx]);
                }
            }
        }
    }

    /// DEBUG: whether primary variables describe a cell without oil, i.e. the
    /// water-only state (Sg state with Sw + Sg = 1).
    static bool hasNoOil_(const PrimaryVariables& pv)
    {
        using GasMeaning = typename PrimaryVariables::GasMeaning;
        using WaterMeaning = typename PrimaryVariables::WaterMeaning;
        return pv.primaryVarsMeaningGas() == GasMeaning::Sg
            && pv.primaryVarsMeaningWater() == WaterMeaning::Sw
            && 1.0 - pv[Indices::waterSwitchIdx] - pv[Indices::compositionSwitchIdx] <= 1e-6;
    }

    /// DEBUG: gas/oil ratio of the oil a cell should hold at the end of the step:
    ///
    ///   Rs = (O_old * Rs_old + dt * sum q_o * Rs_up) / (O_old + dt * sum q_o)
    ///
    /// with O_old the oil in place at the start of the step (surface volume, using
    /// the current bo, as it only sets a weight) and q_o the oil inflow over each face
    /// from the current iterate's intensive quantities. A cell with free gas at the
    /// start of the step has saturated oil, so nothing is computed for it. The O_old
    /// term is left out if OPM_DEBUG_KEEP_RS_NO_OLD is set. Sets event to RsSat if
    /// there is neither oil in place nor inflow, or without per-face neighbour
    /// information (non-TPFA).
    Scalar appearingOilRs_(const unsigned globalDofIdx, KeepRsEvent& event) const
    {
        static const bool useOldOil = std::getenv("OPM_DEBUG_KEEP_RS_NO_OLD") == nullptr;

        using GasMeaning = typename PrimaryVariables::GasMeaning;
        using WaterMeaning = typename PrimaryVariables::WaterMeaning;
        const auto& model = this->model();
        const auto& stepStart = model.solution(/*timeIdx=*/1)[globalDofIdx];
        const bool startInRsState = stepStart.primaryVarsMeaningGas() == GasMeaning::Rs
            && stepStart.primaryVarsMeaningWater() == WaterMeaning::Sw;
        if (!startInRsState && !hasNoOil_(stepStart)) {
            event = KeepRsEvent::RsSat;
            return 0.0;
        }

        if constexpr (requires { model.linearizer().getNeighborInfo(); }) {
            using RateVector = GetPropType<TypeTag, Properties::RateVector>;
            using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
            constexpr auto oilIdx = FluidSystem::oilPhaseIdx;
            constexpr auto rsIdx = Indices::compositionSwitchIdx;
            const unsigned oilEq = FluidSystem::canonicalToActiveCompIdx(
                FluidSystem::solventComponentIndex(oilIdx));
            const auto& intQuants = model.intensiveQuantities(globalDofIdx, /*timeIdx=*/0);
            const Scalar dt = this->simulator_.timeStepSize();

            Scalar oil = 0.0;
            Scalar gas = 0.0;
            if (useOldOil && startInRsState) {
                const Scalar soOld = std::max(1.0 - stepStart[Indices::waterSwitchIdx], 0.0);
                oil = soOld * getValue(intQuants.fluidState().invB(oilIdx))
                    * this->problem().referencePorosity(globalDofIdx, /*timeIdx=*/0)
                    * model.dofTotalVolume(globalDofIdx);
                gas = oil * stepStart[rsIdx];
            }
            const bool haveOldOil = oil > 0.0;

            bool haveInflow = false;
            for (const auto& nb : model.linearizer().getNeighborInfo()[globalDofIdx]) {
                const auto& intQuantsNb = model.intensiveQuantities(nb.neighbor, /*timeIdx=*/0);
                RateVector flux(0.0);
                RateVector darcy(0.0);
                LocalResidual::computeFlux(flux, darcy, globalDofIdx, nb.neighbor, intQuants,
                                           intQuantsNb, nb.res_nbinfo,
                                           this->problem().moduleParams());
                // > 0: out of the cell
                const Scalar qO = getValue(flux[oilEq]) * nb.res_nbinfo.faceArea;
                if (qO < 0.0) {
                    haveInflow = true;
                    oil -= dt * qO;
                    gas -= dt * qO * getValue(intQuantsNb.fluidState().Rs());
                }
            }
            if (oil > 0.0) {
                event = !haveInflow ? KeepRsEvent::OwnOil
                    : haveOldOil ? KeepRsEvent::Mixed : KeepRsEvent::Inflow;
                return gas / oil;
            }
        }
        event = KeepRsEvent::RsSat;
        return 0.0;
    }

    static constexpr Scalar noRs = -1.0;

    /// DEBUG: limit the Rs change per Newton iteration in the Rs state to
    /// OPM_DEBUG_RS_CHOP (in the model's Rs unit, sm3/sm3). Off if not set.
    void chopRsUpdate_(PrimaryVariables& nextValue,
                       const PrimaryVariables& beforeUpdate) const
    {
        static const Scalar maxDrs = [] {
            const char* env = std::getenv("OPM_DEBUG_RS_CHOP");
            return env ? static_cast<Scalar>(std::atof(env)) : Scalar{0.0};
        }();
        if (maxDrs <= 0.0
            || beforeUpdate.primaryVarsMeaningGas() != PrimaryVariables::GasMeaning::Rs) {
            return;
        }
        constexpr auto rsIdx = Indices::compositionSwitchIdx;
        const Scalar drs = std::clamp(nextValue[rsIdx] - beforeUpdate[rsIdx], -maxDrs, maxDrs);
        nextValue[rsIdx] = std::max(beforeUpdate[rsIdx] + drs, Scalar{0.0});
    }

    int numPriVarsSwitched_{};

    BlackoilNewtonParams<Scalar> bparams_{};

    // keep track of cells where the primary variable meaning has changed
    // to detect and hinder oscillations
    std::vector<bool> wasSwitched_{};
    std::vector<bool> wasSwitchedAtStepStart_{};

    // DEBUG: what keepRsOnSwitch_() last did, for the worst-CNV print
    std::vector<KeepRsEvent> keepRsEvent_{};
    // DEBUG: Rs for oil appearing in this update, noRs elsewhere
    std::vector<Scalar> appearingRs_{};
};

} // namespace Opm

#endif // OPM_BLACK_OIL_NEWTHON_METHOD_HPP
