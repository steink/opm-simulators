/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_PROD_GROUP_TREE_NODE_HEADER_INCLUDED
#define OPM_PROD_GROUP_TREE_NODE_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>

#include <array>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace Opm {

/// Represents the control status of a node in the production group tree.
/// Mirrors Matlab mode_cnt values.
enum class ProdNodeCtrlStatus {
    Undetermined,        ///< Top-node: not yet resolved (Matlab: UNDETERMINED)
    GroupControlled,     ///< Node is controlled by its parent group (Matlab: GRUP)
    IndividualControlled,///< Node has hit an individual limit and is now fixed (Matlab: INDIVIDUAL)
    NoGroupChildren      ///< Group node whose GroupControlled children are all resolved (Matlab: NONE)
};

/// Node type: well leaf vs. group interior node.
enum class ProdNodeType {
    Well,
    Group
};

/// One node in the production group tree used by the balancing algorithm.
/// The algorithm always uses 3-component [oil, water, gas] rate vectors;
/// missing phases are set to zero.
/// Rates follow the sign convention: negative = production.
template<class Scalar>
struct ProdGroupTreeNode {
    // ---- Identity --------------------------------------------------------
    std::string name;
    ProdNodeType type{ProdNodeType::Well};

    // ---- Tree topology ---------------------------------------------------
    std::string parent;                ///< parent node name; empty for FIELD
    std::vector<std::string> children; ///< child node names

    // ---- Group control status --------------------------------------------
    bool allowGroupControl{true};      ///< well: isAvailableForGroupControl()
                                       ///< group: productionGroupControlAvailable()

    ProdNodeCtrlStatus ctrlStatus{ProdNodeCtrlStatus::Undetermined};

    /// Which individual control mode is currently active / most limiting.
    Well::ProducerCMode activeIndividualCtrl{Well::ProducerCMode::CMODE_UNDEFINED};

    /// Preferred control for top-level UNDETERMINED nodes
    Well::ProducerCMode preferredCtrl{Well::ProducerCMode::CMODE_UNDEFINED};

    // ---- Rates -----------------------------------------------------------
    /// [oil, water, gas] surface rates; negative = production.
    std::array<Scalar, 3> rates{};

    // ---- Individual limits -----------------------------------------------
    /// Maps each active limit type to the scalar value the algorithm compares
    /// against the corresponding projection of the node's rates.
    ///  ORAT  → oil limit (compared to -rates[0])
    ///  WRAT  → water limit (compared to -rates[1])
    ///  GRAT  → gas limit (compared to -rates[2])
    ///  LRAT  → liquid limit (compared to -(rates[0]+rates[1]))
    ///  RESV  → reservoir-volume limit (compared to sum_p -rates[p]*resv_coeff[p])
    ///  BHP   → equivalent total rate derived from IPR at the BHP limit
    ///          (for wells only; groups do not have a BHP limit)
    std::map<Well::ProducerCMode, Scalar> individualLimits;

    // ---- Guide rates (for wells only, used to split group targets) -------
    /// [oil, water, gas] guide rates; zero for groups.
    std::array<Scalar, 3> guideRates{};

    /// Fixed phase fractions of the well's current rates (for wells only).
    /// fractions[i] = -rates[i] / total_rate, where total_rate = sum(-rates).
    std::array<Scalar, 3> wellRateFractions{};

    // ---- Parametrization (updated during the algorithm) ------------------
    /// The rate change per unit alpha-step: delta_rates = alpha * linearTerm.
    std::array<Scalar, 3> linearTerm{};

    /// Distance (in alpha units) to the next individual limit.
    Scalar alphaToNextLimit{std::numeric_limits<Scalar>::max()};

    /// Which control mode will be hit at alphaToNextLimit.
    Well::ProducerCMode nextLimitCtrl{Well::ProducerCMode::CMODE_UNDEFINED};

    // ---- Group target (set by parent during target assignment) -----------
    struct GroupTarget {
        std::string groupName;
        Group::ProductionCMode ctrlMode{Group::ProductionCMode::NONE};
        Scalar value{0};
        Scalar guideRate{0};
        Scalar guideRateRatio{1};
    };
    GroupTarget groupTarget;

    // ---- RESV conversion coefficients (wells only) -----------------------
    /// convert_coeff[i] for each active phase so that
    ///   resv_rate = sum_p -rates[p] * resvCoeff[p].
    /// Stored in canonical [oil, water, gas] order; zero for inactive phases.
    std::array<Scalar, 3> resvCoeff{};
};

} // namespace Opm

#endif // OPM_PROD_GROUP_TREE_NODE_HEADER_INCLUDED
