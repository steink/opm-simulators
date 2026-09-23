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


#ifndef OPM_AUTODIFF_VFPHELPERS_HPP_
#define OPM_AUTODIFF_VFPHELPERS_HPP_

#include <functional>
#include <map>
#include <vector>
#include <optional>

/**
 * This file contains a set of helper functions used by VFPProd / VFPInj.
 */
namespace Opm {

class VFPInjTable;
class VFPProdTable;

namespace detail {

/**
 * An "ADB-like" structure with a single value and a set of derivatives
 */
template<class Scalar>
struct VFPEvaluation
{
    VFPEvaluation() : value(0.0), dthp(0.0), dwfr(0.0), dgfr(0.0), dalq(0.0), dflo(0.0) {};
    Scalar value;
    Scalar dthp;
    Scalar dwfr;
    Scalar dgfr;
    Scalar dalq;
    Scalar dflo;
};

template<class Scalar>
VFPEvaluation<Scalar> operator+(VFPEvaluation<Scalar> lhs, const VFPEvaluation<Scalar>& rhs);
template<class Scalar>
VFPEvaluation<Scalar> operator-(VFPEvaluation<Scalar> lhs, const VFPEvaluation<Scalar>& rhs);
template<class Scalar>
VFPEvaluation<Scalar> operator*(Scalar lhs, const VFPEvaluation<Scalar>& rhs);

/**
 * What a slope-limited lookup (VFPHelpers::interpolateWithSlopeLimit(),
 * VFPProdProperties::bhp_with_slope_limit()) had to do to keep the tubing
 * curve flatter than the well's own IPR.
 *
 * A production tubing curve is U-shaped in FLO: falling on the low-rate,
 * liquid-loading branch, then rising once friction dominates. Where it falls
 * faster than the well's own IPR, the two can cross more than once and the
 * well has no stable operating point there. A slope-limited lookup replaces
 * that stretch with the linear extrapolation of the first interval flat
 * enough to guarantee a single crossing.
 */
enum class SlopeLimit
{
    //! Left alone: this interval is already flatter than the IPR, so the
    //! result is exactly what the plain lookup would have given.
    Unflattened,
    //! Flattened, with the steep run reaching the first FLO interval -- the
    //! ordinary U-shaped case, where the construction is continuous in FLO.
    Bridged,
    //! Flattened, but the steep run starts above the first interval: some
    //! interval below this one is flat, so the curve is not U-shaped and the
    //! flattened curve has a jump somewhere below this FLO. Still returned,
    //! since the local value is as good as in the Bridged case, but the
    //! caller should not trust a solve that converged here.
    NonMonotone,
    //! No interval up to the end of the FLO axis is flat enough; the slope is
    //! pinned at the limit itself. The table and this well's IPR cannot be
    //! reconciled anywhere on the axis.
    Clamped,
};

/**
 * A slope-limited lookup's value and derivatives, plus what it had to do.
 */
template<class Scalar>
struct SlopeLimitedEvaluation
{
    VFPEvaluation<Scalar> evaluation{};
    SlopeLimit limit{SlopeLimit::Unflattened};
};

/**
 * Helper struct for linear interpolation
 */
template<class Scalar>
struct InterpData
{
    InterpData() : ind_{0, 0}, inv_dist_(0.0), factor_(0.0) {}
    int ind_[2]; //[First element greater than or equal to value, Last element smaller than or equal to value]
    Scalar inv_dist_; // 1 / distance between the two end points of the segment. Used to calculate derivatives and uses 1.0 / 0.0 = 0.0 as a convention
    Scalar factor_; // Interpolation factor
};

/**
 * Computes the flo parameter according to the flo_type_
 * for production tables
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getFlo(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

/**
 * Computes the flo parameter according to the flo_type_
 * for injection tables
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getFlo(const VFPInjTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

/**
 * Computes the wfr parameter according to the wfr_type_
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getWFR(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

/**
 * Computes the gfr parameter according to the gfr_type_
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getGFR(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

/**
 * Returns the table from the map if found, or throws an exception
 */
template <typename T>
const T& getTable(const std::map<int, std::reference_wrapper<const T>>& tables, int table_id);

/**
 * Check whether we have a table with the table number
 */
template <typename T>
bool hasTable(const std::map<int, std::reference_wrapper<const T>>& tables, int table_id) {
    const auto entry = tables.find(table_id);
    return (entry != tables.end() );
}

/**
 * Returns the type variable for FLO/GFR/WFR for production tables
 */
template <typename TYPE, typename TABLE>
TYPE getType(const TABLE& table);

}

template<class Scalar>
class VFPHelpers {
public:
    /**
     * Helper function to find indices etc. for linear interpolation and extrapolation
     *  @param value_in Value to find in values
     *  @param values Sorted list of values to search for value in.
     *  @return Data required to find the interpolated value
     */
    static detail::InterpData<Scalar> findInterpData(const Scalar value_in,
                                                     const std::vector<double>& values);

    /**
     * Helper function which interpolates data using the indices etc. given in the inputs.
     */
    static detail::VFPEvaluation<Scalar> interpolate(const VFPProdTable& table,
                                                     const detail::InterpData<Scalar>& flo_i,
                                                     const detail::InterpData<Scalar>& thp_i,
                                                     const detail::InterpData<Scalar>& wfr_i,
                                                     const detail::InterpData<Scalar>& gfr_i,
                                                     const detail::InterpData<Scalar>& alq_i);

    /**
     * This basically models interpolate(VFPProdTable::array_type, ...)
     * which performs 5D interpolation, but here for the 2D case only
     */
    static detail::VFPEvaluation<Scalar> interpolate(const VFPInjTable& table,
                                                     const detail::InterpData<Scalar>& flo_i,
                                                     const detail::InterpData<Scalar>& thp_i);

    /**
     * interpolate() with the curve's own slope in FLO limited to \p max_slope
     * -- see detail::SlopeLimit for what that means and why.
     *
     * \p max_slope is d(bhp)/d(FLO) of the well's own IPR, oriented along the
     * table's own (positive, increasing) FLO axis and so negative for a
     * producer, with whatever safety margin the caller wants already folded
     * in: a curve exactly parallel to the IPR still admits no unique
     * crossing, so the caller should pass a slightly *less* negative value
     * than the bare IPR slope.
     *
     * Where the interval containing \p flo is already flatter than that, the
     * result is exactly interpolate()'s, flagged Unflattened. Where it is
     * not, the walk moves up the FLO axis to the first interval that is, and
     * returns that interval's linear extrapolation back to \p flo -- value
     * and all five partials from the same line, so the pair stays consistent
     * (unlike clipping the derivative alone, which leaves the linearisation
     * describing a different function than the value does).
     *
     * \p flo is the query FLO oriented along the table's own axis (positive),
     * i.e. already sign-flipped from OPM's negative-for-production rates, and
     * must be the same value \p flo_i was found for.
     */
    static detail::SlopeLimitedEvaluation<Scalar>
    interpolateWithSlopeLimit(const VFPProdTable& table,
                              const Scalar flo,
                              const detail::InterpData<Scalar>& flo_i,
                              const detail::InterpData<Scalar>& thp_i,
                              const detail::InterpData<Scalar>& wfr_i,
                              const detail::InterpData<Scalar>& gfr_i,
                              const detail::InterpData<Scalar>& alq_i,
                              const Scalar max_slope);

    static detail::VFPEvaluation<Scalar> bhp(const VFPProdTable& table,
                                             const Scalar aqua,
                                             const Scalar liquid,
                                             const Scalar vapour,
                                             const Scalar thp,
                                             const Scalar alq,
                                             const Scalar explicit_wfr,
                                             const Scalar explicit_gfr,
                                             const bool   use_vfpexplicit);

    static detail::VFPEvaluation<Scalar> bhp(const VFPInjTable& table,
                                             const Scalar aqua,
                                             const Scalar liquid,
                                             const Scalar vapour,
                                             const Scalar thp);

    /**
     * This function finds the value of THP given a specific BHP.
     * Essentially:
     *   Given the function f(thp_array(x)) = bhp_array(x), which is piecewise linear,
     *   find thp so that f(thp) = bhp.
     */
    static Scalar findTHP(const std::vector<Scalar>& bhp_array,
                          const std::vector<double>& thp_array,
                          Scalar bhp,
                          const bool find_largest = true);

    /**
    * Get (flo, bhp) at minimum bhp for given thp,wfr,gfr,alq
    */
    static std::pair<Scalar, Scalar>
    getMinimumBHPCoordinate(const VFPProdTable& table,
                            const Scalar thp,
                            const Scalar wfr,
                            const Scalar gfr,
                            const Scalar alq);

    /**
    * Get (flo, bhp) at largest occuring stable vfp/ipr-intersection
    * if it exists
    */
    static std::optional<std::pair<Scalar, Scalar>>
    intersectWithIPR(const VFPProdTable& table,
                     const Scalar thp,
                     const Scalar wfr,
                     const Scalar gfr,
                     const Scalar alq,
                     const Scalar ipr_a,
                     const Scalar ipr_b,
                     const std::function<Scalar(const Scalar)>& adjust_bhp);
};

} // namespace

#endif /* OPM_AUTODIFF_VFPHELPERS_HPP_ */
