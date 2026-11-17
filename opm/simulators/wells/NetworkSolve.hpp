/*
  Copyright 2026 Equinor ASA.

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
#ifndef OPM_NETWORK_SOLVE_HEADER_INCLUDED
#define OPM_NETWORK_SOLVE_HEADER_INCLUDED

#include <fmt/format.h>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace Opm::NetworkSolve {

constexpr int NoTable = 9999;   // GNETINJE's "no table": pressure passes through

struct Node
{
    std::string name;
    int parent = -1;            // -1 only for the terminal
    int vfp_table = NoTable;
    /// NEFAC: what this node passes on of what it collects.
    double efficiency = 1.0;
};

template<class Scalar>
struct Result
{
    bool converged = false;
    int iterations = 0;
    std::vector<Scalar> node_pressure;   // every node, terminal included
    std::vector<Scalar> well_rate;       // per well, in the order they were added

    /// Why it stopped, for the caller to report. A converged solve leaves these
    /// at the values that satisfied the test.
    Scalar residual = 0.0;          // max norm of the scaled residual
    bool controls_moving = false;   // an active-set change on the last iteration
    /// One letter per well per iteration for the last few iterations, so a
    /// cycling active set can be read off (SystemBase::controlLetter()).
    std::string control_trace;
    /// Iterations on which some well changed control. A solve that has to move
    /// the active set a few times is working; one that keeps moving it is not.
    int switches = 0;
    /// Every well's phase rates and bhp at the solution.
    std::vector<std::array<Scalar, 3>> well_phase_rates;
    std::vector<Scalar> well_bhp;
    /// A solve that did not converge: its best iterate, the one with the
    /// smallest max norm of the scaled residual over the whole solve, the
    /// start included -- never worse than the start by that measure.
    Scalar best_residual = 0.0;
    std::vector<Scalar> best_node_pressure;
    std::vector<Scalar> best_well_rate;
    std::vector<std::array<Scalar, 3>> best_well_phase_rates;
    std::vector<Scalar> best_well_bhp;
    /// TEMPORARY diagnostics of a solve that did not converge: what the
    /// system reports on its last kDiagnosticIterations iterations and its
    /// final state (SystemBase::describeIteration(), describeState()).
    std::vector<std::string> diagnostics;
};

/// TEMPORARY: how many of the last iterations of a failing solve are described.
inline constexpr int kDiagnosticIterations = 10;

/// Dense square system, solved by Dune. The networks this solves have tens of
/// unknowns, so a dense direct solve is the whole story. Wrapped only to keep
/// (i,j) indexing and to answer "singular" with false instead of an exception,
/// so that the caller can report a failed solve rather than throw out of the
/// well model.
template<class Scalar>
class DenseMatrix
{
public:
    explicit DenseMatrix(const int n) : a_(n, n, Scalar{0}) {}

    Scalar& operator()(const int i, const int j) { return a_[i][j]; }
    Scalar operator()(const int i, const int j) const { return a_[i][j]; }

    /// Solves A y = b. False if A is singular to working precision.
    bool solve(const std::vector<Scalar>& b, std::vector<Scalar>& y) const
    {
        const auto n = a_.N();
        Dune::DynamicVector<Scalar> rhs(n), x(n, Scalar{0});
        std::copy(b.begin(), b.end(), rhs.begin());
        try {
            a_.solve(x, rhs);
        } catch (const Dune::FMatrixError&) {
            return false;
        }
        y.assign(x.begin(), x.end());
        return true;
    }

private:
    Dune::DynamicMatrix<Scalar> a_;
};

/// What solve() needs of a network system: its residual and Jacobian in the
/// unknowns it chooses (node pressures and whatever else), the controls it
/// picks from an iterate, and its results.
template<class Scalar>
class SystemBase
{
public:
    using State = std::vector<Scalar>;
    using ScalarType = Scalar;

    virtual ~SystemBase() = default;

    virtual int size() const = 0;
    virtual int numWells() const = 0;
    /// Magnitude of unknown i, for sizing a difference and scaling a residual.
    virtual Scalar columnScale(const int i) const = 0;

    /// The whole unknown vector to start from, given a guess at the node pressures.
    virtual State start(const State& node_pressure) const = 0;
    virtual State residual(const State& x) const = 0;
    virtual DenseMatrix<Scalar> jacobian(const State& x) const = 0;
    /// False when jacobian() is not assembled from the table derivatives, and the
    /// solver should difference the residual instead.
    virtual bool usesAnalyticJacobian() const = 0;

    /// Pick each well's control from the iterate; true if any of them moved.
    virtual bool updateControls(const State& x) = 0;
    /// One letter per well, so a cycling active set can be read off the trace.
    virtual char controlLetter(const int w) const = 0;

    /// Shorten a Newton step so it stays where the tables are defined.
    virtual State limitStep(const State& x, const State& dx) const = 0;

    virtual State pressures(const State& x) const = 0;
    virtual State wellRates(const State& x) const = 0;
    /// Each well's phase rates and bhp.
    virtual std::vector<std::array<Scalar, 3>> wellPhaseRates(const State& x) const = 0;
    virtual State wellBhps(const State& x) const = 0;

    /// TEMPORARY diagnostics, for a solve close to its iteration limit: lines
    /// describing iteration \p it -- its worst residual rows and its Newton
    /// step, before (\p dx_raw) and after (\p dx_limited) limitStep() -- and
    /// the final state. Empty by default.
    virtual std::vector<std::string> describeIteration(const int /*it*/, const State& /*x*/,
                                                       const State& /*r*/,
                                                       const State& /*dx_raw*/,
                                                       const State& /*dx_limited*/,
                                                       const bool /*controls_moved*/) const
    {
        return {};
    }
    virtual std::vector<std::string> describeState(const State& /*x*/) const { return {}; }
};

/// Take the Newton step as it comes.
struct FullStep
{
    template<class Sys, class State>
    State accept(const Sys&, const State& x, const State&, const State& dx) const
    {
        State next(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            next[i] = x[i] + dx[i];
        }
        return next;
    }
};

/// Backtrack until the residual norm drops. Useful where the residual is not
/// smooth (table knots, flattened stretches of a tubing curve).
struct LineSearch
{
    int max_halvings = 12;

    template<class Sys, class State>
    State accept(const Sys& system, const State& x, const State& r, const State& dx) const
    {
        auto norm2 = [](const State& v) {
            double s = 0.0;
            for (const auto e : v) {
                s += e * e;
            }
            return std::sqrt(s);
        };
        const double f0 = norm2(r);
        double lambda = 1.0;
        State trial(x.size());
        for (int k = 0; k < max_halvings; ++k) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                trial[i] = x[i] + lambda * dx[i];
            }
            if (norm2(system.residual(trial)) < f0) {
                return trial;
            }
            lambda *= 0.5;
        }
        return trial;
    }
};

/// Convergence settings for solve(). No defaults: a caller states what it wants.
template<class Scalar>
struct Parameters
{
    /// Max norm of the scaled residual at which the system is converged.
    Scalar tolerance;
    int max_iterations;
    /// From this iteration on, steps go through LineSearch instead of the
    /// caller's globalisation (0: never). Full steps converge most solves; a
    /// line search late in the solve breaks the cycles they can fall into
    /// without stalling the solves they converge.
    int line_search_from = 0;
};

/// Solve a network system by Newton-Raphson, letting the system choose its
/// wells' controls from each iterate (SystemBase::updateControls()). Takes the
/// node pressures to start from; returns the converged pressures and rates, or
/// a Result with converged false, the reason it stopped and its best iterate.
template<class Sys, class Globalisation>
Result<typename Sys::ScalarType>
solve(Sys& system,
      const std::vector<typename Sys::ScalarType>& node_pressure_guess,
      const Parameters<typename Sys::ScalarType> params,
      Globalisation globalisation)
{
    const auto tolerance = params.tolerance;
    const int max_iterations = params.max_iterations;
    using Scalar = typename Sys::ScalarType;
    auto x = system.start(node_pressure_guess);
    const int n = system.size();
    Result<Scalar> last;
    std::vector<std::string> trace;
    auto joined = [&trace] {
        std::string out;
        for (const auto& e : trace) {
            out += (out.empty() ? "" : " ") + e;
        }
        return out;
    };

    int switches = 0;
    auto best_x = x;
    Scalar best_worst = std::numeric_limits<Scalar>::infinity();
    // TEMPORARY: described iterations of a solve that may be about to fail.
    std::vector<std::string> diagnostics;
    const int describe_from = max_iterations - kDiagnosticIterations + 1;
    auto describe = [&](const int it, const auto& x_it, const auto& r_it, const auto& dx_raw,
                        const auto& dx_limited, const bool moved) {
        if (it < describe_from && it > 2) {   // DEBUG: also the first two, where it started
            return;
        }
        auto lines = system.describeIteration(it, x_it, r_it, dx_raw, dx_limited, moved);
        diagnostics.insert(diagnostics.end(), lines.begin(), lines.end());
    };
    // An unconverged result: where it stopped, and its best iterate.
    auto unconverged = [&](const auto& x_end) {
        last.node_pressure = system.pressures(x_end);
        last.well_rate = system.wellRates(x_end);
        last.best_residual = best_worst;
        last.best_node_pressure = system.pressures(best_x);
        last.best_well_rate = system.wellRates(best_x);
        last.best_well_phase_rates = system.wellPhaseRates(best_x);
        last.best_well_bhp = system.wellBhps(best_x);
        last.diagnostics = std::move(diagnostics);
        const auto state = system.describeState(x_end);
        last.diagnostics.insert(last.diagnostics.end(), state.begin(), state.end());
        return last;
    };
    for (int it = 1; it <= max_iterations; ++it) {
        const bool controls_moved = system.updateControls(x);
        switches += controls_moved ? 1 : 0;
        const auto r = system.residual(x);

        // std::max(worst, std::abs(e)) would keep worst's old value for a NaN
        // e (every comparison with NaN is false), letting a NaN residual entry
        // pass the convergence test. Any non-finite entry fails it instead.
        Scalar worst = 0.0;
        for (const auto e : r) {
            if (!std::isfinite(e)) {
                worst = std::numeric_limits<Scalar>::infinity();
                break;
            }
            worst = std::max(worst, std::abs(e));
        }
        if (worst < best_worst) {
            best_worst = worst;
            best_x = x;
        }
        {   // remember the active set, so a cycle can be seen in the report
            std::string set;
            for (int w = 0; w < system.numWells(); ++w) {
                set += system.controlLetter(w);
            }
            trace.push_back(set);
            if (trace.size() > 8) {
                trace.erase(trace.begin());
            }
        }
        if (worst < tolerance && !controls_moved) {
            Result<Scalar> done{true, it, system.pressures(x), system.wellRates(x), worst,
                                false, {}, switches};
            done.well_phase_rates = system.wellPhaseRates(x);
            done.well_bhp = system.wellBhps(x);
            return done;
        }
        last = {false, it, {}, {}, worst, controls_moved, joined(), switches};

        DenseMatrix<Scalar> J = system.usesAnalyticJacobian()
            ? system.jacobian(x)
            : [&] {
                  DenseMatrix<Scalar> fd(n);
                  for (int j = 0; j < n; ++j) {
                      // The difference step goes through limitStep() like the
                      // accepted step: a column at a hard bound (a THP well's
                      // bhp at its shut-in bhp) would otherwise be perturbed
                      // past it, where the residual is not defined. Forward
                      // first, backward if the forward step is clamped away.
                      const Scalar h_nominal = 1e-2 * system.columnScale(j);
                      std::vector<Scalar> unit_dx(n, Scalar{0});
                      unit_dx[j] = h_nominal;
                      auto limited = system.limitStep(x, unit_dx);
                      if (std::abs(limited[j]) < Scalar{1e-6} * std::abs(h_nominal)) {
                          unit_dx[j] = -h_nominal;
                          limited = system.limitStep(x, unit_dx);
                      }
                      const Scalar h = limited[j];
                      auto shifted = x;
                      for (int k = 0; k < n; ++k) {
                          shifted[k] += limited[k];
                      }
                      const auto rj = system.residual(shifted);
                      for (int i = 0; i < n; ++i) {
                          fd(i, j) = (h != Scalar{0}) ? (rj[i] - r[i]) / h : Scalar{0};
                      }
                  }
                  return fd;
              }();

        std::vector<Scalar> negative(n), dx;
        for (int i = 0; i < n; ++i) {
            negative[i] = -r[i];
        }
        if (!J.solve(negative, dx)) {
            describe(describe_from, x, r, std::vector<Scalar>{}, std::vector<Scalar>{}, controls_moved);
            diagnostics.push_back(fmt::format("iteration {}: singular Jacobian", it));
            return unconverged(x);
        }

        const auto dx_raw = dx;
        dx = system.limitStep(x, dx);
        describe(it, x, r, dx_raw, dx, controls_moved);
        // The residual jumps when a control switches, and that jump is not a
        // failure to make progress. Letting a globalisation veto it stalls the
        // active set instead of resolving it.
        if (controls_moved) {
            for (int i = 0; i < n; ++i) {
                x[i] += dx[i];
            }
        } else if (params.line_search_from > 0 && it >= params.line_search_from) {
            x = LineSearch{}.accept(system, x, r, dx);
        } else {
            x = globalisation.accept(system, x, r, dx);
        }
    }
    last.iterations = max_iterations + 1;
    return unconverged(x);
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_SOLVE_HEADER_INCLUDED
