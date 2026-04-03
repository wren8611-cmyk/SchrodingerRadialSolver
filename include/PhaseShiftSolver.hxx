#pragma once

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "SpecialFunctions.hxx"
#include "WaveFuncResult.hxx"

class ShortRangeIntPhaseShiftSolver final {
public:
    // ===== solve at a certain radial point ===== //
    [[nodiscard]] static double solve(const std::vector<double>& r,
                                      const std::vector<double>& u,
                                      int l,
                                      double k,
                                      int match_index)
    {
        validate_inputs(r, u, l, k, match_index);

        constexpr double x_tol = 1e-12;
        constexpr double u_tol = 1e-14;

        const double h = r[1] - r[0];
        const double x = k * r[match_index];

        if (std::abs(x) < x_tol)                throw std::runtime_error("ShortRangeIntPhaseShiftSolver: matching point too close to origin.");
        if (std::abs(u[match_index]) < u_tol)   throw std::runtime_error("ShortRangeIntPhaseShiftSolver: wavefunction too small at match point.");

        const double uprime = (u[match_index + 1] - u[match_index - 1]) / (2.0 * h);
        const double L = uprime / u[match_index];

        const double jl  = specialfunc::riccati_j(l, x);
        const double nl  = specialfunc::riccati_n(l, x);
        const double jlp = specialfunc::riccati_j_derv(l, x);
        const double nlp = specialfunc::riccati_n_derv(l, x);

        const double num = L * jl - k * jlp;
        const double den = L * nl - k * nlp;

        return std::atan2(num, den);
    }

    [[nodiscard]] static double solve(const WaveFuncResult& wf, int match_index) {
        return solve(wf.r, wf.u, wf.l, wf.k, match_index);
    }

    // ===== solve at end point ===== //
    [[nodiscard]] static double solve(const std::vector<double>& r,
                                      const std::vector<double>& u,
                                      int l,
                                      double k)
    {
        return solve(r, u, l, k, static_cast<int>(r.size()) - 2);
    }

    [[nodiscard]] static double solve(const WaveFuncResult& wf) {
        return solve(wf.r, wf.u, wf.l, wf.k);
    }

private:
    static void validate_inputs(const std::vector<double>& r,
                                const std::vector<double>& u,
                                int l,
                                double k,
                                int match_index)
    {
        constexpr double grid_tol = 1e-08;

        if (r.size() != u.size())   throw std::invalid_argument("ShortRangeIntPhaseShiftSolver: r and u must have the same size.");
        if (r.size() < 3)           throw std::invalid_argument("ShortRangeIntPhaseShiftSolver: need at least 3 grid points.");
        if (l < 0)                  throw std::invalid_argument("ShortRangeIntPhaseShiftSolver: l must be >= 0.");
        if (k <= 0.0)               throw std::invalid_argument("ShortRangeIntPhaseShiftSolver: k must be > 0.");

        if (match_index <= 0 || match_index >= static_cast<int>(r.size()) - 1) throw std::invalid_argument("ShortRangeIntPhaseShiftSolver: invalid match_index.");

        const double h = r[1] - r[0];
        if (!(h > 0.0) || !std::isfinite(h))    throw std::invalid_argument("ShortRangeIntPhaseShiftSolver: invalid radial grid.");

        for (std::size_t i = 1; i + 1 < r.size(); ++i) {
            const double hi = r[i + 1] - r[i];
            if (!std::isfinite(hi) || std::abs(hi - h) > grid_tol * std::abs(h))    throw std::invalid_argument("ShortRangeIntPhaseShiftSolver: radial grid must be uniform.");
        }
    }
};

using SRIntPShift = ShortRangeIntPhaseShiftSolver;  // alias