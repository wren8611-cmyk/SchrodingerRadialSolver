#pragma once

#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "WaveFuncResult.hxx"

class NumerovWaveFuncSolver final {
public:
    using PotentialFunction = std::function<double(double)>;

    NumerovWaveFuncSolver(PotentialFunction potential,
                          int l,
                          double mu,
                          double r_max,
                          int npoints)
        : potential_(std::move(potential)), l_(l), mu_(mu), r_max_(r_max), npoints_(npoints)
    {
        if (!potential_)    throw std::invalid_argument("NumerovWaveFuncSolver: potential must be valid.");
        if (l_ < 0)         throw std::invalid_argument("NumerovWaveFuncSolver: l must be >= 0.");
        if (mu_ <= 0.0)     throw std::invalid_argument("NumerovWaveFuncSolver: mu must be > 0.");
        if (r_max_ <= 0.0)  throw std::invalid_argument("NumerovWaveFuncSolver: r_max must be > 0.");
        if (npoints_ < 3)   throw std::invalid_argument("NumerovWaveFuncSolver: npoints must be >= 3.");

        constexpr double r_min = 1e-14; // internal offset to avoid r = 0 singularity
        h_ = (r_max_ - r_min) / static_cast<double>(npoints_ - 1);

        if (!(h_ > 0.0) || !std::isfinite(h_))  throw std::invalid_argument("NumerovWaveFuncSolver: invalid grid spacing.");

        r_.resize(npoints_);
        for (int i = 0; i < npoints_; ++i)      r_[i] = r_min + static_cast<double>(i) * h_;
    }

    // ===== Numerov Solver (wrapped) ===== //
    [[nodiscard]] WaveFuncResult solve(double k) const {
        if (k <= 0.0)   throw std::invalid_argument("NumerovWaveFuncSolver: k must be > 0.");
        
        return WaveFuncResult{k, energy(k), l_, mu_, r_, integrate(k)};
    }

    // ===== accessors ===== //
    [[nodiscard]] const std::vector<double>& grid() const noexcept  { return r_; }
    [[nodiscard]] double step_size() const noexcept                 { return h_; }

private:
    // ===== Numerov kernel ===== //
    [[nodiscard]] double effective_q(double r, double k) const {
        const double V = potential_(r);
        if (!std::isfinite(V)) throw std::runtime_error("NumerovWaveFuncSolver: potential is not finite at r = " + std::to_string(r) + ".");

        const double centrifugal = static_cast<double>(l_) * static_cast<double>(l_ + 1) / (r * r);

        return 2.0 * mu_ * V - k * k + centrifugal;
    }

    // ===== Numerov Solver (Core) ===== //
    [[nodiscard]] std::vector<double> integrate(double k) const {
        std::vector<double> u(npoints_, 0.0);
        std::vector<double> q(npoints_, 0.0);

        for (int i = 0; i < npoints_; ++i) {
            q[i] = effective_q(r_[i], k);
            if (!std::isfinite(q[i]))   throw std::runtime_error("NumerovWaveFuncSolver: q(r) not finite at grid index " + std::to_string(i) + ", r = " + std::to_string(r_[i]) + ".");
        }

        // Regular near-origin behavior: u(r) = r^(l+1)
        u[0] = std::pow(r_[0], static_cast<double>(l_ + 1));
        u[1] = std::pow(r_[1], static_cast<double>(l_ + 1));

        const double h2 = h_ * h_;
        const double c = h2 / 12.0;

        constexpr double denom_threshold = 1e-14;
        constexpr double overflow_threshold = 1e100;

        // Numerov algorithm
        for (int i = 1; i < npoints_ - 1; ++i) {
            const double denom = 1.0 - c * q[i + 1];
            if (!std::isfinite(denom) || std::abs(denom) < denom_threshold) {
                throw std::runtime_error("NumerovWaveFuncSolver: unstable Numerov denominator at index " + std::to_string(i) + ", r = " + std::to_string(r_[i]) + ".");
            }

            const double next = (2.0 * (1.0 + 5.0 * c * q[i]) * u[i] - (1.0 - c * q[i - 1]) * u[i - 1]) / denom;

            if (!std::isfinite(next) || std::abs(next) > overflow_threshold) {
                throw std::runtime_error("NumerovWaveFuncSolver: overflow at index " + std::to_string(i + 1) + ", r = " + std::to_string(r_[i + 1]) + ".");
            }

            u[i + 1] = next;
        }

        return u;
    }

    // ===== Energy (non-rela) ===== //
    [[nodiscard]] double energy(double k) const noexcept { return (k * k) / (2.0 * mu_); }

private:
    PotentialFunction potential_;
    int l_ = 0;
    double mu_ = 1.0;
    double r_max_ = 0.0;
    int npoints_ = 0;
    double h_ = 0.0;
    std::vector<double> r_;
};

using NumerovSolver  = NumerovWaveFuncSolver;   // alias