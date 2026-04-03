#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <vector>

#include "include/ShortRangeIntWaveFuncResult.hxx"
#include "include/NumerovWaveFuncSolver.hxx"

// g++ -O3 test_solver.cxx

int main() {
    
    // ===== set potential ===== //
    // V(r) = -V0 exp(-(r/a)^2)
    constexpr double V0 = 5.0;
    constexpr double a  = 1.0;

    auto gaussian = [V0, a](double r) -> double {
        const double x = r / a;
        return -V0 * std::exp(-x * x);
    };

    // ===== phys vals ===== //
    constexpr int       l  = 0;     // partial wave
    constexpr double    k  = 1.0;   // momentum
    constexpr double    mu = 1.0;   // reduced mass

    // ===== grid conf ===== //
    constexpr double rmax = 20.0;
    constexpr int npoints = 20000;

    // ===== solver =====//
    NumerovSolver solver(gaussian, l, mu, rmax, npoints);
    const auto wf = solver.solve(k); // solve -> WaveFuncResult

    std::cout << std::setprecision(16);
    std::cout << "=== Numerov scattering test (gaussian potential) ===\n";
    std::cout << "l       = " << wf.l << "\n";
    std::cout << "k       = " << wf.k << "\n";
    std::cout << "mu      = " << wf.mu << "\n";
    std::cout << "E       = " << wf.E << "\n";
    std::cout << "rmax    = " << rmax << "\n";
    std::cout << "N       = " << npoints << "\n";
    std::cout << "dr      = " << solver.step_size() << "\n";
    std::cout << "V0      = " << V0 << "\n";
    std::cout << "a       = " << a << "\n\n";

    // ===== phase shift ===== //
    SRIntWFResult swf(wf);    // WaveFuncResult -> SRIntWaveFuncResult
    const double delta_boundary = swf.phase_shift();
    const double delta_         = swf.phase_shift(npoints-200);

    std::cout << "Phase shift:\n";
    std::cout << "delta = " << delta_boundary << " rad\n";
    std::cout << "delta = " << delta_         << " rad\n";

    return 0;

}