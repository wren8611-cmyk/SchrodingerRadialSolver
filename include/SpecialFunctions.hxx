// special_functions.hpp
#pragma once

#include <cmath>
#include <stdexcept>

#include <boost/math/special_functions/bessel.hpp>

namespace specialfunc {

inline double riccati_j(int l, double x) {
    if (l < 0)              throw std::invalid_argument("special::riccati_j: l must be >= 0.");
    if (!std::isfinite(x))  throw std::invalid_argument("special::riccati_j: x must be finite.");

    return x * boost::math::sph_bessel(l, x);
}

inline double riccati_n(int l, double x) {
    if (l < 0)                  throw std::invalid_argument("special::riccati_n: l must be >= 0.");
    if (!std::isfinite(x))      throw std::invalid_argument("special::riccati_n: x must be finite.");
    if (std::abs(x) < 1e-14)    throw std::runtime_error("special::riccati_n: singular at x = 0.");
    
    return x * boost::math::sph_neumann(l, x);
}

inline double riccati_j_derv(int l, double x) {
    if (l < 0)                  throw std::invalid_argument("special::riccati_j_derv: l must be >= 0.");
    if (!std::isfinite(x))      throw std::invalid_argument("special::riccati_j_derv: x must be finite.");
    if (l == 0)                 return std::cos(x);
    if (std::abs(x) < 1e-14)    throw std::runtime_error("special::riccati_j_derv: unstable near x = 0.");

    return riccati_j(l - 1, x) - (static_cast<double>(l) / x) * riccati_j(l, x);
}

inline double riccati_n_derv(int l, double x) {
    if (l < 0)                  throw std::invalid_argument("special::riccati_n_derv: l must be >= 0.");
    if (!std::isfinite(x))      throw std::invalid_argument("special::riccati_n_derv: x must be finite.");
    if (l == 0)                 return std::sin(x);
    if (std::abs(x) < 1e-14)    throw std::runtime_error("special::riccati_n_derv: singular near x = 0.");

    return riccati_n(l - 1, x) - (static_cast<double>(l) / x) * riccati_n(l, x);
}

} // namespace specialfunc