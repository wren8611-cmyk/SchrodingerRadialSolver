#pragma once

#include <vector>

struct WaveFuncResult {
    double k    = 0.0;
    double E    = 0.0;
    int l       = 0;
    double mu   = 0.0;
    std::vector<double> r;
    std::vector<double> u;
};