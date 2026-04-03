#pragma once

#include <utility>

#include "WaveFuncResult.hxx"
#include "PhaseShiftSolver.hxx"

struct ShortRangeIntWaveFuncResult final : public WaveFuncResult {
    ShortRangeIntWaveFuncResult() = default;
    explicit ShortRangeIntWaveFuncResult(const WaveFuncResult& wf)      : WaveFuncResult(wf) {}
    explicit ShortRangeIntWaveFuncResult(WaveFuncResult&& wf) noexcept  : WaveFuncResult(std::move(wf)) {}

    [[nodiscard]] double phase_shift(int match_index) const {
        return SRIntPShift::solve(r, u, l, k, match_index);
    }

    [[nodiscard]] double phase_shift() const {
        return SRIntPShift::solve(r, u, l, k);
    }
};

using SRIntWFResult = ShortRangeIntWaveFuncResult;    // short alias