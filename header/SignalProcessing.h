#pragma once
// header/SignalProcessing.h
// DSP utilities and physics analysis functions.
// No ROOT graphics — these functions are reusable in other macros.

#include "Config.h"
#include "ButterworthFilter.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

// ════════════════════════════════════════════════════════════
//  BASELINE CORRECTION  (median estimator)
//  Subtracts the median of samples in [t_base_start, t_base_end) ns.
//  Returns the uncorrected signal with a warning if the window is empty.
// ════════════════════════════════════════════════════════════
static std::vector<double> correctBaseline(const std::vector<double>& time,
                                           const std::vector<double>& amp,
                                           double t_base_start,
                                           double t_base_end) {
    std::vector<double> pre;
    for (size_t i = 0; i < time.size(); ++i)
        if (time[i] >= t_base_start && time[i] < t_base_end)
            pre.push_back(amp[i]);
    if (pre.empty()) {
        std::cerr << "  [WARN] No samples in baseline window ["
                  << t_base_start << ", " << t_base_end
                  << ") ns — baseline NOT corrected.\n";
        return amp;
    }
    std::vector<double> tmp = pre;
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
    double offset = tmp[tmp.size()/2];
    std::vector<double> out = amp;
    for (auto& a : out) a -= offset;
    return out;
}

// ════════════════════════════════════════════════════════════
//  FIND FIRST PHYSICAL PEAK of -dN/dV  (1st p.e.)
//
//  Run-length rise+fall algorithm.
//  A real p.e. peak = at least MIN_RISE consecutive rising steps on the
//  smoothed derivative immediately followed by at least MIN_FALL
//  consecutive falling steps.  This is the universal "bump" signature.
//  Working only above REFMAX_THR=8 mV skips the near-zero arming spike.
// ════════════════════════════════════════════════════════════
static int findFirstPeak(const std::vector<double>& thr,
                         const std::vector<double>& der) {
    int n = (int)thr.size();
    if (n < 5) return -1;

    const double REFMAX_THR = 8.0;
    const int    MIN_RISE   = 2;
    const int    MIN_FALL   = 2;

    // 9-point box smooth
    std::vector<double> sm(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0; int cnt = 0;
        for (int d = -4; d <= 4; ++d) {
            int j = i + d;
            if (j >= 0 && j < n) { sum += der[j]; ++cnt; }
        }
        sm[i] = sum / cnt;
    }

    // risingBefore[i]: consecutive steps sm went UP ending at i
    std::vector<int> risingBefore(n, 0);
    for (int i = 1; i < n; ++i)
        risingBefore[i] = (sm[i] > sm[i-1]) ? risingBefore[i-1] + 1 : 0;

    // fallingAfter[i]: consecutive steps sm goes DOWN starting from i
    std::vector<int> fallingAfter(n, 0);
    for (int i = n - 2; i >= 0; --i)
        fallingAfter[i] = (sm[i] > sm[i+1]) ? fallingAfter[i+1] + 1 : 0;

    // Spike-immune reference max
    double globalMax = 0;
    for (int i = 0; i < n; ++i)
        if (thr[i] > REFMAX_THR && sm[i] > globalMax) globalMax = sm[i];
    if (globalMax <= 0) return -1;

    // Pass 1: MIN_RISE + MIN_FALL
    int iPeak = -1;
    for (int i = 0; i < n; ++i) {
        if (thr[i] <= REFMAX_THR) continue;
        if (sm[i] < globalMax * 0.03) continue;
        if (risingBefore[i] >= MIN_RISE && fallingAfter[i] >= MIN_FALL) {
            iPeak = i; break;
        }
    }

    // Pass 2: relax to 1 + 1
    if (iPeak < 0) {
        for (int i = 0; i < n; ++i) {
            if (thr[i] <= REFMAX_THR) continue;
            if (sm[i] < globalMax * 0.03) continue;
            if (risingBefore[i] >= 1 && fallingAfter[i] >= 1) {
                iPeak = i; break;
            }
        }
    }

    // Hard fallback: absolute max above REFMAX_THR
    if (iPeak < 0) {
        std::cout << "  [WARN] findFirstPeak: run-length failed"
                  << " -- using absolute max above " << REFMAX_THR << " mV.\n";
        for (int i = 0; i < n; ++i)
            if (thr[i] > REFMAX_THR && (iPeak < 0 || sm[i] > sm[iPeak]))
                iPeak = i;
    }

    return iPeak;
}
// ════════════════════════════════════════════════════════════
//  LASER TRIGGER TIME
//  First crossing above laser_thr on the laser channel,
//  with linear interpolation for sub-sample precision.
//  Baseline is corrected over [0, BASELINE_END) ns.
//  Returns -999 if no crossing is found.
// ════════════════════════════════════════════════════════════
static double laserTriggerTime(const Double_t* t_laser,
                               const Double_t* a_laser,
                               int N,
                               double laser_thr = 10.0) {
    std::vector<double> pre;
    for (int j = 0; j < N; ++j)
        if (t_laser[j] < BASELINE_END) pre.push_back(a_laser[j]);
    double offset = 0.0;
    if (!pre.empty()) {
        std::vector<double> tmp = pre;
        std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
        offset = tmp[tmp.size()/2];
    }
    for (int j = 1; j < N; ++j) {
        double v0 = a_laser[j-1] - offset;
        double v1 = a_laser[j]   - offset;
        if (v1 > laser_thr && v0 <= laser_thr)
            return t_laser[j-1]
                   + (laser_thr - v0) * (t_laser[j] - t_laser[j-1]) / (v1 - v0);
    }
    return -999.0;
}

// ════════════════════════════════════════════════════════════
//  TRIGGER WINDOW INDICES
//  Converts [t_start, t_end] ns to sample indices j_start, j_end.
//  A half-sample epsilon avoids floating-point boundary misses.
// ════════════════════════════════════════════════════════════
static void triggerWindowIndices(const Double_t* t1, int N,
                                 double t_start, double t_end,
                                 int& j_start, int& j_end,
                                 const std::string& tag = "") {
    double eps = 0.5 * (t1[1] - t1[0]);
    j_start = 0; j_end = N - 1;
    for (int j = 0; j < N; ++j)
        if (t1[j] >= t_start - eps) { j_start = j; break; }
    for (int j = N - 1; j >= 0; --j)
        if (t1[j] <= t_end + eps)   { j_end   = j; break; }
    if (j_start >= j_end) {
        std::cerr << "[WARN] Trigger window [" << t_start << ", " << t_end
                  << "] ns invalid" << (tag.empty() ? "" : " for " + tag)
                  << " — using full waveform.\n";
        j_start = 0; j_end = N - 1;
    }
}
