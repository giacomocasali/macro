#pragma once
// SignalProcessing.h
// DSP utilities used by the event loop:
//   correctBaseline()        — subtract median of pre-signal samples
//   laserTriggerTime()       — first threshold crossing on laser channel
//   triggerWindowIndices()   — convert [t_start,t_end] ns to sample indices
//   findFirstPeak()          — internal: first physical peak in -dN/dV (used by Calibration.h)

#include "Config.h"
#include "ButterworthFilter.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

// Subtract the median of samples in [t_base_start, t_base_end) ns.
// Returns the signal unchanged (with a warning) if the window is empty.
//
// BASELINE QUALITY CHECK (new):
//   If max_rms > 0, computes the RMS of baseline samples around the median.
//   If RMS > max_rms → the baseline is contaminated (dark count, afterpulse,
//   or signal leaking into the baseline window).
//   In that case *baseline_ok is set to false and the caller should skip
//   this waveform entirely.
//
//   Typical values:
//     - Without LP filter: max_rms = 3.0–5.0 mV (RF noise is ~1-2 mV RMS)
//     - With LP filter:    max_rms = 1.5–3.0 mV (filtered noise is lower)
//     - Negative or zero:  no check (backward compatible)
//
//   The RMS is computed AFTER median subtraction (= around zero), so it
//   measures the spread of the baseline, not its absolute level.
//
static std::vector<double> correctBaseline(const std::vector<double>& time,
                                           const std::vector<double>& amp,
                                           double t_base_start,
                                           double t_base_end,
                                           double max_rms    = -1.0,
                                           bool*  baseline_ok = nullptr) {
    if (baseline_ok) *baseline_ok = true;  // optimistic default

    std::vector<double> pre;
    for (size_t i = 0; i < time.size(); ++i)
        if (time[i] >= t_base_start && time[i] < t_base_end)
            pre.push_back(amp[i]);
    if (pre.empty()) {
        std::cerr << "  [WARN] No samples in baseline window ["
                  << t_base_start << ", " << t_base_end
                  << ") ns — baseline NOT corrected.\n";
        if (baseline_ok) *baseline_ok = false;
        return amp;
    }

    // Median
    std::vector<double> tmp = pre;
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
    double offset = tmp[tmp.size()/2];

    // RMS check: compute standard deviation of baseline samples around median
    if (max_rms > 0.0 && baseline_ok) {
        double sum2 = 0.0;
        for (double v : pre) {
            double d = v - offset;
            sum2 += d * d;
        }
        double rms = std::sqrt(sum2 / pre.size());
        if (rms > max_rms) {
            *baseline_ok = false;
            // Non correggiamo nemmeno — la waveform va scartata
            return amp;
        }
    }

    std::vector<double> out = amp;
    for (auto& a : out) a -= offset;
    return out;
}

// Find the first physical peak in -dN/dV (= 1 p.e. candidate).
// Run-length rise+fall: requires MIN_RISE consecutive rising steps then MIN_FALL
// falling steps on the 9-point smoothed derivative, above REFMAX_THR=8 mV.
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

    // Hard fallback: absolute max above REFMAX_THR (silent)
    if (iPeak < 0) {
        for (int i = 0; i < n; ++i)
            if (thr[i] > REFMAX_THR && (iPeak < 0 || sm[i] > sm[iPeak]))
                iPeak = i;
    }

    return iPeak;
}
// First threshold crossing on the laser channel with linear interpolation.
// Baseline corrected over [0, BASELINE_END). Returns -999 if not found.
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

// Convert [t_start, t_end] ns to sample indices j_start, j_end.
// Half-sample epsilon avoids floating-point boundary misses.
// Falls back to full waveform with a warning if window is invalid.
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
