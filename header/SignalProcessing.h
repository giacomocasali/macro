#pragma once
// header/SignalProcessing.h
// DSP utilities and physics analysis functions.
// No ROOT graphics — these functions are reusable in other macros.

#include "Config.h"
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
//  4th-ORDER BUTTERWORTH LOW-PASS FILTER  —  zero-phase (filtfilt)
//
//  Equivalent to scipy.signal.butter(4, fc/(fs/2)) + filtfilt.
//  Implementation:
//    1. Compute SOS (second-order sections) via bilinear transform
//       with frequency pre-warping — identical to scipy.
//    2. Apply each biquad section forward then backward (filtfilt).
//       Forward-backward cancels phase shift → zero group delay.
//    3. Boundary padding (odd-extension, length 3*(order+1)=15 samples)
//       mirrors scipy's default "odd" pad to suppress edge transients.
//
//  fc  = cut-off frequency (MHz)
//  fs  = sampling frequency (MHz)
// ════════════════════════════════════════════════════════════

// Apply one biquad section y[n] = b0*x[n]+b1*x[n-1]+b2*x[n-2]
//                                        -a1*y[n-1]-a2*y[n-2]
// Direct form II transposed (same as scipy lfilter).
static std::vector<double> applyBiquad(const std::vector<double>& x,
                                       double b0, double b1, double b2,
                                       double a1, double a2) {
    size_t n = x.size();
    std::vector<double> y(n, 0.0);
    double s1 = 0.0, s2 = 0.0;   // state variables
    for (size_t i = 0; i < n; ++i) {
        double xi = x[i];
        y[i] = b0*xi + s1;
        s1   = b1*xi - a1*y[i] + s2;
        s2   = b2*xi - a2*y[i];
    }
    return y;
}

static std::vector<double> butterworthLowPass(const std::vector<double>& x,
                                              double fc, double fs) {
    // ── Compute SOS coefficients (bilinear transform, pre-warped) ──
    // Matches scipy butter(4, fc/(fs/2), output='sos') exactly.
    const int N_SECTIONS = 2;   // order/2
    double wc = 2.0 * fs * std::tan(M_PI * fc / fs);
    double c  = 2.0 * fs;

    double b0[N_SECTIONS], b1[N_SECTIONS], b2[N_SECTIONS];
    double a1[N_SECTIONS], a2[N_SECTIONS];

    for (int k = 0; k < N_SECTIONS; ++k) {
        // Poles in reverse order (k=0 → section 1 of scipy, k=1 → section 0)
        double phi = M_PI/2.0 + M_PI*(2.0*(N_SECTIONS-1-k)+1.0)/(2.0*4.0);
        double pRe = wc * std::cos(phi);
        double pIm = wc * std::sin(phi);
        double pm2 = pRe*pRe + pIm*pIm;
        double d0  =  c*c - 2.0*pRe*c + pm2;
        double d1  =  2.0*pm2 - 2.0*c*c;
        double d2  =  c*c + 2.0*pRe*c + pm2;
        double g   =  pm2 / d0;
        b0[k] = g;  b1[k] = 2.0*g;  b2[k] = g;
        a1[k] = d1/d0;
        a2[k] = d2/d0;
    }

    // ── Odd-extension padding (mirrors scipy filtfilt default) ──────
    // Pad length = 3 * (order+1) = 15 samples on each side.
    // Odd extension: pad[i] = 2*x[0] - x[padlen-i]  (left)
    //                pad[i] = 2*x[n-1] - x[n-1-(i-n+padlen)] (right)
    int padlen = 15;
    int n  = (int)x.size();
    int np = n + 2*padlen;
    std::vector<double> xp(np);
    for (int i = 0; i < padlen; ++i)
        xp[i] = 2.0*x[0] - x[padlen - i];          // left odd extension
    for (int i = 0; i < n; ++i)
        xp[padlen + i] = x[i];                       // original signal
    for (int i = 0; i < padlen; ++i)
        xp[padlen + n + i] = 2.0*x[n-1] - x[n-1-1-i]; // right odd extension

    // ── Forward pass through all sections ────────────────────────────
    std::vector<double> y = xp;
    for (int s = 0; s < N_SECTIONS; ++s)
        y = applyBiquad(y, b0[s], b1[s], b2[s], a1[s], a2[s]);

    // ── Reverse ──────────────────────────────────────────────────────
    std::reverse(y.begin(), y.end());

    // ── Backward pass through all sections ───────────────────────────
    for (int s = 0; s < N_SECTIONS; ++s)
        y = applyBiquad(y, b0[s], b1[s], b2[s], a1[s], a2[s]);

    // ── Reverse back and remove padding ──────────────────────────────
    std::reverse(y.begin(), y.end());
    return std::vector<double>(y.begin() + padlen, y.begin() + padlen + n);
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
