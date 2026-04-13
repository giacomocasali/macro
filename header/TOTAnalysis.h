#pragma once
// header/TOTAnalysis.h
// TOT computation and event selection by number of p.e.
//
// v2 fixes (validated in sipm_waveform_diag v3):
//   - Confirmation timeout on falling edge (100 samples = 20 ns at 5 GS/s)
//     Rejects multi-p.e. tails that slowly decay through thr_lo.
//   - Post-fall check uses threshold (not edge_limit) so it works
//     even when edge_thr_frac=100 (no filter mode).
//   - End-of-trace check: rejects events where the SiPM tail is still
//     above thr_lo at the end of the acquisition window.

#include "Config.h"
#include "SignalProcessing.h"
#include "Calibration.h"
#include <vector>
#include <cmath>
#include <algorithm>

// TOT event: all computed quantities for one waveform
struct TOTEvent {
    double tot;      // Time Over Threshold [ns]
    double delta_t;  // t_LET - t_laser [ns]
    double amp_max;  // peak amplitude of the waveform [mV]
    int    n_pe;     // estimated number of p.e. from calibration
};

// Median of amp[i0 .. i0+n_pts-1] via nth_element.
static double edgeMedian(const std::vector<double>& amp,
                         int i0, int n_pts)
{
    int sz  = (int)amp.size();
    int i1  = std::max(0, i0);
    int i2  = std::min(sz, i0 + n_pts);
    if (i2 <= i1) return 0.0;
    std::vector<double> tmp(amp.begin() + i1, amp.begin() + i2);
    int mid = (int)tmp.size() / 2;
    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
    return tmp[mid];
}

// ════════════════════════════════════════════════════════════
//  TOT computation with hysteresis on both edges.
//  Returns {t_rise, t_fall} at threshold.
//  Returns {-1, -1} in any rejection case.
//
//  Parameters:
//    n_edge              = 100   samples at waveform edges for stability check
//    edge_thr_frac       = 0.5   edge median must be < frac*threshold
//    hyst_frac           = 0.5   rise arming: signal must go below thr*hyst_frac
//    pre_check_n         = 10    samples before t_rise for pre-crossing check
//    pre_check_frac      = 0.3   pre-crossing median < frac*threshold
//    min_fall_samp       = 3     consecutive samples below threshold to confirm fall
//    post_check_n        = 50    samples after j_fall_start for post-fall check.
//                                0 = check all remaining samples to j_end.
//    check_trailing_edge = false  whether to apply the trailing edge-stability check
//    confirm_timeout     = 100   max samples between fall candidate and confirmation.
//                                At 5 GS/s, 100 = 20 ns. A 1 p.e. signal drops from
//                                threshold to thr_lo in ~2-5 ns. A 3+ p.e. tail takes
//                                80+ ns — the timeout rejects it.
//                                Set <= 0 to disable (backward compatible).
//
//  Rejection conditions:
//  1. Leading edge-stability: median of first n_edge samples >= edge_limit.
//  2. Trailing edge-stability (only if check_trailing_edge=true).
//  3. Rising-edge hysteresis: arms below thr*hyst_frac, fires above thr.
//  4. Pre-crossing check: median of pre_check_n samples before t_rise
//     must be < pre_check_frac*threshold.
//  5. Falling-edge: t_fall accepted only after min_fall_samp consecutive
//     samples below threshold AND confirmation that the signal reaches
//     thr_lo = threshold*hyst_frac within confirm_timeout samples.
//  6. Post-fall baseline check: median of post_check_n samples after
//     t_fall must be < threshold (the LET level).
//  7. End-of-trace check: median of last 50 samples must be < thr_lo.
//  8. t_rise or t_fall not found within [j_start, j_end].
// ════════════════════════════════════════════════════════════
static std::pair<double,double> computeTOT(const std::vector<double>& time,
                                           const std::vector<double>& amp,
                                           double threshold,
                                           int j_start, int j_end,
                                           int    n_edge              = 100,
                                           double edge_thr_frac       = 0.5,
                                           double hyst_frac           = 0.5,
                                           int    pre_check_n         = 10,
                                           double pre_check_frac      = 0.3,
                                           int    min_fall_samp       = 3,
                                           int    post_check_n        = 50,
                                           bool   check_trailing_edge = false,
                                           int    confirm_timeout     = 100)
{
    // ── 1. Edge-stability checks ─────────────────────────────
    const double edge_limit = edge_thr_frac * threshold;
    if (edgeMedian(amp, 0, n_edge) >= edge_limit)
        return {-1.0, -1.0};
    if (check_trailing_edge &&
        edgeMedian(amp, (int)amp.size() - n_edge, n_edge) >= edge_limit)
        return {-1.0, -1.0};

    // Hysteresis levels (same for rising and falling)
    const double thr_lo = threshold * hyst_frac;

    // ── 2+3. Rising-edge search with hysteresis + pre-check ─
    double t_rise = -1.0, t_fall = -1.0;
    bool armed_rise = false;

    for (int j = j_start; j <= j_end; ++j) {
        if (!armed_rise) {
            if (amp[j] < thr_lo) armed_rise = true;
        } else {
            if (amp[j] >= threshold) {
                double tr = time[j-1] + (threshold - amp[j-1])
                            * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
                int pre_start = std::max(j_start, j - pre_check_n);
                int pre_n     = j - pre_start;
                if (pre_n > 0 &&
                    edgeMedian(amp, pre_start, pre_n) >= pre_check_frac * threshold) {
                    armed_rise = false;
                    continue;
                }
                t_rise = tr;
                break;
            }
        }
    }
    if (t_rise < 0) return {-1.0, -1.0};

    // ── 4. Falling-edge search with hysteresis + timeout ─────
    // Logic:
    //   1. Find min_fall_samp consecutive samples below threshold
    //      → candidate t_fall (interpolated at first crossing)
    //   2. CONFIRM: signal must reach thr_lo within confirm_timeout
    //      samples. If it rises above threshold first → false falling.
    //      If timeout expires → multi-p.e. tail, reject candidate.
    //
    // Without LP filter, RF noise oscillates the SiPM tail
    // above/below threshold. The hysteresis confirmation prevents
    // false falling edges. The timeout prevents accepting multi-p.e.
    // tails that slowly decay through thr_lo after 80+ ns.
    {
        int  consec_below = 0;
        int  j_fall_candidate = -1;

        enum { SEARCHING, CONFIRMING } state = SEARCHING;
        double t_fall_candidate = -1.0;
        int    j_confirm_start = -1;

        for (int j = j_start; j <= j_end; ++j) {
            if (time[j] <= t_rise) continue;

            if (state == SEARCHING) {
                if (amp[j] < threshold) {
                    if (consec_below == 0) j_fall_candidate = j;
                    ++consec_below;
                    if (consec_below >= min_fall_samp) {
                        int jf = j_fall_candidate;
                        if (jf > 0)
                            t_fall_candidate = time[jf-1] + (threshold - amp[jf-1])
                                * (time[jf] - time[jf-1]) / (amp[jf] - amp[jf-1]);
                        else
                            t_fall_candidate = time[jf];
                        j_confirm_start = j;
                        state = CONFIRMING;
                    }
                } else {
                    consec_below = 0;
                    j_fall_candidate = -1;
                }
            } else {
                // CONFIRMING: signal must reach thr_lo
                if (amp[j] <= thr_lo) {
                    // Confirmed: genuine return to baseline
                    t_fall = t_fall_candidate;
                    break;
                }
                if (amp[j] >= threshold) {
                    // Rose above threshold without confirmation → false falling
                    state = SEARCHING;
                    consec_below = 0;
                    j_fall_candidate = -1;
                    t_fall_candidate = -1.0;
                }
                // Timeout: confirmation taking too long → multi-p.e. tail
                if (confirm_timeout > 0 && (j - j_confirm_start) > confirm_timeout) {
                    state = SEARCHING;
                    consec_below = 0;
                    j_fall_candidate = -1;
                    t_fall_candidate = -1.0;
                }
            }
        }
    }
    if (t_fall < 0) return {-1.0, -1.0};

    // ── 5. Post-fall baseline check ──────────────────────────
    // Uses threshold (the LET level) as the post-fall limit.
    // After the falling edge, the median of the next post_check_n
    // samples must be BELOW threshold. This catches false falling
    // edges and pile-up events.
    //
    // Note: we use threshold here (not edge_limit) because
    // edge_limit = edge_thr_frac * threshold, and when
    // edge_thr_frac=100 (no filter), edge_limit=1130 mV which
    // would disable the check entirely.
    {
        int jf_post = 0;
        for (int j = j_start; j <= j_end; ++j)
            if (time[j] >= t_fall) { jf_post = j; break; }

        int post_limit = (post_check_n > 0)
                         ? std::min(jf_post + post_check_n, j_end + 1)
                         : j_end + 1;
        int post_n = post_limit - jf_post;
        if (post_n > 0 && edgeMedian(amp, jf_post, post_n) >= threshold)
            return {-1.0, -1.0};
    }

    // ── 6. End-of-trace check ────────────────────────────────
    // If the signal is still above thr_lo at the end of the
    // acquisition window, the SiPM has not recovered and the
    // TOT measurement is unreliable.
    {
        const int tail_n = 50;  // last 50 samples (~10 ns at 5 GS/s)
        int tail_start = std::max(j_start, j_end + 1 - tail_n);
        int tail_len   = j_end + 1 - tail_start;
        if (tail_len > 0 && edgeMedian(amp, tail_start, tail_len) >= thr_lo)
            return {-1.0, -1.0};
    }

    return {t_rise, t_fall};
}

// n_pe = round((amp_max - cal.q) / cal.m). Returns -1 if unclassifiable.
static int estimatePE(double amp_max, const CalibResult& cal,
                      double tolerance = 0.4) {
    if (!cal.ok || cal.m <= 0) return -1;
    double n_raw = (amp_max - cal.q) / cal.m;
    int    n_pe  = (int)std::round(n_raw);
    if (n_pe < 0) return 0;
    // Accept only if within tolerance*m from the nearest integer peak
    if (std::abs(n_raw - n_pe) > tolerance) return -1;
    return n_pe;
}

// Pre-laser quiet check: scans waveform from BASELINE_END to t_laser.
// Returns false if any sample exceeds quiet_thr (dark count before laser).
static bool preLaserQuiet(const std::vector<double>& time,
                          const std::vector<double>& amp,
                          double t_laser,
                          double quiet_thr,
                          double margin_ns = 1.0)
{
    double t_check_end = t_laser - margin_ns;
    if (t_check_end <= BASELINE_END) return true;

    for (size_t j = 0; j < time.size(); ++j) {
        if (time[j] < BASELINE_END) continue;
        if (time[j] > t_check_end) break;
        if (amp[j] > quiet_thr) return false;
    }
    return true;
}

// Legacy event loop (single file). Prefer collectTOTEvents_fileByFile in VbiasAnalysis_v2.h.
static std::vector<TOTEvent> collectTOTEvents(
        TTree* treeCh1,
        TTree* treeLaser,
        double cutoff_MHz,
        double fs_MHz,
        int j_trig_start, int j_trig_end,
        double let_thr,
        const CalibResult& cal,
        const std::string& tag) {

    std::vector<TOTEvent> events;
    const int N = 1024;
    Double_t t1[N], a1[N], tL[N], aL[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);
    treeLaser->SetBranchAddress("time",      tL);
    treeLaser->SetBranchAddress("amplitude", aL);

    Long64_t nEntries = treeCh1->GetEntries();
    long nNoLaser=0, nNoTOT=0;

    std::cout << "  Processing " << nEntries << " events...\n";

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        treeLaser->GetEntry(i);

        double t_laser = laserTriggerTime(tL, aL, N, cal.laser_thr);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        std::vector<double> v_t(t1,t1+N), v_a(a1,a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, BASELINE_START, BASELINE_END),
            cutoff_MHz, fs_MHz);

        double amp_max = *std::max_element(af.begin() + j_trig_start,
                                            af.begin() + j_trig_end + 1);

        auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                                            j_trig_start, j_trig_end);
        if (t_rise < 0 || t_fall < 0) { ++nNoTOT; continue; }

        double tot     = t_fall - t_rise;
        double delta_t = t_rise - t_laser;
        if (tot <= 0 || tot >= 100.0) continue;

        int n_pe = estimatePE(amp_max, cal);
        events.push_back({tot, delta_t, amp_max, n_pe});
    }

    std::cout << "  Collected: " << events.size()
              << "  |  No laser: " << nNoLaser
              << "  |  No TOT: "   << nNoTOT << "\n";
    return events;
}
