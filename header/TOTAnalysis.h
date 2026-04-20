#pragma once
// TOTAnalysis.h — v3
// Modifiche rispetto a v2:
//   - TOTEvent ora include delta_t_lo e delta_t_hi
//   - computeTOT_rise_only: calcola solo t_rise a soglia arbitraria,
//     senza filtri di rejection, per eventi gia' validati.

#include "Config.h"
#include "SignalProcessing.h"
#include "Calibration.h"
#include <vector>
#include <cmath>
#include <algorithm>

static constexpr double TRISE_INVALID = -999.0;

// ── TOT event ────────────────────────────────────────────────────────────
struct TOTEvent {
    double tot;         // Time Over Threshold [ns]
    double delta_t;     // t_rise(LET) - t_laser [ns]
    double amp_max;     // peak amplitude [mV]
    int    n_pe;        // estimated p.e.
    // Tempi di arrivo alle soglie derivate (stessa validazione, soglia diversa)
    double delta_t_lo;  // t_rise(0.2*LET) - t_laser [ns]  (TRISE_INVALID se non trovato)
    double delta_t_hi;  // t_rise(0.8*LET) - t_laser [ns]  (TRISE_INVALID se non trovato)
};

// ── edgeMedian ───────────────────────────────────────────────────────────
static double edgeMedian(const std::vector<double>& amp, int i0, int n_pts)
{
    int sz = (int)amp.size();
    int i1 = std::max(0, i0);
    int i2 = std::min(sz, i0 + n_pts);
    if (i2 <= i1) return 0.0;
    std::vector<double> tmp(amp.begin() + i1, amp.begin() + i2);
    int mid = (int)tmp.size() / 2;
    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
    return tmp[mid];
}

// ════════════════════════════════════════════════════════════
//  computeTOT — invariato rispetto a v2.
//  Tutti i filtri di rejection attivi.
// ════════════════════════════════════════════════════════════
static std::pair<double,double> computeTOT(
        const std::vector<double>& time,
        const std::vector<double>& amp,
        double threshold,
        int    j_start,
        int    j_end,
        int    n_edge           = 100,
        double edge_thr_frac    = 1000.0,
        double hyst_frac        = 0.5,
        int    pre_check_n      = 10,
        double pre_check_frac   = 0.3,
        int    post_check_n     = 3,
        int    confirm_window   = 50,
        bool   check_trailing   = false,
        int    min_below        = 10)
{
    if (j_start < 0 || j_end >= (int)amp.size() || j_start >= j_end)
        return {-1.0, -1.0};

    if (edge_thr_frac < 999.0) {
        double edge_limit = edge_thr_frac * threshold;
        if (edgeMedian(amp, j_start, n_edge) >= edge_limit) return {-1.0,-1.0};
        if (check_trailing &&
            edgeMedian(amp, j_end - n_edge + 1, n_edge) >= edge_limit)
            return {-1.0,-1.0};
    }

    double thr_lo = threshold * hyst_frac;
    bool   armed  = false;
    double t_rise = -1.0;
    int    j_rise = -1;

    for (int j = j_start; j <= j_end; ++j) {
        if (!armed && amp[j] < thr_lo) { armed = true; continue; }
        if (!armed) continue;
        if (amp[j] >= threshold) {
            if (pre_check_n > 0) {
                double pre = edgeMedian(amp, j - pre_check_n, pre_check_n);
                if (pre >= pre_check_frac * threshold) { armed = false; continue; }
            }
            if (j > j_start && amp[j-1] < threshold)
                t_rise = time[j-1] + (threshold - amp[j-1])
                         * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
            else
                t_rise = time[j];
            j_rise = j;
            break;
        }
    }
    if (t_rise < 0) return {-1.0, -1.0};

    double t_fall = -1.0;
    for (int j = j_rise + 1; j <= j_end; ++j) {
        if (amp[j] < threshold) {
            double t_cand = (j > 0 && amp[j-1] >= threshold)
                ? time[j-1] + (threshold - amp[j-1])
                  * (time[j] - time[j-1]) / (amp[j] - amp[j-1])
                : time[j];
            int j_win_end = std::min(j + confirm_window, j_end);
            int n_sub = 0;
            for (int k = j; k <= j_win_end; ++k)
                if (amp[k] < threshold) ++n_sub;
            if (n_sub >= min_below) { t_fall = t_cand; break; }
            j = j_win_end;
        }
    }
    if (t_fall < 0) return {-1.0, -1.0};

    if (post_check_n > 0) {
        int jf = j_rise;
        while (jf <= j_end && time[jf] < t_fall) ++jf;
        if (edgeMedian(amp, jf, post_check_n) >= threshold)
            return {-1.0, -1.0};
    }

    {
        int tail_n     = 50;
        int tail_start = std::max(j_start, j_end + 1 - tail_n);
        int tail_len   = j_end + 1 - tail_start;
        if (tail_len > 0 && edgeMedian(amp, tail_start, tail_len) >= thr_lo)
            return {-1.0, -1.0};
    }

    return {t_rise, t_fall};
}

// ════════════════════════════════════════════════════════════
//  computeTOT_rise_only
//  Calcola SOLO t_rise a soglia arbitraria, senza filtri.
//  Usato per delta_t_lo e delta_t_hi su eventi gia' validati.
//  L'arming e' rilassato (0.3x invece di 0.5x) perche' a soglie
//  basse il segnale potrebbe non scendere abbastanza tra un
//  fotone e l'altro.
// ════════════════════════════════════════════════════════════
static double computeTOT_rise_only(
        const std::vector<double>& time,
        const std::vector<double>& amp,
        double threshold,
        int j_start,
        int j_end)
{
    if (j_start < 0 || j_end >= (int)amp.size() || j_start >= j_end)
        return TRISE_INVALID;
    if (threshold <= 0.0) return TRISE_INVALID;

    double thr_arm = threshold * 0.3;
    bool   armed   = false;

    for (int j = j_start; j <= j_end; ++j) {
        if (!armed) {
            if (amp[j] < thr_arm) armed = true;
            continue;
        }
        if (amp[j] >= threshold) {
            if (j > j_start && amp[j-1] < threshold)
                return time[j-1] + (threshold - amp[j-1])
                       * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
            return time[j];
        }
    }
    return TRISE_INVALID;
}

// ── estimatePE ───────────────────────────────────────────────────────────
// n_pe = round((amp_max - cal.q) / cal.m). Returns -1 if unclassifiable.
static int estimatePE(double amp_max, const CalibResult& cal,
                      double tolerance = 0.4) {
    if (!cal.ok || cal.m <= 0) return -1;
    double n_raw = (amp_max - cal.q) / cal.m;
    int    n_pe  = (int)std::round(n_raw);
    if (n_pe < 0) return 0;
    if (std::abs(n_raw - n_pe) > tolerance) return -1;
    return n_pe;
}

// ── preLaserQuiet ────────────────────────────────────────────────────────
// Scansiona il waveform da BASELINE_END a (t_laser - margin_ns).
// Ritorna false se qualsiasi campione supera quiet_thr:
// indica dark count o afterpulse prima del laser.
static bool preLaserQuiet(const std::vector<double>& time,
                           const std::vector<double>& amp,
                           double t_laser,
                           double quiet_thr,
                           double margin_ns = 1.0)
{
    double t_check_end = t_laser - margin_ns;
    if (t_check_end <= BASELINE_END) return true;
    for (size_t j = 0; j < time.size(); ++j) {
        if (time[j] < BASELINE_END)    continue;
        if (time[j] > t_check_end)     break;
        if (amp[j]  > quiet_thr)       return false;
    }
    return true;
}
