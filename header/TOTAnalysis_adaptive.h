#pragma once
// TOTAnalysis_adaptive.h
// Parametri adattivi per soglie basse

#include "Config.h"
#include "SignalProcessing.h"
#include "Calibration.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

// Parametri adattivi basati su threshold
struct AdaptiveParams {
    double edge_thr_frac;
    double hyst_frac;
    double pre_check_frac;
    int confirm_window;
    int min_below;
    double baseline_max_rms;
};

// Calcola parametri adattivi — rilassa filtri per soglie basse
static AdaptiveParams getAdaptiveParams(double threshold_mV, bool use_filter) {
    AdaptiveParams p;
    
    // Soglia normalizzata: 10 mV = bassa, 100+ mV = alta
    // Fattore di rilassamento: [0.3, 1.0]
    double relax = std::max(0.3, std::min(1.0, threshold_mV / 30.0));
    
    // ── Edge stability ──
    // Bassa soglia → baseline noise diventa significativo
    // Standard: 0.5*thr  →  Adaptive: 0.2*thr per thr=10 mV
    if (use_filter) {
        p.edge_thr_frac = 0.3 + 0.2 * relax;  // [0.3, 0.5]
    } else {
        p.edge_thr_frac = 100.0;  // disabilitato
    }
    
    // ── Hysteresis ──
    // Bassa soglia → segnale piccolo, isteresi stretta rigetta tutto
    // Standard: 0.5*thr  →  Adaptive: 0.25*thr per thr=10 mV
    p.hyst_frac = 0.25 + 0.25 * relax;  // [0.25, 0.5]
    
    // ── Pre-check ──
    // Baseline noise può superare 0.3*thr per thr bassa
    // Standard: 0.3*thr  →  Adaptive: 0.15*thr per thr=10 mV
    p.pre_check_frac = 0.15 + 0.15 * relax;  // [0.15, 0.3]
    
    // ── Falling edge confirmation ──
    // Segnale piccolo → più sensibile a glitch, ma serve conferma più corta
    // Standard: 50 samples, min 10  →  Adaptive: 30 samples, min 6
    if (relax < 0.5) {
        p.confirm_window = 30;
        p.min_below = 6;
    } else {
        p.confirm_window = 50;
        p.min_below = 10;
    }
    
    // ── Baseline RMS ──
    // Soglia bassa → accetta baseline più rumorosa
    // Standard: 4.0 mV no-filter, 2.0 mV filter
    if (use_filter) {
        p.baseline_max_rms = 1.5 + 0.5 * relax;  // [1.5, 2.0]
    } else {
        p.baseline_max_rms = 3.0 + 1.0 * relax;  // [3.0, 4.0]
    }
    
    return p;
}


// computeTOT con parametri adattivi
static std::pair<double,double> computeTOT_adaptive(
    const std::vector<double>& time,
    const std::vector<double>& amp,
    double threshold,
    int j_start, int j_end,
    bool use_filter = true)
{
    AdaptiveParams p = getAdaptiveParams(threshold, use_filter);
    
    const double edge_limit = p.edge_thr_frac * threshold;
    if (edgeMedian(amp, 0, 100) >= edge_limit)
        return {-1.0, -1.0};
    
    const double thr_lo = threshold * p.hyst_frac;
    double t_rise = -1.0, t_fall = -1.0;
    bool armed_rise = false;
    
    // Rising edge con isteresi adattiva
    for (int j = j_start; j <= j_end; ++j) {
        if (!armed_rise) {
            if (amp[j] < thr_lo) armed_rise = true;
        } else {
            if (amp[j] >= threshold) {
                double tr = time[j-1] + (threshold - amp[j-1])
                            * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
                
                // Pre-check adattivo
                int pre_start = std::max(j_start, j - 10);
                int pre_n = j - pre_start;
                if (pre_n > 0 &&
                    edgeMedian(amp, pre_start, pre_n) >= p.pre_check_frac * threshold) {
                    armed_rise = false;
                    continue;
                }
                t_rise = tr;
                break;
            }
        }
    }
    if (t_rise < 0) return {-1.0, -1.0};
    
    // Falling edge con conferma adattiva
    for (int j = j_start; j <= j_end; ++j) {
        if (time[j] <= t_rise) continue;
        if (amp[j] >= threshold) continue;
        
        int j_cand = j;
        double t_cand = (j_cand > 0)
            ? time[j_cand-1] + (threshold - amp[j_cand-1])
              * (time[j_cand] - time[j_cand-1]) / (amp[j_cand] - amp[j_cand-1])
            : time[j_cand];
        
        int j_win_end = std::min(j_cand + p.confirm_window, j_end);
        int n_sub = 0;
        for (int k = j_cand; k <= j_win_end; ++k)
            if (amp[k] < threshold) ++n_sub;
        
        if (n_sub >= p.min_below) {
            t_fall = t_cand;
            break;
        }
        j = j_win_end;
    }
    if (t_fall < 0) return {-1.0, -1.0};
    
    // Post-fall check
    {
        int jf_post = 0;
        for (int j = j_start; j <= j_end; ++j)
            if (time[j] >= t_fall) { jf_post = j; break; }
        
        int post_limit = std::min(jf_post + 50, j_end + 1);
        int post_n = post_limit - jf_post;
        if (post_n > 0 && edgeMedian(amp, jf_post, post_n) >= threshold)
            return {-1.0, -1.0};
    }
    
    // End-of-trace check
    {
        const int tail_n = 50;
        int tail_start = std::max(j_start, j_end + 1 - tail_n);
        int tail_len = j_end + 1 - tail_start;
        if (tail_len > 0 && edgeMedian(amp, tail_start, tail_len) >= thr_lo)
            return {-1.0, -1.0};
    }
    
    return {t_rise, t_fall};
}


// Stampa parametri correnti
static void printAdaptiveParams(double threshold_mV, bool use_filter) {
    AdaptiveParams p = getAdaptiveParams(threshold_mV, use_filter);
    
    std::cout << "\n┌─── ADAPTIVE PARAMETERS ───────────────────────┐\n";
    std::cout << "│ Threshold:         " << std::setw(6) << std::fixed << std::setprecision(2)
              << threshold_mV << " mV                    │\n";
    std::cout << "│ Filter:            " << (use_filter ? "ON " : "OFF") << "                           │\n";
    std::cout << "├───────────────────────────────────────────────┤\n";
    std::cout << "│ edge_thr_frac:     " << std::setw(6) << std::setprecision(3)
              << p.edge_thr_frac << "  (x threshold)         │\n";
    std::cout << "│ hyst_frac:         " << std::setw(6) << p.hyst_frac
              << "  (x threshold)         │\n";
    std::cout << "│ pre_check_frac:    " << std::setw(6) << p.pre_check_frac
              << "  (x threshold)         │\n";
    std::cout << "│ confirm_window:    " << std::setw(6) << p.confirm_window
              << "  samples               │\n";
    std::cout << "│ min_below:         " << std::setw(6) << p.min_below
              << "  samples               │\n";
    std::cout << "│ baseline_max_rms:  " << std::setw(6) << std::setprecision(2)
              << p.baseline_max_rms << " mV                    │\n";
    std::cout << "└───────────────────────────────────────────────┘\n\n";
}
