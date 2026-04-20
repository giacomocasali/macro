#pragma once
// TOTAnalysis_loose.h
// VERSIONE LARGA: filtri rilassati per massimizzare accettanza
// Mantiene integrità fisica minima, accetta più eventi borderline

#include "Config.h"
#include "SignalProcessing.h"
#include "Calibration.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

// Parametri rilassati per massima accettanza
struct LooseParams {
    double edge_thr_frac;    // Standard 0.5 → Loose 1.0 (più permissivo)
    double hyst_frac;        // Standard 0.5 → Loose 0.8 (arma più facilmente)
    double pre_check_frac;   // Standard 0.3 → Loose 0.6 (più permissivo)
    int    confirm_window;   // Standard 50 → Loose 30
    int    min_below;        // Standard 10 → Loose 3 (meno rigido)
    double baseline_max_rms; // Standard 4.0 → Loose 8.0 (accetta baseline sporca)
    int    post_check_n;     // Standard 50 → Loose 20 (check più corto)
    int    tail_n;           // Standard 50 → Loose 20 (check più corto)
};

// Parametri fissi rilassati
static LooseParams getLooseParams(bool use_filter) {
    LooseParams p;
    
    // Edge leading: molto permissivo
    // No filter: disabilitato completamente
    // Filter: 1.0*thr invece di 0.5*thr
    p.edge_thr_frac = use_filter ? 1.0 : 100.0;
    
    // Hysteresis: più facile armare
    // Standard richiede amp < 0.5*thr per armare
    // Loose richiede amp < 0.8*thr (quasi sempre vero)
    p.hyst_frac = 0.8;
    
    // Pre-check: più permissivo
    // Standard rigetta se mediana pre > 0.3*thr
    // Loose rigetta solo se > 0.6*thr
    p.pre_check_frac = 0.6;
    
    // Falling edge: conferma più corta
    p.confirm_window = 30;
    p.min_below = 3;  // solo 3 samples sotto thr invece di 10
    
    // Baseline RMS: accetta baseline più rumorosa
    p.baseline_max_rms = use_filter ? 5.0 : 8.0;
    
    // Post-checks: window più corta
    p.post_check_n = 20;
    p.tail_n = 20;
    
    return p;
}

// computeTOT LOOSE — parametri rilassati
static std::pair<double,double> computeTOT_loose(
    const std::vector<double>& time,
    const std::vector<double>& amp,
    double threshold,
    int j_start, int j_end,
    bool use_filter = true)
{
    LooseParams p = getLooseParams(use_filter);
    
    // Edge stability
    const double edge_limit = p.edge_thr_frac * threshold;
    if (edgeMedian(amp, 0, 100) >= edge_limit)
        return {-1.0, -1.0};
    
    const double thr_lo = threshold * p.hyst_frac;
    double t_rise = -1.0, t_fall = -1.0;
    bool armed_rise = false;
    
    // Rising edge con isteresi larga
    for (int j = j_start; j <= j_end; ++j) {
        if (!armed_rise) {
            if (amp[j] < thr_lo) armed_rise = true;
        } else {
            if (amp[j] >= threshold) {
                double tr = time[j-1] + (threshold - amp[j-1])
                            * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
                
                // Pre-check rilassato
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
    
    // Falling edge con conferma corta
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
    
    // Post-fall check (window corta)
    {
        int jf_post = 0;
        for (int j = j_start; j <= j_end; ++j)
            if (time[j] >= t_fall) { jf_post = j; break; }
        
        int post_limit = std::min(jf_post + p.post_check_n, j_end + 1);
        int post_n = post_limit - jf_post;
        if (post_n > 0 && edgeMedian(amp, jf_post, post_n) >= threshold)
            return {-1.0, -1.0};
    }
    
    // End-of-trace check (window corta)
    {
        int tail_start = std::max(j_start, j_end + 1 - p.tail_n);
        int tail_len = j_end + 1 - tail_start;
        if (tail_len > 0 && edgeMedian(amp, tail_start, tail_len) >= thr_lo)
            return {-1.0, -1.0};
    }
    
    return {t_rise, t_fall};
}

// Stampa parametri
static void printLooseParams(bool use_filter) {
    LooseParams p = getLooseParams(use_filter);
    std::cout << "\n+--- LOOSE PARAMETERS (wide acceptance) ---+\n";
    std::cout << "| Filter:            " << (use_filter ? "ON " : "OFF") << "                |\n";
    std::cout << "| edge_thr_frac:     " << std::setw(6) << p.edge_thr_frac << " (x threshold) |\n";
    std::cout << "| hyst_frac:         " << std::setw(6) << p.hyst_frac << " (x threshold) |\n";
    std::cout << "| pre_check_frac:    " << std::setw(6) << p.pre_check_frac << " (x threshold) |\n";
    std::cout << "| confirm_window:    " << std::setw(6) << p.confirm_window << " samples       |\n";
    std::cout << "| min_below:         " << std::setw(6) << p.min_below << " samples       |\n";
    std::cout << "| baseline_max_rms:  " << std::setw(6) << p.baseline_max_rms << " mV            |\n";
    std::cout << "+-------------------------------------------+\n\n";
}
