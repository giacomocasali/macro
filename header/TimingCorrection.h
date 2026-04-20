#pragma once
// TimingCorrection.h
// Time-walk correction via TH2::FitSlicesY (ROOT FitResiduals approach).
//
// Modello fisico time walk SiPM:
//   Δt(TOT) = p0 + p1 / TOT^p2        (power law, fisicamente motivato)
//
// Procedura:
//   1. Riempie TH2D (TOT vs Δt) con sola finestra picco (±3σ iterativo)
//   2. FitSlicesY → TH1D con mean(Δt) per slice TOT (ROOT lo fa con Gauss)
//   3. Fit power law su quei punti via TGraphErrors
//   4. Sottrae evento per evento: Δt_corr = Δt - (f(TOT) - f(TOT_ref))

#include "TOTAnalysis.h"
#include "OutputManager.h"
#include <vector>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TObjArray.h>

enum class TWMethod { NONE, FIT_RESIDUALS };

static TWMethod askTimeWalkMethod() {
    std::cout << "\n--- Time-walk correction ---\n";
    while (true) {
        std::cout << "  Apply time-walk correction (FitSlicesY power law)? [y/n]: " << std::flush;
        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cerr << "\n[ERROR] stdin closed - using TWMethod::NONE.\n";
            return TWMethod::NONE;
        }
        auto b = line.find_first_not_of(" \t\r\n");
        if (b == std::string::npos) continue;
        char c = static_cast<char>(std::tolower(static_cast<unsigned char>(line[b])));
        if (c == 'y') return TWMethod::FIT_RESIDUALS;
        if (c == 'n') return TWMethod::NONE;
        std::cerr << "  [!] Enter y or n.\n";
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  DEPRECATED: correctTimeWalkFitResiduals
//  This function is no longer used. The active time-walk correction path is
//  computeTimeWalkFromCache() in ChunkedHistoFill.h, which uses an exponential
//  model f(TOT) = p0 + p1*exp(-TOT/tau) on streaming cache data.
//  This power-law path (f(TOT) = p0 + p1*TOT^-p2) loads full TH2D in RAM and
//  is kept only for backward compatibility. It should not be called.
// ══════════════════════════════════════════════════════════════════════════
static std::vector<TOTEvent> correctTimeWalkFitResiduals(
        std::vector<TOTEvent>& events,
        double fit_lo, double fit_hi,
        const std::string& tag,
        OutCtx& ctx,
        int n_bins_tot = 30)
{
    std::cerr << "\n  ════════════════════════════════════════════════════════\n"
              << "  [TW] ERROR: FitResiduals power-law path is DEPRECATED.\n"
              << "  [TW] The active correction uses computeTimeWalkFromCache()\n"
              << "  [TW] from ChunkedHistoFill.h with exponential model.\n"
              << "  [TW] This function should not be called.\n"
              << "  ════════════════════════════════════════════════════════\n\n";
    return events;  // no-op
}

// ── Dispatcher ───────────────────────────────────────────────────────────
static std::vector<TOTEvent> applyTimeWalkCorrection(
        std::vector<TOTEvent>& events,
        TWMethod method,
        double fit_lo, double fit_hi,
        const std::string& tag,
        OutCtx& ctx)
{
    switch (method) {
        case TWMethod::FIT_RESIDUALS:
            return correctTimeWalkFitResiduals(events, fit_lo, fit_hi, tag, ctx);
        default:
            return events;
    }
}
