#pragma once
// TimingCorrection.h
// Time-walk correction for Delta_t.
//
// Time walk: a signal barely above the LET threshold crosses it later than a
// large signal, so Delta_t depends on signal amplitude (or TOT).
// We correct event-by-event by subtracting a fitted model of this dependence.
//
// Method A — Empirical (recommended for low stats):
//   Bin events by TOT, compute Mean(Δt) per bin, fit with pol2.
//   Reference point = TOT minimum (least time walk = closest to true photon time).
//   Subtract: delta_t_corr = delta_t - [fTW(tot) - fTW(tot_min)]
//
// Method B — Per p.e.:
//   Compute Mean(Δt) separately for each n_pe bin, subtract the per-bin mean.
//   Reference point = global mean of per-bin means (preserves absolute scale).

#include "TOTAnalysis.h"
#include "OutputManager.h"
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPaveText.h>

enum class TWMethod { NONE, EMPIRICAL, AMPLITUDE };

static TWMethod askTimeWalkMethod() {
    std::cout << "\n--- Time walk correction ---\n"
              << "  [0] None\n"
              << "  [1] Empirical: fit Mean(Δt) vs TOT, subtract per event\n"
              << "  [2] Amplitude: subtract mean Δt per p.e. bin\n"
              << "  Choice: ";
    int choice = 0; std::cin >> choice;
    switch (choice) {
        case 1: return TWMethod::EMPIRICAL;
        case 2: return TWMethod::AMPLITUDE;
        default: return TWMethod::NONE;
    }
}

// ════════════════════════════════════════════════════════════
//  Method A — Empirical
// ════════════════════════════════════════════════════════════
static std::vector<TOTEvent> correctTimeWalkEmpirical(
        std::vector<TOTEvent>& events,
        double fit_lo, double fit_hi,
        const std::string& tag,
        OutCtx& ctx,
        int n_bins = 20)
{
    if (events.empty()) return events;

    double tot_min = 1e9, tot_max = -1e9;
    for (auto& e : events) {
        if (e.tot < tot_min) tot_min = e.tot;
        if (e.tot > tot_max) tot_max = e.tot;
    }
    double bin_w = (tot_max - tot_min) / n_bins;
    if (bin_w <= 0) return events;

    std::vector<double> vTOT, vMean, vMeanErr;
    for (int b = 0; b < n_bins; ++b) {
        double lo = tot_min + b * bin_w, hi = lo + bin_w;
        double sum = 0, sum2 = 0; int cnt = 0;
        for (auto& e : events) {
            if (e.tot < lo || e.tot >= hi) continue;
            if (e.delta_t < fit_lo || e.delta_t > fit_hi) continue;
            sum += e.delta_t; sum2 += e.delta_t * e.delta_t; ++cnt;
        }
        if (cnt < 5) continue;
        double mean = sum / cnt;
        double rms  = std::sqrt(std::max(sum2/cnt - mean*mean, 0.0));
        vTOT.push_back(0.5*(lo+hi));
        vMean.push_back(mean);
        vMeanErr.push_back(rms / std::sqrt(cnt));
    }
    if (vTOT.size() < 3) {
        std::cout << "  [TWA] Not enough points for fit — skipping.\n";
        return events;
    }

    TGraphErrors* grMean = new TGraphErrors(
        vTOT.size(), vTOT.data(), vMean.data(), nullptr, vMeanErr.data());

    TF1* fTW = new TF1(Form("fTW_emp_%s", tag.c_str()), "pol2", vTOT[0], tot_max);
    grMean->Fit(fTW, "RQ");
    std::cout << std::scientific << std::setprecision(4)
              << "  TW empirical: Dt(TOT) = "
              << fTW->GetParameter(0) << " + "
              << fTW->GetParameter(1) << "*TOT + "
              << fTW->GetParameter(2) << "*TOT^2\n"
              << std::defaultfloat;

    TCanvas* cTW = new TCanvas(Form("cTW_emp_%s", tag.c_str()),
        Form("Time walk (empirical) — %s", tag.c_str()), 900, 500);
    cTW->SetGrid();
    grMean->SetMarkerStyle(20); grMean->SetMarkerColor(kAzure+1);
    grMean->SetTitle(Form("Mean #Deltat vs TOT   %s;TOT (ns);Mean #Deltat (ns)",
                          tag.c_str()));
    grMean->Draw("AP");
    fTW->SetLineColor(kRed+1); fTW->Draw("same");
    cTW->Update(); cTW->Modified();
    cTW->SaveAs(ctx.png(Form("tw_empirical_%s.png", tag.c_str())).c_str());

    // Subtract: reference = fit value at tot_min (minimum time walk)
    std::vector<TOTEvent> corrected = events;
    double ref = fTW->Eval(tot_min);
    for (auto& e : corrected)
        e.delta_t -= (fTW->Eval(e.tot) - ref);
    return corrected;
}

// ════════════════════════════════════════════════════════════
//  Method B — Per p.e. bin
// ════════════════════════════════════════════════════════════
static std::vector<TOTEvent> correctTimeWalkAmplitude(
        std::vector<TOTEvent>& events,
        double fit_lo, double fit_hi,
        const std::string& tag,
        OutCtx& ctx)
{
    if (events.empty()) return events;

    std::map<int, std::vector<double>> byPE;
    for (auto& e : events) {
        if (e.n_pe < 0) continue;
        if (e.delta_t < fit_lo || e.delta_t > fit_hi) continue;
        byPE[e.n_pe].push_back(e.delta_t);
    }

    std::map<int, double> meanByPE;
    for (auto& [npe, vals] : byPE) {
        if (vals.size() < 5) continue;
        double sum = 0;
        for (double v : vals) sum += v;
        meanByPE[npe] = sum / vals.size();
    }
    if (meanByPE.empty()) {
        std::cerr << "  TW amplitude: not enough events per p.e. bin — skipped.\n";
        return events;
    }

    double ref = 0; int nref = 0;
    for (auto& [npe, m] : meanByPE) { ref += m; ++nref; }
    ref /= nref;

    std::vector<double> vPE, vMean;
    for (auto& [npe, m] : meanByPE) { vPE.push_back(npe); vMean.push_back(m); }
    TGraphErrors* grPE = new TGraphErrors(
        vPE.size(), vPE.data(), vMean.data(), nullptr, nullptr);
    TCanvas* cTW = new TCanvas(Form("cTW_amp_%s", tag.c_str()),
        Form("Time walk (amplitude) — %s", tag.c_str()), 900, 500);
    cTW->SetGrid();
    grPE->SetMarkerStyle(20); grPE->SetMarkerColor(kOrange+7);
    grPE->SetTitle(Form("Mean #Deltat per p.e.   %s;p.e.;Mean #Deltat (ns)",
                        tag.c_str()));
    grPE->Draw("AP");
    cTW->Update(); cTW->Modified();
    cTW->SaveAs(ctx.png(Form("tw_amplitude_%s.png", tag.c_str())).c_str());

    std::vector<TOTEvent> corrected = events;
    for (auto& e : corrected) {
        if (e.n_pe < 0) continue;
        auto it = meanByPE.find(e.n_pe);
        if (it == meanByPE.end()) continue;
        e.delta_t -= (it->second - ref);
    }
    return corrected;
}

// ── Dispatcher ───────────────────────────────────────────────
static std::vector<TOTEvent> applyTimeWalkCorrection(
        std::vector<TOTEvent>& events,
        TWMethod method,
        double fit_lo, double fit_hi,
        const std::string& tag,
        OutCtx& ctx)
{
    switch (method) {
        case TWMethod::EMPIRICAL:
            return correctTimeWalkEmpirical(events, fit_lo, fit_hi, tag, ctx);
        case TWMethod::AMPLITUDE:
            return correctTimeWalkAmplitude(events, fit_lo, fit_hi, tag, ctx);
        default:
            return events;
    }
}
