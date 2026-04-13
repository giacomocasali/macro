#pragma once
// header/VbiasSummary.h
// Two plotting functions used only by the v4 multi-run macro:
//   drawMultiLETOverlay()  -- fine-binned overlay of all LET thresholds
//   drawVbiasSummary()     -- sigma(Delta_t) vs Vbias for each LET

#include "Config.h"
#include "OutputManager.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

// ============================================================
//  MULTI-LET OVERLAY
//  Projects each TH2D onto the Delta_t axis with fine binning
//  (0.02 ns/bin) so the peaks are properly resolved, then
//  normalises each histogram to its own maximum within [fit_lo,
//  fit_hi] for a fair visual comparison between thresholds.
// ============================================================
static void drawMultiLETOverlay(
        const std::vector<double>& fracs_pe,
        const std::vector<TH2D*>&  maps,
        const std::string&         tag,
        double fit_lo, double fit_hi,
        OutCtx& ctx)
{
    if (maps.empty()) return;
    static const int cols[] = {kAzure+1, kRed+1, kGreen+2,
                                kOrange+7, kMagenta+1, kCyan+2};
    const int nc = (int)(sizeof(cols) / sizeof(cols[0]));

    int nBins = std::max(400, (int)std::round((fit_hi - fit_lo) / 0.02));

    TCanvas* cOvl = new TCanvas(
        Form("cMultiLET_%s", tag.c_str()),
        Form("Multi-LET #Deltat overlay -- %s", tag.c_str()), 950, 600);
    cOvl->SetGrid();
    cOvl->SetLeftMargin(PAD_LEFT);    cOvl->SetRightMargin(PAD_RIGHT);
    cOvl->SetBottomMargin(PAD_BOTTOM); cOvl->SetTopMargin(PAD_TOP);

    TLegend* leg = new TLegend(0.62, 0.62, 0.93, 0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetFillColor(0);
    leg->SetTextFont(42);  leg->SetTextSize(0.036);

    bool   first     = true;
    double globalMax = 0.0;
    std::vector<TH1D*> hists;

    // Pass 1: build fine-binned projections and normalise
    for (int i = 0; i < (int)maps.size(); ++i) {
        if (!maps[i]) { hists.push_back(nullptr); continue; }

        TH1D* hP = new TH1D(
            Form("hOvlFine_%s_%d", tag.c_str(), i),
            Form("Multi-LET #Deltat overlay (%s);#Deltat (ns);Normalised counts",
                 tag.c_str()),
            nBins, fit_lo - 1.0, fit_hi + 1.0);
        hP->SetDirectory(nullptr);

        // Fill via ProjectionY (true 0.25 ns/bin resolution of the TH2D),
        // then redistribute uniformly into the finer hP bins.
        // Using GetBinCenter(by) directly would fake 0.02 ns/bin resolution
        // while the actual data sits on a 0.25 ns grid.
        {
            TH1D* hPY = maps[i]->ProjectionY(
                Form("hOvlPY_%s_%d", tag.c_str(), i), 1, maps[i]->GetNbinsX());
            hPY->SetDirectory(nullptr);
            double wPY   = hPY->GetBinWidth(1);
            double wProj = hP->GetBinWidth(1);
            int    nSub  = std::max(1, (int)std::ceil(wPY / wProj));
            for (int by = 1; by <= hPY->GetNbinsX(); ++by) {
                double val = hPY->GetBinContent(by);
                if (val <= 0) continue;
                double xc   = hPY->GetBinCenter(by);
                double step = wPY / nSub;
                double frac = val / nSub;
                for (int s = 0; s < nSub; ++s)
                    hP->Fill(xc - wPY/2.0 + (s + 0.5)*step, frac);
            }
            delete hPY;
        }

        // Normalise to peak within fit window
        int    b1 = hP->FindBin(fit_lo), b2 = hP->FindBin(fit_hi);
        double pk = 0;
        for (int b = b1; b <= b2; ++b)
            if (hP->GetBinContent(b) > pk) pk = hP->GetBinContent(b);
        if (pk > 0) hP->Scale(1.0 / pk);

        hP->SetLineColor(cols[i % nc]);
        hP->SetLineWidth(2);
        hP->GetXaxis()->SetRangeUser(fit_lo - 3.0, fit_hi + 3.0);
        hists.push_back(hP);
        if (hP->GetMaximum() > globalMax) globalMax = hP->GetMaximum();
    }

    // Pass 2: draw
    for (int i = 0; i < (int)hists.size(); ++i) {
        if (!hists[i]) continue;
        hists[i]->SetMaximum(globalMax * 1.15);
        if (first) { hists[i]->Draw("HIST"); first = false; }
        else         hists[i]->Draw("HIST SAME");
        leg->AddEntry(hists[i], Form("LET = %.2f p.e.", fracs_pe[i]), "l");
    }
    if (!first) leg->Draw();

    cOvl->Update(); cOvl->Modified();
    ctx.savePNG(cOvl, Form("tot_multiLET_%s.png", tag.c_str()));
}

// ============================================================
//  CROSS-VBIAS SUMMARY
//  One TGraphErrors per LET: sigma(Delta_t) [ns] vs Vbias [V].
//  Only produced when >= 2 Vbias values have been analysed.
// ============================================================
static void drawVbiasSummary(
        const std::map<int, std::map<double, std::pair<double,double>>>& sigmaResults,
        OutCtx& ctx)
{
    if (sigmaResults.size() < 2) return;

    // Collect the union of all LET values present
    std::vector<double> allLETs;
    for (auto& [vb, letmap] : sigmaResults)
        for (auto& [let, sv] : letmap)
            if (std::find(allLETs.begin(), allLETs.end(), let) == allLETs.end())
                allLETs.push_back(let);
    std::sort(allLETs.begin(), allLETs.end());

    static const int cols[] = {kAzure+1, kRed+1, kGreen+2,
                                kOrange+7, kMagenta+1, kCyan+2};
    const int nc = (int)(sizeof(cols) / sizeof(cols[0]));

    TCanvas* cSum = new TCanvas("cVbiasSummary",
        "#sigma_{#Deltat} vs Vbias (all LET)", 900, 600);
    cSum->SetGrid();
    cSum->SetLeftMargin(PAD_LEFT);    cSum->SetRightMargin(PAD_RIGHT);
    cSum->SetBottomMargin(PAD_BOTTOM); cSum->SetTopMargin(PAD_TOP);

    TLegend* leg = new TLegend(0.62, 0.62, 0.93, 0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetFillColor(0);
    leg->SetTextFont(42);  leg->SetTextSize(0.036);

    bool first = true;
    for (int li = 0; li < (int)allLETs.size(); ++li) {
        double let = allLETs[li];
        std::vector<double> vVb, vSig, vSigErr;
        for (auto& [vb, letmap] : sigmaResults) {
            auto it = letmap.find(let);
            if (it == letmap.end()) continue;
            vVb    .push_back((double)vb);
            vSig   .push_back(it->second.first);
            vSigErr.push_back(it->second.second);
        }
        if (vVb.empty()) continue;

        TGraphErrors* gr = new TGraphErrors(
            (int)vVb.size(), vVb.data(), vSig.data(),
            nullptr, vSigErr.data());
        gr->SetTitle(";Vbias (V);#sigma_{#Deltat} (ns)");
        gr->SetMarkerStyle(20 + li % 10);
        gr->SetMarkerColor(cols[li % nc]);
        gr->SetLineColor  (cols[li % nc]);
        gr->SetLineWidth(2);
        if (first) { gr->Draw("APL"); first = false; }
        else         gr->Draw("PL SAME");
        leg->AddEntry(gr, Form("LET=%.2f p.e.", let), "pl");
    }
    if (!first) leg->Draw();
    cSum->Update(); cSum->Modified();
    ctx.savePNG(cSum, "summary_sigma_vs_vbias.png");
    std::cout << "\n  [SUMMARY] sigma vs Vbias saved.\n";
}
