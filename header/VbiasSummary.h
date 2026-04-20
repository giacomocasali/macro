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
//  RESOLUTION DB WRITER + σ vs LET canvas
//  Writes sigmaResults to TTree "resolution_scan" consumed by
//  plot_resolution_vs_threshold.cpp.
//  Also draws the canvas: one curve per Vbias, X=LET, Y=sigma.
// ============================================================
#include <numeric>
#include <TFile.h>
#include <TTree.h>
static void saveResolutionDB(
        const std::map<int, std::map<double, std::pair<double,double>>>& sigmaResults,
        OutCtx& ctx)
{
    if (sigmaResults.empty()) return;

    // Write ROOT DB — try ../file_root/ first (matches plot script), fallback to rootDir
    std::string dbDir  = ctx.rootDir;
    std::string dbPath = dbDir + "/../file_root/resolution_scan_results.root";
    gSystem->mkdir((dbDir + "/../file_root").c_str(), true);
    TFile* fDB = TFile::Open(dbPath.c_str(), "RECREATE");
    if (!fDB || fDB->IsZombie()) {
        dbPath = dbDir + "/resolution_scan_results.root";
        fDB = TFile::Open(dbPath.c_str(), "RECREATE");
    }
    if (!fDB || fDB->IsZombie()) {
        std::cerr << "  [ResDB] Cannot write: " << dbPath << "\n";
        delete fDB; return;
    }
    Int_t db_vbias = 0; Double_t db_thr=0, db_sigma=0, db_err=0;
    TTree* tDB = new TTree("resolution_scan", "sigma vs LET");
    tDB->Branch("vbias",        &db_vbias, "vbias/I");
    tDB->Branch("threshold_pe", &db_thr,   "threshold_pe/D");
    tDB->Branch("sigma_ns",     &db_sigma, "sigma_ns/D");
    tDB->Branch("sigma_err_ns", &db_err,   "sigma_err_ns/D");
    for (auto& [vb, letmap] : sigmaResults) {
        db_vbias = vb;
        for (auto& [let, sv] : letmap) {
            db_thr=let; db_sigma=sv.first; db_err=sv.second; tDB->Fill();
        }
    }
    fDB->Write(); fDB->Close(); delete fDB;
    std::cout << "  [ResDB] Saved: " << dbPath << "\n";

    // Canvas: σ vs LET, one TGraphErrors per Vbias
    static const int cols[] = {kRed+1, kAzure+1, kGreen+2, kOrange+7, kMagenta+1, kCyan+2};
    const int nc = (int)(sizeof(cols)/sizeof(cols[0]));

    TCanvas* cRL = new TCanvas("cResVsLET","Timing resolution vs p.e. threshold",1000,750);
    cRL->SetGrid();
    cRL->SetLeftMargin(0.14f); cRL->SetBottomMargin(0.12f);
    cRL->SetRightMargin(PAD_RIGHT); cRL->SetTopMargin(PAD_TOP);

    TLegend* leg = new TLegend(0.62,0.68,0.90,0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetFillColor(0);
    leg->SetTextFont(42); leg->SetTextSize(0.036);

    bool first=true; int ci=0;
    double yMin=1e9, yMax=-1e9;
    // first pass: get y range
    for (auto& [vb, letmap] : sigmaResults)
        for (auto& [let, sv] : letmap) {
            yMin=std::min(yMin, sv.first-sv.second);
            yMax=std::max(yMax, sv.first+sv.second);
        }

    for (auto& [vb, letmap] : sigmaResults) {
        std::vector<std::pair<double,std::pair<double,double>>> pts(letmap.begin(), letmap.end());
        std::sort(pts.begin(), pts.end());
        std::vector<double> vX,vY,vEX,vEY;
        for (auto& [let,sv] : pts) {
            vX.push_back(let); vY.push_back(sv.first);
            vEX.push_back(0.0); vEY.push_back(sv.second);
        }
        TGraphErrors* gr = new TGraphErrors(
            (int)vX.size(), vX.data(), vY.data(), vEX.data(), vEY.data());
        gr->SetName(Form("gr_vbias_%d",vb));
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
        gr->SetMarkerColor(cols[ci%nc]); gr->SetLineColor(cols[ci%nc]); gr->SetLineWidth(2);
        if (first) {
            gr->SetTitle("Timing resolution vs p.e. threshold;"
                         "Threshold (p.e.);#sigma_{#Delta t} (ns)");
            gr->Draw("APL");
            double pad=std::max(0.01*(yMax-yMin),0.005);
            gr->GetYaxis()->SetRangeUser(yMin-pad, yMax+pad);
            gr->GetYaxis()->SetTitleOffset(1.4);
            first=false;
        } else { gr->Draw("PL SAME"); }
        leg->AddEntry(gr, Form("Vbias %d V",vb), "lp");
        ++ci;
    }
    if (!first) leg->Draw();
    cRL->Update(); cRL->Modified();
    ctx.savePNG(cRL, "resolution_vs_threshold.png");
    std::cout << "  [ResDB] Canvas: resolution_vs_threshold.png\n";
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
