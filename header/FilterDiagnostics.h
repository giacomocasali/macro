#pragma once
// FilterDiagnostics.h
// Produces two diagnostic canvases to help choose the LP filter cutoff
// and verify the laser trigger threshold.
//
// Canvas 1 — filter_diagnostics_<tag>.png  (two panels):
//   LEFT  : raw vs filtered waveform, zoomed to trigger window ±20 ns.
//           Trigger window in orange, 1 p.e. level in green.
//           Good cutoff: filtered edge smooth, amplitude similar to raw.
//   RIGHT : one-sided PSD (log scale). Cutoff = blue dashed line.
//           Good cutoff: line at the knee between signal roll-off and flat noise.
//   Terminal: power fraction below cutoff.  <80%=too low  85-98%=good  >99%=too high
//   peakLoss is measured with the zero-phase filter to avoid artefacts from
//   causal group delay (the event loop uses causal filtering as it should).
//
// Canvas 2 — laser_diagnostics_<tag>.png  (two panels):
//   LEFT  : t_laser distribution over 1000 events.
//   RIGHT : laser peak amplitude distribution. Orange = 10 mV (old default).
//           Green = recommended threshold (20% of median amplitude).
//   cal.laser_thr is updated to the recommended value.

#include "Config.h"
#include "OutputManager.h"
#include "SignalProcessing.h"
#include "ButterworthFilter.h"
#include "Calibration.h"
#include <TComplex.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TVirtualFFT.h>

// ============================================================
static void drawFilterDiagnostics(
        const std::string& firstFile,
        double             cutoff_MHz,
        double             fs_MHz,
        double             t_trig_start,
        double             t_trig_end,
        CalibResult&       cal,          // non-const: we update laser_thr
        const std::string& tag,
        OutCtx&            ctx,
        int                max_search = 2000)
{
    TFile* f = TFile::Open(firstFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[FilterDiag] Cannot open " << firstFile << "\n";
        return;
    }
    TTree* tr = (TTree*)f->Get("ch1");
    if (!tr) { f->Close(); return; }

    const int N = 1024;
    Double_t t_raw[N], a_raw[N];
    tr->SetBranchAddress("time",      t_raw);
    tr->SetBranchAddress("amplitude", a_raw);

    // Threshold: 0.8 p.e. above baseline.  Fallback 5 mV if no calibration.
    double sig_thr = (cal.ok && cal.m > 0) ? cal.q + 0.8 * cal.m : 5.0;

    // --- search for a waveform with signal ----------------------
    Long64_t nSearch    = std::min((Long64_t)max_search, tr->GetEntries());
    Long64_t foundEntry = 0;
    bool     found      = false;

    std::cout << "  [FilterDiag] Scanning first " << nSearch
              << " entries for a signal waveform"
              << " (thr=" << std::fixed << std::setprecision(1)
              << sig_thr << " mV)..." << std::flush;

    for (Long64_t i = 0; i < nSearch; ++i) {
        tr->GetEntry(i);
        std::vector<double> vt(t_raw, t_raw+N), va(a_raw, a_raw+N);
        std::vector<double> vc = correctBaseline(vt, va,
                                                  BASELINE_START, BASELINE_END);
        // Check peak inside trigger window only
        double pk = -1e9;
        for (int j = 0; j < N; ++j)
            if (t_raw[j] >= t_trig_start && t_raw[j] <= t_trig_end)
                if (vc[j] > pk) pk = vc[j];
        if (pk >= sig_thr) {
            foundEntry = i;
            found      = true;
            break;
        }
    }

    if (found)
        std::cout << " found at entry " << foundEntry << ".\n";
    else {
        std::cout << " not found. Using entry 0 (noise waveform).\n";
        foundEntry = 0;
    }

    // --- read the chosen waveform -------------------------------
    tr->GetEntry(foundEntry);
    f->Close();

    std::vector<double> vt(t_raw, t_raw+N), va(a_raw, a_raw+N);
    std::vector<double> va_corr = correctBaseline(vt, va,
                                                   BASELINE_START, BASELINE_END);
    // Causal filter for the waveform display (consistent with event loop)
    std::vector<double> va_filt = butterworthLowPass(va_corr, cutoff_MHz, fs_MHz);
    // Zero-phase filter for peakLoss only: avoids causal group-delay artefact
    // (the causal filter shifts the peak forward, making it appear lower)
    std::vector<double> va_filt_zp = butterworthLowPassZP(va_corr, cutoff_MHz, fs_MHz);

    double peakRaw  = *std::max_element(va_corr.begin(),    va_corr.end());
    double peakFilt = *std::max_element(va_filt_zp.begin(), va_filt_zp.end());
    double peakLoss = (peakRaw > 0)
                      ? (peakRaw - peakFilt) / peakRaw * 100.0 : 0.0;

    // --- Canvas -------------------------------------------------
    TCanvas* cFD = new TCanvas(
        Form("cFilterDiag_%s", tag.c_str()),
        Form("Filter diagnostics  cutoff=%.0f MHz  %s",
             cutoff_MHz, tag.c_str()),
        1400, 580);
    cFD->Divide(2, 1);

    // ---- LEFT: waveform ----------------------------------------
    cFD->cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(PAD_LEFT);     gPad->SetRightMargin(0.05f);
    gPad->SetBottomMargin(PAD_BOTTOM); gPad->SetTopMargin(PAD_TOP);

    TGraph* grRaw  = new TGraph(N, t_raw, va_corr.data());
    TGraph* grFilt = new TGraph(N, t_raw, va_filt.data());
    grRaw ->SetLineColor(kGray+1);  grRaw ->SetLineWidth(1);
    grFilt->SetLineColor(kAzure+1); grFilt->SetLineWidth(2);

    grRaw->SetTitle(
        Form("Waveform entry %lld  (peak loss %.1f%%)  %s"
             ";Time (ns);Amplitude (mV)",
             foundEntry, peakLoss, tag.c_str()));
    grRaw ->Draw("AL");
    grFilt->Draw("L same");

    // Zoom to trigger window +/- 20 ns
    double xlo_w = std::max((double)t_raw[0],   t_trig_start - 20.0);
    double xhi_w = std::min((double)t_raw[N-1], t_trig_end   + 20.0);
    grRaw->GetXaxis()->SetRangeUser(xlo_w, xhi_w);

    // Trigger window
    double ylo = grRaw->GetYaxis()->GetXmin();
    double yhi = grRaw->GetYaxis()->GetXmax();
    TLine* lS = new TLine(t_trig_start, ylo, t_trig_start, yhi);
    TLine* lE = new TLine(t_trig_end,   ylo, t_trig_end,   yhi);
    for (TLine* l : {lS, lE}) {
        l->SetLineColor(kOrange+7); l->SetLineStyle(2); l->SetLineWidth(2);
        l->Draw("same");
    }

    // 1 p.e. level
    TLine* l1pe = nullptr;
    if (cal.ok && cal.m > 0) {
        double y1pe = cal.q + cal.m;
        l1pe = new TLine(xlo_w, y1pe, xhi_w, y1pe);
        l1pe->SetLineColor(kGreen+2); l1pe->SetLineStyle(3); l1pe->SetLineWidth(1);
        l1pe->Draw("same");
    }

    TLegend* legW = new TLegend(0.14, 0.72, 0.62, 0.88);
    legW->SetBorderSize(1); legW->SetFillColor(0);
    legW->SetTextFont(42);  legW->SetTextSize(0.034);
    legW->AddEntry(grRaw,  "Raw (baseline corrected)", "l");
    legW->AddEntry(grFilt, Form("Filtered  %.0f MHz", cutoff_MHz), "l");
    legW->AddEntry(lS,     "Trigger window", "l");
    if (l1pe && cal.ok)
        legW->AddEntry(l1pe, Form("1 p.e. = %.1f mV", cal.q+cal.m), "l");
    legW->Draw();

    // ---- RIGHT: PSD --------------------------------------------
    cFD->cd(2);
    gPad->SetGrid(); gPad->SetLogy(1);
    gPad->SetLeftMargin(PAD_LEFT);     gPad->SetRightMargin(0.05f);
    gPad->SetBottomMargin(PAD_BOTTOM); gPad->SetTopMargin(PAD_TOP);

    int nFFT = N;
    TVirtualFFT* fft = TVirtualFFT::FFT(1, &nFFT, "R2C M");
    for (int i = 0; i < N; ++i) fft->SetPoint(i, va_corr[i]);
    fft->Transform();

    int    nFreq  = N/2 + 1;
    double df     = fs_MHz / N;
    double totPow = 0, cutPow = 0;
    std::vector<double> vFreq(nFreq), vPSD(nFreq);
    for (int k = 0; k < nFreq; ++k) {
        Double_t re, im;
        fft->GetPointComplex(k, re, im);
        double pwr = (re*re + im*im) / ((double)N*N);
        if (k > 0 && k < nFreq-1) pwr *= 2.0;
        vFreq[k] = k * df;
        vPSD [k] = (pwr > 0) ? pwr : 1e-30;
        totPow  += vPSD[k];
        if (vFreq[k] <= cutoff_MHz) cutPow += vPSD[k];
    }
    delete fft;
    double frac = (totPow > 0) ? cutPow / totPow * 100.0 : 0.0;

    // Noise floor: median of the upper-frequency half
    std::vector<double> upperHalf(vPSD.begin() + nFreq/2, vPSD.end());
    std::nth_element(upperHalf.begin(),
                     upperHalf.begin() + upperHalf.size()/2,
                     upperHalf.end());
    double noiseFloor = upperHalf[upperHalf.size()/2];

    TGraph* grPSD = new TGraph(nFreq, vFreq.data(), vPSD.data());
    grPSD->SetLineColor(kRed+1); grPSD->SetLineWidth(2);
    grPSD->SetTitle(Form("Power Spectral Density  %s"
                         ";Frequency (MHz);PSD (arb.)", tag.c_str()));
    grPSD->Draw("AL");

    double psdMin = *std::min_element(vPSD.begin(), vPSD.end());
    double psdMax = *std::max_element(vPSD.begin(), vPSD.end());
    TLine* lCut = new TLine(cutoff_MHz, psdMin, cutoff_MHz, psdMax);
    lCut->SetLineColor(kAzure+1); lCut->SetLineStyle(2); lCut->SetLineWidth(2);
    lCut->Draw("same");

    TLine* lNF = new TLine(0, noiseFloor, fs_MHz/2.0, noiseFloor);
    lNF->SetLineColor(kGray+2); lNF->SetLineStyle(4); lNF->SetLineWidth(1);
    lNF->Draw("same");

    const char* verdict =
        (frac < 80.0) ? "WARNING: too low (cuts signal)"  :
        (frac > 99.0) ? "WARNING: too high (noise passes)" :
                        "Cutoff looks reasonable";

    TPaveText* pt = new TPaveText(0.36, 0.72, 0.93, 0.88, "NDC");
    pt->SetBorderSize(1); pt->SetFillColor(0);
    pt->SetTextFont(42);  pt->SetTextSize(0.034);
    pt->AddText(Form("Cutoff = %.0f MHz", cutoff_MHz));
    pt->AddText(Form("Power below cutoff: %.1f%%", frac));
    pt->AddText(Form("Peak amplitude loss: %.1f%%", peakLoss));
    pt->AddText(verdict);
    pt->Draw();

    TLegend* legF = new TLegend(0.36, 0.58, 0.93, 0.71);
    legF->SetBorderSize(1); legF->SetFillColor(0);
    legF->SetTextFont(42);  legF->SetTextSize(0.032);
    legF->AddEntry(grPSD, "PSD (signal waveform)", "l");
    legF->AddEntry(lCut,  Form("Cutoff %.0f MHz", cutoff_MHz), "l");
    legF->AddEntry(lNF,   "Noise floor estimate", "l");
    legF->Draw();

    cFD->Update(); cFD->Modified();
    ctx.savePNG(cFD, Form("filter_diagnostics_%s.png", tag.c_str()));

    // --- terminal summary ---------------------------------------
    std::cout << "  [FilterDiag] Cutoff=" << cutoff_MHz << " MHz"
              << "  power_below=" << std::fixed << std::setprecision(1)
              << frac << "%"
              << "  peak_loss=" << peakLoss << "%"
              << "\n  --> " << verdict << "\n";

    // ── Laser diagnostics ────────────────────────────────────
    // Pass 1: estimate median laser amplitude from first 500 events.
    // Pass 2: collect t_laser over 1000 events to check timing stability.
    // Sets cal.laser_thr = 20% of median amplitude (robust threshold).
    {
        TFile* fL = TFile::Open(firstFile.c_str(), "READ");
        if (fL && !fL->IsZombie()) {
            TTree* trL = (TTree*)fL->Get("laser");
            if (trL) {
                Double_t tL[N], aL[N];
                trL->SetBranchAddress("time",      tL);
                trL->SetBranchAddress("amplitude", aL);

                // Pass 1: estimate laser amplitude from first 500 events
                std::vector<double> laserPeaks;
                Long64_t nScan = std::min((Long64_t)500, trL->GetEntries());
                for (Long64_t i = 0; i < nScan; ++i) {
                    trL->GetEntry(i);
                    // Baseline correct laser channel
                    std::vector<double> pre;
                    for (int j = 0; j < N; ++j)
                        if (tL[j] < BASELINE_END) pre.push_back(aL[j]);
                    double lOffset = 0;
                    if (!pre.empty()) {
                        std::vector<double> tmp = pre;
                        std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
                        lOffset = tmp[tmp.size()/2];
                    }
                    double pk = -1e9;
                    for (int j = 0; j < N; ++j)
                        if (aL[j] - lOffset > pk) pk = aL[j] - lOffset;
                    if (pk > 0) laserPeaks.push_back(pk);
                }

                double medLaserAmp = 0;
                if (!laserPeaks.empty()) {
                    std::nth_element(laserPeaks.begin(),
                                     laserPeaks.begin() + laserPeaks.size()/2,
                                     laserPeaks.end());
                    medLaserAmp = laserPeaks[laserPeaks.size()/2];
                }

                // Recommended threshold: 20% of median laser amplitude
                // (robust against noise, well below signal peak)
                double recThr = medLaserAmp * 0.20;

                // Update cal so the event loop uses the correct threshold
                cal.laser_thr = recThr;

                std::cout << "  [LaserDiag] Median laser amplitude: "
                          << std::fixed << std::setprecision(1)
                          << medLaserAmp << " mV\n"
                          << "  [LaserDiag] Threshold set to:      "
                          << recThr << " mV (20% of peak)\n";

                // Pass 2: collect t_laser over first 1000 events
                // using the CORRECT threshold (recThr) to show real jitter
                TH1D* hLT = new TH1D(Form("hLaserTime_%s", tag.c_str()),
                    Form("Laser trigger time distribution   %s"
                         ";t_{laser} (ns);Events", tag.c_str()),
                    512, 0.0, 204.8);
                hLT->SetDirectory(nullptr);

                Long64_t nDiag = std::min((Long64_t)1000, trL->GetEntries());
                int nFound = 0, nMissed = 0;
                for (Long64_t i = 0; i < nDiag; ++i) {
                    trL->GetEntry(i);
                    double t = laserTriggerTime(tL, aL, N, recThr);  // use recThr
                    if (t > -900.0) { hLT->Fill(t); ++nFound; }
                    else              ++nMissed;
                }

                int bpk = hLT->GetMaximumBin();
                double tPk = hLT->GetBinCenter(bpk);
                double tRMS = hLT->GetRMS();

                std::cout << "  [LaserDiag] t_laser peak: "
                          << std::setprecision(2) << tPk << " ns"
                          << "   RMS=" << tRMS << " ns"
                          << "   found=" << nFound
                          << "  missed=" << nMissed << "\n";
                if (tRMS > 1.0)
                    std::cout << "  [LaserDiag] NOTE: t_laser RMS=" << tRMS
                              << " ns (measured with thr=" << recThr << " mV)\n";
                else
                    std::cout << "  [LaserDiag] t_laser jitter looks good (<1 ns).\n";

                // Draw laser timing canvas
                TCanvas* cLD = new TCanvas(
                    Form("cLaserDiag_%s", tag.c_str()),
                    Form("Laser diagnostics   %s", tag.c_str()), 900, 450);
                cLD->Divide(2, 1);

                // Left: laser trigger time distribution
                cLD->cd(1); gPad->SetGrid();
                gPad->SetLeftMargin(PAD_LEFT); gPad->SetRightMargin(0.05f);
                gPad->SetBottomMargin(PAD_BOTTOM); gPad->SetTopMargin(PAD_TOP);
                hLT->SetLineColor(kAzure+1); hLT->SetLineWidth(2);
                // Zoom around peak
                hLT->GetXaxis()->SetRangeUser(tPk - 20.0, tPk + 20.0);
                hLT->Draw("HIST");
                {
                    TPaveText* pt2 = new TPaveText(0.14, 0.75, 0.60, 0.88, "NDC");
                    pt2->SetBorderSize(1); pt2->SetFillColor(0);
                    pt2->SetTextFont(42);  pt2->SetTextSize(0.036);
                    pt2->AddText(Form("Peak at %.2f ns", tPk));
                    pt2->AddText(Form("RMS = %.2f ns", tRMS));
                    pt2->AddText(Form("Found %d / %d events", nFound, (int)nDiag));
                    pt2->Draw();
                }

                // Right: laser amplitude distribution
                TH1D* hLA = new TH1D(Form("hLaserAmp_%s", tag.c_str()),
                    Form("Laser peak amplitude   %s;Amplitude (mV);Events",
                         tag.c_str()),
                    100, 0, medLaserAmp * 2.5);
                hLA->SetDirectory(nullptr);
                for (double pk : laserPeaks) hLA->Fill(pk);

                cLD->cd(2); gPad->SetGrid();
                gPad->SetLeftMargin(PAD_LEFT); gPad->SetRightMargin(0.05f);
                gPad->SetBottomMargin(PAD_BOTTOM); gPad->SetTopMargin(PAD_TOP);
                hLA->SetLineColor(kRed+1); hLA->SetLineWidth(2);
                hLA->Draw("HIST");
                {
                    TLine* lThr = new TLine(10.0, 0, 10.0, hLA->GetMaximum());
                    lThr->SetLineColor(kOrange+7); lThr->SetLineStyle(2);
                    lThr->SetLineWidth(2); lThr->Draw("same");
                    TLine* lRec = new TLine(recThr, 0, recThr, hLA->GetMaximum());
                    lRec->SetLineColor(kGreen+2); lRec->SetLineStyle(2);
                    lRec->SetLineWidth(2); lRec->Draw("same");
                    TPaveText* pt3 = new TPaveText(0.50, 0.72, 0.93, 0.88, "NDC");
                    pt3->SetBorderSize(1); pt3->SetFillColor(0);
                    pt3->SetTextFont(42);  pt3->SetTextSize(0.034);
                    pt3->AddText(Form("Median amp = %.1f mV", medLaserAmp));
                    pt3->AddText(Form("Current thr = 10.0 mV (orange)"));
                    pt3->AddText(Form("Rec. thr = %.1f mV (green)", recThr));
                    pt3->Draw();
                }

                cLD->Update(); cLD->Modified();
                ctx.savePNG(cLD, Form("laser_diagnostics_%s.png", tag.c_str()));
            }
            fL->Close();
        }
    }
}
