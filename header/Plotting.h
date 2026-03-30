#pragma once
// header/Plotting.h

#include "Config.h"
#include "SignalProcessing.h"
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>

// ════════════════════════════════════════════════════════════
//  DRAW SCAN CANVAS
//
//  Single pad, log Y, shared X axis.
//  N(V_thr) [blue] on the left log axis.
//  -dN/dV   [orange] mapped onto the SAME log axis via a linear
//  scale factor, drawn with Draw("PL same").
//  A TGaxis on the right shows the true -dN/dV values.
//  The Gaussian fit [green] is drawn scaled the same way.
//
//  Legend (compact):
//    Trigger: (t_start, t_end) ns
//    mu  ±  err  mV     [fit only]
//    sigma ± err  mV    [fit only]
//    0.5/1/1.5 p.e. lines [fit only, colour via line symbol only]
// ════════════════════════════════════════════════════════════
static void drawScanCanvas(const FileInfo& info,
                           const std::vector<double>& thresholds,
                           const std::vector<double>& counts,
                           const std::vector<double>& thr_der_plot,
                           const std::vector<double>& deriv_plot,
                           const std::vector<double>& err_der_plot,
                           TF1*   fitGaus,
                           double fitMean,    double fitMeanErr,
                           double fitSigma,   double fitSigmaErr,
                           double fitAmp,
                           double fitChi2,    int    fitNdf,
                           double t_trig_start, double t_trig_end) {

    int n_pts      = (int)thresholds.size();
    int n_der_plot = (int)thr_der_plot.size();

    // ── Poisson errors ────────────────────────────────────────
    std::vector<double> x_err_zero(n_pts, 0.0);
    std::vector<double> errors(n_pts);
    for (int k = 0; k < n_pts; ++k)
        errors[k] = (counts[k] > 0) ? std::sqrt(counts[k]) : 0.0;

    // ── Auto X range ──────────────────────────────────────────
    double scan_ymax = *std::max_element(counts.begin(), counts.end());
    double xlo = -10.0;   // fixed wide range — adjust by hand if needed
    double xhi = 120.0;

    // ── Physical range of derivative (ignore near-zero spike) ─
    const double REFMAX_THR_PLOT = 8.0;
    double der_max = 0, der_min_pos = 1e30;
    for (int k = 0; k < n_der_plot; ++k) {
        if (thr_der_plot[k] <= REFMAX_THR_PLOT) continue;
        if (deriv_plot[k] > der_max) der_max = deriv_plot[k];
        if (deriv_plot[k] > 0 && deriv_plot[k] < der_min_pos) der_min_pos = deriv_plot[k];
    }
    if (der_max     <= 0)       der_max     = 1.0;
    if (der_min_pos > der_max)  der_min_pos = der_max * 0.001;

    // ── Shared log Y range ────────────────────────────────────
    // Both curves are displayed on the same log Y axis (Counts).
    // We map the derivative range [der_min_pos, der_max] linearly onto
    // the log-counts range [ylo, yhi] by choosing a scale factor such
    // that der_max maps to scan_ymax * 0.8 (80% of scan peak).
    // Scale: counts_equiv = deriv * derScale
    // The right TGaxis shows the true derivative labels.
    double ylo = 0.5;
    double yhi = 1.0e7;

    double derScale    = (scan_ymax * 0.80) / der_max;
    double der_ylo_map = der_min_pos * derScale * 0.5;
    double der_yhi_map = der_max     * derScale * 1.2;
    // Clamp mapped range inside [ylo, yhi]
    der_ylo_map = std::max(der_ylo_map, ylo);
    der_yhi_map = std::min(der_yhi_map, yhi);

    // Build scaled derivative (clip to ylo for log safety)
    const double LOG_FLOOR = ylo;
    std::vector<double> deriv_scaled(n_der_plot);
    for (int k = 0; k < n_der_plot; ++k) {
        double v = deriv_plot[k] * derScale;
        deriv_scaled[k] = (v > LOG_FLOOR) ? v : LOG_FLOOR;
    }

    // Scaled fit for drawing
    TF1* fitScaled = nullptr;
    if (fitGaus) {
        double sc = derScale;
        TF1*   fg = fitGaus;
        std::string sn = "fitSc_" + info.tag;
        fitScaled = new TF1(sn.c_str(),
            [fg, sc](double* x, double*) {
                double v = fg->Eval(x[0]) * sc;
                return (v > 0.1) ? v : 0.1;
            }, fitGaus->GetXmin(), fitGaus->GetXmax(), 0);
        fitScaled->SetNpx(1000);
        fitScaled->SetLineColor(kGreen+2);
        fitScaled->SetLineWidth(2);
    }

    // ── TGraphErrors ──────────────────────────────────────────
    TGraphErrors* grScan = new TGraphErrors(n_pts,
        thresholds.data(), counts.data(), x_err_zero.data(), errors.data());
    grScan->SetTitle(Form(
        "Threshold scan   V_{bias} = %.0f V,  filter = %.0f"
        ";Threshold (mV);Counts",
        info.vbias, info.filter));
    grScan->SetLineColor(kAzure+1); grScan->SetMarkerColor(kAzure+1);
    grScan->SetMarkerStyle(20); grScan->SetMarkerSize(0.4); grScan->SetLineWidth(2);

    std::vector<double> xerr_z(n_der_plot, 0.0), yerr_z(n_der_plot, 0.0);
    TGraphErrors* grDeriv = new TGraphErrors(n_der_plot,
        thr_der_plot.data(), deriv_scaled.data(), xerr_z.data(), yerr_z.data());
    grDeriv->SetLineColor(kOrange+7); grDeriv->SetMarkerColor(kOrange+7);
    grDeriv->SetMarkerStyle(20); grDeriv->SetMarkerSize(0.4); grDeriv->SetLineWidth(2);

    // ── Canvas: single pad ────────────────────────────────────
    std::string cname = "c_" + info.tag;
    TCanvas* c = new TCanvas(cname.c_str(),
        ("Threshold scan -- " + info.tag).c_str(), 800, 800);
    std::string pn = "pad_" + info.tag;
    // Extra right margin for the TGaxis label
    TPad* pad = new TPad(pn.c_str(), "main", 0.0, 0.0, 1.0, 1.0);
    pad->SetLeftMargin(PAD_LEFT);
    pad->SetRightMargin(PAD_RIGHT);
    pad->SetTopMargin(PAD_TOP);
    pad->SetBottomMargin(PAD_BOTTOM);
    pad->SetTicks(1,1);
    pad->SetLogy();
    pad->Draw(); pad->cd();

    grScan->Draw("APL");
    grScan->GetXaxis()->SetRangeUser(xlo, xhi);
    grScan->GetYaxis()->SetRangeUser(ylo, yhi);
    grScan->GetXaxis()->SetTitleSize(0.050); grScan->GetXaxis()->SetLabelSize(0.045);
    grScan->GetYaxis()->SetTitleSize(0.050); grScan->GetYaxis()->SetLabelSize(0.045);
    grScan->GetYaxis()->SetTitleOffset(1.10);

    // grDeriv and fitScaled kept in memory but not drawn

    // ── p.e. reference lines ──────────────────────────────────
    bool fitOK = (fitGaus && fitMean > 0);
    const double thr_val[3]  = { 0.5*fitMean,    fitMean,    1.5*fitMean    };
    const double thr_err[3]  = { 0.5*fitMeanErr, fitMeanErr, 1.5*fitMeanErr };
    const char*  thr_name[3] = { "0.5 p.e.",     "1 p.e.",   "1.5 p.e."    };
    const int    thr_col[3]  = { kCyan+1, kGreen+2, kMagenta+1 };
    if (fitOK) {
        // Use TGraph with 2 points spanning a huge Y range.
        // Unlike TLine, TGraph is re-clipped at every redraw (including
        // after interactive zoom), so it always fills the full Y axis.
        double y_lo = 1e-3, y_hi = 1e10;
        for (int r = 0; r < 3; ++r) {
            if (thr_val[r] < xlo || thr_val[r] > xhi) continue;
            double gx[2] = { thr_val[r], thr_val[r] };
            double gy[2] = { y_lo, y_hi };
            TGraph* gl = new TGraph(2, gx, gy);
            gl->SetLineColor(thr_col[r]);
            gl->SetLineStyle(2);
            gl->SetLineWidth(2);
            gl->Draw("L same");
        }
    }

    // ── Compact legend ────────────────────────────────────────
    // Legend: top-right corner by default
    double lx2  = 1.0 - PAD_RIGHT - 0.01;
    double lx1  = lx2 - 0.36;
    double ly2  = 1.0 - PAD_TOP - 0.02;
    int nRows   = 1 + (fitOK ? 4 : 0);
    double ly1  = ly2 - nRows * 0.048;

    TLegend* leg = new TLegend(lx1, ly1, lx2, ly2);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);
    leg->SetMargin(0.08);            // tighten left padding

    leg->AddEntry((TObject*)nullptr,
        Form("Trigger: (%.1f, %.1f) ns", t_trig_start, t_trig_end), "");
    if (fitOK) {
        // Significant figures on the error (standard physics convention):
        //   - 1 sig fig, EXCEPT when the leading digit is 1 → keep 2 sig figs
        //     (rounding 1.x to 1 would introduce a ~40% relative error on the error).
        // The value is then reported with the same number of decimal places as the error.
        auto nDecErr = [](double err) -> int {
            if (err <= 0) return 3;
            int mag = (int)std::floor(std::log10(std::abs(err))); // e.g. 0.018 -> -2
            // leading digit
            double leading = std::abs(err) / std::pow(10.0, mag);  // 1.0 <= leading < 10
            int sigfigs = (leading < 2.0) ? 2 : 1;  // 1 sig fig, 2 if leading digit is 1
            int dec = -mag + (sigfigs - 1);
            if (dec < 0) dec = 0;
            if (dec > 6) dec = 6;
            return dec;
        };
        auto fmtVal = [&nDecErr](double val, double err) -> TString {
            int dec = nDecErr(err);
            return TString(Form(Form("%%.%df", dec), val));
        };
        auto fmtErr = [&nDecErr](double err) -> TString {
            int dec = nDecErr(err);
            return TString(Form(Form("%%.%df", dec), err));
        };
        leg->AddEntry((TObject*)nullptr,
            Form("#mu = (%s #pm %s) mV",
                 fmtVal(fitMean,  fitMeanErr).Data(),
                 fmtErr(fitMeanErr).Data()), "");
        leg->AddEntry((TObject*)nullptr,
            Form("#sigma = (%s #pm %s) mV",
                 fmtVal(fitSigma, fitSigmaErr).Data(),
                 fmtErr(fitSigmaErr).Data()), "");
        if (fitNdf > 0)
            leg->AddEntry((TObject*)nullptr,
                Form("#chi^{2}/ndf = %.1f / %d", fitChi2, fitNdf), "");
    }
    leg->Draw();

    // ── Warning box when fit failed ───────────────────────────
    if (!fitOK) {
        TPaveText* warn = new TPaveText(0.12, 0.30, 0.88, 0.55, "NDC");
        warn->SetBorderSize(2); warn->SetFillColor(kYellow);
        warn->SetFillStyle(1001); warn->SetTextFont(62);
        warn->SetTextSize(0.038); warn->SetTextColor(kRed+1);
        warn->SetTextAlign(22);
        warn->AddText("FIT NOT PERFORMED");
        warn->AddText("No narrow p.e. peak found in -dN/dV.");
        warn->AddText("Try a narrower trigger window around the signal peak.");
        warn->Draw("same");
    }

    c->Update(); c->Modified();
    c->SaveAs(("scan_" + info.tag + ".png").c_str());
}

// ════════════════════════════════════════════════════════════
//  LASER TIMING CANVAS
//  Requires a "laser" TTree in rawFile.
//  Top pad  : t_peak_SiPM − t_laser distribution (full waveform).
//  Bottom pad: peak amplitude inside the trigger window.
//  Saves laser_timing_<tag>.png.
// ════════════════════════════════════════════════════════════
static void makeLaserCanvas(TFile*             rawFile,
                            const FileInfo&    info,
                            double             cutoff_MHz,
                            double             t_trig_start,
                            double             t_trig_end,
                            const std::string& tag) {

    TTree* treeLaser = (TTree*)rawFile->Get("laser");
    if (!treeLaser) {
        std::cout << "  [INFO] No 'laser' tree in "
                  << info.path << " — timing canvas skipped.\n";
        return;
    }
    TTree* treeCh1 = (TTree*)rawFile->Get("ch1");
    if (!treeCh1) return;

    const int N = 1024;
    Double_t t1[N], a1[N], tL[N], aL[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);
    treeLaser->SetBranchAddress("time",      tL);
    treeLaser->SetBranchAddress("amplitude", aL);
    treeCh1->GetEntry(0);
    double fs_MHz = 1000.0 / (t1[1] - t1[0]);

    int j_trig_start, j_trig_end;
    triggerWindowIndices(t1, N, t_trig_start, t_trig_end,
                         j_trig_start, j_trig_end, tag);

    double wf_duration = t1[N-1] - t1[0];
    int    timing_bins = (int)(wf_duration / 0.4) + 1;

    TH1D* hTiming = new TH1D(
        ("hTiming_" + tag).c_str(),
        Form("SiPM peak time relative to laser   V_{bias}=%.0fV, filter=%.0f"
             ";t_{peak}^{SiPM} - t_{laser}  (ns);Events", info.vbias, info.filter),
        timing_bins, t1[0], t1[N-1]);
    hTiming->SetDirectory(nullptr);
    hTiming->SetLineColor(kAzure+1); hTiming->SetLineWidth(2);
    hTiming->SetFillColorAlpha(kAzure+1, 0.25);

    TH1D* hAmpTrig = new TH1D(
        ("hAmpTrig_" + tag).c_str(),
        Form("Peak amplitude in trigger window [%.1f, %.1f] ns   "
             "V_{bias}=%.0fV, filter=%.0f;Peak amplitude  (mV);Events",
             t_trig_start, t_trig_end, info.vbias, info.filter),
        200, -5.0, MAX_THR);
    hAmpTrig->SetDirectory(nullptr);
    hAmpTrig->SetLineColor(kRed+1); hAmpTrig->SetLineWidth(2);
    hAmpTrig->SetFillColorAlpha(kRed+1, 0.20);

    Long64_t nEntries = treeCh1->GetEntries();
    long nNoLaser = 0, nInWindow = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 1000 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  [laser canvas] " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        treeLaser->GetEntry(i);

        double t_laser = laserTriggerTime(tL, aL, N, 10.0);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, BASELINE_START, BASELINE_END),
            cutoff_MHz, fs_MHz);

        auto   itMax  = std::max_element(af.begin(), af.end());
        int    i_peak = (int)std::distance(af.begin(), itMax);
        hTiming->Fill(t1[i_peak] - t_laser);

        double A_win = -1e9;
        for (int j = j_trig_start; j <= j_trig_end; ++j)
            if (af[j] > A_win) A_win = af[j];
        if (A_win > -1e8) { hAmpTrig->Fill(A_win); ++nInWindow; }
    }
    std::cout << "\r  [laser canvas] " << nEntries << " / " << nEntries << " — done.\n";
    std::cout << "  No laser trigger: " << nNoLaser
              << " / " << nEntries << "\n"
              << "  In trigger window: " << nInWindow
              << " / " << nEntries << "\n";

    // ── Canvas ───────────────────────────────────────────────
    std::string lcname = "cLaser_" + tag;
    TCanvas* cL = new TCanvas(lcname.c_str(),
                              ("Laser timing — " + tag).c_str(), 950, 1000);

    TPad* lpad1 = new TPad("lpad1", "timing",    0.0, 0.5, 1.0, 1.0);
    TPad* lpad2 = new TPad("lpad2", "amplitude", 0.0, 0.0, 1.0, 0.5);
    for (TPad* p : {lpad1, lpad2}) {
        p->SetLeftMargin(PAD_LEFT);  p->SetRightMargin(PAD_RIGHT);
        p->SetTopMargin(PAD_TOP);    p->SetBottomMargin(PAD_BOTTOM);
        p->SetGrid(); p->SetTicks(1,1); p->Draw();
    }

    lpad1->cd();
    hTiming->GetXaxis()->SetTitleSize(0.055); hTiming->GetXaxis()->SetLabelSize(0.050);
    hTiming->GetYaxis()->SetTitleSize(0.055); hTiming->GetYaxis()->SetLabelSize(0.050);
    hTiming->GetYaxis()->SetTitleOffset(1.1);
    hTiming->Draw("HIST");

    double t_ymax = hTiming->GetMaximum();
    for (double x : {t_trig_start, t_trig_end}) {
        TLine* ll = new TLine(x, 0, x, t_ymax);
        ll->SetLineColor(kOrange+7); ll->SetLineStyle(2); ll->SetLineWidth(2);
        ll->Draw("same");
    }
    TPaveText* ptTrig = new TPaveText(0.14, 0.78, 0.55, 0.88, "NDC");
    ptTrig->SetBorderSize(1); ptTrig->SetFillColor(0); ptTrig->SetFillStyle(1001);
    ptTrig->SetTextFont(42);  ptTrig->SetTextSize(0.038); ptTrig->SetTextAlign(12);
    ptTrig->AddText(Form("Trigger window: [%.1f, %.1f] ns  (orange dashes)",
                         t_trig_start, t_trig_end));
    ptTrig->AddText(Form("No laser trigger: %ld / %lld events", nNoLaser, nEntries));
    ptTrig->Draw();

    lpad2->cd();
    hAmpTrig->GetXaxis()->SetTitleSize(0.055); hAmpTrig->GetXaxis()->SetLabelSize(0.050);
    hAmpTrig->GetYaxis()->SetTitleSize(0.055); hAmpTrig->GetYaxis()->SetLabelSize(0.050);
    hAmpTrig->GetYaxis()->SetTitleOffset(1.1);
    hAmpTrig->Draw("HIST");

    TPaveText* ptAmp = new TPaveText(0.55, 0.78, 0.94, 0.88, "NDC");
    ptAmp->SetBorderSize(1); ptAmp->SetFillColor(0); ptAmp->SetFillStyle(1001);
    ptAmp->SetTextFont(42);  ptAmp->SetTextSize(0.038); ptAmp->SetTextAlign(12);
    ptAmp->AddText(Form("Events in window: %ld / %lld  (%.1f%%)",
                        nInWindow, nEntries,
                        100.0 * nInWindow / (double)nEntries));
    ptAmp->Draw();

    cL->Update(); cL->Modified();
    cL->SaveAs(("laser_timing_" + tag + ".png").c_str());
    std::cout << "  Saved: laser_timing_" << tag << ".png\n";
}
