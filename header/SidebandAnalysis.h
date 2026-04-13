#pragma once
// header/SidebandAnalysis.h
// S+B fit (Gauss + flat background) and sideband subtraction on Δt.
// SB_GAP=1 ns separates sidebands from signal window.
// ============================================================

#include "TOTAnalysis.h"
#include "OutputManager.h"
#include "Config.h"

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TMath.h>

namespace SB {
    constexpr double DEFAULT_SB_WIDTH = 8.0;   // sideband width per side [ns]
    constexpr double DEFAULT_EXT      = 12.0;  // wide window extension [ns]
    constexpr double BIN_WIDTH        =  0.2;  // bin width [ns/bin]
    constexpr double SB_GAP           =  1.0;  // guard gap between signal and sidebands [ns]
    constexpr int    MIN_EVENTS_FIT   =   10;
}

// Result struct filled by fitSplusB()
struct SBResult {
    double mu        = 0.0;
    double muErr     = 0.0;
    double sigma     = 0.0;
    double sigmaErr  = 0.0;
    double N_sig     = 0.0;
    double N_bkg     = 0.0;
    double chi2ndf   = -1.0;
    bool   fit_ok    = false;

    double bkg_per_bin     = 0.0;
    double bkg_err_per_bin = 0.0;

    double mu_sub       = 0.0;
    double sigma_sub    = 0.0;
    double sigmaErr_sub = 0.0;
    bool   sub_ok       = false;

    // h_wide: raw histogram on the wide window (for pad1 plot)
    // h_sub:  subtracted histogram on the signal window (for pad2 plot)
    TH1D*  h_wide = nullptr;
    TH1D*  h_sub  = nullptr;

    std::string sigma_source;
};

static TH1D* sbBuildHisto(const std::vector<TOTEvent>& events,
                           double lo, double hi,
                           const std::string& name,
                           const std::string& title = "")
{
    int nBins = std::max(5, (int)std::round((hi - lo) / SB::BIN_WIDTH));
    TH1D* h = new TH1D(name.c_str(),
                        (title.empty() ? name : title).c_str(),
                        nBins, lo, hi);
    h->SetDirectory(nullptr);
    h->Sumw2();
    for (auto& e : events)
        if (e.delta_t >= lo && e.delta_t < hi)
            h->Fill(e.delta_t);
    return h;
}

// Internal: build background-subtracted histogram from sidebands.
static TH1D* sidebandSubtract(const std::vector<TOTEvent>& events,
                               double fit_lo, double fit_hi,
                               double sb_width,
                               const std::string& tag,
                               double& bkg_per_bin_out,
                               double& bkg_err_per_bin_out)
{
    const double gap = SB::SB_GAP;  // guard gap: keeps signal tails out of sidebands
    TH1D* hSL = sbBuildHisto(events, fit_lo - gap - sb_width, fit_lo - gap, "hSBL_" + tag);
    TH1D* hSR = sbBuildHisto(events, fit_hi + gap, fit_hi + gap + sb_width, "hSBR_" + tag);

    double nL = hSL->Integral(), nR = hSR->Integral();
    int    bL = hSL->GetNbinsX(), bR = hSR->GetNbinsX();

    double bpb = 0.0, berr = 0.0;
    if (nL > 0 && nR > 0) {
        bpb  = 0.5 * (nL/bL + nR/bR);
        berr = 0.5 * std::sqrt(nL/(double)(bL*bL) + nR/(double)(bR*bR));
    } else if (nL > 0) {
        bpb  = nL/bL; berr = std::sqrt(nL)/bL;
    } else if (nR > 0) {
        bpb  = nR/bR; berr = std::sqrt(nR)/bR;
    } else {
        std::cout << "  [SB] WARNING: empty sidebands for " << tag << "\n";
        bkg_per_bin_out = bkg_err_per_bin_out = 0.0;
        delete hSL; delete hSR;
        return nullptr;
    }

    bkg_per_bin_out     = bpb;
    bkg_err_per_bin_out = berr;
    std::cout << "  [SB] Background estimate: " << std::fixed << std::setprecision(3)
              << bpb << " +/- " << berr << " ev/bin"
              << "  (L=" << (int)nL << " R=" << (int)nR << ")\n";

    TH1D* hSig = sbBuildHisto(events, fit_lo, fit_hi, "hSigRaw_" + tag);
    TH1D* hSub = (TH1D*)hSig->Clone(("hSub_" + tag).c_str());
    hSub->SetDirectory(nullptr);
    hSub->SetTitle(Form("Sideband subtracted -- %s;#Deltat (ns);Events - Bkg",
                         tag.c_str()));
    for (int b = 1; b <= hSub->GetNbinsX(); ++b) {
        double raw = hSub->GetBinContent(b);
        hSub->SetBinContent(b, raw - bpb);
        hSub->SetBinError  (b, std::sqrt(raw + bpb + berr*berr));
    }
    delete hSL; delete hSR; delete hSig;
    return hSub;
}

// Main: S+B fit on Δt distribution. Returns SBResult.
static SBResult fitSplusB(const std::vector<TOTEvent>& events,
                           double fit_lo, double fit_hi,
                           double sb_width = SB::DEFAULT_SB_WIDTH,
                           double ext      = SB::DEFAULT_EXT,
                           const std::string& tag = "")
{
    SBResult res;
    if (events.empty()) return res;

    double wide_lo = fit_lo - ext;
    double wide_hi = fit_hi + ext;

    // h_wide: histogram on the wide window, used for the S+B fit and for the plot
    res.h_wide = sbBuildHisto(events, wide_lo, wide_hi, "hWide_" + tag,
        Form("S+B fit -- %s;#Deltat (ns);Events", tag.c_str()));

    if (res.h_wide->Integral() < SB::MIN_EVENTS_FIT) {
        std::cout << "  [SB] Too few events (N="
                  << (int)res.h_wide->Integral() << ") for " << tag << "\n";
        return res;
    }

    // ── 1. Peak finder in the signal window ─────────────────
    int bLo = res.h_wide->FindBin(fit_lo + 1e-6);
    int bHi = res.h_wide->FindBin(fit_hi - 1e-6);
    int bPk = bLo;
    double hMax = 0.0;
    for (int b = bLo; b <= bHi; ++b)
        if (res.h_wide->GetBinContent(b) > hMax) {
            hMax = res.h_wide->GetBinContent(b); bPk = b;
        }
    double mu0 = res.h_wide->GetBinCenter(bPk);

    // Initial sigma: weighted local mean within +/-(fit_hi-fit_lo)/4
    double wHalf = (fit_hi - fit_lo) / 4.0;
    double sw = 0, swx = 0, swx2 = 0;
    for (int b = res.h_wide->FindBin(mu0 - wHalf);
             b <= res.h_wide->FindBin(mu0 + wHalf); ++b) {
        double c = res.h_wide->GetBinContent(b), x = res.h_wide->GetBinCenter(b);
        sw += c; swx += c*x; swx2 += c*x*x;
    }
    double sigma0 = 2.0;
    if (sw > 0) {
        double m = swx/sw, v = swx2/sw - m*m;
        if (v > 0) sigma0 = std::max(0.3, std::min(std::sqrt(v), (fit_hi-fit_lo)/2.0));
        mu0 = m;
    }

    // ── 2. Estimate C0 from the MEAN of bins in the wide window tails ──
    // KEY FIX: constraining C to the level measured outside the peak
    // prevents the Gaussian from widening to absorb the flat background.
    double C0 = 0.0;
    {
        double sum = 0.0; int cnt = 0;
        for (int b = 1; b < bLo; ++b) {
            sum += res.h_wide->GetBinContent(b); ++cnt;
        }
        for (int b = bHi+1; b <= res.h_wide->GetNbinsX(); ++b) {
            sum += res.h_wide->GetBinContent(b); ++cnt;
        }
        if (cnt > 0) C0 = std::max(0.0, sum / cnt);
    }
    double A0 = std::max(1.0, hMax - C0);

    // ── 3. S+B fit on wide window ────────────────────────
    TF1* fSB = new TF1(("fSB_" + tag).c_str(),
                        "gaus(0) + pol0(3)", wide_lo, wide_hi);
    fSB->SetParameters(A0, mu0, sigma0, C0);
    fSB->SetParNames("Amplitude", "Mean", "Sigma", "Background");
    fSB->SetParLimits(0, 0.0,  1e8);
    fSB->SetParLimits(2, 0.1,  fit_hi - fit_lo);

    // Constrain C to [0.5*C0, 2.0*C0] if C0 is significant (> 0.5 ev/bin).
    // This is the main fix: C cannot collapse to zero.
    if (C0 > 0.5)
        fSB->SetParLimits(3, 0.5*C0, 2.0*C0);
    else
        fSB->SetParLimits(3, 0.0, 1e8);

    // Likelihood for low statistics, chi2 fallback
    int status = res.h_wide->Fit(fSB, "RQLN");
    if (status != 0 && status != 4000)
        status = res.h_wide->Fit(fSB, "RQWN");

    res.fit_ok = (status == 0 || status == 4000);

    if (res.fit_ok) {
        res.mu       = fSB->GetParameter(1);
        res.muErr    = fSB->GetParError(1);
        res.sigma    = std::abs(fSB->GetParameter(2));
        res.sigmaErr = fSB->GetParError(2);
        double A_fit = fSB->GetParameter(0);
        double C_fit = fSB->GetParameter(3);
        double binW  = res.h_wide->GetBinWidth(1);

        double sq2   = std::sqrt(2.0);
        double normL = TMath::Erf((fit_lo - res.mu)/(res.sigma*sq2));
        double normR = TMath::Erf((fit_hi - res.mu)/(res.sigma*sq2));
        // N_sig = integral of gaus(A,mu,sigma) from fit_lo to fit_hi
        // = A * sigma * sqrt(2π) * [Φ(hi)-Φ(lo)]  (counts, not density — no /binW)
        res.N_sig = A_fit * res.sigma * std::sqrt(2.0*TMath::Pi())
                    * 0.5*(normR - normL);
        res.N_bkg = C_fit * (int)std::round((fit_hi - fit_lo)/binW);

        double ndf  = fSB->GetNDF();
        res.chi2ndf = (ndf > 0) ? fSB->GetChisquare()/ndf : -1.0;
        double SNR  = (res.N_bkg > 0.1) ? res.N_sig/res.N_bkg : 99.0;

        std::cout << "\n  [SB] S+B fit result -- " << tag << ":\n"
                  << std::fixed << std::setprecision(3)
                  << "    μ     = " << res.mu    << " ± " << res.muErr    << " ns\n"
                  << "    σ     = " << res.sigma << " ± " << res.sigmaErr << " ns\n"
                  << "    N_sig = " << std::setprecision(0) << res.N_sig
                  << "   N_bkg = " << res.N_bkg
                  << "   S/B = "   << std::setprecision(2) << SNR << "\n"
                  << "    χ²/ndf = " << std::setprecision(2) << res.chi2ndf << "\n";

        if (res.sigma < 0.1 || res.sigma > (fit_hi - fit_lo)) {
            std::cout << "  [SB] WARNING: sigma out of physical range -- S+B fit discarded\n";
            res.fit_ok = false;
        }
        if (std::abs(res.mu - mu0) > (fit_hi - fit_lo)) {
            std::cout << "  [SB] WARNING: mu drifted too far from peak -- S+B fit discarded\n";
            res.fit_ok = false;
        }
    } else {
        std::cout << "  [SB] WARNING: S+B fit did not converge for " << tag << "\n";
    }

    // ── 4. Sideband subtraction ────────────────────────────
    res.h_sub = sidebandSubtract(events, fit_lo, fit_hi, sb_width, tag,
                                  res.bkg_per_bin, res.bkg_err_per_bin);

    // ── 5. Gaussian fit on the subtracted histogram ─────────
    if (res.h_sub && res.h_sub->Integral() > SB::MIN_EVENTS_FIT) {
        // Starting point: use S+B fit result if available, else initial estimates
        double mu_s    = res.fit_ok ? res.mu    : mu0;
        double sigma_s = res.fit_ok ? res.sigma : sigma0;
        // Try first with a narrow window +/-2 sigma around the centre
        double g_lo = std::max(fit_lo, mu_s - 2.0*sigma_s);
        double g_hi = std::min(fit_hi, mu_s + 2.0*sigma_s);

        TF1* fG = new TF1(("fGSub_" + tag).c_str(), "gaus", g_lo, g_hi);
        fG->SetParameters(res.h_sub->GetMaximum(), mu_s, sigma_s);
        fG->SetParLimits(2, 0.1, fit_hi - fit_lo);

        int sFit = res.h_sub->Fit(fG, "RQWN");
        res.sub_ok = (sFit == 0 || sFit == 4000);
        if (res.sub_ok) {
            res.sigma_sub    = std::abs(fG->GetParameter(2));
            res.sigmaErr_sub = fG->GetParError(2);
            res.mu_sub       = fG->GetParameter(1);
            if (res.sigma_sub < 0.1 || res.sigma_sub > (fit_hi - fit_lo))
                res.sub_ok = false;
            if (res.sub_ok)
                std::cout << "  [SB] Subtracted fit: sigma = "
                          << std::fixed << std::setprecision(3)
                          << res.sigma_sub << " ± " << res.sigmaErr_sub << " ns\n";
        }
        res.h_sub->GetListOfFunctions()->Add(
            (TF1*)fG->Clone(("fGSubDraw_" + tag).c_str()));
        delete fG;
    }

        // Attach cloned S+B fit to h_wide for the plot
        if (res.fit_ok)
        res.h_wide->GetListOfFunctions()->Add(
            (TF1*)fSB->Clone(("fSBDraw_" + tag).c_str()));

    delete fSB;
    return res;
}

// Draw S+B result canvas (2 pads: wide window + subtracted).
static void drawSidebandResult(SBResult& res,
                                const std::string& tag,
                                double fit_lo, double fit_hi,
                                double sb_width,
                                OutCtx& ctx)
{
    if (!res.h_wide) {
        std::cout << "  [SB] drawSidebandResult: no histogram for " << tag << "\n";
        return;
    }

    TCanvas* cSB = new TCanvas(
        ("cSB_" + tag).c_str(),
        ("Sideband Analysis -- " + tag).c_str(),
        1000, 750);
    cSB->SetFillColor(0);
    // Prevent ROOT from auto-deleting this canvas before savePNG completes
    cSB->SetBit(TObject::kCanDelete, false);

    TPad* pad1 = new TPad(("pSB1_" + tag).c_str(), "", 0.0, 0.32, 1.0, 1.0);
    TPad* pad2 = new TPad(("pSB2_" + tag).c_str(), "", 0.0, 0.00, 1.0, 0.34);
    pad1->SetLeftMargin(PAD_LEFT);  pad1->SetRightMargin(PAD_RIGHT);
    pad1->SetTopMargin(PAD_TOP);    pad1->SetBottomMargin(0.015);
    pad2->SetLeftMargin(PAD_LEFT);  pad2->SetRightMargin(PAD_RIGHT);
    pad2->SetTopMargin(0.025);      pad2->SetBottomMargin(0.32);
    pad1->SetGrid(); pad2->SetGrid();
    pad1->SetFillColor(0); pad2->SetFillColor(0);
    cSB->cd(); pad1->Draw(); pad2->Draw();

    // ── PAD 1: h_wide + S+B fit + sideband boxes ─────────────
    pad1->cd();

    // X range is defined by h_wide which is built on [wide_lo, wide_hi]
    double wide_lo = fit_lo - SB::DEFAULT_EXT;
    double wide_hi = fit_hi + SB::DEFAULT_EXT;
    double yMax = res.h_wide->GetMaximum() * 1.30;
    if (yMax <= 0) yMax = 1.0;

    // Draw order: histogram first (establishes axis frame), then colored boxes,
    // then redraw histogram on top. Boxes before histogram fills the whole pad.
    res.h_wide->SetMaximum(yMax);
    res.h_wide->SetLineColor(kBlack); res.h_wide->SetLineWidth(2);
    res.h_wide->GetXaxis()->SetLabelSize(0.0);
    res.h_wide->GetXaxis()->SetTitleSize(0.0);
    res.h_wide->GetYaxis()->SetTitle("Events");
    res.h_wide->GetYaxis()->SetTitleSize(0.065);
    res.h_wide->GetYaxis()->SetTitleOffset(0.90);
    res.h_wide->GetYaxis()->SetLabelSize(0.060);
    res.h_wide->Draw("HIST");

    double sbL_lo = std::max(wide_lo, fit_lo - sb_width);
    double sbR_hi = std::min(wide_hi, fit_hi + sb_width);
    TBox* bxL = new TBox(sbL_lo, 0.0, fit_lo, yMax);
    TBox* bxR = new TBox(fit_hi, 0.0, sbR_hi, yMax);
    TBox* bxS = new TBox(fit_lo, 0.0, fit_hi, yMax);
    bxL->SetFillColorAlpha(kOrange,  0.35); bxL->SetLineWidth(0);
    bxR->SetFillColorAlpha(kOrange,  0.35); bxR->SetLineWidth(0);
    bxS->SetFillColorAlpha(kAzure-9, 0.20); bxS->SetLineWidth(0);
    bxS->Draw(); bxL->Draw(); bxR->Draw();
    res.h_wide->Draw("HIST same");
    res.h_wide->Draw("E same");

    // Horizontal line at estimated background level
    if (res.bkg_per_bin > 0) {
        TLine* lB = new TLine(wide_lo, res.bkg_per_bin, wide_hi, res.bkg_per_bin);
        lB->SetLineColor(kRed+1); lB->SetLineStyle(2); lB->SetLineWidth(2);
        lB->Draw("same");
    }

    // Fit S+B
    TF1* fDraw = res.fit_ok
        ? (TF1*)res.h_wide->GetListOfFunctions()->FindObject(("fSBDraw_" + tag).c_str())
        : nullptr;
    if (fDraw) {
        fDraw->SetLineColor(kRed); fDraw->SetLineWidth(2);
        fDraw->SetRange(wide_lo, wide_hi);
        fDraw->Draw("same");
    }

    // Info box
    {
        TPaveText* pt = new TPaveText(0.55, 0.56, 0.94, 0.88, "NDC");
        pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
        pt->SetTextFont(42);  pt->SetTextSize(0.052);
        pt->AddText(("S+B fit  " + tag).c_str());
        if (res.fit_ok) {
            pt->AddText(Form("#mu = %.2f #pm %.2f ns",   res.mu,    res.muErr));
            pt->AddText(Form("#sigma = %.2f #pm %.2f ns",res.sigma, res.sigmaErr));
            double SNR = (res.N_bkg > 0.1) ? res.N_sig/res.N_bkg : 99.0;
            pt->AddText(Form("S=%d  B=%d  S/B=%.1f",
                             (int)res.N_sig, (int)res.N_bkg, SNR));
            pt->AddText(Form("#chi^{2}/ndf = %.2f", res.chi2ndf));
        } else {
            pt->AddText("S+B fit not converged");
            if (res.bkg_per_bin > 0)
                pt->AddText(Form("Bkg: %.2f #pm %.2f ev/bin",
                                  res.bkg_per_bin, res.bkg_err_per_bin));
        }
        pt->Draw();
    }
    {
        TLegend* leg = new TLegend(0.14, 0.74, 0.48, 0.88);
        leg->SetBorderSize(1); leg->SetFillColor(0); leg->SetTextSize(0.046);
        leg->AddEntry(res.h_wide, "Raw #Deltat (wide window)", "l");
        if (fDraw) leg->AddEntry(fDraw, "Gauss + flat bkg", "l");
        leg->AddEntry(bxL, "Sideband", "f");
        leg->AddEntry(bxS, "Signal window", "f");
        leg->Draw();
    }

    // ── PAD 2: subtracted histogram + Gaussian fit ───────────
    pad2->cd();

    if (res.h_sub) {
        res.h_sub->SetLineColor(kAzure+1); res.h_sub->SetLineWidth(2);
        res.h_sub->SetMarkerStyle(20); res.h_sub->SetMarkerSize(0.5);
        res.h_sub->SetMarkerColor(kAzure+1);
        res.h_sub->GetXaxis()->SetTitle("#Deltat (ns)");
        res.h_sub->GetXaxis()->SetTitleSize(0.11);
        res.h_sub->GetXaxis()->SetLabelSize(0.10);
        res.h_sub->GetYaxis()->SetTitle("Subtracted");
        res.h_sub->GetYaxis()->SetTitleSize(0.10);
        res.h_sub->GetYaxis()->SetLabelSize(0.09);
        res.h_sub->GetYaxis()->SetTitleOffset(0.50);
        res.h_sub->Draw("E");

        TLine* l0 = new TLine(fit_lo, 0.0, fit_hi, 0.0);
        l0->SetLineColor(kGray+1); l0->SetLineStyle(2); l0->Draw("same");

        TF1* fGD = (TF1*)res.h_sub->GetListOfFunctions()->FindObject(
            ("fGSubDraw_" + tag).c_str());
        if (fGD && res.sub_ok) {
            fGD->SetLineColor(kRed+1); fGD->SetLineWidth(2);
            fGD->Draw("same");
            TPaveText* pt2 = new TPaveText(0.55, 0.55, 0.94, 0.92, "NDC");
            pt2->SetBorderSize(1); pt2->SetFillColor(0); pt2->SetFillStyle(1001);
            pt2->SetTextFont(42); pt2->SetTextSize(0.090);
            pt2->AddText("Fit subt.");
            pt2->AddText(Form("#sigma_{sub} = %.2f #pm %.2f ns",
                               res.sigma_sub, res.sigmaErr_sub));
            pt2->Draw();
        }
    } else {
        TPaveText* pt = new TPaveText(0.15, 0.25, 0.85, 0.75, "NDC");
        pt->SetBorderSize(0); pt->SetFillColor(0);
        pt->SetTextFont(42); pt->SetTextSize(0.12);
        pt->AddText("Sidebands not available");
        pt->Draw();
    }

    cSB->cd();
    cSB->Update(); cSB->Modified();
    ctx.savePNG(cSB, Form("sideband_%s.png", tag.c_str()));
}

static void cleanupSBResult(SBResult& res)
{
    delete res.h_wide; res.h_wide = nullptr;
    delete res.h_sub;  res.h_sub  = nullptr;
}
