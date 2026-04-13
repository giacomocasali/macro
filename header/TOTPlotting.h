#pragma once
// TOTPlotting.h
// Canvas functions for the TOT analysis:
//   fillTH2D()             — fill TOT vs Δt map from event vector
//   drawTOTMap()           — 2D COLZ canvas (log Z)
//   drawGlobalProjection() — Δt projection with Gaussian + q-Gaussian fits
//   drawSlicesAndTrend()   — TOT slices with q-Gaussian fit, σ and mean trend
//   drawByPE()             — per-p.e. Δt distributions with Gaussian fits
//   drawDeltaTHistogram()  — 1D Δt histogram with counts

#include "TOTAnalysis.h"
#include "Calibration.h"
#include "Config.h"
#include "OutputManager.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>

// Asymmetric q-Gaussian for TOT slice fits.
// par: [0]=A [1]=mean [2]=sigma [3]=q1 [4]=q2
// q1 controls the left tail, q2 the right tail. For q→1 reduces to a Gaussian.
// In fitSliceTOT() q1 is fixed to 1 (symmetric left side).
static Double_t qGaussAsymTOT(Double_t* x, Double_t* par) {
    double xx=x[0], A=par[0], mean=par[1], sigma=par[2], q1=par[3], q2=par[4];
    if (sigma<=0) return 0;
    if (xx<=mean) {
        double arg=1.0-(1.0-q1)*(1.0/(3.0-q1))*std::pow((xx-mean)/sigma,2);
        return (arg<=0)?0:A*std::pow(arg,1.0/(1.0-q1));
    } else {
        double arg=1.0-(1.0-q2)*(1.0/(3.0-q2))*std::pow((xx-mean)/sigma,2);
        return (arg<=0)?0:A*std::pow(arg,1.0/(1.0-q2));
    }
}

static std::tuple<double,double,double,double,TF1*>
fitSliceTOT(TH1D* h, int idx, const std::string& tag,
            double fit_lo, double fit_hi) {
    if (!h) return {0,0,0,0,nullptr};
    int b1 = h->FindBin(fit_lo), b2 = h->FindBin(fit_hi);
    double entries_win = h->Integral(b1, b2);
    if (entries_win < 10) return {0,0,0,0,nullptr};

    int maxBin = b1;
    for (int b = b1; b <= b2; ++b)
        if (h->GetBinContent(b) > h->GetBinContent(maxBin)) maxBin = b;
    double center = h->GetBinCenter(maxBin);

    // Sigma estimate: weighted mean in ±5 ns around peak.
    // Narrow window avoids flat tails inflating sigEst and pulling the fit off-peak.
    double sigEst = 0.5;
    {
        double sw=0, swx=0, swx2=0;
        int bl = h->FindBin(center - 5.0), br = h->FindBin(center + 5.0);
        bl = std::max(bl, b1); br = std::min(br, b2);
        for (int b=bl; b<=br; ++b) {
            double c=h->GetBinContent(b), x=h->GetBinCenter(b);
            sw+=c; swx+=c*x; swx2+=c*x*x;
        }
        if (sw>0) {
            double m=swx/sw;
            double var=swx2/sw - m*m;
            sigEst = std::max(std::sqrt(std::max(var,0.0)), 0.05);
            center = m;  // refine center with weighted mean
        }
    }

    std::string fname = Form("fSliceTOT_%s_%d", tag.c_str(), idx);
    TF1* f = new TF1(fname.c_str(), qGaussAsymTOT,
                     std::max(fit_lo, center-3.0*sigEst),
                     std::min(fit_hi, center+3.0*sigEst), 5);
    f->SetParameters(h->GetBinContent(maxBin), center, sigEst, 1.0, 1.2);
    f->FixParameter(3, 1.0);
    h->Fit(f, "RMQN");
    return {f->GetParameter(1), f->GetParError(1),
            std::abs(f->GetParameter(2)), f->GetParError(2), f};
}

// ════════════════════════════════════════════════════════════
//  RIEMPIMENTO TH2D dagli eventi raccolti
// ════════════════════════════════════════════════════════════
static TH2D* fillTH2D(const std::vector<TOTEvent>& events,
                       const std::string& tag,
                       double frac_pe,
                       bool corrected = false) {
    std::string suffix = corrected ? "_corr" : "";
    std::string label  = corrected ? " (corrected)" : "";
    // X: TOT 0-150 ns (300 bins = 0.5 ns/bin)
    // Y: Δt  -50..200 ns (1000 bins = 0.25 ns/bin)
    TH2D* h = new TH2D(
        Form("h2D%s_%s", suffix.c_str(), tag.c_str()),
        Form("TOT vs #Delta t%s   LET=%.2f p.e.   %s"
             ";TOT (ns);#Delta t (ns)",
             label.c_str(), frac_pe, tag.c_str()),
        300, 0.0, 150.0, 1000, -50.0, 200.0);
    h->SetDirectory(nullptr);
    for (auto& e : events)
        h->Fill(e.tot, e.delta_t);
    return h;
}

// TOT vs Δt 2D map — log Z colour scale
static void drawTOTMap(TH2D* h2D, const std::string& tag,
                       OutCtx& ctx, bool corrected = false) {
    std::string suffix = corrected ? "_corr" : "";
    TCanvas* c = new TCanvas(Form("c2D%s_%s", suffix.c_str(), tag.c_str()),
        Form("TOT map%s -- %s", corrected?" (corr)":"", tag.c_str()), 900, 700);
    c->SetRightMargin(0.15); c->SetLeftMargin(PAD_LEFT);
    c->SetBottomMargin(PAD_BOTTOM); c->SetTopMargin(PAD_TOP);
    c->SetLogz(); c->SetGrid();
    h2D->Draw("COLZ");
    c->Update(); c->Modified();
    ctx.savePNG(c, Form("tot_map%s_%s.png", suffix.c_str(), tag.c_str()));
    // delete c;  // FIX: canvas resta aperta per gApplication->Run()
}

// q-Gaussian for global projection fits. Left side (q1): exact Gaussian (q1→1).
// Right side (q2): heavier tail to describe afterpulse/residual time-walk.
// par: [0]=A [1]=mean [2]=sigma [3]=q1 [4]=q2
static Double_t qGaussProj(Double_t* x, Double_t* par) {
    double xx=x[0], A=par[0], mean=par[1], sigma=par[2];
    double q1=par[3], q2=par[4];
    if (sigma<=0) return 0;
    const double eps=1e-6;
    if (xx<=mean) {
        if (std::abs(q1-1.0)<eps) {
            return A*std::exp(-0.5*std::pow((xx-mean)/sigma,2));
        }
        double arg=1.0-(1.0-q1)*(1.0/(3.0-q1))*std::pow((xx-mean)/sigma,2);
        return (arg<=0)?0:A*std::pow(arg,1.0/(1.0-q1));
    } else {
        double arg=1.0-(1.0-q2)*(1.0/(3.0-q2))*std::pow((xx-mean)/sigma,2);
        return (arg<=0)?0:A*std::pow(arg,1.0/(1.0-q2));
    }
}

// Δt projection with two fits:
//   Left panel  : Gaussian over ±3σ around peak
//   Right panel : q-Gaussian over full [fit_lo, fit_hi] — describes right tail
// If eventsPtr is provided, fills the histogram directly from delta_t values
// (0.02 ns/bin resolution). Otherwise falls back to ProjectionY of the TH2D.
static void drawGlobalProjection(TH2D* h2D,
                                  const std::string& tag,
                                  double fit_lo, double fit_hi,
                                  OutCtx& ctx,
                                  bool corrected = false,
                                  const std::vector<TOTEvent>* eventsPtr = nullptr) {
    std::string suffix = corrected ? "_corr" : "";
    // Proiezione su Y con binning fine (target 0.02 ns/bin).
    //
    // Bug fix (bug 15): il vecchio metodo riempiva questo istogramma usando
    // i centri bin del TH2D (risoluzione 0.25 ns/bin), dando l'illusione di
    // avere 0.02 ns/bin ma con risoluzione effettiva limitata a 0.25 ns.
    // Il fit gaussiano su quell'istogramma non aveva quindi migliore precisione
    // di ProjectionY(), contraddicendo i commenti nel codice.
    //
    // Fix: se il vettore eventi grezzi è disponibile (eventsPtr != nullptr),
    // riempi l'istogramma direttamente dai delta_t originali — risoluzione
    // reale 0.02 ns/bin.  Altrimenti usa il fallback dal TH2D (retrocompatibile
    // con chiamate che non passano il vettore).
    int nBinsProj = std::max(400, (int)std::round((fit_hi - fit_lo) / 0.02));
    TH1D* hProj = new TH1D(
        Form("hProj%s_%s", suffix.c_str(), tag.c_str()),
        Form("Global #Deltat projection%s   %s;#Deltat (ns);Events",
             corrected?" (corrected)":"", tag.c_str()),
        nBinsProj, fit_lo - 1.0, fit_hi + 1.0);
    hProj->SetDirectory(nullptr);

    if (eventsPtr && !eventsPtr->empty()) {
        // Percorso corretto: riempimento diretto dai delta_t grezzi → 0.02 ns/bin reale
        for (auto& e : *eventsPtr)
            hProj->Fill(e.delta_t);
    } else {
        // Fallback quando il vettore eventi non è disponibile:
        // usa ProjectionY() del TH2D (risoluzione 0.25 ns/bin).
        // NON si usa più il riempimento via GetBinCenter(by) che fingeva 0.02 ns/bin
        // ma aveva risoluzione reale 0.25 ns (bug 15).
        TH1D* hPY = h2D->ProjectionY(
            Form("hProjYtmp%s_%s", suffix.c_str(), tag.c_str()),
            1, h2D->GetNbinsX());
        hPY->SetDirectory(nullptr);
        // Ricopia nel hProj (con binning adattivo) interpolando i bin di hPY
        for (int by = 1; by <= hPY->GetNbinsX(); ++by) {
            double val = hPY->GetBinContent(by);
            if (val > 0) {
                double xc = hPY->GetBinCenter(by);
                // Distribuisce uniformemente il contenuto del bin
                // su N sottobins di hProj per evitare picchi a dente di sega.
                // N = ceil(hPY->GetBinWidth / hProj->GetBinWidth)
                double wPY   = hPY->GetBinWidth(by);
                double wProj = hProj->GetBinWidth(1);
                int    nSub  = std::max(1, (int)std::ceil(wPY / wProj));
                double step  = wPY / nSub;
                double frac  = val / nSub;
                for (int s = 0; s < nSub; ++s)
                    hProj->Fill(xc - wPY/2.0 + (s + 0.5)*step, frac);
            }
        }
        delete hPY;
    }

    // ── Trova il picco reale: bin di massimo in [fit_lo, fit_hi] ──
    int b1 = hProj->FindBin(fit_lo), b2 = hProj->FindBin(fit_hi);
    int maxBin = b1;
    for (int b=b1; b<=b2; ++b)
        if (hProj->GetBinContent(b) > hProj->GetBinContent(maxBin)) maxBin=b;
    double center = hProj->GetBinCenter(maxBin);
    double A0     = hProj->GetBinContent(maxBin);

    // Stima σ con media pesata in una finestra STRETTA ±10 ns intorno al picco.
    // Usare tutta la finestra [fit_lo,fit_hi] con code piatte darebbe roughSig
    // enormemente sovrastimato, mandando il fit fuori dal picco.
    double wHalf = 10.0;
    double sumW=0, sumWX=0, sumWX2=0;
    int bw1 = hProj->FindBin(center - wHalf);
    int bw2 = hProj->FindBin(center + wHalf);
    bw1 = std::max(bw1, b1); bw2 = std::min(bw2, b2);
    for (int b=bw1; b<=bw2; ++b) {
        double cnt=hProj->GetBinContent(b), xc=hProj->GetBinCenter(b);
        sumW+=cnt; sumWX+=cnt*xc; sumWX2+=cnt*xc*xc;
    }
    double roughSig = 1.0;  // default conservativo: 1 ns
    if (sumW>0) {
        double m=sumWX/sumW;
        double v=sumWX2/sumW - m*m;
        roughSig = std::max(std::sqrt(std::max(v,0.0)), 0.05);
    }
    // Update center with local weighted mean (more precise than bin max alone)
    if (sumW>0) center = sumWX/sumW;

    // ── FIT 1: Gaussiana su ±3σ intorno al picco ───────────
    // La finestra si adatta automaticamente alla larghezza reale del picco.
    double fr_lo = std::max(fit_lo, center - 3.0*roughSig);
    double fr_hi = std::min(fit_hi, center + 3.0*roughSig);
    TF1* fG = new TF1(Form("fGauss%s_%s", suffix.c_str(), tag.c_str()),
                       "gaus", fr_lo, fr_hi);
    fG->SetParameters(A0, center, roughSig);
    hProj->Fit(fG, "RQN");
    // Se il fit gaussiano ha spostato mu fuori dalla finestra ±5σ, riprova
    // with initial parameters closer to the peak (prevents fit divergence)
    if (std::abs(fG->GetParameter(1) - center) > 5.0*roughSig) {
        fG->SetParameters(A0, center, roughSig);
        fG->SetParLimits(1, center - 5.0*roughSig, center + 5.0*roughSig);
        hProj->Fit(fG, "RQNB");
    }

    double mu_g    = fG->GetParameter(1);
    double sig_g   = std::abs(fG->GetParameter(2));
    double muErr_g = fG->GetParError(1);
    double sgErr_g = fG->GetParError(2);
    double chi2_g  = fG->GetChisquare();
    int    ndf_g   = fG->GetNDF();

    // ── FIT 2: Q-gaussiana asimmetrica su tutta la finestra ─
    // q1 fisso a 1 (coda sinistra gaussiana, fisicamente motivato:
    // non ci sono motivi per avere una coda sinistra pesante).
    // q2 libero: descrive la coda destra da afterpulse/time walk residuo.
    TF1* fQ = new TF1(Form("fQGauss%s_%s", suffix.c_str(), tag.c_str()),
                       qGaussProj, fit_lo, fit_hi, 5);
    fQ->SetParNames("A","#mu","#sigma","q_{1}","q_{2}");
    fQ->SetParameters(A0, mu_g, sig_g, 1.0, 1.3);
    fQ->FixParameter(3, 1.0);           // q1=1: coda sinistra gaussiana
    fQ->SetParLimits(4, 1.0, 2.0);      // q2 in [1,2]: heavier right tail
    fQ->SetNpx(1000);
    hProj->Fit(fQ, "RQMN");

    double mu_q    = fQ->GetParameter(1);
    double sig_q   = std::abs(fQ->GetParameter(2));
    double muErr_q = fQ->GetParError(1);
    double sgErr_q = fQ->GetParError(2);
    double q2_val  = fQ->GetParameter(4);
    double chi2_q  = fQ->GetChisquare();
    int    ndf_q   = fQ->GetNDF();

    // ── Stampa a terminale — solo parametri fisici ──────────
    std::cout << "  " << (corrected ? "Corrected" : "Raw") << " Δt projection:\n"
              << "    Gaussian:   μ = " << mu_g  << " ± " << muErr_g  << " ns"
              << "   σ = " << sig_g  << " ± " << sgErr_g  << " ns"
              << "   χ²/ndf = " << (ndf_g>0 ? chi2_g/ndf_g : -1.) << "\n"
              << "    q-Gaussian: μ = " << mu_q  << " ± " << muErr_q  << " ns"
              << "   σ = " << sig_q  << " ± " << sgErr_q  << " ns"
              << "   q₂ = " << q2_val
              << "   χ²/ndf = " << (ndf_q>0 ? chi2_q/ndf_q : -1.) << "\n";

    // ── Canvas con due pannelli ─────────────────────────────
    // Zoom window: centre on the fitted peak, show ±max(8*sigma, 3 ns)
    // so the peak always fills most of the canvas regardless of fit_lo/fit_hi.
    double zoom_half = std::max(8.0 * sig_g, 3.0);
    double zoom_lo   = std::max(fit_lo, mu_g - zoom_half);
    double zoom_hi   = std::min(fit_hi, mu_g + zoom_half);

    TCanvas* cP = new TCanvas(
        Form("cProj%s_%s", suffix.c_str(), tag.c_str()),
        Form("Global projection%s — %s", corrected?" (corr)":"", tag.c_str()),
        1400, 600);
    cP->Divide(2, 1);

    // Left panel: Gaussian
    cP->cd(1);
    gPad->SetGrid(); gPad->SetLeftMargin(PAD_LEFT);
    gPad->SetBottomMargin(PAD_BOTTOM); gPad->SetTopMargin(PAD_TOP);
    TH1D* hL = (TH1D*)hProj->Clone(Form("hProjL%s_%s",suffix.c_str(),tag.c_str()));
    hL->SetDirectory(nullptr);
    hL->SetLineColor(kAzure+1); hL->SetLineWidth(2);
    hL->GetXaxis()->SetRangeUser(zoom_lo, zoom_hi);
    hL->Draw("HIST");
    fG->SetLineColor(kRed+1); fG->SetLineWidth(2);
    fG->Draw("same");
    {
        TPaveText* pt=new TPaveText(0.52,0.60,0.95,0.88,"NDC");
        pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
        pt->SetTextFont(42); pt->SetTextSize(0.038);
        pt->AddText("Gaussian fit");
        pt->AddText(Form("#mu = %.3f #pm %.3f ns",   mu_g, muErr_g));
        pt->AddText(Form("#sigma = %.3f #pm %.3f ns", sig_g, sgErr_g));
        pt->AddText(Form("#chi^{2}/ndf = %.0f / %d",  chi2_g, ndf_g));
        if (corrected) pt->AddText("Time walk corrected");
        pt->Draw();
        delete pt;
    }

    // Right panel: q-Gaussian
    cP->cd(2);
    gPad->SetGrid(); gPad->SetLeftMargin(PAD_LEFT);
    gPad->SetBottomMargin(PAD_BOTTOM); gPad->SetTopMargin(PAD_TOP);
    TH1D* hR = (TH1D*)hProj->Clone(Form("hProjR%s_%s",suffix.c_str(),tag.c_str()));
    hR->SetDirectory(nullptr);
    hR->SetLineColor(kAzure+1); hR->SetLineWidth(2);
    hR->GetXaxis()->SetRangeUser(zoom_lo, zoom_hi);
    hR->Draw("HIST");
    fQ->SetLineColor(kGreen+2); fQ->SetLineWidth(2);
    fQ->Draw("same");
    {
        TPaveText* pt=new TPaveText(0.52,0.55,0.95,0.88,"NDC");
        pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
        pt->SetTextFont(42); pt->SetTextSize(0.038);
        pt->AddText("q-Gaussian fit");
        pt->AddText(Form("#mu = %.3f #pm %.3f ns",   mu_q, muErr_q));
        pt->AddText(Form("#sigma = %.3f #pm %.3f ns", sig_q, sgErr_q));
        pt->AddText(Form("q_{2} = %.3f", q2_val));
        pt->AddText(Form("#chi^{2}/ndf = %.0f / %d",  chi2_q, ndf_q));
        if (corrected) pt->AddText("Time walk corrected");
        pt->Draw();
        delete pt;
    }

    cP->Update(); cP->Modified();
    ctx.savePNG(cP, Form("tot_projection%s_%s.png", suffix.c_str(), tag.c_str()));
    delete fG;
    delete fQ;
    delete hProj;
    delete hL;
    delete hR;
    // delete cP;  // FIX: canvas resta aperta
} //6 equal-width TOT bins, q-Gaussian fit per slice.
// Trend canvas: σ(Δt) vs TOT and Mean(Δt) vs TOT.
static void drawSlicesAndTrend(TH2D* h2D,
                                const std::string& tag,
                                double fit_lo, double fit_hi,
                                double frac_pe,
                                OutCtx& ctx,
                                bool corrected = false) {
    std::string suffix = corrected ? "_corr" : "";
    const int N_SLICES = 6;
    double tot_max = 0.0;
    for (int bx=h2D->GetNbinsX(); bx>=1; --bx) {
        double s=0;
        for (int by=1; by<=h2D->GetNbinsY(); ++by)
            s+=h2D->GetBinContent(bx,by);
        if (s>5) { tot_max=h2D->GetXaxis()->GetBinCenter(bx); break; }
    }
    if (tot_max<=0) tot_max=20.0;
    double slice_w = tot_max / N_SLICES;

    TCanvas* cSl = new TCanvas(
        Form("cSl%s_%s", suffix.c_str(), tag.c_str()),
        Form("TOT slices%s — %s", corrected?" (corr)":"", tag.c_str()),
        1200, 800);
    cSl->Divide(3, 2);

    TGraphErrors* grRes  = new TGraphErrors();
    TGraphErrors* grMean = new TGraphErrors();
    int pts=0;

    for (int s=0; s<N_SLICES; ++s) {
        double xlo_s=s*slice_w, xhi_s=xlo_s+slice_w;
        int bx1=h2D->GetXaxis()->FindBin(xlo_s);
        int bx2=h2D->GetXaxis()->FindBin(xhi_s)-1;
        std::string hname=Form("hsl%s_%s_%d", suffix.c_str(), tag.c_str(), s);
        TH1D* hsl=h2D->ProjectionY(hname.c_str(), bx1, bx2);
        hsl->SetDirectory(nullptr);
        hsl->SetTitle(Form("TOT #in [%.1f, %.1f] ns%s;#Deltat (ns);Counts",
                           xlo_s, xhi_s, corrected?" (corr)":""));
        hsl->GetXaxis()->UnZoom(); hsl->GetYaxis()->UnZoom();
        cSl->cd(s+1);
        auto [mean,meanErr,sigma,sigmaErr,fFit] =
            fitSliceTOT(hsl, s, tag+suffix, fit_lo, fit_hi);
        if (fFit) {
            // Quality cut: reject fits with unphysical σ or mean outside window
            bool fitOk = (sigma > 0.01 && sigma < 5.0 &&
                          mean  > fit_lo && mean < fit_hi &&
                          sigmaErr > 0 && sigmaErr < sigma);
            if (fitOk) {
                hsl->Draw("HIST");
                fFit->SetLineColor(kGreen+2); fFit->Draw("SAME");
                double tc=0.5*(xlo_s+xhi_s);
                grRes->SetPoint(pts,tc,sigma);   grRes->SetPointError(pts,0,sigmaErr);
                grMean->SetPoint(pts,tc,mean);   grMean->SetPointError(pts,0,meanErr);
                ++pts;
            }
            delete fFit;
        }
        delete hsl;
    }
    cSl->Update(); cSl->Modified();
    ctx.savePNG(cSl, Form("tot_slices%s_%s.png", suffix.c_str(), tag.c_str()));
    // delete cSl;  // FIX: canvas resta aperta

    if (pts<2) return;

    TCanvas* cTrend = new TCanvas(
        Form("cTrend%s_%s", suffix.c_str(), tag.c_str()),
        Form("TOT trend%s — %s", corrected?" (corr)":"", tag.c_str()),
        1200, 500);
    cTrend->Divide(2,1);

    cTrend->cd(1); gPad->SetGrid();
    grRes->SetMarkerStyle(20); grRes->SetMarkerColor(kAzure+1);
    grRes->SetLineColor(kAzure+1);
    grRes->SetTitle(Form("#sigma_{#Deltat} vs TOT%s   LET=%.2f p.e."
        ";TOT (ns);#sigma_{#Deltat} (ns)",
        corrected?" (corr)":"", frac_pe));
    grRes->Draw("AP");
    if (pts>=3) {
        TF1* fS=new TF1(Form("fS%s_%s",suffix.c_str(),tag.c_str()),
            "sqrt([0]*[0]/x+[1]*[1])",
            grRes->GetX()[0]*0.8, grRes->GetX()[pts-1]*1.2);
        fS->SetParameters(1.0,0.1); fS->SetLineColor(kRed+1);
        grRes->Fit(fS,"RQ");
        if (!corrected)
            std::cout << "  σ(TOT) fit:   S = " << fS->GetParameter(0)
                      << " ns·√ns   C = " << fS->GetParameter(1) << " ns\n";
        delete fS;
    }

    cTrend->cd(2); gPad->SetGrid();
    grMean->SetMarkerStyle(20); grMean->SetMarkerColor(kOrange+7);
    grMean->SetLineColor(kOrange+7);
    grMean->SetTitle(Form("Mean #Deltat vs TOT%s   LET=%.2f p.e."
        ";TOT (ns);Mean #Deltat (ns)",
        corrected?" (corr)":"", frac_pe));
    grMean->Draw("AP");

    cTrend->Update(); cTrend->Modified();
    ctx.savePNG(cTrend, Form("tot_trend%s_%s.png", suffix.c_str(), tag.c_str()));
    delete grRes;
    delete grMean;
    // delete cTrend;  // FIX: canvas resta aperta
}

// Per-p.e. Δt distributions with Gaussian fits. Also produces σ vs n_pe trend.
static void drawByPE(const std::vector<TOTEvent>& events,
                     const std::string& tag,
                     double fit_lo, double fit_hi,
                     OutCtx& ctx,
                     int max_pe = 6) {
    std::map<int, std::vector<double>> byPE;
    for (auto& e : events) {
        if (e.n_pe < 0 || e.n_pe > max_pe) continue;
        byPE[e.n_pe].push_back(e.delta_t);
    }
    if (byPE.empty()) return;

    int npe_found = (int)byPE.size();
    int ncols = std::min(npe_found, 3);
    int nrows = (npe_found + ncols - 1) / ncols;

    TCanvas* cPE = new TCanvas(Form("cPE_%s", tag.c_str()),
        Form("Δt by p.e. — %s", tag.c_str()),
        400*ncols, 350*nrows);
    cPE->Divide(ncols, nrows);

    TGraphErrors* grSig = new TGraphErrors();
    int pts=0;

    int pad=1;
    for (auto& [npe, vals] : byPE) {
        // Fine binning: 0.05 ns/bin so that a σ~0.15 ns peak has ~3 bins
        int nBinsPE = std::max(100, (int)std::round((fit_hi - fit_lo) / 0.05));
        TH1D* h = new TH1D(Form("hPE%d_%s", npe, tag.c_str()),
            Form("%d p.e.;#Deltat (ns);Events", npe),
            nBinsPE, fit_lo, fit_hi);
        h->SetDirectory(nullptr);
        for (double v : vals) h->Fill(v);

        cPE->cd(pad++);
        gPad->SetGrid();
        h->SetLineColor(kAzure+1); h->SetLineWidth(2);
        h->Draw("HIST");

        // Fit gaussiano sul picco reale (±5 ns intorno al massimo)
        if (h->GetEntries() >= 10) {
            int bmax = h->GetMaximumBin();
            double ctr = h->GetBinCenter(bmax);
            // Weighted-mean sigma in ±5 ns around peak
            double sg = 0.5;
            {
                double sw=0,swx=0,swx2=0;
                int bl=h->FindBin(ctr-5.0), br=h->FindBin(ctr+5.0);
                for(int b=bl;b<=br;++b){
                    double c=h->GetBinContent(b),x=h->GetBinCenter(b);
                    sw+=c; swx+=c*x; swx2+=c*x*x;
                }
                if(sw>0){double m=swx/sw; sg=std::max(std::sqrt(std::max(swx2/sw-m*m,0.0)),0.05);}
            }
            TF1* fg=new TF1(Form("fg%d_%s",npe,tag.c_str()),
                "gaus",
                std::max(fit_lo, ctr-3.0*sg),
                std::min(fit_hi, ctr+3.0*sg));
            fg->SetParameters(h->GetMaximum(),ctr,sg);
            fg->SetParLimits(1, ctr-5.0*sg, ctr+5.0*sg);
            h->Fit(fg,"RQB");
            fg->SetLineColor(kRed+1); fg->Draw("same");

            double sigma=std::abs(fg->GetParameter(2));
            double sgErr=fg->GetParError(2);
            if (sigma>0) {
                grSig->SetPoint(pts,npe,sigma);
                grSig->SetPointError(pts,0,sgErr);
                ++pts;
            }

            TPaveText* pt=new TPaveText(0.55,0.72,0.93,0.88,"NDC");
            pt->SetBorderSize(1); pt->SetFillColor(0);
            pt->SetTextFont(42); pt->SetTextSize(0.040);
            pt->AddText(Form("#sigma = %.2f #pm %.2f ns",sigma,sgErr));
            pt->Draw();
            delete pt;
            delete fg;
        }
        delete h;
    }
    cPE->Update(); cPE->Modified();
    ctx.savePNG(cPE, Form("tot_byPE_%s.png", tag.c_str()));
    // delete cPE;  // FIX: canvas resta aperta

    // σ(Δt) vs n_pe trend
    if (pts >= 2) {
        TCanvas* cSPE=new TCanvas(Form("cSigPE_%s",tag.c_str()),
            Form("#sigma vs p.e. — %s",tag.c_str()), 700, 500);
        cSPE->SetGrid();
        grSig->SetMarkerStyle(20); grSig->SetMarkerColor(kAzure+1);
        grSig->SetTitle(Form("#sigma_{#Deltat} vs p.e.   %s"
            ";p.e. number;#sigma_{#Deltat} (ns)", tag.c_str()));
        grSig->Draw("AP");
        if (pts>=3) {
            TF1* fS=new TF1(Form("fSPE_%s",tag.c_str()),
                "sqrt([0]*[0]/x+[1]*[1])",0.5,max_pe+0.5);
            fS->SetParameters(1.0,0.1);
            fS->SetLineColor(kRed+1);
            grSig->Fit(fS,"RQ");
            std::cout << "  σ(n_pe) fit:  S = " << fS->GetParameter(0)
                      << " ns   C = " << fS->GetParameter(1) << " ns\n";
            delete fS;
        }
        cSPE->Update(); cSPE->Modified();
        ctx.savePNG(cSPE, Form("tot_sigma_vs_pe_%s.png", tag.c_str()));
        delete grSig;
        // delete cSPE;  // FIX: canvas resta aperta
    }
}

// 1D Δt histogram with counts
static void drawDeltaTHistogram(const std::vector<TOTEvent>& events,
                                const std::string& tag,
                                double fit_lo, double fit_hi,
                                OutCtx& ctx,
                                bool corrected = false) {
    std::string suffix = corrected ? "_corr" : "";
    // Fine binning: 0.02 ns/bin resolution
    int nBins = std::max(400, (int)std::round((fit_hi - fit_lo) / 0.02));
    TH1D* h = new TH1D(
        Form("hDeltaT%s_%s", suffix.c_str(), tag.c_str()),
        Form("#Deltat histogram%s   %s;#Deltat (ns);Counts",
             corrected ? " (corrected)" : "", tag.c_str()),
        nBins, fit_lo - 1.0, fit_hi + 1.0);
    h->SetDirectory(nullptr);
    for (auto& e : events)
        h->Fill(e.delta_t);

    TCanvas* c = new TCanvas(
        Form("cDeltaT%s_%s", suffix.c_str(), tag.c_str()),
        Form("#Deltat histogram%s — %s", corrected ? " (corr)" : "", tag.c_str()),
        900, 600);
    c->SetGrid(); c->SetLeftMargin(PAD_LEFT);
    c->SetBottomMargin(PAD_BOTTOM); c->SetTopMargin(PAD_TOP);
    h->SetLineColor(kAzure+1); h->SetLineWidth(2);
    h->Draw("HIST");
    c->Update(); c->Modified();
    ctx.savePNG(c, Form("delta_t_histogram%s_%s.png", suffix.c_str(), tag.c_str()));
    delete h;
    // delete c;  // FIX: canvas resta aperta
}

// saveTOTRoot is replaced by ctx.saveRoot() in OutputManager.h
