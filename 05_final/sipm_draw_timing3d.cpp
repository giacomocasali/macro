/**
 * sipm_draw_timing3d.cpp  — v8  (FIXED)
 * =========================================
 * Post-process di sipm_pos_scan(). Legge map_results_*.root e produce:
 *   mu_3d_<tag>.png        SURF2 3D  μ(Δt) [ns]
 *   mu_2d_<tag>.png        COLZ 2D  μ(Δt) [ns]
 *   sigma_3d_<tag>.png     SURF2 3D  σ(Δt) [ns]
 *   sigma_2d_<tag>.png     COLZ 2D  σ(Δt) [ns]
 *   pdet_<tag>.png         COLZ 2D  P_det [%]
 *   delay_map_2d_<tag>.png COLZ 2D  Δμ(ps)
 *   delay_map_x_<tag>.png  TGraphErrors Δμ vs Δx
 *   delay_map_y_<tag>.png  TGraphErrors Δμ vs Δy
 *
 * FILTRI — identici a drawMap2D in sipm_pos_scan.cpp:
 *   mu_ok=1, sigma < MAD adattivo (mediana+10*MAD), N_acc >= MIN_N_ACC.
 *   Outlier (|delay-med|>3*MAD) esclusi da scala e fit delay.
 *
 * Compile:  .L sipm_draw_timing3d.cpp+
 * Run:      sipm_draw_timing3d()
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "../header/Config.h"
#include "../header/InputHelpers.h"
#include "../header/OutputManager.h"

// ── Constants ─────────────────────────────────────────────────────────────────
static constexpr double C_LIGHT_MMPS  = 0.2998;
static constexpr long   MIN_N_ACC     = 5000;   // min events to include bin
static constexpr double SIGMA_PRECUT  = 0.5;    // initial pre-cut
static constexpr double MU_MAXDEV_NS  = 0.5;    // |mu - median| > 0.5 ns = outlier

static const int VBIAS_COLORS[] = {
    kRed+1, kOrange+7, kGreen+2, kAzure+1, kViolet+1, kCyan+2, kYellow+3 };
static int vbiasColor(int idx) { return VBIAS_COLORS[idx % 7]; }

// ────────────────────────────────────────────────────────────────────────────
static double computeAdaptiveSigmaThreshold(const std::vector<struct MapPt>& pts);

// ── MapPt ─────────────────────────────────────────────────────────────────────
struct MapPt {
    double x, y;
    int    vbias;
    double mu, mu_err;
    double sigma, sigma_err;
    double n_sig, p_det;
    long   n_acc, n_crossing, n_laser;
    double chi2ndf;
    int    mu_ok;
};

// ────────────────────────────────────────────────────────────────────────────
static double computeAdaptiveSigmaThreshold(const std::vector<MapPt>& pts)
{
    std::vector<double> sigVals;
    for (const auto& p : pts)
        if (p.n_acc >= MIN_N_ACC && p.sigma > 0 && p.sigma < SIGMA_PRECUT)
            sigVals.push_back(p.sigma);
    if (sigVals.size() < 3) return SIGMA_PRECUT;
    std::sort(sigVals.begin(), sigVals.end());
    double med = sigVals[sigVals.size()/2];
    std::vector<double> ad; for (double s:sigVals) ad.push_back(std::abs(s-med));
    std::sort(ad.begin(), ad.end());
    double mad = ad[ad.size()/2];
    double thr = med + 10.0*mad;
    std::cout << "  [AdaptiveSigma] median=" << med*1000 << " ps  MAD=" << mad*1000
              << " ps  threshold=" << thr*1000 << " ps\n";
    return thr;
}

static double computeMuMedian(const std::vector<MapPt>& pts, double sigma_thr)
{
    std::vector<double> mv;
    for (const auto& p : pts)
        if (p.n_acc >= MIN_N_ACC && p.mu_ok && p.sigma > 0
            && p.sigma < sigma_thr && std::isfinite(p.mu))
            mv.push_back(p.mu);
    if (mv.empty()) return 0.0;
    std::sort(mv.begin(), mv.end());
    return mv[mv.size()/2];
}

// ── Filters ───────────────────────────────────────────────────────────────────
static bool pointValid(const MapPt& p) {
    return p.n_acc >= MIN_N_ACC && std::isfinite(p.sigma) && p.sigma > 0;
}
static bool pointGoodMu(const MapPt& p, double sigma_thr, double mu_med) {
    return pointValid(p) && p.mu_ok && p.sigma < sigma_thr
           && std::isfinite(p.mu)
           && std::abs(p.mu - mu_med) <= MU_MAXDEV_NS;
}
static bool pointGoodSigma(const MapPt& p, double sigma_thr) {
    return pointValid(p) && p.sigma < sigma_thr;
}

// ── Robust Z Scale MAD ────────────────────────────────────────────────────────
static std::pair<double,double> robustZRange(TH2D* h, double n_sigma=3.0,
                                              double min_half=0.0)
{
    std::vector<double> vals;
    for (int bx=1; bx<=h->GetNbinsX(); ++bx)
        for (int by=1; by<=h->GetNbinsY(); ++by) {
            double v = h->GetBinContent(bx,by);
            if (v < -1e8) continue;   // ignore empty sentinel bins
            vals.push_back(v);
        }
    if (vals.empty()) return {0.0,1.0};
    std::sort(vals.begin(),vals.end());
    double med = vals[vals.size()/2];
    std::vector<double> ad; for (double v:vals) ad.push_back(std::abs(v-med));
    std::sort(ad.begin(),ad.end());
    double mad = ad[ad.size()/2];
    double half = std::max(min_half, n_sigma*mad)*1.05;
    return {med-half, med+half};
}

// ── Axes Style ────────────────────────────────────────────────────────────────
static void styleAxes(TH2D* h)
{
    h->GetXaxis()->SetTitleSize(0.048f); h->GetYaxis()->SetTitleSize(0.048f);
    h->GetXaxis()->SetLabelSize(0.038f); h->GetYaxis()->SetLabelSize(0.038f);
    h->GetXaxis()->SetTitleOffset(1.1f); h->GetYaxis()->SetTitleOffset(1.35f);
    h->GetZaxis()->SetTitleSize(0.042f); h->GetZaxis()->SetLabelSize(0.034f);
    h->GetZaxis()->SetTitleOffset(1.65f);
}

// ── Reversed Pastel Rainbow Palette ───────────────────────────────────────────
static void setRainbowPalette()
{
    const int N = 8;
    double stops[] = {0.00, 0.14, 0.28, 0.43, 0.57, 0.71, 0.86, 1.00};
    double red[]   = {0.95, 0.95, 0.90, 0.55, 0.45, 0.50, 0.55, 0.70};
    double green[] = {0.55, 0.72, 0.90, 0.88, 0.85, 0.70, 0.55, 0.50};
    double blue[]  = {0.50, 0.45, 0.45, 0.50, 0.88, 0.92, 0.95, 0.92};
    TColor::CreateGradientColorTable(N, stops, red, green, blue, 128);
    gStyle->SetNumberContours(128);
}

// ── buildEdges ────────────────────────────────────────────────────────────────
static std::vector<double> buildEdges3d(const std::vector<double>& v)
{
    std::vector<double> e;
    if (v.size()==1){ e.push_back(v[0]-5); e.push_back(v[0]+5); return e; }
    e.reserve(v.size()+1);
    e.push_back(v[0]-0.5*(v[1]-v[0]));
    for (size_t i=1;i<v.size();++i) e.push_back(0.5*(v[i-1]+v[i]));
    e.push_back(v.back()+0.5*(v.back()-v[v.size()-2]));
    return e;
}

// ── Empty bin sentinel (avoid 0.0 being plotted) ─────────────────────────────
static void initEmptyBins(TH2D* h) {
    for (int bx=1; bx<=h->GetNbinsX(); ++bx)
        for (int by=1; by<=h->GetNbinsY(); ++by)
            h->SetBinContent(bx, by, -1e9);
}

// ── buildH2D (with sentinel initialisation) ──────────────────────────────────
static TH2D* buildH2D_mu(const std::vector<MapPt>& pts, int vb, double fp,
                          double sigma_thr, double mu_med)
{
    std::set<double> xs,ys;
    for (auto& p:pts) if (pointGoodMu(p,sigma_thr,mu_med)){xs.insert(p.x);ys.insert(p.y);}
    if (xs.empty()) return nullptr;
    std::vector<double> xv(xs.begin(),xs.end()),yv(ys.begin(),ys.end());
    auto xE=buildEdges3d(xv),yE=buildEdges3d(yv);
    TH2D* h=new TH2D(Form("hMu_v%d_f%.2f",vb,fp),
        Form("V_{bias}=%d V  LET=%.2f p.e.;#Deltax (mm);#Deltay (mm);#mu(#Deltat) (ns)",vb,fp),
        (int)xv.size(),xE.data(),(int)yv.size(),yE.data());
    h->SetDirectory(nullptr);
    initEmptyBins(h);                                      // <-- NEW
    for (auto& p:pts) if (pointGoodMu(p,sigma_thr,mu_med))
        h->SetBinContent(h->GetXaxis()->FindBin(p.x),h->GetYaxis()->FindBin(p.y),p.mu);
    return h;
}

static TH2D* buildH2D_sigma(const std::vector<MapPt>& pts, int vb, double fp,
                              double sigma_thr)
{
    std::set<double> xs,ys;
    for (auto& p:pts) if (pointGoodSigma(p,sigma_thr)){xs.insert(p.x);ys.insert(p.y);}
    if (xs.empty()) return nullptr;
    std::vector<double> xv(xs.begin(),xs.end()),yv(ys.begin(),ys.end());
    auto xE=buildEdges3d(xv),yE=buildEdges3d(yv);
    TH2D* h=new TH2D(Form("hSig_v%d_f%.2f",vb,fp),
        Form("V_{bias}=%d V  LET=%.2f p.e.;#Deltax (mm);#Deltay (mm);#sigma(#Deltat) (ns)",vb,fp),
        (int)xv.size(),xE.data(),(int)yv.size(),yE.data());
    h->SetDirectory(nullptr);
    initEmptyBins(h);                                      // <-- NEW
    for (auto& p:pts) if (pointGoodSigma(p,sigma_thr))
        h->SetBinContent(h->GetXaxis()->FindBin(p.x),h->GetYaxis()->FindBin(p.y),p.sigma);
    return h;
}

static TH2D* buildH2D_pdet(const std::vector<MapPt>& pts, int vb, double fp)
{
    std::set<double> xs,ys;
    for (auto& p:pts) if (pointValid(p)){xs.insert(p.x);ys.insert(p.y);}
    if (xs.empty()) return nullptr;
    std::vector<double> xv(xs.begin(),xs.end()),yv(ys.begin(),ys.end());
    auto xE=buildEdges3d(xv),yE=buildEdges3d(yv);
    TH2D* h=new TH2D(Form("hPdet_v%d_f%.2f",vb,fp),
        Form("V_{bias}=%d V  LET=%.2f p.e.;#Deltax (mm);#Deltay (mm);P_{det} (%%)",vb,fp),
        (int)xv.size(),xE.data(),(int)yv.size(),yE.data());
    h->SetDirectory(nullptr);
    initEmptyBins(h);                                      // <-- NEW
    for (auto& p:pts) if (pointValid(p) && p.p_det>=0)
        h->SetBinContent(h->GetXaxis()->FindBin(p.x),h->GetYaxis()->FindBin(p.y),p.p_det*100.0);
    return h;
}

// ════════════════════════════════════════════════════════════════════════════
//  drawSurf3D
// ════════════════════════════════════════════════════════════════════════════
static void drawSurf3D(TH2D* h, double zLo, double zHi,
                        const std::string& fname, OutCtx& ctx, int /*vb*/)
{
    if (!h) return;
    std::string fullTitle = h->GetTitle();
    std::string mainTitle = fullTitle;
    std::string zTitle    = "";
    size_t sc1 = fullTitle.find(';');
    if (sc1 != std::string::npos) {
        mainTitle = fullTitle.substr(0, sc1);
        size_t sc3 = fullTitle.rfind(';');
        if (sc3 != sc1) zTitle = fullTitle.substr(sc3+1);
    }

    setRainbowPalette();
    TCanvas* c = new TCanvas(Form("cS3d_%s",h->GetName()),mainTitle.c_str(),960,760);
    c->SetLeftMargin(0.10); c->SetRightMargin(0.14);
    c->SetBottomMargin(0.10); c->SetTopMargin(0.10);

    TH2D* hS = static_cast<TH2D*>(h->Clone(Form("hS3_%s",h->GetName())));
    hS->SetDirectory(nullptr);
    hS->SetTitle(mainTitle.c_str());
    
    // Push invalid/empty bins way below zLo so they are ignored by SURF2
    for (int bx=1; bx<=hS->GetNbinsX(); ++bx) {
        for (int by=1; by<=hS->GetNbinsY(); ++by) {
            double v = hS->GetBinContent(bx,by);
            if (v < -1e8 || std::isnan(v)) {
                hS->SetBinContent(bx, by, zLo - 1000.0);
            }
        }
    }

    hS->SetMinimum(zLo);
    hS->SetMaximum(zHi);
    hS->GetZaxis()->SetRangeUser(zLo, zHi);
    hS->GetXaxis()->SetTitle("#Deltax (mm)");
    hS->GetYaxis()->SetTitle("#Deltay (mm)");
    if (!zTitle.empty()) hS->GetZaxis()->SetTitle(zTitle.c_str());
    hS->GetXaxis()->SetTitleOffset(2.0f);
    hS->GetYaxis()->SetTitleOffset(2.0f);
    hS->GetZaxis()->SetTitleOffset(1.8f);
    hS->GetXaxis()->SetTitleSize(0.040f);
    hS->GetYaxis()->SetTitleSize(0.040f);
    hS->GetZaxis()->SetTitleSize(0.040f);
    hS->GetXaxis()->SetLabelSize(0.030f);
    hS->GetYaxis()->SetLabelSize(0.030f);
    hS->GetZaxis()->SetLabelSize(0.030f);
    hS->Draw("SURF2Z");
    c->Update(); c->Modified();
    ctx.savePNG(c, fname);
    delete hS;
}

// ════════════════════════════════════════════════════════════════════════════
//  drawColz2D
// ════════════════════════════════════════════════════════════════════════════
static void drawColz2D(TH2D* h, double zLo, double zHi,
                        const std::string& fname, OutCtx& ctx,
                        bool showText=false, bool withMarker=true)
{
    if (!h) return;
    setRainbowPalette();
    TCanvas* c = new TCanvas(Form("cC2d_%s",h->GetName()),"",840,720);
    c->SetRightMargin(0.20); c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.10);
    styleAxes(h);
    
    // Push empty bins (-1e9) below the real minimum so they become white
    for (int bx=1; bx<=h->GetNbinsX(); ++bx) {
        for (int by=1; by<=h->GetNbinsY(); ++by) {
            double v = h->GetBinContent(bx,by);
            if (v < -1e8 || std::isnan(v)) {
                h->SetBinContent(bx, by, zLo - 1000.0);
            }
        }
    }

    h->SetMinimum(zLo);
    h->SetMaximum(zHi);
    h->GetZaxis()->SetRangeUser(zLo, zHi);
    h->SetContour(128);
    h->Draw(showText ? "COLZ TEXT" : "COLZ");
    if (withMarker) {
        TMarker* mc = new TMarker(0.0, 0.0, 5);
        mc->SetMarkerColor(kBlack); mc->SetMarkerSize(2.0); mc->Draw();
    }
    c->Update(); c->Modified();
    ctx.savePNG(c, fname);
}

// ════════════════════════════════════════════════════════════════════════════
//  loadMapResults
// ════════════════════════════════════════════════════════════════════════════
static std::vector<MapPt> loadMapResults(const std::string& path, int vb_def)
{
    std::vector<MapPt> pts;
    TFile* f=TFile::Open(path.c_str(),"READ");
    if (!f||f->IsZombie()){delete f;return pts;}
    TTree* t=static_cast<TTree*>(f->Get("map"));
    if (!t){f->Close();delete f;return pts;}
    Double_t x=0,y=0,mu=0,mue=0,sig=0,sige=0,nsig=0,pdet=-1,chi2=-1;
    Long64_t nacc=0,ncross=0,nlas=0;
    Int_t vb=vb_def,mu_ok=1;
    t->SetBranchAddress("x",&x); t->SetBranchAddress("y",&y);
    t->SetBranchAddress("mu",&mu); t->SetBranchAddress("sigma",&sig);
    t->SetBranchAddress("n_sig",&nsig); t->SetBranchAddress("n_acc",&nacc);
    t->SetBranchAddress("n_laser",&nlas); t->SetBranchAddress("p_det",&pdet);
    if (t->GetBranch("mu_err"))    t->SetBranchAddress("mu_err",&mue);
    if (t->GetBranch("mu_ok"))     t->SetBranchAddress("mu_ok",&mu_ok);
    if (t->GetBranch("sigma_err")) t->SetBranchAddress("sigma_err",&sige);
    if (t->GetBranch("n_crossing"))t->SetBranchAddress("n_crossing",&ncross);
    if (t->GetBranch("vbias"))     t->SetBranchAddress("vbias",&vb);
    if (t->GetBranch("chi2ndf"))   t->SetBranchAddress("chi2ndf",&chi2);
    const Long64_t N=t->GetEntries();
    pts.reserve((size_t)N);
    for (Long64_t i=0;i<N;++i){
        t->GetEntry(i);
        pts.push_back({x,y,vb,mu,mue,sig,sige,nsig,pdet,
                       (long)nacc,(long)ncross,(long)nlas,chi2,mu_ok});
    }
    f->Close(); delete f;
    std::cout<<"  Loaded "<<pts.size()<<" points from "
             <<path.substr(path.find_last_of("/\\")+1)<<"\n";
    return pts;
}

// ════════════════════════════════════════════════════════════════════════════
//  makeGraph2D
// ════════════════════════════════════════════════════════════════════════════
static TGraph2D* makeGraph2D(TH2D* h, int vb, int col)
{
    if (!h) return nullptr;
    int np=0;
    for (int bx=1;bx<=h->GetNbinsX();++bx)
        for (int by=1;by<=h->GetNbinsY();++by) {
            double v=h->GetBinContent(bx,by);
            if (v > -1e8 && !std::isnan(v)) ++np;
        }
    if (np==0) return nullptr;
    TGraph2D* g=new TGraph2D(np);
    g->SetName(Form("g3d_v%d",vb));
    g->SetMarkerColor(col); g->SetLineColor(col);
    g->SetMarkerStyle(20); g->SetMarkerSize(1.4);
    int ip=0;
    for (int bx=1;bx<=h->GetNbinsX();++bx)
        for (int by=1;by<=h->GetNbinsY();++by){
            double v=h->GetBinContent(bx,by);
            if (v < -1e8 || std::isnan(v)) continue;
            g->SetPoint(ip++,h->GetXaxis()->GetBinCenter(bx),
                             h->GetYaxis()->GetBinCenter(by),v);
        }
    return g;
}

// ════════════════════════════════════════════════════════════════════════════
//  drawDelayMap
// ════════════════════════════════════════════════════════════════════════════
static void drawDelayMap(const std::vector<MapPt>& pts,
                          int vbias, double frac_pe,
                          OutCtx& ctx, const std::string& tag,
                          double sigma_thr, double mu_med)
{
    struct DelayPt { double x,y,delay_ps,err_ps; bool outlier; };

    double mu_ref=0,mu_ref_err=0,minDist=1e18;
    bool found=false;
    for (auto& p:pts){
        if (!pointGoodMu(p,sigma_thr,mu_med)) continue;
        double d=std::hypot(p.x,p.y);
        if (d<minDist){minDist=d;mu_ref=p.mu;mu_ref_err=p.mu_err;found=true;}
    }
    if (!found){std::cout<<"  [DelayMap] no valid points, skip.\n";return;}
    std::cout<<"  [DelayMap] mu_ref="<<std::fixed<<std::setprecision(4)
             <<mu_ref<<" ns  dist_center="<<std::setprecision(1)<<minDist<<" mm\n";

    std::vector<DelayPt> dpts;
    for (auto& p:pts){
        if (!pointGoodMu(p,sigma_thr,mu_med)) continue;
        if (std::isnan(p.mu)) continue;
        double e=std::hypot(p.mu_err,mu_ref_err)*1000.0;
        double d_ps = (p.mu-mu_ref)*1000.0;
        if (std::isnan(d_ps)) continue;
        dpts.push_back({p.x,p.y,d_ps,e,false});
    }
    if (dpts.size()<2) return;

    std::vector<double> dv; for (auto& d:dpts) dv.push_back(d.delay_ps);
    std::sort(dv.begin(),dv.end());
    double med=dv[dv.size()/2];
    std::vector<double> ad; for (double v:dv) ad.push_back(std::abs(v-med));
    std::sort(ad.begin(),ad.end());
    double mad=ad[ad.size()/2];
    double thresh=3.0*mad;
    for (auto& d:dpts) d.outlier=(thresh>0 && std::abs(d.delay_ps-med)>thresh);

    double maxGood=0.0, minGood=1e18;
    for (auto& d:dpts) if (!d.outlier){
        maxGood=std::max(maxGood,d.delay_ps);
        minGood=std::min(minGood,d.delay_ps);
    }
    
    double zLo = (minGood < -0.5*maxGood) ? minGood*1.15 : 0.0;
    double zHi = std::max(5.0, maxGood*1.15);

    std::cout<<"  [DelayMap] med="<<std::fixed<<std::setprecision(1)
             <<med<<" ps  MAD="<<mad<<" ps  range=["<<zLo<<","<<zHi<<"] ps\n";

    std::set<double> uxs,uys;
    for (auto& d:dpts) if (!d.outlier){uxs.insert(d.x);uys.insert(d.y);}
    if (uxs.empty()) return;
    std::vector<double> xv(uxs.begin(),uxs.end()),yv(uys.begin(),uys.end());
    auto xE=buildEdges3d(xv),yE=buildEdges3d(yv);
    double xLo=xv.front()-1,xHi=xv.back()+1;
    double yLo=yv.front()-1,yHi=yv.back()+1;

    // ── Canvas 1: map 2D COLZ ───────────────────────────────────────────────
    {
        TH2D* hD=new TH2D(Form("hDlay_v%d_f%.2f",vbias,frac_pe),
            Form("V_{bias}=%d V  LET=%.2f p.e.  #mu_{ref}=%.4f ns;"
                 "#Deltax (mm);#Deltay (mm);#Delta#mu (ps)",vbias,frac_pe,mu_ref),
            (int)xv.size(),xE.data(),(int)yv.size(),yE.data());
        hD->SetDirectory(nullptr);
        initEmptyBins(hD);                                 // <-- NEW
        
        for (auto& d:dpts){
            if (d.outlier) continue;
            int bx=hD->GetXaxis()->FindBin(d.x);
            int by=hD->GetYaxis()->FindBin(d.y);
            hD->SetBinContent(bx, by, d.delay_ps);
        }

        // replace old drawColz2D_fromZero with unified function
        drawColz2D(hD, zLo, zHi, "delay_map_2d_"+tag+".png", ctx, false, false);

        // ── Canvas 3D delay ──────────────────────────────────────────────────
        {
            setRainbowPalette();
            TH2D* hD3 = static_cast<TH2D*>(hD->Clone(Form("hDlay3d_%s",tag.c_str())));
            hD3->SetDirectory(nullptr);
            hD3->SetTitle(Form("V_{bias}=%d V  LET=%.2f p.e.  #mu_{ref}=%.4f ns",
                               vbias, frac_pe, mu_ref));
            hD3->GetZaxis()->SetTitle("#Delta#mu (ps)");
            
            for (int bx=1; bx<=hD3->GetNbinsX(); ++bx) {
                for (int by=1; by<=hD3->GetNbinsY(); ++by) {
                    double v = hD3->GetBinContent(bx,by);
                    if (v < -1e8 || std::isnan(v)) hD3->SetBinContent(bx, by, zLo - 1000.0);
                }
            }

            hD3->SetMinimum(zLo);                         // <-- use zLo, not 1e-9
            hD3->SetMaximum(zHi);
            hD3->GetZaxis()->SetRangeUser(zLo, zHi);
            
            TCanvas* c3 = new TCanvas(Form("cDlay3d_%s",tag.c_str()),
                Form("V_{bias}=%d V  LET=%.2f p.e. — delay map 3D",vbias,frac_pe),
                960, 760);
            c3->SetLeftMargin(0.10); c3->SetRightMargin(0.14);
            c3->SetBottomMargin(0.10); c3->SetTopMargin(0.10);
            hD3->GetXaxis()->SetTitle("#Deltax (mm)");
            hD3->GetYaxis()->SetTitle("#Deltay (mm)");
            hD3->GetXaxis()->SetTitleOffset(2.0f); hD3->GetYaxis()->SetTitleOffset(2.0f);
            hD3->GetZaxis()->SetTitleOffset(1.8f);
            hD3->GetXaxis()->SetTitleSize(0.040f); hD3->GetYaxis()->SetTitleSize(0.040f);
            hD3->GetZaxis()->SetTitleSize(0.040f);
            hD3->GetXaxis()->SetLabelSize(0.030f); hD3->GetYaxis()->SetLabelSize(0.030f);
            hD3->GetZaxis()->SetLabelSize(0.030f);
            hD3->Draw("SURF2Z");
            c3->Update(); c3->Modified();
            ctx.savePNG(c3, "delay_map_3d_"+tag+".png");
            delete hD3;
        }
        delete hD;
    }

    // ── Linear Fit ────────────────────────────────────────────────────────────
    std::vector<double> gx_x,gx_y,gx_ex,gx_ey;
    std::vector<double> gy_x,gy_y,gy_ex,gy_ey;
    for (auto& d:dpts){
        if (d.outlier) continue;
        gx_x.push_back(d.x); gx_y.push_back(d.delay_ps);
        gx_ex.push_back(0);   gx_ey.push_back(d.err_ps);
        gy_x.push_back(d.y); gy_y.push_back(d.delay_ps);
        gy_ex.push_back(0);   gy_ey.push_back(d.err_ps);
    }
    if (gx_x.empty()) return;

    TGraphErrors grX((int)gx_x.size(),gx_x.data(),gx_y.data(),gx_ex.data(),gx_ey.data());
    grX.SetMarkerStyle(20); grX.SetMarkerSize(1.4);
    grX.SetMarkerColor(kAzure+1); grX.SetLineColor(kAzure+1); grX.SetLineWidth(2);

    TGraphErrors grY((int)gy_x.size(),gy_x.data(),gy_y.data(),gy_ex.data(),gy_ey.data());
    grY.SetMarkerStyle(20); grY.SetMarkerSize(1.4);
    grY.SetMarkerColor(kOrange+7); grY.SetLineColor(kOrange+7); grY.SetLineWidth(2);

    TF1 fLY(Form("fLY_%s",tag.c_str()),"pol1",yLo,yHi);
    fLY.SetLineColor(kRed+1); fLY.SetLineWidth(2);
    grY.Fit(&fLY,"RQ");

    double sY=fLY.GetParameter(1),seY=fLY.GetParError(1);
    double thY=std::asin(std::min(1.0,std::abs(sY)*C_LIGHT_MMPS))*180.0/TMath::Pi();
    std::cout<<"  [DelayMap] slope_y="<<std::setprecision(3)
             <<sY<<"±"<<seY<<" ps/mm  θ_y="<<thY<<"°\n";

    std::vector<double> gxL_x,gxL_y,gxL_ex,gxL_ey;
    std::vector<double> gxR_x,gxR_y,gxR_ex,gxR_ey;
    for (auto& d:dpts){
        if (d.outlier) continue;
        if (d.x < 0) {
            gxL_x.push_back(d.x); gxL_y.push_back(d.delay_ps);
            gxL_ex.push_back(0);  gxL_ey.push_back(d.err_ps);
        } else {
            gxR_x.push_back(d.x); gxR_y.push_back(d.delay_ps);
            gxR_ex.push_back(0);  gxR_ey.push_back(d.err_ps);
        }
    }

    double sXL=0,seXL=0,thXL=0;
    double sXR=0,seXR=0,thXR=0;
    TF1* fLXL = nullptr;
    TF1* fLXR = nullptr;

    if (gxL_x.size() >= 2) {
        fLXL = new TF1(Form("fLXL_%s",tag.c_str()),"pol1",xLo,-0.5);
        fLXL->SetLineColor(kBlue+1); fLXL->SetLineWidth(2);
        TGraphErrors grXL((int)gxL_x.size(),gxL_x.data(),gxL_y.data(),
                           gxL_ex.data(),gxL_ey.data());
        grXL.Fit(fLXL,"RQ");
        sXL  = fLXL->GetParameter(1);
        seXL = fLXL->GetParError(1);
        thXL = std::asin(std::min(1.0,std::abs(sXL)*C_LIGHT_MMPS))*180.0/TMath::Pi();
    }
    if (gxR_x.size() >= 2) {
        fLXR = new TF1(Form("fLXR_%s",tag.c_str()),"pol1",0.5,xHi);
        fLXR->SetLineColor(kRed+1); fLXR->SetLineWidth(2);
        TGraphErrors grXR((int)gxR_x.size(),gxR_x.data(),gxR_y.data(),
                           gxR_ex.data(),gxR_ey.data());
        grXR.Fit(fLXR,"RQ");
        sXR  = fLXR->GetParameter(1);
        seXR = fLXR->GetParError(1);
        thXR = std::asin(std::min(1.0,std::abs(sXR)*C_LIGHT_MMPS))*180.0/TMath::Pi();
    }
    std::cout<<"  [DelayMap] slope_x_left ="<<std::fixed<<std::setprecision(3)
             <<sXL<<"±"<<seXL<<" ps/mm  θ="<<std::setprecision(2)<<thXL<<"°\n"
             <<"  [DelayMap] slope_x_right="<<std::setprecision(3)
             <<sXR<<"±"<<seXR<<" ps/mm  θ="<<std::setprecision(2)<<thXR<<"°\n";

    double yPlo=std::min(0.0,minGood*1.2)-2.0;
    double yPhi=maxGood*1.2+2.0;

    // ── Canvas delay vs Δx ────────────────────────────────────────────────────
    {
        TGraphErrors grX((int)gx_x.size(),gx_x.data(),gx_y.data(),
                          gx_ex.data(),gx_ey.data());
        grX.SetMarkerStyle(20); grX.SetMarkerSize(1.4);
        grX.SetMarkerColor(kAzure+1); grX.SetLineColor(kAzure+1); grX.SetLineWidth(2);
        grX.SetTitle(Form("V_{bias}=%d V  LET=%.2f p.e.;"
                          "#Deltax (mm);#Delta#mu (ps)", vbias, frac_pe));

        TCanvas* cx = new TCanvas(Form("cDMx_%s",tag.c_str()),
            Form("V_{bias}=%d V  LET=%.2f p.e. — delay vs #Deltax",vbias,frac_pe),
            900,640);
        cx->SetGrid();
        cx->SetLeftMargin(0.13); cx->SetRightMargin(0.06);
        cx->SetBottomMargin(0.14); cx->SetTopMargin(0.12);
        grX.GetYaxis()->SetRangeUser(yPlo,yPhi);
        grX.GetYaxis()->SetTitleSize(0.052f); grX.GetYaxis()->SetTitleOffset(1.10f);
        grX.GetYaxis()->SetLabelSize(0.044f);
        grX.GetXaxis()->SetTitleSize(0.052f); grX.GetXaxis()->SetTitleOffset(1.00f);
        grX.GetXaxis()->SetLabelSize(0.044f);
        grX.Draw("AP");
        if (fLXL) fLXL->Draw("same");
        if (fLXR) fLXR->Draw("same");
        TLine lZ(xLo,0.0,xHi,0.0);
        lZ.SetLineStyle(2); lZ.SetLineColor(kGray+2); lZ.SetLineWidth(1); lZ.Draw();
        TLine lC(0.0,yPlo,0.0,yPhi);
        lC.SetLineStyle(3); lC.SetLineColor(kGray+1); lC.SetLineWidth(1); lC.Draw();
        TLatex lt; lt.SetNDC(); lt.SetTextFont(42); lt.SetTextSize(0.038);
        lt.SetTextColor(kBlue+1);
        lt.DrawLatex(0.14,0.88,Form("slope (x<0) = %.2f #pm %.2f ps/mm  (#theta=%.2f#circ)",
                                     sXL,seXL,thXL));
        lt.SetTextColor(kRed+1);
        lt.DrawLatex(0.14,0.83,Form("slope (x>0) = %.2f #pm %.2f ps/mm  (#theta=%.2f#circ)",
                                     sXR,seXR,thXR));
        lt.SetTextSize(0.034); lt.SetTextColor(kGray+1);
        lt.DrawLatex(0.14,0.78,Form("#mu_{ref} = %.4f ns",mu_ref));
        lt.SetTextColor(kBlack);
        cx->Update(); cx->Modified();
        ctx.savePNG(cx,"delay_map_x_"+tag+".png");
    }
    if (fLXL) delete fLXL;
    if (fLXR) delete fLXR;

    // ── Canvas delay vs Δy ────────────────────────────────────────────────────
    {
        grY.SetTitle(Form("V_{bias}=%d V  LET=%.2f p.e.;"
                          "#Deltay (mm);#Delta#mu (ps)", vbias, frac_pe));
        TCanvas* cy = new TCanvas(Form("cDMy_%s",tag.c_str()),
            Form("V_{bias}=%d V  LET=%.2f p.e. — delay vs #Deltay",vbias,frac_pe),
            900,640);
        cy->SetGrid();
        cy->SetLeftMargin(0.13); cy->SetRightMargin(0.06);
        cy->SetBottomMargin(0.14); cy->SetTopMargin(0.12);
        grY.GetYaxis()->SetRangeUser(yPlo,yPhi);
        grY.GetYaxis()->SetTitleSize(0.052f); grY.GetYaxis()->SetTitleOffset(1.10f);
        grY.GetYaxis()->SetLabelSize(0.044f);
        grY.GetXaxis()->SetTitleSize(0.052f); grY.GetXaxis()->SetTitleOffset(1.00f);
        grY.GetXaxis()->SetLabelSize(0.044f);
        grY.Draw("AP");
        fLY.Draw("same");
        TLine lZy(yLo,0.0,yHi,0.0);
        lZy.SetLineStyle(2); lZy.SetLineColor(kGray+2); lZy.SetLineWidth(1); lZy.Draw();
        TLatex lt; lt.SetNDC(); lt.SetTextFont(42); lt.SetTextSize(0.044);
        lt.DrawLatex(0.16,0.88,Form("slope = %.2f #pm %.2f ps/mm",sY,seY));
        lt.SetTextSize(0.036); lt.SetTextColor(kGray+1);
        lt.DrawLatex(0.16,0.82,Form("#theta = %.2f#circ   (c = 0.300 mm/ps)",thY));
        lt.DrawLatex(0.16,0.77,Form("#mu_{ref} = %.4f ns",mu_ref));
        lt.SetTextColor(kBlack);
        cy->Update(); cy->Modified();
        ctx.savePNG(cy,"delay_map_y_"+tag+".png");
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════════════════════
void sipm_draw_timing3d()
{
    g_data_dir_override="";
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(64);

    std::cout<<"\n+==========================================================+\n"
             <<"|     SiPM TIMING 3D — post-process sipm_pos_scan output  |\n"
             <<"+==========================================================+\n\n";

    std::string dataDir;
    {
        std::string root=DATA_DIR;
        while (!root.empty()&&root.back()=='/') root.pop_back();
        const size_t sl=root.find_last_of("/\\");
        if (sl!=std::string::npos) root=root.substr(0,sl);

        std::vector<std::string> subs;
        void* dp=gSystem->OpenDirectory(root.c_str());
        if (dp){
            const char* ent=nullptr;
            while ((ent=gSystem->GetDirEntry(dp))!=nullptr){
                std::string s(ent); if (s=="."||s=="..") continue;
                FileStat_t st;
                if (gSystem->GetPathInfo((root+"/"+s).c_str(),st)==0&&R_ISDIR(st.fMode))
                    subs.push_back(s);
            }
            gSystem->FreeDirectory(dp);
        }
        std::sort(subs.begin(),subs.end());

        if (!subs.empty()){
            std::cout<<"  Folders in "<<root<<":\n";
            for (size_t i=0;i<subs.size();++i)
                std::cout<<"    ["<<(i+1)<<"] "<<subs[i]<<(i==0?"   <-- default":"")<<"\n";
            const std::string l=readLineOrEmpty(Form("\n  Choose [n] or path  [ENTER = 1]: "));
            if (l.empty()) dataDir=root+"/"+subs[0];
            else {
                try{ size_t pos=0; int idx=std::stoi(l,&pos);
                     if (pos==l.size()&&idx>=1&&idx<=(int)subs.size())
                         dataDir=root+"/"+subs[idx-1];
                     else dataDir=l;
                } catch(...){dataDir=l;}
            }
        } else dataDir=readLine("  Data folder path: ");

        while (!dataDir.empty()&&dataDir.back()=='/') dataDir.pop_back();
        if (gSystem->AccessPathName(dataDir.c_str())){
            std::cerr<<"[ERR] Not accessible: "<<dataDir<<"\n"; return;
        }
        g_data_dir_override=dataDir;
        std::cout<<"  --> "<<dataDir<<"\n\n";
    }

    struct ResultFile{ std::string path; int vbias; double frac_pe; std::string mode; };
    std::vector<ResultFile> available;
    {
        static const std::regex re(
            R"(^map_results_vbias(\d+)_let([\d\.]+?)((?:_loose)?(?:_pdet_(?:crossing|accepted))?)\.root$)");
        void* dp=gSystem->OpenDirectory(dataDir.c_str());
        if (!dp){std::cerr<<"[ERR] Cannot open dir.\n";return;}
        const char* ent=nullptr;
        while ((ent=gSystem->GetDirEntry(dp))!=nullptr){
            std::string fn(ent); std::smatch m;
            if (!std::regex_match(fn,m,re)) continue;
            try{
                int vb=std::stoi(m[1].str()); double fr=std::stod(m[2].str());
                std::string suf=m[3].str(),mode="";
                if (suf.find("_pdet_crossing")!=std::string::npos) mode="crossing";
                else if (suf.find("_pdet_accepted")!=std::string::npos) mode="accepted";
                available.push_back({dataDir+"/"+fn,vb,fr,mode});
            } catch(...){}
        }
        gSystem->FreeDirectory(dp);
    }
    if (available.empty()){
        std::cerr<<"[ERR] No map_results_*.root in "<<dataDir<<"\n"; return;
    }
    std::sort(available.begin(),available.end(),
        [](const ResultFile& a,const ResultFile& b){
            return a.vbias!=b.vbias?a.vbias<b.vbias:a.frac_pe<b.frac_pe;});

    std::cout<<"  Available map_results files:\n";
    for (size_t i=0;i<available.size();++i)
        std::cout<<"    ["<<(i+1)<<"] vbias="<<available[i].vbias
                 <<"  let="<<std::fixed<<std::setprecision(2)<<available[i].frac_pe
                 <<(!available[i].mode.empty()?"  ["+available[i].mode+"]":"")<<"\n";

    std::vector<int> chosen;
    {
        const std::string l=readLine("\n  Indices to plot (e.g. 1 2 3) or 'all': ");
        if (l=="all"||l=="ALL"){
            for (size_t i=0;i<available.size();++i) chosen.push_back((int)i);
        } else {
            std::stringstream ss(l); int v=0;
            while (ss>>v) if (v>=1&&v<=(int)available.size()) chosen.push_back(v-1);
        }
    }
    if (chosen.empty()){std::cerr<<"[ERR] No selection.\n";return;}

    double cx=90.0,cy=0.0;
    {
        const std::string l=readLine("\n  Nominal center (x y) [default 90 0]: ");
        std::istringstream ss(l); double a=0,b=0;
        if (ss>>a>>b){cx=a;cy=b;}
    }
    std::cout<<"  --> Center = ("<<cx<<", "<<cy<<")\n\n";

    OutCtx ctx=createOutputDirs("timing3d");

    std::map<double,std::map<int,std::vector<MapPt>>> byFracByVbias;
    for (int ci:chosen){
        auto pts=loadMapResults(available[ci].path,available[ci].vbias);
        if (!pts.empty()) byFracByVbias[available[ci].frac_pe][available[ci].vbias]=pts;
    }

    for (auto& [frac_pe,vbMap]:byFracByVbias){
        for (auto& [vbias,pts]:vbMap){
            std::cout<<"\n--- Vbias="<<vbias
                     <<"  frac_pe="<<std::fixed<<std::setprecision(2)<<frac_pe<<" ---\n";

            std::vector<MapPt> pS=pts;
            for (auto& p:pS){p.x-=cx;p.y-=cy;}

            const double sigma_thr = computeAdaptiveSigmaThreshold(pS);
            const double mu_med    = computeMuMedian(pS, sigma_thr);
            std::cout << "  [Filter] sigma_thr=" << sigma_thr*1000
                      << " ps  mu_med=" << std::fixed << std::setprecision(4)
                      << mu_med << " ns\n";

            const std::string tag=Form("vbias%d_let%.2f",vbias,frac_pe);

            TH2D* hMu=buildH2D_mu(pS,vbias,frac_pe,sigma_thr,mu_med);
            if (!hMu){std::cout<<"  [SKIP] no valid points.\n";continue;}
            auto [muLo,muHi]=robustZRange(hMu,3.0,0.010);
            drawSurf3D(hMu,muLo,muHi,"mu_3d_"+tag+".png",ctx,vbias);
            drawColz2D(hMu,muLo,muHi,"mu_2d_"+tag+".png",ctx);
            delete hMu;

            TH2D* hSig=buildH2D_sigma(pS,vbias,frac_pe,sigma_thr);
            if (hSig){
                auto [sLo,sHi]=robustZRange(hSig,3.0,0.005);
                drawSurf3D(hSig,sLo,sHi,"sigma_3d_"+tag+".png",ctx,vbias);
                drawColz2D(hSig,sLo,sHi,"sigma_2d_"+tag+".png",ctx);
                delete hSig;
            }

            // ── P_det ─────────────────────────────────────────────────────────
            {
                static constexpr double MIN_PDET_DISPLAY = 0.1; // % under this becomes white
                TH2D* hPdet=buildH2D_pdet(pS,vbias,frac_pe);
                if (hPdet){
                    // bin below threshold become empty (sentinel) → white
                    for (int bx=1;bx<=hPdet->GetNbinsX();++bx)
                        for (int by=1;by<=hPdet->GetNbinsY();++by){
                            double v=hPdet->GetBinContent(bx,by);
                            if (v > -1e8 && v < MIN_PDET_DISPLAY)
                                hPdet->SetBinContent(bx,by,-1e9);
                        }
                    double pMax=hPdet->GetMaximum();
                    drawColz2D(hPdet, 0.0, std::min(100.0,pMax*1.10+1.0),
                               "pdet_"+tag+".png", ctx, true, true);
                    delete hPdet;
                }
            }

            // ── Delay map ─────────────────────────────────────────────────────
            drawDelayMap(pS,vbias,frac_pe,ctx,tag,sigma_thr,mu_med);
        }

        // ── Multi-Vbias Overlay ───────────────────────────────────────────────
        if (vbMap.size()>1){
            std::map<int,TH2D*> hMuMap;
            for (auto& [vb,ptsRaw]:vbMap){
                std::vector<MapPt> pO=ptsRaw;
                for (auto& p:pO){p.x-=cx;p.y-=cy;}
                double st=computeAdaptiveSigmaThreshold(pO);
                double mm=computeMuMedian(pO,st);
                TH2D* h=buildH2D_mu(pO,vb,frac_pe,st,mm);
                if (h) hMuMap[vb]=h;
            }
            std::vector<double> allV;
            for (auto& [v,h]:hMuMap)
                for (int bx=1;bx<=h->GetNbinsX();++bx)
                    for (int by=1;by<=h->GetNbinsY();++by){
                        double val=h->GetBinContent(bx,by);
                        if (val > -1e8 && !std::isnan(val)) allV.push_back(val);
                    }
            double gzLo=-1,gzHi=1;
            if (!allV.empty()){
                std::sort(allV.begin(),allV.end());
                double gm=allV[allV.size()/2];
                std::vector<double> ga; for (double v:allV) ga.push_back(std::abs(v-gm));
                std::sort(ga.begin(),ga.end());
                double gmad=ga[ga.size()/2];
                gzLo=gm-3.0*gmad*1.05; gzHi=gm+3.0*gmad*1.05;
            }
            
            TCanvas* cO=new TCanvas(Form("cOvl_let%.2f",frac_pe),"",1100,750);
            cO->SetRightMargin(0.20); 
            TLegend* leg=new TLegend(0.82,0.65,0.99,0.95);
            leg->SetBorderSize(1); leg->SetFillColor(0); leg->SetTextSize(0.025);
            bool first=true; int oi=0;
            std::vector<TGraph2D*> alive;
            for (auto& [vb,h]:hMuMap){
                TGraph2D* g=makeGraph2D(h,vb,vbiasColor(oi++));
                if (!g) continue;
                alive.push_back(g);
                if (first){
                    g->SetTitle(Form("LET=%.2f pe;#Deltax (mm);#Deltay (mm);#mu(#Deltat) (ns)",frac_pe));
                    g->GetZaxis()->SetRangeUser(gzLo,gzHi);
                    g->Draw("P0"); first=false;
                } else g->Draw("P0 SAME");
                leg->AddEntry(g,Form("V_{bias} = %d V",vb),"p");
            }
            if (!first) leg->Draw();
            cO->Update();
            ctx.savePNG(cO,Form("overlay3d_let%.2f.png",frac_pe));
            for (auto* g:alive) delete g;
            for (auto& [v,h]:hMuMap) delete h;
        }
    }

    std::cout<<"\n+==========================================================+\n"
             <<"|  DONE — canvas in: "<<ctx.pngDir<<"\n"
             <<"+==========================================================+\n";
}