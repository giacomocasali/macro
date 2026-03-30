// read_all_scans.cpp  —  v4.0
// Detects all scan_vbiasXX_filterYY.root files and overlays them on two canvases.
//
// Improvements over v3.0:
//   [NEW-1]  Derivative curves normalised to their own maximum → all visible
//   [NEW-2]  X range trimmed to last threshold with meaningful signal
//   [NEW-3]  Scan canvas uses log Y scale to show all three curves together
//   [NEW-4]  Lines only (no markers) → cleaner overlay with many curves
//   [NEW-5]  Derivative Y zoom per-file aware, not driven by the tallest file

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TLine.h>

constexpr double MIN_POSITIVE_THR = 5.0;
constexpr double MIN_THR_DISPLAY  = -10.0;

constexpr float PAD_LEFT   = 0.13f;
constexpr float PAD_RIGHT  = 0.04f;
constexpr float PAD_TOP    = 0.08f;
constexpr float PAD_BOTTOM = 0.12f;

// ════════════════════════════════════════════════════════════
//  FILE PARSING  scan_vbias{V}_filter{F}.root
// ════════════════════════════════════════════════════════════
struct FileInfo { std::string path; double vbias, filter; };

static bool parseFilename(const std::string& fname, FileInfo& info) {
    const std::string pre = "scan_vbias";
    size_t p0 = fname.find(pre);
    if (p0 == std::string::npos) return false;
    size_t vS = p0 + pre.size(), vE = fname.find("_filter", vS);
    if (vE == std::string::npos) return false;
    size_t lS = vE + 7, lE = fname.find(".root", lS);
    if (lE == std::string::npos) return false;
    try {
        info.vbias  = std::stod(fname.substr(vS, vE - vS));
        info.filter = std::stod(fname.substr(lS, lE - lS));
    } catch (...) { return false; }
    info.path = fname;
    return true;
}

// ════════════════════════════════════════════════════════════
//  COLOUR PALETTE
// ════════════════════════════════════════════════════════════
static int fileColor(int idx) {
    static const int pal[] = {
        kAzure+1, kOrange+7, kGreen+2, kMagenta+1,
        kRed+1,   kCyan+2,   kViolet+1, kYellow+3,
    };
    return pal[idx % (int)(sizeof(pal)/sizeof(pal[0]))];
}

// ════════════════════════════════════════════════════════════
//  LOAD SCAN → TGraph (counts only, no y-errors for overlay)
// ════════════════════════════════════════════════════════════
static TGraph* loadScan(TFile* f) {
    TTree* t = (TTree*)f->Get("scan");
    if (!t) return nullptr;
    double thr, cnt;
    t->SetBranchAddress("threshold", &thr);
    t->SetBranchAddress("counts",    &cnt);
    Long64_t n = t->GetEntries();
    std::vector<double> vx(n), vy(n);
    for (Long64_t i = 0; i < n; ++i) { t->GetEntry(i); vx[i]=thr; vy[i]=cnt; }
    return new TGraph(n, vx.data(), vy.data());
}

// ════════════════════════════════════════════════════════════
//  LOAD DERIV → TGraph, NORMALISED to its maximum in
//  the physical region (thr > MIN_POSITIVE_THR)
// ════════════════════════════════════════════════════════════
static TGraph* loadDerivNorm(TFile* f) {
    TTree* t = (TTree*)f->Get("deriv");
    if (!t) return nullptr;
    double thr, der;
    t->SetBranchAddress("threshold", &thr);
    t->SetBranchAddress("deriv",     &der);
    Long64_t n = t->GetEntries();
    std::vector<double> vx(n), vy(n);
    for (Long64_t i = 0; i < n; ++i) { t->GetEntry(i); vx[i]=thr; vy[i]=der; }

    // Find maximum in physical region for normalisation
    double peak = 0;
    for (int i = 0; i < (int)n; ++i)
        if (vx[i] > MIN_POSITIVE_THR && vy[i] > peak) peak = vy[i];
    if (peak <= 0) return new TGraph(n, vx.data(), vy.data());
    for (auto& v : vy) v /= peak;
    return new TGraph(n, vx.data(), vy.data());
}

// ════════════════════════════════════════════════════════════
//  TRIM X RANGE
//  Right edge = last threshold where ANY file has counts > 0.5%
//  of that file's maximum — avoids trailing flat zero regions.
// ════════════════════════════════════════════════════════════
static double trimmedXMax(const std::vector<TGraph*>& gs) {
    double xmax = MIN_THR_DISPLAY;
    for (auto* g : gs) {
        if (!g) continue;
        // Find this graph's max Y
        double ymax = 0;
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y; g->GetPoint(i, x, y);
            if (y > ymax) ymax = y;
        }
        double thresh = ymax * 0.005; // 0.5% of peak
        // Walk from the right until we find a point above thresh
        for (int i = g->GetN()-1; i >= 0; --i) {
            double x, y; g->GetPoint(i, x, y);
            if (y > thresh) { if (x > xmax) xmax = x; break; }
        }
    }
    return xmax;
}

// ════════════════════════════════════════════════════════════
//  AXIS STYLE
// ════════════════════════════════════════════════════════════
static void axisStyle(TGraph* g) {
    if (!g) return;
    g->GetXaxis()->SetTitleSize(0.055); g->GetXaxis()->SetLabelSize(0.048);
    g->GetYaxis()->SetTitleSize(0.055); g->GetYaxis()->SetLabelSize(0.048);
    g->GetYaxis()->SetTitleOffset(1.05);
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════
void read_all_scans() {
    gStyle->SetOptStat(0);

    // ── Find files ───────────────────────────────────────────
    void* dh = gSystem->OpenDirectory(".");
    if (!dh) { std::cerr << "Cannot open current directory.\n"; return; }
    std::vector<FileInfo> files;
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dh)) != nullptr) {
        std::string fname(entry);
        if (fname.find("scan_vbias") == std::string::npos) continue;
        if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
        FileInfo fi;
        if (parseFilename(fname, fi)) {
            files.push_back(fi);
            std::cout << "Found: " << fname
                      << "  vbias=" << fi.vbias
                      << " V,  filter=" << fi.filter << "\n";
        }
    }
    gSystem->FreeDirectory(dh);
    if (files.empty()) {
        std::cerr << "No scan_vbiasXX_filterYY.root files found.\n"; return;
    }
    std::sort(files.begin(), files.end(), [](const FileInfo& a, const FileInfo& b){
        return (a.vbias != b.vbias) ? a.vbias < b.vbias : a.filter < b.filter;
    });
    int nf = (int)files.size();
    std::cout << "\n  " << nf << " file(s) found.\n\n";

    // ── Canvases ─────────────────────────────────────────────
    TCanvas* cScan  = new TCanvas("cScan",  "Threshold scan overlay",    1200, 700);
    TCanvas* cDeriv = new TCanvas("cDeriv", "Derivative overlay (norm)", 1200, 700);
    for (TCanvas* cv : {cScan, cDeriv}) {
        cv->SetGrid(); cv->SetTicks(1,1);
        cv->SetLeftMargin(PAD_LEFT); cv->SetRightMargin(PAD_RIGHT);
        cv->SetTopMargin(PAD_TOP);   cv->SetBottomMargin(PAD_BOTTOM);
    }

    // ── Legend ───────────────────────────────────────────────
    double legH      = std::min(0.07 * nf + 0.04, 0.55);
    double legTextSz = std::max(0.024, 0.040 - 0.002 * nf);
    double legX2 = 1.0 - PAD_RIGHT - 0.01;
    double legX1 = legX2 - 0.30;
    double legY2 = 1.0 - PAD_TOP  - 0.02;
    double legY1 = legY2 - legH;
    TLegend* legS = new TLegend(legX1, legY1, legX2, legY2);
    TLegend* legD = new TLegend(legX1, legY1, legX2, legY2);
    for (auto* leg : {legS, legD}) {
        leg->SetBorderSize(1); leg->SetFillStyle(1001);
        leg->SetFillColor(0);  leg->SetTextSize(legTextSz);
    }

    std::vector<TGraph*> gScans, gDerivs;
    bool firstS = true, firstD = true;

    for (int fi = 0; fi < nf; ++fi) {
        auto& info = files[fi];
        TFile* f = TFile::Open(info.path.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "  [SKIP] " << info.path << "\n"; continue;
        }

        int col = fileColor(fi);
        std::string label = Form("V_{bias}=%.0fV, filter=%.0f",
                                 info.vbias, info.filter);

        // ── Scan (log Y scale) ───────────────────────────────
        TGraph* gS = loadScan(f);
        if (gS) {
            gS->SetTitle(";Threshold (mV);Counts");
            gS->SetLineColor(col); gS->SetLineWidth(2);
            cScan->cd();
            gS->Draw(firstS ? "AL" : "L same");
            if (firstS) { axisStyle(gS); cScan->SetLogy(1); firstS = false; }
            legS->AddEntry(gS, label.c_str(), "l");
            gScans.push_back(gS);
        }

        // ── Normalised derivative ────────────────────────────
        TGraph* gD = loadDerivNorm(f);
        if (gD) {
            gD->SetTitle(";Threshold (mV);-dN/dV  (normalised to peak)");
            gD->SetLineColor(col); gD->SetLineWidth(2);
            cDeriv->cd();
            gD->Draw(firstD ? "AL" : "L same");
            if (firstD) { axisStyle(gD); firstD = false; }
            legD->AddEntry(gD, label.c_str(), "l");
            gDerivs.push_back(gD);
        }
        // Keep file open (graphs reference tree data)
    }

    // ── Adaptive X range: trim trailing zeros ────────────────
    double xmax_data = trimmedXMax(gScans);
    double xRange    = xmax_data - MIN_THR_DISPLAY;
    double xlo = MIN_THR_DISPLAY - xRange * 0.02;
    double xhi = xmax_data       + xRange * 0.03;

    for (auto* g : gScans)  if (g) { g->GetXaxis()->SetLimits(xlo, xhi); }
    for (auto* g : gDerivs) if (g) { g->GetXaxis()->SetLimits(xlo, xhi); }

    // ── Deriv Y range: [-0.15, 1.15] after normalisation ────
    // Small negative margin shows the undershoot, 1.15 gives headroom above 1
    for (auto* g : gDerivs)
        if (g) g->GetYaxis()->SetRangeUser(-0.20, 1.20);

    // ── Legends and save ─────────────────────────────────────
    cScan ->cd(); legS->Draw();
    cDeriv->cd(); legD->Draw();

    cScan ->Update(); cScan ->Modified();
    cDeriv->Update(); cDeriv->Modified();

    cScan ->SaveAs("all_scans.png");
    cDeriv->SaveAs("all_derivs.png");
    std::cout << "Saved: all_scans.png,  all_derivs.png\n";
}
