/**
 * laser_profile.cpp  —  Analisi profilo fascio laser per SiPM pos-scan
 *
 * Usa OutputManager.h, Config.h, InputHelpers.h dal progetto.
 *
 * LOGICA DEFAULT (modificabile a runtime):
 *   CROSSING → P_det, mappe 2D, profilo radiale, slice X/Y
 *   ACCEPTED → sigma(Δt), mu(Δt), delta_mu map (qualità timing)
 *
 * PARAMETRI ESTRATTI:
 *   Da CROSSING: P_det(x,y), sigma_beam (fit gaussiano), profilo radiale
 *   Da ACCEPTED: sigma_t (risoluzione temporale), Δμ(x,y), mu_centro
 *
 * BUG FIX rispetto alla versione precedente:
 *   [1] P0 fit gaussiano: limite superiore 100% (non 105%)
 *   [2] Slice Y: fit solo per y>=0 (coda, non gaussiana intera) — usa half-gauss
 *   [3] Profilo radiale: fit con P0 libero ma bounded [0,100]
 *   [4] p_det letto come 0-1 dal tree, moltiplicato x100 correttamente
 *   [5] n_cross vs n_laser usato per errore binomiale (non n_acc)
 *   [6] Delta_mu: usa solo punti con P_det_crossing > soglia (default 5%)
 *   [7] Outlier Δμ: taglio |Δμ| > 0.5 ns (fit falliti)
 *   [8] Centro fisso a (cx,cy) da input, non calcolato dalla media pesata
 *
 * Compile:  .L laser_profile.cpp+
 * Run:      laser_profile()
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/InputHelpers.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TView.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// ════════════════════════════════════════════════════════════════════════════
//  Struttura punto griglia
// ════════════════════════════════════════════════════════════════════════════
struct MP {
    double x, y;
    double mu, mu_err;
    double sigma, sigma_err;
    double p_det, p_det_err;   // %, 0-100
    long   n_acc, n_cross, n_laser;
    int    vbias;
    double frac;
    bool   ok;    // true se p_det > pdet_cut
};

// ════════════════════════════════════════════════════════════════════════════
//  Carica file ROOT (crossing o accepted)
// ════════════════════════════════════════════════════════════════════════════
static std::vector<MP> loadMap(const std::string& path, int vbias, double frac,
                                double pdet_cut_pct = 5.0)
{
    std::vector<MP> pts;
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "  [WARN] cannot open: " << path << "\n";
        return pts;
    }
    TTree* t = (TTree*)f->Get("map");
    if (!t) { f->Close(); delete f; return pts; }

    Double_t x=0, y=0, mu=0, mu_err=0, sigma=0, sigma_err=0, p_det=0;
    Long64_t n_acc=0, n_cross=0, n_laser=0;

    t->SetBranchAddress("x",     &x);
    t->SetBranchAddress("y",     &y);
    t->SetBranchAddress("mu",    &mu);
    t->SetBranchAddress("sigma", &sigma);
    t->SetBranchAddress("p_det", &p_det);

    bool hme = (t->GetBranch("mu_err")    != nullptr);
    bool hse = (t->GetBranch("sigma_err") != nullptr);
    bool hna = (t->GetBranch("n_acc")     != nullptr);
    bool hnc = (t->GetBranch("n_cross") != nullptr ||
                t->GetBranch("n_crossing") != nullptr);
    bool hnl = (t->GetBranch("n_laser")   != nullptr);

    if (hme) t->SetBranchAddress("mu_err",    &mu_err);
    if (hse) t->SetBranchAddress("sigma_err", &sigma_err);
    if (hna) t->SetBranchAddress("n_acc",     &n_acc);
    if (hnc) {
        if (t->GetBranch("n_crossing"))
            t->SetBranchAddress("n_crossing", &n_cross);
        else
            t->SetBranchAddress("n_cross", &n_cross);
    }
    if (hnl) t->SetBranchAddress("n_laser", &n_laser);

    Long64_t N = t->GetEntries();
    pts.reserve((size_t)N);

    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);

        // p_det salvato come frazione 0-1 nel tree
        double p_pct = p_det * 100.0;
        // denominatore per errore binomiale: n_cross se disponibile, else n_laser
        double den = hnc ? (double)n_cross :
                     hnl ? (double)n_laser : 1000.0;
        double p_frac = p_pct / 100.0;
        double err_pct = (den > 0) ? 100.0 * std::sqrt(p_frac * (1.0 - p_frac) / den) : 2.0;

        MP mp;
        mp.x         = x;
        mp.y         = y;
        mp.mu        = mu;
        mp.mu_err    = hme ? mu_err    : 0.0;
        mp.sigma     = sigma;
        mp.sigma_err = hse ? sigma_err : 0.0;
        mp.p_det     = p_pct;
        mp.p_det_err = err_pct;
        mp.n_acc     = hna ? (long)n_acc   : 0L;
        mp.n_cross   = hnc ? (long)n_cross : 0L;
        mp.n_laser   = hnl ? (long)n_laser : 0L;
        mp.vbias     = vbias;
        mp.frac      = frac;
        mp.ok        = (p_pct >= pdet_cut_pct);
        pts.push_back(mp);
    }
    f->Close(); delete f;
    std::cout << "  Loaded " << pts.size() << " pts  [" << path << "]\n";
    return pts;
}

// ════════════════════════════════════════════════════════════════════════════
//  Scopri cartelle dati disponibili nella base dati
// ════════════════════════════════════════════════════════════════════════════
static std::string chooseDataDir()
{
    // Elenca sottocartelle di DATA_BASE
    std::vector<std::string> dirs;
    void* dh = gSystem->OpenDirectory(DATA_BASE.c_str());
    if (dh) {
        const char* entry;
        while ((entry = gSystem->GetDirEntry(dh)) != nullptr) {
            std::string name(entry);
            if (name == "." || name == "..") continue;
            std::string full = DATA_BASE + "/" + name;
            // controlla che sia una directory
            Long_t id=0, flags=0, modtime=0; Long64_t size=0;
            if (gSystem->GetPathInfo(full.c_str(), &id, &size, &flags, &modtime) == 0
                && (flags & 2)) {   // bit 1 = directory
                dirs.push_back(full);
            }
        }
        gSystem->FreeDirectory(dh);
    }
    std::sort(dirs.begin(), dirs.end());

    if (dirs.empty()) {
        std::cerr << "  [WARN] Nessuna sottocartella in " << DATA_BASE
                  << ". Uso DATA_DIR default.\n";
        return DATA_DIR;
    }

    std::cout << "\n  Cartelle dati disponibili:\n";
    for (size_t i = 0; i < dirs.size(); ++i) {
        // conta i file map_results in questa dir
        int nmaps = 0;
        void* dd = gSystem->OpenDirectory(dirs[i].c_str());
        if (dd) {
            const char* e2;
            while ((e2 = gSystem->GetDirEntry(dd)) != nullptr) {
                std::string fn(e2);
                if (fn.find("map_results_") == 0 &&
                    fn.size() > 5 &&
                    fn.substr(fn.size()-5) == ".root") ++nmaps;
            }
            gSystem->FreeDirectory(dd);
        }
        std::cout << "  [" << i << "]  " << dirs[i]
                  << "  (" << nmaps << " map files)\n";
    }

    // Default: cartella con più file map_results (tipicamente giacomo_serpentina)
    size_t best = 0; int bestN = -1;
    for (size_t i = 0; i < dirs.size(); ++i) {
        void* dd = gSystem->OpenDirectory(dirs[i].c_str());
        int nm = 0;
        if (dd) {
            const char* e2;
            while ((e2 = gSystem->GetDirEntry(dd)) != nullptr) {
                std::string fn(e2);
                if (fn.find("map_results_") == 0 &&
                    fn.size() > 5 &&
                    fn.substr(fn.size()-5) == ".root") ++nm;
            }
            gSystem->FreeDirectory(dd);
        }
        if (nm > bestN) { bestN = nm; best = i; }
    }

    std::cout << "\n  Cartella (INVIO=" << best << "): ";
    std::string sel; std::getline(std::cin, sel);
    if (sel.empty()) return dirs[best];
    try {
        size_t idx = std::stoul(sel);
        if (idx < dirs.size()) return dirs[idx];
    } catch (...) {}
    std::cerr << "  [WARN] Selezione non valida, uso default.\n";
    return dirs[best];
}

// ════════════════════════════════════════════════════════════════════════════
//  Scopri file map_results in una directory
// ════════════════════════════════════════════════════════════════════════════
struct MapFileInfo {
    std::string path;
    int    vbias;
    double frac;
    std::string mode;   // "crossing" o "accepted"
};

static std::vector<MapFileInfo> findMapFiles(const std::string& dir,
                                              const std::string& mode)
{
    std::vector<MapFileInfo> files;
    const std::string suf = "_pdet_" + mode + ".root";
    const std::string pre = "map_results_vbias";

    void* dh = gSystem->OpenDirectory(dir.c_str());
    if (!dh) return files;
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dh)) != nullptr) {
        std::string fn(entry);
        if (fn.size() < pre.size() + suf.size()) continue;
        if (fn.substr(0, pre.size()) != pre) continue;
        if (fn.substr(fn.size() - suf.size()) != suf) continue;

        // parse vbias e frac
        std::string mid = fn.substr(pre.size(), fn.size() - pre.size() - suf.size());
        // mid = "52_let0.60" oppure "52_let0.60_loose"
        size_t pos = mid.find("_let");
        if (pos == std::string::npos) continue;
        MapFileInfo fi;
        try {
            fi.vbias = std::stoi(mid.substr(0, pos));
            // frac: tutto dopo "_let" fino a '_' o fine stringa
            std::string fstr = mid.substr(pos + 4);
            size_t us = fstr.find('_');
            if (us != std::string::npos) fstr = fstr.substr(0, us);
            fi.frac = std::stod(fstr);
        } catch (...) { continue; }
        fi.path = dir + "/" + fn;
        fi.mode = mode;
        files.push_back(fi);
    }
    gSystem->FreeDirectory(dh);
    std::sort(files.begin(), files.end(), [](const MapFileInfo& a, const MapFileInfo& b){
        return a.vbias < b.vbias || (a.vbias == b.vbias && a.frac < b.frac);
    });
    return files;
}

// ════════════════════════════════════════════════════════════════════════════
//  Profilo radiale binned (media punti in annelli)
// ════════════════════════════════════════════════════════════════════════════
static TGraphErrors* radProfile(const std::vector<MP>& pts,
                                 double cx, double cy,
                                 const char* name,
                                 double pdet_cut = 5.0)
{
    const double DR = 8.0;
    std::map<int, std::vector<double>> bins;
    for (const auto& p : pts) {
        if (p.p_det < pdet_cut) continue;
        double r = std::hypot(p.x - cx, p.y - cy);
        bins[(int)(r / DR)].push_back(p.p_det);
    }
    std::vector<double> vR, vP, vRe, vPe;
    for (auto& [bin, v] : bins) {
        double rm = (bin + 0.5) * DR;
        double pm = 0;
        for (double d : v) pm += d;
        pm /= (double)v.size();
        double pv = 0;
        for (double d : v) pv += (d - pm) * (d - pm);
        double se = (v.size() > 1) ?
            std::sqrt(pv / (v.size() - 1)) / std::sqrt((double)v.size()) : 2.0;
        vR.push_back(rm); vP.push_back(pm);
        vRe.push_back(DR * 0.5); vPe.push_back(se);
    }
    if (vR.empty()) return nullptr;
    TGraphErrors* gr = new TGraphErrors((int)vR.size(),
        vR.data(), vP.data(), vRe.data(), vPe.data());
    gr->SetName(name);
    return gr;
}

// ════════════════════════════════════════════════════════════════════════════
//  Mappa sigma_t(x,y) da dati accepted
// ════════════════════════════════════════════════════════════════════════════
static void drawSigmaMap(const std::vector<MP>& pts, int vb, double fr,
                          OutCtx& ctx)
{
    std::set<double> xs, ys;
    for (const auto& p : pts) { xs.insert(p.x); ys.insert(p.y); }
    int nx = (int)xs.size(), ny = (int)ys.size();
    if (nx < 2 || ny < 2) return;
    double xmin=*xs.begin(), xmax=*xs.rbegin();
    double ymin=*ys.begin(), ymax=*ys.rbegin();
    double dx=(*std::next(xs.begin())-xmin), dy=(*std::next(ys.begin())-ymin);

    // sigma in ns → converti in ps per leggibilità
    TH2D* h = new TH2D(Form("hSig_%d_%.2f", vb, fr),
        Form("#sigma_{t}(x,y)   V_{bias}=%dV   LET=%.2f p.e.;"
             "#Deltax (mm);#Deltay (mm);#sigma_{t} (ns)", vb, fr),
        nx, xmin-dx/2, xmax+dx/2,
        ny, ymin-dy/2, ymax+dy/2);
    h->SetDirectory(nullptr);

    double sig_sum = 0, sig_n = 0;
    double sig_min = 1e9, sig_max = 0;
    for (const auto& p : pts) {
        if (!p.ok || p.sigma <= 0) continue;
        h->Fill(p.x, p.y, p.sigma);
        sig_sum += p.sigma; ++sig_n;
        sig_min = std::min(sig_min, p.sigma);
        sig_max = std::max(sig_max, p.sigma);
    }
    double sig_mean = (sig_n > 0) ? sig_sum / sig_n : 0;

    h->GetZaxis()->SetRangeUser(sig_min * 0.95, sig_max * 1.05);

    TCanvas* c = new TCanvas(Form("cSig_%d", vb), "", 820, 700);
    c->SetRightMargin(0.17); c->SetLeftMargin(0.10);
    c->SetBottomMargin(0.10); c->SetTopMargin(0.09);
    gStyle->SetPalette(kBird);
    h->Draw("COLZ");

    for (const auto& p : pts) {
        if (!p.ok || p.sigma <= 0) continue;
        TLatex* lt = new TLatex(p.x, p.y, Form("%.3f", p.sigma));
        lt->SetTextSize(0.018); lt->SetTextAlign(22); lt->Draw("same");
    }

    TPaveText* info = new TPaveText(0.01, 0.01, 0.75, 0.07, "NDC");
    info->SetBorderSize(1); info->SetFillColor(0); info->SetTextSize(0.026);
    info->AddText(Form("#bar{#sigma}=%.3f ns   #sigma_{min}=%.3f ns   #sigma_{max}=%.3f ns   "
                       "Solo P_{det,cross}>5%%", sig_mean, sig_min, sig_max));
    info->Draw();
    c->Update(); c->Modified();
    ctx.savePNG(c, Form("sigma_map_v%d_let%.2f.png", vb, fr));
}

// ════════════════════════════════════════════════════════════════════════════
//  Mappa P_det 2D (da crossing)
// ════════════════════════════════════════════════════════════════════════════
static void drawPdetMap(const std::vector<MP>& pts, int vb, double fr,
                         double cx, double cy, OutCtx& ctx)
{
    std::set<double> xs, ys;
    for (const auto& p : pts) { xs.insert(p.x); ys.insert(p.y); }
    int nx = (int)xs.size(), ny = (int)ys.size();
    if (nx < 2 || ny < 2) return;
    double xmin=*xs.begin(), xmax=*xs.rbegin();
    double ymin=*ys.begin(), ymax=*ys.rbegin();
    double dx=(*std::next(xs.begin())-xmin), dy=(*std::next(ys.begin())-ymin);

    TH2D* h = new TH2D(Form("hPd_%d_%.2f", vb, fr),
        Form("P_{det}   V_{bias}=%dV   LET=%.2f p.e.;"
             "#Deltax (mm);#Deltay (mm)", vb, fr),
        nx, xmin-dx/2, xmax+dx/2,
        ny, ymin-dy/2, ymax+dy/2);
    h->SetDirectory(nullptr);
    h->GetZaxis()->SetTitle("P_{det} (%)");
    h->GetZaxis()->SetRangeUser(0, 100);

    for (const auto& p : pts) h->Fill(p.x, p.y, p.p_det);

    TCanvas* c = new TCanvas(Form("cPd_%d", vb), "", 820, 700);
    c->SetRightMargin(0.16); c->SetLeftMargin(0.10);
    c->SetBottomMargin(0.10); c->SetTopMargin(0.09);
    gStyle->SetPalette(kRainBow);
    h->Draw("COLZ");

    for (const auto& p : pts) {
        TLatex* lt = new TLatex(p.x, p.y, Form("%.0f", p.p_det));
        lt->SetTextSize(0.020); lt->SetTextAlign(22); lt->Draw("same");
    }
    TMarker* mk = new TMarker(cx, cy, 5);
    mk->SetMarkerColor(kWhite); mk->SetMarkerSize(2.5); mk->Draw("same");

    c->Update(); c->Modified();
    ctx.savePNG(c, Form("laser_map2D_v%d_let%.2f.png", vb, fr));
}

// ════════════════════════════════════════════════════════════════════════════
//  Slice X e Y (gaussiana per X, mezza-gaussiana per Y — centro sul bordo)
// ════════════════════════════════════════════════════════════════════════════
static void drawSlices(const std::vector<MP>& pts, int vb, double fr,
                        double cx, double cy, OutCtx& ctx)
{
    // ── Slice X (riga y più vicina a cy) ─────────────────────────────────
    double dyb = 1e9;
    for (const auto& p : pts)
        if (std::abs(p.y - cy) < dyb) dyb = std::abs(p.y - cy);
    std::vector<double> vX, vPx, vPxe;
    for (const auto& p : pts)
        if (std::abs(p.y - cy) <= dyb + 0.5) {
            vX.push_back(p.x - cx);
            vPx.push_back(p.p_det);
            vPxe.push_back(p.p_det_err);
        }
    { // sort per x
        std::vector<size_t> idx(vX.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&vX](size_t a, size_t b){ return vX[a] < vX[b]; });
        std::vector<double> tx, tp, te;
        for (size_t k : idx) { tx.push_back(vX[k]); tp.push_back(vPx[k]); te.push_back(vPxe[k]); }
        vX = tx; vPx = tp; vPxe = te;
    }

    // ── Slice Y (colonna x più vicina a cx) ──────────────────────────────
    double dxb = 1e9;
    for (const auto& p : pts)
        if (std::abs(p.x - cx) < dxb) dxb = std::abs(p.x - cx);
    std::vector<double> vY, vPy, vPye;
    for (const auto& p : pts)
        if (std::abs(p.x - cx) <= dxb + 0.5) {
            vY.push_back(p.y - cy);
            vPy.push_back(p.p_det);
            vPye.push_back(p.p_det_err);
        }
    {
        std::vector<size_t> idx(vY.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&vY](size_t a, size_t b){ return vY[a] < vY[b]; });
        std::vector<double> ty, tp, te;
        for (size_t k : idx) { ty.push_back(vY[k]); tp.push_back(vPy[k]); te.push_back(vPye[k]); }
        vY = ty; vPy = tp; vPye = te;
    }

    TCanvas* cSl = new TCanvas(Form("cSl_%d", vb), "", 1100, 500);
    cSl->Divide(2, 1);

    // ── Fit X: gaussiana simmetrica, P0 bounded [0, 100] ─────────────────
    TF1* fGx = new TF1(Form("fGx_%d", vb),
                        "[0]*exp(-x*x/(2*[1]*[1]))", -80, 80);
    fGx->SetParameters(95, 25);
    fGx->SetParLimits(0,  0.0, 100.0);   // BUG FIX [1]: P0 max=100%
    fGx->SetParLimits(1,  5.0,  80.0);
    fGx->SetLineColor(kRed+1); fGx->SetLineWidth(2); fGx->SetLineStyle(2);

    // ── Fit Y: mezza-gaussiana (solo lato destro, centro a bordo griglia)
    // Usiamo f(y) = P0 * exp(-y²/(2σ²)) con y>=0
    // Il fit ha senso solo se il centro è fuori griglia (y_c = 0 = bordo)
    // BUG FIX [2]: non forziamo gaussiana intera su dati troncati
    TF1* fGy = new TF1(Form("fGy_%d", vb),
                        "[0]*exp(-x*x/(2*[1]*[1]))", 0, 120);
    fGy->SetParameters(98, 50);
    fGy->SetParLimits(0,  0.0, 100.0);   // P0 bounded
    fGy->SetParLimits(1, 10.0, 150.0);   // sigma ampio: vediamo solo la coda
    fGy->SetLineColor(kRed+1); fGy->SetLineWidth(2); fGy->SetLineStyle(2);

    cSl->cd(1); gPad->SetGrid();
    if (!vX.empty()) {
        TGraphErrors* grX = new TGraphErrors((int)vX.size(),
            vX.data(), vPx.data(), nullptr, vPxe.data());
        grX->SetTitle(Form("Slice X   V_{bias}=%dV;x - x_{c} (mm);P_{det} (%%)", vb));
        grX->SetMarkerStyle(20); grX->SetMarkerColor(kBlue+1); grX->SetLineColor(kBlue+1);
        grX->GetYaxis()->SetRangeUser(0, 110);
        grX->Draw("AP");
        grX->Fit(fGx, "RQ", "", -80, 80);
        TLatex* lx = new TLatex(0.15, 0.84,
            Form("#sigma_{X}=%.1f mm   P_{0}=%.1f%%",
                 std::abs(fGx->GetParameter(1)), fGx->GetParameter(0)));
        lx->SetNDC(); lx->SetTextSize(0.044); lx->Draw();
    }

    cSl->cd(2); gPad->SetGrid();
    if (!vY.empty()) {
        TGraphErrors* grY = new TGraphErrors((int)vY.size(),
            vY.data(), vPy.data(), nullptr, vPye.data());
        grY->SetTitle(Form("Slice Y   V_{bias}=%dV   [semi-gauss, centro a bordo];"
                           "y - y_{c} (mm);P_{det} (%%)", vb));
        grY->SetMarkerStyle(21); grY->SetMarkerColor(kGreen+2); grY->SetLineColor(kGreen+2);
        grY->GetYaxis()->SetRangeUser(0, 110);
        grY->Draw("AP");
        grY->Fit(fGy, "RQ", "", 0, 110);
        TLatex* ly = new TLatex(0.15, 0.84,
            Form("#sigma_{Y}=%.1f mm   P_{0}=%.1f%%   [solo coda]",
                 std::abs(fGy->GetParameter(1)), fGy->GetParameter(0)));
        ly->SetNDC(); ly->SetTextSize(0.040); ly->Draw();
        // Nota explicativa
        TLatex* note = new TLatex(0.15, 0.76,
            "#font[12]{Il centro del fascio e' sul bordo inferiore della griglia}");
        note->SetNDC(); note->SetTextSize(0.032); note->SetTextColor(kGray+2); note->Draw();
    }
    cSl->Update(); cSl->Modified();
    ctx.savePNG(cSl, Form("laser_slice_v%d_let%.2f.png", vb, fr));
}

// ════════════════════════════════════════════════════════════════════════════
//  Mappa Delta-mu (timing shift rispetto al centro)
// ════════════════════════════════════════════════════════════════════════════
static void drawDeltaMu(const std::vector<MP>& pts_acc,
                         const std::vector<MP>& pts_cross,
                         int vb, double fr,
                         double cx, double cy,
                         double pdet_cut,
                         OutCtx& ctx)
{
    // mu_centro: media pesata per p_det dei punti r < 15 mm (da crossing)
    double mu_c = 0, wsum = 0;
    for (const auto& p : pts_cross) {
        double r = std::hypot(p.x - cx, p.y - cy);
        if (r < 15.0 && p.p_det >= pdet_cut) {
            mu_c  += p.mu * p.p_det;
            wsum  += p.p_det;
        }
    }
    if (wsum > 0) {
        mu_c /= wsum;
    } else {
        // fallback: punto più vicino al centro (da pts_acc o pts_cross)
        double rmin = 1e9;
        for (const auto& p : pts_acc) {
            double r = std::hypot(p.x - cx, p.y - cy);
            if (r < rmin) { rmin = r; mu_c = p.mu; }
        }
    }
    std::cout << "  [Delta-mu] Vbias=" << vb << "V  mu_centro=" << std::fixed
              << std::setprecision(3) << mu_c << " ns\n";

    // Usa pts_acc per il timing (più pulito), pts_cross per la maschera statistica
    // Costruisci mappa p_det_crossing per lookup veloce
    std::map<std::pair<int,int>, double> pdet_cross_map;
    for (const auto& p : pts_cross) {
        pdet_cross_map[{(int)p.x, (int)p.y}] = p.p_det;
    }

    std::set<double> xs, ys;
    for (const auto& p : pts_acc) { xs.insert(p.x); ys.insert(p.y); }
    int nx = (int)xs.size(), ny = (int)ys.size();
    if (nx < 2 || ny < 2) return;
    double xmin=*xs.begin(), xmax=*xs.rbegin();
    double ymin=*ys.begin(), ymax=*ys.rbegin();
    double dx=(*std::next(xs.begin())-xmin), dy=(*std::next(ys.begin())-ymin);

    TH2D* hMu = new TH2D(Form("hMu_%d_%.2f", vb, fr),
        Form("#Delta#mu(x,y) = #mu(x,y)-#mu_{centro}"
             "   V_{bias}=%dV   LET=%.2f p.e.;"
             "#Deltax (mm);#Deltay (mm);#Delta#mu (ns)", vb, fr),
        nx, xmin-dx/2, xmax+dx/2,
        ny, ymin-dy/2, ymax+dy/2);
    hMu->SetDirectory(nullptr);

    int n_used = 0, n_skipped = 0;
    for (const auto& p : pts_acc) {
        // Maschera: usa p_det crossing per decidere se c'è abbastanza statistica
        double pdc = 0;
        auto it = pdet_cross_map.find({(int)p.x, (int)p.y});
        if (it != pdet_cross_map.end()) pdc = it->second;
        else pdc = p.p_det;   // fallback se non c'è crossing

        if (pdc < pdet_cut) { ++n_skipped; continue; }

        double dmu = p.mu - mu_c;
        // BUG FIX [7]: taglia outlier (fit falliti danno dmu >> 0.5 ns)
        if (std::abs(dmu) > 0.5 || p.n_acc < 50) { ++n_skipped; continue; }

        hMu->Fill(p.x, p.y, dmu);
        ++n_used;
    }
    std::cout << "  [Delta-mu] Usati=" << n_used << "  Scartati=" << n_skipped
              << "  (taglio: P_det_cross<" << pdet_cut << "% o |Δμ|>0.5ns)\n";

    // ── Canvas 2D ────────────────────────────────────────────────────────
    TCanvas* cMu2 = new TCanvas(Form("cMu2_%d", vb), "", 860, 700);
    cMu2->SetRightMargin(0.17); cMu2->SetLeftMargin(0.10);
    cMu2->SetBottomMargin(0.10); cMu2->SetTopMargin(0.09);
    gStyle->SetPalette(kTemperatureMap);
    hMu->Draw("COLZ");

    for (const auto& p : pts_acc) {
        double pdc = 0;
        auto it = pdet_cross_map.find({(int)p.x, (int)p.y});
        if (it != pdet_cross_map.end()) pdc = it->second;
        else pdc = p.p_det;
        if (pdc < pdet_cut) continue;
        double dmu = p.mu - mu_c;
        if (std::abs(dmu) > 0.5) continue;
        TLatex* lt = new TLatex(p.x, p.y, Form("%.3f", dmu));
        lt->SetTextSize(0.019); lt->SetTextAlign(22); lt->Draw("same");
    }
    TMarker* mk = new TMarker(cx, cy, 5);
    mk->SetMarkerColor(kBlack); mk->SetMarkerSize(2); mk->Draw("same");

    TPaveText* ptm = new TPaveText(0.01, 0.01, 0.75, 0.07, "NDC");
    ptm->SetBorderSize(1); ptm->SetFillColor(0); ptm->SetTextSize(0.024);
    ptm->AddText(Form("#mu_{centro}=%.3f ns   Timing da ACCEPTED, maschera da CROSSING"
                      "   P_{det,cross}>%.0f%%   |#Delta#mu|<0.5 ns", mu_c, pdet_cut));
    ptm->Draw();
    cMu2->Update(); cMu2->Modified();
    ctx.savePNG(cMu2, Form("delta_mu_map2D_v%d_let%.2f.png", vb, fr));

    // ── Canvas 3D ────────────────────────────────────────────────────────
    TGraph2D* g3 = new TGraph2D();
    g3->SetName(Form("g3mu_%d", vb));
    int np = 0;
    for (const auto& p : pts_acc) {
        double pdc = 0;
        auto it = pdet_cross_map.find({(int)p.x, (int)p.y});
        if (it != pdet_cross_map.end()) pdc = it->second;
        else pdc = p.p_det;
        if (pdc < pdet_cut) continue;
        double dmu = p.mu - mu_c;
        if (std::abs(dmu) > 0.5) continue;
        g3->SetPoint(np++, p.x - cx, p.y - cy, dmu);
    }
    g3->SetTitle(Form("#Delta#mu(x,y)   V_{bias}=%dV   LET=%.2f p.e.;"
                      "#Deltax (mm);#Deltay (mm);#Delta#mu (ns)", vb, fr));

    TCanvas* cMu3 = new TCanvas(Form("cMu3_%d", vb), "", 900, 750);
    cMu3->SetLeftMargin(0.05); cMu3->SetRightMargin(0.05);
    cMu3->SetBottomMargin(0.05); cMu3->SetTopMargin(0.08);
    gStyle->SetPalette(kTemperatureMap);
    g3->Draw("surf1");
    gPad->Update();
    if (gPad->GetView()) {
        gPad->GetView()->RotateView(30, 20);
        gPad->Modified(); gPad->Update();
    }
    TPaveText* pt3 = new TPaveText(0.05, 0.01, 0.95, 0.06, "NDC");
    pt3->SetBorderSize(0); pt3->SetFillColor(0); pt3->SetTextSize(0.028);
    pt3->AddText("Asimmetria sin-dx = laser inclinato.   Sella simmetrica = ritardo geometrico d/c.");
    pt3->Draw();
    cMu3->Update(); cMu3->Modified();
    ctx.savePNG(cMu3, Form("delta_mu_3D_v%d_let%.2f.png", vb, fr));
}

// ════════════════════════════════════════════════════════════════════════════
//  Grafico sigma_t vs Vbias (media sulla griglia, da accepted)
// ════════════════════════════════════════════════════════════════════════════
static void drawSigmaVsVbias(
    const std::vector<std::vector<MP>>& allAcc,
    const std::vector<int>& vbiases,
    const std::vector<double>& fracs,
    OutCtx& ctx)
{
    if (allAcc.empty()) return;

    // raggruppa per frac
    std::map<double, std::vector<std::pair<int,double>>> byFrac; // frac -> [(vbias, sigma_mean)]
    for (size_t i = 0; i < allAcc.size(); ++i) {
        double sum = 0; int n = 0;
        for (const auto& p : allAcc[i]) {
            if (p.ok && p.sigma > 0 && p.sigma < 0.5) {
                sum += p.sigma; ++n;
            }
        }
        if (n > 0) byFrac[fracs[i]].push_back({vbiases[i], sum / n});
    }
    if (byFrac.empty()) return;

    TCanvas* cSV = new TCanvas("cSigVbias", "", 850, 600);
    cSV->SetGrid(); cSV->SetLeftMargin(0.12); cSV->SetBottomMargin(0.12);
    TLegend* leg = new TLegend(0.55, 0.60, 0.93, 0.88);
    leg->SetBorderSize(1); leg->SetFillColor(0); leg->SetTextSize(0.030);

    const int cols[] = {kBlue+1, kRed+1, kGreen+2, kOrange+7, kMagenta+1};
    bool first = true; int ci = 0;
    for (auto& [frac, pts] : byFrac) {
        std::sort(pts.begin(), pts.end());
        int n = (int)pts.size();
        std::vector<double> vb(n), vs(n);
        for (int k = 0; k < n; ++k) { vb[k] = pts[k].first; vs[k] = pts[k].second * 1000.; }
        TGraph* gr = new TGraph(n, vb.data(), vs.data());
        int col = cols[ci++ % 5];
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.3);
        gr->SetMarkerColor(col); gr->SetLineColor(col); gr->SetLineWidth(2);
        gr->SetTitle(";V_{bias} (V);#bar{#sigma}_{t} (ps)");
        if (first) {
            gr->GetYaxis()->SetRangeUser(0, 250);
            gr->Draw("ALP"); first = false;
        } else gr->Draw("LP same");
        leg->AddEntry(gr, Form("LET=%.2f p.e.", frac), "lp");
    }
    leg->Draw();

    TPaveText* inf = new TPaveText(0.12, 0.82, 0.52, 0.92, "NDC");
    inf->SetBorderSize(1); inf->SetFillColor(0); inf->SetTextSize(0.028);
    inf->AddText("Media #sigma_{t} su griglia (accepted, P_{det,cross}>5%)");
    inf->Draw();
    cSV->Update(); cSV->Modified();
    ctx.savePNG(cSV, "sigma_vs_vbias.png");
}

// ════════════════════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════════════════════
void laser_profile()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    gStyle->SetNumberContours(99);

    std::cout << "\n+============================================================+\n"
              << "|  LASER PROFILE — SiPM position scan analysis              |\n"
              << "|  CROSSING → P_det, profilo fascio                         |\n"
              << "|  ACCEPTED → sigma_t, delta_mu (timing)                    |\n"
              << "+============================================================+\n\n";

    // ── 1. Scegli cartella dati ──────────────────────────────────────────
    std::string dataDir = chooseDataDir();
    std::cout << "  Usando: " << dataDir << "\n";

    // ── 2. Scopri file crossing e accepted ───────────────────────────────
    auto filesCross = findMapFiles(dataDir, "crossing");
    auto filesAcc   = findMapFiles(dataDir, "accepted");

    if (filesCross.empty() && filesAcc.empty()) {
        std::cerr << "  Nessun file map_results trovato in " << dataDir << "\n";
        return;
    }

    // Mostra file disponibili
    std::set<std::pair<int,double>> avail_vb_frac;
    for (const auto& f : filesCross) avail_vb_frac.insert({f.vbias, f.frac});
    for (const auto& f : filesAcc)   avail_vb_frac.insert({f.vbias, f.frac});

    std::vector<std::pair<int,double>> combos(avail_vb_frac.begin(), avail_vb_frac.end());
    std::cout << "\n  Combinazioni (Vbias, LET) disponibili:\n";
    for (size_t i = 0; i < combos.size(); ++i)
        std::cout << "  [" << i << "]  Vbias=" << combos[i].first
                  << "V  LET=" << std::fixed << std::setprecision(2)
                  << combos[i].second << " pe\n";

    std::cout << "\n  Seleziona indici (INVIO=tutti): ";
    std::string sel; std::getline(std::cin, sel);
    std::vector<size_t> chosen;
    if (sel.empty()) {
        for (size_t i = 0; i < combos.size(); ++i) chosen.push_back(i);
    } else {
        std::istringstream ss(sel);
        size_t v;
        while (ss >> v) if (v < combos.size()) chosen.push_back(v);
    }
    if (chosen.empty()) { std::cerr << "  Nessuna selezione.\n"; return; }

    // ── 3. Parametri analisi ─────────────────────────────────────────────
    std::cout << "\n  Centro fascio x_c (INVIO=90): ";
    std::string sx; std::getline(std::cin, sx);
    double cx = sx.empty() ? 90.0 : std::stod(sx);

    std::cout << "  Centro fascio y_c (INVIO=0): ";
    std::string sy; std::getline(std::cin, sy);
    double cy_val = sy.empty() ? 0.0 : std::stod(sy);

    std::cout << "  Soglia P_det per maschera statistica % (INVIO=5): ";
    std::string sp; std::getline(std::cin, sp);
    double pdet_cut = sp.empty() ? 5.0 : std::stod(sp);

    // ── 4. Output dirs ───────────────────────────────────────────────────
    OutCtx ctx = createOutputDirs("laser_profile");


    // ── 5. Loop su combinazioni selezionate ──────────────────────────────
    // Strutture per sigma_vs_vbias
    std::vector<std::vector<MP>> allAcc;
    std::vector<int>    vbAcc;
    std::vector<double> frAcc;

    // Profilo radiale: accumula tutti i crossing
    TCanvas* cRad = new TCanvas("cRadial", "Profilo radiale P_{det}", 900, 650);
    cRad->SetGrid(); cRad->SetLeftMargin(0.12);
    cRad->SetBottomMargin(0.12); cRad->SetTopMargin(0.10);
    TLegend* legRad = new TLegend(0.55, 0.15, 0.93, 0.52);
    legRad->SetBorderSize(1); legRad->SetFillColor(0); legRad->SetTextSize(0.026);

    const int COLS[] = {kBlue+1, kRed+1, kGreen+2, kOrange+7,
                        kMagenta+1, kCyan+2, kBlack, kViolet+1};
    bool firstRad = true;
    TF1* fGrad = new TF1("fGrad", "[0]*exp(-x*x/(2.*[1]*[1]))", 0, 100);
    fGrad->SetLineWidth(2); fGrad->SetLineStyle(2);

    for (size_t ci2 = 0; ci2 < chosen.size(); ++ci2) {
        auto [vb, fr] = combos[chosen[ci2]];
        int col = COLS[ci2 % 8];

        // ── Trova file crossing e accepted per questo (vb, fr) ────────────
        std::string pathCross, pathAcc;
        for (const auto& f : filesCross)
            if (f.vbias == vb && std::abs(f.frac - fr) < 0.001) { pathCross = f.path; break; }
        for (const auto& f : filesAcc)
            if (f.vbias == vb && std::abs(f.frac - fr) < 0.001) { pathAcc = f.path; break; }

        std::vector<MP> ptsCross, ptsAcc;
        if (!pathCross.empty()) ptsCross = loadMap(pathCross, vb, fr, pdet_cut);
        if (!pathAcc.empty())   ptsAcc   = loadMap(pathAcc,   vb, fr, pdet_cut);

        // Se manca uno dei due, usa l'altro anche per il timing
        if (ptsAcc.empty() && !ptsCross.empty())   ptsAcc   = ptsCross;
        if (ptsCross.empty() && !ptsAcc.empty()) ptsCross = ptsAcc;
        if (ptsCross.empty()) continue;

        // ── Profilo radiale (crossing) ────────────────────────────────────
        TGraphErrors* grRad = radProfile(ptsCross, cx, cy_val,
                                          Form("grRad_%d", vb), pdet_cut);
        if (grRad) {
            grRad->SetMarkerStyle(20 + (int)(ci2 % 8));
            grRad->SetMarkerColor(col); grRad->SetLineColor(col); grRad->SetMarkerSize(1.2);
            cRad->cd();
            if (firstRad) {
                grRad->SetTitle(Form("Profilo radiale P_{det}   LET=%.2f p.e.;"
                                     "r dal centro (mm);P_{det} (%%)", fr));
                grRad->GetYaxis()->SetRangeUser(0, 110);
                grRad->GetXaxis()->SetRangeUser(0, 80);
                grRad->Draw("AP"); firstRad = false;
            } else grRad->Draw("P same");
            fGrad->SetLineColor(col);
            fGrad->SetParameters(95, 30);
            fGrad->SetParLimits(0,  0.0, 100.0);  // BUG FIX [3]
            fGrad->SetParLimits(1,  5.0,  80.0);
            grRad->Fit(fGrad, "RQ", "", 0, 75);
            TF1* fc = (TF1*)fGrad->Clone(Form("fcRad_%zu", ci2));
            fc->SetRange(0, 80); fc->Draw("same");
            legRad->AddEntry(grRad,
                Form("V_{bias}=%dV  #sigma=%.1fmm  P_{0}=%.1f%%",
                     vb, std::abs(fGrad->GetParameter(1)),
                     fGrad->GetParameter(0)), "pe");
        }

        // ── Mappa P_det 2D (crossing) ─────────────────────────────────────
        drawPdetMap(ptsCross, vb, fr, cx, cy_val, ctx);

        // ── Slice X/Y (crossing) ─────────────────────────────────────────
        drawSlices(ptsCross, vb, fr, cx, cy_val, ctx);

        // ── Mappa sigma_t 2D (accepted) ───────────────────────────────────
        if (!ptsAcc.empty()) {
            drawSigmaMap(ptsAcc, vb, fr, ctx);
            allAcc.push_back(ptsAcc);
            vbAcc.push_back(vb); frAcc.push_back(fr);
        }

        // ── Delta-mu map (accepted + maschera crossing) ───────────────────
        if (!ptsAcc.empty())
            drawDeltaMu(ptsAcc, ptsCross, vb, fr, cx, cy_val, pdet_cut, ctx);
    }

    // ── 6. Profilo radiale finale ────────────────────────────────────────
    if (!firstRad) {
        cRad->cd();
        legRad->Draw();
        TPaveText* inf = new TPaveText(0.13, 0.82, 0.52, 0.92, "NDC");
        inf->SetBorderSize(1); inf->SetFillColor(0); inf->SetTextSize(0.026);
        inf->AddText("Fit: P(r) = P_{0} exp(-r^{2}/2#sigma^{2})");
        inf->AddText("Diffusore 20#circ   P_{0}#leq100%   [da CROSSING]");
        inf->Draw();
        cRad->Update(); cRad->Modified();
        ctx.savePNG(cRad, "laser_radial_profile.png");
    }

    // ── 7. Sigma_t vs Vbias ──────────────────────────────────────────────
    drawSigmaVsVbias(allAcc, vbAcc, frAcc, ctx);

    // ── 8. Riepilogo ─────────────────────────────────────────────────────
    std::cout << "\n+============================================================+\n"
              << "|  DONE                                                      |\n"
              << "|  Output: " << ctx.pngDir << "\n"
              << "|                                                            |\n"
              << "|  Plot prodotti:                                            |\n"
              << "|   laser_radial_profile.png   — P_det(r) crossing          |\n"
              << "|   laser_map2D_v*.png         — mappa P_det crossing       |\n"
              << "|   laser_slice_v*.png         — slice X/Y crossing         |\n"
              << "|   sigma_map_v*.png           — sigma_t(x,y) accepted      |\n"
              << "|   sigma_vs_vbias.png         — sigma_t media vs Vbias     |\n"
              << "|   delta_mu_map2D_v*.png      — ritardo temporale 2D       |\n"
              << "|   delta_mu_3D_v*.png         — ritardo temporale 3D       |\n"
              << "+============================================================+\n";

    ctx.reopenAllCanvases(20);
}
