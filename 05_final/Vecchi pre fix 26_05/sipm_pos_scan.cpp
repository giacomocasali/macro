/**
 * sipm_pos_scan.cpp
 * =================
 * Scansione TOT su griglia di posizioni laser (x,y) per ogni Vbias × LET.
 *
 * SOSTITUISCE:
 *   - sipm_analisis_pro.cpp  (versione semplice, S+B only)
 *   - sipm_timing_map_v2.cpp (versione completa MA con logica duplicata)
 *
 * PRINCIPIO DRY:
 *   Riusa direttamente le funzioni gia' validate in VbiasAnalysis_v2.h:
 *     - collectTOTEvents_fileByFile   (build cache)
 *     - analyseOneLET_chunked         (S+B fit, sigma, MultiPE, time-walk)
 *   Per ogni posizione = un file raw "data_x_<X>_y_<Y>_vbias_<V>.root".
 *   Su quel file gira la STESSA pipeline che gira su tot_analysis ma
 *   con runMap a un solo elemento.
 *
 * INPUT FILE:
 *   <dataDir>/data_x_<X>_y_<Y>_vbias_<V>.root      (uno per posizione)
 *   <dataDir>/calib_vbias<V>_cut<C>mhz.root        (calibrazione per vbias)
 *
 * PARAMETRI INTERATTIVI:
 *   - cartella dati (subdir auto-listate)
 *   - lista Vbias
 *   - lista frac_pe (LET in unita' p.e.)
 *   - finestra fit [fit_lo, fit_hi] (ns)
 *   - time-walk on/off
 *   - filtro LP on/off
 *   - modalita' ORIGINAL / LOOSE
 *
 * OUTPUT:
 *   <dataDir>/map_results_vbias<V>_let<F>[_loose].root  con TTree "map":
 *     branches: x, y, mu, mu_err, mu_ok, sigma, sigma_err,
 *               n_sig, n_acc, n_crossing, n_laser, p_det, chi2ndf
 *   Canvas PNG: map_mu_*.png, map_sigma_*.png, map_pdet_*.png
 *   Compatibile con sipm_draw_timing3d() per il post-process 3D.
 *   mu_ok=1: S+B fit converged; mu_ok=0: fallback (mu non affidabile).
 *
 * Compile:  .L sipm_pos_scan.cpp+
 * Run:      sipm_pos_scan()
 */

// ── Standard library ────────────────────────────────────────────────────────
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ── ROOT ─────────────────────────────────────────────────────────────────────
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

// ── Project headers ──────────────────────────────────────────────────────────
// Include order: Config.h MUST be first (defines g_analysis_mode, g_data_dir_override
// as inline globals). VbiasAnalysis_v2.h pulls in everything else we need
// (collectTOTEvents_fileByFile, analyseOneLET_chunked, eventCachePath, etc.)
#include "../header/Config.h"
#include "../header/InputHelpers.h"
#include "../header/OutputManager.h"
#include "../header/CalibIO.h"
#include "../header/SignalProcessing.h"     // triggerWindowIndices
#include "../header/EventCache.h"           // eventCachePath
#include "../header/ChunkedHistoFill.h"     // countCacheEvents
#include "../header/SidebandAnalysis.h"     // SBResult, fitSplusB, cleanupSBResult
#include "../header/TimingCorrection.h"     // TWMethod, askTimeWalkMethod
#include "../header/VbiasAnalysis_v2.h"     // collectTOTEvents_fileByFile,
                                            // analyseOneLET_chunked

// ════════════════════════════════════════════════════════════════════════════
//  Struct: file della griglia
// ════════════════════════════════════════════════════════════════════════════
struct PositionFile {
    double      x = 0.0;
    double      y = 0.0;
    int         vbias = 0;
    std::string path;                  // path run-0 (per fs_MHz)
    std::map<int,std::string> runs;    // FIX multirun: tutti i run {N→path}
};

// ════════════════════════════════════════════════════════════════════════════
//  P_det computation mode
// ════════════════════════════════════════════════════════════════════════════
enum class PdetMode {
    ACCEPTED,   // N_sig / N_laser  (after all quality filters)
    CROSSING    // N_crossing / N_laser  (threshold crossed, no quality filters)
};

// ════════════════════════════════════════════════════════════════════════════
//  Struct: result for one position
// ════════════════════════════════════════════════════════════════════════════
struct MapPoint {
    double x = 0.0, y = 0.0;
    int    vbias = 0;
    double frac_pe = 0.0;

    double mu        = 0.0,  mu_err    = 0.0;
    double sigma     = 0.0,  sigma_err = 0.0;
    double n_sig     = 0.0;
    long   n_acc     = 0;
    long   n_crossing = 0;
    long   n_laser   = 0;
    double p_det     = -1.0;  // computed according to PdetMode
    double chi2ndf   = -1.0;

    bool   ok = false;
    bool   mu_reliable = false;  // true only when S+B fit_ok or sub_ok succeeded
    std::string sigma_source;   // "MultiPE", "SB", "fallback", "none"
};

// ════════════════════════════════════════════════════════════════════════════
//  parsePositionFile: parsa "data_x_<X>_y_<Y>_vbias_<V>.root"
//  X,Y interi (con segno), V intero positivo.
// ════════════════════════════════════════════════════════════════════════════
// FIX multirun: riconosce sia data_x_X_y_Y_vbias_V.root
//               che     data_x_X_y_Y_vbias_V_run_N.root
static bool parsePositionFile(const std::string& fname,
                               double& x, double& y, int& vbias, int& run)
{
    static const std::regex reRun(
        R"(^data_x_([\-\d]+)_y_([\-\d]+)_vbias_(\d+)_run_(\d+)\.root$)");
    static const std::regex reNoRun(
        R"(^data_x_([\-\d]+)_y_([\-\d]+)_vbias_(\d+)\.root$)");
    std::smatch m;
    if (std::regex_match(fname, m, reRun)) {
        try {
            x     = std::stod(m[1].str());
            y     = std::stod(m[2].str());
            vbias = std::stoi(m[3].str());
            run   = std::stoi(m[4].str());
        } catch (...) { return false; }
        return true;
    }
    if (std::regex_match(fname, m, reNoRun)) {
        try {
            x     = std::stod(m[1].str());
            y     = std::stod(m[2].str());
            vbias = std::stoi(m[3].str());
            run   = 0;
        } catch (...) { return false; }
        return true;
    }
    return false;
}

// ════════════════════════════════════════════════════════════════════════════
//  readLine / readDouble: input helpers (loop su EOF)
// ════════════════════════════════════════════════════════════════════════════
// ════════════════════════════════════════════════════════════════════════════
//  readLongParamLocal: legge TParameter<long> con fallback a <int>
// ════════════════════════════════════════════════════════════════════════════
static long readLongParamLocal(TFile* f, const char* name)
{
    if (!f || f->IsZombie()) return -1L;
    {
        auto* p = dynamic_cast<TParameter<long>*>(f->Get(name));
        if (p) return p->GetVal();
    }
    {
        auto* p = dynamic_cast<TParameter<int>*>(f->Get(name));
        if (p) return static_cast<long>(p->GetVal());
    }
    return -1L;
}

// ════════════════════════════════════════════════════════════════════════════
//  buildEdges: bordi bin di TH2D centrati sui valori unici
// ════════════════════════════════════════════════════════════════════════════
static std::vector<double> buildEdges(const std::vector<double>& vals)
{
    std::vector<double> edges;
    if (vals.empty()) return edges;
    edges.reserve(vals.size() + 1);
    // Bordo iniziale
    if (vals.size() == 1) {
        edges.push_back(vals[0] - 5.0);
        edges.push_back(vals[0] + 5.0);
        return edges;
    }
    edges.push_back(vals[0] - 0.5 * (vals[1] - vals[0]));
    for (size_t i = 1; i < vals.size(); ++i)
        edges.push_back(0.5 * (vals[i-1] + vals[i]));
    edges.push_back(vals.back() + 0.5 * (vals.back() - vals[vals.size()-2]));
    return edges;
}

// ════════════════════════════════════════════════════════════════════════════
//  drawMap2D: disegna TH2D mu/sigma/p_det e salva PNG via OutCtx
//  FIX: shift coordinate (cx,cy)->0, unita' mm, margini migliorati
//  FIX: outlier sigma > MAD-adattivo o |mu - mu_median| > MU_MAXDEV → bin vuoto
//  FIX: mediana mu calcolata solo su punti con mu_reliable=true
// ════════════════════════════════════════════════════════════════════════════
static void drawMap2D(const std::vector<MapPoint>& pts,
                      const std::string& quantity,
                      const std::string& zTitle,
                      int vbias, double frac_pe,
                      OutCtx& ctx,
                      const std::string& fname,
                      double cx = 0.0, double cy = 0.0,
                      bool drawCentre = true)
{
    // Soglie outlier
    constexpr double MU_MAXDEV   = 0.5;   // ns — |mu - mediana| > 0.5 ns = fit fallito
    // Soglia sigma: mediana + 10*MAD (robusto, si adatta ai dati)
    double sigma_threshold = 0.5; // fallback se non ci sono abbastanza punti
    {
        std::vector<double> sigVals;
        for (const auto& p : pts)
            if (p.ok && p.sigma > 0 && p.sigma < 0.5)
                sigVals.push_back(p.sigma);
        if (sigVals.size() >= 3) {
            std::sort(sigVals.begin(), sigVals.end());
            const double med = sigVals[sigVals.size() / 2];
            std::vector<double> absdev;
            for (double s : sigVals) absdev.push_back(std::abs(s - med));
            std::sort(absdev.begin(), absdev.end());
            const double mad = absdev[absdev.size() / 2];
            sigma_threshold = med + 10.0 * mad;
            std::cout << "  [drawMap2D:sigma] median=" << med*1000
                      << " ps  MAD=" << mad*1000
                      << " ps  threshold=" << sigma_threshold*1000 << " ps\n";
        }
    }

    // Calcola mediana mu SOLO dai punti con mu_reliable=true (S+B converged)
    // FIX: includere outlier mu_reliable=false (es. mu=30ns) inquina la mediana
    // e può spostare la soglia MU_MAXDEV in modo da mascherare punti buoni.
    std::vector<double> muVals;
    for (const auto& p : pts)
        if (p.ok && p.mu_reliable && p.sigma > 0 && p.sigma < 0.5)
            muVals.push_back(p.mu);
    double mu_median = 0.0;
    if (!muVals.empty()) {
        std::sort(muVals.begin(), muVals.end());
        mu_median = muVals[muVals.size() / 2];
    }

    // Raccoglie posizioni uniche shiftate: (x-cx, y-cy)
    // Usa TUTTE le posizioni ok per la griglia (anche outlier — bin vuoto)
    std::set<double> uxs, uys;
    for (const auto& p : pts) if (p.ok) {
        uxs.insert(p.x - cx);
        uys.insert(p.y - cy);
    }
    if (uxs.empty() || uys.empty()) return;

    std::vector<double> xv(uxs.begin(), uxs.end());
    std::vector<double> yv(uys.begin(), uys.end());
    std::vector<double> xEdges = buildEdges(xv);
    std::vector<double> yEdges = buildEdges(yv);
    const int nx = static_cast<int>(xv.size());
    const int ny = static_cast<int>(yv.size());

    TH2D* h = new TH2D(
        Form("h_%s_v%d_let%.2f", quantity.c_str(), vbias, frac_pe),
        Form("V_{bias}=%d V  LET=%.2f p.e.;#Deltax (mm);#Deltay (mm);%s",
             vbias, frac_pe, zTitle.c_str()),
        nx, xEdges.data(), ny, yEdges.data());
    h->SetDirectory(nullptr);

    int nFilled = 0, nMasked = 0;
    for (const auto& p : pts) {
        if (!p.ok) continue;

        // Filtro outlier per sigma e mu
        bool bad = false;
        if (quantity == "sigma" && p.sigma > sigma_threshold) bad = true;
        if (quantity == "mu"    && (std::abs(p.mu - mu_median) > MU_MAXDEV)) bad = true;
        if (bad) { ++nMasked; continue; }  // lascia bin vuoto = NaN in COLZ

        double val = 0.0;
        if      (quantity == "mu")         val = p.mu;
        else if (quantity == "sigma")      val = p.sigma;
        else if (quantity == "p_det")      val = (p.p_det >= 0 ? p.p_det * 100.0 : 0.0);
        else if (quantity == "n_acc")      val = static_cast<double>(p.n_acc);
        else if (quantity == "n_crossing") val = static_cast<double>(p.n_crossing);
        const int bx = h->GetXaxis()->FindBin(p.x - cx);
        const int by = h->GetYaxis()->FindBin(p.y - cy);
        h->SetBinContent(bx, by, val);
        h->SetBinError  (bx, by,
            (quantity == "mu")    ? p.mu_err :
            (quantity == "sigma") ? p.sigma_err : 0.0);
        ++nFilled;
    }
    if (nMasked > 0)
        std::cout << "  [drawMap2D:" << quantity << "] masked " << nMasked
                  << " outlier(s)  filled=" << nFilled << "\n";

    TCanvas* c = new TCanvas(
        Form("c_%s_v%d_let%.2f", quantity.c_str(), vbias, frac_pe),
        Form("%s map  Vbias=%d  LET=%.2f pe", zTitle.c_str(), vbias, frac_pe),
        960, 800);
    c->SetRightMargin(0.18f);
    c->SetLeftMargin(0.14f);
    c->SetBottomMargin(0.13f);
    c->SetTopMargin(0.10f);
    gStyle->SetPalette(kBird);
    h->SetContour(64);
    h->GetXaxis()->SetTitleSize(0.045f);
    h->GetYaxis()->SetTitleSize(0.045f);
    h->GetXaxis()->SetLabelSize(0.035f);
    h->GetYaxis()->SetLabelSize(0.035f);
    h->GetXaxis()->SetTitleOffset(1.1f);
    h->GetYaxis()->SetTitleOffset(1.3f);
    h->GetZaxis()->SetTitleSize(0.040f);
    h->GetZaxis()->SetLabelSize(0.032f);

    // FIX range Z automatico per tutte le quantità: evita scala da 0
    {
        double zMin = 1e9, zMax = -1e9;
        for (int bx = 1; bx <= h->GetNbinsX(); ++bx)
            for (int by = 1; by <= h->GetNbinsY(); ++by) {
                double v = h->GetBinContent(bx, by);
                if (v == 0.0) continue;
                if (v < zMin) zMin = v;
                if (v > zMax) zMax = v;
            }
        if (zMin < zMax) {
            double pad = std::max(0.01, (zMax - zMin) * 0.05);
            h->GetZaxis()->SetRangeUser(zMin - pad, zMax + pad);
        }
    }

    h->Draw("COLZ");

    // Centro marker ora è sempre a (0,0) dopo lo shift
    if (drawCentre) {
        TMarker* mc = new TMarker(0.0, 0.0, 5);
        mc->SetMarkerColor(kRed);
        mc->SetMarkerSize(2.0);
        mc->Draw();
    }

    c->Update(); c->Modified();
    ctx.savePNG(c, fname);   // savePNG fa delete del canvas
    delete h;
}

// ════════════════════════════════════════════════════════════════════════════
//  analyzePosition: per UNA posizione (x,y,vbias,frac_pe)
//    1. costruisce cache se manca (chiama collectTOTEvents_fileByFile)
//    2. esegue analisi (chiama analyseOneLET_chunked)
//    3. legge mu, N_sig dal fit S+B (analyseOneLET non li espone)
//    4. ritorna MapPoint
// ════════════════════════════════════════════════════════════════════════════
static MapPoint analyzePosition(const PositionFile&     pf,
                                 double                  frac_pe,
                                 double                  cutoff_MHz,
                                 const CalibResult&      cal,
                                 double                  fit_lo,
                                 double                  fit_hi,
                                 TWMethod                tw_method,
                                 bool                    use_filter,
                                 const std::string&      cacheSubdir,
                                 OutCtx&                 ctx,
                                 bool                    save_canvas_per_pos,
                                 PdetMode                pdet_mode)
{
    MapPoint mp;
    mp.x       = pf.x;
    mp.y       = pf.y;
    mp.vbias   = pf.vbias;
    mp.frac_pe = frac_pe;

    const double let_thr = cal.q + frac_pe * cal.m;

    // ── 1. fs_MHz e j_start/j_end dal file raw ───────────────────────────────
    double fs_MHz = 0.0;
    int    j_start = 0, j_end = 1023;
    {
        TFile* f0 = TFile::Open(pf.path.c_str(), "READ");
        if (!f0 || f0->IsZombie()) { delete f0; return mp; }
        TTree* tr = static_cast<TTree*>(f0->Get("ch1"));
        if (tr && tr->GetEntries() >= 2 && tr->GetBranch("time")) {
            const int N0 = 1024;
            Double_t tb[N0] = {};
            tr->SetBranchAddress("time", tb);
            tr->GetEntry(0);
            if (tb[1] - tb[0] > 0)
                fs_MHz = 1000.0 / (tb[1] - tb[0]);
            triggerWindowIndices(tb, N0,
                cal.t_trig_start, cal.t_trig_end,
                j_start, j_end,
                Form("vbias%d_x%.0f_y%.0f", pf.vbias, pf.x, pf.y));
        }
        f0->Close(); delete f0;
    }
    if (fs_MHz <= 0.0) {
        std::cout << "  [pos] invalid fs_MHz, skipping.\n";
        return mp;
    }

    // ── 2. cache path (per posizione) ────────────────────────────────────────
    // Suffix include x, y, calib hash e modalita' analisi.
    const std::string cacheSuffix = Form(
        "_pos_js%d_je%d_m%d_q%d_mode%d_x%d_y%d",
        j_start, j_end,
        (int)std::lround(cal.m * 10.0),
        (int)std::lround(cal.q * 10.0),
        g_analysis_mode,
        (int)std::round(pf.x), (int)std::round(pf.y));

    const std::string cachePath = eventCachePath(
        pf.vbias, frac_pe, cutoff_MHz, cal.laser_thr,
        cacheSubdir, use_filter, cacheSuffix);

    const std::string ltag = Form(
        "v%d_x%.0f_y%.0f_let%.2f",
        pf.vbias, pf.x, pf.y, frac_pe);

    // ── 3. Costruisce cache — FIX multirun: usa tutti i run disponibili ────
    const bool cached = (countCacheEvents(cachePath) > 0);
    if (!cached) {
        std::map<int, std::string> runMap = pf.runs;
        if (runMap.empty()) runMap[0] = pf.path;

        std::cout << "  [pos x=" << pf.x << " y=" << pf.y
                  << "] building cache from " << runMap.size() << " run(s)\n";

        collectTOTEvents_fileByFile(
            runMap, cutoff_MHz, fs_MHz,
            j_start, j_end, let_thr, cal,
            ltag, cacheSubdir,
            pf.vbias, frac_pe, use_filter,
            cacheSuffix,
            nullptr, nullptr);
    } else {
        std::cout << "  [pos x=" << pf.x << " y=" << pf.y
                  << "] cache hit (" << pf.runs.size() << " run(s))\n";
    }

    // ── 4. Analizza cache → sigma (MultiPE primario, SB fallback) ───────────
    //    analyseOneLET_chunked salva pure i canvas via ctx.savePNG.
    //    Se save_canvas_per_pos=false, viene ancora chiamato ma i canvas
    //    si accumulano in pngDir — per griglie grandi puo' essere lento.
    OneLETResult res = analyseOneLET_chunked(
        cachePath, frac_pe, let_thr,
        fit_lo, fit_hi, tw_method,
        false,         // do_pe_analysis: NO per griglia (troppi plot)
        ltag, pf.vbias, ctx);

    if (!res.valid) {
        std::cout << "  [pos x=" << pf.x << " y=" << pf.y
                  << "] FIT FAILED — skipping position.\n";
        if (!res.h2d_tmpfile.empty())
            gSystem->Unlink(res.h2d_tmpfile.c_str());
        return mp;
    }

    mp.sigma     = res.sigma;
    mp.sigma_err = res.sigmaErr;
    mp.ok        = true;
    mp.sigma_source = "MultiPE/SB (chunked)";

    // ── 5. mu, N_sig from S+B fit on cached events ──────────────────────────
    //    analyseOneLET_chunked does not expose mu and N_sig directly.
    //    Load events (up to 500k) and rerun fitSplusB with the same
    //    afterpulse cut (tot < tot_1pe_max).
    //    n_crossing is read from the cache TParameter (set by
    //    collectTOTEvents_fileByFile): it counts events where t_rise >= 0
    //    (threshold crossed) BEFORE any quality filter.
    {
        TFile* fc = TFile::Open(cachePath.c_str(), "READ");
        if (fc && !fc->IsZombie()) {
            mp.n_acc      = readLongParamLocal(fc, "n_accepted");
            mp.n_crossing = readLongParamLocal(fc, "n_crossing");
            // FIX bug #23 + FIX crossing-denominatore:
            // Gerarchia denominatore (dalla piu' corretta alla meno):
            //   1. n_laser_found  — laser trovato, NESSUN altro filtro
            //                       (nuovo, identico a n_wf di timing.cpp)
            //   2. n_valid_for_prob — laser trovato + baseline ok (vecchio fix #23)
            //   3. n_laser_tot    — totale grezzo (fallback per cache molto vecchie)
            mp.n_laser = readLongParamLocal(fc, "n_laser_found");
            if (mp.n_laser <= 0)
                mp.n_laser = readLongParamLocal(fc, "n_valid_for_prob");
            if (mp.n_laser <= 0)
                mp.n_laser = readLongParamLocal(fc, "n_laser_tot");
            fc->Close(); delete fc;
        }
    }

    // Stima tot_1pe_max con stessa logica di analyseOneLET (PASS 1)
    double tot_1pe_max = 0.0;
    {
        TH1D* hTp = new TH1D(Form("hTp_local_%s", ltag.c_str()),
            "", 1500, 0, 150);
        hTp->SetDirectory(nullptr);
        {
            ChunkedFiller p(cachePath);
            p.addTH1D(hTp, [](TH1D* h, const CacheEvent& e){
                h->Fill(e.tot); });
            p.run();
        }
        int nB = hTp->GetNbinsX();
        if (hTp->GetEntries() > 0) {
            int bPeak = 1;
            for (int b = 2; b <= nB; ++b)
                if (hTp->GetBinContent(b) > hTp->GetBinContent(bPeak)) bPeak = b;
            double peakVal = hTp->GetBinContent(bPeak);
            double peakTOT = hTp->GetBinCenter(bPeak);
            int bLo = std::max(bPeak + 1, hTp->FindBin(peakTOT + 2.0 + 1e-6));
            int bHi = std::min(nB - 1,    hTp->FindBin(35.0 - 1e-6));
            auto smooth3 = [&](int b) {
                double s = hTp->GetBinContent(b);
                if (b > 1)   s += hTp->GetBinContent(b-1);
                if (b < nB)  s += hTp->GetBinContent(b+1);
                return s / 3.0;
            };
            int bMin = -1; double minVal = 1e18;
            for (int b = bLo; b <= bHi; ++b) {
                double v = smooth3(b);
                if (v < minVal) { minVal = v; bMin = b; }
            }
            if (bMin > 0 && minVal < peakVal * 0.15) {
                tot_1pe_max = hTp->GetBinCenter(bMin);
            } else {
                double total = hTp->Integral(), cumul = 0;
                for (int b = 1; b <= nB; ++b) {
                    cumul += hTp->GetBinContent(b);
                    if (cumul >= total * 0.15) {
                        tot_1pe_max = hTp->GetBinCenter(b); break;
                    }
                }
            }
            tot_1pe_max = std::max(5.0, tot_1pe_max);
        }
        delete hTp;
    }

    // Carica eventi e fit S+B per mu, N_sig (con afterpulse cut)
    std::vector<TOTEvent> evts = loadEventsSmall(cachePath, 500000);
    evts.erase(std::remove_if(evts.begin(), evts.end(),
        [tot_1pe_max](const TOTEvent& e){ return e.tot >= tot_1pe_max; }),
        evts.end());

    // ── Peak finder adattivo per finestra S+B ──────────────────────────────
    // STRATEGIA: cerca il picco in una finestra stretta attorno al picco GIA'
    // trovato da analyseOneLET (res.mu se disponibile) oppure in [fit_lo,fit_hi].
    // Poi restringe la finestra di fit a [tpk-5, tpk+5] ns.
    //
    // FIX rispetto alla versione precedente: il picco veniva cercato su tutto
    // [fit_lo, fit_hi] = [0, 204.6] → per posizioni rumorose con pochi eventi
    // reali il massimo globale era rumore (es. x=130y=80: tpk=40ns invece di
    // 46.75ns). Ora si usa come ancora il mu di analyseOneLET (res.sigma e'
    // valido solo se res.valid=true, e in quel caso res contiene gia' il picco).
    // Come ancora usiamo il picco da analyseOneLET via la proiezione delta_t
    // dell'istogramma filtrato per tot < tot_1pe_max, gia' costruito sopra.
    constexpr double SB_FIT_HALF_WIN = 5.0;
    double sb_fit_lo = fit_lo;
    double sb_fit_hi = fit_hi;
    double tpk_found = -1.0;  // picco trovato, usato anche nel fallback sotto
    {
        TH1D* hPk = new TH1D(Form("hPkSB_%s", ltag.c_str()),
                              "", 1500, -50.0, 200.0);
        hPk->SetDirectory(nullptr);
        for (const auto& e : evts)
            hPk->Fill(e.delta_t);

        // FIX: cerca il massimo in una finestra STRETTA se abbiamo gia' un'ancora
        // dall'analisi precedente (la prima S+B su finestra larga di analyseOneLET
        // non e' disponibile qui, ma il mu di MultiPEFit si — pero' non viene
        // esposto da OneLETResult). Usiamo quindi la finestra [40, 55] ns come
        // ancora fisica: il segnale laser e' sempre in questo range per questo
        // setup (46.75 ns osservato). Se il picco in questa finestra ristretta
        // ha almeno il 30% del massimo globale, lo usiamo; altrimenti cerchiamo
        // sull'intera [fit_lo, fit_hi].
        const double ANCHOR_LO = 40.0;  // ns — limite inferiore finestra fisica segnale
        const double ANCHOR_HI = 55.0;  // ns — limite superiore finestra fisica segnale

        // Cerca max nella finestra fisica ristretta
        int bAncLo = std::max(1, hPk->FindBin(ANCHOR_LO + 1e-6));
        int bAncHi = std::min(hPk->GetNbinsX(), hPk->FindBin(ANCHOR_HI - 1e-6));
        int bMaxAnc = bAncLo;
        for (int b = bAncLo + 1; b <= bAncHi; ++b)
            if (hPk->GetBinContent(b) > hPk->GetBinContent(bMaxAnc)) bMaxAnc = b;
        double peakInAnchor = hPk->GetBinContent(bMaxAnc);
        double globalMax    = hPk->GetMaximum();

        if (peakInAnchor >= 0.30 * globalMax && peakInAnchor > 0) {
            // Segnale trovato nella finestra fisica: usa questo picco
            tpk_found = hPk->GetBinCenter(bMaxAnc);
            std::cout << "  [SBAdapt] peak in [" << ANCHOR_LO << "," << ANCHOR_HI
                      << "] @ " << std::fixed << std::setprecision(2) << tpk_found
                      << " ns  (frac=" << std::setprecision(2)
                      << peakInAnchor / globalMax << ")\n";
        } else {
            // Segnale debole o assente nella finestra fisica: cerca globalmente
            int bLo2 = std::max(1, hPk->FindBin(fit_lo + 1e-6));
            int bHi2 = std::min(hPk->GetNbinsX(), hPk->FindBin(fit_hi - 1e-6));
            int bMax2 = bLo2;
            for (int b = bLo2 + 1; b <= bHi2; ++b)
                if (hPk->GetBinContent(b) > hPk->GetBinContent(bMax2)) bMax2 = b;
            tpk_found = hPk->GetBinCenter(bMax2);
            std::cout << "  [SBAdapt] WARNING: weak signal in anchor window"
                      << " (frac=" << std::setprecision(2)
                      << peakInAnchor / std::max(1.0, globalMax)
                      << ") — using global peak @ "
                      << std::fixed << std::setprecision(2) << tpk_found << " ns\n";
        }
        delete hPk;

        if (evts.size() >= 5) {
            sb_fit_lo = std::max(fit_lo, tpk_found - SB_FIT_HALF_WIN);
            sb_fit_hi = std::min(fit_hi, tpk_found + SB_FIT_HALF_WIN);
            std::cout << "  [SBAdapt] fit window: ["
                      << sb_fit_lo << ", " << sb_fit_hi << "] ns\n";
        }
    }

    SBResult sb = fitSplusB(evts, sb_fit_lo, sb_fit_hi,
                             SB::DEFAULT_SB_WIDTH, SB::DEFAULT_EXT, ltag);

    if (sb.fit_ok) {
        mp.mu          = sb.mu;
        mp.mu_err      = sb.muErr;
        mp.n_sig       = sb.N_sig;
        mp.chi2ndf     = sb.chi2ndf;
        mp.mu_reliable = true;
    } else if (sb.sub_ok) {
        mp.mu          = sb.mu_sub;
        mp.mu_err      = sb.sigmaErr_sub;  // FIX: era 0.0 — usa errore del fit subtracted
        // FIX: n_sig non viene popolato nel path sub_ok di SidebandAnalysis.
        // Stima N_sig = N_tot_in_window - N_bkg_stimato.
        // N_bkg = bkg_per_bin * n_bins_nella_finestra_signal.
        {
            int nBinsSig = std::max(1, (int)std::round((sb_fit_hi - sb_fit_lo) / SB::BIN_WIDTH));
            double n_bkg_est = sb.bkg_per_bin * nBinsSig;
            long n_in_win = 0;
            for (const auto& e : evts)
                if (e.delta_t >= sb_fit_lo && e.delta_t < sb_fit_hi) ++n_in_win;
            mp.n_sig = std::max(0.0, static_cast<double>(n_in_win) - n_bkg_est);
        }
        mp.sigma_source = "MultiPE + SBsub mu";
        mp.mu_reliable  = true;
    } else {
        // Nessun fit converge: mu dalla media degli eventi nella finestra,
        // n_sig NON stimabile in modo affidabile → usa n_acc come upper bound
        // ma segnala che il punto ha mu NON affidabile (mu_reliable=false).
        double sw = 0, swx = 0;
        for (const auto& e : evts) {
            if (e.delta_t >= sb_fit_lo && e.delta_t <= sb_fit_hi) {
                sw += 1.0; swx += e.delta_t;
            }
        }
        double tpk_fallback = (sw > 0) ? swx / sw : tpk_found;
        mp.mu          = tpk_fallback;
        mp.mu_err      = SB_FIT_HALF_WIN;
        mp.n_sig       = static_cast<double>(mp.n_acc);  // upper bound, non affidabile
        mp.mu_reliable = false;
        mp.sigma_source = mp.sigma_source + " (mu=fallback,unreliable)";
        std::cout << "  [SBAdapt] WARNING: S+B fallback mu="
                  << std::fixed << std::setprecision(3) << mp.mu
                  << " ns  mu_reliable=false\n";
    }
    cleanupSBResult(sb);

    // ── 6. P_det computed according to PdetMode ─────────────────────────────
    //   ACCEPTED:  N_sig / N_laser   (events surviving all quality filters)
    //   CROSSING:  N_crossing / N_laser  (events where threshold was crossed,
    //              no afterpulse / time-window / quality cuts applied)
    //   N_crossing >= N_accepted  by construction.
    //
    //   FIX: se mu_reliable=false (S+B fallito del tutto), n_sig = n_acc che e'
    //   un upper bound non affidabile. In quel caso usa sempre CROSSING come
    //   stima piu' robusta, indipendentemente dalla modalita' scelta.
    if (pdet_mode == PdetMode::CROSSING || !mp.mu_reliable) {
        mp.p_det = (mp.n_laser > 0 && mp.n_crossing >= 0)
                   ? std::min(1.0, mp.n_crossing / static_cast<double>(mp.n_laser))
                   : -1.0;
        if (!mp.mu_reliable && pdet_mode == PdetMode::ACCEPTED)
            std::cout << "  [P_det] WARNING: S+B fallito, uso CROSSING come fallback P_det\n";
    } else {
        mp.p_det = (mp.n_laser > 0 && mp.n_sig >= 0)
                   ? std::min(1.0, mp.n_sig / static_cast<double>(mp.n_laser))
                   : -1.0;
    }

    // ── 7. Cleanup tmp file from analyseOneLET ──────────────────────────────
    if (!res.h2d_tmpfile.empty())
        gSystem->Unlink(res.h2d_tmpfile.c_str());

    std::cout << "  [pos x=" << std::fixed << std::setprecision(0) << pf.x
              << " y=" << pf.y << "]"
              << std::setprecision(3)
              << "  mu=" << mp.mu << " ns"
              << (mp.mu_reliable ? "" : " [UNRELIABLE]")
              << "  sigma=" << mp.sigma << " ns"
              << "  N_acc=" << mp.n_acc
              << "  N_cross=" << mp.n_crossing
              << "  P_det=" << std::setprecision(1)
              << (mp.p_det >= 0 ? mp.p_det * 100.0 : -1.0) << "%"
              << "  [" << (pdet_mode == PdetMode::CROSSING ? "CROSSING" : "ACCEPTED")
              << (!mp.mu_reliable && pdet_mode == PdetMode::ACCEPTED ? "->CROSSING_FALLBACK" : "")
              << "]\n";

    return mp;
}

// ════════════════════════════════════════════════════════════════════════════
//  saveMapResults: salva TTree "map" su file ROOT
//  Schema compatibile con sipm_draw_timing3d (legge questo file).
// ════════════════════════════════════════════════════════════════════════════
static void saveMapResults(const std::vector<MapPoint>& pts,
                            int vbias, double frac_pe,
                            const std::string& dataDir,
                            PdetMode pdet_mode)
{
    const std::string modeTag =
        (pdet_mode == PdetMode::CROSSING) ? "_pdet_crossing" : "_pdet_accepted";
    const std::string resPath = dataDir +
        Form("/map_results_vbias%d_let%.2f%s%s.root",
             vbias, frac_pe,
             (g_analysis_mode == 1 ? "_loose" : ""),
             modeTag.c_str());

    TFile* f = TFile::Open(resPath.c_str(), "RECREATE");
    if (!f || f->IsZombie()) {
        std::cerr << "  [MAP] cannot write " << resPath << "\n";
        delete f; return;
    }

    Double_t r_x, r_y, r_mu, r_mue, r_sig, r_sige, r_nsig, r_pdet, r_chi2;
    Long64_t r_nacc, r_ncross, r_nlas;
    Int_t    r_vbias, r_mu_ok;

    TTree* t = new TTree("map",
        Form("Position scan vbias%d let%.2fpe%s P_det=%s",
             vbias, frac_pe,
             (g_analysis_mode == 1 ? " LOOSE" : ""),
             (pdet_mode == PdetMode::CROSSING ? "CROSSING" : "ACCEPTED")));

    t->Branch("x",          &r_x,     "x/D");
    t->Branch("y",          &r_y,     "y/D");
    t->Branch("vbias",      &r_vbias, "vbias/I");
    t->Branch("mu",         &r_mu,    "mu/D");
    t->Branch("mu_err",     &r_mue,   "mu_err/D");
    t->Branch("mu_ok",      &r_mu_ok, "mu_ok/I");   // 1=S+B converged, 0=fallback
    t->Branch("sigma",      &r_sig,   "sigma/D");
    t->Branch("sigma_err",  &r_sige,  "sigma_err/D");
    t->Branch("n_sig",      &r_nsig,  "n_sig/D");
    t->Branch("n_acc",      &r_nacc,  "n_acc/L");
    t->Branch("n_crossing", &r_ncross,"n_crossing/L");
    t->Branch("n_laser",    &r_nlas,  "n_laser/L");
    t->Branch("p_det",      &r_pdet,  "p_det/D");
    t->Branch("chi2ndf",    &r_chi2,  "chi2ndf/D");

    int nOk = 0;
    for (const auto& p : pts) {
        if (!p.ok) continue;
        r_x      = p.x;        r_y       = p.y;
        r_vbias  = p.vbias;
        r_mu     = p.mu;       r_mue     = p.mu_err;
        r_mu_ok  = p.mu_reliable ? 1 : 0;
        r_sig    = p.sigma;    r_sige    = p.sigma_err;
        r_nsig   = p.n_sig;
        r_nacc   = static_cast<Long64_t>(p.n_acc);
        r_ncross = static_cast<Long64_t>(p.n_crossing);
        r_nlas   = static_cast<Long64_t>(p.n_laser);
        r_pdet   = p.p_det;
        r_chi2   = p.chi2ndf;
        t->Fill();
        ++nOk;
    }
    t->Write();
    f->Close(); delete f;
    std::cout << "  [MAP] " << nOk << " points saved -> " << resPath << "\n";
}

// ════════════════════════════════════════════════════════════════════════════
//  saveMapTable: salva tabella TXT + canvas grafica con tutti i dati per posizione
//  Coordinate shiftate: (x-cx, y-cy), unita' mm
//  TXT contiene SEMPRE sia crossing che accepted
// ════════════════════════════════════════════════════════════════════════════
static void saveMapTable(const std::vector<MapPoint>& pts,
                          int vbias, double frac_pe,
                          const std::string& dataDir,
                          OutCtx& ctx,
                          PdetMode /*pdet_mode*/,   // non usato — sempre entrambi
                          double cx, double cy)
{
    const std::string tag = Form("vbias%d_let%.2f%s",
        vbias, frac_pe,
        (g_analysis_mode == 1 ? "_loose" : ""));

    // ── 1. TXT — contiene tutti i valori inclusi crossing e accepted ─────────
    const std::string txtPath = dataDir + "/map_table_" + tag + ".txt";
    {
        std::ofstream ofs(txtPath);
        if (!ofs.is_open()) {
            std::cerr << "  [TABLE] Cannot write TXT: " << txtPath << "\n";
        } else {
            ofs << "# SiPM Position Scan  Vbias=" << vbias
                << " V  LET=" << std::fixed << std::setprecision(2) << frac_pe
                << " p.e.\n"
                << "# Centre: (" << cx << ", " << cy << ")  ->  shift to (0,0)\n"
                << "# Units: mm\n"
                << "# mu_ok: 1=S+B fit converged, 0=fallback (mu unreliable)\n"
                << "# Pdet_acc uses N_sig from S+B fit; if mu_ok=0 uses N_crossing instead\n#\n"
                << std::setw(10) << "dx(mm)"
                << std::setw(10) << "dy(mm)"
                << std::setw(12) << "mu(ns)"
                << std::setw(12) << "mu_err(ns)"
                << std::setw(6)  << "mu_ok"
                << std::setw(12) << "sigma(ns)"
                << std::setw(12) << "sig_err(ns)"
                << std::setw(10) << "N_acc"
                << std::setw(12) << "N_cross"
                << std::setw(10) << "N_laser"
                << std::setw(12) << "Pdet_cross%"
                << std::setw(12) << "Pdet_acc%"
                << "\n" << std::string(130, '-') << "\n";
            for (const auto& p : pts) {
                if (!p.ok) continue;
                const double pdet_cross = (p.n_laser > 0)
                    ? 100.0 * p.n_crossing / p.n_laser : -1.0;
                // Pdet_acc: usa n_sig se mu affidabile, altrimenti crossing
                const double pdet_acc = (p.n_laser > 0)
                    ? (p.mu_reliable
                        ? 100.0 * p.n_sig / p.n_laser
                        : pdet_cross)  // fallback a crossing se S+B fallito
                    : -1.0;
                ofs << std::fixed
                    << std::setw(10) << std::setprecision(1) << (p.x - cx)
                    << std::setw(10) << std::setprecision(1) << (p.y - cy)
                    << std::setw(12) << std::setprecision(4) << p.mu
                    << std::setw(12) << std::setprecision(4) << p.mu_err
                    << std::setw(6)  << (p.mu_reliable ? 1 : 0)
                    << std::setw(12) << std::setprecision(4) << p.sigma
                    << std::setw(12) << std::setprecision(4) << p.sigma_err
                    << std::setw(10) << p.n_acc
                    << std::setw(12) << p.n_crossing
                    << std::setw(10) << p.n_laser
                    << std::setw(12) << std::setprecision(1) << pdet_cross
                    << std::setw(12) << std::setprecision(1) << pdet_acc
                    << "\n";
            }
            ofs.close();
            std::cout << "  [TABLE] TXT -> " << txtPath << "\n";
        }
    }

    // ── 2. Canvas grafica ─────────────────────────────────────────────────────
    std::vector<const MapPoint*> valid;
    for (const auto& p : pts) if (p.ok) valid.push_back(&p);
    if (valid.empty()) return;

    const int nRows = static_cast<int>(valid.size());
    // Altezza pixel: header ~60px + titolo ~40px + ogni riga 22px + margine 40px
    const int canH  = 140 + nRows * 24;
    const int canW  = 1300;

    TCanvas* ct = new TCanvas(
        Form("ct_%s", tag.c_str()),
        Form("Summary table  Vbias=%d  LET=%.2f", vbias, frac_pe),
        canW, canH);
    ct->SetFillColor(0);
    ct->cd();

    // NDC altezza riga proporzionale alla canvas
    const double pH = 1.0 / (nRows + 5.0);

    // Titolo
    TLatex tit;
    tit.SetNDC(); tit.SetTextSize(std::min(0.040, pH * 0.8));
    tit.SetTextFont(62);
    tit.DrawLatex(0.03, 0.97,
        Form("Vbias=%d V   LET=%.2f p.e.   Centre=(%.0f,%.0f)   CROSSING + ACCEPTED",
             vbias, frac_pe, cx, cy));

    // Header colonne
    const char* headers[] = {
        "#Deltax(mm)", "#Deltay(mm)",
        "#mu(ns)", "#mu_{err}",
        "#sigma(ns)", "#sigma_{err}",
        "N_{acc}", "N_{cross}", "N_{laser}",
        "P_{cross}%", "P_{acc}%"
    };
    const int nCols = 11;
    const double colX[] = {
        0.02, 0.10, 0.18, 0.27, 0.36, 0.45,
        0.54, 0.63, 0.72, 0.82, 0.91
    };

    double yH = 0.93 - pH * 0.3;
    TLatex hdr;
    hdr.SetNDC(); hdr.SetTextSize(std::min(0.028, pH * 0.7)); hdr.SetTextFont(62);
    for (int c = 0; c < nCols; ++c)
        hdr.DrawLatex(colX[c], yH, headers[c]);

    TLine ln; ln.SetNDC(); ln.SetLineWidth(1);
    ln.DrawLineNDC(0.01, yH - pH * 0.15, 0.99, yH - pH * 0.15);

    TLatex val;
    val.SetNDC(); val.SetTextSize(std::min(0.025, pH * 0.65)); val.SetTextFont(42);

    double yRow = yH - pH * 0.15 - pH * 0.85;
    for (int r = 0; r < nRows; ++r) {
        const MapPoint& p = *valid[r];
        if (r % 2 == 0) {
            TPave* bg = new TPave(0.01, yRow - pH * 0.1,
                                  0.99, yRow + pH * 0.85, 0, "NDC");
            bg->SetFillColor(18); bg->SetLineColor(0); bg->Draw();
        }
        // Riga in rosso se mu non affidabile (S+B fallito)
        val.SetTextColor(p.mu_reliable ? kBlack : kRed+1);

        const double pdc = (p.n_laser > 0) ? 100.0 * p.n_crossing / p.n_laser : -1.0;
        // pda: usa n_sig se mu affidabile, altrimenti crossing (stesso criterio del TXT)
        const double pda = (p.n_laser > 0)
            ? (p.mu_reliable ? 100.0 * p.n_sig / p.n_laser : pdc)
            : -1.0;

        // Mu con asterisco se non affidabile
        const std::string muStr = p.mu_reliable
            ? Form("%.4f", p.mu)
            : Form("%.4f*", p.mu);

        val.DrawLatex(colX[0],  yRow, Form("%.1f",  p.x - cx));
        val.DrawLatex(colX[1],  yRow, Form("%.1f",  p.y - cy));
        val.DrawLatex(colX[2],  yRow, muStr.c_str());
        val.DrawLatex(colX[3],  yRow, Form("%.4f",  p.mu_err));
        val.DrawLatex(colX[4],  yRow, Form("%.4f",  p.sigma));
        val.DrawLatex(colX[5],  yRow, Form("%.4f",  p.sigma_err));
        val.DrawLatex(colX[6],  yRow, Form("%ld",   p.n_acc));
        val.DrawLatex(colX[7],  yRow, Form("%ld",   p.n_crossing));
        val.DrawLatex(colX[8],  yRow, Form("%ld",   p.n_laser));
        val.DrawLatex(colX[9],  yRow, Form("%.1f",  pdc));
        val.DrawLatex(colX[10], yRow, Form("%.1f",  pda));
        yRow -= pH;
    }
    val.SetTextColor(kBlack);  // reset colore

    ct->Update(); ct->Modified();
    ctx.savePNG(ct, "map_table_" + tag + ".png");
}

// ════════════════════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════════════════════
void sipm_pos_scan()
{
    // Reset any override left by a previous macro in the same ROOT session.
    g_data_dir_override = "";

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::cout << "\n+==========================================================+\n"
              << "|  SiPM POSITION SCAN — grid (x,y) x Vbias x LET          |\n"
              << "|  reuses collectTOTEvents_fileByFile + analyseOneLET      |\n"
              << "+==========================================================+\n\n";

    try {

    // ── 0. DATA DIRECTORY ───────────────────────────────────────────────────
    //   Derives the project data/ root from DATA_DIR (strips the last path
    //   segment, e.g. .../data/data_filter_3 -> .../data), lists all
    //   subdirectories, and proposes data_filter_3 as default.
    //   Env override: SIPM_DATA_DIR skips the selection entirely.
    std::string dataDir;
    {
        // 1) explicit env override: skip selection entirely
        if (const char* envFull = std::getenv("SIPM_DATA_DIR");
            envFull && envFull[0] != '\0') {
            dataDir = envFull;
        } else {
            // 2) derive root = parent of DATA_DIR
            std::string root = DATA_DIR;
            while (!root.empty() && root.back() == '/') root.pop_back();
            const size_t sl = root.find_last_of("/\\");
            if (sl != std::string::npos) root = root.substr(0, sl);

            std::vector<std::string> subs;
            void* dp = gSystem->OpenDirectory(root.c_str());
            if (dp) {
                const char* ent = nullptr;
                while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
                    std::string s(ent);
                    if (s == "." || s == "..") continue;
                    FileStat_t st;
                    if (gSystem->GetPathInfo((root + "/" + s).c_str(), st) == 0 &&
                        R_ISDIR(st.fMode)) subs.push_back(s);
                }
                gSystem->FreeDirectory(dp);
            }
            std::sort(subs.begin(), subs.end());

            // 3) preselect giacomo_serpentina as default, otherwise first entry
            int defIdx = -1;
            const std::string preferred = "giacomo_serpentina";
            for (size_t i = 0; i < subs.size(); ++i)
                if (subs[i] == preferred) { defIdx = static_cast<int>(i) + 1; break; }
            if (defIdx < 0 && !subs.empty()) defIdx = 1;

            if (!subs.empty()) {
                std::cout << "  Folders in " << root << ":\n";
                for (size_t i = 0; i < subs.size(); ++i) {
                    const bool isDef = (static_cast<int>(i) + 1 == defIdx);
                    std::cout << "    [" << (i+1) << "] " << subs[i]
                              << (isDef ? "   <-- default" : "") << "\n";
                }
                const std::string l = readLineOrEmpty(
                    Form("\n  Choose [n] or path  [ENTER = %d]: ", defIdx));
                if (l.empty()) {
                    dataDir = root + "/" + subs[defIdx - 1];
                } else {
                    try {
                        size_t pos = 0;
                        const int idx = std::stoi(l, &pos);
                        if (pos == l.size() && idx >= 1 &&
                            idx <= static_cast<int>(subs.size()))
                            dataDir = root + "/" + subs[idx-1];
                        else
                            dataDir = l;
                    } catch (...) { dataDir = l; }
                }
            } else {
                std::cerr << "  [WARN] Could not list parent of DATA_DIR: "
                          << root << "\n";
                dataDir = readLine("  Data folder path: ");
            }
        }

        while (!dataDir.empty() && dataDir.back() == '/') dataDir.pop_back();
        if (gSystem->AccessPathName(dataDir.c_str())) {
            std::cerr << "[ERR] Not accessible: " << dataDir << "\n";
            return;
        }
        g_data_dir_override = dataDir;
        std::cout << "  --> " << dataDir << "\n\n";
    }

    // ── 1. SCAN flat files — FIX multirun: raggruppa run per posizione ──────
    std::vector<PositionFile> allFiles;
    {
        std::map<std::tuple<double,double,int>, PositionFile> posMap;

        void* dp = gSystem->OpenDirectory(dataDir.c_str());
        if (!dp) { std::cerr << "[ERR] Cannot open directory.\n"; return; }
        const char* ent = nullptr;
        while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
            std::string fn(ent);
            double x = 0, y = 0; int vb = 0, run = 0;
            if (!parsePositionFile(fn, x, y, vb, run)) continue;
            auto key = std::make_tuple(x, y, vb);
            auto& pf = posMap[key];
            pf.x = x; pf.y = y; pf.vbias = vb;
            pf.runs[run] = dataDir + "/" + fn;
            if (run == 0 || pf.path.empty())
                pf.path = dataDir + "/" + fn;
        }
        gSystem->FreeDirectory(dp);

        for (auto& [key, pf] : posMap) {
            std::cout << "  pos x=" << pf.x << " y=" << pf.y
                      << " vbias=" << pf.vbias
                      << " runs=" << pf.runs.size() << "\n";
            allFiles.push_back(pf);
        }
    }
    if (allFiles.empty()) {
        std::cerr << "[ERR] No data_x_*_y_*_vbias_*.root files found in:\n"
                  << "      " << dataDir << "\n";
        return;
    }

    std::set<int> availVbias;
    for (const auto& f : allFiles) availVbias.insert(f.vbias);
    std::cout << "  Files found: " << allFiles.size()
              << "  Available Vbias:";
    for (int v : availVbias) std::cout << " " << v;
    std::cout << "\n\n";

    // ── 2. ANALYSIS MODE ────────────────────────────────────────────────────
    {
        std::cout << "  Analysis mode: 0=ORIGINAL (strict), 1=LOOSE (permissive)\n";
        while (true) {
            const std::string l = readLine("  Choose [0/1]: ");
            if (l == "0") { g_analysis_mode = 0; break; }
            if (l == "1") { g_analysis_mode = 1; break; }
        }
        std::cout << "  --> "
                  << (g_analysis_mode == 0 ? "ORIGINAL" : "LOOSE") << "\n\n";
    }

    // ── 3. VBIAS LIST ───────────────────────────────────────────────────────
    std::vector<int> vbiasList;
    {
        const std::string l = readLine("  Vbias list (e.g. 55 56): ");
        std::stringstream ss(l);
        int v = 0;
        while (ss >> v) {
            if (availVbias.count(v)) vbiasList.push_back(v);
            else std::cerr << "  [WARN] Vbias=" << v << " not in folder, skipping.\n";
        }
    }
    if (vbiasList.empty()) {
        std::cerr << "[ERR] No valid Vbias selected.\n"; return;
    }

    // ── 4. FRAC_PE LIST (range or list) ─────────────────────────────────────
    std::vector<double> fracList;
    {
        std::cout << "\n  frac_pe — range (lo hi step) or space-separated list:\n";
        const std::string l = readLine("  > ");
        std::vector<double> nums;
        { std::stringstream ss(l); double v = 0; while (ss >> v) nums.push_back(v); }
        const bool isRange =
            nums.size() == 3 && nums[0] > 0 && nums[1] > nums[0] &&
            nums[2] > 0 && nums[2] <= (nums[1] - nums[0]) + 1e-9;
        if (isRange) {
            const int n = static_cast<int>(std::round((nums[1] - nums[0]) / nums[2]));
            for (int i = 0; i <= n; ++i) {
                const double v = nums[0] + i * nums[2];
                if (v <= nums[1] + 1e-9) fracList.push_back(v);
            }
        } else {
            for (double v : nums) if (v > 0) fracList.push_back(v);
        }
        if (fracList.empty()) fracList.push_back(1.0);
        std::cout << "  --> " << fracList.size() << " threshold(s)\n";
    }

    // ── 5. CUTOFF / FIT WINDOW / OPTIONS ────────────────────────────────────
    double cutoff_global = 0.0;
    {
        const std::string l = readLine("\n  Cutoff MHz [0 = auto-detect / no filter]: ");
        try { cutoff_global = std::stod(l); } catch (...) {}
    }
    const bool use_filter = (cutoff_global > 0.0);

    const double fit_lo = readDouble("\n  Fit window start [ns]: ");
    const double fit_hi = readDouble("  Fit window end   [ns]: ");
    if (fit_lo >= fit_hi) {
        std::cerr << "[ERR] fit_lo >= fit_hi.\n"; return;
    }

    const TWMethod tw_method = askTimeWalkMethod();

    const bool save_canvas_per_pos = readYN(
        "\n  Save canvas for each position? [y/n] (slow, many files): ");

    // ── 5b. P_DET MODE ──────────────────────────────────────────────────────
    //   ACCEPTED  = N_sig / N_laser
    //               Events that pass all quality filters (afterpulse cut,
    //               time-walk correction, fit time-window). Conservative.
    //
    //   CROSSING  = N_crossing / N_laser
    //               Events where the waveform crossed the LET threshold at
    //               least once (t_rise >= 0), regardless of quality filters.
    //               Always >= ACCEPTED by construction. Less conservative,
    //               closer to the physical detection probability.
    PdetMode pdet_mode;
    {
        std::cout << "\n  P_det mode:\n"
                  << "    [0] ACCEPTED   — N_sig / N_laser  (after all quality filters)\n"
                  << "    [1] CROSSING   — N_crossing / N_laser  (threshold crossed, no filters)\n";
        while (true) {
            const std::string l = readLine("  Choose [0/1]: ");
            if (l == "0") { pdet_mode = PdetMode::ACCEPTED;  break; }
            if (l == "1") { pdet_mode = PdetMode::CROSSING;  break; }
            std::cerr << "  [!] Enter 0 or 1.\n";
        }
        std::cout << "  --> P_det mode: "
                  << (pdet_mode == PdetMode::CROSSING ? "CROSSING" : "ACCEPTED")
                  << "\n\n";
    }

    // ── 6. OUTPUT CONTEXT ───────────────────────────────────────────────────
    std::string runTag = "pos_scan";
    for (int v : vbiasList) runTag += "_v" + std::to_string(v);
    runTag += (g_analysis_mode == 1) ? "_loose" : "_orig";
    runTag += (pdet_mode == PdetMode::CROSSING) ? "_crossing" : "_accepted";
    OutCtx ctx = createOutputDirs(runTag);

    // ── 6b. CENTRO NOMINALE (per shift coordinate e marker) ─────────────────
    double map_cx = 0.0, map_cy = 0.0;
    {
        // Suggerisci il centro come mediana delle coordinate disponibili
        std::set<double> allX, allY;
        for (const auto& f : allFiles) { allX.insert(f.x); allY.insert(f.y); }
        if (!allX.empty()) {
            auto it = allX.begin();
            std::advance(it, allX.size() / 2);
            map_cx = *it;
        }
        if (!allY.empty()) {
            auto it = allY.begin();
            std::advance(it, allY.size() / 2);
            map_cy = *it;
        }
        std::cout << "\n  Centro nominale (x y) per shift coordinate [default "
                  << map_cx << " " << map_cy << "]: ";
        const std::string l = readLineOrEmpty("");
        if (!l.empty()) {
            std::istringstream ss(l);
            double a = 0, b = 0;
            if (ss >> a >> b) { map_cx = a; map_cy = b; }
        }
        std::cout << "  --> Centro = (" << map_cx << ", " << map_cy << ")\n\n";
    }

    // ── 7. LOOP PRINCIPALE (Vbias × frac_pe × posizioni) ────────────────────
    const auto t0 = std::chrono::steady_clock::now();

    for (int vbias : vbiasList) {

        // ── 7a. Carica calibrazione ──────────────────────────────────────────
        double cutoff_cal = cutoff_global;
        if (cutoff_cal == 0.0) {
            // auto-detect: prende il cutoff piu' grande disponibile
            const std::string pre = "calib_vbias" + std::to_string(vbias) + "_cut";
            const std::string suf = "mhz.root";
            std::vector<double> fc;
            void* dp = gSystem->OpenDirectory(dataDir.c_str());
            if (dp) {
                const char* ent = nullptr;
                while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
                    std::string fn(ent);
                    if (fn.size() < pre.size() + suf.size()) continue;
                    if (fn.substr(0, pre.size()) != pre) continue;
                    if (fn.substr(fn.size() - suf.size()) != suf) continue;
                    const std::string mid = fn.substr(
                        pre.size(), fn.size() - pre.size() - suf.size());
                    try { fc.push_back(std::stod(mid)); } catch (...) {}
                }
                gSystem->FreeDirectory(dp);
            }
            std::sort(fc.begin(), fc.end());
            if (!fc.empty()) cutoff_cal = fc.back();
        }

        CalibResult cal;
        if (!loadCalibration(cal, vbias, cutoff_cal, dataDir)) {
            std::cerr << "[WARN] Vbias=" << vbias
                      << ": calibration not found (cutoff="
                      << cutoff_cal << " MHz), skipping.\n";
            continue;
        }
        if (cal.m < GAIN_MIN || cal.m > GAIN_MAX) {
            std::cerr << "[WARN] Vbias=" << vbias << ": gain " << cal.m
                      << " out of range [" << GAIN_MIN << ", "
                      << GAIN_MAX << "], skipping.\n";
            continue;
        }
        // FIX CRITICO: sovrascrive sempre t_trig_start/end con i valori
        // passati dall'utente da terminale. I file .root di calibrazione
        // vecchi contengono trig_start=95, trig_end=125 (hardcoded bug),
        // che porterebbero triggerWindowIndices a j_start~475, j_end~640
        // escludendo completamente il segnale laser a ~46.77 ns.
        // Il parse ha gia' chiesto fit_lo/fit_hi all'utente — quelli
        // hanno precedenza assoluta su qualsiasi valore nel file .root.
        if (cal.t_trig_start != fit_lo || cal.t_trig_end != fit_hi) {
            std::cout << "  [CalibFix] Overriding trig window from file ["
                      << cal.t_trig_start << ", " << cal.t_trig_end
                      << "] ns with user input ["
                      << fit_lo << ", " << fit_hi << "] ns\n";
            cal.t_trig_start = fit_lo;
            cal.t_trig_end   = fit_hi;
        }

        std::cout << "\n+==========================================================+\n"
                  << "|  Vbias=" << vbias
                  << " V  gain=" << cal.m << " mV/pe"
                  << "  laser_thr=" << cal.laser_thr << " mV"
                  << "  cutoff=" << cutoff_cal << " MHz\n"
                  << "+==========================================================+\n";

        // Files for this Vbias, sorted by x then y
        std::vector<PositionFile> posFiles;
        for (const auto& f : allFiles)
            if (f.vbias == vbias) posFiles.push_back(f);
        std::sort(posFiles.begin(), posFiles.end(),
            [](const PositionFile& a, const PositionFile& b){
                if (a.x != b.x) return a.x < b.x;
                return a.y < b.y;
            });

        for (double frac_pe : fracList) {

            const double let_thr = cal.q + frac_pe * cal.m;
            std::cout << "\n+--- LET=" << frac_pe << " pe"
                      << "  thr=" << std::fixed << std::setprecision(2)
                      << let_thr << " mV"
                      << "  positions=" << posFiles.size()
                      << "  P_det="
                      << (pdet_mode == PdetMode::CROSSING ? "CROSSING" : "ACCEPTED")
                      << " ---\n";

            // Cache subdir for this (vbias, frac_pe) combination
            const std::string cacheSubdir = dataDir
                + "/cache_vbias" + std::to_string(vbias)
                + Form("_let%.2f", frac_pe)
                + (use_filter
                    ? Form("_filt%d", (int)std::round(cutoff_cal))
                    : "_nofilt")
                + (g_analysis_mode == 1 ? "_loose" : "");
            gSystem->mkdir(cacheSubdir.c_str(), true);

            // Process each position
            std::vector<MapPoint> mapPts;
            mapPts.reserve(posFiles.size());
            int nOk = 0, nFail = 0;
            const auto tStart = std::chrono::steady_clock::now();

            for (const auto& pf : posFiles) {
                if (gROOT->IsInterrupted()) goto end_outer;

                std::cout << "\n  [pos x=" << std::setw(4) << (int)pf.x
                          << " y=" << std::setw(4) << (int)pf.y << "]\n";

                MapPoint mp = analyzePosition(
                    pf, frac_pe, cutoff_cal, cal,
                    fit_lo, fit_hi, tw_method, use_filter,
                    cacheSubdir, ctx, save_canvas_per_pos,
                    pdet_mode);

                mapPts.push_back(mp);
                if (mp.ok) ++nOk; else ++nFail;
            }

            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - tStart).count();
            std::cout << "\n  Positions OK=" << nOk << "  FAIL=" << nFail
                      << "  elapsed=" << std::fixed << std::setprecision(1)
                      << elapsed << " s\n";

            // Save results
            saveMapResults(mapPts, vbias, frac_pe, dataDir, pdet_mode);

            // Draw 2D maps
            const std::string tag = Form("vbias%d_let%.2f%s_%s",
                vbias, frac_pe,
                (g_analysis_mode == 1 ? "_loose" : ""),
                (pdet_mode == PdetMode::CROSSING ? "crossing" : "accepted"));

            // ── 7 canvas fisse — indipendenti da pdet_mode ──────────────────
            const std::string tagBase = Form("vbias%d_let%.2f%s",
                vbias, frac_pe,
                (g_analysis_mode == 1 ? "_loose" : ""));

            // 1. mu
            drawMap2D(mapPts, "mu", "#mu(#Delta t) (ns)",
                      vbias, frac_pe, ctx,
                      Form("map_mu_%s.png", tagBase.c_str()),
                      map_cx, map_cy);
            // 2. sigma
            drawMap2D(mapPts, "sigma", "#sigma(#Delta t) (ns)",
                      vbias, frac_pe, ctx,
                      Form("map_sigma_%s.png", tagBase.c_str()),
                      map_cx, map_cy);
            // 3. P_det CROSSING — usa n_crossing/n_laser
            //    Override temporaneo p_det con crossing per questa canvas
            {
                std::vector<MapPoint> ptsCross = mapPts;
                for (auto& mp : ptsCross)
                    if (mp.ok && mp.n_laser > 0)
                        mp.p_det = static_cast<double>(mp.n_crossing) / mp.n_laser;
                drawMap2D(ptsCross, "p_det", "P_{det} crossing (%)",
                          vbias, frac_pe, ctx,
                          Form("map_pdet_crossing_%s.png", tagBase.c_str()),
                          map_cx, map_cy);
                drawMap2D(ptsCross, "n_crossing", "N_{crossing}",
                          vbias, frac_pe, ctx,
                          Form("map_ncross_crossing_%s.png", tagBase.c_str()),
                          map_cx, map_cy);
            }
            // 5. P_det ACCEPTED — usa n_sig/n_laser (solo punti con S+B ok)
            //    FIX: era n_acc/n_laser che e' segnale+fondo, non solo segnale
            {
                std::vector<MapPoint> ptsAcc = mapPts;
                for (auto& mp : ptsAcc)
                    if (mp.ok && mp.n_laser > 0)
                        mp.p_det = mp.mu_reliable
                            ? std::min(1.0, mp.n_sig / static_cast<double>(mp.n_laser))
                            : static_cast<double>(mp.n_crossing) / mp.n_laser;
                drawMap2D(ptsAcc, "p_det", "P_{det} accepted (%)",
                          vbias, frac_pe, ctx,
                          Form("map_pdet_accepted_%s.png", tagBase.c_str()),
                          map_cx, map_cy);
                drawMap2D(ptsAcc, "n_acc", "N_{accepted}",
                          vbias, frac_pe, ctx,
                          Form("map_nacc_accepted_%s.png", tagBase.c_str()),
                          map_cx, map_cy);
            }
            // 7. Tabella TXT + canvas grafica
            saveMapTable(mapPts, vbias, frac_pe, dataDir, ctx, pdet_mode, map_cx, map_cy);

            // Summary table
            int nUnreliable = 0;
            for (const auto& mp : mapPts) if (mp.ok && !mp.mu_reliable) ++nUnreliable;
            std::cout << "\n+-- SUMMARY  Vbias=" << vbias
                      << " V  frac=" << frac_pe << " pe"
                      << "  P_det=" << (pdet_mode == PdetMode::CROSSING ? "CROSSING" : "ACCEPTED")
                      << "  centre=(" << map_cx << "," << map_cy << ")"
                      << (nUnreliable > 0 ? Form("  [*=%d mu unreliable]", nUnreliable) : "")
                      << " --\n"
                      << std::setw(10) << "dx(mm)"
                      << std::setw(10) << "dy(mm)"
                      << std::setw(12) << "mu(ns)"
                      << std::setw(11) << "sig(ns)"
                      << std::setw(10) << "N_acc"
                      << std::setw(10) << "N_cross"
                      << std::setw(9)  << "P_det%"
                      << std::setw(7)  << "mu_ok" << "\n"
                      << std::string(79, '-') << "\n";
            for (const auto& mp : mapPts) {
                if (!mp.ok) continue;
                std::cout << std::fixed << std::setprecision(1)
                          << std::setw(10) << (mp.x - map_cx)
                          << std::setw(10) << (mp.y - map_cy)
                          << std::setprecision(3)
                          << std::setw(12) << mp.mu
                          << std::setw(11) << mp.sigma
                          << std::setw(10) << mp.n_acc
                          << std::setw(10) << mp.n_crossing
                          << std::setprecision(1)
                          << std::setw(9)
                          << (mp.p_det >= 0 ? mp.p_det * 100.0 : -1.0)
                          << std::setw(7) << (mp.mu_reliable ? "OK" : "*FAIL")
                          << "\n";
            }
        } // end frac_pe loop
    } // end vbias loop

end_outer:
    std::cout << "\n+==========================================================+\n"
              << "|  POSITION SCAN DONE  "
              << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - t0).count()
              << " s\n"
              << "|  P_det mode : "
              << (pdet_mode == PdetMode::CROSSING ? "CROSSING" : "ACCEPTED") << "\n"
              << "|  Results    : " << dataDir << "/map_results_*.root\n"
              << "|  Canvas     : " << ctx.pngDir << "\n"
              << "|  Next step  : sipm_draw_timing3d()\n"
              << "+==========================================================+\n";

    ctx.reopenAllCanvases(20);

    } catch (const std::runtime_error& e) {
        std::cerr << "\n[ABORT] " << e.what() << "\n";
    }
}
