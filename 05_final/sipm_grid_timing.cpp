/**
 * sipm_grid_timing.cpp – ROBUSTO e COMPLETO
 * Analisi spaziale del tempo di arrivo con TUTTE le logiche di selezione
 * fisica usate in tot_analysis (baseline RMS, preLaserQuiet, edge check,
 * isteresi, pre‑check, falling‑edge confirmation, post‑check, tail check).
 * Scritto da zero, senza toccare gli header.
 *
 * Uso:
 *   .L sipm_grid_timing.cpp+
 *   sipm_grid_timing("percorso/dati", 0.5, 200, true, 95, 125)
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/CalibIO.h"
#include "../header/ProbabilityAnalysis.h"

static constexpr int N_SAMPLES = 1024;

// ── Copia locale delle funzioni di selezione ──────────────────────────────

static double edgeMedian(const std::vector<double>& amp, int i0, int n_pts) {
    int sz = (int)amp.size();
    int i1 = std::max(0, i0);
    int i2 = std::min(sz, i0 + n_pts);
    if (i2 <= i1) return 0.0;
    std::vector<double> tmp(amp.begin() + i1, amp.begin() + i2);
    int mid = (int)tmp.size() / 2;
    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
    return tmp[mid];
}

static std::vector<double> correctBaselineRMS(const std::vector<double>& time,
                                              const std::vector<double>& amp,
                                              double t_base_start,
                                              double t_base_end,
                                              double max_rms,
                                              bool* baseline_ok) {
    *baseline_ok = true;
    std::vector<double> pre;
    for (size_t i = 0; i < time.size(); ++i)
        if (time[i] >= t_base_start && time[i] < t_base_end)
            pre.push_back(amp[i]);
    if (pre.empty()) { *baseline_ok = false; return amp; }
    std::vector<double> tmp = pre;
    size_t mid = tmp.size() / 2;
    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
    double offset = tmp[mid];
    if (max_rms > 0.0) {
        double sum2 = 0.0;
        for (double v : pre) sum2 += (v - offset) * (v - offset);
        double rms = std::sqrt(sum2 / pre.size());
        if (rms > max_rms) { *baseline_ok = false; return amp; }
    }
    std::vector<double> out = amp;
    for (auto& a : out) a -= offset;
    return out;
}

static std::pair<double,double> computeTOT_local(
        const std::vector<double>& time,
        const std::vector<double>& amp,
        double threshold,
        int    j_start,
        int    j_end,
        int    n_edge           = 100,
        double edge_thr_frac    = 1000.0,
        double hyst_frac        = 0.5,
        int    pre_check_n      = 10,
        double pre_check_frac   = 0.3,
        int    post_check_n     = 30,
        int    confirm_window   = 50,
        int    min_below        = 10)
{
    if (j_start < 0 || j_end >= (int)amp.size() || j_start >= j_end)
        return {-1.0, -1.0};

    if (edge_thr_frac < 999.0) {
        double edge_limit = edge_thr_frac * threshold;
        if (edgeMedian(amp, j_start, n_edge) >= edge_limit) return {-1.0,-1.0};
    }

    double thr_lo = threshold * hyst_frac;
    bool   armed  = true;
    double t_rise = -1.0;
    int    j_rise = -1;

    for (int j = j_start; j <= j_end; ++j) {
        if (!armed) {
            if (amp[j] < thr_lo) armed = true;
            continue;
        }
        if (amp[j] >= threshold) {
            if (pre_check_n > 0) {
                double pre = edgeMedian(amp, j - pre_check_n, pre_check_n);
                if (pre >= pre_check_frac * threshold) { armed = false; continue; }
            }
            if (j > j_start && amp[j-1] < threshold)
                t_rise = time[j-1] + (threshold - amp[j-1])
                         * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
            else
                t_rise = time[j];
            j_rise = j;
            break;
        }
    }
    if (t_rise < 0) return {-1.0, -1.0};

    double t_fall = -1.0;
    for (int j = j_rise + 1; j <= j_end; ++j) {
        if (amp[j] < threshold) {
            double t_cand = (j > 0 && amp[j-1] >= threshold)
                ? time[j-1] + (threshold - amp[j-1])
                  * (time[j] - time[j-1]) / (amp[j] - amp[j-1])
                : time[j];
            int j_win_end = std::min(j + confirm_window, j_end);
            int n_sub = 0;
            for (int k = j; k <= j_win_end; ++k)
                if (amp[k] < threshold) ++n_sub;
            if (n_sub >= min_below) { t_fall = t_cand; break; }
            j = j_win_end;
        }
    }
    if (t_fall < 0) return {-1.0, -1.0};

    if (post_check_n > 0) {
        int jf = j_rise;
        while (jf <= j_end && time[jf] < t_fall) ++jf;
        if (jf <= j_end && edgeMedian(amp, jf, post_check_n) >= threshold)
            return {-1.0, -1.0};
    }

    int tail_n     = 50;
    int tail_start = std::max(j_start, j_end + 1 - tail_n);
    int tail_len   = j_end + 1 - tail_start;
    if (tail_len > 0 && edgeMedian(amp, tail_start, tail_len) >= thr_lo)
        return {-1.0, -1.0};

    return {t_rise, t_fall};
}

static bool preLaserQuiet(const std::vector<double>& time,
                           const std::vector<double>& amp,
                           double t_laser, double quiet_thr,
                           double margin_ns = 1.0)
{
    double t_check_end = t_laser - margin_ns;
    if (t_check_end <= BASELINE_END) return true;
    for (size_t j = 0; j < time.size(); ++j) {
        if (time[j] < BASELINE_END) continue;
        if (time[j] > t_check_end) break;
        if (amp[j]  > quiet_thr)   return false;
    }
    return true;
}

// ═════════════════════════════════════════════════════════════════
//  Processa un singolo file con TUTTA la selezione
// ═════════════════════════════════════════════════════════════════
static void processOneFile_robust(const std::string& path,
                                  double cutoff_MHz,
                                  double fs_MHz,
                                  int    j_trig_start,
                                  int    j_trig_end,
                                  double let_thr,
                                  double laser_thr,
                                  double baseline_max_rms,
                                  double& mu_out, double& sigma_out,
                                  long& n_events_out)
{
    mu_out = sigma_out = 0.0; n_events_out = 0;
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return; }
    TTree* tCh1   = (TTree*)f->Get("ch1");
    TTree* tLaser = (TTree*)f->Get("laser");
    if (!tCh1 || !tLaser) { f->Close(); delete f; return; }

    Double_t t1[N_SAMPLES], a1[N_SAMPLES], tL[N_SAMPLES], aL[N_SAMPLES];
    tCh1  ->SetBranchAddress("time", t1);
    tCh1  ->SetBranchAddress("amplitude", a1);
    tLaser->SetBranchAddress("time", tL);
    tLaser->SetBranchAddress("amplitude", aL);
    tCh1->SetCacheSize(2 * 1024 * 1024);
    tLaser->SetCacheSize(2 * 1024 * 1024);

    Long64_t nFile  = tCh1->GetEntries();
    Long64_t nLaser = tLaser->GetEntries();
    Long64_t nEv    = std::min(nFile, nLaser);

    std::vector<double> pre_dt;
    for (Long64_t i = 0; i < nEv && pre_dt.size() < 2000; ++i) {
        tCh1->GetEntry(i); tLaser->GetEntry(i);
        double t_laser = laserTriggerTime(tL, aL, N_SAMPLES, laser_thr);
        if (t_laser < -900) continue;
        std::vector<double> vt(t1, t1+N_SAMPLES), va(a1, a1+N_SAMPLES);
        bool bl_ok;
        auto af = correctBaselineRMS(vt, va, BASELINE_START, BASELINE_END,
                                      baseline_max_rms, &bl_ok);
        if (!bl_ok) continue;
        if (cutoff_MHz > 0) af = butterworthLowPass(af, cutoff_MHz, fs_MHz);
        if (!preLaserQuiet(vt, af, t_laser, let_thr)) continue;
        auto tot = computeTOT_local(vt, af, let_thr, j_trig_start, j_trig_end,
                                    100, 1000.0, 0.5, 10, 0.3, 30, 50, 10);
        if (tot.first < 0 || tot.second < 0) continue;
        double dt = tot.first - t_laser;
        pre_dt.push_back(dt);
    }
    if (pre_dt.empty()) { f->Close(); delete f; return; }
    std::nth_element(pre_dt.begin(), pre_dt.begin() + pre_dt.size()/2, pre_dt.end());
    double median_dt = pre_dt[pre_dt.size()/2];
    double hlo = median_dt - 20.0, hhi = median_dt + 20.0;
    if (hlo < -50) hlo = -50; if (hhi > 150) hhi = 150;

    TH1D* h = new TH1D("h_dt", ";#Deltat (ns);Events", 800, hlo, hhi);
    h->SetDirectory(nullptr); h->Sumw2();

    for (Long64_t i = 0; i < nEv; ++i) {
        tCh1->GetEntry(i); tLaser->GetEntry(i);
        double t_laser = laserTriggerTime(tL, aL, N_SAMPLES, laser_thr);
        if (t_laser < -900) continue;
        std::vector<double> vt2(t1, t1+N_SAMPLES), va2(a1, a1+N_SAMPLES);
        bool bl2;
        auto af2 = correctBaselineRMS(vt2, va2, BASELINE_START, BASELINE_END,
                                       baseline_max_rms, &bl2);
        if (!bl2) continue;
        if (cutoff_MHz > 0) af2 = butterworthLowPass(af2, cutoff_MHz, fs_MHz);
        if (!preLaserQuiet(vt2, af2, t_laser, let_thr)) continue;
        auto tot2 = computeTOT_local(vt2, af2, let_thr, j_trig_start, j_trig_end,
                                     100, 1000.0, 0.5, 10, 0.3, 30, 50, 10);
        if (tot2.first < 0 || tot2.second < 0) continue;
        double dt = tot2.first - t_laser;
        if (dt >= hlo && dt <= hhi) h->Fill(dt);
    }
    f->Close(); delete f;

    n_events_out = (long)h->GetEntries();
    if (n_events_out < 30) { delete h; return; }

    int bMax = h->GetMaximumBin();
    double peak_pos = h->GetBinCenter(bMax);
    double peak_val = h->GetBinContent(bMax);
    double sig_est  = h->GetRMS(); if (sig_est < 0.1) sig_est = 0.5;

    TF1* fg = new TF1("fg", "gaus", peak_pos - 2.5*sig_est, peak_pos + 2.5*sig_est);
    fg->SetParameters(peak_val, peak_pos, sig_est);
    fg->SetParLimits(2, 0.05, 10.0);
    int status = h->Fit(fg, "QNR");
    if (status == 0 || status == 4000) {
        mu_out    = fg->GetParameter(1);
        sigma_out = std::abs(fg->GetParameter(2));
    } else {
        double sw=0,swx=0,swx2=0;
        for (int b=1; b<=h->GetNbinsX(); ++b) {
            double cnt=h->GetBinContent(b), xc=h->GetBinCenter(b);
            sw+=cnt; swx+=cnt*xc; swx2+=cnt*xc*xc;
        }
        if (sw>0) { mu_out=swx/sw; sigma_out=std::sqrt(std::max(swx2/sw - mu_out*mu_out,0.0)); }
    }
    delete fg; delete h;
}

// ═════════════════════════════════════════════════════════════════
//  MAIN
// ═════════════════════════════════════════════════════════════════
void sipm_grid_timing(std::string dataDir         = "",
                      double frac_pe              = 0.5,
                      double cutoff_MHz           = 0.0,
                      bool   use_filter           = false,
                      double t_trig_start         = 0.0,
                      double t_trig_end           = 204.6)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    if (dataDir.empty()) {
        std::cout << "Data directory (default " << DATA_DIR << "): ";
        std::getline(std::cin, dataDir);
        if (dataDir.empty()) dataDir = DATA_DIR;
    }
    while (!dataDir.empty() && dataDir.back() == '/') dataDir.pop_back();
    g_data_dir_override = dataDir;

    // Raccoglie tutti i file data_x_…_y_…_vbias_….root
    std::vector<std::pair<std::string, std::tuple<double,double,int>>> posFiles;
    void* dp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dp) { std::cerr << "Directory non trovata\n"; return; }
    const char* entry;
    std::regex re(R"(^data_x_([\-\d]+)_y_([\-\d]+)_vbias_(\d+)\.root$)");
    while ((entry = gSystem->GetDirEntry(dp))) {
        std::string fn(entry);
        std::smatch m;
        if (std::regex_match(fn, m, re))
            posFiles.push_back({dataDir+"/"+fn,
                                {std::stod(m[1]), std::stod(m[2]), std::stoi(m[3])}});
    }
    gSystem->FreeDirectory(dp);
    if (posFiles.empty()) { std::cerr << "Nessun file data_x_... trovato\n"; return; }

    // Vbias unici
    std::set<int> vbiasSet;
    for (auto& [path, pos] : posFiles) vbiasSet.insert(std::get<2>(pos));

    OutCtx ctx = createOutputDirs("grid_timing");
    double fs_MHz = 5000.0;

    for (int vbias : vbiasSet) {
        // Carica calibrazione
        CalibResult cal;
        if (!loadCalibration(cal, vbias, cutoff_MHz, dataDir)) {
            void* d2 = gSystem->OpenDirectory(dataDir.c_str());
            if (d2) {
                std::regex calre("calib_vbias" + std::to_string(vbias) + "_cut(\\d+)mhz\\.root");
                const char* e2;
                while ((e2 = gSystem->GetDirEntry(d2))) {
                    std::string cf(e2);
                    std::smatch cm;
                    if (std::regex_match(cf, cm, calre)) {
                        double cc = std::stod(cm[1]);
                        if (loadCalibration(cal, vbias, cc, dataDir)) break;
                    }
                }
                gSystem->FreeDirectory(d2);
            }
        }
        if (!cal.ok) {
            std::cerr << "Nessuna calibrazione per vbias " << vbias << "\n";
            continue;
        }
        double let_thr = cal.q + frac_pe * cal.m;
        double baseline_max_rms = use_filter ? 2.0 : 4.0;

        // Determina j_start e j_end dal primo file valido per questo vbias
        int j_start = 0, j_end = N_SAMPLES - 1;
        bool found_trigger_idx = false;
        for (auto& [path, pos] : posFiles) {
            if (std::get<2>(pos) != vbias) continue;
            TFile* ff = TFile::Open(path.c_str());
            if (!ff) continue;
            TTree* tt = (TTree*)ff->Get("ch1");
            if (tt && tt->GetEntries() > 0) {
                Double_t t_arr[N_SAMPLES];
                tt->SetBranchAddress("time", t_arr);
                tt->GetEntry(0);
                triggerWindowIndices(t_arr, N_SAMPLES, t_trig_start, t_trig_end,
                                     j_start, j_end, std::to_string(vbias));
                ff->Close(); delete ff;
                found_trigger_idx = true;
                break;
            }
            ff->Close(); delete ff;
        }
        if (!found_trigger_idx) {
            std::cerr << "  [WARN] Impossibile determinare indici trigger per Vbias "
                      << vbias << " – salto.\n";
            continue;
        }

        // Colleziona le coordinate x e y per questo vbias
        std::set<double> xs, ys;
        for (auto& [path, pos] : posFiles)
            if (std::get<2>(pos) == vbias) {
                xs.insert(std::get<0>(pos));
                ys.insert(std::get<1>(pos));
            }
        std::vector<double> xv(xs.begin(), xs.end()), yv(ys.begin(), ys.end());
        int nx = xv.size(), ny = yv.size();
        if (nx == 0 || ny == 0) continue;

        auto edges = [](const std::vector<double>& v) {
            std::vector<double> e;
            e.reserve(v.size()+1);
            for (size_t i=0; i<v.size(); ++i) {
                double lo = (i==0) ? v[0]-5.0 : 0.5*(v[i-1]+v[i]);
                e.push_back(lo);
            }
            e.push_back(v.back()+5.0);
            return e;
        };
        auto xEdges = edges(xv), yEdges = edges(yv);

        TH2D* hMu = new TH2D(Form("hMu_v%d",vbias),
            Form("<#Deltat> map  Vbias=%d V;x (mm);y (mm);<#Deltat> (ns)",vbias),
            nx, xEdges.data(), ny, yEdges.data());
        TH2D* hSigma = new TH2D(Form("hSig_v%d",vbias),
            Form("#sigma map  Vbias=%d V;x (mm);y (mm);#sigma (ns)",vbias),
            nx, xEdges.data(), ny, yEdges.data());
        TH2D* hNev = new TH2D(Form("hNev_v%d",vbias),
            Form("N_{acc} map  Vbias=%d V;x (mm);y (mm);N",vbias),
            nx, xEdges.data(), ny, yEdges.data());
        hMu->SetDirectory(0); hSigma->SetDirectory(0); hNev->SetDirectory(0);

        int nDone = 0;
        std::cout << "Vbias = " << vbias << " V, LET = " << frac_pe << " p.e.\n";
        for (auto& [path, pos] : posFiles) {
            if (std::get<2>(pos) != vbias) continue;
            double x = std::get<0>(pos), y = std::get<1>(pos);
            double mu, sigma; long nEv;
            processOneFile_robust(path, cutoff_MHz, fs_MHz,
                                  j_start, j_end, let_thr, cal.laser_thr,
                                  baseline_max_rms, mu, sigma, nEv);
            if (nEv < 30) {
                std::cout << "  pos ("<<x<<","<<y<<") scartata, N="<<nEv<<"\n";
                continue;
            }
            int bx = hMu->GetXaxis()->FindBin(x);
            int by = hMu->GetYaxis()->FindBin(y);
            hMu  ->SetBinContent(bx, by, mu);
            hSigma->SetBinContent(bx, by, sigma);
            hNev  ->SetBinContent(bx, by, (double)nEv);
            std::cout << "  ("<<x<<","<<y<<") mu="<<mu<<" sigma="<<sigma<<" N="<<nEv<<"\n";
            ++nDone;
        }

        for (auto* h : {hMu, hSigma, hNev}) {
            TCanvas* c = new TCanvas(h->GetName(), h->GetTitle(), 800, 600);
            c->SetRightMargin(0.15);
            h->Draw("COLZ TEXT");
            ctx.savePNG(c, std::string(h->GetName()) + ".png");
        }
        delete hMu; delete hSigma; delete hNev;
    }
    std::cout << "\nGrafici salvati in " << ctx.pngDir << "\n";

    // ── Analisi probabilità (metodo crossing, senza filtri TOT) ──────────────
    // Costruisce cal_map da tutti i Vbias già caricati nel loop precedente
    {
        std::map<int, CalibResult> cal_map;
        for (int vb : vbiasSet) {
            CalibResult cal;
            if (loadCalibration(cal, vb, cutoff_MHz, dataDir)) {
                cal_map[vb] = cal;
            } else {
                // Prova a trovare qualsiasi cutoff disponibile
                void* d3 = gSystem->OpenDirectory(dataDir.c_str());
                if (d3) {
                    std::regex calre("calib_vbias" + std::to_string(vb) + "_cut(\\d+)mhz\\.root");
                    const char* e3;
                    while ((e3 = gSystem->GetDirEntry(d3))) {
                        std::string cf(e3);
                        std::smatch cm;
                        if (std::regex_match(cf, cm, calre)) {
                            double cc = std::stod(cm[1]);
                            if (loadCalibration(cal, vb, cc, dataDir)) {
                                cal_map[vb] = cal;
                                break;
                            }
                        }
                    }
                    gSystem->FreeDirectory(d3);
                }
            }
        }

        runProbabilityAnalysis(posFiles, vbiasSet, cal_map,
                               frac_pe, cutoff_MHz, fs_MHz, use_filter,
                               t_trig_start, t_trig_end, ctx);
    }
}