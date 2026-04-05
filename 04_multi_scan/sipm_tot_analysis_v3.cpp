// sipm_tot_analysis_v3.cpp
// Analisi TOT completa con calibrazione tramite threshold scan.
//
// Output:
//   macro/canvas/YYYYMMDD_HHMMSS/    ← tutti i PNG in ordine logico
//   macro/file_root/YYYYMMDD_HHMMSS/ ← tutti i ROOT
//
// Ordine numerico dei PNG:
//   01  calibration_<tag>.png
//   ── per ogni file × ogni LET ─────────────────────────────
//   02  tot_map_<tag>.png             mappa 2D raw
//   03  tot_projection_<tag>.png      proiezione Δt raw
//   04  tw_empirical/amplitude_...    time walk diagnostico
//   05  tot_map_corr_<tag>.png        mappa 2D corretta
//   06  tot_projection_corr_<tag>.png proiezione Δt corretta
//   07  tot_slices_<tag>.png          slice in TOT raw
//   08  tot_trend_<tag>.png           σ e Mean vs TOT raw
//   09  tot_slices_corr_<tag>.png     slice in TOT corrette
//   10  tot_trend_corr_<tag>.png      σ e Mean vs TOT corrette
//   11  tot_byPE_<tag>.png            Δt per p.e.
//   12  tot_sigma_vs_pe_<tag>.png     σ vs n_pe
//   ── al termine ───────────────────────────────────────────
//   13  tot_multiLET_<tag>.png        overlay multi-LET

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/Calibration.h"
#include "../header/TOTAnalysis.h"
#include "../header/TimingCorrection.h"
#include "../header/TOTPlotting.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

// ════════════════════════════════════════════════════════════
//  OVERLAY MULTI-LET  (PNG 15)
// ════════════════════════════════════════════════════════════
static void drawMultiLETOverlay(
        const std::vector<double>& fracs_pe,
        const std::vector<TH2D*>&  maps,
        const std::string&         tag,
        double fit_lo, double fit_hi,
        OutCtx& ctx) {

    if (maps.empty()) return;
    static const int cols[] = {kAzure+1, kRed+1, kGreen+2,
                                kOrange+7, kMagenta+1, kCyan+2};
    int nc = (int)(sizeof(cols)/sizeof(cols[0]));

    TCanvas* cOvl = new TCanvas(Form("cMultiLET_%s", tag.c_str()),
        Form("Multi-LET #Deltat overlay — %s", tag.c_str()), 900, 600);
    cOvl->SetGrid();
    cOvl->SetLeftMargin(PAD_LEFT);   cOvl->SetRightMargin(PAD_RIGHT);
    cOvl->SetBottomMargin(PAD_BOTTOM); cOvl->SetTopMargin(PAD_TOP);

    TLegend* leg = new TLegend(0.62, 0.62, 0.93, 0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetFillColor(0);
    leg->SetTextFont(42);  leg->SetTextSize(0.036);

    bool first = true;
    for (int i = 0; i < (int)maps.size(); ++i) {
        if (!maps[i]) continue;
        TH1D* hProj = maps[i]->ProjectionY(Form("hOvl_%s_%d",tag.c_str(),i));
        hProj->SetDirectory(nullptr);
        int b1=hProj->FindBin(fit_lo), b2=hProj->FindBin(fit_hi);
        double maxVal=0;
        for (int b=b1; b<=b2; ++b)
            if (hProj->GetBinContent(b)>maxVal) maxVal=hProj->GetBinContent(b);
        if (maxVal>0) hProj->Scale(1.0/maxVal);
        hProj->SetLineColor(cols[i%nc]); hProj->SetLineWidth(2);
        hProj->GetXaxis()->SetRangeUser(fit_lo-5.0, fit_hi+5.0);
        hProj->SetTitle(Form("Multi-LET #Deltat overlay (corrected)   %s"
                             ";#Deltat (ns);Normalised counts", tag.c_str()));
        if (first) { hProj->Draw("HIST"); first=false; }
        else        hProj->Draw("HIST same");
        leg->AddEntry(hProj, Form("LET = %.2f p.e.", fracs_pe[i]), "l");
    }
    leg->Draw();
    cOvl->Update(); cOvl->Modified();
    ctx.savePNG(cOvl, Form("tot_multiLET_%s.png", tag.c_str()));
}

// ════════════════════════════════════════════════════════════
//  PROCESS ONE LET  (PNG 04–14, ROOT)
// ════════════════════════════════════════════════════════════
static TH2D* processOneLET(const FileInfo&    info,
                             double             cutoff_MHz,
                             double             t_trig_start,
                             double             t_trig_end,
                             double             fit_lo,
                             double             fit_hi,
                             double             frac_pe,
                             const CalibResult& cal,
                             TWMethod           tw_method,
                             bool               do_pe_analysis,
                             OutCtx&            ctx) {

    TFile* file = TFile::Open(info.path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[SKIP] Cannot open: " << info.path << "\n"; return nullptr;
    }
    TTree* treeCh1   = (TTree*)file->Get("ch1");
    TTree* treeLaser = (TTree*)file->Get("laser");
    if (!treeCh1 || !treeLaser) {
        std::cerr << "[SKIP] Missing trees in: " << info.path << "\n";
        file->Close(); return nullptr;
    }

    const int N = 1024; Double_t t1[N];
    treeCh1->SetBranchAddress("time", t1); treeCh1->GetEntry(0);
    double fs_MHz = 1000.0 / (t1[1] - t1[0]);

    int j_trig_start, j_trig_end;
    triggerWindowIndices(t1, N, t_trig_start, t_trig_end,
                         j_trig_start, j_trig_end, info.tag);

    double let_thr = cal.q + frac_pe * cal.m;
    std::cout << "\n  LET = " << frac_pe << " p.e."
              << "   threshold = " << let_thr << " mV\n";

    auto events = collectTOTEvents(treeCh1, treeLaser,
                                   cutoff_MHz, fs_MHz,
                                   j_trig_start, j_trig_end,
                                   let_thr, cal, info.tag);
    file->Close();

    if (events.empty()) { std::cout << "  No events — skipping.\n"; return nullptr; }

    std::string ltag = Form("%s_let%.2fpe", info.tag.c_str(), frac_pe);

    // 04: mappa raw
    TH2D* h2D = fillTH2D(events, ltag, frac_pe, false);
    drawTOTMap(h2D, ltag, ctx, false);

    // 05: proiezione raw
    drawGlobalProjection(h2D, ltag, fit_lo, fit_hi, ctx, false);

    // 06 (time walk) + 07 mappa corr + 08 proiezione corr
    auto events_corr = applyTimeWalkCorrection(
        events, tw_method, fit_lo, fit_hi, ltag, ctx);

    TH2D* h2D_corr = nullptr;
    if (tw_method != TWMethod::NONE) {
        h2D_corr = fillTH2D(events_corr, ltag, frac_pe, true);
        drawTOTMap(h2D_corr, ltag, ctx, true);
        drawGlobalProjection(h2D_corr, ltag, fit_lo, fit_hi, ctx, true);
    }

    // 09-10: slice e trend raw
    drawSlicesAndTrend(h2D, ltag, fit_lo, fit_hi, frac_pe, ctx, false);

    // 11-12: slice e trend corretti
    if (tw_method != TWMethod::NONE)
        drawSlicesAndTrend(h2D_corr, ltag, fit_lo, fit_hi, frac_pe, ctx, true);

    // 13-14: per p.e.
    if (do_pe_analysis)
        drawByPE((tw_method != TWMethod::NONE) ? events_corr : events,
                 ltag, fit_lo, fit_hi, ctx);

    // ROOT
    ctx.saveRoot(h2D, h2D_corr, ltag, frac_pe);

    return (h2D_corr != nullptr) ? h2D_corr : h2D;
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════
void sipm_tot_analysis_v3() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // Crea subito le cartelle di output
    OutCtx ctx = createOutputDirs();

    std::vector<FileInfo> files = findInputFiles();
    if (files.empty()) {
        std::cerr << "[ERROR] No data.vbias_{V}_L.root files found.\n"; return;
    }
    printFileList(files);

    int calIdx = selectOneFile(files,
        "\nSelect calibration file index (brightest): ");
    std::vector<int> anaIdx = selectMultipleFiles(files,
        "Select analysis file indices (e.g. \"1\" or \"0,1\" or \"all\"): ");

    double cutoff_MHz;
    std::cout << "\nLow-pass filter cut-off [MHz]: ";
    std::cin >> cutoff_MHz;

    double t_wf_min=0, t_wf_max=0;
    {
        TFile* f0 = TFile::Open(files[calIdx].path.c_str(),"READ");
        if (f0 && !f0->IsZombie()) {
            TTree* tr=(TTree*)f0->Get("ch1");
            if (tr) {
                const int N0=1024; Double_t t0[N0];
                tr->SetBranchAddress("time",t0); tr->GetEntry(0);
                t_wf_min=t0[0]; t_wf_max=t0[N0-1];
            }
            f0->Close();
        }
    }

    double t_trig_start, t_trig_end;
    std::cout << "\n--- Trigger window ---\n"
              << "  Waveform: [" << t_wf_min << ", " << t_wf_max << "] ns\n"
              << "  Baseline: [" << BASELINE_START << ", " << BASELINE_END << ") ns\n"
              << "  Start [ns]: "; std::cin >> t_trig_start;
    std::cout << "  End   [ns]: "; std::cin >> t_trig_end;
    if (t_trig_start >= t_trig_end) {
        std::cerr << "[ERROR] Invalid trigger window.\n"; return;
    }

    double fit_lo, fit_hi;
    std::cout << "\n--- Fit window (Δt peak) ---\n"
              << "  Start [ns]: "; std::cin >> fit_lo;
    std::cout << "  End   [ns]: "; std::cin >> fit_hi;
    if (fit_lo >= fit_hi) { fit_lo=t_trig_start; fit_hi=t_trig_end; }

    std::vector<double> fracs_pe;
    std::cout << "\n--- LET thresholds (p.e., space-separated) ---\n  > ";
    {
        std::string line; std::cin.ignore(); std::getline(std::cin, line);
        std::stringstream ss(line); double v;
        while (ss >> v) if (v>0) fracs_pe.push_back(v);
    }
    if (fracs_pe.empty()) fracs_pe.push_back(1.0);

    TWMethod tw_method = askTimeWalkMethod();

    char ans_pe=0;
    while (ans_pe!='y' && ans_pe!='n') {
        std::cout << "\nPer-p.e. timing analysis? [y/n]: "; std::cin >> ans_pe;
    }
    bool do_pe_analysis = (ans_pe=='y');

    std::cout << "\n╔══════════════════════════════════════╗\n"
              << "║         Settings summary             ║\n"
              << "╚══════════════════════════════════════╝\n"
              << "  LP cut-off      : " << cutoff_MHz << " MHz\n"
              << "  Trigger window  : [" << t_trig_start << ", " << t_trig_end << "] ns\n"
              << "  Fit window      : [" << fit_lo << ", " << fit_hi << "] ns\n"
              << "  LET thresholds  : ";
    for (double f : fracs_pe) std::cout << f << " ";
    std::cout << "p.e.\n"
              << "  Time walk corr  : "
              << (tw_method==TWMethod::EMPIRICAL ? "Empirical" :
                  tw_method==TWMethod::AMPLITUDE ? "Amplitude" : "None")
              << "\n  Per-p.e. analysis: " << (do_pe_analysis?"yes":"no") << "\n\n";

    // ── Calibrazione  →  PNG 01, 02, 03 ─────────────────────
    std::cout << "══════ Calibration ══════\n";
    TFile* fCal = TFile::Open(files[calIdx].path.c_str(), "READ");
    if (!fCal || fCal->IsZombie()) {
        std::cerr << "[ERROR] Cannot open calibration file.\n"; return;
    }
    TTree* tCal = (TTree*)fCal->Get("ch1");
    if (!tCal) { fCal->Close(); std::cerr << "[ERROR] No ch1 tree.\n"; return; }
    const int NC=1024; Double_t tc[NC];
    tCal->SetBranchAddress("time",tc); tCal->GetEntry(0);
    double fs_cal = 1000.0/(tc[1]-tc[0]);

    CalibResult cal = calibrateSpectrum(tCal, cutoff_MHz, fs_cal,
                                         files[calIdx].tag, ctx);
    fCal->Close();
    if (!cal.ok) { std::cerr << "[ERROR] Calibration failed.\n"; return; }

    // ── Analisi  →  PNG 04…15 per LET ───────────────────────
    for (int idx : anaIdx) {
        std::cout << "\n══════ " << files[idx].tag << " ══════\n";
        std::vector<TH2D*> maps;
        for (double frac : fracs_pe) {
            TH2D* h = processOneLET(files[idx], cutoff_MHz,
                                    t_trig_start, t_trig_end,
                                    fit_lo, fit_hi,
                                    frac, cal, tw_method, do_pe_analysis, ctx);
            maps.push_back(h);
        }
        if (fracs_pe.size() > 1)
            drawMultiLETOverlay(fracs_pe, maps, files[idx].tag,
                                fit_lo, fit_hi, ctx);
    }

    std::cout << "\n=== Done ===\n"
              << "  PNG  → " << ctx.pngDir  << "\n"
              << "  ROOT → " << ctx.rootDir << "\n";
}
