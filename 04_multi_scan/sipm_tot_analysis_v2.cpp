// sipm_tot_analysis_v2.cpp
// Analisi TOT completa con calibrazione incrociata.
//
// Funzionalità:
//   1. Calibrazione spettro p.e. da file luminoso (cross-calibrazione)
//   2. Raccolta eventi TOT con stima numero p.e.
//   3. Correzione time walk (empirica o per ampiezza — scelta operatore)
//   4. Mappa 2D TOT vs Δt (raw e corretta)
//   5. Proiezione globale sull'asse Δt con fit gaussiano
//   6. Slice in TOT con fit q-gaussiano
//   7. Distribuzioni Δt per numero di p.e.
//   8. Analisi multi-LET: lista di soglie, canvas overlay
//
// Header:
//   Config.h          — parametri globali, FileInfo, scansione file
//   SignalProcessing.h — baseline, filtro, trigger laser
//   Calibration.h     — calibrateSpectrum, selezione file
//   TOTAnalysis.h     — computeTOT, collectTOTEvents, estimatePE
//   TimingCorrection.h — correzione time walk (metodo A e B)
//   TOTPlotting.h     — tutte le canvas

#include "../header/Config.h"
#include "../header/SignalProcessing.h"
#include "../header/Calibration.h"
#include "../header/TOTAnalysis.h"
#include "../header/TimingCorrection.h"
#include "../header/TOTPlotting.h"

#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>

// ════════════════════════════════════════════════════════════
//  CANVAS OVERLAY MULTI-LET
//  Sovrappone le mappe 2D (proiezione sull'asse Δt) per diversi
//  valori di LET, per confrontare la risoluzione temporale.
// ════════════════════════════════════════════════════════════
static void drawMultiLETOverlay(
        const std::vector<double>& fracs_pe,
        const std::vector<TH2D*>& maps,
        const std::string& tag,
        double fit_lo, double fit_hi) {

    if (maps.empty()) return;

    // Palette colori
    static const int cols[] = {kAzure+1, kRed+1, kGreen+2,
                                kOrange+7, kMagenta+1, kCyan+2};
    int nc = sizeof(cols)/sizeof(cols[0]);

    TCanvas* cOvl = new TCanvas(Form("cMultiLET_%s", tag.c_str()),
        Form("Multi-LET overlay — %s", tag.c_str()), 900, 600);
    cOvl->SetGrid(); cOvl->SetLeftMargin(PAD_LEFT);
    cOvl->SetBottomMargin(PAD_BOTTOM); cOvl->SetTopMargin(PAD_TOP);

    TLegend* leg = new TLegend(0.65, 0.65, 0.93, 0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001);
    leg->SetTextFont(42); leg->SetTextSize(0.036);

    bool first = true;
    for (int i = 0; i < (int)maps.size(); ++i) {
        if (!maps[i]) continue;
        TH1D* hProj = maps[i]->ProjectionY(
            Form("hOvl_%s_%d", tag.c_str(), i));
        hProj->SetDirectory(nullptr);

        // Normalizza al massimo nella finestra di fit
        int b1=hProj->FindBin(fit_lo), b2=hProj->FindBin(fit_hi);
        double maxVal=0;
        for (int b=b1; b<=b2; ++b)
            if (hProj->GetBinContent(b)>maxVal) maxVal=hProj->GetBinContent(b);
        if (maxVal>0) hProj->Scale(1.0/maxVal);

        hProj->SetLineColor(cols[i % nc]);
        hProj->SetLineWidth(2);
        hProj->GetXaxis()->SetRangeUser(fit_lo - 10, fit_hi + 10);
        hProj->SetTitle(Form("Multi-LET #Deltat overlay   %s"
                             ";#Deltat (ns);Normalised counts", tag.c_str()));

        if (first) { hProj->Draw("HIST"); first = false; }
        else        hProj->Draw("HIST same");

        leg->AddEntry(hProj, Form("LET = %.2f p.e.", fracs_pe[i]), "l");
    }
    leg->Draw();
    cOvl->Update(); cOvl->Modified();
    cOvl->SaveAs(Form("tot_multiLET_%s.png", tag.c_str()));
}

// ════════════════════════════════════════════════════════════
//  PROCESS A SINGLE FILE per un dato valore di LET
// ════════════════════════════════════════════════════════════
static TH2D* processOneLET(const FileInfo& info,
                             double cutoff_MHz,
                             double t_trig_start,
                             double t_trig_end,
                             double fit_lo,
                             double fit_hi,
                             double frac_pe,
                             const CalibResult& cal,
                             TWMethod tw_method,
                             bool do_pe_analysis) {

    TFile* file = TFile::Open(info.path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[SKIP] Cannot open: " << info.path << "\n";
        return nullptr;
    }
    TTree* treeCh1   = (TTree*)file->Get("ch1");
    TTree* treeLaser = (TTree*)file->Get("laser");
    if (!treeCh1 || !treeLaser) {
        std::cerr << "[SKIP] Missing trees in: " << info.path << "\n";
        file->Close(); return nullptr;
    }

    const int N = 1024;
    Double_t t1[N];
    treeCh1->SetBranchAddress("time", t1);
    treeCh1->GetEntry(0);
    double fs_MHz = 1000.0 / (t1[1] - t1[0]);

    int j_trig_start, j_trig_end;
    triggerWindowIndices(t1, N, t_trig_start, t_trig_end,
                         j_trig_start, j_trig_end, info.tag);

    double let_thr = cal.q + frac_pe * cal.m;
    std::cout << "\n[" << info.tag << "]  LET=" << frac_pe
              << " p.e.  threshold=" << let_thr << " mV\n"
              << "  Trigger window : [" << t1[j_trig_start]
              << ", " << t1[j_trig_end] << "] ns\n";

    // Raccolta eventi
    auto events = collectTOTEvents(treeCh1, treeLaser,
                                   cutoff_MHz, fs_MHz,
                                   j_trig_start, j_trig_end,
                                   let_thr, cal, info.tag);
    file->Close();

    if (events.empty()) {
        std::cout << "  [WARN] No events collected.\n";
        return nullptr;
    }

    // Tag specifico per questo LET
    std::string ltag = Form("%s_let%.2fpe", info.tag.c_str(), frac_pe);

    // Mappa 2D raw
    TH2D* h2D = fillTH2D(events, ltag, frac_pe, false);

    // Canvas raw
    drawTOTMap(h2D, ltag, false);
    drawGlobalProjection(h2D, ltag, fit_lo, fit_hi, false);
    drawSlicesAndTrend(h2D, ltag, fit_lo, fit_hi, frac_pe, false);

    // Correzione time walk
    auto events_corr = applyTimeWalkCorrection(
        events, tw_method, fit_lo, fit_hi, ltag);

    TH2D* h2D_corr = nullptr;
    if (tw_method != TWMethod::NONE) {
        h2D_corr = fillTH2D(events_corr, ltag, frac_pe, true);
        drawTOTMap(h2D_corr, ltag, true);
        drawGlobalProjection(h2D_corr, ltag, fit_lo, fit_hi, true);
        drawSlicesAndTrend(h2D_corr, ltag, fit_lo, fit_hi, frac_pe, true);
    }

    // Analisi per p.e.
    if (do_pe_analysis)
        drawByPE(tw_method != TWMethod::NONE ? events_corr : events,
                 ltag, fit_lo, fit_hi);

    // Salva ROOT
    saveTOTRoot(h2D, h2D_corr, ltag, frac_pe);

    return h2D;  // restituisce la mappa raw per l'overlay multi-LET
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════
void sipm_tot_analysis_v2() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // Scansiona file
    std::vector<FileInfo> files = findInputFiles();
    if (files.empty()) {
        std::cerr << "No data.vbias_{V}_L.root files found.\n"; return;
    }

    printFileList(files);

    int calIdx = selectOneFile(files,
        "\nSelect calibration file index (brightest): ");
    std::vector<int> anaIdx = selectMultipleFiles(files,
        "Select analysis file indices (e.g. \"1\" or \"0,1\" or \"all\"): ");

    // Parametri
    double cutoff_MHz;
    std::cout << "\nLow-pass filter cut-off [MHz]: ";
    std::cin >> cutoff_MHz;

    // Asse temporale
    double t_wf_min=0, t_wf_max=0;
    {
        TFile* f0=TFile::Open(files[calIdx].path.c_str(),"READ");
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

    // Finestra trigger
    double t_trig_start, t_trig_end;
    std::cout << "\n--- Trigger window (for TOT search) ---\n"
              << "  Waveform range  : [" << t_wf_min << ", " << t_wf_max << "] ns\n"
              << "  Baseline (fixed): [" << BASELINE_START << ", "
                                         << BASELINE_END << ") ns\n"
              << "  Start [ns]: ";
    std::cin >> t_trig_start;
    std::cout << "  End   [ns]: ";
    std::cin >> t_trig_end;
    if (t_trig_start >= t_trig_end) {
        std::cerr << "[ERROR] Invalid trigger window.\n"; return;
    }

    // Finestra fit slice
    double fit_lo, fit_hi;
    std::cout << "\n--- Fit window (signal peak in Δt) ---\n"
              << "  (Use laser_timing canvas to identify peak position)\n"
              << "  Start [ns]: ";
    std::cin >> fit_lo;
    std::cout << "  End   [ns]: ";
    std::cin >> fit_hi;
    if (fit_lo >= fit_hi) {
        std::cerr << "[WARN] Invalid fit window — using trigger window.\n";
        fit_lo = t_trig_start; fit_hi = t_trig_end;
    }

    // Lista soglie LET
    std::vector<double> fracs_pe;
    std::cout << "\n--- LET thresholds ---\n"
              << "  Enter p.e. values separated by spaces, then Enter\n"
              << "  (e.g.: 0.5 1.0 1.5)\n  > ";
    {
        std::string line;
        std::cin.ignore();
        std::getline(std::cin, line);
        std::stringstream ss(line);
        double v;
        while (ss >> v) if (v > 0) fracs_pe.push_back(v);
    }
    if (fracs_pe.empty()) {
        std::cout << "  No valid values — using 0.5 p.e.\n";
        fracs_pe.push_back(0.5);
    }

    // Correzione time walk
    TWMethod tw_method = askTimeWalkMethod();

    // Analisi per p.e.?
    char ans_pe = 0;
    while (ans_pe != 'y' && ans_pe != 'n') {
        std::cout << "\nProduce per-p.e. timing analysis? [y/n]: ";
        std::cin >> ans_pe;
    }
    bool do_pe_analysis = (ans_pe == 'y');

    // Calibrazione
    TFile* fCal = TFile::Open(files[calIdx].path.c_str(), "READ");
    if (!fCal || fCal->IsZombie()) {
        std::cerr << "[ERROR] Cannot open calibration file.\n"; return;
    }
    TTree* tCal = (TTree*)fCal->Get("ch1");
    const int NC=1024; Double_t tc[NC];
    tCal->SetBranchAddress("time",tc); tCal->GetEntry(0);
    double fs_cal = 1000.0/(tc[1]-tc[0]);

    CalibResult cal = calibrateSpectrum(tCal, cutoff_MHz, fs_cal,
                                         files[calIdx].tag);
    fCal->Close();

    if (!cal.ok) {
        std::cerr << "[ERROR] Calibration failed.\n"; return;
    }

    std::cout << "\nSettings summary:\n"
              << "  LP cut-off      : " << cutoff_MHz    << " MHz\n"
              << "  Trigger window  : [" << t_trig_start << ", "
                                         << t_trig_end   << "] ns\n"
              << "  Fit window      : [" << fit_lo       << ", "
                                         << fit_hi       << "] ns\n"
              << "  LET fractions   : ";
    for (double f : fracs_pe) std::cout << f << " ";
    std::cout << "p.e.\n"
              << "  Time walk corr  : "
              << (tw_method==TWMethod::EMPIRICAL ? "Empirical" :
                  tw_method==TWMethod::AMPLITUDE ? "Amplitude" : "None")
              << "\n\n";

    // Loop sui file di analisi
    for (int idx : anaIdx) {
        std::cout << "\n══════ " << files[idx].tag << " ══════\n";

        std::vector<TH2D*> maps;

        for (double frac : fracs_pe) {
            TH2D* h = processOneLET(files[idx], cutoff_MHz,
                                    t_trig_start, t_trig_end,
                                    fit_lo, fit_hi,
                                    frac, cal, tw_method, do_pe_analysis);
            maps.push_back(h);
        }

        // Overlay multi-LET (solo se più di una soglia)
        if (fracs_pe.size() > 1)
            drawMultiLETOverlay(fracs_pe, maps, files[idx].tag,
                                fit_lo, fit_hi);
    }

    std::cout << "\n=== All files processed. ===\n";
}
