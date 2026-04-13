#pragma once
// VbiasAnalysis.h
// Full analysis pipeline for a single Vbias value.
//
// collectTOTEvents_v4() — event loop over a TChain with progress bar
// peakFinderFromCache() — quick Δt peak estimate from event cache
// processOneVbias()     — orchestrates everything: file scan, calibration,
//                         event loop, plots, sideband fit, time-walk correction
//
// RAM strategy:
//   - The event loop writes TOTEvent directly to the cache file one event at
//     a time via a streaming TTree (no in-memory accumulation).
//   - All canvases are closed after each LET iteration (OutCtx::closeAllCanvases).
//   - WaveformPlotter and PersistencePlotter are disabled by default.
//   - The TChain is deleted after each Vbias is done.

#include "Config.h"
#include "OutputManager.h"
#include "SignalProcessing.h"
#include "Calibration.h"
#include "CalibIO.h"
#include "EventCache.h"
#include "TOTAnalysis.h"
#include "TimingCorrection.h"
#include "TOTPlotting.h"
#include "VbiasSummary.h"
#include "FilterDiagnostics.h"
#include "ProgressBar.h"
#include "SidebandAnalysis.h"
#include "WaveformPlotter.h"
#include "PersistencePlot.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TSystem.h>
#include <TROOT.h>

// ════════════════════════════════════════════════════════════
//  collectTOTEvents_v4
//  Event loop su TChain. Scrive direttamente su file ROOT (streaming),
//  senza tenere in RAM il vettore completo degli eventi.
//
//  RAM: O(1) per evento — solo i branch buffer del TTree attivo.
//  Il file cache viene scritto evento per evento: nessun accumulo in memoria.
// ════════════════════════════════════════════════════════════
static void collectTOTEvents_v4(
        TTree*              treeCh1,
        TTree*              treeLaser,
        double              cutoff_MHz,
        double              fs_MHz,
        int                 j_trig_start,
        int                 j_trig_end,
        double              let_thr,
        const CalibResult&  cal,
        const std::string&  tag,
        const std::string&  dataDir,
        int                 vbias,
        double              frac_pe,
        bool                use_filter,
        WaveformPlotter*    wfPlotter    = nullptr,
        WaveformPlotter*    wfPlotterNeg = nullptr,
        PersistencePlotter* persPlotter  = nullptr)
{
    const int N = 1024;
    Double_t t1[N], a1[N], tL[N], aL[N];
    treeCh1  ->SetBranchAddress("time",      t1);
    treeCh1  ->SetBranchAddress("amplitude", a1);
    treeLaser->SetBranchAddress("time",      tL);
    treeLaser->SetBranchAddress("amplitude", aL);

    Long64_t nEntries = treeCh1->GetEntries();
    long nNoLaser = 0, nNoTOT = 0, nAccepted = 0;

    // ── Apri il file di output e crea il TTree in streaming ──────
    // Scriviamo direttamente sul file cache finale, evento per evento.
    // Nessun vettore in RAM: RAM = O(1) per evento.
    std::string cachePath = eventCachePath(vbias, frac_pe, cutoff_MHz,
                                            cal.laser_thr, dataDir, use_filter);
    TFile* fOut = new TFile(cachePath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "[collectTOTEvents_v4] Cannot create cache: " << cachePath << "\n";
        delete fOut;
        return;
    }

    Double_t tot_br, delta_t_br, amp_max_br, laser_thr_br = cal.laser_thr;
    Int_t    n_pe_br;
    TTree* tOut = new TTree("events",
        Form("TOTEvents vbias%d let%.2fpe cut%.0fMHz lthr%.4fmV",
             vbias, frac_pe, cutoff_MHz, cal.laser_thr));
    tOut->Branch("tot",             &tot_br,       "tot/D");
    tOut->Branch("delta_t",         &delta_t_br,   "delta_t/D");
    tOut->Branch("amp_max",         &amp_max_br,   "amp_max/D");
    tOut->Branch("n_pe",            &n_pe_br,      "n_pe/I");
    tOut->Branch("laser_thr_saved", &laser_thr_br, "laser_thr_saved/D");
    // Imposta AutoSave ogni 50k eventi per evitare di perdere dati su crash
    tOut->SetAutoSave(50000);

    ProgressBar bar(nEntries,
        "TOT collect  tag=" + tag +
        "  thr=" + std::to_string((int)std::round(let_thr)) + " mV" +
        (use_filter ? "" : "  [NO FILTER]"));

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (bar.update(i, nNoLaser, nNoTOT, nAccepted)) break;

        treeCh1  ->GetEntry(i);
        treeLaser->GetEntry(i);

        double t_laser = laserTriggerTime(tL, aL, N, cal.laser_thr);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        std::vector<double> v_t(t1, t1 + N), v_a(a1, a1 + N);
        std::vector<double> af = correctBaseline(v_t, v_a,
                                                 BASELINE_START, BASELINE_END);
        if (use_filter)
            af = butterworthLowPass(af, cutoff_MHz, fs_MHz);

        // amp_max nel solo trigger window — evita bias da afterpulse
        double amp_max = *std::max_element(af.begin() + j_trig_start,
                                            af.begin() + j_trig_end + 1);

        auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                                            j_trig_start, j_trig_end,
                                            100, 1000.0, 0.5, 10, 0.3, 3, 50, false);
        if (t_rise < 0 || t_fall < 0) { ++nNoTOT; continue; }

        double tot     = t_fall - t_rise;
        double delta_t = t_rise - t_laser;
        if (tot <= 0.0 || tot >= 150.0) continue;

        int n_pe = estimatePE(amp_max, cal);

        // Scrivi direttamente su disco — nessuna copia in RAM
        tot_br     = tot;
        delta_t_br = delta_t;
        amp_max_br = amp_max;
        n_pe_br    = n_pe;
        tOut->Fill();
        ++nAccepted;

        if (wfPlotter    && wfPlotter->collecting())
            wfPlotter->collect(i, v_t, af, t_laser, t_rise, t_fall, let_thr, n_pe);
        if (wfPlotterNeg && wfPlotterNeg->collecting())
            wfPlotterNeg->collect(i, v_t, af, t_laser, t_rise, t_fall, let_thr, n_pe);
        if (persPlotter)
            persPlotter->collect(v_t, af, t_laser, t_rise, let_thr);
    }
    bar.done();

    fOut->Write();
    fOut->Close();
    delete fOut;

    std::cout << "  [EventCache] Saved " << nAccepted
              << " events -> " << cachePath << "\n";
}

// ════════════════════════════════════════════════════════════
//  peakFinderFromCache
//  Carica la cache (se esiste) e restituisce la posizione del picco Δt
//  per eventi con TOT piccolo (15° percentile).
//  Ritorna {peak_ns, peak-5, peak+5} come finestra di fit suggerita.
//  Ritorna {-1,-1,-1} se la cache non esiste ancora.
// ════════════════════════════════════════════════════════════
static std::tuple<double,double,double> peakFinderFromCache(
        int vbias, double frac_pe, double cutoff_MHz, double laser_thr,
        const std::string& dataDir = DATA_DIR)
{
    std::vector<TOTEvent> events;
    if (!loadEvents(events, vbias, frac_pe, cutoff_MHz, laser_thr, dataDir))
        return {-1.0, -1.0, -1.0};
    if (events.empty()) return {-1.0, -1.0, -1.0};

    std::vector<double> tots;
    for (auto& e : events) tots.push_back(e.tot);
    std::sort(tots.begin(), tots.end());
    double tot_cut = tots[std::min((size_t)(tots.size()*0.15), tots.size()-1)];
    tot_cut = std::max(2.0, std::min(tot_cut, 15.0));

    int nB = std::max(200, (int)std::round(150.0 / 0.1));
    TH1D* hPk = new TH1D(Form("hPFCache_v%d_%.2fpe", vbias, frac_pe),
                          "", nB, -50.0, 200.0);
    hPk->SetDirectory(nullptr);
    for (auto& e : events)
        if (e.tot > 0 && e.tot < tot_cut) hPk->Fill(e.delta_t);

    double tpk = hPk->GetBinCenter(hPk->GetMaximumBin());
    delete hPk;
    return {tpk, tpk - 5.0, tpk + 5.0};
}

// ════════════════════════════════════════════════════════════
//  processOneVbias
//  Pipeline completa per un Vbias:
//    1. Cerca file data.vbias_<V>_run_<N>.root in DATA_DIR
//    2. Legge sampling rate dal primo file
//    3. Usa calibrazione passata come argomento (già scelta dal main)
//    4. Costruisce due TChain paralleli (ch1 + laser), ordinati per run
//    5. Per ogni LET threshold:
//       a. Carica o raccoglie eventi (streaming su disco)
//       b. Peak-finder diagnostico
//       c. TOT map, proiezione Δt globale, slices TOT
//       d. Correzione time walk (se richiesta)
//       e. Sideband subtraction + S+B fit
//    6. Multi-LET overlay
//  Ritorna mappa: frac_pe → {sigma [ns], sigmaErr [ns]}
// ════════════════════════════════════════════════════════════
static std::map<double, std::pair<double,double>> processOneVbias(
        int                        vbias,
        const std::vector<double>& fracs_pe,
        double                     cutoff_MHz,
        double                     t_trig_start,
        double                     t_trig_end,
        double                     fit_lo,
        double                     fit_hi,
        TWMethod                   tw_method,
        bool                       do_pe_analysis,
        bool                       use_filter,
        OutCtx&                    ctx,
        const CalibResult&         cal_in)   // <-- calibrazione già caricata dal main
{
    std::map<double, std::pair<double,double>> sigmaOut;

    // 1. Trova i file di run
    const std::string dataDir = DATA_DIR;
    const std::string pattern = "data.vbias_" + std::to_string(vbias) + "_run_";
    std::map<int, std::string> foundRuns;
    {
        void* dirp = gSystem->OpenDirectory(dataDir.c_str());
        if (!dirp) {
            std::cerr << "[ERROR] Cannot open: " << dataDir << "\n";
            return sigmaOut;
        }
        const char* entry;
        while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
            std::string fname(entry);
            if (fname.find(pattern) == std::string::npos) continue;
            if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
            size_t pos = fname.find("_run_");
            if (pos == std::string::npos) continue;
            try {
                std::string sub = fname.substr(pos + 5);
                size_t dot = sub.find(".root");
                if (dot != std::string::npos) sub = sub.substr(0, dot);
                foundRuns[std::stoi(sub)] = dataDir + "/" + fname;
            } catch (...) {}
        }
        gSystem->FreeDirectory(dirp);
    }
    if (foundRuns.empty()) {
        std::cerr << "[WARN] No files for Vbias=" << vbias << "\n";
        return sigmaOut;
    }
    std::cout << "\n+-- Vbias=" << vbias << " V  " << foundRuns.size() << " run(s)\n";
    for (auto& [r, path] : foundRuns)
        std::cout << "    run " << std::setw(2) << r << " : " << path << "\n";

    // 2. Sampling rate
    double fs_MHz = 0.0;
    {
        TFile* f0 = TFile::Open(foundRuns.begin()->second.c_str(), "READ");
        if (!f0 || f0->IsZombie()) {
            std::cerr << "[ERROR] Cannot open first run file.\n"; return sigmaOut;
        }
        TTree* tr = (TTree*)f0->Get("ch1");
        if (tr) {
            const int N0 = 1024; Double_t tb[N0];
            tr->SetBranchAddress("time", tb); tr->GetEntry(0);
            fs_MHz = 1000.0 / (tb[1] - tb[0]);
        }
        f0->Close();
    }
    if (fs_MHz <= 0.0) {
        std::cerr << "[ERROR] Cannot read sampling rate.\n"; return sigmaOut;
    }
    std::cout << "  Sampling rate: " << fs_MHz << " MHz\n";

    // 3. Usa la calibrazione passata dal main (già scelta interattivamente).
    //    Esegui solo FilterDiagnostics per aggiornare laser_thr in memoria
    //    (non viene risalvata su disco — il main ha già caricato il valore corretto).
    CalibResult cal = cal_in;  // copia locale per aggiornare laser_thr
    const std::string calTag = "vbias" + std::to_string(vbias);

    drawFilterDiagnostics(foundRuns.begin()->second, cutoff_MHz, fs_MHz,
                          t_trig_start, t_trig_end, cal, calTag, ctx);

    if (!cal.ok) {
        std::cerr << "[ERROR] Calibration not valid for Vbias=" << vbias << "\n";
        return sigmaOut;
    }
    std::cout << "  Gain=" << cal.m << " mV/p.e.  Offset=" << cal.q << " mV\n";

    // 4. TChain — stesso ordine per ch1 e laser garantisce allineamento 1:1
    TChain* chainCh1   = new TChain("ch1");
    TChain* chainLaser = new TChain("laser");
    for (auto& [r, path] : foundRuns) {
        chainCh1  ->Add(path.c_str());
        chainLaser->Add(path.c_str());
    }
    Long64_t totEvt = chainCh1->GetEntries();
    std::cout << "  Total events: " << totEvt << "\n";
    if (totEvt != chainLaser->GetEntries()) {
        std::cerr << "[ERROR] ch1 and laser chain sizes differ.\n";
        delete chainCh1; delete chainLaser; return sigmaOut;
    }

    // 5. Indici del trigger window (calcolati una volta dal primo file)
    int j_start = 0, j_end = 1023;
    {
        TFile* f0 = TFile::Open(foundRuns.begin()->second.c_str(), "READ");
        if (f0 && !f0->IsZombie()) {
            TTree* tr = (TTree*)f0->Get("ch1");
            if (tr) {
                const int N0 = 1024; Double_t tb[N0];
                tr->SetBranchAddress("time", tb); tr->GetEntry(0);
                triggerWindowIndices(tb, N0, t_trig_start, t_trig_end,
                                     j_start, j_end, calTag);
            }
            f0->Close();
        }
    }
    std::cout << "  Trigger indices: [" << j_start << ", " << j_end << "]\n";

    // 6. Loop sui threshold LET
    const double sb_width = SB::DEFAULT_SB_WIDTH;
    const double sb_ext   = SB::DEFAULT_EXT;
    std::vector<TH2D*> maps_for_overlay;

    // Diagnostica waveform e persistence disabilitate di default per risparmiare RAM.
    // Per abilitarle: metti true e ricompila.
    const bool doWaveformPlots = false;
    const bool doPersistence   = false;

    for (double frac : fracs_pe) {
        if (gROOT->IsInterrupted()) break;

        double let_thr = cal.q + frac * cal.m;
        std::cout << "\n  +-- LET=" << frac << " p.e.  thr="
                  << std::fixed << std::setprecision(2) << let_thr << " mV\n";
        const std::string ltag = Form("vbias%d_let%.2fpe", vbias, frac);

        // ── Prova la cache. Se assente, esegui il loop in streaming ──
        std::vector<TOTEvent> events;
        bool fromCache = loadEvents(events, vbias, frac, cutoff_MHz,
                                    cal.laser_thr, dataDir, use_filter);

        WaveformPlotter    wfPlotter   (10.0, 100);
        WaveformPlotter    wfPlotterNeg(10.0, 100, -30.0, 0.0);
        PersistencePlotter persPlotter (let_thr);

        if (!fromCache) {
            // Streaming: scrive direttamente su disco senza accumulare in RAM
            collectTOTEvents_v4(chainCh1, chainLaser,
                cutoff_MHz, fs_MHz, j_start, j_end, let_thr, cal, ltag,
                dataDir, vbias, frac, use_filter,
                doWaveformPlots ? &wfPlotter    : nullptr,
                doWaveformPlots ? &wfPlotterNeg : nullptr,
                doPersistence   ? &persPlotter  : nullptr);
            // Carica in RAM solo dopo che il file è chiuso e completo
            loadEvents(events, vbias, frac, cutoff_MHz,
                       cal.laser_thr, dataDir, use_filter);
        } else {
            // Cache hit: opzionale rescan per persistence (prima N eventi)
            if (doPersistence) {
                std::cout << "  [Cache hit] Re-scanning for persistence map...\n";
                const Long64_t PERSIST_RESCAN_MAX = 50000;
                const int N = 1024;
                Double_t t1[N], a1[N], tL[N], aL[N];
                chainCh1  ->SetBranchAddress("time",      t1);
                chainCh1  ->SetBranchAddress("amplitude", a1);
                chainLaser->SetBranchAddress("time",      tL);
                chainLaser->SetBranchAddress("amplitude", aL);
                Long64_t nEntries = std::min(chainCh1->GetEntries(), PERSIST_RESCAN_MAX);
                for (Long64_t i = 0; i < nEntries; ++i) {
                    if (i % 5000 == 0) {
                        if (gROOT->IsInterrupted()) break;
                        gSystem->ProcessEvents();
                    }
                    chainCh1->GetEntry(i); chainLaser->GetEntry(i);
                    double t_laser = laserTriggerTime(tL, aL, N, cal.laser_thr);
                    if (t_laser < -900.0) continue;
                    std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
                    std::vector<double> af = correctBaseline(v_t, v_a,
                                                             BASELINE_START, BASELINE_END);
                    if (use_filter) af = butterworthLowPass(af, cutoff_MHz, fs_MHz);
                    // Usa gli stessi parametri del loop principale
                    auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                                                        j_start, j_end,
                                                        100, 1000.0, 0.5, 10, 0.3, 3, 50, false);
                    if (t_rise < 0) continue;
                    persPlotter.collect(v_t, af, t_laser, t_rise, let_thr);
                }
                std::cout << "  [Persistence] " << persPlotter.nWaveforms()
                          << " waveforms (from first " << nEntries << " events).\n";
            }
        }

        if (doWaveformPlots && wfPlotter.size()    > 0) wfPlotter.draw(ltag, ctx);
        if (doWaveformPlots && wfPlotterNeg.size() > 0) wfPlotterNeg.draw(ltag + "_negDt", ctx);
        if (doPersistence && persPlotter.nWaveforms() > 0) persPlotter.draw(ltag, ctx);

        if (events.empty()) {
            std::cout << "  +-- No events, skipping.\n";
            maps_for_overlay.push_back(nullptr);
            continue;
        }

        // ── Peak-finder diagnostico ──────────────────────────────
        {
            std::vector<double> tots;
            for (auto& e : events) tots.push_back(e.tot);
            std::sort(tots.begin(), tots.end());
            double tot_cut = tots[std::min((size_t)(tots.size()*0.15), tots.size()-1)];
            tot_cut = std::max(2.0, std::min(tot_cut, 15.0));

            int   nB  = std::max(200, (int)std::round(150.0 / 0.1));
            TH1D* hPk = new TH1D(Form("hPeakFinder_%s", ltag.c_str()),
                Form("Peak finder (TOT < %.1f ns)  %s;#Deltat (ns);Events",
                     tot_cut, ltag.c_str()), nB, -50.0, 200.0);
            hPk->SetDirectory(nullptr);
            int nSmall = 0;
            for (auto& e : events)
                if (e.tot > 0 && e.tot < tot_cut) { hPk->Fill(e.delta_t); ++nSmall; }

            double tpk = hPk->GetBinCenter(hPk->GetMaximumBin());
            double hpk = hPk->GetMaximum();
            std::cout << "  [PeakFinder] TOT<" << tot_cut << " ns: " << nSmall
                      << " events  Peak at " << tpk << " ns\n"
                      << "  [PeakFinder] Suggested window: ["
                      << tpk-5.0 << ", " << tpk+5.0 << "] ns\n";

            TCanvas* cPk = new TCanvas(Form("cPeakFinder_%s", ltag.c_str()),
                Form("Peak finder -- %s", ltag.c_str()), 900, 500);
            cPk->SetGrid();
            cPk->SetLeftMargin(PAD_LEFT); cPk->SetRightMargin(PAD_RIGHT);
            cPk->SetBottomMargin(PAD_BOTTOM); cPk->SetTopMargin(PAD_TOP);
            hPk->SetLineColor(kAzure+1); hPk->SetLineWidth(2);
            hPk->GetXaxis()->SetRangeUser(tpk - 30.0, tpk + 30.0);
            hPk->Draw("HIST");
            TLine* lPk = new TLine(tpk, 0, tpk, hpk);
            lPk->SetLineColor(kRed+1); lPk->SetLineStyle(2); lPk->SetLineWidth(2);
            lPk->Draw("same");
            {
                TPaveText* pt = new TPaveText(0.55, 0.72, 0.93, 0.88, "NDC");
                pt->SetBorderSize(1); pt->SetFillColor(0);
                pt->SetTextFont(42); pt->SetTextSize(0.036);
                pt->AddText(Form("TOT < %.1f ns  (%d events)", tot_cut, nSmall));
                pt->AddText(Form("Peak at #Deltat = %.2f ns", tpk));
                pt->AddText(Form("Suggested: [%.1f, %.1f] ns", tpk-5.0, tpk+5.0));
                pt->Draw();
                delete pt;
            }
            cPk->Update(); cPk->Modified();
            ctx.savePNG(cPk, Form("peak_finder_%s.png", ltag.c_str()));
            delete hPk;
            delete lPk;
            delete cPk;
        }

        // ── Plot standard ────────────────────────────────────────
        TH2D* h2D = fillTH2D(events, ltag, frac, false);
        drawTOTMap(h2D, ltag, ctx, false);
        drawGlobalProjection(h2D, ltag, fit_lo, fit_hi, ctx, false, &events);
        drawDeltaTHistogram(events, ltag, fit_lo, fit_hi, ctx, false);

        // ── Correzione time walk ─────────────────────────────────
        auto events_corr = applyTimeWalkCorrection(
            events, tw_method, fit_lo, fit_hi, ltag, ctx);

        TH2D* h2D_final = h2D;
        if (tw_method != TWMethod::NONE) {
            TH2D* h2D_corr = fillTH2D(events_corr, ltag, frac, true);
            drawTOTMap(h2D_corr, ltag, ctx, true);
            drawGlobalProjection(h2D_corr, ltag, fit_lo, fit_hi, ctx, true, &events_corr);
            drawDeltaTHistogram(events_corr, ltag, fit_lo, fit_hi, ctx, true);
            if (do_pe_analysis) drawByPE(events_corr, ltag, fit_lo, fit_hi, ctx);
            ctx.saveRoot(h2D, h2D_corr, ltag, frac);
            h2D_final = h2D_corr;
        } else {
            if (do_pe_analysis) drawByPE(events, ltag, fit_lo, fit_hi, ctx);
            ctx.saveRoot(h2D, nullptr, ltag, frac);
        }

        // ── Sideband fit ─────────────────────────────────────────
        {
            const auto& evSrc = (tw_method != TWMethod::NONE) ? events_corr : events;

            const double tot_cut_sb = 0.0;
            std::vector<TOTEvent> evSB;
            evSB.reserve(evSrc.size());
            for (auto& e : evSrc)
                if (tot_cut_sb <= 0.0 || e.tot < tot_cut_sb)
                    evSB.push_back(e);

            double sb_fit_lo = fit_lo, sb_fit_hi = fit_hi;
            {
                double win    = fit_hi - fit_lo;
                double centre = 0.5 * (fit_lo + fit_hi);
                int    nB     = std::max(50, (int)std::round(win / 0.2));
                TH1D*  hPk    = new TH1D(Form("hSBpk_%s", ltag.c_str()),
                                         "", nB, fit_lo, fit_hi);
                hPk->SetDirectory(nullptr);
                for (auto& e : evSB)
                    if (e.delta_t >= fit_lo && e.delta_t < fit_hi)
                        hPk->Fill(e.delta_t);

                int bMax = hPk->FindBin(centre - 2.0);
                for (int b = bMax; b <= hPk->GetNbinsX(); ++b)
                    if (hPk->GetBinContent(b) > hPk->GetBinContent(bMax)) bMax = b;
                double rough_pk = hPk->GetBinCenter(bMax);
                double sw = 0, swx = 0;
                for (int b = 1; b <= hPk->GetNbinsX(); ++b) {
                    double x = hPk->GetBinCenter(b), c = hPk->GetBinContent(b);
                    if (std::abs(x - rough_pk) <= 4.0) { sw += c; swx += c*x; }
                }
                double tpk = (sw > 0) ? swx/sw : rough_pk;
                delete hPk;
                if (std::abs(tpk - centre) > 0.3) {
                    sb_fit_lo = tpk - win/2.0;
                    sb_fit_hi = tpk + win/2.0;
                    std::cout << "  [SB] Window re-centred: ["
                              << std::fixed << std::setprecision(1)
                              << sb_fit_lo << ", " << sb_fit_hi << "] ns\n";
                }
            }

            SBResult sbr = fitSplusB(evSB, sb_fit_lo, sb_fit_hi, sb_width, sb_ext, ltag);
            drawSidebandResult(sbr, ltag, sb_fit_lo, sb_fit_hi, sb_width, ctx);

            double sigma = 0.0, sigmaErr = 0.0;
            if (sbr.fit_ok && sbr.sigma > 0.1 && sbr.sigma < (sb_fit_hi - sb_fit_lo)) {
                sigma = sbr.sigma; sigmaErr = sbr.sigmaErr;
                sbr.sigma_source = "S+B fit";
            } else if (sbr.sub_ok && sbr.sigma_sub > 0.1
                       && sbr.sigma_sub < (sb_fit_hi - sb_fit_lo)) {
                sigma = sbr.sigma_sub; sigmaErr = sbr.sigmaErr_sub;
                sbr.sigma_source = "subtracted gaussian";
            } else {
                int nBP = std::max(400, (int)std::round((sb_fit_hi - sb_fit_lo) / 0.02));
                TH1D* hF = new TH1D(Form("hSumFit_%s", ltag.c_str()), "",
                                     nBP, sb_fit_lo - 1.0, sb_fit_hi + 1.0);
                hF->SetDirectory(nullptr);
                for (auto& e : evSB)
                    if (e.delta_t >= sb_fit_lo && e.delta_t <= sb_fit_hi)
                        hF->Fill(e.delta_t);
                double ctr = hF->GetBinCenter(hF->GetMaximumBin());
                double sg0 = std::max(hF->GetRMS(), 0.1);
                TF1* fG = new TF1(Form("fSumG_%s", ltag.c_str()), "gaus",
                    std::max(sb_fit_lo, ctr - 2.0*sg0),
                    std::min(sb_fit_hi, ctr + 2.0*sg0));
                fG->SetParameters(hF->GetMaximum(), ctr, sg0);
                hF->Fit(fG, "RQN");
                sigma = std::abs(fG->GetParameter(2));
                sigmaErr = fG->GetParError(2);
                sbr.sigma_source = "legacy gaussian";
                std::cout << "  [SB] WARNING: S+B failed, using plain Gaussian.\n";
                delete hF; delete fG;
            }

            std::cout << "  [SB] sigma = " << std::fixed << std::setprecision(3)
                      << sigma << " +/- " << sigmaErr << " ns"
                      << "  (" << sbr.sigma_source << ")\n";

            if (sigma > 0.0 && sigma < (sb_fit_hi - sb_fit_lo))
                sigmaOut[frac] = {sigma, sigmaErr};
            cleanupSBResult(sbr);
        }

        maps_for_overlay.push_back(h2D_final);
        std::cout << "  +-- LET " << frac << " p.e. done.\n";

        // Libera tutte le canvas accumulate durante questa iterazione LET.
        // ROOT mantiene ogni TCanvas in gROOT->GetListOfCanvases() finché
        // non viene chiusa — con ~20 canvas per LET × N Vbias questa è
        // la principale fonte di crescita RAM durante run lunghi.
        OutCtx::closeAllCanvases();
    }

    // 7. Multi-LET overlay
    if (fracs_pe.size() > 1 && !maps_for_overlay.empty())
        drawMultiLETOverlay(fracs_pe, maps_for_overlay, calTag, fit_lo, fit_hi, ctx);

    delete chainCh1;
    delete chainLaser;
    return sigmaOut;
}
