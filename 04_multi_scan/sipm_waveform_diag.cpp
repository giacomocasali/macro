/**
 * sipm_waveform_diag.cpp
 * ======================
 * Visual diagnostics: shows waveforms with Δt < 0 and Δt > 0,
 * plus a persistence plot of all accepted waveforms.
 *
 * Processes file-by-file (low RAM).
 * Produces 3 canvases:
 *   1. 10 waveforms with Δt < 0 (why negative?)
 *   2. 10 waveforms with Δt > 0 (reference)
 *   3. Persistence plot (all accepted waveforms, aligned on t_laser)
 *
 * v3 changes:
 *   - Fixed post-fall check: uses threshold*0.5 as post-fall limit
 *     instead of edge_limit (which was disabled at 1130 mV when
 *     edge_thr_frac=100, letting long-tail events through).
 *   - Pre-laser quiet check to reject dark counts before laser.
 *   - All output in English.
 *   - Interactive calibration selection (supports cut0).
 *   - Canvas remain open and interactive (no gApplication->Run block).
 *
 * Compile: .L sipm_waveform_diag.cpp+
 * Run:     sipm_waveform_diag()
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/TOTAnalysis.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>

// ── Saved waveform structure ────────────────────────────────
struct SavedWF {
    std::vector<double> time, amp;    // processed waveform (filtered or raw)
    std::vector<double> amp_raw;      // raw waveform (for comparison overlay)
    double t_laser, t_rise, t_fall;
    double delta_t, tot, amp_max;
    double let_thr;
    int    n_pe;
    Long64_t entry;
    int    run;
};

// preLaserQuiet() is now in TOTAnalysis.h



// ── Draw a grid of waveforms ────────────────────────────────
static void drawWaveformGrid(const std::vector<SavedWF>& wfs,
                              const std::string& title,
                              const std::string& filename,
                              OutCtx& ctx)
{
    if (wfs.empty()) { std::cout << "  [WF] No waveforms for " << title << "\n"; return; }

    int n = std::min((int)wfs.size(), 10);
    int COLS = 5, ROWS = 2;

    TCanvas* c = new TCanvas(Form("cWF_%s", filename.c_str()),
        title.c_str(), 350*COLS, 280*ROWS);
    c->SetFillColor(0);
    c->Divide(COLS, ROWS, 0.002, 0.002);

    for (int k = 0; k < n; ++k) {
        const SavedWF& w = wfs[k];
        c->cd(k+1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.03);
        gPad->SetBottomMargin(0.18);
        gPad->SetTopMargin(0.15);
        gPad->SetGrid();
        gPad->SetFillColor(0);

        double x_lo = w.t_laser - 30.0;
        double x_hi = w.t_laser + 80.0;
        if (!w.time.empty()) {
            x_lo = std::max(x_lo, w.time.front());
            x_hi = std::min(x_hi, w.time.back());
        }

        double amp_max_zoom = 0;
        for (size_t j = 0; j < w.time.size(); ++j)
            if (w.time[j] >= x_lo && w.time[j] <= x_hi)
                if (w.amp[j] > amp_max_zoom) amp_max_zoom = w.amp[j];
        amp_max_zoom = std::max(amp_max_zoom, w.let_thr * 2.0);
        double y_lo = -10.0, y_hi = amp_max_zoom * 1.4;

        std::vector<double> gx, gy, grx, gry;
        for (size_t j = 0; j < w.time.size(); ++j) {
            if (w.time[j] >= x_lo && w.time[j] <= x_hi) {
                gx.push_back(w.time[j]); gy.push_back(w.amp[j]);
                if (j < w.amp_raw.size()) { grx.push_back(w.time[j]); gry.push_back(w.amp_raw[j]); }
            }
        }

        // Raw (grey)
        if (!grx.empty()) {
            TGraph* grRaw = new TGraph(grx.size(), grx.data(), gry.data());
            grRaw->SetTitle("");
            grRaw->SetLineColor(kGray);
            grRaw->SetLineWidth(1);
            grRaw->GetXaxis()->SetTitle("t (ns)");
            grRaw->GetYaxis()->SetTitle("A (mV)");
            grRaw->GetXaxis()->SetTitleSize(0.07);
            grRaw->GetYaxis()->SetTitleSize(0.07);
            grRaw->GetXaxis()->SetLabelSize(0.06);
            grRaw->GetYaxis()->SetLabelSize(0.06);
            grRaw->GetXaxis()->SetTitleOffset(0.9);
            grRaw->GetYaxis()->SetTitleOffset(0.85);
            grRaw->GetXaxis()->SetLimits(x_lo, x_hi);
            grRaw->GetYaxis()->SetRangeUser(y_lo, y_hi);
            grRaw->Draw("AL");
        }

        // Filtered/processed (blue)
        TGraph* gr = new TGraph(gx.size(), gx.data(), gy.data());
        gr->SetTitle("");
        gr->SetLineColor(kAzure+1);
        gr->SetLineWidth(2);
        if (grx.empty()) {
            gr->GetXaxis()->SetTitle("t (ns)");
            gr->GetYaxis()->SetTitle("A (mV)");
            gr->GetXaxis()->SetTitleSize(0.07);
            gr->GetYaxis()->SetTitleSize(0.07);
            gr->GetXaxis()->SetLabelSize(0.06);
            gr->GetYaxis()->SetLabelSize(0.06);
            gr->GetXaxis()->SetLimits(x_lo, x_hi);
            gr->GetYaxis()->SetRangeUser(y_lo, y_hi);
            gr->Draw("AL");
        } else {
            gr->Draw("L same");
        }

        // LET threshold (blue dashed)
        TLine* lThr = new TLine(x_lo, w.let_thr, x_hi, w.let_thr);
        lThr->SetLineColor(kBlue); lThr->SetLineStyle(2); lThr->SetLineWidth(1);
        lThr->Draw("same");

        // t_laser (green)
        if (w.t_laser > x_lo && w.t_laser < x_hi) {
            TLine* lLas = new TLine(w.t_laser, y_lo, w.t_laser, y_hi*0.85);
            lLas->SetLineColor(kGreen+2); lLas->SetLineStyle(2); lLas->SetLineWidth(2);
            lLas->Draw("same");
        }

        // t_rise (red)
        if (w.t_rise > x_lo && w.t_rise < x_hi) {
            TLine* lR = new TLine(w.t_rise, y_lo, w.t_rise, y_hi*0.85);
            lR->SetLineColor(kRed+1); lR->SetLineStyle(2); lR->SetLineWidth(2);
            lR->Draw("same");
        }

        // t_fall (orange)
        if (w.t_fall > x_lo && w.t_fall < x_hi) {
            TLine* lF = new TLine(w.t_fall, y_lo, w.t_fall, y_hi*0.7);
            lF->SetLineColor(kOrange+7); lF->SetLineStyle(2); lF->SetLineWidth(1);
            lF->Draw("same");
        }

        // Info box
        TPaveText* pt = new TPaveText(0.01, 0.86, 0.99, 0.99, "NDC");
        pt->SetBorderSize(0); pt->SetFillColor(0); pt->SetFillStyle(0);
        pt->SetTextFont(42); pt->SetTextSize(0.065); pt->SetTextAlign(22);
        pt->AddText(Form("run%d ev%lld  #Deltat=%.1f ns  TOT=%.1f ns  %dpe",
                         w.run, w.entry, w.delta_t, w.tot, w.n_pe));
        pt->Draw();
    }

    c->cd();
    c->Update(); c->Modified();
    ctx.savePNG(c, Form("%s.png", filename.c_str()));
    std::cout << "  [WF] " << n << " waveforms saved: " << filename << ".png\n";
}


void sipm_waveform_diag()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();
    const std::string dataDir = DATA_DIR;
    const int N = 1024;

    std::cout << "\n+==========================================================+\n"
              << "|    WAVEFORM DIAGNOSTICS  v3                               |\n"
              << "+==========================================================+\n"
              << "  Data: " << dataDir << "\n\n";

    // ── Parameters ────────────────────────────────────────────
    int vbias;
    std::cout << "Vbias [V]: " << std::flush; std::cin >> vbias;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // ── Find all available calibrations (including cut0) ─────
    std::vector<double> cutoffs;
    {
        std::string prefix = "calib_vbias" + std::to_string(vbias) + "_cut";
        std::string suffix = "mhz.root";
        void* dirp = gSystem->OpenDirectory(dataDir.c_str());
        if (dirp) {
            const char* entry;
            while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
                std::string fname(entry);
                if (fname.size() < prefix.size()+suffix.size()) continue;
                if (fname.substr(0,prefix.size()) != prefix) continue;
                if (fname.substr(fname.size()-suffix.size()) != suffix) continue;
                std::string mid = fname.substr(prefix.size(), fname.size()-prefix.size()-suffix.size());
                try {
                    double co = std::stod(mid);
                    if (co >= 0) cutoffs.push_back(co);
                } catch(...) {}
            }
            gSystem->FreeDirectory(dirp);
        }
    }
    std::sort(cutoffs.begin(), cutoffs.end());
    if (cutoffs.empty()) {
        std::cerr << "No calibration found for Vbias=" << vbias << ".\n";
        return;
    }

    // ── Interactive calibration selection ─────────────────────
    double cutoff_MHz = cutoffs.back();
    if (cutoffs.size() > 1) {
        std::cout << "  Available calibrations: ";
        for (double co : cutoffs) std::cout << (int)co << " ";
        std::cout << "MHz\n";
        std::cout << "  Which cutoff? [" << (int)cutoff_MHz << "]: " << std::flush;
        std::string line;
        std::getline(std::cin, line);
        if (!line.empty()) {
            try {
                double inp = std::stod(line);
                double best = cutoff_MHz, bestD = 1e9;
                for (double co : cutoffs)
                    if (std::abs(co - inp) < bestD) { bestD = std::abs(co - inp); best = co; }
                cutoff_MHz = best;
            } catch(...) {}
        }
    }
    std::cout << "  Using calibration: cut" << (int)cutoff_MHz << "mhz\n";

    CalibResult cal;
    if (!loadCalibration(cal, vbias, cutoff_MHz, dataDir)) {
        std::cerr << "Failed to load calibration.\n";
        return;
    }
    std::cout << "  Gain=" << cal.m << " mV/p.e.  Offset=" << cal.q << "\n";

    double frac_pe;
    std::cout << "LET threshold [p.e.]: " << std::flush; std::cin >> frac_pe;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    double let_thr = cal.q + frac_pe * cal.m;
    std::cout << "  LET threshold = " << std::fixed << std::setprecision(2) << let_thr << " mV\n";

    char filt_ans = 0;
    while (filt_ans != 'y' && filt_ans != 'n') {
        std::cout << "Apply LP filter? [y/n]: " << std::flush; std::cin >> filt_ans; }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    bool use_filter = (filt_ans == 'y');

    double filter_cutoff = cutoff_MHz;
    if (use_filter && cutoff_MHz <= 0) {
        std::cout << "  Loaded calibration has no filter (cut0).\n"
                  << "  LP cutoff to apply [MHz]: " << std::flush;
        std::cin >> filter_cutoff;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (filter_cutoff <= 0) {
            std::cerr << "  Invalid cutoff, disabling filter.\n";
            use_filter = false;
        }
    }

    char quiet_ans = 0;
    while (quiet_ans != 'y' && quiet_ans != 'n') {
        std::cout << "Pre-laser quiet check (reject dark counts before laser)? [y/n]: "
                  << std::flush;
        std::cin >> quiet_ans;
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    bool do_quiet_check = (quiet_ans == 'y');

    std::cout << "\n  -- Configuration --\n"
              << "  Vbias:           " << vbias << " V\n"
              << "  Calibration:     cut" << (int)cutoff_MHz << " MHz"
              << "  (gain=" << cal.m << " mV/p.e.)\n"
              << "  LET:             " << frac_pe << " p.e. = "
              << let_thr << " mV\n"
              << "  LP filter:       " << (use_filter ? "ON" : "OFF");
    if (use_filter) std::cout << " @ " << filter_cutoff << " MHz";
    std::cout << "\n"
              << "  Pre-laser check: " << (do_quiet_check ? "ON" : "OFF")
              << " (thr=" << let_thr << " mV in ["
              << BASELINE_END << ", t_laser-1] ns)\n\n";

    // ── Find run files ───────────────────────────────────────
    std::string run_pattern = "data.vbias_" + std::to_string(vbias) + "_run_";
    std::map<int,std::string> foundRuns;
    {
        void* dirp = gSystem->OpenDirectory(dataDir.c_str());
        if (dirp) {
            const char* entry;
            while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
                std::string fname(entry);
                if (fname.find(run_pattern) == std::string::npos) continue;
                if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
                size_t pos = fname.find("_run_"); if (pos == std::string::npos) continue;
                try { std::string sub = fname.substr(pos+5); size_t dot = sub.find(".root");
                    if (dot != std::string::npos) sub = sub.substr(0,dot);
                    foundRuns[std::stoi(sub)] = dataDir + "/" + fname; } catch(...) {}
            }
            gSystem->FreeDirectory(dirp);
        }
    }

    if (foundRuns.empty()) {
        std::cerr << "No run files found for Vbias=" << vbias << ".\n";
        return;
    }

    // ── Estimate laser threshold at 5% ───────────────────────
    double laser_thr = cal.laser_thr;
    if (laser_thr <= 0) laser_thr = 10.0;
    {
        TFile* f = TFile::Open(foundRuns.begin()->second.c_str(), "READ");
        if (f && !f->IsZombie()) {
            TTree* tr = (TTree*)f->Get("laser");
            if (tr) {
                Double_t tL[N], aL[N];
                tr->SetBranchAddress("time", tL);
                tr->SetBranchAddress("amplitude", aL);
                std::vector<double> peaks;
                Long64_t nS = std::min((Long64_t)500, tr->GetEntries());
                for (Long64_t i = 0; i < nS; ++i) {
                    tr->GetEntry(i);
                    std::vector<double> pre;
                    for (int j = 0; j < N; ++j) if (tL[j] < BASELINE_END) pre.push_back(aL[j]);
                    double off = 0;
                    if (!pre.empty()) { auto tmp = pre; std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end()); off = tmp[tmp.size()/2]; }
                    double pk = -1e9;
                    for (int j = 0; j < N; ++j) if (aL[j]-off > pk) pk = aL[j]-off;
                    if (pk > 0) peaks.push_back(pk);
                }
                if (!peaks.empty()) {
                    std::nth_element(peaks.begin(), peaks.begin()+peaks.size()/2, peaks.end());
                    double med = peaks[peaks.size()/2];
                    laser_thr = std::max(10.0, std::min(med * 0.05, 50.0));
                }
            }
            f->Close(); delete f;
        }
    }
    std::cout << "  Laser threshold = " << std::fixed << std::setprecision(1)
              << laser_thr << " mV (5%)\n";

    // Sampling rate
    double fs_MHz = 5000.0;
    {
        TFile* f = TFile::Open(foundRuns.begin()->second.c_str(), "READ");
        if (f && !f->IsZombie()) {
            TTree* tr = (TTree*)f->Get("ch1");
            if (tr) { Double_t tb[N]; tr->SetBranchAddress("time", tb); tr->GetEntry(0);
                fs_MHz = 1000.0 / (tb[1] - tb[0]); }
            f->Close(); delete f;
        }
    }
    std::cout << "  fs = " << fs_MHz << " MHz\n";

    const double edge_thr_frac = use_filter ? 0.5 : 100.0;

    // ── Collectors ───────────────────────────────────────────
    std::vector<SavedWF> wf_neg;
    std::vector<SavedWF> wf_pos;
    const int MAX_WF = 10;

    TH2D* hPersist = new TH2D("hPersist",
        Form("Persistence (all accepted, aligned on t_{laser})   Vbias=%dV  LET=%.2f p.e.  %s"
             ";t - t_{laser} (ns);Amplitude (mV)",
             vbias, frac_pe, use_filter ? Form("LP=%dMHz",(int)filter_cutoff) : "NO FILT"),
        1024, -50, 200, 440, -20, 200);
    hPersist->SetDirectory(nullptr);

    long totalEvt = 0, nAccepted = 0, nNoLaser = 0, nNoTOT = 0;
    long nNegDt = 0, nPosDt = 0, nDirtyBL = 0, nPreLaserReject = 0;
    const double bl_max_rms = use_filter ? 2.0 : 4.0;

    // ── File-by-file event loop ──────────────────────────────
    std::cout << "\n  Processing " << foundRuns.size() << " run files...\n";

    for (auto& [runNum, path] : foundRuns) {
        TFile* fIn = TFile::Open(path.c_str(), "READ");
        if (!fIn || fIn->IsZombie()) { delete fIn; continue; }

        TTree* treeCh1 = (TTree*)fIn->Get("ch1");
        TTree* treeLaser = (TTree*)fIn->Get("laser");
        if (!treeCh1 || !treeLaser) { fIn->Close(); delete fIn; continue; }

        Double_t t1[N], a1[N], tL2[N], aL2[N];
        treeCh1->SetBranchAddress("time", t1);
        treeCh1->SetBranchAddress("amplitude", a1);
        treeLaser->SetBranchAddress("time", tL2);
        treeLaser->SetBranchAddress("amplitude", aL2);
        treeCh1->SetCacheSize(2*1024*1024);
        treeLaser->SetCacheSize(2*1024*1024);

        Long64_t nFile = treeCh1->GetEntries();
        std::cout << "  Run " << runNum << ": " << nFile << " events" << std::flush;

        for (Long64_t i = 0; i < nFile; ++i) {
            if (i % 10000 == 0 && gROOT->IsInterrupted()) break;
            ++totalEvt;

            treeCh1->GetEntry(i);
            treeLaser->GetEntry(i);

            double t_laser = laserTriggerTime(tL2, aL2, N, laser_thr);
            if (t_laser < -900.0) { ++nNoLaser; continue; }

            std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
            bool bl_ok = true;
            std::vector<double> v_bc = correctBaseline(v_t, v_a,
                BASELINE_START, BASELINE_END, bl_max_rms, &bl_ok);
            if (!bl_ok) { ++nDirtyBL; continue; }
            std::vector<double> af = use_filter
                ? butterworthLowPass(v_bc, filter_cutoff, fs_MHz)
                : v_bc;

            // Pre-laser quiet check
            if (do_quiet_check && !preLaserQuiet(v_t, af, t_laser, let_thr)) {
                ++nPreLaserReject;
                continue;
            }

            double amp_max = *std::max_element(af.begin(), af.end());

            // Use computeTOT from TOTAnalysis.h (with all v2 fixes)
            auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                0, N-1, 100, edge_thr_frac, 0.5, 10, 0.3, 3, 50, false, 100);
            if (t_rise < 0 || t_fall < 0) { ++nNoTOT; continue; }

            double tot = t_fall - t_rise;
            double delta_t = t_rise - t_laser;
            if (tot <= 0 || tot >= 150) continue;

            int n_pe = estimatePE(amp_max, cal);
            ++nAccepted;

            // Persistence: all accepted waveforms
            for (int j = 0; j < N; ++j)
                hPersist->Fill(v_t[j] - t_laser, af[j]);

            if (delta_t < 0) {
                ++nNegDt;
                if ((int)wf_neg.size() < MAX_WF) {
                    wf_neg.push_back({v_t, af, v_bc, t_laser, t_rise, t_fall,
                                      delta_t, tot, amp_max, let_thr, n_pe, i, runNum});
                }
            } else {
                ++nPosDt;
                if ((int)wf_pos.size() < MAX_WF) {
                    wf_pos.push_back({v_t, af, v_bc, t_laser, t_rise, t_fall,
                                      delta_t, tot, amp_max, let_thr, n_pe, i, runNum});
                }
            }
        }

        treeCh1->SetCacheSize(0);
        treeLaser->SetCacheSize(0);
        fIn->Close(); delete fIn;
        std::cout << "  acc=" << nAccepted << "\n";
    }

    std::cout << "\n  -- Results --\n"
              << "  Total:           " << totalEvt << "\n"
              << "  DirtyBL:         " << nDirtyBL << "\n"
              << "  PreLaserReject:  " << nPreLaserReject << "\n"
              << "  NoTOT:           " << nNoTOT << "\n"
              << "  NoLaser:         " << nNoLaser << "\n"
              << "  Accepted:        " << nAccepted << "\n"
              << "  dt < 0:          " << nNegDt << " (" << std::fixed << std::setprecision(1)
              << (nAccepted > 0 ? 100.0*nNegDt/nAccepted : 0) << "%)\n"
              << "  dt > 0:          " << nPosDt << "\n\n";

    // ── 1. Waveforms with dt < 0 ────────────────────────────
    drawWaveformGrid(wf_neg,
        Form("Waveforms with #Deltat < 0  (Vbias=%dV  LET=%.2f p.e.  laser_thr=%.1f mV)",
             vbias, frac_pe, laser_thr),
        Form("wf_neg_dt_vbias%d_let%.2fpe", vbias, frac_pe),
        ctx);

    // ── 2. Waveforms with dt > 0 ────────────────────────────
    drawWaveformGrid(wf_pos,
        Form("Waveforms with #Deltat > 0  (Vbias=%dV  LET=%.2f p.e.  laser_thr=%.1f mV)",
             vbias, frac_pe, laser_thr),
        Form("wf_pos_dt_vbias%d_let%.2fpe", vbias, frac_pe),
        ctx);

    // ── 3. Persistence plot (all accepted) ───────────────────
    if (nAccepted > 0) {
        TCanvas* cP = new TCanvas("cPersist",
            Form("Persistence (all accepted)  Vbias=%dV  LET=%.2f p.e.  %s",
                 vbias, frac_pe,
                 do_quiet_check ? "quiet=ON" : "quiet=OFF"),
            1100, 600);
        cP->SetLeftMargin(PAD_LEFT);
        cP->SetRightMargin(0.13f);
        cP->SetBottomMargin(PAD_BOTTOM);
        cP->SetTopMargin(PAD_TOP);
        cP->SetGrid();
        cP->SetLogz(1);

        gStyle->SetPalette(kBird);
        hPersist->SetContour(100);
        hPersist->GetXaxis()->SetTitleSize(0.048f);
        hPersist->GetYaxis()->SetTitleSize(0.048f);
        hPersist->GetZaxis()->SetTitle("Counts / bin");
        hPersist->Draw("COLZ");

        TLine* lLas = new TLine(0, -20, 0, 200);
        lLas->SetLineColor(kGreen+2); lLas->SetLineStyle(2); lLas->SetLineWidth(2);
        lLas->Draw("same");

        TLine* lThr = new TLine(-50, let_thr, 200, let_thr);
        lThr->SetLineColor(kRed+1); lThr->SetLineStyle(2); lThr->SetLineWidth(2);
        lThr->Draw("same");

        TPaveText* pt = new TPaveText(0.15, 0.66, 0.50, 0.88, "NDC");
        pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextFont(42); pt->SetTextSize(0.034);
        pt->AddText(Form("N accepted = %ld", nAccepted));
        pt->AddText(Form("#Deltat<0 = %ld (%.1f%%)", nNegDt,
                         nAccepted > 0 ? 100.0*nNegDt/nAccepted : 0.));
        if (do_quiet_check)
            pt->AddText(Form("Pre-laser rejected = %ld", nPreLaserReject));
        pt->AddText(Form("laser thr = %.1f mV", laser_thr));
        pt->AddText(Form("LET thr = %.1f mV", let_thr));
        pt->AddText(Form("post-fall limit = %.1f mV (= LET)", let_thr));
        pt->AddText("confirm timeout = 100 samp (20 ns)");
        pt->AddText(Form("end-of-trace limit = %.1f mV (last 10 ns)", let_thr * 0.5));
        pt->Draw();

        TLegend* leg = new TLegend(0.55, 0.78, 0.85, 0.88);
        leg->SetBorderSize(1); leg->SetFillColor(0); leg->SetTextSize(0.033);
        leg->AddEntry(lLas, "t_{laser} = 0", "l");
        leg->AddEntry(lThr, Form("LET = %.1f mV", let_thr), "l");
        leg->Draw();

        cP->Update(); cP->Modified();
        ctx.savePNG(cP, Form("persistence_vbias%d_let%.2fpe.png", vbias, frac_pe));
    } else {
        delete hPersist;
    }

    std::cout << "\n+==========================================================+\n"
              << "|  DONE  PNG: " << ctx.pngDir << "\n"
              << "+==========================================================+\n"
              << "\n  Canvases are open. Use the ROOT prompt to inspect them.\n"
              << "  Type .q to exit.\n";
}
