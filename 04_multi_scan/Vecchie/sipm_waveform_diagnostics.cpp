/**
 * sipm_waveform_diagnostics.cpp
 * ==============================
 * Waveform-level diagnostic tool for SiPM TOT analysis.
 *
 * After applying the Butterworth low-pass filter and computing
 * t_laser, t_rise, t_fall for each event, produces four diagnostic
 * canvas sets:
 *
 *  Canvas 1 — ALL waveforms aligned on t_laser (persistence + grid)
 *             Δt histogram of all events with a valid t_laser.
 *
 *  Canvas 2 — Waveforms WITH a threshold crossing (t_rise found),
 *             regardless of whether t_fall is also found.
 *             Δt histogram (t_rise - t_laser).
 *
 *  Canvas 3 — Waveforms with VALID TOT (both t_rise AND t_fall found
 *             inside the acquisition window).
 *             Δt histogram.
 *
 *  Canvas 4 — Waveforms with t_rise found but NO valid t_fall
 *             (crossing but truncated / ringing / no return below thr).
 *             Δt histogram.
 *
 * Each canvas set = one persistence TH2D (time vs amplitude, log Z)
 * + one Δt TH1D, side by side on the same canvas.
 *
 * Waveforms are aligned on t_laser (x axis = t - t_laser).
 * A separate grid canvas shows up to MAX_GRID individual waveforms
 * per category (5x4 layout, same as WaveformPlotter).
 *
 * Compile: .L sipm_waveform_diagnostics.cpp+
 * Run:     sipm_waveform_diagnostics()
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/CalibIO.h"
#include "../header/TOTAnalysis.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TAxis.h>

// ============================================================
//  Per-category accumulator
//  Holds a persistence TH2D and a Δt TH1D.
//  Also stores up to MAX_GRID waveforms for the grid canvas.
// ============================================================
static const int MAX_GRID = 100;   // max individual waveforms per category

struct DiagCategory {
    std::string name;           // short label, e.g. "all", "crossing", etc.
    std::string title;          // long description for canvas title

    // Persistence map: x = t - t_laser [ns], y = amplitude [mV]
    TH2D* hPers  = nullptr;
    // Delta_t histogram: x = t_rise - t_laser [ns]  (or t_laser jitter for cat 1)
    TH1D* hDt    = nullptr;

    // Individual waveforms for grid canvas (aligned on t_laser)
    struct WF {
        std::vector<double> t_rel;  // t - t_laser
        std::vector<double> amp;
        double delta_t = 0;         // t_rise - t_laser  (-999 if no crossing)
        double tot     = 0;         // t_fall - t_rise   (-1 if not valid)
    };
    std::vector<WF> grid_wf;

    long n_events = 0;

    void init(const std::string& n, const std::string& ttl,
              double t_lo, double t_hi, int t_bins,
              double amp_lo, double amp_hi, int amp_bins,
              double dt_lo, double dt_hi, int dt_bins)
    {
        name  = n;
        title = ttl;
        hPers = new TH2D(
            Form("hPers_%s", n.c_str()), "",
            t_bins,   t_lo,   t_hi,
            amp_bins, amp_lo, amp_hi);
        hPers->SetDirectory(nullptr);

        hDt = new TH1D(
            Form("hDt_%s", n.c_str()), "",
            dt_bins, dt_lo, dt_hi);
        hDt->SetDirectory(nullptr);
    }

    void fillPers(const std::vector<double>& t_rel_v,
                  const std::vector<double>& amp_v)
    {
        int ns = (int)std::min(t_rel_v.size(), amp_v.size());
        for (int j = 0; j < ns; ++j)
            hPers->Fill(t_rel_v[j], amp_v[j]);
    }

    void fillDt(double dt) { if (hDt) hDt->Fill(dt); }

    void collectGrid(const std::vector<double>& t_rel_v,
                     const std::vector<double>& amp_v,
                     double delta_t, double tot)
    {
        if ((int)grid_wf.size() >= MAX_GRID) return;
        WF w;
        w.t_rel   = t_rel_v;
        w.amp     = amp_v;
        w.delta_t = delta_t;
        w.tot     = tot;
        grid_wf.push_back(std::move(w));
    }
};

// ============================================================
//  Draw one summary canvas: persistence (left) + Δt histo (right)
// ============================================================
static void drawSummaryCanvas(DiagCategory& cat,
                               double let_thr,
                               OutCtx& ctx)
{
    if (cat.n_events == 0) {
        std::cout << "  [Diag] Category '" << cat.name
                  << "': no events — skipping.\n";
        return;
    }

    cat.hPers->SetTitle(
        Form("Persistence (aligned t_{laser})   %s   N=%ld"
             ";t #minus t_{laser} (ns);Amplitude (mV)",
             cat.title.c_str(), cat.n_events));
    cat.hDt->SetTitle(
        Form("#Deltat distribution   %s"
             ";t_{rise} #minus t_{laser} (ns);Events",
             cat.title.c_str()));

    TCanvas* cSum = new TCanvas(
        Form("cSummary_%s", cat.name.c_str()),
        Form("Diagnostics: %s", cat.title.c_str()),
        1600, 600);
    cSum->Divide(2, 1);

    // Left: persistence
    cSum->cd(1);
    gPad->SetGrid();
    gPad->SetLogz(1);
    gPad->SetLeftMargin(PAD_LEFT);
    gPad->SetRightMargin(0.14f);
    gPad->SetBottomMargin(PAD_BOTTOM);
    gPad->SetTopMargin(PAD_TOP);
    gStyle->SetPalette(kBird);
    cat.hPers->SetContour(100);
    cat.hPers->GetXaxis()->SetTitleSize(0.048f);
    cat.hPers->GetYaxis()->SetTitleSize(0.048f);
    cat.hPers->GetXaxis()->SetLabelSize(0.040f);
    cat.hPers->GetYaxis()->SetLabelSize(0.040f);
    cat.hPers->GetZaxis()->SetTitle("Counts/bin");
    cat.hPers->Draw("COLZ");

    // Marker lines on persistence
    double t_lo = cat.hPers->GetXaxis()->GetXmin();
    double t_hi = cat.hPers->GetXaxis()->GetXmax();
    double a_lo = cat.hPers->GetYaxis()->GetXmin();
    double a_hi = cat.hPers->GetYaxis()->GetXmax();

    // t_laser at 0 (green)
    TLine* lL = new TLine(0, a_lo, 0, a_hi);
    lL->SetLineColor(kGreen+2); lL->SetLineStyle(2); lL->SetLineWidth(2);
    lL->Draw("same");

    // LET threshold (blue horizontal)
    TLine* lT = new TLine(t_lo, let_thr, t_hi, let_thr);
    lT->SetLineColor(kAzure+1); lT->SetLineStyle(2); lT->SetLineWidth(2);
    lT->Draw("same");

    // Median Δt (red) if hDt has entries
    if (cat.hDt->GetEntries() > 0) {
        double median = 0;
        double total  = cat.hDt->Integral();
        double cumul  = 0;
        for (int b = 1; b <= cat.hDt->GetNbinsX(); ++b) {
            cumul += cat.hDt->GetBinContent(b);
            if (cumul >= total * 0.5) {
                median = cat.hDt->GetBinCenter(b);
                break;
            }
        }
        TLine* lM = new TLine(median, a_lo, median, a_hi);
        lM->SetLineColor(kRed+1); lM->SetLineStyle(2); lM->SetLineWidth(2);
        lM->Draw("same");
    }

    // Right: Δt histogram
    cSum->cd(2);
    gPad->SetGrid();
    gPad->SetLeftMargin(PAD_LEFT);
    gPad->SetRightMargin(PAD_RIGHT);
    gPad->SetBottomMargin(PAD_BOTTOM);
    gPad->SetTopMargin(PAD_TOP);
    cat.hDt->SetLineColor(kAzure+1);
    cat.hDt->SetLineWidth(2);
    cat.hDt->SetFillColor(kAzure-9);
    cat.hDt->GetXaxis()->SetTitleSize(0.048f);
    cat.hDt->GetYaxis()->SetTitleSize(0.048f);
    cat.hDt->Draw("HIST");

    // Stats box
    TPaveText* pt = new TPaveText(0.55, 0.72, 0.94, 0.88, "NDC");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
    pt->SetTextFont(42); pt->SetTextSize(0.038);
    pt->AddText(Form("N = %ld", cat.n_events));
    if (cat.hDt->GetEntries() > 0) {
        pt->AddText(Form("Mean  = %.2f ns", cat.hDt->GetMean()));
        pt->AddText(Form("RMS   = %.2f ns", cat.hDt->GetRMS()));
    }
    pt->Draw();

    cSum->Update(); cSum->Modified();
    ctx.savePNG(cSum, Form("wfdiag_summary_%s.png", cat.name.c_str()));
}

// ============================================================
//  Draw grid of individual waveforms (5x4 per page)
// ============================================================
static void drawGridCanvas(DiagCategory& cat,
                            double let_thr,
                            OutCtx& ctx)
{
    if (cat.grid_wf.empty()) return;

    const int COLS = 5, ROWS = 4, PER_PAGE = COLS * ROWS;
    int n_pages = ((int)cat.grid_wf.size() + PER_PAGE - 1) / PER_PAGE;

    std::cout << "  [Diag] Grid '" << cat.name << "': "
              << cat.grid_wf.size() << " waveforms, "
              << n_pages << " page(s)\n";

    for (int pg = 0; pg < n_pages; ++pg) {
        int first = pg * PER_PAGE;
        int last  = std::min(first + PER_PAGE, (int)cat.grid_wf.size());
        int count = last - first;

        TCanvas* cG = new TCanvas(
            Form("cGrid_%s_p%d", cat.name.c_str(), pg),
            Form("Waveforms: %s  (page %d/%d)",
                 cat.title.c_str(), pg+1, n_pages),
            300*COLS, 220*ROWS);
        cG->SetFillColor(0);
        cG->Divide(COLS, ROWS, 0.002, 0.002);

        for (int k = 0; k < count; ++k) {
            const auto& wf = cat.grid_wf[first + k];
            cG->cd(k + 1);
            gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03);
            gPad->SetBottomMargin(0.18); gPad->SetTopMargin(0.15);
            gPad->SetGrid(); gPad->SetFillColor(0);

            if (wf.t_rel.empty()) continue;

            // Y range
            double amp_max = *std::max_element(wf.amp.begin(), wf.amp.end());
            double amp_min = *std::min_element(wf.amp.begin(), wf.amp.end());
            double y_lo = amp_min - 5.0;
            double y_hi = std::max(amp_max * 1.25, let_thr * 2.5);

            TGraph* gr = new TGraph((int)wf.t_rel.size(),
                                    wf.t_rel.data(), wf.amp.data());
            gr->SetTitle("");
            gr->SetLineColor(kAzure+1); gr->SetLineWidth(1);
            gr->GetXaxis()->SetTitle("t#minust_{laser} (ns)");
            gr->GetYaxis()->SetTitle("A (mV)");
            gr->GetXaxis()->SetTitleSize(0.08); gr->GetYaxis()->SetTitleSize(0.08);
            gr->GetXaxis()->SetLabelSize(0.07); gr->GetYaxis()->SetLabelSize(0.07);
            gr->GetXaxis()->SetTitleOffset(0.9); gr->GetYaxis()->SetTitleOffset(0.85);
            gr->GetXaxis()->SetLimits(wf.t_rel.front(), wf.t_rel.back());
            gr->GetYaxis()->SetRangeUser(y_lo, y_hi);
            gr->Draw("AL");

            // LET threshold (blue horizontal)
            TLine* lT = new TLine(wf.t_rel.front(), let_thr,
                                   wf.t_rel.back(),  let_thr);
            lT->SetLineColor(kAzure+1); lT->SetLineStyle(2);
            lT->SetLineWidth(1); lT->Draw("same");

            // t_laser at 0 (green)
            TLine* lL = new TLine(0, y_lo, 0, y_hi * 0.85);
            lL->SetLineColor(kGreen+2); lL->SetLineStyle(2);
            lL->SetLineWidth(1); lL->Draw("same");

            // t_rise (red) — if crossing found
            if (wf.delta_t > -990) {
                TLine* lR = new TLine(wf.delta_t, y_lo,
                                       wf.delta_t, y_hi * 0.85);
                lR->SetLineColor(kRed+1); lR->SetLineStyle(2);
                lR->SetLineWidth(1); lR->Draw("same");
            }

            // t_fall (orange) — if TOT valid
            if (wf.tot > 0 && wf.delta_t > -990) {
                double t_fall_rel = wf.delta_t + wf.tot;
                TLine* lF = new TLine(t_fall_rel, y_lo,
                                       t_fall_rel, y_hi * 0.70);
                lF->SetLineColor(kOrange+7); lF->SetLineStyle(2);
                lF->SetLineWidth(1); lF->Draw("same");
            }

            // Title box
            TPaveText* pt = new TPaveText(0.01, 0.86, 0.99, 0.99, "NDC");
            pt->SetBorderSize(0); pt->SetFillColor(0); pt->SetFillStyle(0);
            pt->SetTextFont(42); pt->SetTextSize(0.072); pt->SetTextAlign(22);
            if (wf.delta_t > -990 && wf.tot > 0)
                pt->AddText(Form("#Deltat=%.1f  TOT=%.1f ns",
                                  wf.delta_t, wf.tot));
            else if (wf.delta_t > -990)
                pt->AddText(Form("#Deltat=%.1f ns  no t_{fall}",
                                  wf.delta_t));
            else
                pt->AddText("no crossing");
            pt->Draw();
        }

        // Blank remaining pads
        for (int k = count; k < PER_PAGE; ++k) {
            cG->cd(k + 1);
            gPad->SetFillColor(kWhite);
        }

        cG->Update(); cG->Modified();
        ctx.savePNG(cG, Form("wfdiag_grid_%s_p%d.png",
                              cat.name.c_str(), pg+1));
    }
}

// ============================================================
//  MAIN
// ============================================================
void sipm_waveform_diagnostics()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();

    std::cout << "\n+==========================================================+\n"
              << "|  SiPM WAVEFORM DIAGNOSTICS\n"
              << "+==========================================================+\n";

    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        std::cout << prompt << std::flush;
        while (line.empty()) std::getline(std::cin, line);
        return line;
    };

    // --- 1. Parameters ------------------------------------------
    int vbias = 0;
    std::cout << "Vbias [V]: " << std::flush;
    std::cin >> vbias;

    double cutoff_MHz = 0;
    std::cout << "Low-pass cutoff [MHz] (0 = no filter): " << std::flush;
    std::cin >> cutoff_MHz;

    const std::string dataDir = DATA_DIR;

    // --- 2. Load calibration ------------------------------------
    // Find available calibrations for this Vbias
    std::vector<double> cutoffs;
    {
        const std::string prefix = "calib_vbias" + std::to_string(vbias) + "_cut";
        const std::string suffix = "mhz.root";
        void* dirp = gSystem->OpenDirectory(dataDir.c_str());
        if (dirp) {
            const char* entry;
            while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
                std::string fname(entry);
                if (fname.size() < prefix.size() + suffix.size()) continue;
                if (fname.substr(0, prefix.size()) != prefix) continue;
                if (fname.substr(fname.size()-suffix.size()) != suffix) continue;
                std::string mid = fname.substr(prefix.size(),
                    fname.size()-prefix.size()-suffix.size());
                try { double co = std::stod(mid); if (co>0) cutoffs.push_back(co); }
                catch (...) {}
            }
            gSystem->FreeDirectory(dirp);
        }
        std::sort(cutoffs.begin(), cutoffs.end());
    }

    CalibResult cal;
    bool calLoaded = false;

    if (!cutoffs.empty()) {
        std::cout << "  Available calibrations: ";
        for (double co : cutoffs) std::cout << (int)co << " MHz  ";
        double chosen = cutoffs.back();
        if (cutoffs.size() > 1) {
            std::cout << "\n  Which cutoff [MHz]? [default "
                      << (int)chosen << "]: " << std::flush;
            std::string line; std::getline(std::cin, line);
            if (!line.empty()) {
                try {
                    double inp = std::stod(line);
                    double best = chosen, bestDiff = 1e9;
                    for (double co : cutoffs)
                        if (std::abs(co-inp) < bestDiff)
                            { bestDiff=std::abs(co-inp); best=co; }
                    chosen = best;
                } catch (...) {}
            }
        } else { std::cout << "\n"; }
        calLoaded = loadCalibration(cal, vbias, chosen, dataDir);
    }

    if (!calLoaded) {
        std::cout << "  [WARN] No calibration loaded — using defaults.\n"
                  << "  laser_thr=10 mV  gain=1  offset=0\n";
        cal.laser_thr  = 10.0;
        cal.m          = 1.0;
        cal.q          = 0.0;
        cal.ok         = true;
    }

    // LET threshold
    double frac_pe = 0.5;
    std::cout << "LET threshold [p.e.]: " << std::flush;
    std::cin >> frac_pe;
    double let_thr = cal.q + frac_pe * cal.m;
    std::cout << "  LET threshold = " << std::fixed << std::setprecision(2)
              << let_thr << " mV\n";

    // Max events to process (safety valve)
    long max_events = 0;
    std::cout << "Max events to process (0 = all): " << std::flush;
    std::cin >> max_events;

    // --- 3. Find and chain run files ----------------------------
    const std::string pattern = "data.vbias_" + std::to_string(vbias) + "_run_";
    std::map<int, std::string> foundRuns;
    {
        void* dirp = gSystem->OpenDirectory(dataDir.c_str());
        if (!dirp) { std::cerr << "[ERROR] Cannot open " << dataDir << "\n"; return; }
        const char* entry;
        while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
            std::string fname(entry);
            if (fname.find(pattern) == std::string::npos) continue;
            if (fname.size()<5 || fname.substr(fname.size()-5)!=".root") continue;
            size_t pos = fname.find("_run_");
            if (pos == std::string::npos) continue;
            try {
                std::string sub = fname.substr(pos+5);
                size_t dot = sub.find(".root");
                if (dot!=std::string::npos) sub=sub.substr(0,dot);
                foundRuns[std::stoi(sub)] = dataDir + "/" + fname;
            } catch (...) {}
        }
        gSystem->FreeDirectory(dirp);
    }
    if (foundRuns.empty()) {
        std::cerr << "[ERROR] No run files for Vbias=" << vbias << "\n";
        return;
    }
    std::cout << "  Found " << foundRuns.size() << " run file(s)\n";

    TChain* chainCh1   = new TChain("ch1");
    TChain* chainLaser = new TChain("laser");
    for (auto& [r, path] : foundRuns) {
        chainCh1  ->Add(path.c_str());
        chainLaser->Add(path.c_str());
    }
    Long64_t nEntries = chainCh1->GetEntries();
    if (max_events > 0 && (Long64_t)max_events < nEntries)
        nEntries = (Long64_t)max_events;
    std::cout << "  Processing " << nEntries << " events\n";

    // Sampling rate
    double fs_MHz = 0;
    {
        TFile* f0 = TFile::Open(foundRuns.begin()->second.c_str(), "READ");
        if (f0 && !f0->IsZombie()) {
            TTree* tr = (TTree*)f0->Get("ch1");
            if (tr) {
                const int N0=1024; Double_t tb[N0];
                tr->SetBranchAddress("time", tb); tr->GetEntry(0);
                fs_MHz = 1000.0/(tb[1]-tb[0]);
            }
            f0->Close();
        }
    }
    if (fs_MHz <= 0) { std::cerr << "[ERROR] Cannot read fs.\n"; return; }
    std::cout << "  Sampling rate: " << fs_MHz << " MHz\n";

    // --- 4. Histogram / persistence parameters ------------------
    // Time axis: full waveform relative to t_laser
    // At 5 GS/s, 1024 samples = 204.6 ns.
    // t_laser is typically at ~25 ns -> range [-30, 180] ns covers all.
    const double T_LO   = -50.0, T_HI   = 200.0;
    const int    T_BINS = 1024;
    const double A_LO   = -20.0, A_HI   = 200.0;
    const int    A_BINS = 440;
    const double DT_LO  = -60.0, DT_HI  = 220.0;
    const int    DT_BINS = 1400;   // 0.2 ns/bin

    // --- 5. Initialise categories -------------------------------
    DiagCategory cat1, cat2, cat3, cat4;
    cat1.init("all",      "All events (t_laser found)",
              T_LO, T_HI, T_BINS, A_LO, A_HI, A_BINS,
              DT_LO, DT_HI, DT_BINS);
    cat2.init("crossing", "Events with threshold crossing (t_rise found)",
              T_LO, T_HI, T_BINS, A_LO, A_HI, A_BINS,
              DT_LO, DT_HI, DT_BINS);
    cat3.init("valid_tot","Events with valid TOT (t_rise AND t_fall found)",
              T_LO, T_HI, T_BINS, A_LO, A_HI, A_BINS,
              DT_LO, DT_HI, DT_BINS);
    cat4.init("no_tfall", "Crossing found but t_fall missing",
              T_LO, T_HI, T_BINS, A_LO, A_HI, A_BINS,
              DT_LO, DT_HI, DT_BINS);

    // --- 6. Branch setup ----------------------------------------
    const int N = 1024;
    Double_t t1[N], a1[N], tL[N], aL[N];
    chainCh1  ->SetBranchAddress("time",      t1);
    chainCh1  ->SetBranchAddress("amplitude", a1);
    chainLaser->SetBranchAddress("time",      tL);
    chainLaser->SetBranchAddress("amplitude", aL);

    long nNoLaser=0, nNoCross=0, nNoFall=0, nValid=0;

    std::cout << "\n  Running event loop...\n";

    // --- 7. Event loop ------------------------------------------
    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 10000 == 0) {
            if (gROOT->IsInterrupted()) break;
            gSystem->ProcessEvents();
            std::cout << "\r  " << i << " / " << nEntries
                      << "  valid=" << nValid
                      << "  noLaser=" << nNoLaser
                      << "  noCross=" << nNoCross
                      << "  noFall=" << nNoFall
                      << "   " << std::flush;
        }

        chainCh1  ->GetEntry(i);
        chainLaser->GetEntry(i);

        // Laser trigger time
        double t_laser = laserTriggerTime(tL, aL, N, cal.laser_thr);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        // Baseline correction + optional filter
        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = correctBaseline(v_t, v_a,
                                                  BASELINE_START, BASELINE_END);
        if (cutoff_MHz > 0)
            af = butterworthLowPass(af, cutoff_MHz, fs_MHz);

        // Relative time axis (aligned on t_laser)
        std::vector<double> t_rel(N);
        for (int j = 0; j < N; ++j) t_rel[j] = v_t[j] - t_laser;

        // --- Category 1: all events with valid t_laser ----------
        cat1.n_events++;
        cat1.fillPers(t_rel, af);
        cat1.fillDt(t_laser);
        cat1.collectGrid(t_rel, af, -999.0, -1.0);

        // Edge-stability check — same logic as computeTOT.
        // Rejects events where the signal was already elevated at
        // acquisition start or had not returned to baseline at the end.
        const int    N_EDGE       = 100;
        const double EDGE_FRAC    = 0.5;
        const double edge_limit   = EDGE_FRAC * let_thr;

        double med_first = edgeMedian(af, 0, N_EDGE);
        double med_last  = edgeMedian(af, N - N_EDGE, N_EDGE);
        bool   edge_ok   = (med_first < edge_limit) && (med_last < edge_limit);

        // --- Find t_rise: first crossing above let_thr ----------
        // Search the full waveform (no trigger window restriction here —
        // this is the diagnostic, we want to see everything)
        double t_rise = -1.0, t_fall = -1.0;
        for (int j = 1; j < N; ++j) {
            if (af[j-1] < let_thr && af[j] >= let_thr) {
                t_rise = v_t[j-1] + (let_thr - af[j-1])
                         * (v_t[j]-v_t[j-1]) / (af[j]-af[j-1]);
                break;
            }
        }

        if (t_rise < 0) {
            ++nNoCross;
            continue;
        }

        // Events that fail the edge check are counted but not entered
        // in cat2/3/4 — they would pollute the crossing/TOT categories.
        if (!edge_ok) {
            ++nNoCross;   // count as "no valid crossing" for summary
            continue;
        }

        double delta_t = t_rise - t_laser;

        // --- Category 2: any crossing ---------------------------
        cat2.n_events++;
        cat2.fillPers(t_rel, af);
        cat2.fillDt(delta_t);

        // --- Find t_fall: first return below let_thr after t_rise
        for (int j = 1; j < N; ++j) {
            if (v_t[j] <= t_rise) continue;
            if (af[j-1] >= let_thr && af[j] < let_thr) {
                t_fall = v_t[j-1] + (let_thr - af[j-1])
                         * (v_t[j]-v_t[j-1]) / (af[j]-af[j-1]);
                break;
            }
        }

        if (t_fall > 0) {
            // --- Category 3: valid TOT --------------------------
            double tot = t_fall - t_rise;
            ++nValid;
            cat3.n_events++;
            cat3.fillPers(t_rel, af);
            cat3.fillDt(delta_t);
            cat3.collectGrid(t_rel, af, delta_t, tot);
            // also add to cat2 grid
            cat2.collectGrid(t_rel, af, delta_t, tot);
        } else {
            // --- Category 4: crossing but no t_fall -------------
            ++nNoFall;
            cat4.n_events++;
            cat4.fillPers(t_rel, af);
            cat4.fillDt(delta_t);
            cat4.collectGrid(t_rel, af, delta_t, -1.0);
            // also add to cat2 grid
            cat2.collectGrid(t_rel, af, delta_t, -1.0);
        }
    }

    std::cout << "\n\n  +----- Event loop summary -----\n"
              << "  | Total processed : " << nEntries   << "\n"
              << "  | No laser        : " << nNoLaser   << "\n"
              << "  | No crossing     : " << nNoCross   << "\n"
              << "  | Valid TOT       : " << nValid     << "\n"
              << "  | Cross, no fall  : " << nNoFall    << "\n"
              << "  +------------------------------\n\n";

    // Fix cat1 Δt axis label (it shows t_laser, not t_rise - t_laser)
    cat1.hDt->GetXaxis()->SetTitle("t_{laser} (ns)  [jitter distribution]");

    // --- 8. Draw canvases ---------------------------------------
    std::cout << "  Drawing summary canvases...\n";
    drawSummaryCanvas(cat1, let_thr, ctx);
    drawSummaryCanvas(cat2, let_thr, ctx);
    drawSummaryCanvas(cat3, let_thr, ctx);
    drawSummaryCanvas(cat4, let_thr, ctx);

    std::cout << "  Drawing grid canvases...\n";
    drawGridCanvas(cat1, let_thr, ctx);
    drawGridCanvas(cat2, let_thr, ctx);
    drawGridCanvas(cat3, let_thr, ctx);
    drawGridCanvas(cat4, let_thr, ctx);

    std::cout << "\n+==========================================================+\n"
              << "|  DONE\n"
              << "|  PNG: " << ctx.pngDir << "\n"
              << "|  All canvases open. Type .q to exit.\n"
              << "+==========================================================+\n";

    if (gApplication) gApplication->Run(kTRUE);
}
