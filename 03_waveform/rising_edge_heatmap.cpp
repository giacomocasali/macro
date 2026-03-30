// rising_edge_heatmap.cpp
//
// Builds 2D heatmaps (TH2D) of CFD-aligned rising edges from waveforms
// classified as signal (SNR + temporal window criteria).
//
// Alignment strategy:
//   1. Laser trigger time t_laser: first crossing > 10 mV on laser channel
//      (linear interpolation, same as sipm_timing_unfiltered_v3)
//   2. CFD time t_CFD: crossing at 50% of local maximum, scanning
//      BACKWARDS from the peak (fix from v3 -- avoids spurious crossings)
//   3. Aligned time: t_aligned = t - t_CFD  (rising edge always at t=0)
//
// Signal selection (both criteria must be satisfied):
//   SNR = A_max / sigma_noise  >=  snr_threshold
//   t_peak (relative to laser) in [t_sig_lo, t_sig_hi]
//
// Canvases produced:
//   1. Heatmap -- raw baseline-corrected rising edges
//   2. Heatmap -- Butterworth LP-filtered rising edges
//   3. Amplitude spectrum (sanity check, log-y)
//
// Usage:
//   root -l 'rising_edge_heatmap.cpp+'
//   rising_edge_heatmap()

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include "../header/ButterworthFilter.h"  // scipy-equivalent Butterworth IIR
inline std::vector<double> butterworthLP(
    const std::vector<double>& x, double fc, double fs, int order=4)
{ return butterworthLowPass(x, fc, fs, order); }

// ═══════════════════════════════════════════════════════════════
// TERMINAL UTILITIES
// ═══════════════════════════════════════════════════════════════

static std::string repeatStr(const std::string& s, size_t n)
{
    std::string o; o.reserve(s.size() * n);
    for (size_t i = 0; i < n; ++i) o += s;
    return o;
}

static void printBox(const std::vector<std::string>& lines)
{
    size_t w = 0;
    for (const auto& l : lines) w = std::max(w, l.size());
    w += 4;
    std::cout << "\n  +" << repeatStr("-", w) << "+" << std::endl;
    for (const auto& l : lines) {
        size_t pad = w - l.size() - 2;
        std::cout << "  |  " << l << std::string(pad, ' ') << "|" << std::endl;
    }
    std::cout << "  +" << repeatStr("-", w) << "+" << std::endl;
}

static double promptDouble(const std::string& msg, double lo, double hi)
{
    double v = lo - 1.0;
    while (v < lo || v > hi) {
        std::cout << "  " << msg << " [" << lo << " .. " << hi << "]: ";
        std::cin >> v;
        if (std::cin.fail() || v < lo || v > hi) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "  [ERROR] Out of range or invalid input.\n";
            v = lo - 1.0;
        }
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return v;
}

// ═══════════════════════════════════════════════════════════════
// BASELINE CORRECTION  (median of pre-signal region)
// Identical logic to correctWaveforms() in Utils.h
// ═══════════════════════════════════════════════════════════════

static std::vector<double> correctBaseline(
    const Double_t* time,
    const Double_t* amp,
    int N,
    double pre_end,
    double& sigma_out)
{
    std::vector<double> pre;
    pre.reserve(N / 4);
    for (int i = 0; i < N; ++i)
        if (time[i] <= pre_end) pre.push_back(amp[i]);
    if (pre.empty()) {
        int nPre = std::max(10, N / 10);
        for (int i = 0; i < nPre; ++i) pre.push_back(amp[i]);
    }

    std::vector<double> tmp = pre;
    size_t mid = tmp.size() / 2;
    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
    double offset = tmp[mid];

    double mean = 0.0;
    for (double v : pre) mean += v;
    mean /= (double)pre.size();
    double var = 0.0;
    for (double v : pre) var += (v - mean) * (v - mean);
    sigma_out = (pre.size() > 1) ? std::sqrt(var / (double)pre.size()) : 1e-9;
    if (sigma_out < 1e-9) sigma_out = 1e-9;

    std::vector<double> out(N);
    for (int i = 0; i < N; ++i) out[i] = amp[i] - offset;
    return out;
}

// ═══════════════════════════════════════════════════════════════
// BUTTERWORTH LOW-PASS  (4th order, biquad cascade, IIR)
// ═══════════════════════════════════════════════════════════════

// butterworthLowPass is provided by ButterworthFilter.h
// To use zero-phase (sosfiltfilt): butterworthLowPassZP(x, fc, fs, 4)
// To use fast pre-computed 500MHz:  butterworthLP_500MHz(x)

// ═══════════════════════════════════════════════════════════════
// LASER TRIGGER TIME
// First crossing above laser_thr with linear interpolation.
// Same logic used in all your timing analyses.
// ═══════════════════════════════════════════════════════════════

static double laserTriggerTime(
    const Double_t* t_laser,
    const std::vector<double>& a_laser,
    int N,
    double laser_thr = 10.0)
{
    // Compute median baseline over the first 30 ns of the laser channel.
    // Without this correction a DC offset on the laser channel shifts the
    // threshold crossing and introduces a systematic error on all delta-t values.
    std::vector<double> pre;
    for (int j = 0; j < N; ++j)
        if (t_laser[j] < 30.0) pre.push_back(a_laser[j]);
    double offset = 0.0;
    if (!pre.empty()) {
        std::vector<double> tmp = pre;
        std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
        offset = tmp[tmp.size() / 2];
    }
    for (int j = 1; j < N; ++j) {
        double v0 = a_laser[j-1] - offset;
        double v1 = a_laser[j]   - offset;
        if (v1 > laser_thr && v0 <= laser_thr) {
            return t_laser[j-1]
                   + (laser_thr - v0)
                   * (t_laser[j] - t_laser[j-1])
                   / (v1 - v0);
        }
    }
    return -999.0;
}

// ═══════════════════════════════════════════════════════════════
// CFD TIME  (50% of local maximum, scanning BACKWARDS from peak)
// The backward scan is the v3 fix: it guarantees we find the
// rising-edge crossing and not a spurious re-crossing on the tail.
// ═══════════════════════════════════════════════════════════════

static double cfdTime(
    const Double_t* t_arr,
    const std::vector<double>& amp,
    int N,
    double cfd_fraction = 0.5)
{
    auto itMax = std::max_element(amp.begin(), amp.end());
    double A_max = *itMax;
    int    i_peak = (int)std::distance(amp.begin(), itMax);

    if (A_max <= 0.0 || i_peak < 1) return -999.0;

    double thr = cfd_fraction * A_max;

    // Scan backwards from peak
    for (int j = i_peak; j >= 1; --j) {
        if (amp[j] >= thr && amp[j-1] < thr) {
            return t_arr[j-1]
                   + (thr - amp[j-1])
                   * (t_arr[j] - t_arr[j-1])
                   / (amp[j] - amp[j-1]);
        }
    }
    return -999.0;
}

// ═══════════════════════════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════════════════════════

void rising_edge_heatmap()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);

    printBox({"Rising Edge Heatmap  --  CFD alignment + laser trigger"});

    // ── Open file ─────────────────────────────────────────────
    TFile* file = TFile::Open("../data/data.vbias_{55}.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Cannot open data.vbias_{55}.root" << std::endl;
        return;
    }
    TTree* treeCh1   = (TTree*)file->Get("ch1");
    TTree* treeLaser = (TTree*)file->Get("laser");
    if (!treeCh1 || !treeLaser) {
        std::cerr << "[ERROR] TTree 'ch1' or 'laser' not found." << std::endl;
        file->Close(); return;
    }
    Long64_t nEntries = treeCh1->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[ERROR] TTree is empty." << std::endl;
        file->Close(); return;
    }

    const int N = 1024;
    Double_t t1[N], a1[N], tL[N], aL[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);
    treeLaser->SetBranchAddress("time",      tL);
    treeLaser->SetBranchAddress("amplitude", aL);

    // Sampling frequency from first entry
    treeCh1->GetEntry(0);
    double dt_ns  = (t1[N-1] - t1[0]) / (N - 1);
    double fs_MHz = 1000.0 / dt_ns;

    std::cout << "\n  Entries     : " << nEntries    << std::endl;
    std::cout << "  dt          : " << dt_ns         << " ns"  << std::endl;
    std::cout << "  fs          : " << fs_MHz        << " MHz" << std::endl;
    std::cout << "  Nyquist     : " << fs_MHz / 2.0  << " MHz" << std::endl;

    // ── User parameters ───────────────────────────────────────
    printBox({"Analysis parameters"});

    std::cout << "\n  [1/6] SNR threshold (A_max / sigma_noise).\n"
              << "        Typical: 3-5 clean data, 5-8 noisy.\n";
    double snr_thr = promptDouble("SNR threshold", 1.0, 100.0);

    std::cout << "\n  [2/6] Signal window lower bound [ns] relative to laser trigger.\n"
              << "        Peak must appear AFTER this time.  (~45 ns)\n";
    double t_sig_lo = promptDouble("t_sig_lo [ns]", 0.0, 500.0);

    std::cout << "\n  [3/6] Signal window upper bound [ns] relative to laser trigger.\n";
    double t_sig_hi = promptDouble("t_sig_hi [ns]", t_sig_lo + 1.0, 1000.0);

    std::cout << "\n  [4/6] Heatmap window: ns BEFORE the CFD point to show.\n"
              << "        Shows the baseline just before the edge.  (5-15 ns)\n";
    double win_pre  = promptDouble("Window pre-CFD [ns]", 1.0, 100.0);

    std::cout << "\n  [5/6] Heatmap window: ns AFTER the CFD point to show.\n"
              << "        Should cover the full rise to the peak.  (10-40 ns)\n";
    double win_post = promptDouble("Window post-CFD [ns]", 1.0, 200.0);

    std::cout << "\n  [6/6] LP cutoff frequency [MHz] for the filtered heatmap.\n";
    double fc_MHz = promptDouble("LP cutoff [MHz]", 1.0, fs_MHz / 2.0 - 1.0);

    // ── Histogram setup ───────────────────────────────────────
    const int    T_BINS = 500;
    const int    A_BINS = 300;
    const double A_LO   = -5.0;
    const double A_HI   = 60.0;

    TH2D* hRaw = new TH2D("hRaw",
        "Rising Edge Heatmap (raw);t - t_{CFD} (ns);amplitude (mV)",
        T_BINS, -win_pre, win_post,
        A_BINS, A_LO, A_HI);

    TH2D* hFilt = new TH2D("hFilt",
        Form("Rising Edge Heatmap (LP %.0f MHz);t - t_{CFD} (ns);amplitude (mV)", fc_MHz),
        T_BINS, -win_pre, win_post,
        A_BINS, A_LO, A_HI);

    TH1D* hAmp = new TH1D("hAmp",
        "Signal Amplitude Spectrum;A_{max} (mV);counts",
        300, A_LO, A_HI);

    // ── Main loop ─────────────────────────────────────────────
    printBox({"Classifying and filling heatmaps..."});

    long nSig = 0, nNoise = 0, nNoLaser = 0, nNoCFD = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 1000 == 0) {
            std::cout << "\r  Progress: " << i << " / " << nEntries
                      << "  (" << (int)(100.0 * i / nEntries) << "%)"
                      << std::flush;
            gSystem->ProcessEvents();
        }

        treeCh1->GetEntry(i);
        treeLaser->GetEntry(i);

        // Baseline correction
        double sigma_ch1 = 0.0, sigma_dummy = 0.0;
        std::vector<double> ac1 = correctBaseline(t1, a1, N, 30.0, sigma_ch1);
        std::vector<double> acL = correctBaseline(tL, aL, N,  5.0, sigma_dummy);

        // Laser trigger
        double t_laser = laserTriggerTime(tL, acL, N, 10.0);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        // Signal classification
        auto   itMax    = std::max_element(ac1.begin(), ac1.end());
        double A_max    = *itMax;
        int    i_peak   = (int)std::distance(ac1.begin(), itMax);
        double t_pk_rel = t1[i_peak] - t_laser;
        double SNR      = A_max / sigma_ch1;

        bool is_signal = (SNR       >= snr_thr)
                      && (t_pk_rel  >= t_sig_lo)
                      && (t_pk_rel  <= t_sig_hi);

        if (!is_signal) { ++nNoise; continue; }
        ++nSig;
        hAmp->Fill(A_max);

        // CFD time (50%, backwards from peak)
        double t_cfd = cfdTime(t1, ac1, N, 0.5);
        if (t_cfd < -900.0) { ++nNoCFD; continue; }

        // Fill raw heatmap
        for (int j = 0; j < N; ++j) {
            double t_al = t1[j] - t_cfd;
            if (t_al < -win_pre || t_al > win_post) continue;
            hRaw->Fill(t_al, ac1[j]);
        }

        // LP filter and fill filtered heatmap
        std::vector<double> af = butterworthLP(ac1, fc_MHz, fs_MHz, 4);
        for (int j = 0; j < N; ++j) {
            double t_al = t1[j] - t_cfd;
            if (t_al < -win_pre || t_al > win_post) continue;
            hFilt->Fill(t_al, af[j]);
        }
    }

    std::cout << "\r  Progress: " << nEntries << " / " << nEntries
              << "  (100%)" << std::endl;

    double fracSig = 100.0 * nSig / (double)nEntries;

    printBox({
        "Results",
        "Total entries  : " + std::to_string(nEntries),
        "Signal         : " + std::to_string(nSig)
            + "  (" + std::to_string((int)fracSig) + "%)",
        "Noise (cut)    : " + std::to_string(nNoise),
        "No laser trig  : " + std::to_string(nNoLaser),
        "No CFD found   : " + std::to_string(nNoCFD),
        "SNR threshold  : " + std::to_string((int)snr_thr),
        "Peak window    : [" + std::to_string((int)t_sig_lo)
            + ", " + std::to_string((int)t_sig_hi) + "] ns"
    });

    // ── Canvas 1: Raw heatmap ─────────────────────────────────
    TCanvas* cRaw = new TCanvas("cRaw",
        "Rising Edge Heatmap -- Raw", 1000, 750);
    cRaw->SetLeftMargin(0.12);
    cRaw->SetRightMargin(0.14);
    cRaw->SetBottomMargin(0.12);
    cRaw->SetLogz();
    hRaw->SetContour(80);
    hRaw->GetXaxis()->SetTitleSize(0.045);
    hRaw->GetYaxis()->SetTitleSize(0.045);
    hRaw->GetYaxis()->SetTitleOffset(1.1);
    hRaw->Draw("COLZ");

    // Dashed vertical line at t=0 (CFD crossing point)
    TLine* lRaw = new TLine(0, A_LO, 0, A_HI);
    lRaw->SetLineColor(kRed+1); lRaw->SetLineStyle(2); lRaw->SetLineWidth(2);
    lRaw->Draw("SAME");

    TPaveText* ptRaw = new TPaveText(0.13, 0.76, 0.58, 0.90, "brNDC");
    ptRaw->SetBorderSize(1); ptRaw->SetFillColor(kWhite); ptRaw->SetFillStyle(1001);
    ptRaw->SetTextFont(42); ptRaw->SetTextSize(0.026); ptRaw->SetTextAlign(12);
    ptRaw->AddText(Form("Signal events: %ld  (%.1f%%)", nSig, fracSig));
    ptRaw->AddText(Form("SNR > %.1f  |  Peak in [%.0f, %.0f] ns (rel. laser)",
        snr_thr, t_sig_lo, t_sig_hi));
    ptRaw->AddText("Aligned on CFD 50%  |  Raw (baseline corrected)");
    ptRaw->Draw();
    cRaw->Update();

    // ── Canvas 2: Filtered heatmap ────────────────────────────
    TCanvas* cFilt = new TCanvas("cFilt",
        Form("Rising Edge Heatmap -- LP %.0f MHz", fc_MHz), 1000, 750);
    cFilt->SetLeftMargin(0.12);
    cFilt->SetRightMargin(0.14);
    cFilt->SetBottomMargin(0.12);
    cFilt->SetLogz();
    hFilt->SetContour(80);
    hFilt->GetXaxis()->SetTitleSize(0.045);
    hFilt->GetYaxis()->SetTitleSize(0.045);
    hFilt->GetYaxis()->SetTitleOffset(1.1);
    hFilt->Draw("COLZ");

    TLine* lFilt = new TLine(0, A_LO, 0, A_HI);
    lFilt->SetLineColor(kRed+1); lFilt->SetLineStyle(2); lFilt->SetLineWidth(2);
    lFilt->Draw("SAME");

    TPaveText* ptFilt = new TPaveText(0.13, 0.76, 0.58, 0.90, "brNDC");
    ptFilt->SetBorderSize(1); ptFilt->SetFillColor(kWhite); ptFilt->SetFillStyle(1001);
    ptFilt->SetTextFont(42); ptFilt->SetTextSize(0.026); ptFilt->SetTextAlign(12);
    ptFilt->AddText(Form("Signal events: %ld  (%.1f%%)", nSig, fracSig));
    ptFilt->AddText(Form("SNR > %.1f  |  Peak in [%.0f, %.0f] ns (rel. laser)",
        snr_thr, t_sig_lo, t_sig_hi));
    ptFilt->AddText(Form("Aligned on CFD 50%%  |  Butterworth LP %.0f MHz", fc_MHz));
    ptFilt->Draw();
    cFilt->Update();

    // ── Canvas 3: Amplitude spectrum ──────────────────────────
    TCanvas* cAmpC = new TCanvas("cAmpC",
        "Signal Amplitude Spectrum", 900, 580);
    cAmpC->SetLeftMargin(0.12); cAmpC->SetLogy();
    hAmp->SetLineColor(kTeal-4);
    hAmp->SetFillColorAlpha(kTeal-4, 0.35);
    hAmp->GetXaxis()->SetTitleSize(0.045);
    hAmp->GetYaxis()->SetTitleSize(0.045);
    hAmp->GetYaxis()->SetTitleOffset(1.1);
    hAmp->Draw("HIST");
    cAmpC->Update();

    // ── Save ──────────────────────────────────────────────────
    cRaw ->SaveAs("reh_01_raw_heatmap.png");
    cFilt->SaveAs("reh_02_filtered_heatmap.png");
    cAmpC->SaveAs("reh_03_amplitude_spectrum.png");

    printBox({
        "Canvases saved:",
        "reh_01  Raw heatmap (CFD-aligned)",
        "reh_02  Filtered heatmap (LP " + std::to_string((int)fc_MHz) + " MHz)",
        "reh_03  Signal amplitude spectrum"
    });

    file->Close();

    // Keep canvases open until the user presses Enter.
    // ProcessEvents() keeps the GUI alive without crashing.
    std::cout << "\n  Press [Enter] to exit and close all canvases..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();

    std::cout << "\n  Done.\n" << std::endl;
}
