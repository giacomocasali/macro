// signal_classifier.cpp
//
// Classifies every waveform in data.vbias_{53}.root (TTree "ch1") as
// "signal" or "noise" using an adaptive per-waveform SNR criterion
// combined with a temporal-window constraint on the peak position.
//
// Classification criteria (BOTH must be satisfied):
//   1.  SNR = A_max / sigma_noise  >=  snr_threshold
//   2.  t_peak  in  [t_sig_lo, t_sig_hi]
//
// Canvases produced:
//   1. SNR distribution          (signal | noise, stacked, log-y)
//   2. A_max distribution        (signal | noise, stacked, log-y)
//   3. A_max vs SNR 2D map       (COLZ, log-z)
//   4. Best-SNR signal vs random noise waveform
//   5. Best-SNR signal + LP filter applied
//   6. Grid of N_SHOW random signal waveforms with LP filter applied
//
// Usage:
//   root -l 'signal_classifier.cpp+'
//   signal_classifier()

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TLine.h>
#include <TRandom3.h>
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

static int promptInt(const std::string& msg, int lo, int hi)
{
    int v = lo - 1;
    while (v < lo || v > hi) {
        std::cout << "  " << msg << " [" << lo << " .. " << hi << "]: ";
        std::cin >> v;
        if (std::cin.fail() || v < lo || v > hi) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "  [ERROR] Out of range or invalid input.\n";
            v = lo - 1;
        }
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return v;
}

// ═══════════════════════════════════════════════════════════════
// WAVEFORM METRICS
// ═══════════════════════════════════════════════════════════════

struct WaveformMetrics {
    double   sigma_noise;
    double   A_max;
    double   t_peak;
    double   SNR;
    bool     is_signal;
    Long64_t entry;
};

static std::vector<double> correctBaseline(
    const Double_t* time,
    const Double_t* amp,
    int N,
    double pre_end,
    double& sigma_noise_out)
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

    double mean_pre = 0.0;
    for (double v : pre) mean_pre += v;
    mean_pre /= (double)pre.size();
    double var = 0.0;
    for (double v : pre) var += (v - mean_pre) * (v - mean_pre);
    sigma_noise_out = (pre.size() > 1)
        ? std::sqrt(var / (double)pre.size()) : 1e-9;
    if (sigma_noise_out < 1e-9) sigma_noise_out = 1e-9;

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
// CLEANUP
// ═══════════════════════════════════════════════════════════════

static void cleanupTempFiles()
{
    std::cout << "\n[INFO] Cleaning up ACLiC temporary files..." << std::endl;
    for (const auto& f : {
            "signal_classifier_cpp.so",
            "signal_classifier_cpp.d",
            "signal_classifier_cpp_ACLiC_dict_rdict.pcm"})
        gSystem->Exec((std::string("rm -f ") + f).c_str());
}

// ═══════════════════════════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════════════════════════

void signal_classifier()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    std::cout << "\n[INFO] Working directory: "
              << gSystem->WorkingDirectory() << std::endl;

    printBox({"Signal Classifier  --  Adaptive SNR + Temporal Window"});

    // ── Open file ─────────────────────────────────────────────
    TFile* file = TFile::Open("../../data/data.vbias_{53}.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Cannot open data.vbias_{53}.root" << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get("ch1");
    if (!tree) {
        std::cerr << "[ERROR] TTree 'ch1' not found." << std::endl;
        file->Close(); return;
    }
    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[ERROR] TTree 'ch1' is empty." << std::endl;
        file->Close(); return;
    }

    const int N = 1024;
    Double_t t_arr[N], a_arr[N];
    tree->SetBranchAddress("time",      t_arr);
    tree->SetBranchAddress("amplitude", a_arr);

    tree->GetEntry(0);
    double dt_ns  = (t_arr[N-1] - t_arr[0]) / (N - 1);
    double fs_MHz = 1000.0 / dt_ns;

    std::cout << "\n  Entries     : " << nEntries    << std::endl;
    std::cout << "  Samples/wfm : " << N             << std::endl;
    std::cout << "  dt          : " << dt_ns         << " ns" << std::endl;
    std::cout << "  fs          : " << fs_MHz        << " MHz" << std::endl;
    std::cout << "  Nyquist     : " << fs_MHz / 2.0  << " MHz" << std::endl;

    // ── User parameters ───────────────────────────────────────
    printBox({"Classification + display parameters"});

    std::cout << "\n  [1/5] SNR threshold  (A_max / sigma_noise).\n"
              << "        Typical: 3-5 clean data, 5-8 noisy data.\n";
    double snr_thr = promptDouble("SNR threshold", 1.0, 100.0);

    std::cout << "\n  [2/5] Signal window lower bound t_sig_lo [ns].\n"
              << "        Peak must appear AFTER this time.  (~35 ns)\n";
    double t_sig_lo = promptDouble("t_sig_lo [ns]", 0.0, 500.0);

    std::cout << "\n  [3/5] Signal window upper bound t_sig_hi [ns].\n";
    double t_sig_hi = promptDouble("t_sig_hi [ns]", t_sig_lo + 1.0, 1000.0);

    std::cout << "\n  [4/5] Low-pass cutoff [MHz] for filtered canvases.\n";
    double fc_MHz = promptDouble("Cutoff [MHz]", 1.0, fs_MHz / 2.0 - 1.0);

    std::cout << "\n  [5/5] How many random signal waveforms in the grid?\n"
              << "        (1-25, suggested 9 for 3x3 or 16 for 4x4)\n";
    int N_SHOW = promptInt("Number of waveforms", 1, 25);

    // ── Classify all waveforms ────────────────────────────────
    printBox({"Classifying " + std::to_string(nEntries) + " waveforms..."});

    std::vector<WaveformMetrics> metrics;
    metrics.reserve(nEntries);

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 1000 == 0) {
            std::cout << "\r  Progress: " << i << " / " << nEntries
                      << "  (" << (int)(100.0 * i / nEntries) << "%)"
                      << std::flush;
            gSystem->ProcessEvents();
        }
        tree->GetEntry(i);

        double sigma = 0.0;
        std::vector<double> ac = correctBaseline(t_arr, a_arr, N, 30.0, sigma);

        auto   itMax = std::max_element(ac.begin(), ac.end());
        double A_max = *itMax;
        int    iMax  = (int)std::distance(ac.begin(), itMax);
        double t_pk  = t_arr[iMax];
        double SNR   = A_max / sigma;

        bool is_sig = (SNR >= snr_thr)
                   && (t_pk >= t_sig_lo)
                   && (t_pk <= t_sig_hi);

        metrics.push_back({sigma, A_max, t_pk, SNR, is_sig, i});
    }
    std::cout << "\r  Progress: " << nEntries << " / " << nEntries
              << "  (100%)" << std::endl;

    // ── Summary ───────────────────────────────────────────────
    long nSig = 0, nNoise = 0;
    for (const auto& m : metrics) { if (m.is_signal) ++nSig; else ++nNoise; }
    double fracSig = 100.0 * nSig / (double)nEntries;

    // Signal entries sorted by SNR descending
    std::vector<size_t> sig_sorted;
    sig_sorted.reserve(nSig);
    for (size_t i = 0; i < metrics.size(); ++i)
        if (metrics[i].is_signal) sig_sorted.push_back(i);
    std::sort(sig_sorted.begin(), sig_sorted.end(),
        [&](size_t a, size_t b){ return metrics[a].SNR > metrics[b].SNR; });

    Long64_t best_entry = sig_sorted.empty() ? -1 : metrics[sig_sorted[0]].entry;
    double   best_snr   = sig_sorted.empty() ? 0.0 : metrics[sig_sorted[0]].SNR;

    // Random noise entry
    std::vector<Long64_t> noise_pool;
    for (const auto& m : metrics) if (!m.is_signal) noise_pool.push_back(m.entry);
    TRandom3 rng(0);
    Long64_t noise_entry = noise_pool.empty() ? -1
        : noise_pool[(Long64_t)(rng.Uniform(0, noise_pool.size()))];
    double noise_snr = 0.0;
    for (const auto& m : metrics)
        if (m.entry == noise_entry) { noise_snr = m.SNR; break; }

    // N_SHOW random signal entries for grid (shuffled from full signal pool)
    std::vector<Long64_t> grid_entries;
    {
        std::vector<Long64_t> pool;
        pool.reserve(nSig);
        for (const auto& idx : sig_sorted) pool.push_back(metrics[idx].entry);
        for (int i = (int)pool.size() - 1; i > 0; --i) {
            int j = (int)(rng.Uniform(0, i + 1));
            std::swap(pool[i], pool[j]);
        }
        int take = std::min((int)pool.size(), N_SHOW);
        for (int i = 0; i < take; ++i) grid_entries.push_back(pool[i]);
    }

    printBox({
        "Results",
        "Total    : " + std::to_string(nEntries),
        "Signal   : " + std::to_string(nSig)
            + "  (" + std::to_string((int)fracSig) + "%)",
        "Noise    : " + std::to_string(nNoise)
            + "  (" + std::to_string((int)(100 - fracSig)) + "%)",
        "SNR thr  : " + std::to_string((int)snr_thr),
        "Window   : [" + std::to_string((int)t_sig_lo)
            + ", " + std::to_string((int)t_sig_hi) + "] ns",
        "Best SNR : entry #" + std::to_string(best_entry)
            + "  (SNR=" + std::to_string((int)best_snr) + ")"
    });

    // ── Helper: reload + baseline-correct one entry ───────────
    auto reloadWaveform = [&](Long64_t idx,
                               std::vector<double>& tvec,
                               std::vector<double>& avec)
    {
        tree->GetEntry(idx);
        double dummy = 0.0;
        tvec.assign(t_arr, t_arr + N);
        avec = correctBaseline(t_arr, a_arr, N, 30.0, dummy);
    };

    // ── Histogram ranges (99th percentile) ───────────────────
    double snr_max = 0.0, amax_max = 0.0;
    {
        std::vector<double> sv, av;
        sv.reserve(nEntries); av.reserve(nEntries);
        for (const auto& m : metrics) { sv.push_back(m.SNR); av.push_back(m.A_max); }
        std::sort(sv.begin(), sv.end()); std::sort(av.begin(), av.end());
        size_t p99 = (size_t)(0.99 * sv.size());
        snr_max  = sv[std::min(p99, sv.size() - 1)];
        amax_max = av[std::min(p99, av.size() - 1)];
    }

    const int NBINS = 150;
    TH1D* hSNR_sig   = new TH1D("hSNR_sig",   "", NBINS, 0, snr_max  * 1.05);
    TH1D* hSNR_noise = new TH1D("hSNR_noise", "", NBINS, 0, snr_max  * 1.05);
    TH1D* hAmp_sig   = new TH1D("hAmp_sig",   "", NBINS, 0, amax_max * 1.05);
    TH1D* hAmp_noise = new TH1D("hAmp_noise", "", NBINS, 0, amax_max * 1.05);
    TH2D* h2D = new TH2D("h2D", "", 200, 0, snr_max*1.05, 200, 0, amax_max*1.05);

    for (const auto& m : metrics) {
        if (m.is_signal) { hSNR_sig->Fill(m.SNR); hAmp_sig->Fill(m.A_max); }
        else             { hSNR_noise->Fill(m.SNR); hAmp_noise->Fill(m.A_max); }
        h2D->Fill(m.SNR, m.A_max);
    }

    // ── CANVAS 1: SNR distribution ────────────────────────────
    TCanvas* cSNR = new TCanvas("cSNR", "SNR Distribution", 950, 580);
    cSNR->SetLeftMargin(0.12); cSNR->SetLogy();
    hSNR_sig->SetTitle(Form(
        "SNR Distribution  (threshold = %.1f);A_{max}/#sigma_{noise};counts", snr_thr));
    hSNR_sig->SetLineColor(kAzure+1);   hSNR_sig->SetFillColorAlpha(kAzure+1, 0.35);
    hSNR_noise->SetLineColor(kGray+1);  hSNR_noise->SetFillColorAlpha(kGray+1, 0.45);
    (hSNR_noise->GetMaximum() > hSNR_sig->GetMaximum())
        ? (hSNR_noise->Draw("HIST"), hSNR_sig->Draw("HIST SAME"))
        : (hSNR_sig->Draw("HIST"),   hSNR_noise->Draw("HIST SAME"));
    TLine* lSNR = new TLine(snr_thr, cSNR->GetUymin(),
                             snr_thr, hSNR_sig->GetMaximum() * 2.0);
    lSNR->SetLineColor(kRed+1); lSNR->SetLineStyle(2); lSNR->SetLineWidth(2);
    lSNR->Draw("SAME");
    TLegend* legSNR = new TLegend(0.55, 0.72, 0.88, 0.88);
    legSNR->SetBorderSize(1); legSNR->SetFillColor(kWhite);
    legSNR->AddEntry(hSNR_sig,   Form("Signal  (%ld)",  nSig),              "f");
    legSNR->AddEntry(hSNR_noise, Form("Noise   (%ld)", nNoise),             "f");
    legSNR->AddEntry(lSNR,       Form("Threshold = %.1f", snr_thr),         "l");
    legSNR->Draw();
    cSNR->Update();

    // ── CANVAS 2: Amplitude distribution ─────────────────────
    TCanvas* cAmp = new TCanvas("cAmp", "Amplitude Distribution", 950, 580);
    cAmp->SetLeftMargin(0.12); cAmp->SetLogy();
    hAmp_sig->SetTitle("A_{max} Distribution;A_{max} (mV);counts");
    hAmp_sig->SetLineColor(kTeal-4);    hAmp_sig->SetFillColorAlpha(kTeal-4, 0.35);
    hAmp_noise->SetLineColor(kGray+1);  hAmp_noise->SetFillColorAlpha(kGray+1, 0.45);
    (hAmp_noise->GetMaximum() > hAmp_sig->GetMaximum())
        ? (hAmp_noise->Draw("HIST"), hAmp_sig->Draw("HIST SAME"))
        : (hAmp_sig->Draw("HIST"),   hAmp_noise->Draw("HIST SAME"));
    TLegend* legAmp = new TLegend(0.55, 0.72, 0.88, 0.88);
    legAmp->SetBorderSize(1); legAmp->SetFillColor(kWhite);
    legAmp->AddEntry(hAmp_sig,   Form("Signal  (%ld)",  nSig),  "f");
    legAmp->AddEntry(hAmp_noise, Form("Noise   (%ld)", nNoise), "f");
    legAmp->Draw();
    cAmp->Update();

    // ── CANVAS 3: 2D scatter A_max vs SNR ────────────────────
    TCanvas* c2D = new TCanvas("c2D", "A_max vs SNR 2D", 950, 700);
    c2D->SetLeftMargin(0.13); c2D->SetRightMargin(0.13); c2D->SetLogz();
    h2D->SetTitle(Form(
        "A_{max} vs SNR  |  threshold = %.1f"
        ";SNR  (A_{max}/#sigma_{noise});A_{max} (mV)", snr_thr));
    h2D->SetContour(80); h2D->GetYaxis()->SetTitleOffset(1.3); h2D->Draw("COLZ");
    TLine* lv = new TLine(snr_thr, 0, snr_thr, amax_max * 1.05);
    lv->SetLineColor(kRed+1); lv->SetLineStyle(2); lv->SetLineWidth(2); lv->Draw("SAME");
    TPaveText* pthr = new TPaveText(0.14, 0.80, 0.52, 0.90, "brNDC");
    pthr->SetBorderSize(1); pthr->SetFillColor(kWhite); pthr->SetFillStyle(1001);
    pthr->SetTextFont(42); pthr->SetTextSize(0.028); pthr->SetTextAlign(12);
    pthr->AddText(Form("SNR cut = %.1f  |  Peak in [%.0f, %.0f] ns",
        snr_thr, t_sig_lo, t_sig_hi));
    pthr->AddText(Form("Signal: %ld (%.1f%%)  |  Noise: %ld (%.1f%%)",
        nSig, fracSig, nNoise, 100.0 - fracSig));
    pthr->Draw(); c2D->Update();

    // ── CANVAS 4: best-SNR signal vs random noise ─────────────
    TCanvas* cEx = new TCanvas("cEx", "Signal vs Noise -- Example", 1000, 600);
    cEx->SetLeftMargin(0.12); cEx->SetRightMargin(0.05); cEx->SetBottomMargin(0.12);
    std::vector<double> t_best, a_best, t_noi, a_noi;
    if (best_entry  >= 0) reloadWaveform(best_entry,  t_best, a_best);
    if (noise_entry >= 0) reloadWaveform(noise_entry, t_noi,  a_noi);
    TGraph* grBest = t_best.empty() ? nullptr
        : new TGraph(N, t_best.data(), a_best.data());
    TGraph* grNoi  = t_noi.empty()  ? nullptr
        : new TGraph(N, t_noi.data(),  a_noi.data());
    if (grBest) { grBest->SetLineColor(kAzure+1); grBest->SetLineWidth(2); }
    if (grNoi)  { grNoi->SetLineColor(kGray+1);   grNoi->SetLineWidth(1);  }
    if (grBest) {
        grBest->SetTitle(";time (ns);amplitude (mV)");
        grBest->GetXaxis()->SetTitleSize(0.045);
        grBest->GetYaxis()->SetTitleSize(0.045);
        grBest->GetYaxis()->SetTitleOffset(1.1);
        grBest->Draw("AL");
    }
    if (grNoi) grNoi->Draw("L SAME");
    TLegend* legEx = new TLegend(0.55, 0.78, 0.93, 0.90);
    legEx->SetBorderSize(1); legEx->SetFillColor(kWhite); legEx->SetTextSize(0.033);
    if (grBest) legEx->AddEntry(grBest,
        Form("Best signal  #%lld  (SNR=%.1f)", best_entry, best_snr), "l");
    if (grNoi)  legEx->AddEntry(grNoi,
        Form("Noise        #%lld  (SNR=%.1f)", noise_entry, noise_snr), "l");
    legEx->Draw();
    TPaveText* ptEx = new TPaveText(0.13, 0.78, 0.50, 0.90, "brNDC");
    ptEx->SetBorderSize(1); ptEx->SetFillColor(kWhite); ptEx->SetFillStyle(1001);
    ptEx->SetTextFont(42); ptEx->SetTextSize(0.030); ptEx->SetTextAlign(12);
    ptEx->AddText(Form("Best signal entry #%lld:  SNR = %.1f", best_entry, best_snr));
    ptEx->AddText(Form("Noise entry       #%lld:  SNR = %.1f", noise_entry, noise_snr));
    ptEx->Draw(); cEx->Update();

    // ── CANVAS 5: best-SNR signal + LP filter ────────────────
    TCanvas* cFilt = nullptr;
    if (!t_best.empty()) {
        std::vector<double> af = butterworthLP(a_best, fc_MHz, fs_MHz, 4);
        TGraph* gR = new TGraph(N, t_best.data(), a_best.data());
        TGraph* gF = new TGraph(N, t_best.data(), af.data());
        gR->SetLineColor(kGray+1);    gR->SetLineWidth(1);
        gF->SetLineColor(kAzure-4);   gF->SetLineWidth(2);
        cFilt = new TCanvas("cFilt",
            Form("Best signal #%lld  LP f_c=%.1f MHz", best_entry, fc_MHz),
            1000, 580);
        cFilt->SetLeftMargin(0.12); cFilt->SetRightMargin(0.05);
        cFilt->SetBottomMargin(0.12);
        gR->SetTitle(Form(
            "Best signal entry #%lld  |  Butterworth LP  f_{c} = %.1f MHz"
            "  (f_{s} = %.0f MHz);time (ns);amplitude (mV)",
            best_entry, fc_MHz, fs_MHz));
        gR->GetXaxis()->SetTitleSize(0.045);
        gR->GetYaxis()->SetTitleSize(0.045);
        gR->GetYaxis()->SetTitleOffset(1.1);
        gR->Draw("AL"); gF->Draw("L SAME");
        TLegend* lgF = new TLegend(0.55, 0.78, 0.93, 0.90);
        lgF->SetBorderSize(1); lgF->SetFillColor(kWhite); lgF->SetTextSize(0.033);
        lgF->AddEntry(gR, "Raw (baseline-corrected)", "l");
        lgF->AddEntry(gF, Form("Filtered  f_{c} = %.1f MHz", fc_MHz), "l");
        lgF->Draw();
        TPaveText* ptF = new TPaveText(0.13, 0.78, 0.52, 0.90, "brNDC");
        ptF->SetBorderSize(1); ptF->SetFillColor(kWhite); ptF->SetFillStyle(1001);
        ptF->SetTextFont(42); ptF->SetTextSize(0.030); ptF->SetTextAlign(12);
        ptF->AddText(Form("Entry #%lld  |  SNR = %.1f  |  Best signal",
            best_entry, best_snr));
        ptF->AddText(Form("f_{s} = %.0f MHz   f_{c} = %.1f MHz", fs_MHz, fc_MHz));
        ptF->AddText("4^{th}-order Butterworth LP");
        ptF->Draw(); cFilt->Update();
    }

    // ── CANVAS 6: grid of N_SHOW random signal waveforms ─────
    TCanvas* cGrid = nullptr;
    if (!grid_entries.empty()) {
        int nShow = (int)grid_entries.size();
        int nCols = (int)std::ceil(std::sqrt((double)nShow));
        int nRows = (int)std::ceil((double)nShow / nCols);
        int cW = std::min(1800, 350 * nCols);
        int cH = std::min(1200, 280 * nRows);

        cGrid = new TCanvas("cGrid",
            Form("%d Random Signal Waveforms  |  LP f_c=%.1f MHz", nShow, fc_MHz),
            cW, cH);
        cGrid->Divide(nCols, nRows, 0.001, 0.001);

        for (int k = 0; k < nShow; ++k) {
            cGrid->cd(k + 1);
            gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03);
            gPad->SetTopMargin(0.12);  gPad->SetBottomMargin(0.14);
            gPad->SetGridx(true);      gPad->SetGridy(true);

            Long64_t eidx = grid_entries[k];
            std::vector<double> tv, av;
            reloadWaveform(eidx, tv, av);
            std::vector<double> af = butterworthLP(av, fc_MHz, fs_MHz, 4);

            double snr_k = 0.0;
            for (const auto& m : metrics)
                if (m.entry == eidx) { snr_k = m.SNR; break; }

            TGraph* gR = new TGraph(N, tv.data(), av.data());
            TGraph* gF = new TGraph(N, tv.data(), af.data());
            gR->SetLineColor(kGray+1);  gR->SetLineWidth(1);
            gF->SetLineColor(kAzure-4); gF->SetLineWidth(2);

            gR->SetTitle(Form("Entry #%lld  SNR=%.1f;t (ns);amp (mV)", eidx, snr_k));
            gR->GetXaxis()->SetTitleSize(0.055);
            gR->GetYaxis()->SetTitleSize(0.055);
            gR->GetXaxis()->SetLabelSize(0.050);
            gR->GetYaxis()->SetLabelSize(0.050);
            gR->GetYaxis()->SetTitleOffset(0.9);
            gR->Draw("AL");
            gF->Draw("L SAME");

            // Legend only on first pad
            if (k == 0) {
                TLegend* lgk = new TLegend(0.35, 0.76, 0.96, 0.96);
                lgk->SetBorderSize(1); lgk->SetFillColor(kWhite);
                lgk->SetTextSize(0.060);
                lgk->AddEntry(gR, "Raw",                         "l");
                lgk->AddEntry(gF, Form("LP %.0f MHz", fc_MHz),   "l");
                lgk->Draw();
            }
        }
        cGrid->Update();
    }

    // ── Save all canvases ─────────────────────────────────────
    std::string wd = gSystem->WorkingDirectory();
    std::cout << "\n[INFO] Saving PNGs to: " << wd << "/" << std::endl;

    cSNR->SaveAs("sc_01_snr_distribution.png");
    cAmp->SaveAs("sc_02_amplitude_distribution.png");
    c2D ->SaveAs("sc_03_amax_vs_snr_2D.png");
    cEx ->SaveAs("sc_04_signal_vs_noise_example.png");
    if (cFilt) cFilt->SaveAs("sc_05_best_signal_filtered.png");
    if (cGrid) cGrid->SaveAs("sc_06_signal_grid_filtered.png");

    printBox({
        "Canvases saved to: " + wd,
        "sc_01  SNR distribution",
        "sc_02  Amplitude distribution",
        "sc_03  A_max vs SNR  (2D map)",
        "sc_04  Best signal vs random noise",
        "sc_05  Best signal + LP filter",
        "sc_06  Grid of " + std::to_string((int)grid_entries.size())
            + " random signals + LP filter"
    });

    file->Close();
    cleanupTempFiles();
    std::cout << "\n  Done.\n" << std::endl;
}
