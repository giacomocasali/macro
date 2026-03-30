// waveform_analysis.cpp
//
// Full two-channel waveform analysis with:
//   - Dual-channel readout (ch1, laser)
//   - Baseline correction (median of pre-signal window)
//   - 4th-order Butterworth low-pass filter
//   - Fixed threshold crossing detection
//   - Peak amplitude, pulse integral, rise time, fall time extraction
//   - Interactive cutoff & threshold selection
//   - CSV export of per-entry results
//   - Sampling-rate guard (expected 5 GS/s)
//
// Usage:
//   root -l 'waveform_analysis.cpp+'
//   then call:  waveform_analysis()

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>
#include <limits>
#include <iomanip>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "../header/ButterworthFilter.h"  // scipy-equivalent Butterworth IIR
#include <TMultiGraph.h>

// ─────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────
static const int    N_SAMPLES        = 1024;
static const double EXPECTED_FS_MHZ  = 5000.0;   // 5 GS/s
static const double FS_TOLERANCE     = 0.05;      // 5% tolerance
static const double PRE_SIGNAL_NS    = 30.0;      // baseline window [ns]

// ─────────────────────────────────────────────────────────────
// Utility: repeat string n times
// ─────────────────────────────────────────────────────────────
static std::string repeatStr(const std::string& s, size_t n)
{
    std::string out;
    out.reserve(s.size() * n);
    for (size_t i = 0; i < n; ++i) out += s;
    return out;
}

// ─────────────────────────────────────────────────────────────
// Pretty ASCII box
// ─────────────────────────────────────────────────────────────
static void printBox(const std::vector<std::string>& lines)
{
    size_t maxLen = 0;
    for (const auto& l : lines)
        if (l.size() > maxLen) maxLen = l.size();

    const size_t width = maxLen + 4;
    std::cout << "\n  +" << repeatStr("-", width) << "+" << std::endl;
    for (const auto& l : lines) {
        size_t pad = width - l.size() - 2;
        std::cout << "  |  " << l << std::string(pad, ' ') << "|" << std::endl;
    }
    std::cout << "  +" << repeatStr("-", width) << "+" << std::endl;
}

// ─────────────────────────────────────────────────────────────
// Cleanup ACLiC temp files
// ─────────────────────────────────────────────────────────────
static void cleanupTempFiles()
{
    std::cout << "\n[INFO] Cleaning up temporary files..." << std::endl;
    const std::vector<std::string> files = {
        "waveform_analysis_cpp.so",
        "waveform_analysis_cpp.d",
        "waveform_analysis_cpp_ACLiC_dict_rdict.pcm"
    };
    for (const auto& f : files) gSystem->Exec(("rm -f " + f).c_str());
    std::cout << "[INFO] Cleanup done." << std::endl;
}

// ─────────────────────────────────────────────────────────────
// Baseline correction: subtract median of pre-signal samples
// ─────────────────────────────────────────────────────────────
static std::vector<double> correctBaseline(
    const std::vector<double>& time,
    const std::vector<double>& amp,
    double pre_signal_end = PRE_SIGNAL_NS)
{
    std::vector<double> pre;
    for (size_t i = 0; i < time.size(); ++i)
        if (time[i] <= pre_signal_end) pre.push_back(amp[i]);

    if (pre.empty()) {
        size_t n = std::max((size_t)10, time.size() / 10);
        for (size_t i = 0; i < n && i < amp.size(); ++i) pre.push_back(amp[i]);
    }

    std::vector<double> tmp = pre;
    size_t mid = tmp.size() / 2;
    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
    double offset = tmp[mid];

    std::vector<double> out = amp;
    for (auto& a : out) a -= offset;
    return out;
}

// ─────────────────────────────────────────────────────────────
// Baseline RMS noise (std-dev of pre-signal window)
// ─────────────────────────────────────────────────────────────
static double baselineRMS(
    const std::vector<double>& time,
    const std::vector<double>& amp_corr,
    double pre_signal_end = PRE_SIGNAL_NS)
{
    std::vector<double> pre;
    for (size_t i = 0; i < time.size(); ++i)
        if (time[i] <= pre_signal_end) pre.push_back(amp_corr[i]);
    if (pre.size() < 2) return 0.0;

    double mean = std::accumulate(pre.begin(), pre.end(), 0.0) / pre.size();
    double sq   = 0.0;
    for (double v : pre) sq += (v - mean) * (v - mean);
    return std::sqrt(sq / (pre.size() - 1));
}

// ─────────────────────────────────────────────────────────────
// 4th-order Butterworth low-pass filter (IIR, biquad cascade)
// ─────────────────────────────────────────────────────────────
// butterworthLowPass is provided by ButterworthFilter.h
// To use zero-phase (sosfiltfilt): butterworthLowPassZP(x, fc, fs, 4)
// To use fast pre-computed 500MHz:  butterworthLP_500MHz(x)

// ─────────────────────────────────────────────────────────────
// Pulse metrics extracted from a filtered waveform
// ─────────────────────────────────────────────────────────────
struct PulseMetrics {
    bool   found        = false;
    double peak_amp_mV  = 0.0;   // amplitude at peak [mV]
    double peak_time_ns = 0.0;   // time of peak [ns]
    double integral_mVns= 0.0;   // trapezoid integral over pulse window [mV·ns]
    double rise_time_ns = 0.0;   // 10%-90% rise time [ns]
    double fall_time_ns = 0.0;   // 90%-10% fall time [ns]
    double thresh_cross_ns = 0.0;// first threshold crossing [ns]
    int    n_crossings  = 0;     // number of threshold crossings (signal activity)
    double noise_rms_mV = 0.0;   // baseline RMS [mV]
    std::string channel;
};

static PulseMetrics extractMetrics(
    const std::vector<double>& time,
    const std::vector<double>& amp_filt,
    const std::vector<double>& amp_raw_corr,
    double threshold_mV,
    const std::string& chName)
{
    PulseMetrics m;
    m.channel      = chName;
    m.noise_rms_mV = baselineRMS(time, amp_raw_corr);

    int N = (int)time.size();
    if (N < 4) return m;

    double dt = time[1] - time[0];

    // ── Peak ──────────────────────────────────────────────────
    // Support both positive and negative pulses
    auto itMax = std::max_element(amp_filt.begin(), amp_filt.end());
    auto itMin = std::min_element(amp_filt.begin(), amp_filt.end());
    double peakPos = *itMax;
    double peakNeg = *itMin;
    bool negPulse  = std::abs(peakNeg) > std::abs(peakPos);
    int iPeak      = negPulse
                     ? (int)(itMin - amp_filt.begin())
                     : (int)(itMax - amp_filt.begin());
    m.peak_amp_mV  = amp_filt[iPeak];
    m.peak_time_ns = time[iPeak];

    // ── Threshold crossings ───────────────────────────────────
    double thr = negPulse ? -std::abs(threshold_mV) : std::abs(threshold_mV);
    bool above  = false;
    bool firstFound = false;
    for (int i = 1; i < N; ++i) {
        bool now = negPulse ? (amp_filt[i] < thr) : (amp_filt[i] > thr);
        if (now && !above) {
            ++m.n_crossings;
            above = true;
            if (!firstFound) {
                m.thresh_cross_ns = time[i];
                firstFound = true;
                m.found = true;
            }
        }
        if (!now) above = false;
    }
    if (!m.found) return m;

    // ── Pulse window: from first crossing to peak + 3× rise window ──
    int iStart = iPeak;
    for (int i = iPeak; i >= 0; --i) {
        bool outside = negPulse ? (amp_filt[i] > thr) : (amp_filt[i] < thr);
        if (outside) { iStart = i; break; }
    }
    int iEnd = std::min(N - 1, iPeak + (iPeak - iStart) * 3);

    // ── Integral (trapezoid) over pulse window ─────────────────
    double integ = 0.0;
    for (int i = iStart; i < iEnd; ++i)
        integ += 0.5 * (amp_filt[i] + amp_filt[i + 1]) * dt;
    m.integral_mVns = integ;

    // ── Rise time (10%–90% of peak) ───────────────────────────
    double p10 = 0.1 * m.peak_amp_mV;
    double p90 = 0.9 * m.peak_amp_mV;
    double t10 = -1.0, t90 = -1.0;
    for (int i = iStart; i < iPeak; ++i) {
        if (t10 < 0 && (negPulse ? amp_filt[i] < p10 : amp_filt[i] > p10)) t10 = time[i];
        if (t90 < 0 && (negPulse ? amp_filt[i] < p90 : amp_filt[i] > p90)) t90 = time[i];
    }
    if (t10 >= 0 && t90 >= 0) m.rise_time_ns = std::abs(t90 - t10);

    // ── Fall time (90%–10% after peak) ────────────────────────
    double f90 = -1.0, f10 = -1.0;
    for (int i = iPeak; i < iEnd; ++i) {
        if (f90 < 0 && (negPulse ? amp_filt[i] > p90 : amp_filt[i] < p90)) f90 = time[i];
        if (f10 < 0 && (negPulse ? amp_filt[i] > p10 : amp_filt[i] < p10)) f10 = time[i];
    }
    if (f90 >= 0 && f10 >= 0) m.fall_time_ns = std::abs(f10 - f90);

    return m;
}

// ─────────────────────────────────────────────────────────────
// CSV result row
// ─────────────────────────────────────────────────────────────
struct ResultRow {
    Long64_t entry;
    int      run;
    double   cutoff_MHz;
    double   threshold_mV;
    PulseMetrics ch1, laser;
};

static void writeCSVHeader(std::ofstream& f)
{
    f << "entry,run,cutoff_MHz,threshold_mV,"
      << "ch1_found,ch1_peak_mV,ch1_peak_time_ns,ch1_integral_mVns,"
      << "ch1_rise_ns,ch1_fall_ns,ch1_thresh_cross_ns,ch1_n_crossings,ch1_noise_rms_mV,"
      << "laser_found,laser_peak_mV,laser_peak_time_ns,laser_integral_mVns,"
      << "laser_rise_ns,laser_fall_ns,laser_thresh_cross_ns,laser_n_crossings,laser_noise_rms_mV\n";
}

static void writeCSVRow(std::ofstream& f, const ResultRow& r)
{
    auto b = [](bool v){ return v ? "1" : "0"; };
    f << r.entry << "," << r.run << ","
      << std::fixed << std::setprecision(3)
      << r.cutoff_MHz    << "," << r.threshold_mV   << ","
      // ch1
      << b(r.ch1.found)  << "," << r.ch1.peak_amp_mV    << "," << r.ch1.peak_time_ns  << ","
      << r.ch1.integral_mVns  << "," << r.ch1.rise_time_ns  << "," << r.ch1.fall_time_ns  << ","
      << r.ch1.thresh_cross_ns << "," << r.ch1.n_crossings  << "," << r.ch1.noise_rms_mV  << ","
      // laser
      << b(r.laser.found)  << "," << r.laser.peak_amp_mV    << "," << r.laser.peak_time_ns  << ","
      << r.laser.integral_mVns  << "," << r.laser.rise_time_ns  << "," << r.laser.fall_time_ns  << ","
      << r.laser.thresh_cross_ns << "," << r.laser.n_crossings  << "," << r.laser.noise_rms_mV  << "\n";
}

// ─────────────────────────────────────────────────────────────
// stdin helpers
// ─────────────────────────────────────────────────────────────
static char askYesNo(const std::string& prompt)
{
    char choice = 0;
    while (choice != 'y' && choice != 'n') {
        std::cout << prompt;
        std::string line;
        std::getline(std::cin, line);
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos != std::string::npos)
            choice = std::tolower((unsigned char)line[pos]);
        if (choice != 'y' && choice != 'n')
            std::cout << "  [ERROR] Please enter 'y' or 'n'." << std::endl;
    }
    return choice;
}

static double promptPositiveDouble(const std::string& msg)
{
    double val = -1.0;
    while (val <= 0.0) {
        std::cout << msg;
        std::cin >> val;
        if (std::cin.fail() || val <= 0.0) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "  [ERROR] Please enter a positive number." << std::endl;
            val = -1.0;
        }
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return val;
}

// ─────────────────────────────────────────────────────────────
// Load one entry from a TTree (static branch arrays)
// ─────────────────────────────────────────────────────────────
static bool loadEntry(
    TTree* tree,
    Long64_t idx,
    std::vector<double>& time,
    std::vector<double>& amp)
{
    static Double_t t_arr[1024], a_arr[1024];
    tree->SetBranchAddress("time",      t_arr);
    tree->SetBranchAddress("amplitude", a_arr);
    tree->GetEntry(idx);
    time.assign(t_arr, t_arr + N_SAMPLES);
    amp .assign(a_arr, a_arr + N_SAMPLES);
    return true;
}

// ─────────────────────────────────────────────────────────────
// Sampling-rate check
// ─────────────────────────────────────────────────────────────
static void checkSamplingRate(const std::vector<double>& time, double& fs_MHz)
{
    double dt_ns = (time.back() - time.front()) / (N_SAMPLES - 1);
    fs_MHz = 1000.0 / dt_ns;
    double ratio = fs_MHz / EXPECTED_FS_MHZ;
    if (std::abs(ratio - 1.0) > FS_TOLERANCE)
        std::cerr << "[WARNING] Measured fs = " << fs_MHz
                  << " MHz, expected " << EXPECTED_FS_MHZ
                  << " MHz. Check digitizer settings." << std::endl;
}

// ─────────────────────────────────────────────────────────────
// PHASE 1 – Preview loop (dual-channel overlay)
// ─────────────────────────────────────────────────────────────
static bool previewLoop(
    TTree* tree1,  TTree* tree2,
    Long64_t nEntries,
    TRandom3& rng,
    std::vector<double>& out_time,
    std::vector<double>& out_amp1_corr,
    std::vector<double>& out_amp2_corr,
    Long64_t& out_entry,
    double& out_fs_MHz)
{
    TCanvas* cPrev = nullptr;
    int attempt = 0;

    while (true) {
        ++attempt;

        std::vector<double> time, amp1, amp2;
        Long64_t idx = (Long64_t)(rng.Uniform(0, (double)nEntries));

        loadEntry(tree1, idx, time, amp1);
        std::vector<double> t2_dummy, amp2_raw;
        loadEntry(tree2, idx, t2_dummy, amp2_raw);

        std::vector<double> amp1_corr = correctBaseline(time, amp1);
        std::vector<double> amp2_corr = correctBaseline(time, amp2_raw);

        double fs_MHz = 0.0;
        checkSamplingRate(time, fs_MHz);

        if (cPrev) { cPrev->Close(); delete cPrev; cPrev = nullptr; }

        TGraph* gr1 = new TGraph(N_SAMPLES, time.data(), amp1_corr.data());
        TGraph* gr2 = new TGraph(N_SAMPLES, time.data(), amp2_corr.data());
        gr1->SetLineColor(kTeal  - 5);  gr1->SetLineWidth(2);
        gr2->SetLineColor(kOrange + 7); gr2->SetLineWidth(2);

        cPrev = new TCanvas("waveform_preview",
            Form("PREVIEW -- Entry #%lld (attempt %d)", idx, attempt),
            1100, 550);
        cPrev->SetLeftMargin(0.11); cPrev->SetRightMargin(0.05);
        cPrev->SetBottomMargin(0.12);
        gStyle->SetOptStat(0);
        gStyle->SetPadGridX(true); gStyle->SetPadGridY(true);

        TMultiGraph* mg = new TMultiGraph();
        mg->Add(gr1, "L");
        mg->Add(gr2, "L");
        mg->SetTitle(Form(
            "PREVIEW -- Entry #%lld  (attempt %d)  |  f_{s} = %.0f MHz"
            ";time (ns);amplitude (mV)", idx, attempt, fs_MHz));
        mg->Draw("A");
        mg->GetXaxis()->SetTitleSize(0.045);
        mg->GetYaxis()->SetTitleSize(0.045);
        mg->GetYaxis()->SetTitleOffset(1.0);

        TLegend* leg = new TLegend(0.55, 0.80, 0.93, 0.90);
        leg->SetBorderSize(1); leg->SetFillColor(kWhite);
        leg->SetTextSize(0.033);
        leg->AddEntry(gr1, "ch1 (baseline-corrected)", "l");
        leg->AddEntry(gr2, "laser (baseline-corrected)", "l");
        leg->Draw();

        TPaveText* info = new TPaveText(0.12, 0.78, 0.50, 0.90, "brNDC");
        info->SetBorderSize(1); info->SetFillColor(kWhite); info->SetFillStyle(1001);
        info->SetTextFont(42);  info->SetTextSize(0.030); info->SetTextAlign(12);
        info->AddText(Form("Entry:  %lld / %lld", idx, nEntries));
        info->AddText(Form("f_{s} = %.0f MHz  (expected %.0f MHz)", fs_MHz, EXPECTED_FS_MHZ));
        info->AddText(Form("Attempt: %d", attempt));
        info->Draw();
        cPrev->Update();

        printBox({
            std::string(Form("PREVIEW  --  Entry #%lld  (attempt %d)", idx, attempt)),
            "Is this waveform good for analysis?"
        });

        if (askYesNo("  Your choice [y/n]: ") == 'y') {
            out_time      = time;
            out_amp1_corr = amp1_corr;
            out_amp2_corr = amp2_corr;
            out_entry     = idx;
            out_fs_MHz    = fs_MHz;
            return true;
        }
    }
}

// ─────────────────────────────────────────────────────────────
// Draw one analysis canvas (two sub-pads, one per channel)
// ─────────────────────────────────────────────────────────────
static TCanvas* drawAnalysisCanvas(
    const std::vector<double>& time,
    const std::vector<double>& amp1_corr,
    const std::vector<double>& amp1_filt,
    const std::vector<double>& amp2_corr,
    const std::vector<double>& amp2_filt,
    const PulseMetrics& m1,
    const PulseMetrics& m2,
    Long64_t entry, Long64_t nEntries,
    double cutoff_MHz, double threshold_mV, double fs_MHz,
    int runIndex)
{
    TString cName = Form("waveform_analysis_run%d", runIndex);
    TCanvas* c = new TCanvas(cName,
        Form("Entry #%lld -- LP @ %.1f MHz / thr %.1f mV  [Run %d]",
             entry, cutoff_MHz, threshold_mV, runIndex),
        1100, 700);
    c->Divide(1, 2, 0.01, 0.01);

    auto drawPad = [&](int pad,
                       const std::vector<double>& raw,
                       const std::vector<double>& filt,
                       const PulseMetrics& m,
                       const char* ch)
    {
        c->cd(pad);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.14);
        gStyle->SetOptStat(0);
        gStyle->SetPadGridX(true);
        gStyle->SetPadGridY(true);

        TGraph* grR = new TGraph(N_SAMPLES, time.data(), raw.data());
        TGraph* grF = new TGraph(N_SAMPLES, time.data(), filt.data());
        grR->SetLineColor(kGray + 1); grR->SetLineWidth(1);
        grF->SetLineColor(pad == 1 ? kAzure - 4 : kOrange + 7); grF->SetLineWidth(2);

        grR->SetTitle(Form(
            "%s  |  Entry #%lld  LP f_{c}=%.1f MHz  thr=%.1f mV  (Run %d)"
            ";time (ns);amplitude (mV)", ch, entry, cutoff_MHz, threshold_mV, runIndex));
        grR->GetXaxis()->SetTitleSize(0.05);
        grR->GetYaxis()->SetTitleSize(0.05);
        grR->GetYaxis()->SetTitleOffset(0.9);
        grR->Draw("AL");
        grF->Draw("L SAME");

        // Threshold line
        double xMin = time.front(), xMax = time.back();
        double thr = m.peak_amp_mV < 0
                     ? -std::abs(threshold_mV) : std::abs(threshold_mV);
        TLine* thrLine = new TLine(xMin, thr, xMax, thr);
        thrLine->SetLineColor(kRed); thrLine->SetLineWidth(1);
        thrLine->SetLineStyle(2);    thrLine->Draw();

        TLegend* leg = new TLegend(0.55, 0.74, 0.94, 0.90);
        leg->SetBorderSize(1); leg->SetFillColor(kWhite); leg->SetTextSize(0.040);
        leg->AddEntry(grR,     "Raw (baseline-corrected)", "l");
        leg->AddEntry(grF,     Form("Filtered f_{c}=%.1f MHz", cutoff_MHz), "l");
        leg->AddEntry(thrLine, Form("Threshold %.1f mV", thr), "l");
        leg->Draw();

        // Info box
        TPaveText* info = new TPaveText(0.12, 0.74, 0.52, 0.90, "brNDC");
        info->SetBorderSize(1); info->SetFillColor(kWhite); info->SetFillStyle(1001);
        info->SetTextFont(42);  info->SetTextSize(0.036); info->SetTextAlign(12);
        info->AddText(Form("Entry %lld/%lld  f_{s}=%.0f MHz", entry, nEntries, fs_MHz));
        if (m.found) {
            info->AddText(Form("Peak=%.2f mV  @%.1f ns", m.peak_amp_mV, m.peak_time_ns));
            info->AddText(Form("Integral=%.1f mV#upointns  Rise=%.2f ns  Fall=%.2f ns",
                               m.integral_mVns, m.rise_time_ns, m.fall_time_ns));
        } else {
            info->AddText("No pulse above threshold");
        }
        info->Draw();

        gPad->Update();
    };

    drawPad(1, amp1_corr, amp1_filt, m1, "ch1");
    drawPad(2, amp2_corr, amp2_filt, m2, "laser");
    c->Update();
    return c;
}

// ─────────────────────────────────────────────────────────────
// Print metrics summary to stdout
// ─────────────────────────────────────────────────────────────
static void printMetrics(const PulseMetrics& m, const std::string& label)
{
    std::cout << "\n  ── " << label << " ──────────────────────────────" << std::endl;
    if (!m.found) {
        std::cout << "  [!] No pulse crossing threshold." << std::endl;
        std::cout << "  Noise RMS: " << std::fixed << std::setprecision(3)
                  << m.noise_rms_mV << " mV" << std::endl;
        return;
    }
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  Peak amplitude : " << m.peak_amp_mV    << " mV" << std::endl;
    std::cout << "  Peak time      : " << m.peak_time_ns   << " ns" << std::endl;
    std::cout << "  Integral       : " << m.integral_mVns  << " mV·ns" << std::endl;
    std::cout << "  Rise time      : " << m.rise_time_ns   << " ns (10%-90%)" << std::endl;
    std::cout << "  Fall time      : " << m.fall_time_ns   << " ns (90%-10%)" << std::endl;
    std::cout << "  Thr crossing   : " << m.thresh_cross_ns<< " ns" << std::endl;
    std::cout << "  N crossings    : " << m.n_crossings    << std::endl;
    std::cout << "  Noise RMS      : " << m.noise_rms_mV   << " mV" << std::endl;
}

// ─────────────────────────────────────────────────────────────
// PHASE 2 – Analysis loop
// ─────────────────────────────────────────────────────────────
static void analysisLoop(
    const std::vector<double>& time,
    const std::vector<double>& amp1_corr,
    const std::vector<double>& amp2_corr,
    Long64_t entry,
    Long64_t nEntries,
    double fs_MHz,
    std::ofstream& csvFile)
{
    int runIndex = 1;

    while (true) {
        // ── User inputs ────────────────────────────────────────
        double cutoff_MHz   = promptPositiveDouble(
            "\n  Enter low-pass cutoff frequency [MHz]: ");
        double threshold_mV = promptPositiveDouble(
            "  Enter threshold (absolute value) [mV]: ");

        if (cutoff_MHz >= fs_MHz / 2.0)
            std::cerr << "[WARNING] Cutoff >= Nyquist ("
                      << fs_MHz / 2.0 << " MHz). Filter may be ineffective." << std::endl;

        // ── Filter both channels ───────────────────────────────
        std::vector<double> amp1_filt = butterworthLowPass(amp1_corr, cutoff_MHz, fs_MHz, 4);
        std::vector<double> amp2_filt = butterworthLowPass(amp2_corr, cutoff_MHz, fs_MHz, 4);

        // ── Extract metrics ────────────────────────────────────
        PulseMetrics m1 = extractMetrics(time, amp1_filt, amp1_corr, threshold_mV, "ch1");
        PulseMetrics m2 = extractMetrics(time, amp2_filt, amp2_corr, threshold_mV, "laser");

        printMetrics(m1, "ch1 results");
        printMetrics(m2, "laser results");

        // ── Draw canvas ────────────────────────────────────────
        TCanvas* c = drawAnalysisCanvas(
            time,
            amp1_corr, amp1_filt,
            amp2_corr, amp2_filt,
            m1, m2,
            entry, nEntries,
            cutoff_MHz, threshold_mV, fs_MHz,
            runIndex);

        // ── Save PNG ───────────────────────────────────────────
        TString fname = Form("waveform_entry%lld_run%d_fc%.0fMHz_thr%.0fmV.png",
                             entry, runIndex, cutoff_MHz, threshold_mV);
        c->SaveAs(fname);
        std::cout << "\n[INFO] Canvas saved as: " << fname << std::endl;

        // ── Write CSV row ──────────────────────────────────────
        ResultRow row;
        row.entry        = entry;
        row.run          = runIndex;
        row.cutoff_MHz   = cutoff_MHz;
        row.threshold_mV = threshold_mV;
        row.ch1          = m1;
        row.laser        = m2;
        writeCSVRow(csvFile, row);
        std::cout << "[INFO] Results appended to CSV." << std::endl;

        // ── Continue? ──────────────────────────────────────────
        printBox({"Run another analysis with a different cutoff / threshold?"});
        if (askYesNo("  Your choice [y/n]: ") == 'n') break;
        ++runIndex;
    }
}

// ─────────────────────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────────────────────
void waveform_analysis()
{
    printBox({
        "Waveform Analysis  --  Interactive Mode",
        "Dual-channel  |  Low-pass filter  |  Threshold  |  Amplitude & Integral"
    });

    // ── Open ROOT file ─────────────────────────────────────────
    TFile* file = TFile::Open("../../data/data.vbias_{53}.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Cannot open data.vbias_{53}.root" << std::endl;
        return;
    }

    // ── Load both TTrees ───────────────────────────────────────
    TTree* tree1 = (TTree*)file->Get("ch1");
    TTree* tree2 = (TTree*)file->Get("laser");
    if (!tree1) { std::cerr << "[ERROR] TTree 'ch1' not found." << std::endl; file->Close(); return; }
    if (!tree2) { std::cerr << "[ERROR] TTree 'laser' not found." << std::endl; file->Close(); return; }

    Long64_t nEntries = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    if (nEntries == 0) { std::cerr << "[ERROR] TTree 'ch1' is empty." << std::endl; file->Close(); return; }
    if (nEntries2 != nEntries)
        std::cerr << "[WARNING] ch1 and laser have different number of entries ("
                  << nEntries << " vs " << nEntries2 << "). Using min." << std::endl;
    nEntries = std::min(nEntries, nEntries2);

    // ── Open CSV output ────────────────────────────────────────
    std::string csvName = "waveform_analysis_results.csv";
    std::ofstream csvFile(csvName, std::ios::out | std::ios::trunc);
    if (!csvFile.is_open()) {
        std::cerr << "[ERROR] Cannot open " << csvName << " for writing." << std::endl;
        file->Close(); return;
    }
    writeCSVHeader(csvFile);
    std::cout << "\n[INFO] Results will be saved to: " << csvName << std::endl;

    // ── PHASE 1: waveform selection ────────────────────────────
    TRandom3 rng(0);
    std::cout << "\n  [PHASE 1]  Dual-channel waveform selection" << std::endl;

    std::vector<double> time, amp1_corr, amp2_corr;
    Long64_t entry;
    double fs_MHz = 0.0;

    if (!previewLoop(tree1, tree2, nEntries, rng,
                     time, amp1_corr, amp2_corr, entry, fs_MHz)) {
        std::cerr << "[ERROR] Preview loop failed." << std::endl;
        csvFile.close(); file->Close(); return;
    }

    // ── PHASE 2: analysis ─────────────────────────────────────
    std::cout << "\n  [PHASE 2]  Analysis  --  Entry #" << entry
              << "  |  f_s = " << fs_MHz << " MHz" << std::endl;

    analysisLoop(time, amp1_corr, amp2_corr,
                 entry, nEntries, fs_MHz, csvFile);

    // ── Cleanup ────────────────────────────────────────────────
    csvFile.close();
    std::cout << "\n[INFO] CSV closed: " << csvName << std::endl;
    file->Close();
    cleanupTempFiles();

    std::cout << "\n  Goodbye! All canvases remain open until you close ROOT.\n" << std::endl;
}