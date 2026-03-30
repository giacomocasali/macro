// waveform_lowpass.cpp
// Usage:
//   root -l 'waveform_lowpass.cpp+'
//   then call:  waveform_lowpass()

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
#include <TGraph.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "ButterworthFilter.h"  // scipy-equivalent Butterworth IIR

// ─────────────────────────────────────────────────────────────
// Repeat a string n times (used for box drawing with ASCII)
// ─────────────────────────────────────────────────────────────
static std::string repeatStr(const std::string& s, size_t n)
{
    std::string out;
    out.reserve(s.size() * n);
    for (size_t i = 0; i < n; ++i) out += s;
    return out;
}

// ─────────────────────────────────────────────────────────────
// Pretty-print a uniform ASCII box with dynamic width
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
// Baseline correction: subtract median of samples before t_pre
// ─────────────────────────────────────────────────────────────
static std::vector<double> correctBaseline(
    const std::vector<double>& time,
    const std::vector<double>& amp,
    double pre_signal_end = 0.0)
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
// 4th-order Butterworth low-pass filter (IIR, biquad cascade)
// ─────────────────────────────────────────────────────────────
// butterworthLowPass is provided by ButterworthFilter.h
// To use zero-phase (sosfiltfilt): butterworthLowPassZP(x, fc, fs, 4)
// To use fast pre-computed 500MHz:  butterworthLP_500MHz(x)

// ─────────────────────────────────────────────────────────────
// Read a yes/no answer from stdin
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

// ─────────────────────────────────────────────────────────────
// Prompt for a positive cutoff frequency [MHz]
// ─────────────────────────────────────────────────────────────
static double promptCutoff()
{
    double fc = -1.0;
    while (fc <= 0.0) {
        std::cout << "\n  Enter low-pass cutoff frequency [MHz]: ";
        std::cin >> fc;
        if (std::cin.fail() || fc <= 0.0) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "  [ERROR] Please enter a positive number." << std::endl;
            fc = -1.0;
        }
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return fc;
}

// ─────────────────────────────────────────────────────────────
// Load a random entry from the TTree
// ─────────────────────────────────────────────────────────────
static bool loadRandomEntry(
    TTree* tree,
    Long64_t nEntries,
    TRandom3& rng,
    int N,
    std::vector<double>& time,
    std::vector<double>& amp,
    Long64_t& entryIndex)
{
    static Double_t t_arr[1024], a_arr[1024];
    tree->SetBranchAddress("time",      t_arr);
    tree->SetBranchAddress("amplitude", a_arr);

    entryIndex = (Long64_t)(rng.Uniform(0, (double)nEntries));
    tree->GetEntry(entryIndex);

    time.assign(t_arr, t_arr + N);
    amp .assign(a_arr, a_arr + N);
    return true;
}

// ─────────────────────────────────────────────────────────────
// PHASE 1 -- Preview loop
// ─────────────────────────────────────────────────────────────
static bool previewLoop(
    TTree* tree,
    Long64_t nEntries,
    TRandom3& rng,
    int N,
    std::vector<double>& out_time,
    std::vector<double>& out_amp_corr,
    Long64_t& out_entry)
{
    TCanvas* cPrev = nullptr;
    int attempt = 0;

    while (true) {
        ++attempt;

        std::vector<double> time, amp;
        Long64_t entry;
        if (!loadRandomEntry(tree, nEntries, rng, N, time, amp, entry))
            return false;

        std::vector<double> amp_corr = correctBaseline(time, amp, 30.0);

        if (cPrev) { cPrev->Close(); cPrev = nullptr; }

        TGraph* grPrev = new TGraph(N, time.data(), amp_corr.data());
        grPrev->SetLineColor(kTeal - 5);
        grPrev->SetLineWidth(2);

        double dt_ns  = (time.back() - time.front()) / (N - 1);
        double fs_MHz = 1000.0 / dt_ns;

        cPrev = new TCanvas("waveform_preview",
            Form("PREVIEW -- Entry #%lld (attempt %d)", entry, attempt),
            1000, 550);
        cPrev->SetLeftMargin(0.12);
        cPrev->SetRightMargin(0.05);
        cPrev->SetBottomMargin(0.12);

        gStyle->SetOptStat(0);
        gStyle->SetPadGridX(true);
        gStyle->SetPadGridY(true);

        grPrev->SetTitle(Form(
            "PREVIEW -- Entry #%lld  (attempt %d)  |  f_{s} = %.0f MHz"
            ";time (ns);amplitude (mV)", entry, attempt, fs_MHz));
        grPrev->GetXaxis()->SetTitleSize(0.045);
        grPrev->GetYaxis()->SetTitleSize(0.045);
        grPrev->GetYaxis()->SetTitleOffset(1.1);
        grPrev->Draw("AL");

        TPaveText* info = new TPaveText(0.13, 0.78, 0.48, 0.90, "brNDC");
        info->SetBorderSize(1);
        info->SetFillColor(kWhite);
        info->SetFillStyle(1001);
        info->SetTextFont(42);
        info->SetTextSize(0.032);
        info->SetTextAlign(12);
        info->AddText(Form("Entry:  %lld / %lld", entry, nEntries));
        info->AddText(Form("f_{s} = %.0f MHz   |   Baseline corrected", fs_MHz));
        info->AddText(Form("Attempt: %d", attempt));
        info->Draw();

        cPrev->Update();

        printBox({
            std::string(Form("PREVIEW  --  Entry #%lld  (attempt %d)", entry, attempt)),
            "Is this waveform good for analysis?"
        });

        char ans = askYesNo("  Your choice [y/n]: ");

        if (ans == 'y') {
            out_time     = time;
            out_amp_corr = amp_corr;
            out_entry    = entry;
            return true;
        }
    }
}

// ─────────────────────────────────────────────────────────────
// PHASE 2 -- Analysis loop
// ─────────────────────────────────────────────────────────────
static void analysisLoop(
    const std::vector<double>& time,
    const std::vector<double>& amp_corr,
    Long64_t entry,
    Long64_t nEntries,
    int N)
{
    double dt_ns  = (time.back() - time.front()) / (N - 1);
    double fs_MHz = 1000.0 / dt_ns;

    int runIndex = 1;

    while (true) {
        double cutoff_MHz = promptCutoff();

        if (cutoff_MHz >= fs_MHz / 2.0)
            std::cerr << "[WARNING] Cutoff >= Nyquist (" << fs_MHz / 2.0
                      << " MHz). Filter will have no effect." << std::endl;

        std::vector<double> amp_filt = butterworthLowPass(amp_corr, cutoff_MHz, fs_MHz, 4);

        TGraph* grRaw  = new TGraph(N, time.data(), amp_corr.data());
        TGraph* grFilt = new TGraph(N, time.data(), amp_filt.data());

        grRaw->SetLineColor(kGray + 1);
        grRaw->SetLineWidth(1);
        grFilt->SetLineColor(kAzure - 4);
        grFilt->SetLineWidth(2);

        TString cName = Form("waveform_analysis_run%d", runIndex);
        TCanvas* c = new TCanvas(cName,
            Form("Entry #%lld -- LP @ %.1f MHz  [Run %d]", entry, cutoff_MHz, runIndex),
            1000, 550);
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.12);

        grRaw->SetTitle(Form(
            "Entry #%lld  |  Butterworth LP  f_{c} = %.1f MHz  (f_{s} = %.0f MHz)"
            ";time (ns);amplitude (mV)", entry, cutoff_MHz, fs_MHz));
        grRaw->GetXaxis()->SetTitleSize(0.045);
        grRaw->GetYaxis()->SetTitleSize(0.045);
        grRaw->GetYaxis()->SetTitleOffset(1.1);

        grRaw->Draw("AL");
        grFilt->Draw("L SAME");

        TLegend* leg = new TLegend(0.55, 0.78, 0.93, 0.90);
        leg->SetBorderSize(1);
        leg->SetFillColor(kWhite);
        leg->SetFillStyle(1001);
        leg->SetTextSize(0.033);
        leg->AddEntry(grRaw,  "Raw (baseline-corrected)", "l");
        leg->AddEntry(grFilt, Form("Filtered  f_{c} = %.1f MHz", cutoff_MHz), "l");
        leg->Draw();

        TPaveText* info = new TPaveText(0.13, 0.78, 0.48, 0.90, "brNDC");
        info->SetBorderSize(1);
        info->SetFillColor(kWhite);
        info->SetFillStyle(1001);
        info->SetTextFont(42);
        info->SetTextSize(0.030);
        info->SetTextAlign(12);
        info->AddText(Form("Entry:  %lld / %lld", entry, nEntries));
        info->AddText(Form("f_{s} = %.0f MHz   f_{c} = %.1f MHz", fs_MHz, cutoff_MHz));
        info->AddText(Form("Run %d  |  4^{th}-order Butterworth LP", runIndex));
        info->Draw();

        c->Update();

        TString fname = Form("waveform_entry%lld_run%d_fc%.0fMHz.png",
                             entry, runIndex, cutoff_MHz);
        c->SaveAs(fname);
        std::cout << "[INFO] Canvas saved as: " << fname << std::endl;

        printBox({"Run another analysis with a different cutoff?"});
        char ans = askYesNo("  Your choice [y/n]: ");
        if (ans == 'n') break;

        ++runIndex;
    }
}

// ─────────────────────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────────────────────
void waveform_lowpass()
{
    printBox({"Waveform Low-Pass Filter  --  Interactive Mode"});

    TFile* file = TFile::Open("data.vbias_{55}.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Cannot open data.vbias_{55}.root" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("ch1");
    if (!tree) {
        std::cerr << "[ERROR] TTree 'ch1' not found." << std::endl;
        file->Close();
        return;
    }

    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[ERROR] TTree 'ch1' is empty." << std::endl;
        file->Close();
        return;
    }

    const int N = 1024;
    TRandom3 rng(0);

    std::cout << "\n  [PHASE 1]  Waveform selection" << std::endl;

    std::vector<double> time, amp_corr;
    Long64_t entry;

    if (!previewLoop(tree, nEntries, rng, N, time, amp_corr, entry)) {
        std::cerr << "[ERROR] Preview loop failed." << std::endl;
        file->Close();
        return;
    }

    std::cout << "\n  [PHASE 2]  Low-pass filter analysis  --  Entry #"
              << entry << std::endl;

    analysisLoop(time, amp_corr, entry, nEntries, N);

    file->Close();

    std::cout << "\n  Done.\n" << std::endl;
}