#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include "ButterworthFilter.h"  // scipy-equivalent Butterworth IIR

static std::vector<double> correctBaseline(const std::vector<double>& time, const std::vector<double>& amp, double pre_signal_end = 30.0) {
    std::vector<double> pre;
    for (size_t i = 0; i < time.size(); ++i) if (time[i] <= pre_signal_end) pre.push_back(amp[i]);
    std::vector<double> tmp = pre;
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
    double offset = tmp[tmp.size()/2];
    std::vector<double> out = amp;
    for (auto& a : out) a -= offset;
    return out;
}

// butterworthLowPass is provided by ButterworthFilter.h
// To use zero-phase (sosfiltfilt): butterworthLowPassZP(x, fc, fs, 4)
// To use fast pre-computed 500MHz:  butterworthLP_500MHz(x)

void threshold_scan2() {
    gStyle->SetOptStat(0);

    TFile* file = TFile::Open("data.vbias_{55}.root", "READ");
    if (!file || file->IsZombie()) { std::cerr << "Errore apertura file!\n"; return; }

    TTree* treeCh1 = (TTree*)file->Get("ch1");
    if (!treeCh1) { std::cerr << "Tree ch1 non trovato!\n"; return; }

    const int N = 1024;
    Double_t t1[N], a1[N];
    treeCh1->SetBranchAddress("time", t1);
    treeCh1->SetBranchAddress("amplitude", a1);

    treeCh1->GetEntry(0);
    double fs_MHz = 1000.0 / (t1[1] - t1[0]);
    std::cout << "Sampling frequency: " << fs_MHz << " MHz\n";

    double cutoff_MHz;
    std::cout << "Frequenza di taglio LP [MHz]: ";
    std::cin >> cutoff_MHz;

    // FIX: virgola mancante corretta
    double min_thr = -100.0, max_thr = 100.0, step = 0.1;
    int n_steps = (int)std::round((max_thr - min_thr) / step);
    std::vector<double> thresholds(n_steps + 1), counts(n_steps + 1, 0.0);
    for (int k = 0; k <= n_steps; ++k) thresholds[k] = min_thr + k * step;

    Long64_t nEntries = treeCh1->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\rProgress: " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);

        std::vector<double> v_time(t1, t1 + N);
        std::vector<double> v_amp(a1, a1 + N);
        std::vector<double> amp_filt = butterworthLowPass(
            correctBaseline(v_time, v_amp), cutoff_MHz, fs_MHz);

        for (int k = 0; k <= n_steps; ++k) {
            double threshold = thresholds[k];
            double sign = (threshold >= 0) ? 1.0 : -1.0;
            bool armed = false, crossed = false;

            for (int j = 0; j < N; ++j) {
                double y = amp_filt[j];
                if (!armed) {
                    if (sign * y < sign * threshold * 0.5) armed = true;
                } else {
                    if (sign * y > sign * threshold) {
                        crossed = true;
                        armed = false;
                    }
                }
            }
            if (crossed) counts[k] += 1.0;
        }
    }

    std::cout << "\nCounts sample: ";
    for (int k = 0; k <= n_steps && k < 10; ++k)
        std::cout << thresholds[k] << "mV=" << counts[k] << "  ";
    std::cout << "\n";

    // Chiudi il file PRIMA di disegnare (il TGraph ha già i dati copiati)
    file->Close();

    TCanvas* cScan = new TCanvas("cScan", "Threshold Scan with Arming Logic", 800, 600);
    cScan->SetGrid();
    TGraph* gr = new TGraph(n_steps + 1, thresholds.data(), counts.data());
    gr->SetTitle("Cumulative Scan (Arming Logic 50%);Threshold (mV);Counts");
    gr->SetLineColor(kAzure + 1);
    gr->SetLineWidth(2);
    gr->Draw("AL");
    cScan->Update();
    cScan->Modified();

    // Salva sempre su file per sicurezza
    cScan->SaveAs("threshold_scan.png");
    cScan->SaveAs("threshold_scan.root");

    std::cout << "\nDONE. Grafico salvato. Chiudi la finestra per uscire.\n";
    // Niente gSystem->Run() → ROOT rimane responsivo
}