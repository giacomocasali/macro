#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

/**
 * Official code language: English [cite: 2026-02-26]
 * Self-contained version to avoid Utils.h dependency errors.
 */

// --- INLINED UTILS ---
static Double_t* setupBranch_dgz(TTree* tree, const char* branchName) {
    Double_t* arr = new Double_t[1024];
    tree->SetBranchAddress(branchName, arr);
    return arr;
}

static std::vector<std::pair<std::vector<double>, std::vector<double>>> correctWaveforms(
    const std::vector<std::pair<std::vector<double>, std::vector<double>>>& waveforms,
    double pre_signal_end = 30.0
) {
    std::vector<std::pair<std::vector<double>, std::vector<double>>> corrected;
    corrected.reserve(waveforms.size());
    for (const auto& entry : waveforms) {
        const std::vector<double>& time = entry.first;
        const std::vector<double>& amp  = entry.second;
        std::vector<double> pre_amp;
        for (size_t i = 0; i < time.size(); ++i)
            if (time[i] < pre_signal_end) pre_amp.push_back(amp[i]);
        // Fallback: if no samples in baseline window, use first 10% of waveform
        if (pre_amp.empty()) {
            size_t nPre = std::max((size_t)10, time.size() / 10);
            for (size_t i = 0; i < nPre && i < amp.size(); ++i)
                pre_amp.push_back(amp[i]);
        }
        double offset = 0.0;
        if (!pre_amp.empty()) {
            std::vector<double> tmp = pre_amp;
            size_t mid = tmp.size() / 2;
            std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
            offset = tmp[mid];  // median, robust to outliers
        }
        std::vector<double> corrected_amp = amp;
        for (auto& a : corrected_amp) a -= offset;
        corrected.emplace_back(time, corrected_amp);
    }
    return corrected;
}

void waveform_2D_alignment() {
    // Apri il file (assicurati che il nome sia corretto per la tua cartella)
    TFile *file = TFile::Open("data.vbias_{54}.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open source file!" << std::endl;
        return;
    }

    TTree *treeCh1 = (TTree*)file->Get("ch1");
    TTree *treeLaser = (TTree*)file->Get("laser");

    if (!treeCh1 || !treeLaser) {
        std::cerr << "Error: Trees (ch1 or laser) not found!" << std::endl;
        return;
    }

    Double_t* timeCh1 = setupBranch_dgz(treeCh1, "time");
    Double_t* ampCh1  = setupBranch_dgz(treeCh1, "amplitude");
    Double_t* timeLaser = setupBranch_dgz(treeLaser, "time");
    Double_t* ampLaser  = setupBranch_dgz(treeLaser, "amplitude");

    Long64_t nEntries = treeCh1->GetEntries();
    
    // TH2D: Asse X = Tempo relativo [ns], Asse Y = Ampiezza corretta [mV]
    TH2D *hCh1Corr = new TH2D("hCh1Corr", "Aligned Waveforms (Baseline Corrected);Time relative to Laser [ns];Amplitude [mV]", 
                              600, -20, 80, 500, -10, 250);

    const double threshold = 50.0; // Soglia trigger per il laser
    const int samples = 1024;

    for (Long64_t i = 0; i < nEntries; i++) {
        treeCh1->GetEntry(i);
        treeLaser->GetEntry(i);

        // Prepariamo i dati per la correzione baseline
        std::vector<std::pair<std::vector<double>, std::vector<double>>> tempCh1, tempLaser;
        std::vector<double> t1(samples), a1(samples), tL(samples), aL(samples);
        for(int s=0; s<samples; s++){
            t1[s]=timeCh1[s]; a1[s]=ampCh1[s];
            tL[s]=timeLaser[s]; aL[s]=ampLaser[s];
        }
        tempCh1.push_back({t1, a1});
        tempLaser.push_back({tL, aL});

        // Applichiamo la sottrazione dell'offset (baseline)
        auto correctedCh1 = correctWaveforms(tempCh1, 30.0);
        auto correctedLaser = correctWaveforms(tempLaser, 30.0);

        // Troviamo il tempo del laser (trigger) per allineare le forme d'onda
        double tLaserTrigger = -999.0;
        for (size_t j = 1; j < correctedLaser[0].first.size(); j++) {
            double v1 = correctedLaser[0].second[j-1];
            double v2 = correctedLaser[0].second[j];
            if (v2 > threshold && v1 <= threshold) {
                double t1_l = correctedLaser[0].first[j-1];
                double t2_l = correctedLaser[0].first[j];
                // Interpolazione lineare per precisione sub-sample
                tLaserTrigger = t1_l + (threshold - v1) * (t2_l - t1_l) / (v2 - v1);
                break; 
            }
        }

        // Se il trigger è valido, riempiamo l'istogramma 2D
        if (tLaserTrigger != -999.0) {
            for (size_t j = 0; j < correctedCh1[0].first.size(); j++) {
                double tRel = correctedCh1[0].first[j] - tLaserTrigger;
                double amp  = correctedCh1[0].second[j];
                hCh1Corr->Fill(tRel, amp);
            }
        }
        
        if (i % 2000 == 0) std::cout << "Event " << i << " processed..." << std::endl;
    }

    // Visualizzazione
    gStyle->SetPalette(kBird);
    TCanvas *c1 = new TCanvas("c1", "Waveform Alignment", 900, 700);
    gPad->SetLogz(); // Scala logaritmica per vedere meglio i segnali rari
    hCh1Corr->Draw("COLZ");
}