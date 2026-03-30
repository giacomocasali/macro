#include "Utils.h"
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TStyle.h>

void pulse_interarrival_time() {
    gStyle->SetOptStat(1111);

    TFile *file = TFile::Open("../../data/calibration.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    TTree *treeCh1 = (TTree*)file->Get("ch1");
    if (!treeCh1) return;

    Double_t* timeCh1 = Utils::setupBranch_dgz(treeCh1, "time");
    Double_t* ampCh1  = Utils::setupBranch_dgz(treeCh1, "amplitude");

    // Range set to [0, 100] as per your latest margins
    TH1D *hDeltaT = new TH1D("hDeltaT", "Inter-arrival Time (First 2 Pulses); #Delta t [ns]; Counts", 80, 66, 67);

    const double threshold = 10.0; 
    const double deadTime  = 12.0; 
    const int samples      = 1024;

    Long64_t nEntries = treeCh1->GetEntries();
    int zeroPulses = 0;
    int singlePulses = 0;
    int multiPulses = 0;

    std::cout << "Processing " << nEntries << " waveforms..." << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        // EMERGENCY STOP: Press Esc or Ctrl+C
        if (i % 500 == 0 && gSystem->ProcessEvents()) break;

        treeCh1->GetEntry(i);

        std::vector<double> arrivalTimes;
        
        for (int j = 1; j < samples; j++) {
            double valPrev = ampCh1[j-1];
            double valCurr = ampCh1[j];

            if (valCurr > threshold && valPrev <= threshold) {
                double tCrossing = timeCh1[j] + (threshold - valCurr) * (timeCh1[j-1] - timeCh1[j]) / (valPrev - valCurr);


                std::cout << "crossing " << tCrossing << "\n";
                
                if (arrivalTimes.empty() || (tCrossing - arrivalTimes.back() > deadTime)) {
                    arrivalTimes.push_back(tCrossing);
                    std::cout << "arrival time " << tCrossing - arrivalTimes.back() << "\n";
                }
            }
            // Stop after finding the first two pulses for performance
            if (arrivalTimes.size() == 2) break; 
        }

        if (arrivalTimes.size() >= 2) {
            hDeltaT->Fill(arrivalTimes[1] - arrivalTimes[0]);
            multiPulses++;
            std::cout << " delta arrival time " << arrivalTimes[1] - arrivalTimes[0] << "\n";
        } else if (arrivalTimes.size() == 1) {
            singlePulses++;
        } else {
            zeroPulses++;
        }

        if (i % 10000 == 0) printf("\rProgress: %.1f%%", (double)i/nEntries*100.0);
    }

    std::cout << "\n\n=== ANALYSIS REPORT ===" << std::endl;
    std::cout << "Total entries:  " << nEntries << std::endl;
    std::cout << "Valid pairs:    " << multiPulses << std::endl;
    std::cout << "Single pulses:  " << singlePulses << std::endl;
    std::cout << "Zero pulses:    " << zeroPulses << std::endl;
    std::cout << "=======================\n" << std::endl;

    TCanvas *c1 = new TCanvas("c1", "Timing Analysis", 900, 900);
    
    // Set Log Scale on Y axis
    c1->SetLogy();
    
    hDeltaT->SetLineColor(kRed+1);
    hDeltaT->SetFillColor(kOrange-9);
    hDeltaT->Draw("HIST");

    delete[] timeCh1; 
    delete[] ampCh1;
}
