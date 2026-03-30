#include "Utils.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLegend.h>

// Official code language: English [cite: 2026-02-26]

Double_t multiGauss(Double_t *x, Double_t *par) {
    int n = (int)par[0];
    double sum = 0;
    for(int i=0; i<n; i++) sum += par[1+3*i] * TMath::Gaus(x[0], par[2+3*i], par[3+3*i]);
    return sum;
}

void sipm_risetime_analysis() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile *file = TFile::Open("../../data/data.vbias_{40}.root", "READ");
    if(!file || file->IsZombie()) return;

    TTree *treeCh1 = (TTree*)file->Get("ch1");
    TTree *treeLaser = (TTree*)file->Get("laser");

    Double_t* t1 = Utils::setupBranch_dgz(treeCh1, "time");
    Double_t* a1 = Utils::setupBranch_dgz(treeCh1, "amplitude");
    Double_t* tL = Utils::setupBranch_dgz(treeLaser, "time");
    Double_t* aL = Utils::setupBranch_dgz(treeLaser, "amplitude");

    // --- STEP 1: SPECTRUM & CALIBRATION ---
    TH1D *hSpec = new TH1D("hSpec", "Amplitude Spectrum;Amplitude (mV);Counts", 600, -5, 60);
    
    std::cout << "Processing Spectrum..." << std::endl;
    for(Long64_t i=0; i<treeCh1->GetEntries(); i++) {
        if(i % 500 == 0 && gSystem->ProcessEvents()) break; // Esc [cite: 2026-03-04]
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);
        
        auto c1 = Utils::correctWaveforms({{{t1, t1+1024}, {a1, a1+1024}}}, 30.0);
        auto cL = Utils::correctWaveforms({{{tL, tL+1024}, {aL, aL+1024}}}, 30.0);
        
        double tTrig = -999;
        for(size_t j=1; j<cL[0].first.size(); j++) {
            if(cL[0].second[j] > 10 && cL[0].second[j-1] <= 10) {
                tTrig = cL[0].first[j-1] + (10-cL[0].second[j-1])*(cL[0].first[j]-cL[0].first[j-1])/(cL[0].second[j]-cL[0].second[j-1]);
                break;
            }
        }
        if(tTrig == -999) continue;

        double maxA = -999;
        for(size_t j=0; j<c1[0].first.size(); j++) {
            double tRel = c1[0].first[j] - tTrig;
            if(tRel >= 45 && tRel <= 55) if(c1[0].second[j] > maxA) maxA = c1[0].second[j];
        }
        if(maxA > -900) hSpec->Fill(maxA);
    }

    TSpectrum s(20);
    int nFound = s.Search(hSpec, 1.5, "goff", 0.005);
    double *xPos = s.GetPositionX();
    std::vector<double> peaks(xPos, xPos + nFound);
    std::sort(peaks.begin(), peaks.end());

    int nPeaks = 8; 
    TF1 *fMulti = new TF1("fMulti", multiGauss, -5, 55, 1+3*nPeaks);
    fMulti->SetNpx(1000);
    fMulti->FixParameter(0, nPeaks);
    for(int i=0; i<nPeaks; i++) {
        double x_est = (i < (int)peaks.size()) ? peaks[i] : (1.5 + 4.3 * i);
        fMulti->SetParameter(1+3*i, hSpec->GetBinContent(hSpec->FindBin(x_est)));
        fMulti->SetParameter(2+3*i, x_est);
        fMulti->SetParLimits(2+3*i, x_est - 1.0, x_est + 1.0); // tightened
        fMulti->SetParameter(3+3*i, 0.6);
        fMulti->SetParLimits(3+3*i, 0.4, 1.1); // sigma limits added
    }
    TVirtualFitter::SetDefaultFitter("Minuit2");
    hSpec->Fit(fMulti, "RQML0"); // Maximum Likelihood, 0 = don't draw

    TGraphErrors *grLin = new TGraphErrors();
    for(int i=0; i<nPeaks; i++) {
        grLin->SetPoint(i, i, fMulti->GetParameter(2+3*i));
        grLin->SetPointError(i, 0, fMulti->GetParError(2+3*i));
    }
    TF1 *fLin = new TF1("fLin", "pol1", -0.5, nPeaks-0.5);
    grLin->Fit(fLin, "RQ0");
    
    double m = fLin->GetParameter(1);
    double q = fLin->GetParameter(0);
    double let = (m / 2.0) + q;

    // --- STEP 2: RISE TIME COLLECTION ---
    TH2D *hRise2D = new TH2D("hRise2D", "Rise time vs amplitude;amplitude (mV);t_{rise} (ns)", 400, -5, 60, 250, -1, 10);
    std::vector<TH1D*> hProjections;
    for(int n=0; n<nPeaks; n++) {
        hProjections.push_back(new TH1D(Form("hProj_%d", n), Form("%d p.e. rise time;#Deltat (ns);counts", n), 100, 0, 10));
    }

    std::cout << "Collecting Rise Time data..." << std::endl;
    for(Long64_t i=0; i<treeCh1->GetEntries(); i++) {
        if(i % 500 == 0 && gSystem->ProcessEvents()) break; 
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);
        auto c1 = Utils::correctWaveforms({{{t1, t1+1024}, {a1, a1+1024}}}, 30.0);
        auto cL = Utils::correctWaveforms({{{tL, tL+1024}, {aL, aL+1024}}}, 30.0);
        
        double tTrig = -999;
        for(size_t j=1; j<cL[0].first.size(); j++) {
            if(cL[0].second[j] > 10 && cL[0].second[j-1] <= 10) {
                tTrig = cL[0].first[j-1] + (10-cL[0].second[j-1])*(cL[0].first[j]-cL[0].first[j-1])/(cL[0].second[j]-cL[0].second[j-1]);
                break;
            }
        }
        if(tTrig == -999) continue;

        const auto& times = c1[0].first; const auto& ampls = c1[0].second;
        double tLET = -999; size_t idxLET = 0;
        for(size_t j=1; j<times.size(); j++) {
            if(ampls[j] > let && ampls[j-1] <= let) {
                tLET = times[j-1] + (let - ampls[j-1])*(times[j]-times[j-1])/(ampls[j]-ampls[j-1]);
                idxLET = j; break;
            }
        }

        if(tLET != -999) {
            double tPeak = -999, aPeak = -999;
            for(size_t k = idxLET; k < ampls.size() - 10; k++) {
                bool isMax = true;
                for(size_t m_pt = 1; m_pt <= 10; m_pt++) if(ampls[k+m_pt] >= ampls[k]) { isMax = false; break; }
                if(isMax) { aPeak = ampls[k]; tPeak = times[k]; break; }
            }
            
            if(tPeak != -999) {
                double dt = tPeak - tLET;
                hRise2D->Fill(aPeak, dt);
                for(int n=0; n<nPeaks; n++) {
                    double low = q + (n - 0.5)*m;
                    double high = q + (n + 0.5)*m;
                    if(aPeak >= low && aPeak < high) {
                        hProjections[n]->Fill(dt);
                        break;
                    }
                }
            }
        }
    }

    // --- CANVAS 1: 2D DENSITY PLOT ---
    TCanvas *can1 = new TCanvas("can1", "2D Rise Time Analysis", 900, 700);
    can1->SetRightMargin(0.15); can1->SetLogz(); can1->SetGrid();
    hRise2D->Draw("COLZ");
    for(int n=0; n<nPeaks; n++) {
        double x_border = q + (n + 0.5)*m;
        if(x_border < 60) {
            TLine *l = new TLine(x_border, -1, x_border, 10);
            l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw();
        }
    }

    // --- CANVAS 2: 1D RAINBOW PROJECTIONS ---
    TCanvas *can2 = new TCanvas("can2", "Rise Time Projections", 900, 700);
    can2->SetGrid();
    TLegend *leg = new TLegend(0.55, 0.45, 0.88, 0.88);
    leg->SetBorderSize(1);
    
    // Rainbow Palette for 0 to 7 p.e.
    std::vector<int> rainbow = {kGray+2, kRed, kOrange+1, kSpring-1, kGreen+2, kAzure+7, kBlue+1, kViolet-3};

    double maxVal = 0;
    for(int n=0; n<nPeaks; n++) if(hProjections[n]->GetMaximum() > maxVal) maxVal = hProjections[n]->GetMaximum();

    for(int n=0; n<nPeaks; n++) {
        if(hProjections[n]->GetEntries() < 10) continue;
        int color = (n < (int)rainbow.size()) ? rainbow[n] : kBlack;
        
        hProjections[n]->SetLineColor(color);
        hProjections[n]->SetLineWidth(3);
        hProjections[n]->GetYaxis()->SetRangeUser(0, maxVal * 1.1);
        hProjections[n]->Draw(n==0 ? "HIST" : "HIST SAME");
        
        if(hProjections[n]->GetEntries() > 20) {
            leg->AddEntry(hProjections[n], Form("%d p.e. (Mean: %.2f ns)", n, hProjections[n]->GetMean()), "l");
        }
    }
    leg->Draw();

    std::cout << "\nAnalysis Complete.\nGain: " << m << " mV/pe\nLET Threshold: " << let << " mV" << std::endl;
}