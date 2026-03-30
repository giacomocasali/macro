#include "Utils.h"
#include "gauss_stuff.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TLine.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TLegend.h>

// Official code language: English [cite: 2026-02-26]

Double_t multiGauss(Double_t *x, Double_t *par) {
    int n = (int)par[0];
    double sum = 0;
    for(int i=0; i<n; i++) {
        sum += par[1+3*i] * TMath::Gaus(x[0], par[2+3*i], par[3+3*i]);
    }
    return sum;
}

void sipm_timing_LET_vs_CFD() {
    gStyle->SetOptFit(1111); // Show fit parameters
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetLabelSize(0.04, "XY");

    TFile *file = TFile::Open("data.vbias_{40}.root", "READ");
    if(!file || file->IsZombie()) return;

    TTree *treeCh1   = (TTree*)file->Get("ch1");
    TTree *treeLaser = (TTree*)file->Get("laser");

    Double_t* t1 = Utils::setupBranch_dgz(treeCh1,   "time");
    Double_t* a1 = Utils::setupBranch_dgz(treeCh1,   "amplitude");
    Double_t* tL = Utils::setupBranch_dgz(treeLaser,  "time");
    Double_t* aL = Utils::setupBranch_dgz(treeLaser,  "amplitude");

    // =========================================================
    // 2. PASS 1: AMPLITUDE SPECTRUM
    // =========================================================
    TH1D *hAmp = new TH1D("hAmp", "SiPM Amplitude Spectrum;Amplitude (mV);Counts", 400, -5, 60);

    for(Long64_t i=0; i<treeCh1->GetEntries(); i++) {
        if(i % 500 == 0 && gSystem->ProcessEvents()) break; // Esc [cite: 2026-03-04]
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);

        auto waveforms = Utils::correctWaveforms({{{t1, t1+1024}, {a1, a1+1024}}}, 30.0, true);
        const auto& as = waveforms[0].second;

        // Simple max find for spectrum
        double maxA = *std::max_element(as.begin(), as.end());
        if(maxA > -2.0) hAmp->Fill(maxA);
    }

    TSpectrum s(15);
    int nFound = s.Search(hAmp, 2, "gauss", 0.01);
    double *xPos = s.GetPositionX();
    std::vector<double> peaks(xPos, xPos + nFound);
    std::sort(peaks.begin(), peaks.end());

    int nPeaks = 8; 
    TF1 *fitAmp = new TF1("fitAmp", multiGauss, -2, 45, 1+3*nPeaks);
    fitAmp->FixParameter(0, nPeaks);
    fitAmp->SetLineColor(kBlue);
    fitAmp->SetLineWidth(2);

    for(int i=0; i<nPeaks; i++) {
        double x_est = (i < (int)peaks.size()) ? peaks[i] : (1.5 + 4.5 * i);
        fitAmp->SetParameter(1+3*i, hAmp->GetBinContent(hAmp->FindBin(x_est)));
        fitAmp->SetParameter(2+3*i, x_est);
        fitAmp->SetParameter(3+3*i, 0.6);
    }
    hAmp->Fit(fitAmp, "RQ");

    TCanvas *cAmp = new TCanvas("cAmp", "Spectrum Analysis", 900, 700);
    cAmp->SetLogy();
    hAmp->SetMinimum(0.5);
    hAmp->Draw("HIST");
    fitAmp->Draw("SAME");

    // =========================================================
    // 3. LINEARITY (cLin.png) - Improved with "draw 3" style
    // =========================================================
    std::vector<double> xPE, yM, yME, xE(nPeaks, 0);
    for(int i=0; i<nPeaks; i++) {
        xPE.push_back(i+1); // Start from 1 p.e.
        yM.push_back(fitAmp->GetParameter(2+3*i));
        yME.push_back(fitAmp->GetParError(2+3*i));
    }
    TCanvas *cLin = new TCanvas("cLin", "Linearity Calibration", 800, 600);
    cLin->SetGrid();
    TGraphErrors *grL = new TGraphErrors(nPeaks, &xPE[0], &yM[0], &xE[0], &yME[0]);
    grL->SetTitle("SiPM Linearity;Photoelectrons (n. p.e.);Peak Position (mV)");
    grL->SetMarkerStyle(21);
    grL->SetMarkerSize(1.2);
    TF1 *linFit = new TF1("linFit", "pol1", 0.5, nPeaks+0.5);
    grL->Fit(linFit, "RQ");
    grL->Draw("AP"); // AP shows axes and points

    double m_lin = linFit->GetParameter(1);
    double q_lin = linFit->GetParameter(0);
    double let_thr = (m_lin * 0.5) + q_lin;

    // =========================================================
    // 4. TIMING ANALYSIS (LET & CFD)
    // =========================================================
    TH2D *h2DLET = new TH2D("h2DLET", "LET Timing;Amplitude (mV);#Delta t (ns)", 150, 0, 50, 400, 43, 47);
    TH2D *h2DCFD = new TH2D("h2DCFD", "CFD Timing;Amplitude (mV);#Delta t (ns)", 150, 0, 50, 400, 43, 47);

    for(Long64_t i=0; i<treeCh1->GetEntries(); i++) {
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);
        auto c1 = Utils::correctWaveforms({{{t1, t1+1024}, {a1, a1+1024}}}, 30.0, true);
        auto cL = Utils::correctWaveforms({{{tL, tL+1024}, {aL, aL+1024}}}, 5., true);

        // Laser Trigger
        double tTrig = -999;
        const auto& alL = cL[0].second;
        for(size_t j=1; j<alL.size(); j++) {
            if(alL[j] > 50 && alL[j-1] <= 50) {
                tTrig = cL[0].first[j-1] + (50 - alL[j-1]) * (cL[0].first[j] - cL[0].first[j-1]) / (alL[j] - alL[j-1]);
                break;
            }
        }
        if(tTrig == -999) continue;

        const auto& as = c1[0].second;
        auto itMax = std::max_element(as.begin(), as.end());
        double maxA = *itMax;
        size_t idxMax = std::distance(as.begin(), itMax);

        // Fill 2D Maps
        double tLET = -999, tCFD = -999;
        for(size_t j=1; j<as.size(); j++) {
            if(tLET == -999 && as[j] > let_thr) 
                tLET = c1[0].first[j-1] + (let_thr - as[j-1]) * (c1[0].first[j]-c1[0].first[j-1])/(as[j]-as[j-1]);
            if(tCFD == -999 && j <= idxMax && as[j] > maxA*0.5)
                tCFD = c1[0].first[j-1] + (maxA*0.5 - as[j-1]) * (c1[0].first[j]-c1[0].first[j-1])/(as[j]-as[j-1]);
        }
        if(tLET != -999) h2DLET->Fill(maxA, tLET - tTrig);
        if(tCFD != -999) h2DCFD->Fill(maxA, tCFD - tTrig);
    }

    // =========================================================
    // 5. TRENDS AND SLICES (Corrected Ranges)
    // =========================================================
    auto process = [&](TH2D* h2D, const char* name, int color) {
        TCanvas *cSl = new TCanvas(Form("cSl_%s", name), Form("Slices %s", name), 1200, 800);
        cSl->Divide(4,2);
        std::vector<double> vp, vm, vme, vs, vse;

        for(int n=1; n<=nPeaks; n++) {
            cSl->cd(n);
            double low = q_lin + (n-0.4)*m_lin;
            double high = q_lin + (n+0.4)*m_lin;
            TH1D *hS = h2D->ProjectionY(Form("hS_%s_%d", name, n), h2D->GetXaxis()->FindBin(low), h2D->GetXaxis()->FindBin(high));
            hS->SetTitle(Form("%s Peak %d p.e.;#Delta t (ns);Counts", name, n));
            
            // ZOOM: Allargato l'intervallo visibile intorno al picco
            double central = hS->GetBinCenter(hS->GetMaximumBin());
            hS->GetXaxis()->SetRangeUser(central - 1.5, central + 1.5);

            TF1* fG = GaussStuff::FitQGauss(hS, Form("f%s_%d", name, n));
            if(fG) {
                vp.push_back(n); vm.push_back(fG->GetParameter(1)); vme.push_back(fG->GetParError(1));
                vs.push_back(fG->GetParameter(2)); vse.push_back(fG->GetParError(2));
                fG->SetLineColor(kGreen+2);
                fG->Draw("SAME");
            }
            hS->Draw("HIST SAME");
        }

        TGraphErrors *gR = new TGraphErrors(vp.size(), &vp[0], &vs[0], nullptr, &vse[0]);
        gR->SetMarkerStyle(20); gR->SetMarkerColor(color);
        gR->SetTitle(Form("%s Resolution Trend;Photoelectrons (n);#sigma_{#Delta t} (ns)", name));
        
        TGraphErrors *gM = new TGraphErrors(vp.size(), &vp[0], &vm[0], nullptr, &vme[0]);
        gM->SetMarkerStyle(24); gM->SetMarkerColor(color);
        gM->SetTitle(Form("%s Mean Trend;Photoelectrons (n);Mean #Delta t (ns)", name));

        return std::make_pair(gM, gR);
    };

    auto resLET = process(h2DLET, "LET", kBlue);
    auto resCFD = process(h2DCFD, "CFD", kRed);

    // Final Canvas (cComp.png)
    TCanvas *cFinal = new TCanvas("cFinal", "Comparison", 1100, 500);
    cFinal->Divide(2,1);
    cFinal->cd(1); gPad->SetGrid();
    resLET.first->Draw("AP"); resCFD.first->Draw("P SAME");
    TLegend *l1 = new TLegend(0.6, 0.7, 0.88, 0.88); l1->AddEntry(resLET.first, "LET Mean", "p"); l1->AddEntry(resCFD.first, "CFD Mean", "p"); l1->Draw();

    cFinal->cd(2); gPad->SetGrid();
    resLET.second->Draw("AP"); resCFD.second->Draw("P SAME");
    TF1 *fRes = new TF1("fRes", "sqrt(([0]*[0])/x + [1]*[1])", 0.5, 8.5);
    fRes->SetLineColor(kBlack); resCFD.second->Fit(fRes, "R");
    TLegend *l2 = new TLegend(0.6, 0.7, 0.88, 0.88); l2->AddEntry(resLET.second, "LET Res", "p"); l2->AddEntry(resCFD.second, "CFD Res", "p"); l2->Draw();

    gSystem->Exec("rm -f *.so *.pcm *.d");
}