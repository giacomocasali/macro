
//questo file è quello iniziale, senza modifiche di Edo, che fa i grafici carini cercando intervalli tra 45 e 55 ns.

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
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
#include <TVirtualFitter.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TLegend.h>

// Official code language: English [cite: 2026-02-26]
// Self-contained: Utils and GaussStuff inlined below

// =========================================================
// FROM Utils.h — only what is needed
// =========================================================

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

        // Try time-based cut first; fall back to first 10% of samples if empty
        std::vector<double> pre_amp;
        for (size_t i = 0; i < time.size(); ++i)
            if (time[i] <= pre_signal_end) pre_amp.push_back(amp[i]);
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
            offset = tmp[mid];
        }
        std::vector<double> corrected_amp = amp;
        for (auto& a : corrected_amp) a -= offset;
        corrected.emplace_back(time, corrected_amp);
    }
    return corrected;
}

// =========================================================
// FROM gauss_stuff.h — only q_gauss_projection + qGaussAsym
// =========================================================

static Double_t qGaussAsym(Double_t* x, Double_t* par) {
    double xx    = x[0];
    double A     = par[0];
    double mean  = par[1];
    double sigma = par[2];
    double q1    = par[3];
    double q2    = par[4];
    if (sigma <= 0) return 0;
    const double eps = 1e-6;
    if (xx <= mean) {
        if (TMath::Abs(q1 - 1.0) < eps) {
            double arg_sq = TMath::Power((xx - mean) / sigma, 2);
            return A * TMath::Exp(-0.5 * arg_sq);
        } else {
            double arg = 1 - (1 - q1) * (1.0 / (3 - q1)) * TMath::Power((xx - mean) / sigma, 2);
            if (arg <= 0) return 0;
            return A * TMath::Power(arg, 1.0 / (1 - q1));
        }
    } else {
        double arg = 1 - (1 - q2) * (1.0 / (3 - q2)) * TMath::Power((xx - mean) / sigma, 2);
        if (arg <= 0) return 0;
        return A * TMath::Power(arg, 1.0 / (1 - q2));
    }
}

static std::tuple<double,double,double,double,TF1*> q_gauss_projection(TH1D* hProj) {
    if (!hProj || hProj->GetEntries() == 0) return {0,0,0,0,nullptr};
    double mean_guess  = hProj->GetMean();
    double sigma_guess = hProj->GetRMS();
    double A_guess     = hProj->GetMaximum();
    TF1* fFit = new TF1("qGaussFit", qGaussAsym,
                        hProj->GetXaxis()->GetXmin(),
                        hProj->GetXaxis()->GetXmax(), 5);
    fFit->SetParameters(A_guess, mean_guess, sigma_guess, 1.0, 1.3);
    fFit->SetParName(0, "A");
    fFit->SetParName(1, "#mu");
    fFit->SetParName(2, "#sigma");
    fFit->SetParName(3, "q_{1}");
    fFit->SetParName(4, "q_{2}");
    fFit->FixParameter(3, 1.0);
    fFit->SetNpx(1000);
    hProj->Fit(fFit, "IMRE");
    double mean      = fFit->GetParameter(1);
    double mean_err  = fFit->GetParError(1);
    double sigma     = fFit->GetParameter(2);
    double sigma_err = fFit->GetParError(2);
    return {mean, mean_err, sigma, sigma_err, fFit};
}

// =========================================================
// MAIN ANALYSIS FUNCTION
// =========================================================

Double_t multiGauss(Double_t *x, Double_t *par) {
    int n = (int)par[0];
    double sum = 0;
    for(int i=0; i<n; i++)
        sum += par[1+3*i] * TMath::Gaus(x[0], par[2+3*i], par[3+3*i]);
    return sum;
}

void sipm_timing_windowed() {
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TFile *file = TFile::Open("../../data/data.vbias_{40}.root", "READ");
    if(!file || file->IsZombie()) return;

    TTree *treeCh1   = (TTree*)file->Get("ch1");
    TTree *treeLaser = (TTree*)file->Get("laser");

    Double_t* t1 = setupBranch_dgz(treeCh1,   "time");
    Double_t* a1 = setupBranch_dgz(treeCh1,   "amplitude");
    Double_t* tL = setupBranch_dgz(treeLaser,  "time");
    Double_t* aL = setupBranch_dgz(treeLaser,  "amplitude");

    // DEBUG — remove after checking
    treeCh1->GetEntry(0); treeLaser->GetEntry(0);
    std::cout << "ch1   time range: " << t1[0] << " to " << t1[1023] << " ns" << std::endl;
    std::cout << "laser time range: " << tL[0] << " to " << tL[1023] << " ns" << std::endl;

    // =========================================================
    // PASS 1 — Build amplitude spectrum & fit to find LET
    // =========================================================
    TH1D *hTemp = new TH1D("hTemp", "Total SiPM Spectrum;amplitude (mV);counts", 300, -5, 60);

    for(Long64_t i=0; i<treeCh1->GetEntries(); i++) {
        if(i % 500 == 0 && gSystem->ProcessEvents()) break; // Esc
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);

        auto c1 = correctWaveforms({{{t1, t1+1024}, {a1, a1+1024}}}, 30.0);
        auto cL = correctWaveforms({{{tL, tL+1024}, {aL, aL+1024}}}, 30.0);

        double tTrig = -999;
        for(size_t j=1; j<cL[0].first.size(); j++) {
            if(cL[0].second[j] > 10 && cL[0].second[j-1] <= 10) {
                tTrig = cL[0].first[j-1] + (10 - cL[0].second[j-1]) * (cL[0].first[j] - cL[0].first[j-1]) / (cL[0].second[j] - cL[0].second[j-1]);
                break;
            }
        }
        if(tTrig == -999) continue;

        double maxA = -999;
        for(size_t j=0; j<c1[0].first.size(); j++) {
            double tRel = c1[0].first[j] - tTrig;
            if(tRel >= 45 && tRel <= 55)
                if(c1[0].second[j] > maxA) maxA = c1[0].second[j];
        }
        if(maxA > -900) hTemp->Fill(maxA);
    }

    TSpectrum s(20);
    int nFound = s.Search(hTemp, 1.5, "goff", 0.005);
    double *xPos = s.GetPositionX();
    std::vector<double> peaks(xPos, xPos + nFound);
    std::sort(peaks.begin(), peaks.end());

    int nPeaks = 8;
    TF1 *fitFunc = new TF1("fitFunc", multiGauss, -5, 55, 1+3*nPeaks);
    fitFunc->SetNpx(1000);
    fitFunc->FixParameter(0, nPeaks);
    for(int i=0; i<nPeaks; i++) {
        double x_est = (i < (int)peaks.size()) ? peaks[i] : (1.5 + 4.3 * i);
        fitFunc->SetParameter(1+3*i, hTemp->GetBinContent(hTemp->FindBin(x_est)));
        fitFunc->SetParameter(2+3*i, x_est);
        fitFunc->SetParLimits(2+3*i, x_est - 1.0, x_est + 1.0);
        fitFunc->SetParameter(3+3*i, 0.6);
        fitFunc->SetParLimits(3+3*i, 0.4, 1.1);
    }
    TVirtualFitter::SetDefaultFitter("Minuit2");
    hTemp->Fit(fitFunc, "RQML");

    // --- CANVAS 1: SPECTRUM ---
    TCanvas *cAmp = new TCanvas("cAmp", "Spectrum Analysis", 900, 900);
    cAmp->SetLeftMargin(0.15);
    cAmp->SetLogy();
    hTemp->SetLineColor(kBlack);
    hTemp->GetYaxis()->SetTitleOffset(1.6);
    hTemp->Draw("HIST");
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("SAME");
    for(int i=0; i<nPeaks; i++) {
        TF1 *gi = new TF1(Form("gi_%d", i), "gaus", -5, 60);
        gi->SetParameters(fitFunc->GetParameter(1+3*i), fitFunc->GetParameter(2+3*i), fitFunc->GetParameter(3+3*i));
        gi->SetLineColor(kAzure-4+i); gi->SetLineStyle(2); gi->Draw("SAME");
    }
    TPaveText *pt = new TPaveText(0.65, 0.68, 0.95, 0.92, "brNDC");
    pt->SetBorderSize(1); pt->SetFillColor(0);
    for(int i=0; i<nPeaks; i++) pt->AddText(Form("Peak %d: %.3f mV", i, fitFunc->GetParameter(2+3*i)));
    pt->Draw();

    // --- CANVAS 2: LINEARITY ---
    std::vector<double> xPE, yMean, yErrMean, ySigma, xErr;
    for(int i=0; i<nPeaks; i++) {
        xPE.push_back(i); xErr.push_back(0);
        yMean.push_back(fitFunc->GetParameter(2 + 3*i));
        yErrMean.push_back(fitFunc->GetParError(2 + 3*i));
        ySigma.push_back(fitFunc->GetParameter(3 + 3*i));
    }
    TCanvas *cLin = new TCanvas("cLin", "Linearity", 900, 900);
    cLin->SetGrid();
    cLin->SetLeftMargin(0.15);
    TGraphErrors *grSigma = new TGraphErrors(xPE.size(), &xPE[0], &yMean[0], &xErr[0], &ySigma[0]);
    grSigma->SetTitle("Linearity;photoelectrons (n);peak position (mV)");
    grSigma->GetYaxis()->SetTitleOffset(1.6);
    grSigma->SetFillColorAlpha(kAzure+1, 0.25); grSigma->Draw("A3");
    TGraphErrors *gr = new TGraphErrors(xPE.size(), &xPE[0], &yMean[0], &xErr[0], &yErrMean[0]);
    gr->SetMarkerStyle(21); gr->SetMarkerSize(1.2); gr->SetMarkerColor(kAzure+2);
    TF1 *linFit = new TF1("linFit", "pol1", -0.5, nPeaks - 0.5);
    linFit->SetLineColor(kRed+1); gr->Fit(linFit, "RQ"); gr->Draw("P SAME");
    TLegend *legLin = new TLegend(0.15, 0.75, 0.45, 0.88);
    legLin->SetBorderSize(1);
    legLin->AddEntry((TObject*)0, Form("m: %.3f #pm %.3f", linFit->GetParameter(1), linFit->GetParError(1)), "");
    legLin->AddEntry((TObject*)0, Form("q: %.3f #pm %.3f", linFit->GetParameter(0), linFit->GetParError(0)), "");
    legLin->Draw();

    double let = (linFit->GetParameter(1) / 2.0) + linFit->GetParameter(0);

    // =========================================================
    // PASS 2 — Timing, 2D map, and projection
    // =========================================================
    TH1D *hTime = new TH1D("hTime", Form("Time projection @ LET = %.2f mV;#Deltat (ns);counts", let), 500, 45, 50);
    TH2D *h2D = new TH2D("h2D", Form("Amp vs t_{LET}  @ LET = %.2f mV;amplitude (mV);t_{LET} - t_{laser} (ns)", let), 200, -5, 60, 200, 45, 55);

    for(Long64_t i=0; i<treeCh1->GetEntries(); i++) {
        if(i % 500 == 0 && gSystem->ProcessEvents()) break; // Esc
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);
        auto c1 = correctWaveforms({{{t1, t1+1024}, {a1, a1+1024}}}, 30.0);
        auto cL = correctWaveforms({{{tL, tL+1024}, {aL, aL+1024}}}, 30.0);

        double tTrig = -999;
        for(size_t j=1; j<cL[0].first.size(); j++) {
            if(cL[0].second[j] > 10 && cL[0].second[j-1] <= 10) {
                tTrig = cL[0].first[j-1] + (10 - cL[0].second[j-1]) * (cL[0].first[j] - cL[0].first[j-1]) / (cL[0].second[j] - cL[0].second[j-1]);
                break;
            }
        }
        if(tTrig == -999) continue;

        double maxA = -999, tLET = -999;
        const auto& times = c1[0].first;
        const auto& ampls = c1[0].second;

        for(size_t j=1; j<times.size(); j++) {
            if(ampls[j] > let && tLET == -999)
                tLET = times[j-1] + (let - ampls[j-1]) * (times[j] - times[j-1]) / (ampls[j] - ampls[j-1]);
            if(times[j]-tTrig >= 45 && times[j]-tTrig <= 55)
                if(ampls[j] > maxA) maxA = ampls[j];
        }

        if(tLET != -999 && maxA > -900) {
            hTime->Fill(tLET - tTrig);
            h2D->Fill(maxA, tLET - tTrig);
        }
    }

    // --- CANVAS 3: TIMING ---
    TCanvas *cTime = new TCanvas("cTime", "Timing Analysis", 900, 900);
    cTime->SetLeftMargin(0.15);
    hTime->GetYaxis()->SetTitleOffset(1.6);
    hTime->SetLineColor(kTeal+3); hTime->SetFillColorAlpha(kTeal, 0.3);
    hTime->Draw();

    // --- CANVAS 4: 2D MAP ---
    TCanvas *c2D = new TCanvas("c2D", "2D Map", 1000, 900);
    c2D->SetRightMargin(0.13);
    c2D->SetLeftMargin(0.15);
    h2D->SetContour(100);
    h2D->GetYaxis()->SetTitleOffset(1.6);
    h2D->Draw("COLZ");
    for(int n=0; n<nPeaks; n++) {
        double x_border = linFit->GetParameter(0) + (n + 0.5) * linFit->GetParameter(1);
        if(x_border < 60) {
            TLine *l = new TLine(x_border, 45, x_border, 50);
            l->SetLineColor(kRed); l->SetLineStyle(2); l->SetLineWidth(2);
            l->Draw("SAME");
        }
    }

    // --- CANVAS 5: AMPLITUDE PROJECTION ---
    TH1D *hProjX = (TH1D*)h2D->ProjectionX("hProjX");
    double x_min_proj = linFit->GetParameter(0) + 0.5 * linFit->GetParameter(1);
    hProjX->GetXaxis()->SetRangeUser(x_min_proj, 60.0);
    hProjX->SetTitle("Amplitude Projection;amplitude (mV);counts");
    TCanvas *cProj = new TCanvas("cProj", "Amplitude Projection", 900, 900);
    cProj->SetGrid();
    cProj->SetLeftMargin(0.15);
    hProjX->GetYaxis()->SetTitleOffset(1.6);
    hProjX->SetLineColor(kTeal+3); hProjX->SetFillColorAlpha(kTeal, 0.3);
    hProjX->Draw("HIST");

    // --- CANVAS 6: OVERLAY ---
    TCanvas *cOverlay = new TCanvas("cOverlay", "Spectrum Overlay Clean", 900, 900);
    cOverlay->SetLogy();
    cOverlay->SetGrid();
    cOverlay->SetLeftMargin(0.15);
    TH1D *hTempC6 = (TH1D*)hTemp->Clone("hTempC6");
    TH1D *hProjC6 = (TH1D*)hProjX->Clone("hProjC6");
    hTempC6->SetTitle("Total spectrum vs timing-selected projection;amplitude (mV);counts");
    hTempC6->GetYaxis()->SetTitleOffset(1.6);
    hTempC6->SetLineColor(kOrange+7);
    hTempC6->SetLineWidth(3);
    hTempC6->Draw("HIST");
    hProjC6->SetLineColor(kSpring-1);
    hProjC6->SetLineWidth(3);
    hProjC6->Draw("HIST SAME");
    TLegend *legOver = new TLegend(0.55, 0.70, 0.88, 0.85);
    legOver->AddEntry(hTempC6, "Total spectrum", "l");
    legOver->AddEntry(hProjC6, "Gated projection", "l");
    legOver->Draw();

    // --- CANVAS 7: TIMING SLICES ---
    TCanvas *cSlices = new TCanvas("cSlices", "Timing Slices per p.e.", 1200, 800);
    cSlices->Divide(4, 2);
    std::vector<double> vPE, vMeanT, vSigmaT, vEMean, vESigma;
    double m_lin = linFit->GetParameter(1);
    double q_lin = linFit->GetParameter(0);

    for(int n=1; n<nPeaks; n++) {
        cSlices->cd(n);
        gPad->SetLeftMargin(0.20);
        double low_edge  = q_lin + (n - 0.5) * m_lin;
        double high_edge = q_lin + (n + 0.5) * m_lin;
        int binLow  = h2D->GetXaxis()->FindBin(low_edge);
        int binHigh = h2D->GetXaxis()->FindBin(high_edge);
        TH1D *hSlice = h2D->ProjectionY(Form("hSlice_%d", n), binLow, binHigh);
        hSlice->SetTitle(Form("Peak %d p.e.;#Deltat (ns);counts", n));
        hSlice->GetYaxis()->SetTitleOffset(2.0);
        hSlice->SetLineColor(kAzure+7);
        hSlice->SetFillColorAlpha(kAzure+7, 0.15);
        hSlice->Draw("HIST");
        auto fitResult = q_gauss_projection(hSlice);
        if(std::get<4>(fitResult)) {
            std::get<4>(fitResult)->SetLineColor(kRed+1);
            std::get<4>(fitResult)->Draw("SAME");
            vPE.push_back(n);
            vMeanT.push_back(std::get<0>(fitResult));
            vEMean.push_back(std::get<1>(fitResult));
            vSigmaT.push_back(std::get<2>(fitResult));
            vESigma.push_back(std::get<3>(fitResult));
        }
    }

    // --- CANVAS 8: MEAN TREND ---
    TCanvas *cMeanTrend = new TCanvas("cMeanTrend", "Mean Timing Trend", 900, 900);
    cMeanTrend->SetGrid();
    cMeanTrend->SetLeftMargin(0.15);
    TGraphErrors *grMean = new TGraphErrors(vPE.size(), &vPE[0], &vMeanT[0], nullptr, &vEMean[0]);
    grMean->SetTitle("Mean #Deltat vs p.e.;photoelectrons (n);#mu_{#Deltat} (ns)");
    grMean->SetMarkerStyle(20);
    grMean->SetMarkerSize(1.2);
    grMean->SetMarkerColor(kAzure+2);
    grMean->GetYaxis()->SetTitleOffset(1.6);
    TF1 *fitMean = new TF1("fitMean", "[0] + [1]/sqrt(x)", 0.5, nPeaks - 0.5);
    fitMean->SetParameters(vMeanT.back(), 0.5);
    fitMean->SetLineColor(kOrange+7);
    grMean->Fit(fitMean, "RQ");
    grMean->Draw("AP");
    TLegend *legMean = new TLegend(0.45, 0.70, 0.88, 0.88);
    legMean->AddEntry(grMean, "Data Points", "pe");
    legMean->AddEntry(fitMean, "Fit: a + b/#sqrt{n}", "l");
    legMean->AddEntry((TObject*)0, Form("a: %.4f #pm %.4f", fitMean->GetParameter(0), fitMean->GetParError(0)), "");
    legMean->AddEntry((TObject*)0, Form("b: %.4f #pm %.4f", fitMean->GetParameter(1), fitMean->GetParError(1)), "");
    legMean->Draw();

    // --- CANVAS 9: RESOLUTION TREND ---
    TCanvas *cResTrend = new TCanvas("cResTrend", "Timing Resolution Trend", 900, 900);
    cResTrend->SetGrid();
    cResTrend->SetLeftMargin(0.15);
    TGraphErrors *grRes = new TGraphErrors(vPE.size(), &vPE[0], &vSigmaT[0], nullptr, &vESigma[0]);
    grRes->SetTitle("Timing resolution vs p.e.;photoelectrons (n);#sigma_{#Deltat} (ns)");
    grRes->SetMarkerStyle(21);
    grRes->SetMarkerSize(1.2);
    grRes->SetMarkerColor(kRed+1);
    grRes->GetYaxis()->SetTitleOffset(1.6);
    TF1 *fitSigma = new TF1("fitSigma", "sqrt(([0]*[0])/x + [1]*[1])", 0.5, nPeaks - 0.5);
    fitSigma->SetParNames("#sigma_{stoc}", "#sigma_{const}");
    fitSigma->SetParameters(0.5, 0.1);
    fitSigma->SetLineColor(kBlack);
    fitSigma->SetLineStyle(2);
    grRes->Fit(fitSigma, "RQ");
    grRes->Draw("AP");
    TLegend *legRes = new TLegend(0.45, 0.70, 0.88, 0.88);
    legRes->AddEntry(grRes, "Exp. Resolution", "pe");
    legRes->AddEntry(fitSigma, "Fit: #sqrt{#sigma_{stoc}^{2}/n + #sigma_{const}^{2}}", "l");
    legRes->AddEntry((TObject*)0, Form("#sigma_{stoc}: %.4f #pm %.4f", fitSigma->GetParameter(0), fitSigma->GetParError(0)), "");
    legRes->AddEntry((TObject*)0, Form("#sigma_{const}: %.4f #pm %.4f", fitSigma->GetParameter(1), fitSigma->GetParError(1)), "");
    legRes->Draw();

    std::cout << "\nAnalysis Complete. LET: " << let << " mV" << std::endl;

    // =========================================================
    // PASS 3 — CFD (Constant Fraction Discriminator) ANALYSIS
    // =========================================================
    TH1D *hTimeCFD = new TH1D("hTimeCFD", "Time projection (CFD 50%);#Deltat_{CFD} (ns);counts", 500, 45, 50);
    TH2D *h2DCFD   = new TH2D("h2DCFD", "Amp vs t_{CFD} (50%);amplitude (mV);t_{CFD} - t_{laser} (ns)", 200, -5, 60, 200, 45, 55);

    for(Long64_t i=0; i<treeCh1->GetEntries(); i++) {
        if(i % 500 == 0 && gSystem->ProcessEvents()) break; // Esc
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);

        auto c1 = correctWaveforms({{{t1, t1+1024}, {a1, a1+1024}}}, 30.0);
        auto cL = correctWaveforms({{{tL, tL+1024}, {aL, aL+1024}}}, 30.0);

        double tTrig = -999;
        for(size_t j=1; j<cL[0].first.size(); j++) {
            if(cL[0].second[j] > 10 && cL[0].second[j-1] <= 10) {
                tTrig = cL[0].first[j-1] + (10 - cL[0].second[j-1]) * (cL[0].first[j] - cL[0].first[j-1]) / (cL[0].second[j] - cL[0].second[j-1]);
                break;
            }
        }
        if(tTrig == -999) continue;

        double maxA = -999;
        int binMax = -1;
        const auto& times = c1[0].first;
        const auto& ampls = c1[0].second;

        for(size_t j=0; j<times.size(); j++) {
            if(times[j]-tTrig >= 45 && times[j]-tTrig <= 55) {
                if(ampls[j] > maxA) { maxA = ampls[j]; binMax = j; }
            }
        }

        if(maxA > 2.0 && binMax > 0) {
            double vCFD = maxA * 0.5;
            double tCFD = -999;
            for(int j=1; j <= binMax; j++) {
                if(ampls[j] >= vCFD && ampls[j-1] < vCFD) {
                    tCFD = times[j-1] + (vCFD - ampls[j-1]) * (times[j] - times[j-1]) / (ampls[j] - ampls[j-1]);
                    break;
                }
            }
            if(tCFD != -999) {
                hTimeCFD->Fill(tCFD - tTrig);
                h2DCFD->Fill(maxA, tCFD - tTrig);
            }
        }
    }

    // --- CANVAS 10: TIMING CFD ---
    TCanvas *cTimeCFD = new TCanvas("cTimeCFD", "Timing Analysis CFD", 900, 900);
    cTimeCFD->SetLeftMargin(0.15);
    hTimeCFD->GetYaxis()->SetTitleOffset(1.6);
    hTimeCFD->SetLineColor(kRed+1); hTimeCFD->SetFillColorAlpha(kRed+1, 0.3);
    hTimeCFD->Draw("HIST");

    // --- CANVAS 11: 2D MAP CFD ---
    TCanvas *c2DCFD = new TCanvas("c2DCFD", "2D Map CFD", 1000, 900);
    c2DCFD->SetRightMargin(0.13); c2DCFD->SetLeftMargin(0.15);
    h2DCFD->SetContour(100);
    h2DCFD->GetYaxis()->SetTitleOffset(1.6);
    h2DCFD->Draw("COLZ");

    // --- CANVAS 12: SLICES CFD ---
    TCanvas *cSlicesCFD = new TCanvas("cSlicesCFD", "Timing Slices CFD", 1200, 800);
    cSlicesCFD->Divide(4, 2);
    std::vector<double> vPE_C, vMeanT_C, vSigmaT_C, vEMean_C, vESigma_C;

    for(int n=1; n<nPeaks; n++) {
        cSlicesCFD->cd(n); gPad->SetLeftMargin(0.20);
        double low = q_lin + (n - 0.5) * m_lin;
        double high = q_lin + (n + 0.5) * m_lin;
        TH1D *hS = h2DCFD->ProjectionY(Form("hS_CFD_%d", n), h2DCFD->GetXaxis()->FindBin(low), h2DCFD->GetXaxis()->FindBin(high));
        hS->SetTitle(Form("Peak %d p.e. (CFD);#Deltat_{CFD} (ns);counts", n));
        hS->GetYaxis()->SetTitleOffset(2.0); hS->SetLineColor(kRed-7); hS->Draw("HIST");
        auto res = q_gauss_projection(hS);
        if(std::get<4>(res)) {
            std::get<4>(res)->SetLineColor(kBlack); std::get<4>(res)->Draw("SAME");
            vPE_C.push_back(n); vMeanT_C.push_back(std::get<0>(res));
            vEMean_C.push_back(std::get<1>(res)); vSigmaT_C.push_back(std::get<2>(res));
            vESigma_C.push_back(std::get<3>(res));
        }
    }

    // --- CANVAS 13: MEAN TREND CFD ---
    TCanvas *cMeanCFD = new TCanvas("cMeanCFD", "Mean Trend CFD", 900, 900);
    cMeanCFD->SetGrid(); cMeanCFD->SetLeftMargin(0.15);
    TGraphErrors *grMeanC = new TGraphErrors(vPE_C.size(), &vPE_C[0], &vMeanT_C[0], nullptr, &vEMean_C[0]);
    grMeanC->SetTitle("Mean #Deltat_{CFD} vs p.e.;photoelectrons (n);#mu_{CFD} (ns)");
    grMeanC->SetMarkerStyle(20); grMeanC->SetMarkerColor(kRed+2);
    grMeanC->GetYaxis()->SetTitleOffset(1.6);
    TF1 *fMC = new TF1("fMC", "[0] + [1]/sqrt(x)", 0.5, nPeaks-0.5);
    grMeanC->Fit(fMC, "RQ"); grMeanC->Draw("AP");

    // --- CANVAS 14: RESOLUTION TREND CFD ---
    TCanvas *cResCFD = new TCanvas("cResCFD", "Resolution Trend CFD", 900, 900);
    cResCFD->SetGrid(); cResCFD->SetLeftMargin(0.15);
    TGraphErrors *grResC = new TGraphErrors(vPE_C.size(), &vPE_C[0], &vSigmaT_C[0], nullptr, &vESigma_C[0]);
    grResC->SetTitle("Resolution CFD vs p.e.;photoelectrons (n);#sigma_{CFD} (ns)");
    grResC->SetMarkerStyle(21); grResC->SetMarkerColor(kRed+2);
    grResC->GetYaxis()->SetTitleOffset(1.6);
    TF1 *fSC = new TF1("fSC", "sqrt(([0]*[0])/x + [1]*[1])", 0.5, nPeaks-0.5);
    fSC->SetLineColor(kRed); grResC->Fit(fSC, "RQ"); grResC->Draw("AP");

    std::cout << "CFD Analysis finished." << std::endl;

    // --- SAVE ALL CANVASES AS PNG ---
    cAmp->SaveAs("plot_01_spectrum.png");
    cLin->SaveAs("plot_02_linearity.png");
    cTime->SaveAs("plot_03_timing_LET.png");
    c2D->SaveAs("plot_04_2Dmap_LET.png");
    cProj->SaveAs("plot_05_amplitude_projection.png");
    cOverlay->SaveAs("plot_06_overlay.png");
    cSlices->SaveAs("plot_07_timing_slices.png");
    cMeanTrend->SaveAs("plot_08_mean_trend.png");
    cResTrend->SaveAs("plot_09_resolution_trend.png");
    cTimeCFD->SaveAs("plot_10_timing_CFD.png");
    c2DCFD->SaveAs("plot_11_2Dmap_CFD.png");
    cSlicesCFD->SaveAs("plot_12_slices_CFD.png");
    cMeanCFD->SaveAs("plot_13_mean_trend_CFD.png");
    cResCFD->SaveAs("plot_14_resolution_CFD.png");
    std::cout << "All canvases saved as PNG." << std::endl;

    gSystem->Exec("find . -maxdepth 1 -type f ! -name '*.root' ! -name '*.cpp' ! -name '*.h' ! -name '*.png' -delete");
}
