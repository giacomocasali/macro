#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <iomanip>
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
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h> // Header aggiunto per correggere l'errore
#include <TGraph.h>

/**
 * Official code language: English [cite: 2026-02-26]
 * Scan automatico del LET con set di Canvas diagnostiche.
 */

// --- UTILS ---
static Double_t* setupBranch_dgz(TTree* tree, const char* branchName) {
    Double_t* arr = new Double_t[1024];
    tree->SetBranchAddress(branchName, arr);
    return arr;
}

static std::vector<std::pair<std::vector<double>, std::vector<double>>> correctWaveforms(
    const std::vector<std::pair<std::vector<double>, std::vector<double>>>& waveforms,
    double pre_signal_end = 30.0) {
    std::vector<std::pair<std::vector<double>, std::vector<double>>> corrected;
    for (const auto& entry : waveforms) {
        const std::vector<double>& time = entry.first;
        const std::vector<double>& amp = entry.second;
        double offset = 0.0;
        std::vector<double> pre_amp;
        for (size_t i = 0; i < (int)time.size(); ++i) if (time[i] <= pre_signal_end) pre_amp.push_back(amp[i]);
        if (!pre_amp.empty()) {
            std::vector<double> tmp = pre_amp;
            size_t mid = tmp.size() / 2;
            std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
            offset = tmp[mid];
        }
        std::vector<double> c_amp = amp;
        for (auto& a : c_amp) a -= offset;
        corrected.emplace_back(time, c_amp);
    }
    return corrected;
}

static Double_t qGaussAsym(Double_t* x, Double_t* par) {
    double xx = x[0], A = par[0], mean = par[1], sigma = par[2], q1 = par[3], q2 = par[4];
    if (sigma <= 0) return 0;
    if (xx <= mean) {
        double arg = 1.0 - (1.0 - q1) * (1.0 / (3.0 - q1)) * TMath::Power((xx - mean) / sigma, 2);
        return (arg <= 0) ? 0 : A * TMath::Power(arg, 1.0 / (1.0 - q1));
    } else {
        double arg = 1.0 - (1.0 - q2) * (1.0 / (3.0 - q2)) * TMath::Power((xx - mean) / sigma, 2);
        return (arg <= 0) ? 0 : A * TMath::Power(arg, 1.0 / (1.0 - q2));
    }
}

static std::tuple<double, double, double, double, TF1*> q_gauss_projection_zoom(TH1D* hProj, int pIdx) {
    if (!hProj || hProj->GetEntries() < 40) return {0,0,0,0, nullptr};
    int maxBin = hProj->GetMaximumBin();
    double center = hProj->GetBinCenter(maxBin);
    TF1* f = new TF1(Form("fZoom_%d", pIdx), qGaussAsym, center - 1.5, center + 1.5, 5);
    f->SetParameters(hProj->GetMaximum(), center, 0.15, 1.0, 1.2);
    f->FixParameter(3, 1.0);
    hProj->Fit(f, "RMQN");
    return {f->GetParameter(1), f->GetParError(1), std::abs(f->GetParameter(2)), f->GetParError(2), f};
}

// --- PROGRAMMA ---
void sipm_LET_scan() {
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);

    TFile *file = TFile::Open("../../data/data.vbias_{40}.root", "READ");
    if (!file) return;

    TTree *tCh1 = (TTree*)file->Get("ch1");
    TTree *tL = (TTree*)file->Get("laser");
    Double_t *time1 = setupBranch_dgz(tCh1, "time"), *amp1 = setupBranch_dgz(tCh1, "amplitude");
    Double_t *timeL = setupBranch_dgz(tL, "time"), *ampL = setupBranch_dgz(tL, "amplitude");

    // 1. CALIBRAZIONE
    TH1D *hCal = new TH1D("hCal", "Calibration;mV;Counts", 500, -5, 100);
    for (Long64_t i = 0; i < tCh1->GetEntries(); ++i) {
        tCh1->GetEntry(i);
        auto c = correctWaveforms({{{time1, time1+1024}, {amp1, amp1+1024}}});
        hCal->Fill(*std::max_element(c[0].second.begin(), c[0].second.end()));
    }
    TSpectrum s(10);
    s.Search(hCal, 2, "goff", 0.05);
    std::vector<double> pk(s.GetPositionX(), s.GetPositionX() + s.GetNPeaks());
    std::sort(pk.begin(), pk.end());
    TGraphErrors *grL = new TGraphErrors();
    for(int i=0; i<(int)pk.size() && i<7; i++) grL->SetPoint(i, i, pk[i]);
    TF1 *lin = new TF1("lin", "pol1", -0.5, 6.5);
    grL->Fit(lin, "RQ");
    double m_lin = lin->GetParameter(1), q_lin = lin->GetParameter(0);

    // 2. CICLO DI ANALISI LET
    std::vector<double> fracs = {0.2, 0.4, 0.6, 0.8};
    struct Res { double f, let, stoc, cst; };
    std::vector<Res> finalResults;

    for (double f : fracs) {
        double currentLET = q_lin + f * m_lin;
        TString folder = Form("scan_f%.1f", f);
        gSystem->Exec("mkdir -p " + folder);

        TH2D *h2D = new TH2D("h2D", Form("Map LET %.2f mV;Amp (mV);#Delta t (ns)", currentLET), 250, -5, 80, 400, -10, 80);
        TCanvas *cOvr = new TCanvas("cOvr", "Overlay", 800, 600);
        
        for (Long64_t i = 0; i < tCh1->GetEntries(); ++i) {
            tCh1->GetEntry(i); tL->GetEntry(i);
            auto c1 = correctWaveforms({{{time1, time1+1024}, {amp1, amp1+1024}}});
            auto cL = correctWaveforms({{{timeL, timeL+1024}, {ampL, ampL+1024}}});
            
            double tTrig = -1, tLET = -1;
            for(size_t j=1; j<1024; j++) if(cL[0].second[j]>10 && cL[0].second[j-1]<=10) {
                tTrig = cL[0].first[j-1] + (10-cL[0].second[j-1])*(cL[0].first[j]-cL[0].first[j-1])/(cL[0].second[j]-cL[0].second[j-1]);
                break;
            }
            if(tTrig < 0) continue;
            for(size_t j=1; j<1024; j++) if(c1[0].second[j]>currentLET && c1[0].second[j-1]<=currentLET) {
                tLET = c1[0].first[j-1] + (currentLET-c1[0].second[j-1])*(c1[0].first[j]-c1[0].first[j-1])/(c1[0].second[j]-c1[0].second[j-1]);
                break;
            }
            if(tLET > 0) h2D->Fill(*std::max_element(c1[0].second.begin(), c1[0].second.end()), tLET - tTrig);
            
            if(i < 15) {
                TGraph *gW = new TGraph(1024, &c1[0].first[0], &c1[0].second[0]);
                gW->SetLineColor(i+1); gW->Draw(i==0?"AL":"L");
                if(i==0) { gW->GetXaxis()->SetRangeUser(40, 100); gW->SetMaximum(60); }
            }
        }
        TLine *lLET = new TLine(40, currentLET, 100, currentLET); 
        lLET->SetLineColor(kRed); lLET->SetLineStyle(2); lLET->Draw();
        cOvr->SaveAs(folder + "/01_overlay.png");

        TCanvas *cSl = new TCanvas("cSl", "Slices", 1200, 800); cSl->Divide(3,2);
        TGraphErrors *grRes = new TGraphErrors();
        int pts = 0;
        for (int n = 1; n <= 6; ++n) {
            cSl->cd(n);
            TH1D *sl = h2D->ProjectionY(Form("sl_%d_f%.1f",n,f), h2D->GetXaxis()->FindBin(q_lin+(n-0.4)*m_lin), h2D->GetXaxis()->FindBin(q_lin+(n+0.4)*m_lin));
            auto fr = q_gauss_projection_zoom(sl, n);
            if (std::get<4>(fr)) {
                sl->Draw(); std::get<4>(fr)->Draw("SAME");
                grRes->SetPoint(pts, n, std::get<2>(fr));
                grRes->SetPointError(pts++, 0, std::get<3>(fr));
            }
        }
        cSl->SaveAs(folder + "/02_slices.png");

        TCanvas *cTrend = new TCanvas("cTrend", "Trend", 1200, 500); cTrend->Divide(2,1);
        cTrend->cd(1); h2D->Draw("COLZ");
        cTrend->cd(2);
        TF1 *fS = new TF1("fS", "sqrt(([0]*[0])/x + [1]*[1])", 0.5, 6.5);
        fS->SetParameters(0.4, 0.05); grRes->Fit(fS, "RQ");
        grRes->SetMarkerStyle(20); grRes->Draw("AP");
        cTrend->SaveAs(folder + "/03_final_trend.png");

        finalResults.push_back({f, currentLET, fS->GetParameter(0), fS->GetParameter(1)});
        delete h2D; delete cOvr; delete cSl; delete cTrend;
    }

    std::cout << "\n--- FINAL TABLE ---" << std::endl;
    printf("%-10s | %-10s | %-12s | %-12s\n", "Frac", "LET (mV)", "Stoc (ns)", "Const (ns)");
    for(auto &r : finalResults) printf("%.1f        | %6.2f     | %.4f       | %.4f\n", r.f, r.let, r.stoc, r.cst);
}