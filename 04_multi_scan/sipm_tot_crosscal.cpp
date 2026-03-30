// sipm_tot_crosscal.cpp
// Analisi TOT con calibrazione incrociata (Approccio 3).
//
// La soglia LET è ricavata calibrando lo spettro p.e. su un file
// con illuminazione sufficiente (es. filter=3 o filter=4), e poi
// applicando la stessa soglia in mV ai file di analisi scelti.
// Questo è fisicamente valido perché il guadagno m [mV/p.e.] dipende
// solo da V_bias, non dall'intensità della luce.
//
//   ampiezza = q + n * m   (fit lineare picchi p.e. — file di calibrazione)
//   let_thr  = q + frac_pe * m   [mV]  — applicato a tutti i file scelti
//
// Header:
//   header/Config.h           — parametri, FileInfo, scansione file
//   header/SignalProcessing.h — baseline, filtro, trigger laser
//   header/Plotting.h         — canvas (incluso per coerenza)

#include "../header/Config.h"
#include "../header/SignalProcessing.h"
#include "../header/Plotting.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>

// ════════════════════════════════════════════════════════════
//  SELEZIONE INTERATTIVA DEI FILE
//  Mostra la lista dei file trovati e chiede all'utente quali
//  usare per calibrazione e quali per l'analisi TOT.
// ════════════════════════════════════════════════════════════
static void printFileList(const std::vector<FileInfo>& files) {
    std::cout << "\nFiles found:\n";
    for (int i = 0; i < (int)files.size(); ++i)
        std::cout << "  [" << i << "] " << files[i].tag << "\n";
}

static int selectOneFile(const std::vector<FileInfo>& files,
                         const std::string& prompt) {
    int idx = -1;
    while (idx < 0 || idx >= (int)files.size()) {
        std::cout << prompt;
        std::cin >> idx;
        if (idx < 0 || idx >= (int)files.size())
            std::cout << "  [ERROR] Invalid index. Enter 0-"
                      << files.size()-1 << ".\n";
    }
    return idx;
}

static std::vector<int> selectMultipleFiles(const std::vector<FileInfo>& files,
                                             const std::string& prompt) {
    std::vector<int> selected;
    std::string line;
    std::cin.ignore();
    while (selected.empty()) {
        std::cout << prompt;
        std::getline(std::cin, line);
        // "all" seleziona tutti
        if (line == "all" || line == "ALL") {
            for (int i = 0; i < (int)files.size(); ++i) selected.push_back(i);
            break;
        }
        // parsing lista separata da virgole: "0,2,3"
        std::stringstream ss(line);
        std::string tok;
        bool ok = true;
        while (std::getline(ss, tok, ',')) {
            try {
                int i = std::stoi(tok);
                if (i < 0 || i >= (int)files.size()) { ok = false; break; }
                selected.push_back(i);
            } catch (...) { ok = false; break; }
        }
        if (!ok || selected.empty()) {
            selected.clear();
            std::cout << "  [ERROR] Enter indices like \"2\" or \"0,2\" or \"all\".\n";
        }
    }
    return selected;
}

// ════════════════════════════════════════════════════════════
//  CALIBRAZIONE SPETTRO p.e. (su file di calibrazione)
// ════════════════════════════════════════════════════════════
static bool calibrateSpectrum(TTree* treeCh1,
                               double cutoff_MHz,
                               double fs_MHz,
                               int j_trig_start, int j_trig_end,
                               const std::string& tag,
                               double& m_cal, double& q_cal) {
    const int N = 1024;
    Double_t t1[N], a1[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);

    TH1D* hSpec = new TH1D(Form("hSpec_%s", tag.c_str()),
        Form("p.e. spectrum   %s;Amplitude (mV);Counts", tag.c_str()),
        500, -5.0, 120.0);
    hSpec->SetDirectory(nullptr);

    Long64_t nEntries = treeCh1->GetEntries();
    std::cout << "  [CAL] Building amplitude spectrum from "
              << tag << " (" << nEntries << " events)...\n";

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 1000 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  [CAL] " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, BASELINE_START, BASELINE_END),
            cutoff_MHz, fs_MHz);
        // Massimo sull'intera waveform — stesso approccio delle macro originali.
        // Non limitato alla finestra trigger per avere lo spettro p.e. completo.
        double ampMax = *std::max_element(af.begin(), af.end());
        if (ampMax > -1e8) hSpec->Fill(ampMax);
    }
    std::cout << "\r  [CAL] " << nEntries << " / " << nEntries << " — done.\n";

    // Soglia 0.05 come nelle macro originali — robusta e stabile
    TSpectrum sp(15);
    int nPeaks = sp.Search(hSpec, 2, "goff", 0.05);
    if (nPeaks < 2) {
        std::cerr << "  [CAL] Not enough peaks (" << nPeaks
                  << ") in calibration file — try a brighter file.\n";
        return false;
    }

    std::vector<double> pkPos(sp.GetPositionX(), sp.GetPositionX() + nPeaks);
    std::sort(pkPos.begin(), pkPos.end());
    int nFit = std::min((int)pkPos.size(), 7);

    TGraphErrors* grCal = new TGraphErrors();
    for (int i = 0; i < nFit; ++i)
        grCal->SetPoint(i, i, pkPos[i]);

    TF1* fLin = new TF1(Form("fLin_%s", tag.c_str()), "pol1", -0.5, nFit-0.5);
    grCal->Fit(fLin, "RQ");
    q_cal = fLin->GetParameter(0);
    m_cal = fLin->GetParameter(1);

    std::cout << "  [CAL] Gain = " << m_cal << " mV/p.e."
              << "   Offset = " << q_cal << " mV\n";

    // Canvas di calibrazione
    TCanvas* cCal = new TCanvas(Form("cCal_%s", tag.c_str()),
        Form("Calibration — %s", tag.c_str()), 1200, 500);
    cCal->Divide(2, 1);

    cCal->cd(1); gPad->SetLogy(); gPad->SetGrid();
    hSpec->SetLineColor(kAzure+1); hSpec->SetLineWidth(2);
    hSpec->Draw("HIST");
    for (int i = 0; i < nFit; ++i) {
        TLine* lp = new TLine(pkPos[i], 0, pkPos[i], hSpec->GetMaximum());
        lp->SetLineColor(kRed+1); lp->SetLineStyle(2); lp->SetLineWidth(2);
        lp->Draw("same");
    }

    cCal->cd(2); gPad->SetGrid();
    grCal->SetMarkerStyle(20); grCal->SetMarkerColor(kAzure+1);
    grCal->SetTitle(Form("Linearity   %s;p.e. number;Peak amplitude (mV)", tag.c_str()));
    grCal->Draw("AP");
    fLin->SetLineColor(kRed+1); fLin->Draw("same");

    TPaveText* pt = new TPaveText(0.15, 0.75, 0.55, 0.88, "NDC");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
    pt->SetTextFont(42); pt->SetTextSize(0.038);
    pt->AddText(Form("Gain   = %.3f mV/p.e.", m_cal));
    pt->AddText(Form("Offset = %.3f mV",      q_cal));
    pt->Draw();

    cCal->Update(); cCal->Modified();
    cCal->SaveAs(Form("tot_calibration_%s.png", tag.c_str()));

    return true;
}

// ════════════════════════════════════════════════════════════
//  CALCOLO TOT con interpolazione lineare
// ════════════════════════════════════════════════════════════
static std::pair<double,double> computeTOT(const std::vector<double>& time,
                                           const std::vector<double>& amp,
                                           double threshold,
                                           int j_start, int j_end) {
    double t_rise = -1.0, t_fall = -1.0;
    for (int j = j_start + 1; j <= j_end; ++j) {
        if (amp[j-1] < threshold && amp[j] >= threshold) {
            t_rise = time[j-1] + (threshold - amp[j-1])
                     * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
            break;
        }
    }
    if (t_rise < 0) return {-1.0, -1.0};
    for (int j = j_start + 1; j <= j_end; ++j) {
        if (time[j] <= t_rise) continue;
        if (amp[j-1] >= threshold && amp[j] < threshold) {
            t_fall = time[j-1] + (threshold - amp[j-1])
                     * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
            break;
        }
    }
    if (t_fall < 0) return {-1.0, -1.0};
    return {t_rise, t_fall};
}

// ════════════════════════════════════════════════════════════
//  FIT Q-GAUSS ASIMMETRICO sulle slice
// ════════════════════════════════════════════════════════════
static Double_t qGaussAsym(Double_t* x, Double_t* par) {
    double xx = x[0], A = par[0], mean = par[1], sigma = par[2];
    double q1 = par[3], q2 = par[4];
    if (sigma <= 0) return 0;
    if (xx <= mean) {
        double arg = 1.0-(1.0-q1)*(1.0/(3.0-q1))*std::pow((xx-mean)/sigma,2);
        return (arg<=0) ? 0 : A*std::pow(arg,1.0/(1.0-q1));
    } else {
        double arg = 1.0-(1.0-q2)*(1.0/(3.0-q2))*std::pow((xx-mean)/sigma,2);
        return (arg<=0) ? 0 : A*std::pow(arg,1.0/(1.0-q2));
    }
}

static std::tuple<double,double,double,double,TF1*>
fitSlice(TH1D* h, int idx, const std::string& tag,
         double fit_lo, double fit_hi) {
    // Restringe lo spettro alla finestra di fit specificata dall'operatore
    // prima di cercare il massimo — evita di fittare rumore fuori finestra
    if (!h || h->GetEntries() < 10) return {0,0,0,0,nullptr};

    // Conta eventi nella finestra di fit
    int b1 = h->FindBin(fit_lo);
    int b2 = h->FindBin(fit_hi);
    double entries_in_window = h->Integral(b1, b2);
    if (entries_in_window < 10) return {0,0,0,0,nullptr};

    // Trova il massimo nella finestra di fit
    int maxBin = b1;
    for (int b = b1; b <= b2; ++b)
        if (h->GetBinContent(b) > h->GetBinContent(maxBin)) maxBin = b;
    double center = h->GetBinCenter(maxBin);

    // Stima sigma dal RMS nella finestra
    double sumW=0, sumWX=0, sumWX2=0;
    for (int b = b1; b <= b2; ++b) {
        double c = h->GetBinContent(b);
        double x = h->GetBinCenter(b);
        sumW += c; sumWX += c*x; sumWX2 += c*x*x;
    }
    double sigEst = 0.5;
    if (sumW > 0) {
        double mean = sumWX/sumW;
        double var  = sumWX2/sumW - mean*mean;
        sigEst = std::max(std::sqrt(std::max(var, 0.0)), 0.3);
    }

    std::string fname = Form("fSlice_%s_%d", tag.c_str(), idx);
    TF1* f = new TF1(fname.c_str(), qGaussAsym,
                     std::max(fit_lo, center-3.0*sigEst),
                     std::min(fit_hi, center+3.0*sigEst), 5);
    f->SetParameters(h->GetBinContent(maxBin), center, sigEst, 1.0, 1.2);
    f->FixParameter(3, 1.0);
    h->Fit(f, "RMQN");
    return {f->GetParameter(1), f->GetParError(1),
            std::abs(f->GetParameter(2)), f->GetParError(2), f};
}

// ════════════════════════════════════════════════════════════
//  PROCESS A SINGLE FILE (con soglia let_thr già nota)
// ════════════════════════════════════════════════════════════
static void processFile(const FileInfo& info,
                        double cutoff_MHz,
                        double t_trig_start,
                        double t_trig_end,
                        double fit_lo,
                        double fit_hi,
                        double let_thr,
                        double frac_pe) {

    TFile* file = TFile::Open(info.path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[SKIP] Cannot open: " << info.path << "\n"; return;
    }
    TTree* treeCh1   = (TTree*)file->Get("ch1");
    TTree* treeLaser = (TTree*)file->Get("laser");
    if (!treeCh1 || !treeLaser) {
        std::cerr << "[SKIP] Missing trees in: " << info.path << "\n";
        file->Close(); return;
    }

    const int N = 1024;
    Double_t t1[N], a1[N], tL[N], aL[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);
    treeLaser->SetBranchAddress("time",      tL);
    treeLaser->SetBranchAddress("amplitude", aL);
    treeCh1->GetEntry(0);
    double fs_MHz = 1000.0 / (t1[1] - t1[0]);

    int j_trig_start, j_trig_end;
    triggerWindowIndices(t1, N, t_trig_start, t_trig_end,
                         j_trig_start, j_trig_end, info.tag);

    std::cout << "\n[" << info.tag << "]  "
              << treeCh1->GetEntries() << " events,  fs=" << fs_MHz << " MHz\n"
              << "  Trigger window : [" << t1[j_trig_start] << ", "
              << t1[j_trig_end] << "] ns\n"
              << "  LET threshold  : " << let_thr << " mV"
              << "  (= " << frac_pe << " p.e., cross-calibrated)\n";

    TH2D* h2D = new TH2D(
        Form("h2D_%s", info.tag.c_str()),
        Form("TOT vs #Delta t   V_{bias}=%.0f V, filter=%.0f, LET=%.2f p.e."
             ";TOT (ns);#Delta t = t_{LET} - t_{laser} (ns)",
             info.vbias, info.filter, frac_pe),
        200, 0.0, 50.0, 400, -30.0, 100.0);
    h2D->SetDirectory(nullptr);

    Long64_t nEntries = treeCh1->GetEntries();
    long nNoLaser = 0, nNoTOT = 0, nFilled = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  [TOT] " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        treeLaser->GetEntry(i);

        double t_laser = laserTriggerTime(tL, aL, N, 10.0);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, BASELINE_START, BASELINE_END),
            cutoff_MHz, fs_MHz);

        auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                                            j_trig_start, j_trig_end);
        if (t_rise < 0 || t_fall < 0) { ++nNoTOT; continue; }

        double tot     = t_fall - t_rise;
        double delta_t = t_rise - t_laser;
        if (tot > 0 && tot < 100.0) { h2D->Fill(tot, delta_t); ++nFilled; }
    }
    std::cout << "\r  [TOT] " << nEntries << " / " << nEntries << " — done.\n"
              << "  Filled: " << nFilled
              << "  |  No laser: " << nNoLaser
              << "  |  No TOT: "   << nNoTOT << "\n";

    // Canvas mappa 2D
    TCanvas* c2D = new TCanvas(Form("c2D_%s", info.tag.c_str()),
        Form("TOT map — %s", info.tag.c_str()), 900, 700);
    c2D->SetRightMargin(0.15); c2D->SetLeftMargin(PAD_LEFT);
    c2D->SetBottomMargin(PAD_BOTTOM); c2D->SetTopMargin(PAD_TOP);
    c2D->SetLogz(); c2D->SetGrid();
    h2D->Draw("COLZ");
    c2D->Update(); c2D->Modified();
    c2D->SaveAs(Form("tot_map_%s.png", info.tag.c_str()));

    std::cout << "  Fit window     : [" << fit_lo << ", " << fit_hi << "] ns\n";

    // Slice
    const int N_SLICES = 6;
    double tot_max = 0.0;
    for (int bx = h2D->GetNbinsX(); bx >= 1; --bx) {
        double s = 0;
        for (int by = 1; by <= h2D->GetNbinsY(); ++by)
            s += h2D->GetBinContent(bx, by);
        if (s > 10) { tot_max = h2D->GetXaxis()->GetBinCenter(bx); break; }
    }
    if (tot_max <= 0) tot_max = 20.0;
    double slice_w = tot_max / N_SLICES;

    TCanvas* cSl = new TCanvas(Form("cSl_%s", info.tag.c_str()),
        Form("TOT slices — %s", info.tag.c_str()), 1200, 800);
    cSl->Divide(3, 2);

    TGraphErrors* grRes  = new TGraphErrors();
    TGraphErrors* grMean = new TGraphErrors();
    int pts = 0;

    for (int s = 0; s < N_SLICES; ++s) {
        double xlo_s = s * slice_w;
        double xhi_s = xlo_s + slice_w;
        int bx1 = h2D->GetXaxis()->FindBin(xlo_s);
        int bx2 = h2D->GetXaxis()->FindBin(xhi_s) - 1;
        std::string hname = Form("hsl_%s_%d", info.tag.c_str(), s);
        TH1D* hsl = h2D->ProjectionY(hname.c_str(), bx1, bx2);
        hsl->SetDirectory(nullptr);
        hsl->SetTitle(Form("TOT #in [%.1f, %.1f] ns;#Delta t (ns);Counts",
                           xlo_s, xhi_s));
        hsl->GetXaxis()->UnZoom(); hsl->GetYaxis()->UnZoom();
        cSl->cd(s+1);
        auto [mean, meanErr, sigma, sigmaErr, fFit] =
            fitSlice(hsl, s, info.tag, fit_lo, fit_hi);
        if (fFit) {
            hsl->Draw("HIST");
            fFit->SetLineColor(kGreen+2); fFit->Draw("SAME");
            double tc = 0.5*(xlo_s+xhi_s);
            grRes->SetPoint(pts, tc, sigma);   grRes->SetPointError(pts, 0, sigmaErr);
            grMean->SetPoint(pts, tc, mean);   grMean->SetPointError(pts, 0, meanErr);
            ++pts;
        }
    }
    cSl->Update(); cSl->Modified();
    cSl->SaveAs(Form("tot_slices_%s.png", info.tag.c_str()));

    // Trend
    if (pts > 1) {
        TCanvas* cTrend = new TCanvas(Form("cTrend_%s", info.tag.c_str()),
            Form("TOT trend — %s", info.tag.c_str()), 1200, 500);
        cTrend->Divide(2, 1);
        cTrend->cd(1); gPad->SetGrid();
        grRes->SetMarkerStyle(20); grRes->SetMarkerColor(kAzure+1);
        grRes->SetLineColor(kAzure+1);
        grRes->SetTitle(Form("#sigma_{#Delta t} vs TOT   LET=%.2f p.e."
            ";TOT (ns);#sigma_{#Delta t} (ns)", frac_pe));
        grRes->Draw("AP");
        if (pts >= 3) {
            TF1* fS = new TF1(Form("fS_%s", info.tag.c_str()),
                "sqrt([0]*[0]/x+[1]*[1])",
                grRes->GetX()[0]*0.8, grRes->GetX()[pts-1]*1.2);
            fS->SetParameters(1.0, 0.1); fS->SetLineColor(kRed+1);
            grRes->Fit(fS, "RQ");
            std::cout << "  Stochastic fit: S=" << fS->GetParameter(0)
                      << "  C=" << fS->GetParameter(1) << " ns\n";
        }
        cTrend->cd(2); gPad->SetGrid();
        grMean->SetMarkerStyle(20); grMean->SetMarkerColor(kOrange+7);
        grMean->SetLineColor(kOrange+7);
        grMean->SetTitle(Form("Mean #Delta t vs TOT   LET=%.2f p.e."
            ";TOT (ns);Mean #Delta t (ns)", frac_pe));
        grMean->Draw("AP");
        cTrend->Update(); cTrend->Modified();
        cTrend->SaveAs(Form("tot_trend_%s.png", info.tag.c_str()));
    }

    // Salva ROOT
    std::string outname = Form("tot_%s_let%.2fpe.root", info.tag.c_str(), frac_pe);
    TFile* fout = new TFile(outname.c_str(), "RECREATE");
    h2D->Write();
    if (pts > 0) { grRes->Write("grSigma"); grMean->Write("grMean"); }
    fout->Write(); fout->Close();
    std::cout << "  Saved: " << outname << "\n";

    file->Close();
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════
void sipm_tot_crosscal() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::vector<FileInfo> files = findInputFiles();
    if (files.empty()) {
        std::cerr << "No data.vbias_{V}_L.root files found.\n"; return;
    }

    // Mostra lista e chiede selezione
    printFileList(files);

    int calIdx = selectOneFile(files,
        "\nSelect calibration file index (brightest, e.g. filter=3): ");

    std::vector<int> anaIdx = selectMultipleFiles(files,
        "Select analysis file indices (e.g. \"2\" or \"0,2\" or \"all\"): ");

    // Parametri
    double cutoff_MHz;
    std::cout << "\nLow-pass filter cut-off [MHz]: ";
    std::cin >> cutoff_MHz;

    // Asse temporale
    double t_wf_min = 0.0, t_wf_max = 0.0;
    {
        TFile* f0 = TFile::Open(files[calIdx].path.c_str(), "READ");
        if (f0 && !f0->IsZombie()) {
            TTree* tr = (TTree*)f0->Get("ch1");
            if (tr) {
                const int N0 = 1024; Double_t t0[N0];
                tr->SetBranchAddress("time", t0); tr->GetEntry(0);
                t_wf_min = t0[0]; t_wf_max = t0[N0-1];
            }
            f0->Close();
        }
    }

    double t_trig_start, t_trig_end;
    std::cout << "\n--- Trigger search window ---\n"
              << "  Waveform range  : [" << t_wf_min << ", " << t_wf_max << "] ns\n"
              << "  Baseline (fixed): [" << BASELINE_START << ", "
                                         << BASELINE_END << ") ns\n"
              << "  Start [ns]: ";
    std::cin >> t_trig_start;
    std::cout << "  End   [ns]: ";
    std::cin >> t_trig_end;
    if (t_trig_start >= t_trig_end) {
        std::cerr << "[ERROR] Invalid trigger window.\n"; return;
    }

    double frac_pe;
    std::cout << "\nLET threshold in p.e. units (e.g. 0.5, 1.0, 1.5): ";
    std::cin >> frac_pe;
    while (frac_pe <= 0.0) {
        std::cout << "  [ERROR] Enter a positive value: ";
        std::cin >> frac_pe;
    }

    // Finestra per il fit delle slice — stretta intorno al picco del segnale
    // (diversa dalla finestra trigger che può essere larga)
    double fit_lo, fit_hi;
    std::cout << "\n--- Fit window for slice analysis ---\n"
              << "  Enter the Δt range [ns] where the signal peak is expected.\n"
              << "  (Use laser_timing canvas to identify the peak position)\n"
              << "  Start [ns]: ";
    std::cin >> fit_lo;
    std::cout << "  End   [ns]: ";
    std::cin >> fit_hi;
    if (fit_lo >= fit_hi) {
        std::cerr << "[ERROR] Invalid fit window — using trigger window.\n";
        fit_lo = t_trig_start; fit_hi = t_trig_end;
    }

    // Calibrazione sul file scelto
    TFile* fCal = TFile::Open(files[calIdx].path.c_str(), "READ");
    if (!fCal || fCal->IsZombie()) {
        std::cerr << "[ERROR] Cannot open calibration file.\n"; return;
    }
    TTree* tCal = (TTree*)fCal->Get("ch1");

    // Indici trigger per la calibrazione
    const int NC = 1024; Double_t tc[NC];
    tCal->SetBranchAddress("time", tc); tCal->GetEntry(0);
    double fs_cal = 1000.0 / (tc[1] - tc[0]);
    int jc_start, jc_end;
    triggerWindowIndices(tc, NC, t_trig_start, t_trig_end,
                         jc_start, jc_end, files[calIdx].tag);

    double m_cal = 0, q_cal = 0;
    if (!calibrateSpectrum(tCal, cutoff_MHz, fs_cal,
                           jc_start, jc_end,
                           files[calIdx].tag, m_cal, q_cal)) {
        fCal->Close(); return;
    }
    fCal->Close();

    double let_thr = q_cal + frac_pe * m_cal;
    std::cout << "\nCalibration result:\n"
              << "  Gain   = " << m_cal  << " mV/p.e.\n"
              << "  Offset = " << q_cal  << " mV\n"
              << "  LET threshold = " << let_thr << " mV  (= "
              << frac_pe << " p.e.)\n"
              << "  Fit window    = [" << fit_lo << ", " << fit_hi << "] ns\n\n";

    // Analisi sui file selezionati
    for (int idx : anaIdx)
        processFile(files[idx], cutoff_MHz, t_trig_start, t_trig_end,
                    fit_lo, fit_hi, let_thr, frac_pe);

    std::cout << "\n=== All files processed. ===\n";
}
