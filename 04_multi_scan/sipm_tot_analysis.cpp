// sipm_tot_analysis.cpp
// Analisi TOT (Time Over Threshold) vs risoluzione temporale.
//
// La soglia LET è ricavata dalla calibrazione dello spettro dei fotoelettroni:
//   ampiezza = q + n * m   (fit lineare dei picchi p.e.)
//   let_thr  = q + frac_pe * m   [mV]  — soglia fissa per tutti gli eventi
//
// Per ogni evento:
//   1. Corregge la baseline e applica il filtro passa basso
//   2. Calcola il TOT come t_discesa - t_salita sulla soglia let_thr
//   3. Misura il tempo di arrivo LET rispetto al trigger laser
//   4. Produce TH2D: TOT (X) vs Delta_t (Y)
//   5. Slice del TH2D con fit q-gaussiano per estrarre sigma(Delta_t) vs TOT
//
// Header:
//   header/Config.h           — parametri, FileInfo, scansione file
//   header/SignalProcessing.h — baseline, filtro, trigger laser
//   header/Plotting.h         — canvas overlay (incluso per coerenza)

#include "../header/Config.h"
#include "../header/SignalProcessing.h"
#include "../header/Plotting.h"

#include <iostream>
#include <vector>
#include <string>
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
//  CALIBRAZIONE SPETTRO p.e.
//  Scansiona tutti gli eventi, riempie lo spettro delle ampiezze
//  massime nella finestra trigger, trova i picchi con TSpectrum,
//  fitta la retta  ampiezza = q + n * m  e restituisce m e q.
//  Produce e salva un canvas di calibrazione.
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
    std::cout << "  [CAL] Building amplitude spectrum (" << nEntries << " events)...\n";

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
        double ampMax = -1e9;
        for (int j = j_trig_start; j <= j_trig_end; ++j)
            if (af[j] > ampMax) ampMax = af[j];
        if (ampMax > -1e8) hSpec->Fill(ampMax);
    }
    std::cout << "\r  [CAL] " << nEntries << " / " << nEntries << " — done.\n";

    // Ricerca picchi con TSpectrum
    TSpectrum sp(15);
    int nPeaks = sp.Search(hSpec, 2, "goff", 0.005);
    if (nPeaks < 2) {
        std::cerr << "  [CAL] Not enough peaks found (" << nPeaks
                  << ") — cannot calibrate. Skipping file.\n";
        return false;
    }

    // Ordina i picchi per posizione in ampiezza
    std::vector<double> pkPos(sp.GetPositionX(), sp.GetPositionX() + nPeaks);
    std::sort(pkPos.begin(), pkPos.end());

    // Fit lineare: x = indice p.e. (0,1,2,...), y = ampiezza del picco
    int nFit = std::min((int)pkPos.size(), 7);
    TGraphErrors* grCal = new TGraphErrors();
    for (int i = 0; i < nFit; ++i)
        grCal->SetPoint(i, i, pkPos[i]);

    TF1* fLin = new TF1(Form("fLin_%s", tag.c_str()), "pol1", -0.5, nFit - 0.5);
    grCal->Fit(fLin, "RQ");
    q_cal = fLin->GetParameter(0);  // intercetta [mV]
    m_cal = fLin->GetParameter(1);  // guadagno [mV/p.e.]

    std::cout << "  [CAL] Gain = " << m_cal << " mV/p.e.   Offset = " << q_cal << " mV\n";

    // Canvas di calibrazione
    TCanvas* cCal = new TCanvas(Form("cCal_%s", tag.c_str()),
        Form("Calibration — %s", tag.c_str()), 1200, 500);
    cCal->Divide(2, 1);

    cCal->cd(1);
    gPad->SetLogy(); gPad->SetGrid();
    hSpec->SetLineColor(kAzure+1); hSpec->SetLineWidth(2);
    hSpec->Draw("HIST");
    for (int i = 0; i < nFit; ++i) {
        TLine* lp = new TLine(pkPos[i], 0, pkPos[i], hSpec->GetMaximum());
        lp->SetLineColor(kRed+1); lp->SetLineStyle(2); lp->SetLineWidth(2);
        lp->Draw("same");
    }

    cCal->cd(2);
    gPad->SetGrid();
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
//  CALCOLO DEL TOT con interpolazione lineare
//  Restituisce {t_rise, t_fall} sulla soglia fissa let_thr.
//  Ritorna {-1, -1} se non trova entrambi i crossing.
// ════════════════════════════════════════════════════════════
static std::pair<double,double> computeTOT(const std::vector<double>& time,
                                           const std::vector<double>& amp,
                                           double threshold,
                                           int j_start, int j_end) {
    double t_rise = -1.0, t_fall = -1.0;

    // Salita: primo crossing in avanti
    for (int j = j_start + 1; j <= j_end; ++j) {
        if (amp[j-1] < threshold && amp[j] >= threshold) {
            t_rise = time[j-1] + (threshold - amp[j-1])
                     * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
            break;
        }
    }
    if (t_rise < 0) return {-1.0, -1.0};

    // Discesa: primo crossing dopo t_rise
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
//  FIT Q-GAUSS ASIMMETRICO sulle slice del TH2D
// ════════════════════════════════════════════════════════════
static Double_t qGaussAsym(Double_t* x, Double_t* par) {
    double xx = x[0], A = par[0], mean = par[1], sigma = par[2];
    double q1 = par[3], q2 = par[4];
    if (sigma <= 0) return 0;
    if (xx <= mean) {
        double arg = 1.0 - (1.0-q1)*(1.0/(3.0-q1))*std::pow((xx-mean)/sigma, 2);
        return (arg <= 0) ? 0 : A * std::pow(arg, 1.0/(1.0-q1));
    } else {
        double arg = 1.0 - (1.0-q2)*(1.0/(3.0-q2))*std::pow((xx-mean)/sigma, 2);
        return (arg <= 0) ? 0 : A * std::pow(arg, 1.0/(1.0-q2));
    }
}

static std::tuple<double,double,double,double,TF1*>
fitSlice(TH1D* h, int idx, const std::string& tag) {
    if (!h || h->GetEntries() < 40) return {0,0,0,0,nullptr};
    int maxBin    = h->GetMaximumBin();
    double center = h->GetBinCenter(maxBin);
    double sigEst = std::max(h->GetRMS(), 0.3);
    std::string fname = Form("fSlice_%s_%d", tag.c_str(), idx);
    TF1* f = new TF1(fname.c_str(), qGaussAsym,
                     center - 3.0*sigEst, center + 3.0*sigEst, 5);
    f->SetParameters(h->GetMaximum(), center, sigEst, 1.0, 1.2);
    f->FixParameter(3, 1.0);
    h->Fit(f, "RMQN");
    return {f->GetParameter(1), f->GetParError(1),
            std::abs(f->GetParameter(2)), f->GetParError(2), f};
}

// ════════════════════════════════════════════════════════════
//  PROCESS A SINGLE FILE
// ════════════════════════════════════════════════════════════
static void processFile(const FileInfo& info,
                        double cutoff_MHz,
                        double t_trig_start,
                        double t_trig_end,
                        double frac_pe) {

    TFile* file = TFile::Open(info.path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[SKIP] Cannot open: " << info.path << "\n"; return;
    }
    TTree* treeCh1   = (TTree*)file->Get("ch1");
    TTree* treeLaser = (TTree*)file->Get("laser");
    if (!treeCh1) {
        std::cerr << "[SKIP] Tree 'ch1' not found in: " << info.path << "\n";
        file->Close(); return;
    }
    if (!treeLaser) {
        std::cerr << "[SKIP] Tree 'laser' not found in: " << info.path
                  << " — laser trigger required.\n";
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

    // Indici della finestra trigger
    int j_trig_start, j_trig_end;
    triggerWindowIndices(t1, N, t_trig_start, t_trig_end,
                         j_trig_start, j_trig_end, info.tag);

    std::cout << "\n[" << info.tag << "]  "
              << treeCh1->GetEntries() << " events,  fs=" << fs_MHz << " MHz\n"
              << "  Trigger window : [" << t1[j_trig_start] << ", "
              << t1[j_trig_end] << "] ns\n";

    // ── Calibrazione spettro p.e. ────────────────────────────
    double m_cal = 0, q_cal = 0;
    if (!calibrateSpectrum(treeCh1, cutoff_MHz, fs_MHz,
                           j_trig_start, j_trig_end,
                           info.tag, m_cal, q_cal)) {
        file->Close(); return;
    }

    // Soglia LET fissa in mV
    double let_thr = q_cal + frac_pe * m_cal;
    std::cout << "  LET threshold  : " << let_thr << " mV"
              << "  (= " << frac_pe << " p.e.)\n";

    // ── Istogramma 2D: TOT (X) vs Delta_t (Y) ───────────────
    TH2D* h2D = new TH2D(
        Form("h2D_%s", info.tag.c_str()),
        Form("TOT vs #Delta t   V_{bias}=%.0f V, filter=%.0f, LET=%.2f p.e."
             ";TOT (ns);#Delta t = t_{LET} - t_{laser} (ns)",
             info.vbias, info.filter, frac_pe),
        200, 0.0, 50.0,
        400, -30.0, 100.0);
    h2D->SetDirectory(nullptr);

    Long64_t nEntries = treeCh1->GetEntries();
    long nNoLaser = 0, nNoTOT = 0, nFilled = 0;

    std::cout << "  [TOT] Processing events...\n";

    // ── Event loop ───────────────────────────────────────────
    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  [TOT] Progress: " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        treeLaser->GetEntry(i);

        // Trigger laser
        double t_laser = laserTriggerTime(tL, aL, N, 10.0);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        // Baseline correction + filtro passa basso
        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, BASELINE_START, BASELINE_END),
            cutoff_MHz, fs_MHz);

        // TOT sulla soglia fissa let_thr
        auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                                            j_trig_start, j_trig_end);
        if (t_rise < 0 || t_fall < 0) { ++nNoTOT; continue; }

        double tot     = t_fall - t_rise;
        double delta_t = t_rise - t_laser;

        if (tot > 0 && tot < 100.0) {
            h2D->Fill(tot, delta_t);
            ++nFilled;
        }
    }
    std::cout << "\r  [TOT] Progress: " << nEntries << " / " << nEntries << " — done.\n"
              << "  Filled: "      << nFilled
              << "  |  No laser: " << nNoLaser
              << "  |  No TOT: "   << nNoTOT << "\n";

    // ── Canvas 1: mappa 2D ───────────────────────────────────
    TCanvas* c2D = new TCanvas(Form("c2D_%s", info.tag.c_str()),
        Form("TOT map — %s", info.tag.c_str()), 900, 700);
    c2D->SetRightMargin(0.15);
    c2D->SetLeftMargin(PAD_LEFT);
    c2D->SetBottomMargin(PAD_BOTTOM);
    c2D->SetTopMargin(PAD_TOP);
    c2D->SetLogz(); c2D->SetGrid();
    h2D->Draw("COLZ");
    c2D->Update(); c2D->Modified();
    c2D->SaveAs(Form("tot_map_%s.png", info.tag.c_str()));

    // ── Canvas 2: slice in TOT ───────────────────────────────
    const int N_SLICES = 6;
    double tot_min = 0.0, tot_max = 0.0;
    for (int bx = h2D->GetNbinsX(); bx >= 1; --bx) {
        double s = 0;
        for (int by = 1; by <= h2D->GetNbinsY(); ++by)
            s += h2D->GetBinContent(bx, by);
        if (s > 10) { tot_max = h2D->GetXaxis()->GetBinCenter(bx); break; }
    }
    if (tot_max <= 0) tot_max = 20.0;
    double slice_w = (tot_max - tot_min) / N_SLICES;

    TCanvas* cSl = new TCanvas(Form("cSl_%s", info.tag.c_str()),
        Form("TOT slices — %s", info.tag.c_str()), 1200, 800);
    cSl->Divide(3, 2);

    TGraphErrors* grRes  = new TGraphErrors();
    TGraphErrors* grMean = new TGraphErrors();
    int pts = 0;

    for (int s = 0; s < N_SLICES; ++s) {
        double xlo_s = tot_min + s * slice_w;
        double xhi_s = xlo_s + slice_w;
        int bx1 = h2D->GetXaxis()->FindBin(xlo_s);
        int bx2 = h2D->GetXaxis()->FindBin(xhi_s) - 1;

        std::string hname = Form("hsl_%s_%d", info.tag.c_str(), s);
        TH1D* hsl = h2D->ProjectionY(hname.c_str(), bx1, bx2);
        hsl->SetDirectory(nullptr);
        hsl->SetTitle(Form("TOT #in [%.1f, %.1f] ns;#Delta t (ns);Counts",
                           xlo_s, xhi_s));
        hsl->GetXaxis()->UnZoom();
        hsl->GetYaxis()->UnZoom();

        cSl->cd(s+1);
        auto [mean, meanErr, sigma, sigmaErr, fFit] =
            fitSlice(hsl, s, info.tag);

        if (fFit) {
            hsl->Draw("HIST");
            fFit->SetLineColor(kGreen+2);
            fFit->Draw("SAME");
            double tot_center = 0.5*(xlo_s + xhi_s);
            grRes->SetPoint(pts,  tot_center, sigma);
            grRes->SetPointError(pts, 0, sigmaErr);
            grMean->SetPoint(pts, tot_center, mean);
            grMean->SetPointError(pts, 0, meanErr);
            ++pts;
        }
    }
    cSl->Update(); cSl->Modified();
    cSl->SaveAs(Form("tot_slices_%s.png", info.tag.c_str()));

    // ── Canvas 3: trend sigma e mean vs TOT ─────────────────
    if (pts > 1) {
        TCanvas* cTrend = new TCanvas(Form("cTrend_%s", info.tag.c_str()),
            Form("TOT trend — %s", info.tag.c_str()), 1200, 500);
        cTrend->Divide(2, 1);

        cTrend->cd(1);
        gPad->SetGrid();
        grRes->SetMarkerStyle(20);
        grRes->SetMarkerColor(kAzure+1);
        grRes->SetLineColor(kAzure+1);
        grRes->SetTitle(Form(
            "#sigma_{#Delta t} vs TOT   LET=%.2f p.e."
            ";TOT (ns);#sigma_{#Delta t} (ns)", frac_pe));
        grRes->Draw("AP");

        if (pts >= 3) {
            TF1* fStoc = new TF1(Form("fStoc_%s", info.tag.c_str()),
                "sqrt([0]*[0]/x + [1]*[1])",
                grRes->GetX()[0] * 0.8,
                grRes->GetX()[pts-1] * 1.2);
            fStoc->SetParameters(1.0, 0.1);
            fStoc->SetLineColor(kRed+1);
            grRes->Fit(fStoc, "RQ");
            std::cout << "  Stochastic fit: S=" << fStoc->GetParameter(0)
                      << " ns*sqrt(ns)   C=" << fStoc->GetParameter(1) << " ns\n";
        }

        cTrend->cd(2);
        gPad->SetGrid();
        grMean->SetMarkerStyle(20);
        grMean->SetMarkerColor(kOrange+7);
        grMean->SetLineColor(kOrange+7);
        grMean->SetTitle(Form(
            "Mean #Delta t vs TOT   LET=%.2f p.e."
            ";TOT (ns);Mean #Delta t (ns)", frac_pe));
        grMean->Draw("AP");

        cTrend->Update(); cTrend->Modified();
        cTrend->SaveAs(Form("tot_trend_%s.png", info.tag.c_str()));
    }

    // ── Salva ROOT file ──────────────────────────────────────
    std::string outname = Form("tot_%s_let%.2fpe.root",
                               info.tag.c_str(), frac_pe);
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
void sipm_tot_analysis() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::vector<FileInfo> files = findInputFiles();
    if (files.empty()) {
        std::cerr << "No data.vbias_{V}_L.root files found.\n"; return;
    }

    double cutoff_MHz;
    std::cout << "\nLow-pass filter cut-off [MHz]: ";
    std::cin >> cutoff_MHz;

    // Asse temporale dal primo file
    double t_wf_min = 0.0, t_wf_max = 0.0;
    {
        TFile* f0 = TFile::Open(files[0].path.c_str(), "READ");
        if (f0 && !f0->IsZombie()) {
            TTree* tr = (TTree*)f0->Get("ch1");
            if (tr) {
                const int N0 = 1024;
                Double_t t0[N0];
                tr->SetBranchAddress("time", t0);
                tr->GetEntry(0);
                t_wf_min = t0[0];
                t_wf_max = t0[N0-1];
            }
            f0->Close();
        }
    }

    double t_trig_start, t_trig_end;
    std::cout << "\n--- Trigger search window ---\n"
              << "  Waveform range  : [" << t_wf_min << ", " << t_wf_max << "] ns\n"
              << "  Baseline (fixed): [" << BASELINE_START << ", "
                                         << BASELINE_END   << ") ns\n"
              << "  Start [ns]: ";
    std::cin >> t_trig_start;
    std::cout << "  End   [ns]: ";
    std::cin >> t_trig_end;
    if (t_trig_start >= t_trig_end) {
        std::cerr << "[ERROR] t_trig_start >= t_trig_end. Aborting.\n"; return;
    }

    // Frazione p.e. per la soglia LET — es. 0.5 = mezza p.e., 1.0 = 1 p.e.
    double frac_pe;
    std::cout << "\nLET threshold in p.e. units (e.g. 0.5 = half p.e., 1.0 = 1 p.e.): ";
    std::cin >> frac_pe;
    while (frac_pe <= 0.0) {
        std::cout << "  [ERROR] Please enter a positive value: ";
        std::cin >> frac_pe;
    }

    std::cout << "\nSettings summary:\n"
              << "  LP cut-off     : " << cutoff_MHz    << " MHz\n"
              << "  Trigger window : [" << t_trig_start << ", "
                                        << t_trig_end   << "] ns\n"
              << "  LET fraction   : " << frac_pe       << " p.e.\n\n";

    for (auto& f : files)
        processFile(f, cutoff_MHz, t_trig_start, t_trig_end, frac_pe);

    std::cout << "\n=== All files processed. ===\n";
}
