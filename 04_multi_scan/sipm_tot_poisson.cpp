// sipm_tot_poisson.cpp
// Analisi TOT con calibrazione tramite fit poissoniano (Approccio 1).
//
// Lo spettro delle ampiezze è fittato con una somma di gaussiane
// pesate da coefficienti poissoniani:
//
//   S(x) = sum_{n=0}^{N} [ e^{-mu} * mu^n / n! ] * G(x; q+n*m, sigma_n)
//   sigma_n = sqrt(sigma_0^2 + n * sigma_1^2)
//
// Questo approccio è fisicamente il più corretto: funziona anche con
// pochissima luce (mu piccolo) perché sfrutta la struttura statistica
// della distribuzione di Poisson per vincolare m e q anche quando
// i picchi p.e. non sono visivamente distinguibili.
//
// Selezione interattiva dei file per l'analisi TOT.
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
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>

// ════════════════════════════════════════════════════════════
//  SELEZIONE INTERATTIVA DEI FILE
// ════════════════════════════════════════════════════════════
static void printFileList(const std::vector<FileInfo>& files) {
    std::cout << "\nFiles found:\n";
    for (int i = 0; i < (int)files.size(); ++i)
        std::cout << "  [" << i << "] " << files[i].tag << "\n";
}

static std::vector<int> selectMultipleFiles(const std::vector<FileInfo>& files,
                                             const std::string& prompt) {
    std::vector<int> selected;
    std::string line;
    std::cin.ignore();
    while (selected.empty()) {
        std::cout << prompt;
        std::getline(std::cin, line);
        if (line == "all" || line == "ALL") {
            for (int i = 0; i < (int)files.size(); ++i) selected.push_back(i);
            break;
        }
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
//  FIT POISSONIANO DELLO SPETTRO p.e.
//
//  Modello: somma di N_PE+1 gaussiane (da n=0 a n=N_PE),
//  ciascuna pesata dal coefficiente di Poisson P(n; mu).
//  Parametri liberi: mu, m [mV/p.e.], q [mV], sigma_0, sigma_1
//
//  sigma_n = sqrt(sigma_0^2 + n * sigma_1^2)
//    sigma_0 = rumore elettronico (larghezza del picco a 0 p.e.)
//    sigma_1 = fluttuazione del guadagno per singolo p.e.
// ════════════════════════════════════════════════════════════
// N_PE_MAX ridotto a 5: con mu piccolo i termini n>5 sono trascurabili
// e ridurre il numero di termini accelera il fit e migliora la stabilità
static const int N_PE_MAX = 5;

static Double_t poissonSpectrum(Double_t* x, Double_t* par) {
    // par[0] = mu      (media Poisson)
    // par[1] = q       (offset [mV])
    // par[2] = m       (guadagno [mV/p.e.])
    // par[3] = sigma_0 (rumore elettronico [mV])
    // par[4] = sigma_1 (fluttuazione guadagno [mV/sqrt(p.e.)])
    // par[5] = N       (normalizzazione)
    double mu      = par[0];
    double q       = par[1];
    double m       = par[2];
    double sigma_0 = par[3];
    double sigma_1 = par[4];
    double norm    = par[5];

    if (mu <= 0 || m <= 0 || sigma_0 <= 0) return 0;

    double sum = 0;
    double logMu = std::log(mu);
    for (int n = 0; n <= N_PE_MAX; ++n) {
        double logP = -mu + n * logMu - TMath::LnGamma(n + 1);
        double P    = std::exp(logP);
        if (P < 1e-10) continue;
        double mean_n  = q + n * m;
        double sigma_n = std::sqrt(sigma_0*sigma_0 + n*sigma_1*sigma_1);
        if (sigma_n <= 0) continue;
        sum += P * TMath::Gaus(x[0], mean_n, sigma_n, true);
    }
    return norm * sum;
}

static bool calibratePoisson(TTree* treeCh1,
                              double cutoff_MHz,
                              double fs_MHz,
                              int j_trig_start, int j_trig_end,
                              const std::string& tag,
                              double& m_cal, double& q_cal) {
    const int N = 1024;
    Double_t t1[N], a1[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);

    TH1D* hSpec = new TH1D(Form("hSpecP_%s", tag.c_str()),
        Form("p.e. spectrum (Poisson fit)   %s;Amplitude (mV);Counts",
             tag.c_str()),
        300, -5.0, 60.0);
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
        // Massimo sull'intera waveform — stesso approccio delle macro originali
        double ampMax = *std::max_element(af.begin(), af.end());
        if (ampMax > -1e8) hSpec->Fill(ampMax);
    }
    std::cout << "\r  [CAL] " << nEntries << " / " << nEntries << " — done.\n";

    // ── Stima iniziale robusta ───────────────────────────────
    // Trova il bin del picco principale (0 p.e.) — è sempre il massimo
    int bin0     = hSpec->GetMaximumBin();
    double q_est = hSpec->GetBinCenter(bin0);  // posizione del picco a 0 p.e.

    // Stima sigma_0 dalla larghezza a metà altezza del picco principale
    double hmax  = hSpec->GetMaximum();
    double sigma0_est = 0.5;  // fallback
    for (int b = bin0; b < hSpec->GetNbinsX(); ++b) {
        if (hSpec->GetBinContent(b) < hmax * 0.5) {
            sigma0_est = hSpec->GetBinCenter(b) - q_est;
            break;
        }
    }
    sigma0_est = std::max(sigma0_est, 0.1);

    // Stima mu dalla frazione di eventi al di sopra di q + 2*sigma0
    // P(n=0) = e^{-mu} → mu = -ln(N0/Ntot)
    double thr_noise = q_est + 2.0 * sigma0_est;
    double n_above   = hSpec->Integral(hSpec->FindBin(thr_noise),
                                        hSpec->GetNbinsX());
    double n_tot     = hSpec->GetEntries();
    double frac_above = std::max(n_above / n_tot, 1e-4);
    double mu_est    = std::max(-std::log(1.0 - frac_above), 0.01);
    mu_est           = std::min(mu_est, 5.0);

    // Stima m dal centroide degli eventi sopra soglia
    double sum_w = 0, sum_wx = 0;
    for (int b = hSpec->FindBin(thr_noise); b <= hSpec->GetNbinsX(); ++b) {
        double c = hSpec->GetBinContent(b);
        double x = hSpec->GetBinCenter(b);
        sum_w += c; sum_wx += c * x;
    }
    double mean_above = (sum_w > 0) ? sum_wx / sum_w : q_est + 5.0;
    double m_est      = std::max(mean_above - q_est, 1.0);

    // sigma_1 fissato a 20% di sigma_0 — mal determinato con poca luce,
    // fissarlo evita che il fit diventi instabile
    double sigma1_est = sigma0_est * 0.2;
    double norm_est   = hSpec->GetMaximum() * hSpec->GetBinWidth(1);

    std::cout << "  [CAL] Initial estimates:"
              << "  mu="      << mu_est
              << "  m="       << m_est      << " mV/p.e."
              << "  q="       << q_est      << " mV"
              << "  sigma_0=" << sigma0_est << " mV\n";

    // Fit range: dal minimo dello spettro a q + (N_PE_MAX+1)*m
    double fitLo = hSpec->GetXaxis()->GetXmin();
    double fitHi = std::min(q_est + (N_PE_MAX + 1.5) * m_est,
                            hSpec->GetXaxis()->GetXmax());

    TF1* fPois = new TF1(Form("fPois_%s", tag.c_str()),
                          poissonSpectrum, fitLo, fitHi, 6);
    fPois->SetParNames("mu", "q", "m", "sigma_0", "sigma_1", "N");
    fPois->SetParameters(mu_est, q_est, m_est, sigma0_est, sigma1_est, norm_est);

    fPois->SetParLimits(0, 1e-4, 10.0);                      // mu
    fPois->SetParLimits(1, q_est - 2.0, q_est + 2.0);        // q vicino al picco 0 p.e.
    fPois->SetParLimits(2, 0.5, 50.0);                        // m
    fPois->SetParLimits(3, 0.05, sigma0_est * 3.0);           // sigma_0
    fPois->FixParameter(4, sigma1_est);                        // sigma_1 fissato
    fPois->SetParLimits(5, 1.0, 1e9);                         // N

    fPois->SetLineColor(kRed+1);
    fPois->SetLineWidth(2);
    fPois->SetNpx(1000);

    // "N" = no draw durante fit, "L" = log-likelihood, più robusto con bin vuoti
    int fitStatus = hSpec->Fit(fPois, "RQLN");
    if (fitStatus != 0)
        std::cout << "  [CAL] WARNING: fit status = " << fitStatus
                  << " — result may be unreliable.\n";

    m_cal = fPois->GetParameter(2);
    q_cal = fPois->GetParameter(1);
    double mu_fit      = fPois->GetParameter(0);
    double sigma0_fit  = fPois->GetParameter(3);
    double sigma1_fit  = fPois->GetParameter(4);
    double chi2_ndf    = (fPois->GetNDF() > 0)
                         ? fPois->GetChisquare() / fPois->GetNDF() : -1;

    std::cout << "  [CAL] Poisson fit result:\n"
              << "    mu      = " << mu_fit     << "\n"
              << "    q       = " << q_cal      << " mV\n"
              << "    m       = " << m_cal      << " mV/p.e.\n"
              << "    sigma_0 = " << sigma0_fit << " mV\n"
              << "    sigma_1 = " << sigma1_fit << " mV/sqrt(p.e.)\n"
              << "    chi2/ndf= " << chi2_ndf   << "\n";

    if (chi2_ndf > 10.0)
        std::cout << "  [CAL] WARNING: chi2/ndf = " << chi2_ndf
                  << " > 10 — fit quality is poor.\n"
                  << "  Consider using sipm_tot_crosscal for a more robust result.\n";

    // Canvas di calibrazione
    TCanvas* cCal = new TCanvas(Form("cCalP_%s", tag.c_str()),
        Form("Poisson calibration — %s", tag.c_str()), 900, 700);
    cCal->SetLeftMargin(PAD_LEFT); cCal->SetBottomMargin(PAD_BOTTOM);
    cCal->SetTopMargin(PAD_TOP);   cCal->SetRightMargin(PAD_RIGHT);
    cCal->SetGrid();

    hSpec->SetLineColor(kAzure+1); hSpec->SetLineWidth(2);
    hSpec->Draw("HIST");
    fPois->Draw("SAME");

    // Linee verticali per ogni picco p.e.
    for (int n = 0; n <= std::min(5, N_PE_MAX); ++n) {
        double xp = q_cal + n * m_cal;
        if (xp < fitLo || xp > fitHi) continue;
        TLine* lp = new TLine(xp, 0, xp, hSpec->GetMaximum()*0.8);
        lp->SetLineColor(kGreen+2); lp->SetLineStyle(2); lp->SetLineWidth(2);
        lp->Draw("same");
    }

    TPaveText* pt = new TPaveText(0.55, 0.65, 0.93, 0.88, "NDC");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
    pt->SetTextFont(42);  pt->SetTextSize(0.033);
    pt->AddText(Form("#mu = %.3f",           mu_fit));
    pt->AddText(Form("Gain  m = %.3f mV/p.e.", m_cal));
    pt->AddText(Form("Offset q = %.3f mV",    q_cal));
    pt->AddText(Form("#sigma_{0} = %.3f mV",  sigma0_fit));
    pt->AddText(Form("#chi^{2}/ndf = %.2f",   chi2_ndf));
    pt->Draw();

    cCal->Update(); cCal->Modified();
    cCal->SaveAs(Form("tot_poisson_cal_%s.png", tag.c_str()));

    return (m_cal > 0);
}

// ════════════════════════════════════════════════════════════
//  CALCOLO TOT con interpolazione lineare
// ════════════════════════════════════════════════════════════
static std::pair<double,double> computeTOT(const std::vector<double>& time,
                                           const std::vector<double>& amp,
                                           double threshold,
                                           int j_start, int j_end) {
    double t_rise = -1.0, t_fall = -1.0;
    for (int j = j_start+1; j <= j_end; ++j) {
        if (amp[j-1] < threshold && amp[j] >= threshold) {
            t_rise = time[j-1] + (threshold-amp[j-1])
                     *(time[j]-time[j-1])/(amp[j]-amp[j-1]);
            break;
        }
    }
    if (t_rise < 0) return {-1.0,-1.0};
    for (int j = j_start+1; j <= j_end; ++j) {
        if (time[j] <= t_rise) continue;
        if (amp[j-1] >= threshold && amp[j] < threshold) {
            t_fall = time[j-1] + (threshold-amp[j-1])
                     *(time[j]-time[j-1])/(amp[j]-amp[j-1]);
            break;
        }
    }
    if (t_fall < 0) return {-1.0,-1.0};
    return {t_rise, t_fall};
}

// ════════════════════════════════════════════════════════════
//  FIT Q-GAUSS ASIMMETRICO sulle slice
// ════════════════════════════════════════════════════════════
static Double_t qGaussAsym(Double_t* x, Double_t* par) {
    double xx=x[0], A=par[0], mean=par[1], sigma=par[2], q1=par[3], q2=par[4];
    if (sigma<=0) return 0;
    if (xx<=mean) {
        double arg=1.0-(1.0-q1)*(1.0/(3.0-q1))*std::pow((xx-mean)/sigma,2);
        return (arg<=0)?0:A*std::pow(arg,1.0/(1.0-q1));
    } else {
        double arg=1.0-(1.0-q2)*(1.0/(3.0-q2))*std::pow((xx-mean)/sigma,2);
        return (arg<=0)?0:A*std::pow(arg,1.0/(1.0-q2));
    }
}

static std::tuple<double,double,double,double,TF1*>
fitSlice(TH1D* h, int idx, const std::string& tag) {
    if (!h || h->GetEntries()<40) return {0,0,0,0,nullptr};
    int maxBin    = h->GetMaximumBin();
    double center = h->GetBinCenter(maxBin);
    double sigEst = std::max(h->GetRMS(), 0.3);
    std::string fname = Form("fSliceP_%s_%d", tag.c_str(), idx);
    TF1* f = new TF1(fname.c_str(), qGaussAsym,
                     center-3.0*sigEst, center+3.0*sigEst, 5);
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
              << "  (= " << frac_pe << " p.e., Poisson-calibrated)\n";

    TH2D* h2D = new TH2D(
        Form("h2DP_%s", info.tag.c_str()),
        Form("TOT vs #Delta t   V_{bias}=%.0f V, filter=%.0f, LET=%.2f p.e."
             ";TOT (ns);#Delta t = t_{LET} - t_{laser} (ns)",
             info.vbias, info.filter, frac_pe),
        200, 0.0, 50.0, 400, -30.0, 100.0);
    h2D->SetDirectory(nullptr);

    Long64_t nEntries = treeCh1->GetEntries();
    long nNoLaser=0, nNoTOT=0, nFilled=0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  [TOT] " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i); treeLaser->GetEntry(i);

        double t_laser = laserTriggerTime(tL, aL, N, 10.0);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        std::vector<double> v_t(t1,t1+N), v_a(a1,a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t,v_a,BASELINE_START,BASELINE_END),
            cutoff_MHz, fs_MHz);

        auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                                            j_trig_start, j_trig_end);
        if (t_rise<0 || t_fall<0) { ++nNoTOT; continue; }

        double tot     = t_fall - t_rise;
        double delta_t = t_rise - t_laser;
        if (tot>0 && tot<100.0) { h2D->Fill(tot, delta_t); ++nFilled; }
    }
    std::cout << "\r  [TOT] " << nEntries << " / " << nEntries << " — done.\n"
              << "  Filled: " << nFilled
              << "  |  No laser: " << nNoLaser
              << "  |  No TOT: "   << nNoTOT << "\n";

    // Canvas mappa 2D
    TCanvas* c2D = new TCanvas(Form("c2DP_%s", info.tag.c_str()),
        Form("TOT map — %s", info.tag.c_str()), 900, 700);
    c2D->SetRightMargin(0.15); c2D->SetLeftMargin(PAD_LEFT);
    c2D->SetBottomMargin(PAD_BOTTOM); c2D->SetTopMargin(PAD_TOP);
    c2D->SetLogz(); c2D->SetGrid();
    h2D->Draw("COLZ");
    c2D->Update(); c2D->Modified();
    c2D->SaveAs(Form("tot_map_poisson_%s.png", info.tag.c_str()));

    // Slice
    const int N_SLICES = 6;
    double tot_max = 0.0;
    for (int bx=h2D->GetNbinsX(); bx>=1; --bx) {
        double s=0;
        for (int by=1; by<=h2D->GetNbinsY(); ++by) s+=h2D->GetBinContent(bx,by);
        if (s>10) { tot_max=h2D->GetXaxis()->GetBinCenter(bx); break; }
    }
    if (tot_max<=0) tot_max=20.0;
    double slice_w = tot_max / N_SLICES;

    TCanvas* cSl = new TCanvas(Form("cSlP_%s", info.tag.c_str()),
        Form("TOT slices — %s", info.tag.c_str()), 1200, 800);
    cSl->Divide(3,2);

    TGraphErrors* grRes  = new TGraphErrors();
    TGraphErrors* grMean = new TGraphErrors();
    int pts=0;

    for (int s=0; s<N_SLICES; ++s) {
        double xlo_s = s*slice_w, xhi_s = xlo_s+slice_w;
        int bx1=h2D->GetXaxis()->FindBin(xlo_s);
        int bx2=h2D->GetXaxis()->FindBin(xhi_s)-1;
        std::string hname=Form("hslP_%s_%d", info.tag.c_str(), s);
        TH1D* hsl=h2D->ProjectionY(hname.c_str(), bx1, bx2);
        hsl->SetDirectory(nullptr);
        hsl->SetTitle(Form("TOT #in [%.1f, %.1f] ns;#Delta t (ns);Counts",
                           xlo_s, xhi_s));
        hsl->GetXaxis()->UnZoom(); hsl->GetYaxis()->UnZoom();
        cSl->cd(s+1);
        auto [mean,meanErr,sigma,sigmaErr,fFit]=fitSlice(hsl,s,info.tag);
        if (fFit) {
            hsl->Draw("HIST");
            fFit->SetLineColor(kGreen+2); fFit->Draw("SAME");
            double tc=0.5*(xlo_s+xhi_s);
            grRes->SetPoint(pts,tc,sigma);   grRes->SetPointError(pts,0,sigmaErr);
            grMean->SetPoint(pts,tc,mean);   grMean->SetPointError(pts,0,meanErr);
            ++pts;
        }
    }
    cSl->Update(); cSl->Modified();
    cSl->SaveAs(Form("tot_slices_poisson_%s.png", info.tag.c_str()));

    if (pts>1) {
        TCanvas* cTrend=new TCanvas(Form("cTrendP_%s",info.tag.c_str()),
            Form("TOT trend — %s",info.tag.c_str()), 1200, 500);
        cTrend->Divide(2,1);
        cTrend->cd(1); gPad->SetGrid();
        grRes->SetMarkerStyle(20); grRes->SetMarkerColor(kAzure+1);
        grRes->SetLineColor(kAzure+1);
        grRes->SetTitle(Form("#sigma_{#Delta t} vs TOT   LET=%.2f p.e."
            ";TOT (ns);#sigma_{#Delta t} (ns)", frac_pe));
        grRes->Draw("AP");
        if (pts>=3) {
            TF1* fS=new TF1(Form("fSP_%s",info.tag.c_str()),
                "sqrt([0]*[0]/x+[1]*[1])",
                grRes->GetX()[0]*0.8, grRes->GetX()[pts-1]*1.2);
            fS->SetParameters(1.0,0.1); fS->SetLineColor(kRed+1);
            grRes->Fit(fS,"RQ");
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
        cTrend->SaveAs(Form("tot_trend_poisson_%s.png", info.tag.c_str()));
    }

    std::string outname=Form("tot_poisson_%s_let%.2fpe.root",
                              info.tag.c_str(), frac_pe);
    TFile* fout=new TFile(outname.c_str(),"RECREATE");
    h2D->Write();
    if (pts>0) { grRes->Write("grSigma"); grMean->Write("grMean"); }
    fout->Write(); fout->Close();
    std::cout << "  Saved: " << outname << "\n";

    file->Close();
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════
void sipm_tot_poisson() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::vector<FileInfo> files = findInputFiles();
    if (files.empty()) {
        std::cerr << "No data.vbias_{V}_L.root files found.\n"; return;
    }

    printFileList(files);
    std::vector<int> anaIdx = selectMultipleFiles(files,
        "Select files for analysis (e.g. \"2\" or \"0,2\" or \"all\"): ");

    double cutoff_MHz;
    std::cout << "\nLow-pass filter cut-off [MHz]: ";
    std::cin >> cutoff_MHz;

    // Asse temporale
    double t_wf_min=0.0, t_wf_max=0.0;
    {
        TFile* f0=TFile::Open(files[anaIdx[0]].path.c_str(),"READ");
        if (f0 && !f0->IsZombie()) {
            TTree* tr=(TTree*)f0->Get("ch1");
            if (tr) {
                const int N0=1024; Double_t t0[N0];
                tr->SetBranchAddress("time",t0); tr->GetEntry(0);
                t_wf_min=t0[0]; t_wf_max=t0[N0-1];
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
    while (frac_pe<=0.0) {
        std::cout << "  [ERROR] Enter a positive value: ";
        std::cin >> frac_pe;
    }

    std::cout << "\nSettings summary:\n"
              << "  LP cut-off     : " << cutoff_MHz    << " MHz\n"
              << "  Trigger window : [" << t_trig_start << ", "
                                        << t_trig_end   << "] ns\n"
              << "  LET fraction   : " << frac_pe       << " p.e.\n\n";

    // Per ogni file selezionato: calibrazione poissoniana + analisi TOT
    for (int idx : anaIdx) {
        TFile* fAna = TFile::Open(files[idx].path.c_str(), "READ");
        if (!fAna || fAna->IsZombie()) {
            std::cerr << "[SKIP] Cannot open: " << files[idx].path << "\n";
            continue;
        }
        TTree* tAna = (TTree*)fAna->Get("ch1");
        if (!tAna) { fAna->Close(); continue; }

        const int NA=1024; Double_t ta[NA];
        tAna->SetBranchAddress("time",ta); tAna->GetEntry(0);
        double fs_a = 1000.0/(ta[1]-ta[0]);
        int ja_start, ja_end;
        triggerWindowIndices(ta, NA, t_trig_start, t_trig_end,
                             ja_start, ja_end, files[idx].tag);

        double m_cal=0, q_cal=0;
        bool ok = calibratePoisson(tAna, cutoff_MHz, fs_a,
                                   ja_start, ja_end,
                                   files[idx].tag, m_cal, q_cal);
        fAna->Close();

        if (!ok) {
            std::cerr << "[SKIP] Poisson calibration failed for "
                      << files[idx].tag << "\n";
            continue;
        }

        double let_thr = q_cal + frac_pe * m_cal;
        std::cout << "\n  LET threshold = " << let_thr << " mV"
                  << "  (= " << frac_pe << " p.e.)\n";

        processFile(files[idx], cutoff_MHz, t_trig_start, t_trig_end,
                    let_thr, frac_pe);
    }

    std::cout << "\n=== All files processed. ===\n";
}
