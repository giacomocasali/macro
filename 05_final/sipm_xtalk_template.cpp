// ════════════════════════════════════════════════════════════════
//  sipm_xtalk_template.cpp
//  ----------------------------------------------------------------
//  Discriminazione xtalk vs 2-pe vere via TEMPLATE FIT su waveform.
//
//  PERCHÉ NON USARE m  o  dV/dt₂/dV/dt₁ ?
//  ----------------------------------------
//  1. m = (Δt₈₀% − Δt₂₀%) / 0.6
//     Stima del rise-time grossolana: 2 punti, dominata da jitter
//     elettronico (~150 ps su Δt → σ_m ≈ 0.07 ns/p.e.). La firma
//     vera dello xtalk ha Δm ≈ 50–80 ps/p.e. → SOTTO il rumore.
//     Risultato sui dati reali: separazione 3.2σ ma sovrapposizione
//     >80%, cut Fisher inutile.
//
//  2. dV/dt₂/dV/dt₁ (rapporto picchi della derivata)
//     Box-smooth k=3 + sampling 200 ps → distrugge le features
//     sub-ns dove xtalk si manifesta. Le distribuzioni low-light e
//     high-light risultano sovrapposte (vedi 07_c3_dvdt_overlay).
//
//  3. Rise-time 10–90% (img 10) ha INFORMAZIONE vera: low-light ha
//     coda lunga, high-light e' stretto. Ma viene buttata via.
//
//  COSA USIAMO INVECE
//  -------------------
//  Template T(τ): forma media normalizzata di un singolo p.e.
//  costruito su low-light selezionando 1pe puliti.
//  Per ogni evento candidato (amp ∈ [1.5, 2.5] p.e.) facciamo:
//
//    fit_1:  m(t) = α·T(t−t0) + b
//    fit_2:  m(t) = α₁·T(t−t0) + α₂·T(t−t0−Δt) + b
//
//  Discriminanti:
//    Δt              → separazione fra i due fotoni
//                      xtalk: ≈ 0  (entro 50–100 ps)
//                      2-pe vere: jitter laser × √2 (>200 ps)
//    α₂/α₁           → rapporto cariche (≈1 per pari, <1 per disp.)
//    Δχ²=χ²₁−χ²₂     → guadagno informativo del modello a 2 fotoni
//
//  La frazione di xtalk si stima fittando dt come somma di:
//    (a) gaussiana centrata su <dt>≈jitter laser
//    (b) gaussiana stretta centrata in 0 (xtalk)
//
//  ----------------------------------------------------------------
//  Dipendenze (header riutilizzati dal progetto):
//     Config.h  CalibIO.h  Calibration.h  EventCache.h
//     SignalProcessing.h  ButterworthFilter.h  ProgressBar.h
//     OutputManager.h
//
//  Header NUOVI (in questo programma):
//     WaveformTemplate.h    XtalkDiscriminator.h    TemplateFit.h
//
//  Compile:
//     .L sipm_xtalk_template.cpp+
//  Run:
//     sipm_xtalk_template(55, 200.0)   // vbias, cutoff
// ════════════════════════════════════════════════════════════════

#include "../header/Config.h"
#include "../header/CalibIO.h"
#include "../header/Calibration.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/ProgressBar.h"

#include "../header/WaveformTemplate.h"
#include "../header/XtalkDiscriminator.h"
#include "../header/TemplateFit.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>

// global mode (riutilizziamo i cache esistenti)
#ifndef G_ANALYSIS_MODE_DEFINED
inline int g_analysis_mode = 0;
#define G_ANALYSIS_MODE_DEFINED
#endif

static constexpr int N_SAMPLES = 1024;

// ─── Utility: trova run files (riusa la logica di xtalk) ─────────
static std::map<int, std::string>
findRunFilesXT(int vbias, const std::string& dataDir) {
    std::map<int,std::string> runs;
    const std::string patBrace = makeRunPattern(vbias);
    const std::string patPlain = RUN_PREFIX + std::to_string(vbias) + "_run_";
    void* dirp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp) return runs;
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        std::string fn(entry);
        bool match = (fn.find(patBrace) != std::string::npos);
        if (!match) {
            size_t p = fn.find(patPlain);
            if (p != std::string::npos) {
                bool start_ok = (p == 0);
                bool end_ok = (p + patPlain.size() < fn.size())
                              && std::isdigit((unsigned char)fn[p + patPlain.size()]);
                if (start_ok && end_ok) match = true;
            }
        }
        if (!match) continue;
        if (fn.size() < 5 || fn.substr(fn.size()-5) != ".root") continue;
        size_t pos = fn.find("_run_");
        if (pos == std::string::npos) continue;
        try {
            std::string sub = fn.substr(pos+5);
            size_t dot = sub.find(".root");
            if (dot != std::string::npos) sub = sub.substr(0, dot);
            runs[std::stoi(sub)] = dataDir + "/" + fn;
        } catch (...) {}
    }
    gSystem->FreeDirectory(dirp);
    return runs;
}

// ─── Stima sampling rate dalla prima entry di un file ────────────
static double estimateSamplingMHz(const std::string& path) {
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return 5000.0; }
    TTree* tr = (TTree*)f->Get("ch1");
    if (!tr || tr->GetEntries() < 1) { f->Close(); delete f; return 5000.0; }
    Double_t tb[N_SAMPLES];
    tr->SetBranchAddress("time", tb);
    tr->GetEntry(0);
    double dt = tb[1] - tb[0];
    f->Close(); delete f;
    if (dt <= 0) return 5000.0;
    return 1000.0 / dt;
}

// ─── Plot template + RMS band ────────────────────────────────────
static void plotTemplate(const xtalk::Template1PE& T, const std::string& outDir,
                          int vbias) {
    if (!T.ok) return;
    int N = (int)T.tau.size();
    std::vector<double> tauv = T.tau, mu = T.mean, sig = T.rms;
    std::vector<double> upper(N), lower(N);
    for (int i = 0; i < N; ++i) {
        upper[i] = mu[i] + sig[i];
        lower[i] = mu[i] - sig[i];
    }

    TCanvas c("c_template", "1pe template", 900, 600);
    c.SetGrid();
    c.SetMargin(PAD_LEFT, PAD_RIGHT, PAD_BOTTOM, PAD_TOP);

    TGraph* gMu = new TGraph(N, tauv.data(), mu.data());
    gMu->SetLineColor(kBlue+1); gMu->SetLineWidth(2);
    gMu->SetTitle(Form("1pe template Vbias=%d V (N=%d);#tau (ns);Amplitude (norm)",
                       vbias, T.N_used));

    TGraph* gUp = new TGraph(N, tauv.data(), upper.data());
    TGraph* gLo = new TGraph(N, tauv.data(), lower.data());
    gUp->SetLineColor(kAzure+7); gUp->SetLineStyle(2);
    gLo->SetLineColor(kAzure+7); gLo->SetLineStyle(2);

    gMu->Draw("AL");
    gMu->GetYaxis()->SetRangeUser(-0.2, 1.3);
    gUp->Draw("L SAME");
    gLo->Draw("L SAME");

    auto* leg = new TLegend(0.55, 0.75, 0.92, 0.88);
    leg->SetTextSize(0.035); leg->SetBorderSize(1); leg->SetFillColor(0);
    leg->AddEntry(gMu, "Mean template", "l");
    leg->AddEntry(gUp, "Mean +/- RMS",  "l");
    leg->Draw();

    c.SaveAs(Form("%s/01_template_vbias%d.png", outDir.c_str(), vbias));
}

// ─── Plot dt distribution low vs high ────────────────────────────
static void plotDtDistribution(const xtalk::LightStats& sLo,
                                const xtalk::LightStats& sHi,
                                const std::string& outDir, int vbias,
                                double min_ratio = 0.10,
                                double min_dchi2 = 5.0) {
    TH1D* hLo = new TH1D(Form("h_dt_lo_v%d", vbias),
        ";#Deltat (ns);Freq. norm.", 100, -0.5, 2.0);
    TH1D* hHi = new TH1D(Form("h_dt_hi_v%d", vbias),
        ";#Deltat (ns);Freq. norm.", 100, -0.5, 2.0);
    hLo->SetDirectory(nullptr); hHi->SetDirectory(nullptr);
    hLo->SetLineColor(kAzure+1);  hLo->SetLineWidth(2);
    hHi->SetLineColor(kRed+1);    hHi->SetLineWidth(2);
    hLo->SetFillColorAlpha(kAzure+1, 0.30); hLo->SetFillStyle(3244);
    hHi->SetFillColorAlpha(kRed+1, 0.30);   hHi->SetFillStyle(3245);

    auto pass = [&](const xtalk::CandidateEvent& e){
        return e.fit_ok && e.ratio >= min_ratio && e.dchi2 >= min_dchi2;
    };
    for (auto& e : sLo.events) if (pass(e)) hLo->Fill(e.dt);
    for (auto& e : sHi.events) if (pass(e)) hHi->Fill(e.dt);

    auto norm = [](TH1D* h){ if (h->Integral()>0) h->Scale(1.0/h->Integral()); };
    norm(hLo); norm(hHi);

    TCanvas c("c_dt", "dt distribution", 900, 600);
    c.SetGrid();
    c.SetMargin(PAD_LEFT, PAD_RIGHT, PAD_BOTTOM, PAD_TOP);
    double ymax = std::max(hLo->GetMaximum(), hHi->GetMaximum()) * 1.25;
    hHi->SetMaximum(ymax); hHi->SetMinimum(0);
    hHi->Draw("HIST");
    hLo->Draw("HIST SAME");

    auto* leg = new TLegend(0.55, 0.72, 0.92, 0.88);
    leg->SetBorderSize(1); leg->SetFillColor(0); leg->SetTextSize(0.032);
    leg->AddEntry(hLo, Form("Low light  N=%d  <dt>=%.3f #pm %.3f ns",
                  sLo.N_fit_ok, sLo.mean_dt, sLo.rms_dt), "lf");
    leg->AddEntry(hHi, Form("High light N=%d  <dt>=%.3f #pm %.3f ns",
                  sHi.N_fit_ok, sHi.mean_dt, sHi.rms_dt), "lf");
    leg->Draw();

    auto* pt = new TPaveText(PAD_LEFT+0.02, 0.78, 0.50, 0.92, "NDC");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextSize(0.032);
    pt->AddText(Form("Vbias = %d V", vbias));
    pt->AddText(Form("Amp [1.5, 2.5] p.e."));
    pt->AddText(Form("CUTS: #alpha_{2}/#alpha_{1}>%.2f, #Delta#chi^{2}>%.0f",
                     min_ratio, min_dchi2));
    pt->Draw();

    c.SaveAs(Form("%s/02_dt_overlay_vbias%d.png", outDir.c_str(), vbias));
}

// ─── Plot ratio_alpha & dchi2 ────────────────────────────────────
static void plotAuxDiscriminants(const xtalk::LightStats& sLo,
                                  const xtalk::LightStats& sHi,
                                  const std::string& outDir, int vbias) {
    TH1D* hRLo = new TH1D(Form("h_rL_v%d", vbias),
        ";#alpha_{2}/#alpha_{1};Freq. norm.", 100, 0.0, 1.05);
    TH1D* hRHi = new TH1D(Form("h_rH_v%d", vbias),
        ";#alpha_{2}/#alpha_{1};Freq. norm.", 100, 0.0, 1.05);
    TH1D* hCLo = new TH1D(Form("h_cL_v%d", vbias),
        ";#chi^{2}_{1} - #chi^{2}_{2};Freq. norm.", 100, -10, 200);
    TH1D* hCHi = new TH1D(Form("h_cH_v%d", vbias),
        ";#chi^{2}_{1} - #chi^{2}_{2};Freq. norm.", 100, -10, 200);

    for (auto* h : {hRLo, hRHi, hCLo, hCHi}) h->SetDirectory(nullptr);
    hRLo->SetLineColor(kAzure+1); hRLo->SetLineWidth(2);
    hRHi->SetLineColor(kRed+1);   hRHi->SetLineWidth(2);
    hCLo->SetLineColor(kAzure+1); hCLo->SetLineWidth(2);
    hCHi->SetLineColor(kRed+1);   hCHi->SetLineWidth(2);

    for (auto& e : sLo.events) if (e.fit_ok) {
        hRLo->Fill(e.ratio); hCLo->Fill(e.dchi2);
    }
    for (auto& e : sHi.events) if (e.fit_ok) {
        hRHi->Fill(e.ratio); hCHi->Fill(e.dchi2);
    }
    auto norm = [](TH1D* h){ if (h->Integral()>0) h->Scale(1.0/h->Integral()); };
    norm(hRLo); norm(hRHi); norm(hCLo); norm(hCHi);

    {
        TCanvas c("c_r", "alpha ratio", 900, 600);
        c.SetGrid();
        c.SetMargin(PAD_LEFT, PAD_RIGHT, PAD_BOTTOM, PAD_TOP);
        double ymax = std::max(hRLo->GetMaximum(), hRHi->GetMaximum()) * 1.25;
        hRHi->SetMaximum(ymax); hRHi->Draw("HIST"); hRLo->Draw("HIST SAME");
        auto* leg = new TLegend(0.15, 0.72, 0.55, 0.88);
        leg->SetTextSize(0.032);
        leg->AddEntry(hRLo, Form("Low  N=%d", sLo.N_fit_ok), "l");
        leg->AddEntry(hRHi, Form("High N=%d", sHi.N_fit_ok), "l");
        leg->Draw();
        c.SaveAs(Form("%s/03_alpha_ratio_vbias%d.png", outDir.c_str(), vbias));
    }
    {
        TCanvas c("c_dchi2", "delta chi2", 900, 600);
        c.SetGrid(); c.SetLogy();
        c.SetMargin(PAD_LEFT, PAD_RIGHT, PAD_BOTTOM, PAD_TOP);
        hCHi->Draw("HIST"); hCLo->Draw("HIST SAME");
        auto* leg = new TLegend(0.55, 0.72, 0.92, 0.88);
        leg->SetTextSize(0.032);
        leg->AddEntry(hCLo, Form("Low  N=%d", sLo.N_fit_ok), "l");
        leg->AddEntry(hCHi, Form("High N=%d", sHi.N_fit_ok), "l");
        leg->Draw();
        c.SaveAs(Form("%s/04_dchi2_vbias%d.png", outDir.c_str(), vbias));
    }
}

// ─── 2D scatter dt vs ratio_alpha ────────────────────────────────
static void plotDtVsRatio(const xtalk::LightStats& s,
                           const std::string& outDir, int vbias,
                           const std::string& tag,
                           double min_dchi2 = 5.0) {
    // Plot raw (tutti)
    TH2D h(Form("h2_dt_r_%s_v%d", tag.c_str(), vbias),
           Form("Delta t vs alpha2/alpha1 -- %s Vbias=%d (raw);#Deltat (ns);#alpha_{2}/#alpha_{1}",
                tag.c_str(), vbias),
           80, -0.3, 2.0, 60, 0.0, 1.05);
    h.SetDirectory(nullptr);
    for (auto& e : s.events)
        if (e.fit_ok) h.Fill(e.dt, e.ratio);

    TCanvas c("c_2d", "scatter", 900, 700);
    c.SetGrid(); c.SetLogz();
    c.SetMargin(PAD_LEFT, 0.13f, PAD_BOTTOM, PAD_TOP);
    h.Draw("COLZ");
    c.SaveAs(Form("%s/05_dt_vs_ratio_%s_vbias%d.png",
                   outDir.c_str(), tag.c_str(), vbias));

    // Plot filtrato (solo good fit doppio)
    TH2D h2(Form("h2_dt_r_clean_%s_v%d", tag.c_str(), vbias),
           Form("Delta t vs alpha2/alpha1 CLEAN -- %s Vbias=%d (#Delta#chi^{2}>%.0f);#Deltat (ns);#alpha_{2}/#alpha_{1}",
                tag.c_str(), vbias, min_dchi2),
           80, -0.3, 2.0, 60, 0.0, 1.05);
    h2.SetDirectory(nullptr);
    for (auto& e : s.events)
        if (e.fit_ok && e.dchi2 >= min_dchi2 && e.ratio > 0.05)
            h2.Fill(e.dt, e.ratio);
    TCanvas c2("c_2d_clean", "scatter clean", 900, 700);
    c2.SetGrid(); c2.SetLogz();
    c2.SetMargin(PAD_LEFT, 0.13f, PAD_BOTTOM, PAD_TOP);
    h2.Draw("COLZ");
    c2.SaveAs(Form("%s/06_dt_vs_ratio_clean_%s_vbias%d.png",
                   outDir.c_str(), tag.c_str(), vbias));
}

// ─── Salva dump TTree con tutti i candidati (per analisi offline) ─
static void saveCandidatesTree(const xtalk::LightStats& s,
                                const std::string& outDir,
                                int vbias, const std::string& tag) {
    std::string path = Form("%s/candidates_%s_vbias%d.root",
                              outDir.c_str(), tag.c_str(), vbias);
    TFile f(path.c_str(), "RECREATE");
    if (f.IsZombie()) return;
    Double_t amp_pe, t_cfd, m_legacy;
    Int_t    fit_ok;
    Double_t dt, dt_err, alpha_a, alpha_b, ratio, sum, dchi2, chi2_1, chi2_2;
    Int_t    ndf_1, ndf_2;
    TTree* tr = new TTree("cand", "xtalk candidates");
    tr->Branch("amp_pe",   &amp_pe);
    tr->Branch("t_cfd",    &t_cfd);
    tr->Branch("m_legacy", &m_legacy);
    tr->Branch("fit_ok",   &fit_ok);
    tr->Branch("dt",       &dt);
    tr->Branch("dt_err",   &dt_err);
    tr->Branch("alpha_a",  &alpha_a);
    tr->Branch("alpha_b",  &alpha_b);
    tr->Branch("ratio",    &ratio);
    tr->Branch("sum",      &sum);
    tr->Branch("dchi2",    &dchi2);
    tr->Branch("chi2_1",   &chi2_1);
    tr->Branch("chi2_2",   &chi2_2);
    tr->Branch("ndf_1",    &ndf_1);
    tr->Branch("ndf_2",    &ndf_2);
    for (auto& e : s.events) {
        amp_pe   = e.amp_pe;
        t_cfd    = e.t_cfd;
        m_legacy = e.m_legacy;
        fit_ok   = e.fit_ok ? 1 : 0;
        dt       = e.dt;       dt_err = e.dt_err;
        alpha_a  = e.alpha_a;  alpha_b = e.alpha_b;
        ratio    = e.ratio;    sum     = e.sum;
        dchi2    = e.dchi2;
        chi2_1   = e.chi2_1;   chi2_2  = e.chi2_2;
        ndf_1    = e.ndf_1;    ndf_2   = e.ndf_2;
        tr->Fill();
    }
    f.Write(); f.Close();
    std::cout << "  [Save] " << path << "\n";
}

// ════════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════════
void sipm_xtalk_template(int vbias = 55, double cutoff_MHz = 0.0,
                          long max_template_events = 5000,
                          long max_candidate_events = 50000) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::cout << "\n=================================================\n"
              << "  SiPM XTALK DISCRIMINATION via TEMPLATE FIT\n"
              << "  Vbias = " << vbias << " V    cutoff = " << cutoff_MHz << " MHz\n"
              << "=================================================\n";

    const std::string dataLow  = DATA_DIR_LOW;
    const std::string dataHigh = DATA_DIR_HIGH;
    std::cout << "  Low  : " << dataLow  << "\n";
    std::cout << "  High : " << dataHigh << "\n";

    // ── 1. carica calibrazione ───────────────────────────────
    CalibResult cal;
    if (!loadCalibration(cal, vbias, cutoff_MHz, dataLow) || !cal.ok) {
        std::cerr << "[ERR] no calibration found\n";
        return;
    }
    std::cout << "  Gain = " << cal.m << " mV/p.e.\n"
              << "  Trig = [" << cal.t_trig_start << ", " << cal.t_trig_end << "] ns\n";

    // ── 2. costruisci template (o ricarica se gia' presente) ──
    auto runsLo = findRunFilesXT(vbias, dataLow);
    auto runsHi = findRunFilesXT(vbias, dataHigh);
    if (runsLo.empty()) {
        std::cerr << "[ERR] no low-light runs\n"; return;
    }

    double fs_MHz = estimateSamplingMHz(runsLo.begin()->second);
    std::cout << "  fs   = " << fs_MHz << " MHz\n\n";

    OutCtx ctx = createOutputDirs("xtalk_template");
    std::string tmpl_path = Form("%s/template_vbias%d.root",
                                  ctx.rootDir.c_str(), vbias);

    xtalk::Template1PE tmpl;
    if (xtalk::loadTemplate(tmpl, tmpl_path)) {
        std::cout << "  [Template] reloaded from " << tmpl_path
                  << " (" << tmpl.N_used << " events)\n";
    } else {
        xtalk::TemplateConfig tcfg;
        tcfg.amp_pe_lo  = 0.85;
        tcfg.amp_pe_hi  = 1.15;
        tcfg.use_filter = false;
        tcfg.cutoff_MHz = cutoff_MHz;
        tcfg.fs_MHz     = fs_MHz;
        tcfg.max_events = max_template_events;
        std::cout << "  [Template] building from low-light...\n";
        tmpl = xtalk::buildTemplate(runsLo, cal, tcfg, N_SAMPLES);
        if (!tmpl.ok) {
            std::cerr << "[ERR] template construction failed\n";
            return;
        }
        xtalk::saveTemplate(tmpl, tmpl_path);
    }
    plotTemplate(tmpl, ctx.pngDir, vbias);

    // ── 3. processa candidati 2-pe in low e high ──────────────
    xtalk::DiscriminatorConfig dcfg;
    dcfg.amp_pe_lo  = 1.5;
    dcfg.amp_pe_hi  = 2.5;
    dcfg.use_filter = false;
    dcfg.cutoff_MHz = cutoff_MHz;
    dcfg.fs_MHz     = fs_MHz;
    dcfg.max_events = max_candidate_events;

    xtalk::LightStats sLo; sLo.label = "low";
    xtalk::LightStats sHi; sHi.label = "high";

    std::cout << "\n--- Processing LOW light ---\n";
    xtalk::processRunsForCandidates(runsLo, cal, tmpl, dcfg, "low",
                                     sLo.events, N_SAMPLES);
    if (!runsHi.empty()) {
        std::cout << "\n--- Processing HIGH light ---\n";
        xtalk::processRunsForCandidates(runsHi, cal, tmpl, dcfg, "high",
                                         sHi.events, N_SAMPLES);
    } else {
        std::cerr << "[WARN] no high-light runs available\n";
    }

    // ── 4. summary stats ──────────────────────────────────────
    xtalk::summarizeStats(sLo);
    xtalk::summarizeStats(sHi);

    std::cout << "\n=== RESULTS (with quality cuts: ratio>0.10, dchi2>5) ===\n";
    std::cout << std::fixed << std::setprecision(4);
    auto report = [](const xtalk::LightStats& s){
        std::cout << "  [" << s.label << "]"
                  << "  N_total=" << s.N_total
                  << "  N_good=" << s.N_fit_ok
                  << "  <dt>=" << s.mean_dt << " ns"
                  << "  rms_dt=" << s.rms_dt << " ns"
                  << "  <dchi2>=" << s.mean_dchi2
                  << "  frac(|dt|<0.15ns)=" << s.frac_short_dt << "\n";
    };
    report(sLo);
    report(sHi);

    // confronto: differenza frac short = stima xtalk
    double xt_excess = sLo.frac_short_dt - sHi.frac_short_dt;
    std::cout << "\n  >>> excess(low - high) at |dt|<0.15ns = "
              << xt_excess * 100.0 << " %\n";
    std::cout << "      (positivo = xtalk in low; negativo = piu' 2pe-veri in high)\n";

    // ── 5. stima frazione xtalk in low light ──────────────────
    double xt_frac = xtalk::estimateXtalkFraction(sLo.events);
    std::cout << "\n  Xtalk fraction in LOW light (likelihood fit): "
              << xt_frac * 100.0 << " %\n";

    // ── 6. plot ───────────────────────────────────────────────
    plotDtDistribution(sLo, sHi, ctx.pngDir, vbias);
    plotAuxDiscriminants(sLo, sHi, ctx.pngDir, vbias);
    plotDtVsRatio(sLo, ctx.pngDir, vbias, "low");
    if (sHi.N_fit_ok > 0)
        plotDtVsRatio(sHi, ctx.pngDir, vbias, "high");

    // ── 7. salva TTree per analisi offline ────────────────────
    saveCandidatesTree(sLo, ctx.rootDir, vbias, "low");
    if (sHi.N_fit_ok > 0)
        saveCandidatesTree(sHi, ctx.rootDir, vbias, "high");

    std::cout << "\nDone. Output in:\n"
              << "  PNG  :  " << ctx.pngDir << "\n"
              << "  ROOT :  " << ctx.rootDir << "\n";
}
