// plot_rising_edge_slope.cpp — versione stabile con parsing robusto
//
// Analisi del fronte di salita SiPM tramite interpolazione lineare.
//
// Uso interattivo:
//   .L plot_rising_edge_slope.cpp+
//   plot_rising_edge_slope()                // scegli da lista
//   plot_rising_edge_slope(55, 0.6, 0.0)    // salta selezione
//   plot_rising_edge_slope(55, 0.6, 0.0, 44.0, 47.0)

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/EventCache.h"
#include "../header/ChunkedHistoFill.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TAxis.h>

#ifndef G_ANALYSIS_MODE_DEFINED
int g_analysis_mode = 0;
#define G_ANALYSIS_MODE_DEFINED
#endif

// ─────────────────────────────────────────────────────────────────────────────
// Elenca tutte le cache v3 e fa scegliere all'utente
static std::string selectCacheInteractively(const std::string& dataDir = DATA_DIR) {
    void* dirp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp) {
        std::cerr << "[ERROR] Cannot open data directory: " << dataDir << "\n";
        return "";
    }

    struct CacheInfo {
        std::string path;
        int         vbias;
        double      let_pe;
        double      cutoff_MHz;
    };
    std::vector<CacheInfo> caches;

    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        std::string fname(entry);
        if (fname.find("events_v3_vbias") != 0) continue;
        if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
        if (fname.find(".tmp") != std::string::npos) continue;

        // Estrai vbias
        size_t p_v = fname.find("vbias");
        size_t p_let = fname.find("_let", p_v);
        if (p_v == std::string::npos || p_let == std::string::npos) continue;
        std::string v_str = fname.substr(p_v+5, p_let - (p_v+5));
        int vbias = 0;
        try { vbias = std::stoi(v_str); } catch(...) { continue; }

        // Estrai LET (p.e.) — robusto: da "_let" fino a "pe"
        size_t let_start = fname.find("_let", p_v);
        if (let_start == std::string::npos) continue;
        let_start += 4;
        size_t let_end = fname.find("pe", let_start);
        if (let_end == std::string::npos) continue;
        std::string let_str = fname.substr(let_start, let_end - let_start);
        double let_pe = 0.0;
        try { let_pe = std::stod(let_str); } catch(...) { continue; }

        // Estrai cutoff (MHz) — robusto: da "_cut" fino a "mhz"
        size_t cut_start = fname.find("_cut");
        if (cut_start == std::string::npos) continue;
        cut_start += 4;
        size_t cut_end = fname.find("mhz", cut_start);
        if (cut_end == std::string::npos) continue;
        std::string cut_str = fname.substr(cut_start, cut_end - cut_start);
        double cutoff = 0.0;
        try { cutoff = std::stod(cut_str); } catch(...) { continue; }

        std::string fullPath = dataDir + "/" + fname;
        if (countCacheEvents(fullPath) == 0) continue;

        caches.push_back({fullPath, vbias, let_pe, cutoff});
    }
    gSystem->FreeDirectory(dirp);

    if (caches.empty()) {
        std::cerr << "[ERROR] Nessuna cache v3 trovata in " << dataDir << "\n";
        return "";
    }

    std::sort(caches.begin(), caches.end(),
        [](const CacheInfo& a, const CacheInfo& b) {
            if (a.vbias != b.vbias) return a.vbias < b.vbias;
            if (std::abs(a.let_pe - b.let_pe) > 1e-4) return a.let_pe < b.let_pe;
            return a.cutoff_MHz < b.cutoff_MHz;
        });

    std::cout << "\n+----------------------------------------------------------+\n"
              << "|  Cache v3 disponibili in:\n|  " << dataDir << "\n"
              << "+----------------------------------------------------------+\n"
              << "|  #  |  Vbias [V]  |  LET [p.e.]  |  Cutoff [MHz]\n"
              << "+-----+--------------+---------------+-----------------+\n";
    for (size_t i = 0; i < caches.size(); ++i) {
        std::cout << "| " << std::setw(3) << i+1 << " | "
                  << std::setw(12) << caches[i].vbias << " | "
                  << std::setw(13) << std::fixed << std::setprecision(2) << caches[i].let_pe << " | "
                  << std::setw(15) << std::setprecision(0) << caches[i].cutoff_MHz << " |\n";
    }
    std::cout << "+-----+--------------+---------------+-----------------+\n";

    int selection = -1;
    while (true) {
        std::cout << "Scegli il numero della cache (1-" << caches.size() << ") o 0 per uscire: ";
        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cout << "\nInput interrotto.\n";
            return "";
        }
        try {
            selection = std::stoi(line);
        } catch(...) {
            std::cout << "  [!] Inserisci un numero valido.\n";
            continue;
        }
        if (selection == 0) return "";
        if (selection >= 1 && selection <= (int)caches.size()) break;
        std::cout << "  [!] Numero fuori range.\n";
    }

    const CacheInfo& chosen = caches[selection-1];
    std::cout << "\n  --> Cache selezionata: Vbias=" << chosen.vbias
              << " V  LET=" << chosen.let_pe << " p.e.  Cutoff=" << chosen.cutoff_MHz << " MHz\n"
              << "  File: " << chosen.path << "\n\n";

    return chosen.path;
}

// ── Fit gaussiano su slice TOT ────────────────────────────────────────────
static std::pair<double,double> fitBinMu(TH2D* h2, int ib,
                                          double fit_lo, double fit_hi)
{
    TH1D* hY = h2->ProjectionY(Form("__hY_%s_%d", h2->GetName(), ib), ib, ib);
    hY->SetDirectory(nullptr);
    if (hY->Integral() < 15) { delete hY; return {std::nan(""), std::nan("")}; }
    int b1 = hY->FindBin(fit_lo+1e-4), b2 = hY->FindBin(fit_hi-1e-4);
    if (b1 >= b2) { delete hY; return {std::nan(""), std::nan("")}; }
    int bPk = b1;
    for (int b = b1+1; b <= b2; ++b)
        if (hY->GetBinContent(b) > hY->GetBinContent(bPk)) bPk = b;
    double mu0 = hY->GetBinCenter(bPk);
    double sg0 = std::max(0.08, hY->GetRMS()*0.4);
    double flo = std::max(fit_lo, mu0-3*sg0), fhi = std::min(fit_hi, mu0+3*sg0);
    if (fhi-flo < 0.05) { delete hY; return {std::nan(""), std::nan("")}; }
    TF1* g = new TF1(Form("__g_%s_%d", h2->GetName(), ib), "gaus", flo, fhi);
    g->SetParameters(hY->GetBinContent(bPk), mu0, sg0);
    int st = hY->Fit(g, "RQN");
    double mu = std::nan(""), err = std::nan("");
    if (st == 0 || st == 4000) {
        double mf = g->GetParameter(1), ef = g->GetParError(1);
        if (std::isfinite(mf) && mf > fit_lo && mf < fit_hi &&
            std::isfinite(ef) && ef > 0.0) { mu = mf; err = ef; }
    }
    delete g; delete hY;
    return {mu, err};
}

// ── Funzione principale ───────────────────────────────────────────────────
void plot_rising_edge_slope(
        int    vbias      = -1,     // -1 = interattivo
        double thr_opt    = -1.0,
        double cutoff_MHz = -1.0,
        double fit_lo     = 44.0,
        double fit_hi     = 47.0,
        int    n_tot_bins = 20,
        double tot_lo     = 5.0,
        double tot_hi     = 50.0)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    const std::string dataDir = DATA_DIR;

    std::string cachePath;
    double pe_lo, pe_hi, delta_pe;
    std::string tag;

    // ── Selezione cache ──────────────────────────────────────────────────
    if (vbias < 0 || thr_opt < 0 || cutoff_MHz < 0) {
        // Modalità interattiva
        cachePath = selectCacheInteractively(dataDir);
        if (cachePath.empty()) return;

        // Estrai parametri dal nome del file selezionato
        std::string fname = cachePath.substr(cachePath.find_last_of("/\\") + 1);
        size_t p_v = fname.find("vbias");
        size_t p_let = fname.find("_let", p_v);
        size_t let_start = p_let + 4;
        size_t let_end = fname.find("pe", let_start);
        size_t cut_start = fname.find("_cut") + 4;
        size_t cut_end = fname.find("mhz", cut_start);

        vbias      = std::stoi(fname.substr(p_v+5, p_let - (p_v+5)));
        thr_opt    = std::stod(fname.substr(let_start, let_end - let_start));
        cutoff_MHz = std::stod(fname.substr(cut_start, cut_end - cut_start));

        pe_lo    = 0.2 * thr_opt;
        pe_hi    = 0.8 * thr_opt;
        delta_pe = pe_hi - pe_lo;
        tag      = Form("vbias%d_thr%.2fpe_cut%.0fMHz", vbias, thr_opt, cutoff_MHz);
    } else {
        // Modalità non interattiva
        pe_lo    = 0.2 * thr_opt;
        pe_hi    = 0.8 * thr_opt;
        delta_pe = pe_hi - pe_lo;
        tag      = Form("vbias%d_thr%.2fpe_cut%.0fMHz", vbias, thr_opt, cutoff_MHz);

        // Cerca la cache più vicina (retrocompatibilità)
        auto findBestCache = [&](int vb, double thr, double cut) -> std::string {
            void* dirp = gSystem->OpenDirectory(dataDir.c_str());
            if (!dirp) return "";
            std::string pat_v = Form("events_v3_vbias%d_let", vb);
            std::string pat_c = Form("_cut%dmhz", (int)std::round(cut));
            double best_d = 1e18;
            std::string best_path;
            const char* e = nullptr;
            while ((e = gSystem->GetDirEntry(dirp)) != nullptr) {
                std::string fn(e);
                if (fn.size() < 5 || fn.substr(fn.size()-5) != ".root") continue;
                if (fn.find(pat_v) == std::string::npos) continue;
                if (fn.find(pat_c) == std::string::npos) continue;
                if (fn.find(".tmp") != std::string::npos) continue;
                size_t p0 = fn.find("_let"), p1 = fn.find("pe_", p0);
                if (p0 == std::string::npos || p1 == std::string::npos) continue;
                double frac = 0.0;
                try { frac = std::stod(fn.substr(p0+4, p1-p0-4)); } catch (...) { continue; }
                if (frac <= 0.0 || frac > 10.0) continue;
                std::string fp = dataDir + "/" + fn;
                if (countCacheEvents(fp) == 0) continue;
                double d = std::abs(frac - thr);
                if (d < best_d) { best_d = d; best_path = fp; }
            }
            gSystem->FreeDirectory(dirp);
            if (best_path.empty())
                std::cerr << "[FindCache] No v3 cache for vbias=" << vb
                          << " thr=" << thr << " cutoff=" << cut << "\n"
                          << "  Rigenera con sipm_tot_analysis aggiornato.\n";
            return best_path;
        };
        cachePath = findBestCache(vbias, thr_opt, cutoff_MHz);
        if (cachePath.empty()) return;
    }

    std::cout << "\n+============================================+\n"
              << "|  Rising Edge Slope Analysis (v3 cache)     |\n"
              << "+============================================+\n"
              << "  Vbias     = " << vbias      << " V\n"
              << "  thr_opt   = " << thr_opt    << " p.e.\n"
              << "  pe_lo     = " << pe_lo      << " p.e. (0.2x)\n"
              << "  pe_hi     = " << pe_hi      << " p.e. (0.8x)\n"
              << "  cutoff    = " << cutoff_MHz << " MHz\n"
              << "  fit [ns]  = [" << fit_lo << ", " << fit_hi << "]\n"
              << "  TOT [ns]  = [" << tot_lo << ", " << tot_hi << "]\n"
              << "  Cache     : " << cachePath << "\n\n";

    // Verifica branch v3
    {
        TFile* fc = TFile::Open(cachePath.c_str(), "READ");
        if (!fc || fc->IsZombie()) { std::cerr << "[ERROR] Cannot open cache.\n"; delete fc; return; }
        TTree* tc = (TTree*)fc->Get("events");
        bool ok = tc && tc->GetBranch("delta_t_lo") && tc->GetBranch("delta_t_hi");
        fc->Close(); delete fc;
        if (!ok) {
            std::cerr << "[ERROR] Cache senza delta_t_lo/hi — rigenera con sipm_tot_analysis v3.\n";
            return;
        }
    }

    OutCtx ctx = createOutputDirs(Form("rising_edge_%s", tag.c_str()));

    // Riempi tre TH2D in un singolo passaggio
    double dt_lo = fit_lo - 8.0, dt_hi = fit_hi + 8.0;
    int    n_dt  = std::max(300, (int)std::round((dt_hi-dt_lo)/0.1));

    TH2D* h2_opt = new TH2D(Form("h2opt_%s",tag.c_str()),
        Form("TOT vs #Deltat  %.2f p.e.  Vbias %d V;TOT (ns);#Deltat (ns)", thr_opt, vbias),
        n_tot_bins, tot_lo, tot_hi, n_dt, dt_lo, dt_hi);
    h2_opt->SetDirectory(nullptr);

    TH2D* h2_lo = new TH2D(Form("h2lo_%s",tag.c_str()),
        Form("TOT vs #Deltat_{lo}  %.2f p.e.;TOT (ns);#Deltat (ns)", pe_lo),
        n_tot_bins, tot_lo, tot_hi, n_dt, dt_lo, dt_hi);
    h2_lo->SetDirectory(nullptr);

    TH2D* h2_hi = new TH2D(Form("h2hi_%s",tag.c_str()),
        Form("TOT vs #Deltat_{hi}  %.2f p.e.;TOT (ns);#Deltat (ns)", pe_hi),
        n_tot_bins, tot_lo, tot_hi, n_dt, dt_lo, dt_hi);
    h2_hi->SetDirectory(nullptr);

    std::cout << "  Filling (single pass)...\n";
    {
        ChunkedFiller filler(cachePath);
        filler.addTH2D(h2_opt, [tot_lo,tot_hi,dt_lo,dt_hi](TH2D* h, const CacheEvent& e){
            if (e.tot>tot_lo && e.tot<tot_hi && e.delta_t>dt_lo && e.delta_t<dt_hi)
                h->Fill(e.tot, e.delta_t);
        });
        filler.addTH2D(h2_lo, [tot_lo,tot_hi,dt_lo,dt_hi](TH2D* h, const CacheEvent& e){
            if (e.tot>tot_lo && e.tot<tot_hi &&
                e.delta_t_lo > TRISE_INVALID+1.0 &&
                e.delta_t_lo > dt_lo && e.delta_t_lo < dt_hi)
                h->Fill(e.tot, e.delta_t_lo);
        });
        filler.addTH2D(h2_hi, [tot_lo,tot_hi,dt_lo,dt_hi](TH2D* h, const CacheEvent& e){
            if (e.tot>tot_lo && e.tot<tot_hi &&
                e.delta_t_hi > TRISE_INVALID+1.0 &&
                e.delta_t_hi > dt_lo && e.delta_t_hi < dt_hi)
                h->Fill(e.tot, e.delta_t_hi);
        });
        Long64_t n = filler.run();
        std::cout << "  [Fill] " << n << " events  h2_opt="
                  << (long)h2_opt->GetEntries() << "  h2_lo="
                  << (long)h2_lo->GetEntries()  << "  h2_hi="
                  << (long)h2_hi->GetEntries()  << "\n\n";
    }

    // Calcola m e q per bin TOT
    struct BinRes { double tc,mlo,mhi,elo,ehi,m,q,em,eq; };
    std::vector<BinRes> res;
    res.reserve(n_tot_bins);

    std::cout << "  Fitting slices...\n";
    for (int ib = 1; ib <= n_tot_bins; ++ib) {
        double tc = h2_lo->GetXaxis()->GetBinCenter(ib);
        auto [mlo,elo] = fitBinMu(h2_lo, ib, fit_lo, fit_hi);
        auto [mhi,ehi] = fitBinMu(h2_hi, ib, fit_lo, fit_hi);
        if (!std::isfinite(mlo) || !std::isfinite(mhi)) {
            std::cout << "    bin " << std::setw(2) << ib
                      << "  TOT=" << std::fixed << std::setprecision(1) << tc << " -> skip\n";
            continue;
        }
        BinRes br;
        br.tc  = tc;
        br.mlo = mlo; br.elo = std::isfinite(elo) ? elo : 0.0;
        br.mhi = mhi; br.ehi = std::isfinite(ehi) ? ehi : 0.0;
        br.m   = (mhi-mlo)/delta_pe;
        br.em  = std::hypot(br.elo, br.ehi)/std::abs(delta_pe);
        br.q   = mlo - br.m*pe_lo;
        br.eq  = std::hypot(br.elo, std::abs(pe_lo)*br.em);
        res.push_back(br);
        std::cout << "    bin " << std::setw(2) << ib
                  << "  TOT=" << std::fixed << std::setprecision(1) << tc
                  << "  m=" << std::setprecision(3) << br.m
                  << " ns/pe  q=" << br.q << " ns\n";
    }

    if (res.empty()) {
        std::cerr << "[ERROR] No valid bins.\n";
        delete h2_opt; delete h2_lo; delete h2_hi; return;
    }

    int np = (int)res.size();
    double hb = (tot_hi-tot_lo)/n_tot_bins*0.5;
    std::vector<double> vt(np),vm(np),vq(np),vet(np,hb),vem(np),veq(np);
    for (int i=0;i<np;++i){ vt[i]=res[i].tc; vm[i]=res[i].m; vq[i]=res[i].q;
                              vem[i]=res[i].em; veq[i]=res[i].eq; }

    TGraphErrors* gr_m = new TGraphErrors(np, vt.data(),vm.data(),vet.data(),vem.data());
    gr_m->SetName(Form("gr_m_%s",tag.c_str()));
    gr_m->SetTitle(Form("Slope m vs TOT  Vbias %d V  thr_{opt}=%.2f p.e.;TOT (ns);m (ns/p.e.)",vbias,thr_opt));
    gr_m->SetMarkerStyle(20); gr_m->SetMarkerSize(1.0);
    gr_m->SetMarkerColor(kAzure+1); gr_m->SetLineColor(kAzure+1); gr_m->SetLineWidth(2);

    TGraphErrors* gr_q = new TGraphErrors(np, vt.data(),vq.data(),vet.data(),veq.data());
    gr_q->SetName(Form("gr_q_%s",tag.c_str()));
    gr_q->SetTitle(Form("Arrival time q vs TOT  Vbias %d V  thr_{opt}=%.2f p.e.;TOT (ns);q (ns)",vbias,thr_opt));
    gr_q->SetMarkerStyle(21); gr_q->SetMarkerSize(1.0);
    gr_q->SetMarkerColor(kRed+1); gr_q->SetLineColor(kRed+1); gr_q->SetLineWidth(2);

    // Canvas 1
    TCanvas* c1 = new TCanvas(Form("c1_%s",tag.c_str()),
                               Form("TOT vs #Deltat  Vbias %d V",vbias), 1300, 550);
    c1->Divide(2,1);

    c1->cd(1);
    gPad->SetLogz(); gPad->SetGrid();
    gPad->SetLeftMargin(PAD_LEFT); gPad->SetBottomMargin(PAD_BOTTOM);
    gPad->SetRightMargin(0.15);   gPad->SetTopMargin(PAD_TOP);
    h2_opt->Draw("COLZ");

    TProfile* plo = h2_lo->ProfileX(Form("plo_%s",tag.c_str()),1,-1,"s");
    TProfile* phi = h2_hi->ProfileX(Form("phi_%s",tag.c_str()),1,-1,"s");
    plo->SetLineColor(kGreen+2);  plo->SetLineWidth(2);
    phi->SetLineColor(kOrange+7); phi->SetLineWidth(2);
    plo->Draw("same"); phi->Draw("same");

    TLegend* leg1 = new TLegend(0.18,0.74,0.58,0.88);
    leg1->SetBorderSize(1); leg1->SetFillColor(0);
    leg1->SetTextFont(42);  leg1->SetTextSize(0.034);
    leg1->AddEntry(plo, Form("#mu @ %.2f p.e. (0.2x)",pe_lo), "l");
    leg1->AddEntry(phi, Form("#mu @ %.2f p.e. (0.8x)",pe_hi), "l");
    leg1->Draw();

    c1->cd(2);
    gPad->SetGrid();
    gPad->SetLeftMargin(PAD_LEFT+0.02); gPad->SetBottomMargin(PAD_BOTTOM); gPad->SetTopMargin(PAD_TOP);
    gr_m->Draw("AP"); gr_m->GetYaxis()->SetTitleOffset(1.5);
    { TPaveText* pt=new TPaveText(0.55,0.74,0.94,0.88,"NDC");
      pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextFont(42); pt->SetTextSize(0.032);
      pt->AddText(Form("Vbias = %d V",vbias));
      pt->AddText(Form("thr_{opt} = %.2f p.e.",thr_opt));
      pt->AddText(Form("lo=%.2f  hi=%.2f p.e.",pe_lo,pe_hi));
      pt->Draw(); }
    c1->Update(); c1->Modified();
    ctx.savePNG(c1, Form("rising_edge_c1_%s.png",tag.c_str()));

    // Canvas 2
    TCanvas* c2 = new TCanvas(Form("c2_%s",tag.c_str()),
                               Form("Arrival time q vs TOT  Vbias %d V",vbias), 900, 650);
    c2->SetGrid();
    c2->SetLeftMargin(PAD_LEFT+0.02); c2->SetBottomMargin(PAD_BOTTOM); c2->SetTopMargin(PAD_TOP);
    gr_q->Draw("AP"); gr_q->GetYaxis()->SetTitleOffset(1.5);

    TF1* fLin = new TF1(Form("fLin_%s",tag.c_str()),"pol1",vt.front()-hb,vt.back()+hb);
    fLin->SetLineColor(kGray+2); fLin->SetLineStyle(2); fLin->SetLineWidth(2);
    gr_q->Fit(fLin,"RQ");

    { TPaveText* pt2=new TPaveText(0.15,0.68,0.68,0.88,"NDC");
      pt2->SetBorderSize(1); pt2->SetFillColor(0); pt2->SetTextFont(42); pt2->SetTextSize(0.032);
      pt2->AddText(Form("Vbias = %d V",vbias));
      pt2->AddText(Form("thr_{opt} = %.2f p.e.",thr_opt));
      pt2->AddText(Form("lo=%.2f  hi=%.2f p.e.",pe_lo,pe_hi));
      pt2->AddText(Form("N bins valid = %d / %d",np,n_tot_bins));
      pt2->AddText(Form("Lin fit: q = %.3f + %.4f #cdot TOT",
                        fLin->GetParameter(0),fLin->GetParameter(1)));
      pt2->Draw(); }
    c2->Update(); c2->Modified();
    ctx.savePNG(c2, Form("rising_edge_c2_%s.png",tag.c_str()));

    std::cout << "\n[OK] " << ctx.pngDir << "\n"
              << "[OK] N valid bins: " << np << " / " << n_tot_bins << "\n";

    delete h2_opt; delete h2_lo; delete h2_hi;
}