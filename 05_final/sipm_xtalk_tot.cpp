// ════════════════════════════════════════════════════════════════
//  sipm_xtalk_tot.cpp
//  ----------------------------------------------------------------
//  Cerca xtalk via TOT usando i file cache prodotti da
//  sipm_tot_analysis (TTree "events" con tot, amp_max, delta_t, n_pe).
//
//  Input automatico: stessi file cache eventCachePath() del modulo
//  EventCache. Niente riprocessamento waveforms.
//
//  Compile:
//     .L sipm_xtalk_tot.cpp+
//  Run:
//     sipm_xtalk_tot(55, 0.60, 200.0)
//                    vbias  let_pe  cutoff_MHz
// ════════════════════════════════════════════════════════════════

#include "../header/Config.h"
#include "../header/CalibIO.h"
#include "../header/Calibration.h"
#include "../header/OutputManager.h"
#include "../header/EventCache.h"

#include "../header/XtalkTOT.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGraphErrors.h>

#ifndef G_ANALYSIS_MODE_DEFINED
inline int g_analysis_mode = 0;
#define G_ANALYSIS_MODE_DEFINED
#endif

// ─── Plot TOT distribution con fit gauss + cuts marker ───────────
static void plotSlice(const std::vector<xtot::CacheEv>& evLo,
                       const std::vector<xtot::CacheEv>& evHi,
                       const xtot::AmpSlice& slice,
                       const xtot::SliceFit& fLo,
                       const xtot::SliceFit& fHi,
                       const std::string& outDir, int vbias) {
    int N = 300;
    double tot_lo = 0, tot_hi = 150;
    TH1D* hLo = new TH1D(Form("h_lo_n%d", slice.n_pe_target),
        Form("TOT distribution Vbias=%d  amp [%.2f, %.2f] pe;TOT (ns);Freq. norm.",
             vbias, slice.amp_lo, slice.amp_hi),
        N, tot_lo, tot_hi);
    TH1D* hHi = new TH1D(Form("h_hi_n%d", slice.n_pe_target),
        ";TOT (ns);Freq. norm.", N, tot_lo, tot_hi);
    hLo->SetDirectory(nullptr); hHi->SetDirectory(nullptr);
    hLo->SetLineColor(kAzure+1); hLo->SetLineWidth(2);
    hHi->SetLineColor(kRed+1);   hHi->SetLineWidth(2);
    hLo->SetFillColorAlpha(kAzure+1, 0.30); hLo->SetFillStyle(3244);
    hHi->SetFillColorAlpha(kRed+1, 0.30);   hHi->SetFillStyle(3245);

    for (auto& e : evLo)
        if (e.amp_pe >= slice.amp_lo && e.amp_pe < slice.amp_hi
            && e.tot > 0.5 && e.tot < tot_hi) hLo->Fill(e.tot);
    for (auto& e : evHi)
        if (e.amp_pe >= slice.amp_lo && e.amp_pe < slice.amp_hi
            && e.tot > 0.5 && e.tot < tot_hi) hHi->Fill(e.tot);
    auto norm = [](TH1D* h){ if (h->Integral()>0) h->Scale(1.0/h->Integral()); };
    norm(hLo); norm(hHi);

    TCanvas c("c", "tot slice", 1000, 600);
    c.SetGrid();
    c.SetMargin(PAD_LEFT, PAD_RIGHT, PAD_BOTTOM, PAD_TOP);
    double ymax = std::max(hLo->GetMaximum(), hHi->GetMaximum()) * 1.30;
    if (ymax <= 0) ymax = 0.05;
    hHi->SetMaximum(ymax); hHi->SetMinimum(0);
    // zoom intorno ai picchi
    double xmin = std::min(fLo.mu_tot, fHi.mu_tot) - 6 * std::max(fLo.sigma_tot, fHi.sigma_tot);
    double xmax = std::max(fLo.mu_tot, fHi.mu_tot) + 6 * std::max(fLo.sigma_tot, fHi.sigma_tot);
    if (xmin < 0) xmin = 0;
    if (xmax > tot_hi) xmax = tot_hi;
    if (xmax - xmin < 5) { xmin = std::max(0.0, xmin - 5); xmax += 5; }
    hHi->GetXaxis()->SetRangeUser(xmin, xmax);
    hHi->Draw("HIST");
    hLo->Draw("HIST SAME");

    // marker cuts left tail
    auto drawCut = [&](double mu, double sig, int col){
        if (sig <= 0) return;
        double x = mu - 2 * sig;
        TLine* ln = new TLine(x, 0, x, ymax * 0.9);
        ln->SetLineColor(col); ln->SetLineStyle(2); ln->SetLineWidth(2);
        ln->Draw();
    };
    drawCut(fLo.mu_tot, fLo.sigma_tot, kAzure+1);
    drawCut(fHi.mu_tot, fHi.sigma_tot, kRed+1);

    auto* leg = new TLegend(0.55, 0.65, 0.92, 0.88);
    leg->SetBorderSize(1); leg->SetFillColor(0); leg->SetTextSize(0.030);
    leg->AddEntry(hLo, Form("Low N=%d  mu=%.2f  sig=%.2f  fL=%.3f",
                  fLo.N_in_slice, fLo.mu_tot, fLo.sigma_tot, fLo.frac_left), "lf");
    leg->AddEntry(hHi, Form("High N=%d  mu=%.2f  sig=%.2f  fL=%.3f",
                  fHi.N_in_slice, fHi.mu_tot, fHi.sigma_tot, fHi.frac_left), "lf");
    leg->Draw();

    auto* pt = new TPaveText(PAD_LEFT+0.02, 0.70, 0.50, 0.88, "NDC");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetTextSize(0.030);
    pt->AddText(Form("Vbias = %d V", vbias));
    pt->AddText(Form("amp slice: %.2f - %.2f pe  (target %d pe)",
                     slice.amp_lo, slice.amp_hi, slice.n_pe_target));
    double exc = fLo.frac_left - fHi.frac_left;
    pt->AddText(Form("excess(low - high) left tail = %.3f", exc));
    pt->Draw();

    c.SaveAs(Form("%s/tot_slice_n%dpe_vbias%d.png",
                   outDir.c_str(), slice.n_pe_target, vbias));
}

// ─── Plot summary: TOT vs amp 2D + medie picco ───────────────────
static void plotTOTvsAmp(const std::vector<xtot::CacheEv>& ev,
                          const std::vector<xtot::SliceFit>& fits,
                          const std::string& outDir, int vbias,
                          const std::string& tag) {
    TH2D h(Form("h2_%s_v%d", tag.c_str(), vbias),
           Form("TOT vs amp -- %s Vbias=%d;amp (pe);TOT (ns)",
                tag.c_str(), vbias),
           80, 0, 4.5, 200, 0, 100);
    h.SetDirectory(nullptr);
    for (auto& e : ev)
        if (e.amp_pe > 0 && e.amp_pe < 4.5 && e.tot > 0.5 && e.tot < 100)
            h.Fill(e.amp_pe, e.tot);

    TCanvas c("c2d","tot vs amp",900,700);
    c.SetGrid(); c.SetLogz();
    c.SetMargin(PAD_LEFT, 0.13f, PAD_BOTTOM, PAD_TOP);
    h.Draw("COLZ");

    // Sovrappongo i mu del fit
    if (!fits.empty()) {
        std::vector<double> xs, ys, exs, eys;
        for (auto& f : fits) {
            if (!f.ok) continue;
            xs.push_back((double)f.n_pe_target);
            ys.push_back(f.mu_tot);
            exs.push_back(0.3);
            eys.push_back(f.sigma_tot);
        }
        if (!xs.empty()) {
            auto* g = new TGraphErrors((int)xs.size(),
                                        xs.data(), ys.data(),
                                        exs.data(), eys.data());
            g->SetMarkerStyle(20); g->SetMarkerColor(kBlack);
            g->SetLineColor(kBlack); g->SetLineWidth(2);
            g->Draw("PZ SAME");
        }
    }
    c.SaveAs(Form("%s/tot_vs_amp_%s_vbias%d.png",
                   outDir.c_str(), tag.c_str(), vbias));
}

// ════════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════════
void sipm_xtalk_tot(int vbias = 55, double frac_pe = 0.60,
                     double cutoff_MHz = 0.0,
                     int nmax = 4, double half_amp = 0.30) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::cout << "\n=================================================\n"
              << "  SiPM XTALK via TOT\n"
              << "  Vbias=" << vbias << "V  let=" << frac_pe << "pe  "
              << "cutoff=" << cutoff_MHz << "MHz\n"
              << "=================================================\n";

    // ── carica calibrazione ──
    CalibResult cal;
    if (!loadCalibration(cal, vbias, cutoff_MHz, DATA_DIR_LOW) || !cal.ok) {
        std::cerr << "[ERR] calibration load failed\n";
        return;
    }
    std::cout << "  gain=" << cal.m << " mV/pe  laser_thr=" << cal.laser_thr << " mV\n";

    // ── path cache esistenti (cerca via prefisso, nome puo' avere suffissi) ──
    std::string pathLo = xtot::findCacheFile(DATA_DIR_LOW, vbias, frac_pe,
                                              cutoff_MHz, cal.laser_thr);
    std::string pathHi = xtot::findCacheFile(DATA_DIR_HIGH, vbias, frac_pe,
                                              cutoff_MHz, cal.laser_thr);
    if (pathLo.empty()) {
        std::cerr << "[ERR] no LOW cache matching prefix in " << DATA_DIR_LOW << "\n";
        return;
    }
    std::cout << "  cache low : " << pathLo << "\n";
    std::cout << "  cache high: " << (pathHi.empty() ? "(not found)" : pathHi) << "\n";

    // ── load ──
    std::vector<xtot::CacheEv> evLo, evHi;
    if (!xtot::loadCache(pathLo, cal, evLo)) {
        std::cerr << "[ERR] cannot load LOW cache. Run sipm_tot_analysis first.\n";
        return;
    }
    bool has_high = (!pathHi.empty()) && xtot::loadCache(pathHi, cal, evHi);
    if (!has_high) {
        std::cerr << "[WARN] HIGH cache not available; only LOW will be analyzed\n";
    }

    // ── slices ──
    auto slices = xtot::defaultSlices(nmax, half_amp);
    std::vector<xtot::SliceFit> fitsLo, fitsHi;

    OutCtx ctx = createOutputDirs("xtalk_tot");

    std::cout << "\n=== SLICE FIT RESULTS ===\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << " npe |   slice (pe)  |              LOW                |              HIGH\n";
    std::cout << "     |               |   N    mu_tot  sig_tot  fL_tail |   N    mu_tot  sig_tot  fL_tail | excess(low-high)\n";
    std::cout << "-----+---------------+---------------------------------+---------------------------------+-----------------\n";

    for (auto& s : slices) {
        auto fL = xtot::fitSlice(evLo, s);
        xtot::SliceFit fH;
        if (has_high) fH = xtot::fitSlice(evHi, s);
        fitsLo.push_back(fL);
        fitsHi.push_back(fH);

        double exc = fL.frac_left - fH.frac_left;
        std::cout << "  " << s.n_pe_target << "  | "
                  << std::setw(5) << s.amp_lo << " - " << std::setw(5) << s.amp_hi << " | "
                  << std::setw(5) << fL.N_in_slice << "  "
                  << std::setw(6) << fL.mu_tot << "   "
                  << std::setw(5) << fL.sigma_tot << "   "
                  << std::setw(5) << fL.frac_left << "  | "
                  << std::setw(5) << fH.N_in_slice << "  "
                  << std::setw(6) << fH.mu_tot << "   "
                  << std::setw(5) << fH.sigma_tot << "   "
                  << std::setw(5) << fH.frac_left << "  |   "
                  << std::setw(7) << exc << "\n";

        // plot per slice
        plotSlice(evLo, evHi, s, fL, fH, ctx.pngDir, vbias);
    }

    // 2D summary
    plotTOTvsAmp(evLo, fitsLo, ctx.pngDir, vbias, "low");
    if (has_high)
        plotTOTvsAmp(evHi, fitsHi, ctx.pngDir, vbias, "high");

    // ── final verdict ──
    std::cout << "\n=== INTERPRETATION ===\n";
    std::cout << "excess(low - high) > 0  on LEFT tail of N-pe peak\n"
              << "  => low light has more events with TOT < mu - 2sig at same amp\n"
              << "  => candidate xtalk signature (paired-cell decay shape)\n\n";
    std::cout << "If excess is consistent across n=1,2,3 pe -> robust xtalk evidence.\n"
              << "If excess ~ 0 or random sign -> no detectable xtalk via TOT.\n";

    std::cout << "\nOutput PNG: " << ctx.pngDir << "\n";
}
