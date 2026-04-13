#pragma once
// WaveformPlotter.h
// Collects and plots filtered waveforms for events with small TOT.
// Purpose: visually inspect signal shapes of events that barely crossed
// the LET threshold — purest single-p.e. laser events, least time walk.
//
// Usage:
//   WaveformPlotter wfp(10.0, 100);   // tot_cut=10 ns, max 100 waveforms
//   // inside event loop:
//   if (wfp.collecting()) wfp.collect(i, v_t, af, t_laser, t_rise, t_fall, let_thr, n_pe);
//   // after loop:
//   wfp.draw("vbias55_let0.50pe", ctx);
//
// Canvas: 5×4 grid per page, X zoomed to [t_laser-20, t_rise+60] ns.
// Markers: green=t_laser, red=t_rise, orange=t_fall, blue=LET threshold.

#include "Config.h"
#include "OutputManager.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>

#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TLine.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TStyle.h>

// One stored waveform
struct WaveformRecord {
    std::vector<double> time;    // sample times [ns]
    std::vector<double> amp;     // filtered amplitude [mV]
    double t_laser  = 0.0;       // laser trigger time [ns]
    double t_rise   = 0.0;       // LET threshold crossing [ns]
    double t_fall   = 0.0;       // falling edge [ns]
    double tot      = 0.0;       // TOT = t_fall - t_rise [ns]
    double delta_t  = 0.0;       // t_rise - t_laser [ns]
    double let_thr  = 0.0;       // LET threshold level [mV]
    int    n_pe     = 0;
    Long64_t entry  = -1;        // TChain entry index
};

class WaveformPlotter {
public:
    // tot_cut:       collect only events with TOT < tot_cut [ns]
    // dt_lo/dt_hi:   optional delta_t window filter (default: disabled)
    // max_waveforms: stop collecting after this many events
    explicit WaveformPlotter(double tot_cut       = 10.0,
                             int    max_waveforms = 100,
                             double dt_lo         = -1e9,
                             double dt_hi         = +1e9)
        : tot_cut_(tot_cut), max_waveforms_(max_waveforms),
          dt_lo_(dt_lo), dt_hi_(dt_hi) {}

    // Returns true if still collecting (not yet full)
    bool collecting() const {
        return (int)records_.size() < max_waveforms_;
    }

    // Call during the event loop for every accepted event.
    // Only stores the waveform if TOT < tot_cut and not yet full.
    void collect(Long64_t            entry,
                 const std::vector<double>& time_v,
                 const std::vector<double>& amp_filtered,
                 double t_laser,
                 double t_rise,
                 double t_fall,
                 double let_thr,
                 int    n_pe)
    {
        double tot     = t_fall - t_rise;
        double delta_t = t_rise - t_laser;
        if (tot <= 0.0 || tot >= tot_cut_) return;
        if (delta_t < dt_lo_ || delta_t > dt_hi_) return;
        if (!collecting()) return;

        WaveformRecord rec;
        rec.time    = time_v;
        rec.amp     = amp_filtered;
        rec.t_laser = t_laser;
        rec.t_rise  = t_rise;
        rec.t_fall  = t_fall;
        rec.tot     = tot;
        rec.delta_t = delta_t;
        rec.let_thr = let_thr;
        rec.n_pe    = n_pe;
        rec.entry   = entry;
        records_.push_back(std::move(rec));
    }

    int size() const { return (int)records_.size(); }

    // Draw all collected waveforms and save PNGs.
    // tag:      label used in canvas/file names (e.g. "vbias55_let0.50pe")
    // ctx:      OutputManager context for numbered PNG saving
    void draw(const std::string& tag, OutCtx& ctx) const
    {
        if (records_.empty()) {
            std::cout << "  [WFPlotter] No waveforms collected for " << tag << "\n";
            return;
        }

        const int COLS  = 5;
        const int ROWS  = 4;
        const int PER_PAGE = COLS * ROWS;  // 20 per canvas
        int n_pages = ((int)records_.size() + PER_PAGE - 1) / PER_PAGE;

        std::cout << "  [WFPlotter] Drawing " << records_.size()
                  << " waveforms (TOT < " << tot_cut_ << " ns)  "
                  << n_pages << " page(s)  tag=" << tag << "\n";

        for (int pg = 0; pg < n_pages; ++pg) {
            int first = pg * PER_PAGE;
            int last  = std::min(first + PER_PAGE, (int)records_.size());
            int count = last - first;

            TCanvas* cWF = new TCanvas(
                Form("cWF_%s_pg%d", tag.c_str(), pg),
                Form("Waveforms TOT < %.0f ns -- %s (page %d/%d)",
                     tot_cut_, tag.c_str(), pg+1, n_pages),
                300 * COLS, 220 * ROWS);
            cWF->SetFillColor(0);
            cWF->Divide(COLS, ROWS, 0.002, 0.002);

            for (int k = 0; k < count; ++k) {
                const WaveformRecord& rec = records_[first + k];
                cWF->cd(k + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.03);
                gPad->SetBottomMargin(0.18);
                gPad->SetTopMargin(0.15);
                gPad->SetGrid();
                gPad->SetFillColor(0);

                // Zoom window: 20 ns before laser, 60 ns after t_rise
                double x_lo = rec.t_laser - 20.0;
                double x_hi = rec.t_rise  + 60.0;
                // Clamp to waveform range
                if (!rec.time.empty()) {
                    x_lo = std::max(x_lo, rec.time.front());
                    x_hi = std::min(x_hi, rec.time.back());
                }

                // Y range: -5 mV to 1.3 * max amplitude in zoom window
                double amp_max_zoom = 0.0;
                for (size_t j = 0; j < rec.time.size(); ++j)
                    if (rec.time[j] >= x_lo && rec.time[j] <= x_hi)
                        if (rec.amp[j] > amp_max_zoom)
                            amp_max_zoom = rec.amp[j];
                amp_max_zoom = std::max(amp_max_zoom, rec.let_thr * 2.0);
                double y_lo = -5.0;
                double y_hi = amp_max_zoom * 1.35;

                // Build TGraph for the zoomed region
                std::vector<double> gx, gy;
                gx.reserve(rec.time.size());
                gy.reserve(rec.amp.size());
                for (size_t j = 0; j < rec.time.size(); ++j) {
                    if (rec.time[j] >= x_lo && rec.time[j] <= x_hi) {
                        gx.push_back(rec.time[j]);
                        gy.push_back(rec.amp[j]);
                    }
                }

                TGraph* gr = new TGraph((int)gx.size(), gx.data(), gy.data());
                gr->SetTitle("");
                gr->SetLineColor(kAzure+1);
                gr->SetLineWidth(1);
                gr->GetXaxis()->SetTitle("t (ns)");
                gr->GetYaxis()->SetTitle("A (mV)");
                gr->GetXaxis()->SetTitleSize(0.08);
                gr->GetYaxis()->SetTitleSize(0.08);
                gr->GetXaxis()->SetLabelSize(0.07);
                gr->GetYaxis()->SetLabelSize(0.07);
                gr->GetXaxis()->SetTitleOffset(0.9);
                gr->GetYaxis()->SetTitleOffset(0.85);
                gr->GetXaxis()->SetLimits(x_lo, x_hi);
                gr->GetYaxis()->SetRangeUser(y_lo, y_hi);
                gr->Draw("AL");

                TLine* lThr = new TLine(x_lo, rec.let_thr, x_hi, rec.let_thr);
                lThr->SetLineColor(kBlue); lThr->SetLineStyle(2); lThr->SetLineWidth(1);
                lThr->Draw("same");

                TLine* lLas = nullptr;
                if (rec.t_laser > x_lo && rec.t_laser < x_hi) {
                    lLas = new TLine(rec.t_laser, y_lo,
                                     rec.t_laser, y_hi * 0.85);
                    lLas->SetLineColor(kGreen+2); lLas->SetLineStyle(2);
                    lLas->SetLineWidth(1); lLas->Draw("same");
                }

                TLine* lRise = nullptr;
                if (rec.t_rise > x_lo && rec.t_rise < x_hi) {
                    lRise = new TLine(rec.t_rise, y_lo,
                                      rec.t_rise, y_hi * 0.85);
                    lRise->SetLineColor(kRed+1); lRise->SetLineStyle(2);
                    lRise->SetLineWidth(1); lRise->Draw("same");
                }

                TLine* lFall = nullptr;
                if (rec.t_fall > x_lo && rec.t_fall < x_hi) {
                    lFall = new TLine(rec.t_fall, y_lo,
                                      rec.t_fall, y_hi * 0.7);
                    lFall->SetLineColor(kOrange+7); lFall->SetLineStyle(2);
                    lFall->SetLineWidth(1); lFall->Draw("same");
                }

                TPaveText* pt = new TPaveText(0.01, 0.86, 0.99, 0.99, "NDC");
                pt->SetBorderSize(0); pt->SetFillColor(0); pt->SetFillStyle(0);
                pt->SetTextFont(42); pt->SetTextSize(0.075);
                pt->SetTextAlign(22);
                pt->AddText(Form("ev %lld  TOT=%.1f ns  #Deltat=%.1f ns  %dpe",
                                 rec.entry, rec.tot, rec.delta_t, rec.n_pe));
                pt->Draw();

                delete pt;
                delete lThr;
                if (lLas) delete lLas;
                if (lRise) delete lRise;
                if (lFall) delete lFall;
                delete gr;
            }

            // Fill empty pads with blank
            for (int k = count; k < PER_PAGE; ++k) {
                cWF->cd(k + 1);
                gPad->SetFillColor(kWhite);
            }

            cWF->cd();
            cWF->Update(); cWF->Modified();
            ctx.savePNG(cWF, Form("waveforms_smallTOT_%s_p%d.png",
                                  tag.c_str(), pg+1));
            // cWF->Close();   // FIX: canvas resta aperta
            // delete cWF;     // FIX
        }

        std::cout << "  [WFPlotter] Done. " << n_pages
                  << " PNG(s) saved.\n";
    }

    void clear() { records_.clear(); }

private:
    double tot_cut_;
    int    max_waveforms_;
    double dt_lo_;
    double dt_hi_;
    std::vector<WaveformRecord> records_;
};
