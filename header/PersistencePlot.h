#pragma once
// PersistencePlot.h
// Digital oscilloscope persistence display.
// Fills a 2D histogram (time-t_laser vs amplitude) with every sample of every
// accepted waveform. Result: bright = many waveforms passed through, dark = rare.
//
// RAM: O(t_bins * amp_bins * 8 bytes) ≈ 5 MB — no raw waveform buffer.
//
// Usage:
//   PersistencePlotter persist(let_thr);
//   // inside event loop: persist.collect(v_t, af, t_laser, t_rise, let_thr);
//   // after loop:        persist.draw("vbias55_let0.50pe", ctx);

#include "Config.h"
#include "OutputManager.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>

#include <TCanvas.h>
#include <TH2D.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TAxis.h>

class PersistencePlotter {
public:
    explicit PersistencePlotter(double let_thr  = 10.0,
                                double t_lo     = -50.0,
                                double t_hi     = 200.0,
                                int    t_bins   = 1024,
                                double amp_lo   = -20.0,
                                double amp_hi   = 200.0,
                                int    amp_bins = 440)
        : let_thr_(let_thr),
          t_lo_(t_lo), t_hi_(t_hi), t_bins_(t_bins),
          amp_lo_(amp_lo), amp_hi_(amp_hi), amp_bins_(amp_bins)
    {}

    ~PersistencePlotter() { delete h_; }

    // Call for every accepted event. Fills the TH2D sample-by-sample.
    // RAM cost: O(1) per event.
    void collect(const std::vector<double>& time_v,
                 const std::vector<double>& amp_v,
                 double t_laser,
                 double t_rise,
                 double /*let_thr_unused*/)
    {
        if (time_v.empty()) return;

        // Lazy histogram creation on the first call.
        // Using (void*)this in the name guarantees uniqueness across
        // multiple PersistencePlotter instances in the same ROOT session.
        if (!h_) {
            h_ = new TH2D(
                Form("hPersist_%p", (void*)this),
                "",
                t_bins_,   t_lo_,   t_hi_,
                amp_bins_, amp_lo_, amp_hi_);
            h_->SetDirectory(nullptr);
        }

        // Fill every sample, shifted to t_laser = 0
        int ns = (int)std::min(time_v.size(), amp_v.size());
        for (int j = 0; j < ns; ++j)
            h_->Fill(time_v[j] - t_laser, amp_v[j]);

        // Store only Δt (one double per event — negligible RAM)
        delta_t_vec_.push_back(t_rise - t_laser);
        n_waveforms_++;
    }

    int nWaveforms() const { return n_waveforms_; }

    // Draw the persistence map and save PNG.
    void draw(const std::string& tag, OutCtx& ctx) const
    {
        if (n_waveforms_ == 0 || !h_) {
            std::cout << "  [Persistence] No waveforms for " << tag << "\n";
            return;
        }

        std::cout << "  [Persistence] Drawing map: "
                  << n_waveforms_ << " waveforms  tag=" << tag << "\n";

        // Median and RMS of Δt
        std::vector<double> dts = delta_t_vec_;
        std::sort(dts.begin(), dts.end());
        double median_dt = dts[dts.size() / 2];
        double mean_dt   = 0.0;
        for (double d : dts) mean_dt += d;
        mean_dt /= (double)dts.size();
        double rms_dt = 0.0;
        for (double d : dts) rms_dt += (d - mean_dt) * (d - mean_dt);
        rms_dt = std::sqrt(rms_dt / (double)dts.size());

        h_->SetTitle(
            Form("Waveform persistence (aligned on t_{laser})   %s"
                 ";t #minus t_{laser} (ns);Amplitude (mV)", tag.c_str()));

        TCanvas* cP = new TCanvas(
            Form("cPersist_%s", tag.c_str()),
            Form("Persistence  %s", tag.c_str()),
            1100, 600);
        cP->SetLeftMargin(PAD_LEFT);
        cP->SetRightMargin(0.13f);
        cP->SetBottomMargin(PAD_BOTTOM);
        cP->SetTopMargin(PAD_TOP);
        cP->SetGrid();

        gStyle->SetPalette(kBird);
        h_->SetContour(100);
        h_->GetXaxis()->SetTitleSize(0.048f);
        h_->GetYaxis()->SetTitleSize(0.048f);
        h_->GetXaxis()->SetLabelSize(0.040f);
        h_->GetYaxis()->SetLabelSize(0.040f);
        h_->GetZaxis()->SetTitle("Counts / bin");
        h_->GetZaxis()->SetTitleSize(0.040f);
        h_->GetZaxis()->SetLabelSize(0.036f);
        cP->SetLogz(1);
        h_->Draw("COLZ");

        // Vertical line at t=0: laser trigger  (green dashed)
        TLine* lLas = new TLine(0.0, amp_lo_, 0.0, amp_hi_);
        lLas->SetLineColor(kGreen+2);
        lLas->SetLineStyle(2);
        lLas->SetLineWidth(2);
        lLas->Draw("same");

        // Vertical line at median Δt: median t_rise  (red dashed)
        TLine* lRise = new TLine(median_dt, amp_lo_, median_dt, amp_hi_);
        lRise->SetLineColor(kRed+1);
        lRise->SetLineStyle(2);
        lRise->SetLineWidth(2);
        lRise->Draw("same");

        // Horizontal line at LET threshold  (blue dashed)
        TLine* lThr = new TLine(t_lo_, let_thr_, t_hi_, let_thr_);
        lThr->SetLineColor(kAzure+1);
        lThr->SetLineStyle(2);
        lThr->SetLineWidth(2);
        lThr->Draw("same");

        // Stats box
        TPaveText* pt = new TPaveText(0.15, 0.70, 0.48, 0.88, "NDC");
        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->SetFillStyle(1001);
        pt->SetTextFont(42);
        pt->SetTextSize(0.036);
        pt->AddText(Form("N waveforms = %d",          n_waveforms_));
        pt->AddText(Form("Median #Deltat = %.2f ns",  median_dt));
        pt->AddText(Form("RMS #Deltat    = %.2f ns",  rms_dt));
        pt->AddText(Form("LET thr = %.1f mV  (blue)", let_thr_));
        pt->Draw();

        cP->Update();
        cP->Modified();
        ctx.savePNG(cP, Form("persistence_%s.png", tag.c_str()));
        // cP->Close();   // FIX: canvas resta aperta
        // delete cP;     // FIX
        delete lLas;
        delete lRise;
        delete lThr;
        delete pt;

        std::cout << "  [Persistence] Done."
                  << "  N=" << n_waveforms_
                  << "  median Δt=" << std::fixed << std::setprecision(2)
                  << median_dt << " ns"
                  << "  RMS=" << rms_dt << " ns\n";
    }

    void clear() {
        delete h_;
        h_ = nullptr;
        n_waveforms_ = 0;
        delta_t_vec_.clear();
    }

private:
    // Histogram parameters
    double let_thr_;
    double t_lo_,   t_hi_;    int t_bins_;
    double amp_lo_, amp_hi_;  int amp_bins_;

    // State
    TH2D*               h_           = nullptr;   // the persistence map
    int                 n_waveforms_ = 0;
    std::vector<double> delta_t_vec_;             // one entry per event
};
