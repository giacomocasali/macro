#pragma once
// TOTBandPersistence.h
// Persistence plot of waveforms whose TOT falls in a given [tot_lo, tot_hi] band.
//
// Usage (inside the event loop of VbiasAnalysis_v2 or any other loop):
//
//   TOTBandPersistence p1(40.0, 43.0, let_thr);   // TOT in [40,43] ns
//   TOTBandPersistence p2(50.0, 60.0, let_thr);   // TOT in [50,60] ns
//
//   // inside event loop (after t_rise, t_fall, t_laser are computed):
//   p1.collect(v_t, af, t_laser, t_rise, t_fall);
//   p2.collect(v_t, af, t_laser, t_rise, t_fall);
//
//   // after loop:
//   p1.draw("vbias55_let0.50pe_tot40_43", ctx);
//   p2.draw("vbias55_let0.50pe_tot50_60", ctx);
//
// RAM: O(t_bins * amp_bins * 8) ≈ 5 MB per instance — no waveform buffering.

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

class TOTBandPersistence {
public:
    // tot_lo, tot_hi : TOT band [ns] — only waveforms with tot_lo <= TOT < tot_hi are accepted
    // let_thr        : LET threshold [mV] — drawn as horizontal dashed line
    // t_lo/t_hi      : time axis range relative to t_laser [ns]
    // amp_lo/amp_hi  : amplitude axis range [mV]
    explicit TOTBandPersistence(double tot_lo,
                                double tot_hi,
                                double let_thr  = 10.0,
                                double t_lo     = -50.0,
                                double t_hi     = 200.0,
                                int    t_bins   = 1024,
                                double amp_lo   = -20.0,
                                double amp_hi   = 200.0,
                                int    amp_bins = 440)
        : tot_lo_(tot_lo), tot_hi_(tot_hi), let_thr_(let_thr),
          t_lo_(t_lo), t_hi_(t_hi), t_bins_(t_bins),
          amp_lo_(amp_lo), amp_hi_(amp_hi), amp_bins_(amp_bins)
    {}

    ~TOTBandPersistence() { delete h_; }

    // Call for every accepted event.
    // v_t, amp_v : full waveform vectors (time [ns], amplitude [mV])
    // t_laser    : laser trigger time [ns]
    // t_rise     : LET threshold crossing time [ns]
    // t_fall     : falling edge time [ns]
    // Returns true if the event was accepted (TOT in band).
    bool collect(const std::vector<double>& v_t,
                 const std::vector<double>& amp_v,
                 double t_laser,
                 double t_rise,
                 double t_fall)
    {
        double tot = t_fall - t_rise;
        if (tot < tot_lo_ || tot >= tot_hi_) return false;

        // Lazy histogram creation
        if (!h_) {
            h_ = new TH2D(
                Form("hTOTBandPers_%.0f_%.0f_%p", tot_lo_, tot_hi_, (void*)this),
                "",
                t_bins_,   t_lo_,   t_hi_,
                amp_bins_, amp_lo_, amp_hi_);
            h_->SetDirectory(nullptr);
        }

        int ns = (int)std::min(v_t.size(), amp_v.size());
        for (int j = 0; j < ns; ++j)
            h_->Fill(v_t[j] - t_laser, amp_v[j]);

        delta_t_vec_.push_back(t_rise - t_laser);
        tot_vec_.push_back(tot);
        n_waveforms_++;
        return true;
    }

    int nWaveforms() const { return n_waveforms_; }

    // Draw and save PNG.
    // tag : label used in canvas/file names (e.g. "vbias55_let0.50pe")
    void draw(const std::string& tag, OutCtx& ctx) const
    {
        if (n_waveforms_ == 0 || !h_) {
            std::cout << "  [TOTBandPers] No waveforms for TOT in ["
                      << tot_lo_ << ", " << tot_hi_ << ") ns  tag=" << tag << "\n";
            return;
        }

        std::cout << "  [TOTBandPers] Drawing: " << n_waveforms_
                  << " waveforms  TOT in [" << tot_lo_ << ", " << tot_hi_ << ") ns"
                  << "  tag=" << tag << "\n";

        // Δt stats
        std::vector<double> dts = delta_t_vec_;
        std::sort(dts.begin(), dts.end());
        double median_dt = dts[dts.size() / 2];
        double mean_dt = 0; for (double d : dts) mean_dt += d; mean_dt /= dts.size();
        double rms_dt = 0;  for (double d : dts) rms_dt += (d-mean_dt)*(d-mean_dt);
        rms_dt = std::sqrt(rms_dt / dts.size());

        // TOT stats
        std::vector<double> tots = tot_vec_;
        std::sort(tots.begin(), tots.end());
        double mean_tot = 0; for (double t : tots) mean_tot += t; mean_tot /= tots.size();
        double rms_tot = 0;  for (double t : tots) rms_tot += (t-mean_tot)*(t-mean_tot);
        rms_tot = std::sqrt(rms_tot / tots.size());

        h_->SetTitle(
            Form("Persistence  TOT #in [%.0f, %.0f) ns   %s"
                 ";t #minus t_{laser} (ns);Amplitude (mV)",
                 tot_lo_, tot_hi_, tag.c_str()));

        TCanvas* c = new TCanvas(
            Form("cTOTBandPers_%.0f_%.0f_%s", tot_lo_, tot_hi_, tag.c_str()),
            Form("Persistence TOT [%.0f,%.0f) ns  %s", tot_lo_, tot_hi_, tag.c_str()),
            1100, 600);
        c->SetLeftMargin(PAD_LEFT);
        c->SetRightMargin(0.13f);
        c->SetBottomMargin(PAD_BOTTOM);
        c->SetTopMargin(PAD_TOP);
        c->SetGrid();
        c->SetLogz(1);

        gStyle->SetPalette(kBird);
        h_->SetContour(100);
        h_->GetXaxis()->SetTitleSize(0.048f);
        h_->GetYaxis()->SetTitleSize(0.048f);
        h_->GetXaxis()->SetLabelSize(0.040f);
        h_->GetYaxis()->SetLabelSize(0.040f);
        h_->GetZaxis()->SetTitle("Counts / bin");
        h_->GetZaxis()->SetTitleSize(0.040f);
        h_->GetZaxis()->SetLabelSize(0.036f);
        h_->Draw("COLZ");

        // t_laser = 0 (green dashed)
        TLine* lLas = new TLine(0.0, amp_lo_, 0.0, amp_hi_);
        lLas->SetLineColor(kGreen+2); lLas->SetLineStyle(2); lLas->SetLineWidth(2);
        lLas->Draw("same");

        // median t_rise (red dashed)
        TLine* lRise = new TLine(median_dt, amp_lo_, median_dt, amp_hi_);
        lRise->SetLineColor(kRed+1); lRise->SetLineStyle(2); lRise->SetLineWidth(2);
        lRise->Draw("same");

        // LET threshold (blue dashed)
        TLine* lThr = new TLine(t_lo_, let_thr_, t_hi_, let_thr_);
        lThr->SetLineColor(kAzure+1); lThr->SetLineStyle(2); lThr->SetLineWidth(2);
        lThr->Draw("same");

        TPaveText* pt = new TPaveText(0.15, 0.62, 0.50, 0.88, "NDC");
        pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
        pt->SetTextFont(42);  pt->SetTextSize(0.034);
        pt->AddText(Form("N waveforms = %d",            n_waveforms_));
        pt->AddText(Form("TOT band: [%.0f, %.0f) ns",   tot_lo_, tot_hi_));
        pt->AddText(Form("Mean TOT = %.2f #pm %.2f ns", mean_tot, rms_tot));
        pt->AddText(Form("Median #Deltat = %.2f ns",    median_dt));
        pt->AddText(Form("RMS #Deltat    = %.2f ns",    rms_dt));
        pt->Draw();

        c->Update(); c->Modified();
        ctx.savePNG(c, Form("persistence_tot%.0f_%.0f_%s.png",
                            tot_lo_, tot_hi_, tag.c_str()));
        delete lLas; delete lRise; delete lThr; delete pt;
    }

    void clear() {
        delete h_; h_ = nullptr;
        n_waveforms_ = 0;
        delta_t_vec_.clear();
        tot_vec_.clear();
    }

private:
    double tot_lo_, tot_hi_;
    double let_thr_;
    double t_lo_, t_hi_;    int t_bins_;
    double amp_lo_, amp_hi_; int amp_bins_;

    TH2D*               h_           = nullptr;
    int                 n_waveforms_ = 0;
    std::vector<double> delta_t_vec_;
    std::vector<double> tot_vec_;
};
