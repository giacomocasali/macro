// final.cpp  —  v4.0
// Automatically detects all data.vbias_{V}_L.root files in the current directory.
//
// Changes from v3.0:
//   [NEW-1]  Baseline window and trigger window are now INDEPENDENT.
//            The operator specifies:
//              - baseline window [t_base_start, t_base_end] ns: region used to
//                compute the median baseline offset. Must be noise-only (pre-signal).
//              - trigger window  [t_trig_start, t_trig_end] ns: region where the
//                discriminator runs. Can be [0,45], [45,55], [0,200] or anything.
//            This fixes the bug where t_trig_start=0 left zero samples for
//            baseline estimation, causing correctBaseline() to fall back to
//            the uncorrected waveform.
//   [NEW-2]  findFirstPeak() rewritten: uses a heavier 5-point smoothing kernel
//            and searches for the true peak maximum (not the rising edge), fixing
//            the misidentification on steep narrow peaks (e.g. lux=3).
//   [NEW-3]  All output text is in English.
//
// Inherited fixes from v2.0 / v3.0:
//   x_err=0, step_der=1 mV, min_positive_thr=5 mV, sigma limits on fit,
//   adaptive fit window, chi2/ndf warnings, TLine bug fix, uniform pad
//   margins via explicit TPad, negative threshold support.

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
#include <TPad.h>
#include <TPaveText.h>
#include "../header/ButterworthFilter.h"  // scipy-equivalent Butterworth IIR
#include <TLegend.h>
#include <TH1D.h>

// ════════════════════════════════════════════════════════════
//  GLOBAL PARAMETERS  —  edit here to change behaviour
// ════════════════════════════════════════════════════════════
static const double MIN_THR          = -10.0;  // minimum threshold (mV)
static const double MAX_THR          = 110.0;  // maximum threshold (mV)
static const double STEP_SCAN        =   0.1;  // threshold scan step (mV)
static const double STEP_DER         =   0.1;  // derivative step for fit/peak finding (mV)
static const double STEP_DER_PLOT    =   0.5;  // derivative step for display only (mV)
static const double MIN_POSITIVE_THR =   5.0;  // spike-rejection threshold (mV)
                                                // fit and Y-zoom only above this
static const double FIT_SIGMA_MAX    =  15.0;  // upper limit on Gaussian sigma (mV)
static const double FIT_SIGMA_MIN    =   0.1;  // lower limit on Gaussian sigma (mV)

// Baseline window: always the first BASELINE_END ns of the waveform.
// This is a quiet noise-only region, fixed by hardware — do not change
// unless your acquisition window changes.
static const double BASELINE_START   =   0.0;  // ns
static const double BASELINE_END     =  30.0;  // ns

// Pad margins — identical for both pads → same physical X scale
static const float  PAD_LEFT   = 0.13f;
static const float  PAD_RIGHT  = 0.04f;
static const float  PAD_TOP    = 0.10f;
static const float  PAD_BOTTOM = 0.13f;

// ════════════════════════════════════════════════════════════
//  BASELINE CORRECTION  (median over baseline window)
//  Uses samples with  t_base_start <= t < t_base_end
// ════════════════════════════════════════════════════════════
static std::vector<double> correctBaseline(const std::vector<double>& time,
                                           const std::vector<double>& amp,
                                           double t_base_start,
                                           double t_base_end) {
    std::vector<double> pre;
    for (size_t i = 0; i < time.size(); ++i)
        if (time[i] >= t_base_start && time[i] < t_base_end)
            pre.push_back(amp[i]);

    if (pre.empty()) {
        std::cerr << "  [WARN] No samples found in baseline window ["
                  << t_base_start << ", " << t_base_end
                  << ") ns — baseline NOT corrected.\n";
        return amp;
    }
    std::vector<double> tmp = pre;
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
    double offset = tmp[tmp.size()/2];
    std::vector<double> out = amp;
    for (auto& a : out) a -= offset;
    return out;
}

// ════════════════════════════════════════════════════════════
//  4th-ORDER BUTTERWORTH LOW-PASS FILTER  (IIR, biquad cascade)
// ════════════════════════════════════════════════════════════
// butterworthLowPass is provided by ButterworthFilter.h
// To use zero-phase (sosfiltfilt): butterworthLowPassZP(x, fc, fs, 4)
// To use fast pre-computed 500MHz:  butterworthLP_500MHz(x)

// ════════════════════════════════════════════════════════════
//  FILE PARSING:  data.vbias_{55}_3.root  →  vbias=55, lux=3
// ════════════════════════════════════════════════════════════
struct FileInfo {
    std::string path;
    double      vbias;
    double      lux;
    std::string tag;   // e.g. "vbias55_lux3"
};

static bool parseFilename(const std::string& path, FileInfo& info) {
    size_t slash = path.find_last_of("/\\");
    std::string fname = (slash == std::string::npos) ? path : path.substr(slash + 1);
    const std::string prefix = "data.vbias_{";
    size_t p0 = fname.find(prefix);
    if (p0 == std::string::npos) return false;
    size_t vStart = p0 + prefix.size();
    size_t vEnd   = fname.find('}', vStart);
    if (vEnd == std::string::npos) return false;
    std::string vStr = fname.substr(vStart, vEnd - vStart);
    size_t lStart = vEnd + 2;
    size_t lEnd   = fname.find(".root", lStart);
    if (lEnd == std::string::npos) return false;
    std::string lStr = fname.substr(lStart, lEnd - lStart);
    try { info.vbias = std::stod(vStr); info.lux = std::stod(lStr); }
    catch (...) { return false; }
    info.path = path;
    info.tag  = "vbias" + vStr + "_lux" + lStr;
    return true;
}

// ════════════════════════════════════════════════════════════
//  FIND FIRST PHYSICAL PEAK of -dN/dV  (1st p.e.)
//  in the region thr > MIN_POSITIVE_THR.
//
//  Algorithm:
//   1. 5-point smoothing (suppresses numerical noise on steep flanks).
//   2. Find global maximum of smoothed derivative (scale reference).
//   3. Scan left-to-right looking for a strict local maximum that:
//        a) exceeds 10% of the global maximum
//        b) was preceded by a genuine rising flank of at least
//           MIN_RISING_POINTS consecutive increasing samples.
//      This requirement (b) prevents latching onto a descending
//      front, which locally looks like a maximum but has no
//      rising flank — exactly the failure mode seen with narrow
//      trigger windows where the scan starts mid-pulse.
//   4. Fallback: absolute maximum of raw derivative (no flank check).
// ════════════════════════════════════════════════════════════
// Minimum rising flank width in mV before accepting a peak.
// With STEP_DER=0.1 mV this means 20 consecutive rising points (2 mV).
// Expressed in mV so it stays meaningful if STEP_DER changes.
static const double MIN_RISING_MV = 2.0;

static int findFirstPeak(const std::vector<double>& thr,
                         const std::vector<double>& der) {
    int n = (int)thr.size();
    // Convert mV threshold to number of derivative points
    int minRisingPts = std::max(2, (int)std::round(MIN_RISING_MV / STEP_DER));
    if (n < minRisingPts + 2) return -1;

    // 5-point smoothing
    std::vector<double> sm(n, 0.0);
    sm[0]   = der[0];
    sm[1]   = (der[0] + der[1] + der[2]) / 3.0;
    sm[n-2] = (der[n-3] + der[n-2] + der[n-1]) / 3.0;
    sm[n-1] = der[n-1];
    for (int i = 2; i < n - 2; ++i)
        sm[i] = (der[i-2] + der[i-1] + der[i] + der[i+1] + der[i+2]) / 5.0;

    // Global maximum in the physical region
    double globalMax = 0;
    for (int i = 0; i < n; ++i)
        if (thr[i] > MIN_POSITIVE_THR && sm[i] > globalMax) globalMax = sm[i];
    if (globalMax <= 0) return -1;

    const double minFraction = 0.10;

    // Count how many consecutive rising steps precede position i
    // A step is "rising" if sm[i] > sm[i-1]
    std::vector<int> risingBefore(n, 0);
    for (int i = 1; i < n; ++i)
        risingBefore[i] = (sm[i] > sm[i-1]) ? risingBefore[i-1] + 1 : 0;

    // First strict local maximum with sufficient rising flank
    int iPeak = -1;
    for (int i = 1; i < n - 1; ++i) {
        if (thr[i] <= MIN_POSITIVE_THR) continue;
        if (sm[i] < globalMax * minFraction) continue;
        // Strict local maximum: higher than both neighbours
        if (sm[i] <= sm[i-1] || sm[i] <= sm[i+1]) continue;
        // Must have been rising for at least minRisingPts steps (= MIN_RISING_MV mV)
        if (risingBefore[i] < minRisingPts) continue;
        iPeak = i;
        break;
    }

    // Fallback: absolute maximum of raw derivative (no flank check)
    if (iPeak < 0) {
        std::cout << "  [WARN] findFirstPeak: no peak with "
                  << minRisingPts << " rising points (" << MIN_RISING_MV
                  << " mV) found — falling back to absolute maximum.\n";
        for (int i = 0; i < n; ++i)
            if (thr[i] > MIN_POSITIVE_THR && (iPeak < 0 || der[i] > der[iPeak]))
                iPeak = i;
    }
    return iPeak;
}

// ════════════════════════════════════════════════════════════
//  LASER TRIGGER TIME
//  First crossing above laser_thr on the laser channel,
//  with linear interpolation for sub-sample precision.
//  Returns -999 if no crossing is found.
//  Laser baseline is corrected over [0, BASELINE_END) ns
//  using the same median estimator as the SiPM channel.
// ════════════════════════════════════════════════════════════
static double laserTriggerTime(const Double_t* t_laser,
                               const Double_t* a_laser,
                               int N,
                               double laser_thr = 10.0) {
    // Median baseline over first BASELINE_END ns
    std::vector<double> pre;
    for (int j = 0; j < N; ++j)
        if (t_laser[j] < BASELINE_END) pre.push_back(a_laser[j]);
    double offset = 0.0;
    if (!pre.empty()) {
        std::vector<double> tmp = pre;
        std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
        offset = tmp[tmp.size()/2];
    }
    for (int j = 1; j < N; ++j) {
        double v0 = a_laser[j-1] - offset;
        double v1 = a_laser[j]   - offset;
        if (v1 > laser_thr && v0 <= laser_thr)
            return t_laser[j-1]
                   + (laser_thr - v0) * (t_laser[j] - t_laser[j-1]) / (v1 - v0);
    }
    return -999.0;
}

// ════════════════════════════════════════════════════════════
//  LASER TIMING CANVAS
//
//  Reads the "laser" TTree from the input file (if present)
//  and produces a dedicated canvas with two pads:
//
//  Top pad — TH1D of  (t_peak_SiPM − t_laser_trigger)
//    The peak time is the time of the sample with the maximum
//    amplitude in the FULL waveform (no trigger-window cut),
//    after baseline correction and LP filtering.
//    This tells you WHERE the SiPM signal falls in absolute
//    time relative to the laser pulse.
//    → Use this to confirm that [0,45] ns is truly dark and
//      that the signal window [45,55] ns is correctly placed.
//
//  Bottom pad — TH1D of the peak amplitude for events where
//    the peak falls inside the trigger window [t_trig_start,
//    t_trig_end].  This is the amplitude spectrum of the
//    selected events, useful as a cross-check of the scan.
//
//  If the "laser" tree is absent the function prints a warning
//  and returns without drawing anything.
// ════════════════════════════════════════════════════════════
static void makeLaserCanvas(TFile*             rawFile,
                            const FileInfo&    info,
                            double             cutoff_MHz,
                            double             t_trig_start,
                            double             t_trig_end,
                            const std::string& tag) {

    TTree* treeLaser = (TTree*)rawFile->Get("laser");
    if (!treeLaser) {
        std::cout << "  [INFO] No 'laser' tree found in "
                  << info.path << " — timing canvas skipped.\n";
        return;
    }
    TTree* treeCh1 = (TTree*)rawFile->Get("ch1");
    if (!treeCh1) return;   // should never happen here

    const int N = 1024;
    Double_t t1[N], a1[N], tL[N], aL[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);
    treeLaser->SetBranchAddress("time",      tL);
    treeLaser->SetBranchAddress("amplitude", aL);
    treeCh1->GetEntry(0);   // read time axis for fs
    double fs_MHz = 1000.0 / (t1[1] - t1[0]);
    double dt_ns  = t1[1] - t1[0];
    double eps    = 0.5 * dt_ns;

    // Locate trigger-window indices (same epsilon logic as main loop)
    int j_trig_start = 0, j_trig_end = N - 1;
    for (int j = 0; j < N; ++j)
        if (t1[j] >= t_trig_start - eps) { j_trig_start = j; break; }
    for (int j = N - 1; j >= 0; --j)
        if (t1[j] <= t_trig_end + eps)   { j_trig_end   = j; break; }

    // Histogram: t_peak_SiPM − t_laser over the full waveform
    // Range covers the whole waveform width with 0.4 ns bins
    double wf_duration = t1[N-1] - t1[0];
    int    timing_bins = (int)(wf_duration / 0.4) + 1;
    TH1D* hTiming = new TH1D(
        ("hTiming_" + tag).c_str(),
        Form("SiPM peak time relative to laser   V_{bias}=%.0fV, lux=%.0f"
             ";t_{peak}^{SiPM} - t_{laser}  (ns);Events",
             info.vbias, info.lux),
        timing_bins, t1[0], t1[N-1]);
    hTiming->SetDirectory(nullptr);   // detach from rawFile — survives Close()
    hTiming->SetLineColor(kAzure+1);
    hTiming->SetLineWidth(2);
    hTiming->SetFillColorAlpha(kAzure+1, 0.25);

    // Histogram: peak amplitude for events inside trigger window
    TH1D* hAmpTrig = new TH1D(
        ("hAmpTrig_" + tag).c_str(),
        Form("Peak amplitude in trigger window [%.1f, %.1f] ns   "
             "V_{bias}=%.0fV, lux=%.0f"
             ";Peak amplitude  (mV);Events",
             t_trig_start, t_trig_end, info.vbias, info.lux),
        200, -5.0, MAX_THR);
    hAmpTrig->SetDirectory(nullptr);  // detach from rawFile — survives Close()
    hAmpTrig->SetLineColor(kRed+1);
    hAmpTrig->SetLineWidth(2);
    hAmpTrig->SetFillColorAlpha(kRed+1, 0.20);

    Long64_t nEntries   = treeCh1->GetEntries();
    long     nNoLaser   = 0;
    long     nInWindow  = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 1000 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  [laser canvas] Progress: "
                      << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->SetBranchAddress("time",      t1);
        treeCh1->SetBranchAddress("amplitude", a1);
        treeLaser->SetBranchAddress("time",      tL);
        treeLaser->SetBranchAddress("amplitude", aL);
        treeCh1->GetEntry(i);
        treeLaser->GetEntry(i);

        // Laser trigger time
        double t_laser = laserTriggerTime(tL, aL, N, 10.0);
        if (t_laser < -900.0) { ++nNoLaser; continue; }

        // SiPM: baseline correction + LP filter (full waveform)
        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, BASELINE_START, BASELINE_END),
            cutoff_MHz, fs_MHz);

        // Peak over FULL waveform → timing histogram
        auto   itMax  = std::max_element(af.begin(), af.end());
        double A_max  = *itMax;
        int    i_peak = (int)std::distance(af.begin(), itMax);
        double t_peak_rel = t1[i_peak] - t_laser;
        hTiming->Fill(t_peak_rel);

        // Peak inside trigger window → amplitude histogram
        double A_win = -1e9;
        for (int j = j_trig_start; j <= j_trig_end; ++j)
            if (af[j] > A_win) A_win = af[j];
        if (A_win > -1e8) {
            hAmpTrig->Fill(A_win);
            ++nInWindow;
        }
    }
    std::cout << "\r  [laser canvas] Progress: "
              << nEntries << " / " << nEntries << " — done.\n";
    std::cout << "  No laser trigger : " << nNoLaser
              << " / " << nEntries << " events\n";
    std::cout << "  In trigger window: " << nInWindow
              << " / " << nEntries << " events\n";

    // ── Canvas ───────────────────────────────────────────────
    std::string lcname = "cLaser_" + tag;
    TCanvas* cL = new TCanvas(lcname.c_str(),
                              ("Laser timing — " + tag).c_str(),
                              950, 1000);

    TPad* lpad1 = new TPad("lpad1", "timing", 0.0, 0.5, 1.0, 1.0);
    lpad1->SetLeftMargin(PAD_LEFT); lpad1->SetRightMargin(PAD_RIGHT);
    lpad1->SetTopMargin(PAD_TOP);   lpad1->SetBottomMargin(PAD_BOTTOM);
    lpad1->SetGrid(); lpad1->SetTicks(1,1);
    lpad1->Draw();

    TPad* lpad2 = new TPad("lpad2", "amplitude", 0.0, 0.0, 1.0, 0.5);
    lpad2->SetLeftMargin(PAD_LEFT); lpad2->SetRightMargin(PAD_RIGHT);
    lpad2->SetTopMargin(PAD_TOP);   lpad2->SetBottomMargin(PAD_BOTTOM);
    lpad2->SetGrid(); lpad2->SetTicks(1,1);
    lpad2->Draw();

    // Top pad: timing distribution
    lpad1->cd();
    hTiming->GetXaxis()->SetTitleSize(0.055);
    hTiming->GetXaxis()->SetLabelSize(0.050);
    hTiming->GetYaxis()->SetTitleSize(0.055);
    hTiming->GetYaxis()->SetLabelSize(0.050);
    hTiming->GetYaxis()->SetTitleOffset(1.1);
    hTiming->Draw("HIST");

    // Draw the trigger window boundaries as vertical dashed lines
    double t_ymax = hTiming->GetMaximum();
    TLine* lTrigLo = new TLine(t_trig_start, 0, t_trig_start, t_ymax);
    TLine* lTrigHi = new TLine(t_trig_end,   0, t_trig_end,   t_ymax);
    for (TLine* ll : {lTrigLo, lTrigHi}) {
        ll->SetLineColor(kOrange+7);
        ll->SetLineStyle(2);
        ll->SetLineWidth(2);
        ll->Draw("same");
    }

    // Label the trigger window
    TPaveText* ptTrig = new TPaveText(0.14, 0.78, 0.55, 0.88, "NDC");
    ptTrig->SetBorderSize(1); ptTrig->SetFillColor(0); ptTrig->SetFillStyle(1001);
    ptTrig->SetTextFont(42);  ptTrig->SetTextSize(0.038); ptTrig->SetTextAlign(12);
    ptTrig->AddText(Form("Trigger window: [%.1f, %.1f] ns  (orange dashes)",
                         t_trig_start, t_trig_end));
    ptTrig->AddText(Form("No laser trigger: %ld / %lld events", nNoLaser, nEntries));
    ptTrig->Draw();

    // Bottom pad: amplitude spectrum inside trigger window
    lpad2->cd();
    hAmpTrig->GetXaxis()->SetTitleSize(0.055);
    hAmpTrig->GetXaxis()->SetLabelSize(0.050);
    hAmpTrig->GetYaxis()->SetTitleSize(0.055);
    hAmpTrig->GetYaxis()->SetLabelSize(0.050);
    hAmpTrig->GetYaxis()->SetTitleOffset(1.1);
    hAmpTrig->Draw("HIST");

    // If the fit succeeded (fitMean passed as 0 means not available),
    // draw a note with the in-window event count
    TPaveText* ptAmp = new TPaveText(0.55, 0.78, 0.94, 0.88, "NDC");
    ptAmp->SetBorderSize(1); ptAmp->SetFillColor(0); ptAmp->SetFillStyle(1001);
    ptAmp->SetTextFont(42);  ptAmp->SetTextSize(0.038); ptAmp->SetTextAlign(12);
    ptAmp->AddText(Form("Events in window: %ld / %lld  (%.1f%%)",
                        nInWindow, nEntries,
                        100.0 * nInWindow / (double)nEntries));
    ptAmp->Draw();

    cL->Update(); cL->Modified();
    cL->SaveAs(("laser_timing_" + tag + ".png").c_str());
    std::cout << "  Saved: laser_timing_" << tag << ".png\n";
}

// ════════════════════════════════════════════════════════════
//  OVERLAY CANVAS — single pad, threshold scan + derivative
//
//  The threshold scan (blue, left Y axis) is drawn first and
//  defines the Y scale shown to the user.
//  The derivative -dN/dV (red) is scaled to "sit on top":
//    deriv_scaled[k] = counts_max * (deriv[k] / deriv_max)
//  so its peaks reach the top of the scan curve.
//  A TGaxis on the right side shows the true derivative scale.
//  The Gaussian fit and reference lines are drawn on top.
// ════════════════════════════════════════════════════════════
static void drawScanOverlay(const FileInfo&            info,
                            const std::vector<double>& thresholds,
                            const std::vector<double>& counts,
                            const std::vector<double>& thr_der_plot,
                            const std::vector<double>& deriv_plot,
                            const std::vector<double>& thr_der,
                            const std::vector<double>& deriv,
                            TF1*                       fitGaus,
                            double                     fitMean,
                            double                     fitMeanErr,
                            double                     fitSigma,
                            double                     fitSigmaErr,
                            double                     fitAmp,
                            double                     fitChi2,
                            int                        fitNdf,
                            double                     t_trig_start,
                            double                     t_trig_end) {

    int n_pts      = (int)thresholds.size();
    int n_der_plot = (int)thr_der_plot.size();

    // ── X range ──────────────────────────────────────────────
    double xhi_data = MIN_THR;
    for (int k = 0; k < n_pts; ++k)
        if (counts[k] > 0 && thresholds[k] > xhi_data) xhi_data = thresholds[k];
    double xRange = xhi_data - MIN_THR;
    double xlo    = MIN_THR  - xRange * 0.02;
    double xhi    = xhi_data + xRange * 0.03;

    // ── Y scale: defined entirely by the threshold scan ──────
    double counts_max = *std::max_element(counts.begin(), counts.end());
    double ylo = -counts_max * 0.05;
    double yhi =  counts_max * 1.15;

    // ── Scale derivative to sit within [ylo, yhi] ────────────
    // Find max of derivative in physical region (avoid noise spike near 0)
    double deriv_max = 0.0;
    for (int k = 0; k < n_der_plot; ++k)
        if (thr_der_plot[k] > MIN_POSITIVE_THR && deriv_plot[k] > deriv_max)
            deriv_max = deriv_plot[k];
    if (deriv_max <= 0) deriv_max = 1.0;

    // Map deriv_max → 90% of counts_max so peaks sit just below the scan top
    double scale = counts_max * 0.90 / deriv_max;

    std::vector<double> deriv_scaled(n_der_plot);
    for (int k = 0; k < n_der_plot; ++k)
        deriv_scaled[k] = deriv_plot[k] * scale;

    // ── Graphs ───────────────────────────────────────────────
    TGraph* grScan = new TGraph(n_pts, thresholds.data(), counts.data());
    grScan->SetLineColor(kAzure+1);
    grScan->SetLineWidth(2);
    grScan->SetMarkerColor(kAzure+1);
    grScan->SetMarkerStyle(20);
    grScan->SetMarkerSize(0.3);

    TGraph* grDer = new TGraph(n_der_plot, thr_der_plot.data(), deriv_scaled.data());
    grDer->SetLineColor(kRed+1);
    grDer->SetLineWidth(2);
    grDer->SetMarkerColor(kRed+1);
    grDer->SetMarkerStyle(20);
    grDer->SetMarkerSize(0.3);

    // ── Canvas ───────────────────────────────────────────────
    std::string cname = "c_overlay_" + info.tag;
    TCanvas* c = new TCanvas(cname.c_str(),
        ("Overlay — " + info.tag).c_str(), 1000, 650);
    c->SetLeftMargin(PAD_LEFT);
    c->SetRightMargin(PAD_RIGHT);
    c->SetTopMargin(PAD_TOP);
    c->SetBottomMargin(PAD_BOTTOM);
    c->SetGrid();
    c->SetTicks(1, 1);

    // ── Draw scan (sets axes) ─────────────────────────────────
    grScan->SetTitle(Form(
        "Threshold scan   V_{bias} = %.0f V,  lux = %.0f"
        ";Threshold (mV);Counts",
        info.vbias, info.lux));
    grScan->Draw("APL");
    grScan->GetXaxis()->SetRangeUser(xlo, xhi);
    grScan->GetYaxis()->SetRangeUser(ylo, yhi);
    grScan->GetXaxis()->SetTitleSize(0.050);
    grScan->GetXaxis()->SetLabelSize(0.045);
    grScan->GetYaxis()->SetTitleSize(0.050);
    grScan->GetYaxis()->SetLabelSize(0.045);
    grScan->GetYaxis()->SetTitleOffset(1.15);

    // ── Draw derivative data on top ───────────────────────────
    grDer->Draw("PL same");

    // ── Gaussian fit rescaled ────────────────────────────────
    if (fitGaus && fitMean > 0) {
        TF1* fitOvl = new TF1(("gaus_ovl_" + info.tag).c_str(),
                               "gaus",
                               fitGaus->GetXmin(), fitGaus->GetXmax());
        fitOvl->SetParameters(fitGaus->GetParameter(0) * scale,
                               fitGaus->GetParameter(1),
                               fitGaus->GetParameter(2));
        fitOvl->SetLineColor(kGreen+2);
        fitOvl->SetLineWidth(2);
        fitOvl->Draw("same");

        // ── Reference lines: always clipped to [ylo, yhi] ───
        // ROOT clips TLine to the pad range automatically when
        // the line endpoints are set exactly to ylo and yhi.
        const double thr_val[3] = { 0.5*fitMean, fitMean, 1.5*fitMean };
        const double thr_err[3] = { 0.5*fitMeanErr, fitMeanErr, 1.5*fitMeanErr };
        const int    thr_col[3] = { kCyan+1, kGreen+2, kMagenta+1 };
        const char*  thr_name[3]= { "0.5 p.e.", "1 p.e.", "1.5 p.e." };

        for (int r = 0; r < 3; ++r) {
            if (thr_val[r] < xlo || thr_val[r] > xhi) continue;
            // Use NDC-y endpoints so the line always spans the full
            // visible Y range regardless of zoom
            TLine* rl = new TLine(thr_val[r], ylo, thr_val[r], yhi);
            rl->SetLineColor(thr_col[r]);
            rl->SetLineStyle(2);
            rl->SetLineWidth(2);
            rl->Draw("same");
        }

        // ── Legend: only μ, σ, and reference lines ──────────
        double lx1 = 1.0 - PAD_RIGHT - 0.36;
        double lx2 = 1.0 - PAD_RIGHT - 0.01;
        double ly2 = 1.0 - PAD_TOP   - 0.02;
        double ly1 = ly2 - 0.30;
        TLegend* leg = new TLegend(lx1, ly1, lx2, ly2);
        leg->SetBorderSize(1);
        leg->SetFillStyle(1001);
        leg->SetFillColor(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.038);

        leg->AddEntry((TObject*)nullptr,
            Form("#mu = %.3f #pm %.3f mV", fitMean, fitMeanErr), "");
        leg->AddEntry((TObject*)nullptr,
            Form("#sigma = %.3f #pm %.3f mV", fitSigma, fitSigmaErr), "");

        for (int r = 0; r < 3; ++r) {
            if (thr_val[r] < xlo || thr_val[r] > xhi) continue;
            TLine* ld = new TLine();
            ld->SetLineColor(thr_col[r]);
            ld->SetLineStyle(2);
            ld->SetLineWidth(2);
            leg->AddEntry(ld,
                Form("%s: %.3f #pm %.3f mV", thr_name[r], thr_val[r], thr_err[r]),
                "l");
        }
        leg->Draw();
    }

    c->Update(); c->Modified();
    c->SaveAs(("overlay_" + info.tag + ".png").c_str());
    std::cout << "  Saved: overlay_" << info.tag << ".png\n";
}

// ════════════════════════════════════════════════════════════
//  PROCESS A SINGLE FILE
//
//  Parameters:
//    info          — file path and metadata
//    cutoff_MHz    — low-pass filter cut-off frequency
//    t_base_start  — start of baseline estimation window (ns)  [NEW-1]
//    t_base_end    — end   of baseline estimation window (ns)  [NEW-1]
//    t_trig_start  — start of trigger search window (ns)       [NEW-1]
//    t_trig_end    — end   of trigger search window (ns)       [NEW-1]
// ════════════════════════════════════════════════════════════
static void processFile(const FileInfo& info,
                        double cutoff_MHz,
                        double t_trig_start,
                        double t_trig_end) {

    // Baseline window is fixed globally — no user input needed
    const double t_base_start = BASELINE_START;
    const double t_base_end   = BASELINE_END;

    TFile* file = TFile::Open(info.path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[SKIP] Cannot open: " << info.path << "\n"; return;
    }
    TTree* treeCh1 = (TTree*)file->Get("ch1");
    if (!treeCh1) {
        std::cerr << "[SKIP] Tree 'ch1' not found in: " << info.path << "\n";
        file->Close(); return;
    }

    const int N = 1024;
    Double_t t1[N], a1[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);
    treeCh1->GetEntry(0);
    double fs_MHz = 1000.0 / (t1[1] - t1[0]);

    // Locate trigger window sample indices.
    // A small epsilon (half a sample interval) is added to t_trig_end to
    // avoid missing the last sample due to floating-point rounding:
    // e.g. with fs=5000 MHz, t1[1000]=200.000...01 ns would be excluded
    // by a strict <= 200.0 comparison even though it is nominally inside
    // the requested window.
    const double dt_ns  = t1[1] - t1[0];          // sample interval (ns)
    const double eps    = 0.5 * dt_ns;             // half-sample tolerance

    int j_trig_start = 0, j_trig_end = N - 1;
    for (int j = 0; j < N; ++j)
        if (t1[j] >= t_trig_start - eps) { j_trig_start = j; break; }
    for (int j = N - 1; j >= 0; --j)
        if (t1[j] <= t_trig_end + eps)   { j_trig_end   = j; break; }
    if (j_trig_start >= j_trig_end) {
        std::cerr << "[WARN] Trigger window [" << t_trig_start << ", "
                  << t_trig_end << "] ns is invalid for " << info.tag
                  << " (waveform spans [" << t1[0] << ", " << t1[N-1]
                  << "] ns). Using full waveform.\n";
        j_trig_start = 0; j_trig_end = N - 1;
    }

    std::cout << "\n[" << info.tag << "]  "
              << treeCh1->GetEntries() << " events,  fs=" << fs_MHz << " MHz\n"
              << "  Trigger  window : [" << t1[j_trig_start] << ", "
              << t1[j_trig_end] << "] ns  (samples "
              << j_trig_start << "-" << j_trig_end << ")\n";

    int n_pts = (int)std::round((MAX_THR - MIN_THR) / STEP_SCAN) + 1;
    std::vector<double> thresholds(n_pts), counts(n_pts, 0.0);
    for (int k = 0; k < n_pts; ++k) thresholds[k] = MIN_THR + k * STEP_SCAN;

    Long64_t nEntries = treeCh1->GetEntries();

    // ── Event loop ──────────────────────────────────────────
    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  Progress: " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);

        // [NEW-1] Baseline estimated from baseline window, independent of trigger
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, t_base_start, t_base_end),
            cutoff_MHz, fs_MHz);

        for (int k = 0; k < n_pts; ++k) {
            double thr = thresholds[k];
            bool crossed = false;

            if (thr > 0) {
                // Hysteresis trigger: armed when signal drops below thr/2,
                // fires when signal rises above thr.
                // Search restricted to trigger window [NEW-1].
                bool armed = false;
                for (int j = j_trig_start; j <= j_trig_end; ++j) {
                    double y = af[j];
                    if (!armed) { if (y < thr * 0.5) armed = true; }
                    else        { if (y > thr) { crossed = true; break; } }
                }
            } else {
                // For thr <= 0: count if the minimum in the trigger window
                // falls below thr (measures noise floor / undershoot).
                double ymin = 0.0;
                for (int j = j_trig_start; j <= j_trig_end; ++j)
                    if (af[j] < ymin) ymin = af[j];
                if (ymin < thr) crossed = true;
            }
            if (crossed) counts[k] += 1.0;
        }
    }
    file->Close();
    std::cout << "\r  Progress: " << nEntries << " / " << nEntries << " — done.\n";

    // ── Poisson errors ──────────────────────────────────────
    // Empty bins have error = 0 (no measured uncertainty)
    std::vector<double> errors(n_pts);
    for (int k = 0; k < n_pts; ++k)
        errors[k] = (counts[k] > 0) ? std::sqrt(counts[k]) : 0.0;

    // ── Centred numerical derivative (fine, for fit and peak finding) ──
    int h = (int)std::round(STEP_DER / STEP_SCAN);
    if (h < 1) h = 1;
    int n_der = n_pts - 2*h;
    std::vector<double> thr_der(n_der), deriv(n_der), err_der(n_der),
                        xerr_der(n_der, 0.0);
    for (int k = h; k < n_pts - h; ++k) {
        int idx      = k - h;
        thr_der[idx] = thresholds[k];
        deriv[idx]   = -(counts[k+h] - counts[k-h]) / (2.0 * STEP_DER);
        double eU    = (counts[k+h] > 0) ? std::sqrt(counts[k+h]) : 0.0;
        double eD    = (counts[k-h] > 0) ? std::sqrt(counts[k-h]) : 0.0;
        err_der[idx] = std::sqrt(eU*eU + eD*eD) / (2.0 * STEP_DER);
    }

    // ── Centred numerical derivative (coarse, for display only) ─
    // Computed at STEP_DER_PLOT = 0.5 mV — smoother curve on the canvas.
    // The fit and peak finder always use the fine derivative above.
    int h_plot = (int)std::round(STEP_DER_PLOT / STEP_SCAN);
    if (h_plot < 1) h_plot = 1;
    int n_der_plot = n_pts - 2*h_plot;
    std::vector<double> thr_der_plot(n_der_plot), deriv_plot(n_der_plot),
                        err_der_plot(n_der_plot), xerr_plot(n_der_plot, 0.0);
    for (int k = h_plot; k < n_pts - h_plot; ++k) {
        int idx           = k - h_plot;
        thr_der_plot[idx] = thresholds[k];
        deriv_plot[idx]   = -(counts[k+h_plot] - counts[k-h_plot]) / (2.0 * STEP_DER_PLOT);
        double eU         = (counts[k+h_plot] > 0) ? std::sqrt(counts[k+h_plot]) : 0.0;
        double eD         = (counts[k-h_plot] > 0) ? std::sqrt(counts[k-h_plot]) : 0.0;
        err_der_plot[idx] = std::sqrt(eU*eU + eD*eD) / (2.0 * STEP_DER_PLOT);
    }

    // ── TGraphErrors ────────────────────────────────────────
    std::vector<double> x_err_zero(n_pts, 0.0);
    TGraphErrors* grScan = new TGraphErrors(n_pts,
        thresholds.data(), counts.data(), x_err_zero.data(), errors.data());
    grScan->SetTitle(Form(
        "Threshold scan   V_{bias} = %.0f V,  lux = %.0f;Threshold (mV);Counts",
        info.vbias, info.lux));
    grScan->SetLineColor(kAzure+1); grScan->SetMarkerColor(kAzure+1);
    grScan->SetMarkerStyle(20); grScan->SetMarkerSize(0.4); grScan->SetLineWidth(2);

    // grDeriv uses the COARSE derivative (STEP_DER_PLOT) for a cleaner display.
    // The fit and peak finder work on the fine derivative (STEP_DER) above.
    TGraphErrors* grDeriv = new TGraphErrors(n_der_plot,
        thr_der_plot.data(), deriv_plot.data(), xerr_plot.data(), err_der_plot.data());
    grDeriv->SetTitle(Form(
        "-dN/dV   V_{bias} = %.0f V,  lux = %.0f;Threshold (mV);-dN/dV  (mV^{-1})",
        info.vbias, info.lux));
    grDeriv->SetLineColor(kRed+1); grDeriv->SetMarkerColor(kRed+1);
    grDeriv->SetMarkerStyle(20); grDeriv->SetMarkerSize(0.4); grDeriv->SetLineWidth(2);

    // ── Gaussian fit on 1st p.e. peak  [NEW-2] ──────────────
    int iPeak = findFirstPeak(thr_der, deriv);
    TF1* fitGaus = nullptr;
    double fitMean=0, fitMeanErr=0, fitSigma=0, fitSigmaErr=0;
    double fitAmp=0, fitAmpErr=0, fitChi2=0;
    int    fitNdf=0;

    if (iPeak >= 0) {
        double center = thr_der[iPeak];
        double amp0   = deriv[iPeak];

        // Estimate sigma from left half-maximum point
        double sigmaEst = 3.0;
        for (int i = iPeak; i > 0; --i) {
            if (thr_der[i] <= MIN_POSITIVE_THR) break;
            if (deriv[i] < amp0 * 0.5) {
                sigmaEst = (center - thr_der[i]) / 1.177; break;
            }
        }
        sigmaEst = std::max(FIT_SIGMA_MIN, std::min(sigmaEst, FIT_SIGMA_MAX));

        // Fit range: centre ± 2*sigmaEst, never below MIN_POSITIVE_THR.
        // Using 2σ (not 3σ) avoids including the valley between adjacent
        // p.e. peaks, which would pull the Gaussian tail below zero.
        double fitLo = std::max(MIN_POSITIVE_THR, center - 2.0 * sigmaEst);
        double fitHi = center + 2.0 * sigmaEst;

        fitGaus = new TF1(("gaus_" + info.tag).c_str(), "gaus", fitLo, fitHi);
        fitGaus->SetParameters(amp0, center, sigmaEst);
        fitGaus->SetParNames("Amplitude", "Mean", "Sigma");
        fitGaus->SetParLimits(1, MIN_POSITIVE_THR, MAX_THR);
        fitGaus->SetParLimits(2, FIT_SIGMA_MIN, FIT_SIGMA_MAX);
        fitGaus->SetLineColor(kGreen+2); fitGaus->SetLineWidth(2);
        grDeriv->Fit(fitGaus, "RQ");

        fitMean     = fitGaus->GetParameter(1);
        fitMeanErr  = fitGaus->GetParError(1);
        fitSigma    = std::abs(fitGaus->GetParameter(2));
        fitSigmaErr = fitGaus->GetParError(2);
        fitAmp      = fitGaus->GetParameter(0);
        fitAmpErr   = fitGaus->GetParError(0);
        fitChi2     = fitGaus->GetChisquare();
        fitNdf      = fitGaus->GetNDF();

        std::cout << "  Fit 1st peak: mean=" << fitMean << " +/- " << fitMeanErr
                  << "  sigma=" << fitSigma << " +/- " << fitSigmaErr
                  << "  chi2/ndf=" << fitChi2 << "/" << fitNdf;
        if (fitNdf > 0) {
            double r = fitChi2 / fitNdf;
            if (r > 5.0) std::cout << "  *** WARNING: chi2/ndf=" << r << " > 5 (poor fit) ***";
            if (r < 0.3) std::cout << "  *** WARNING: chi2/ndf=" << r << " < 0.3 (over-free?) ***";
        }
        std::cout << "\n";
    } else {
        std::cerr << "  [WARN] First peak not found for " << info.tag << "\n";
    }

    // ════════════════════════════════════════════════════════
    //  CANVAS — two explicit TPads with identical margins
    //  pad1 (scan)  : top half    (y: 0.5 → 1.0)
    //  pad2 (deriv) : bottom half (y: 0.0 → 0.5)
    //  Identical left/right/top/bottom margins → same physical
    //  data-area width → visually consistent X scales.
    // ════════════════════════════════════════════════════════
    std::string cname = "c_" + info.tag;
    TCanvas* c = new TCanvas(cname.c_str(),
                             ("Threshold scan — " + info.tag).c_str(),
                             950, 1000);

    TPad* pad1 = new TPad("pad1", "scan",  0.0, 0.5, 1.0, 1.0);
    pad1->SetLeftMargin(PAD_LEFT);   pad1->SetRightMargin(PAD_RIGHT);
    pad1->SetTopMargin(PAD_TOP);     pad1->SetBottomMargin(PAD_BOTTOM);
    pad1->SetGrid(); pad1->SetTicks(1,1);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "deriv", 0.0, 0.0, 1.0, 0.5);
    pad2->SetLeftMargin(PAD_LEFT);   pad2->SetRightMargin(PAD_RIGHT);
    pad2->SetTopMargin(PAD_TOP);     pad2->SetBottomMargin(PAD_BOTTOM);
    pad2->SetGrid(); pad2->SetTicks(1,1);
    pad2->Draw();

    // Common X range: from MIN_THR to last bin with counts > 0
    double xhi_data = MIN_THR;
    for (int k = 0; k < n_pts; ++k)
        if (counts[k] > 0 && thresholds[k] > xhi_data) xhi_data = thresholds[k];
    double xRange = xhi_data - MIN_THR;
    double xlo    = MIN_THR  - xRange * 0.02;
    double xhi    = xhi_data + xRange * 0.03;

    // ── Pad 1: threshold scan ────────────────────────────────
    pad1->cd();
    grScan->Draw("APL");
    grScan->GetXaxis()->SetRangeUser(xlo, xhi);
    grScan->GetXaxis()->SetTitleSize(0.055); grScan->GetXaxis()->SetLabelSize(0.050);
    grScan->GetYaxis()->SetTitleSize(0.055); grScan->GetYaxis()->SetLabelSize(0.050);
    grScan->GetYaxis()->SetTitleOffset(1.1);

    // Vertical reference lines on the scan pad + legend with values.
    // Reference positions derived from fitMean (1st p.e. threshold):
    //   0.5 pe : 0.5 × μ  →  σ = 0.5 × σ_μ  (linear error propagation)
    //   1   pe : μ         →  σ = σ_μ
    //   1.5 pe : 1.5 × μ  →  σ = 1.5 × σ_μ  (linear error propagation)
    if (fitGaus && fitMean > 0) {
        double scan_ymax = *std::max_element(counts.begin(), counts.end());
        double scan_ymin = 0.0;

        // Three reference thresholds and their propagated errors
        const double thr_val[3]  = { 0.5*fitMean,       fitMean,      1.5*fitMean      };
        const double thr_err[3]  = { 0.5*fitMeanErr,    fitMeanErr,   1.5*fitMeanErr   };
        const char*  thr_name[3] = { "0.5 p.e.",        "1 p.e.",     "1.5 p.e."       };
        // Distinct colours: cyan, green, magenta — all visible on blue scan curve
        const int    thr_col[3]  = { kCyan+1,           kGreen+2,     kMagenta+1       };
        const int    ref_style   = 2;    // standard dashes — clearly visible
        const float  ref_width   = 2.5f; // thick enough to stand out

        TLine* refLines[3];
        for (int r = 0; r < 3; ++r) {
            if (thr_val[r] < xlo || thr_val[r] > xhi) continue;
            refLines[r] = new TLine(thr_val[r], scan_ymin, thr_val[r], scan_ymax);
            refLines[r]->SetLineColor(thr_col[r]);
            refLines[r]->SetLineStyle(ref_style);
            refLines[r]->SetLineWidth(ref_width);
            refLines[r]->Draw("same");
        }

        // Legend in the top-right corner of the scan pad
        // Height sized to fit 4 entries (μ header + 3 thresholds)
        double lx1 = 1.0 - PAD_RIGHT - 0.42;
        double lx2 = 1.0 - PAD_RIGHT - 0.01;
        double ly2 = 1.0 - PAD_TOP   - 0.02;
        double ly1 = ly2 - 0.28;
        TLegend* scanLeg = new TLegend(lx1, ly1, lx2, ly2);
        scanLeg->SetBorderSize(1);
        scanLeg->SetFillStyle(1001);
        scanLeg->SetFillColor(0);
        scanLeg->SetTextFont(42);
        scanLeg->SetTextSize(0.040);

        // Header: fit mean with its error
        scanLeg->AddEntry((TObject*)nullptr,
            Form("#mu = %.3f #pm %.3f mV", fitMean, fitMeanErr), "");

        // One entry per reference line: name, value ± error
        for (int r = 0; r < 3; ++r) {
            // Create a dummy TLine just to get the coloured line in the legend
            TLine* leg_dummy = new TLine();
            leg_dummy->SetLineColor(thr_col[r]);
            leg_dummy->SetLineStyle(ref_style);
            leg_dummy->SetLineWidth(ref_width);
            scanLeg->AddEntry(leg_dummy,
                Form("%s: %.3f #pm %.3f mV", thr_name[r], thr_val[r], thr_err[r]),
                "l");
        }
        scanLeg->Draw();
    }

    // ── Pad 2: derivative ────────────────────────────────────
    pad2->cd();
    grDeriv->Draw("APL");
    grDeriv->GetXaxis()->SetRangeUser(xlo, xhi);   // same X range as pad1

    // Y zoom: exclude |thr| <= MIN_POSITIVE_THR (spike region)
    // Uses the COARSE derivative (same as grDeriv) for consistent scaling
    double ymax_pos=0, ymin_neg=0, ymin_display=0;
    for (int k = 0; k < n_der_plot; ++k) {
        if (std::abs(thr_der_plot[k]) <= MIN_POSITIVE_THR) continue;
        if (deriv_plot[k] > ymax_pos) ymax_pos = deriv_plot[k];
        if (deriv_plot[k] < ymin_neg) ymin_neg = deriv_plot[k];
    }
    if (ymax_pos > 0) {
        ymin_display = std::min(ymin_neg * 1.3, -ymax_pos * 0.15);
        grDeriv->GetYaxis()->SetRangeUser(ymin_display, ymax_pos * 1.3);
    }
    grDeriv->GetXaxis()->SetTitleSize(0.055); grDeriv->GetXaxis()->SetLabelSize(0.050);
    grDeriv->GetYaxis()->SetTitleSize(0.055); grDeriv->GetYaxis()->SetLabelSize(0.050);
    grDeriv->GetYaxis()->SetTitleOffset(1.1);

    if (fitGaus) {
        fitGaus->Draw("same");

        // Vertical dashed line at mu, from ymin_display to peak amplitude
        TLine* vline = new TLine(fitMean, ymin_display, fitMean, fitAmp);
        vline->SetLineColor(kGreen+2); vline->SetLineStyle(2); vline->SetLineWidth(2);
        vline->Draw("same");

        // Same three reference lines as the scan pad (no legend — colours speak)
        // 0.5 p.e. = cyan, 1 p.e. = green (coincides with mu line above),
        // 1.5 p.e. = magenta
        if (fitMean > 0) {
            const double thr_val[3] = { 0.5*fitMean, fitMean, 1.5*fitMean };
            const int    thr_col[3] = { kCyan+1, kGreen+2, kMagenta+1 };
            for (int r = 0; r < 3; ++r) {
                if (thr_val[r] < xlo || thr_val[r] > xhi) continue;
                TLine* rl = new TLine(thr_val[r], ymin_display,
                                      thr_val[r], ymax_pos * 1.3);
                rl->SetLineColor(thr_col[r]);
                rl->SetLineStyle(2);
                rl->SetLineWidth(2.5);
                rl->Draw("same");
            }
        }

        // Result box: position on the less-crowded half
        double relPos = (fitMean - xlo) / (xhi - xlo);
        double bx1 = (relPos > 0.5) ? 0.14 : 0.55;
        double bx2 = bx1 + 0.40;
        double by2 = 0.88, by1 = by2 - 0.30;

        TPaveText* pave = new TPaveText(bx1, by1, bx2, by2, "NDC");
        pave->SetBorderSize(1); pave->SetFillColor(0); pave->SetFillStyle(1001);
        pave->SetTextFont(42);  pave->SetTextSize(0.038); pave->SetTextAlign(12);
        pave->AddText(Form("Trigger: [%.1f, %.1f] ns", t_trig_start, t_trig_end));
        pave->AddText(Form("#mu = %.3f #pm %.3f mV", fitMean, fitMeanErr));
        pave->AddText(Form("#sigma = %.3f #pm %.3f mV", fitSigma, fitSigmaErr));
        pave->AddText(Form("#chi^{2}/ndf = %.1f / %d", fitChi2, fitNdf));
        pave->Draw("same");
    }

    c->Update(); c->Modified();
    c->SaveAs(("scan_" + info.tag + ".png").c_str());

    // ── Overlay canvas: scan + derivative on a single pad ────
    drawScanOverlay(info, thresholds, counts,
                    thr_der_plot, deriv_plot,
                    thr_der, deriv,
                    fitGaus,
                    fitMean, fitMeanErr,
                    fitSigma, fitSigmaErr,
                    fitAmp,
                    fitChi2, fitNdf,
                    t_trig_start, t_trig_end);

    // ── Save to ROOT file ────────────────────────────────────
    std::string outname = "scan_" + info.tag + ".root";
    TFile* fout = new TFile(outname.c_str(), "RECREATE");

    double b_vbias=info.vbias, b_lux=info.lux, b_cutoff=cutoff_MHz;
    double b_tbase_start=t_base_start, b_tbase_end=t_base_end;
    double b_ttrig_start=t_trig_start, b_ttrig_end=t_trig_end;
    double b_thr, b_count, b_err_count, b_deriv_val, b_err_deriv;

    TTree* tScan = new TTree("scan", "Threshold scan");
    tScan->Branch("vbias",        &b_vbias,       "vbias/D");
    tScan->Branch("lux",          &b_lux,         "lux/D");
    tScan->Branch("cutoff_MHz",   &b_cutoff,      "cutoff_MHz/D");
    tScan->Branch("t_base_start", &b_tbase_start, "t_base_start/D");
    tScan->Branch("t_base_end",   &b_tbase_end,   "t_base_end/D");
    tScan->Branch("t_trig_start", &b_ttrig_start, "t_trig_start/D");
    tScan->Branch("t_trig_end",   &b_ttrig_end,   "t_trig_end/D");
    tScan->Branch("threshold",    &b_thr,         "threshold/D");
    tScan->Branch("counts",       &b_count,       "counts/D");
    tScan->Branch("err_counts",   &b_err_count,   "err_counts/D");
    for (int k = 0; k < n_pts; ++k) {
        b_thr=thresholds[k]; b_count=counts[k]; b_err_count=errors[k];
        tScan->Fill();
    }

    TTree* tDeriv = new TTree("deriv", "Derivative -dN/dV (fine step for analysis)");
    double b_step_der = STEP_DER;
    tDeriv->Branch("vbias",        &b_vbias,       "vbias/D");
    tDeriv->Branch("lux",          &b_lux,         "lux/D");
    tDeriv->Branch("cutoff_MHz",   &b_cutoff,      "cutoff_MHz/D");
    tDeriv->Branch("step_der_mV",  &b_step_der,    "step_der_mV/D");
    tDeriv->Branch("t_base_start", &b_tbase_start, "t_base_start/D");
    tDeriv->Branch("t_base_end",   &b_tbase_end,   "t_base_end/D");
    tDeriv->Branch("t_trig_start", &b_ttrig_start, "t_trig_start/D");
    tDeriv->Branch("t_trig_end",   &b_ttrig_end,   "t_trig_end/D");
    tDeriv->Branch("threshold",    &b_thr,         "threshold/D");
    tDeriv->Branch("deriv",        &b_deriv_val,   "deriv/D");
    tDeriv->Branch("err_deriv",    &b_err_deriv,   "err_deriv/D");
    for (int k = 0; k < n_der; ++k) {
        b_thr=thr_der[k]; b_deriv_val=deriv[k]; b_err_deriv=err_der[k];
        tDeriv->Fill();
    }

    TTree* tFit = new TTree("fit_results", "Gaussian fit on 1st p.e. peak");
    double b_mean, b_mean_err, b_sigma, b_sigma_err, b_amp, b_amp_err, b_chi2;
    int    b_ndf;
    tFit->Branch("vbias",        &b_vbias,       "vbias/D");
    tFit->Branch("lux",          &b_lux,         "lux/D");
    tFit->Branch("t_base_start", &b_tbase_start, "t_base_start/D");
    tFit->Branch("t_base_end",   &b_tbase_end,   "t_base_end/D");
    tFit->Branch("t_trig_start", &b_ttrig_start, "t_trig_start/D");
    tFit->Branch("t_trig_end",   &b_ttrig_end,   "t_trig_end/D");
    tFit->Branch("mean",         &b_mean,        "mean/D");
    tFit->Branch("mean_err",     &b_mean_err,    "mean_err/D");
    tFit->Branch("sigma",        &b_sigma,       "sigma/D");
    tFit->Branch("sigma_err",    &b_sigma_err,   "sigma_err/D");
    tFit->Branch("amplitude",    &b_amp,         "amplitude/D");
    tFit->Branch("amp_err",      &b_amp_err,     "amp_err/D");
    tFit->Branch("chi2",         &b_chi2,        "chi2/D");
    tFit->Branch("ndf",          &b_ndf,         "ndf/I");
    b_mean=fitMean; b_mean_err=fitMeanErr; b_sigma=fitSigma; b_sigma_err=fitSigmaErr;
    b_amp=fitAmp;   b_amp_err=fitAmpErr;   b_chi2=fitChi2;   b_ndf=fitNdf;
    tFit->Fill();

    fout->Write(); fout->Close();
    std::cout << "  Saved: " << outname << "\n";

    // ── Laser timing canvas (requires 'laser' tree in raw file) ──
    // Reopen the raw input file — it was closed after the event loop.
    TFile* rawForLaser = TFile::Open(info.path.c_str(), "READ");
    if (rawForLaser && !rawForLaser->IsZombie()) {
        makeLaserCanvas(rawForLaser, info, cutoff_MHz,
                        t_trig_start, t_trig_end, info.tag);
        rawForLaser->Close();
    }
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════
void sipm_threshold_scan_full() {
    gStyle->SetOptStat(0);

    // ── Scan current directory for input files ───────────────
    void* dirHandle = gSystem->OpenDirectory("../../data");
    if (!dirHandle) { std::cerr << "Cannot open current directory.\n"; return; }

    std::vector<FileInfo> files;
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirHandle)) != nullptr) {
        std::string fname(entry);
        if (fname.find("data.vbias_{") == std::string::npos) continue;
        if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
        FileInfo info;
        if (parseFilename("../../data/" + fname, info)) {
            files.push_back(info);
            std::cout << "Found: " << fname
                      << "  →  vbias=" << info.vbias
                      << " V,  lux=" << info.lux << "\n";
        } else {
            std::cerr << "Cannot parse filename: " << fname << "\n";
        }
    }
    gSystem->FreeDirectory(dirHandle);

    if (files.empty()) {
        std::cerr << "No data.vbias_{V}_L.root files found.\n"; return;
    }
    std::sort(files.begin(), files.end(), [](const FileInfo& a, const FileInfo& b){
        return (a.vbias != b.vbias) ? a.vbias < b.vbias : a.lux < b.lux;
    });

    // ── Interactive input ────────────────────────────────────
    double cutoff_MHz;
    std::cout << "\nLow-pass filter cut-off [MHz]: ";
    std::cin >> cutoff_MHz;

    // Read the time axis of the first file to inform the operator
    // of the actual waveform range before asking for the trigger window.
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
                t_wf_max = t0[N0 - 1];
            }
            f0->Close();
        }
    }

    // Trigger window — the only thing the operator needs to specify.
    // Baseline is always computed over [BASELINE_START, BASELINE_END) ns.
    double t_trig_start, t_trig_end;
    std::cout << "\n--- Trigger search window ---\n";
    std::cout << "  Waveform range  : [" << t_wf_min << ", " << t_wf_max << "] ns\n";
    std::cout << "  Baseline (fixed): [" << BASELINE_START << ", " << BASELINE_END << ") ns\n";
    std::cout << "  Start [ns]: ";
    std::cin >> t_trig_start;
    std::cout << "  End   [ns]: ";
    std::cin >> t_trig_end;
    if (t_trig_start >= t_trig_end) {
        std::cerr << "[ERROR] t_trig_start >= t_trig_end. Aborting.\n"; return;
    }

    std::cout << "\nSettings summary:\n"
              << "  LP cut-off      : " << cutoff_MHz    << " MHz\n"
              << "  Trigger window  : [" << t_trig_start << ", "
                                         << t_trig_end   << "] ns\n\n";

    for (auto& f : files)
        processFile(f, cutoff_MHz, t_trig_start, t_trig_end);

    std::cout << "\n=== All files processed. ===\n";
}
