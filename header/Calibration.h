#pragma once
// Calibration.h
// p.e. calibration via threshold scan — TOT-compatible method.
//
// Physics: vary the discriminator threshold and count events that exceed it.
// N(threshold) is the cumulative amplitude distribution.
// -dN/dV shows discrete peaks at multiples of the p.e. charge.
//
// Exported:
//   calibrateSpectrum()   — single-file calibration, saves canvas
//   drawCalibOverlay()    — overlay of multiple files
//   printFileList()       — interactive helpers
//   selectOneFile()
//   selectMultipleFiles()

#include "Config.h"
#include "SignalProcessing.h"
#include "OutputManager.h"
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TSystem.h>

// ── Interactive file helpers ─────────────────────────────────
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
            std::cout << "  [ERROR] Enter 0-" << files.size()-1 << ".\n";
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
            std::cout << "  [ERROR] Enter e.g. \"2\" or \"0,2\" or \"all\".\n";
        }
    }
    return selected;
}

// ── Result structs ───────────────────────────────────────────

// Calibration result — passed to all downstream analysis functions.
struct CalibResult {
    double m            = 0;      // gain   [mV/p.e.]
    double q            = 0;      // offset [mV]  (threshold for n p.e. = q + n*m)
    bool   ok           = false;
    double laser_thr    = 10.0;   // laser trigger threshold [mV] — set by FilterDiagnostics
    double t_trig_start = 95.0;   // trigger window start [ns]
    double t_trig_end   = 125.0;  // trigger window end   [ns]
    double cutoff_MHz   = 200.0;  // LP filter cutoff used during calibration
};

// Raw scan data — filled by calibrateSpectrum() when scanOut != nullptr.
// Needed only by drawCalibOverlay(); not required for normal analysis.
struct CalibScanData {
    std::string         tag;
    std::vector<double> thresholds;  // scan thresholds [mV]
    std::vector<double> counts;      // N(threshold): events above each threshold
    std::vector<double> thr_der;     // derivative x-axis [mV]
    std::vector<double> deriv;       // -dN/dV [mV⁻¹]
    std::vector<double> pkPos;       // accepted peak positions [mV]
    double              gain   = 0;
    double              offset = 0;
};

static int calibColor(int idx) {
    static const int pal[] = {
        kAzure+1, kOrange+7, kGreen+2, kMagenta+1,
        kRed+1,   kCyan+2,   kViolet+1, kYellow+3
    };
    return pal[idx % (int)(sizeof(pal)/sizeof(pal[0]))];
}

// ════════════════════════════════════════════════════════════
//  calibrateSpectrum()
//  Performs a threshold scan on treeCh1, fits the p.e. spectrum,
//  and returns gain [mV/p.e.] and offset [mV].
//
//  Algorithm:
//    1. Estimate baseline noise RMS on first 200 events.
//       ARM_LEVEL = 2*noise_rms is the fixed arming threshold for the
//       discriminator (independent of the scanned threshold).
//    2. Event loop: for each threshold, arm below ARM_LEVEL, fire above threshold.
//    3. Compute centred derivative -dN/dV (= amplitude spectrum).
//    4. Find peaks: 9-point smooth, local maxima, deduplication.
//    5. Estimate gain_raw from the gap between 1 p.e. and 2 p.e. peaks.
//    6. Assign n_pe = round((pos - offset_guess)/gain_raw) to each peak.
//       Peaks that land >0.4 p.e. from an integer are discarded as spurious
//       (optical crosstalk, afterpulse, noise).
//    7. Linear fit: peak_position vs n_pe → final gain and offset.
//
//  If scanOut != nullptr, it is filled for use by drawCalibOverlay().
// ════════════════════════════════════════════════════════════
static CalibResult calibrateSpectrum(TTree*             treeCh1,
                                      double             cutoff_MHz,
                                      double             fs_MHz,
                                      const std::string& tag,
                                      OutCtx&            ctx,
                                      CalibScanData*     scanOut = nullptr) {
    CalibResult res{0, 0, false};
    const int N = 1024;
    Double_t t1[N], a1[N];
    treeCh1->SetBranchAddress("time",      t1);
    treeCh1->SetBranchAddress("amplitude", a1);
    treeCh1->GetEntry(0);

    const double SCAN_MIN  = -10.0;  // mV — negative end shows pedestal in plot
    const double SCAN_MAX  = 110.0;  // mV
    const double SCAN_STEP =   0.2;  // mV
    const double DER_STEP  =   1.0;  // mV — centred-difference step

    int n_pts = (int)std::round((SCAN_MAX - SCAN_MIN) / SCAN_STEP) + 1;
    std::vector<double> thresholds(n_pts), counts(n_pts, 0.0);
    for (int k = 0; k < n_pts; ++k)
        thresholds[k] = SCAN_MIN + k * SCAN_STEP;

    Long64_t nEntries = treeCh1->GetEntries();
    std::cout << "  Calibrating " << tag << " (" << nEntries << " events)...\n";

    // ── Step 1: noise RMS → fixed arming level ───────────────
    // ARM_LEVEL is fixed per-file, not proportional to threshold.
    // This avoids the pathological case where thr/2 falls inside the noise
    // band (for low thresholds) or below zero (for negative thresholds).
    const double ARM_SIGMA = 2.0;
    double noise_rms = 1.0;
    {
        double sum2 = 0.0; int cnt = 0;
        Long64_t nNoise = std::min(nEntries, (Long64_t)200);
        for (Long64_t i = 0; i < nNoise; ++i) {
            treeCh1->GetEntry(i);
            std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
            std::vector<double> bc = correctBaseline(v_t, v_a,
                                                      BASELINE_START, BASELINE_END);
            for (int j = 0; j < N; ++j)
                if (v_t[j] >= BASELINE_START && v_t[j] < BASELINE_END) {
                    sum2 += bc[j] * bc[j]; ++cnt;
                }
        }
        if (cnt > 0) noise_rms = std::sqrt(sum2 / cnt);
        noise_rms = std::max(noise_rms, 0.2);
    }
    const double ARM_LEVEL = ARM_SIGMA * noise_rms;
    std::cout << "  [CAL] Noise RMS = " << std::fixed << std::setprecision(2)
              << noise_rms << " mV   arm level = " << ARM_LEVEL << " mV\n";

    // ── Step 2: threshold scan ───────────────────────────────
    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = correctBaseline(v_t, v_a, BASELINE_START, BASELINE_END);
        if (cutoff_MHz > 0)
            af = butterworthLowPass(af, cutoff_MHz, fs_MHz);

        for (int k = 0; k < n_pts; ++k) {
            double thr = thresholds[k];
            bool armed = false, crossed = false;
            for (int j = 0; j < N; ++j) {
                double y = af[j];
                if (!armed) { if (y < ARM_LEVEL) armed = true; }
                else        { if (y > thr) { crossed = true; break; } }
            }
            if (crossed) counts[k] += 1.0;
        }
    }

    // ── Step 3: centred derivative -dN/dV ───────────────────
    int h = (int)std::round(DER_STEP / SCAN_STEP);
    if (h < 1) h = 1;
    int n_der = n_pts - 2*h;
    std::vector<double> thr_der(n_der), deriv(n_der), err_der(n_der);
    for (int k = h; k < n_pts - h; ++k) {
        int idx      = k - h;
        thr_der[idx] = thresholds[k];
        deriv[idx]   = -(counts[k+h] - counts[k-h]) / (2.0 * DER_STEP);
        double eU    = (counts[k+h] > 0) ? std::sqrt(counts[k+h]) : 0.0;
        double eD    = (counts[k-h] > 0) ? std::sqrt(counts[k-h]) : 0.0;
        err_der[idx] = std::sqrt(eU*eU + eD*eD) / (2.0 * DER_STEP);
    }

    // ── Step 4: find peaks ───────────────────────────────────
    const double MIN_THR_MV = 6.0;   // ignore peaks below this [mV]
    const double MIN_REL    = 0.03;  // peak must be >= 3% of global max
    const double MIN_SEP_MV = 3.0;  // minimum separation between peaks [mV]

    // 9-point box smooth
    std::vector<double> sm(n_der, 0.0);
    for (int i = 0; i < n_der; ++i) {
        double s = 0; int c = 0;
        for (int d = -4; d <= 4; ++d) {
            int j = i + d;
            if (j >= 0 && j < n_der) { s += deriv[j]; ++c; }
        }
        sm[i] = s / c;
    }

    double globalMax = 0;
    for (int i = 0; i < n_der; ++i)
        if (thr_der[i] > MIN_THR_MV && sm[i] > globalMax) globalMax = sm[i];
    if (globalMax <= 0) {
        std::cerr << "[CAL] ERROR: no peaks above threshold.\n"; return res;
    }

    std::vector<std::pair<double,double>> rawPeaks;
    for (int i = 1; i < n_der - 1; ++i) {
        if (thr_der[i] <= MIN_THR_MV) continue;
        if (sm[i] < globalMax * MIN_REL) continue;
        if (sm[i] > sm[i-1] && sm[i] >= sm[i+1])
            rawPeaks.push_back({thr_der[i], sm[i]});
    }
    if (rawPeaks.empty()) {
        std::cerr << "[CAL] ERROR: no local maxima found.\n"; return res;
    }

    // Deduplication: if two peaks are closer than MIN_SEP_MV, keep the taller one
    std::sort(rawPeaks.begin(), rawPeaks.end());
    std::vector<std::pair<double,double>> dedupPeaks;
    dedupPeaks.push_back(rawPeaks[0]);
    for (size_t i = 1; i < rawPeaks.size(); ++i) {
        if (rawPeaks[i].first - dedupPeaks.back().first < MIN_SEP_MV) {
            if (rawPeaks[i].second > dedupPeaks.back().second)
                dedupPeaks.back() = rawPeaks[i];
        } else {
            dedupPeaks.push_back(rawPeaks[i]);
        }
    }

    std::cout << "  [CAL] Local maxima (" << dedupPeaks.size() << "): ";
    for (auto& [pos, ht] : dedupPeaks)
        std::cout << std::fixed << std::setprecision(1) << pos << "mV ";
    std::cout << "\n";

    // Dominant peak = tallest (assumed to be 1 p.e. at low light level)
    int dominantIdx = 0;
    for (int i = 1; i < (int)dedupPeaks.size(); ++i)
        if (dedupPeaks[i].second > dedupPeaks[dominantIdx].second)
            dominantIdx = i;
    double dominantPos = dedupPeaks[dominantIdx].first;
    std::cout << "  [CAL] Dominant peak (1 p.e. candidate): "
              << std::fixed << std::setprecision(2) << dominantPos << " mV\n";

    // ── Step 5+6: assign n_pe via gain_raw estimate ──────────
    //
    // gain_raw = gap between 1 p.e. and 2 p.e. peaks.
    // n_pe = round((pos - offset_guess) / gain_raw).
    // Peaks that do not land within 0.4 p.e. of an integer are discarded.
    // This prevents a spurious peak between 2 p.e. and 3 p.e. from shifting
    // all higher-n_pe assignments by +1 (which would give a wrong gain).
    double gain_raw = dominantPos;  // fallback if only one peak
    for (auto& [pos, ht] : dedupPeaks) {
        if (pos > dominantPos + MIN_SEP_MV) {
            gain_raw = pos - dominantPos; break;
        }
    }
    if (gain_raw < 3.0 || gain_raw > 60.0) {
        std::cerr << "[CAL] WARNING: raw gain " << gain_raw
                  << " mV out of [3,60] — using fallback\n";
        gain_raw = dominantPos;
    }
    double offset_guess = dominantPos - gain_raw;
    std::cout << "  [CAL] gain_raw = " << std::fixed << std::setprecision(2)
              << gain_raw << " mV/p.e.   offset_guess = " << offset_guess << " mV\n";

    struct PeakWithNpe { double pos; double ht; int npe; };
    std::vector<PeakWithNpe> assignedPeaks;
    for (auto& [pos, ht] : dedupPeaks) {
        double n_raw = (pos - offset_guess) / gain_raw;
        int    npe   = (int)std::round(n_raw);
        if (npe < 1) continue;
        if (std::abs(n_raw - npe) > 0.4) continue;
        assignedPeaks.push_back({pos, ht, npe});
    }

    // Remove duplicates: same n_pe → keep taller peak
    std::sort(assignedPeaks.begin(), assignedPeaks.end(),
              [](const PeakWithNpe& a, const PeakWithNpe& b){ return a.npe < b.npe; });
    {
        std::vector<PeakWithNpe> dedup2;
        for (auto& pk : assignedPeaks) {
            if (!dedup2.empty() && dedup2.back().npe == pk.npe) {
                if (pk.ht > dedup2.back().ht) dedup2.back() = pk;
            } else {
                dedup2.push_back(pk);
            }
        }
        assignedPeaks = std::move(dedup2);
    }

    if (assignedPeaks.empty()) {
        std::cerr << "[CAL] ERROR: no peaks with valid n_pe assignment.\n"; return res;
    }

    int nFit = std::min((int)assignedPeaks.size(), 7);
    std::vector<double> pkPos;
    std::vector<int>    npeList;
    for (int i = 0; i < nFit; ++i) {
        pkPos  .push_back(assignedPeaks[i].pos);
        npeList.push_back(assignedPeaks[i].npe);
    }

    // Iterative outlier removal: at most 3 passes, dominant (n_pe=1) protected
    if (nFit >= 3) {
        for (int iter = 0; iter < 3 && nFit >= 3; ++iter) {
            TGraphErrors grCheck;
            for (int i = 0; i < nFit; ++i)
                grCheck.SetPoint(i, (double)npeList[i], pkPos[i]);
            TF1 fCheck("fCheck_tmp2", "pol1", 0.5, npeList.back()+0.5);
            grCheck.Fit(&fCheck, "RQ");
            double m0 = fCheck.GetParameter(1), q0 = fCheck.GetParameter(0);
            double maxRes = 0; int maxIdx = -1;
            for (int i = 0; i < nFit; ++i) {
                if (npeList[i] == 1) continue;
                double ri = std::abs(pkPos[i] - (q0 + npeList[i]*m0));
                if (ri > maxRes) { maxRes = ri; maxIdx = i; }
            }
            if (maxRes < 3.0 || maxIdx < 0) break;
            pkPos  .erase(pkPos  .begin() + maxIdx);
            npeList.erase(npeList.begin() + maxIdx);
            nFit = (int)pkPos.size();
        }
    }

    std::cout << "  [CAL] Peaks used for fit (" << nFit << "): ";
    for (int i = 0; i < nFit; ++i)
        std::cout << std::fixed << std::setprecision(1) << pkPos[i]
                  << "(n=" << npeList[i] << ") ";
    std::cout << "mV\n";

    // ── Step 7: linear fit peak_position vs n_pe ────────────
    TGraphErrors* grCal = new TGraphErrors();
    for (int i = 0; i < nFit; ++i)
        grCal->SetPoint(i, (double)npeList[i], pkPos[i]);
    double npe_max = (double)npeList[nFit-1];
    TF1* fLin = new TF1(Form("fLin_%s", tag.c_str()), "pol1", 0.5, npe_max+0.5);
    if (nFit == 1) { fLin->FixParameter(0, 0.0); fLin->SetParameter(1, pkPos[0]); }
    grCal->Fit(fLin, "RQ");
    res.q  = fLin->GetParameter(0);
    res.m  = fLin->GetParameter(1);
    res.ok = (res.m > 0);

    const double GAIN_MIN = 5.0, GAIN_MAX = 50.0;
    if (res.m < GAIN_MIN || res.m > GAIN_MAX) {
        std::cerr << "[CAL] WARNING: gain = " << res.m
                  << " mV/p.e. outside [" << GAIN_MIN << ", " << GAIN_MAX << "]\n";
        res.ok = false;
    }
    if (std::abs(res.q) > res.m * 0.5)
        std::cerr << "[CAL] WARNING: offset = " << res.q << " mV is large (>gain/2).\n";

    std::cout << "  Calibration:  gain = " << res.m << " mV/p.e."
              << "   offset = " << res.q << " mV"
              << "   (" << nFit << " peaks)\n";

    // ── Canvas ───────────────────────────────────────────────
    std::vector<double> xz(n_pts, 0.0), xzd(n_der, 0.0);
    TGraphErrors* grScan = new TGraphErrors(n_pts,
        thresholds.data(), counts.data(), xz.data(), xz.data());
    grScan->SetTitle(Form("Threshold scan   %s;Threshold (mV);Counts", tag.c_str()));
    grScan->SetLineColor(kAzure+1); grScan->SetLineWidth(2);

    TGraphErrors* grDer = new TGraphErrors(n_der,
        thr_der.data(), deriv.data(), xzd.data(), err_der.data());
    grDer->SetTitle(Form("-dN/dV   %s;Threshold (mV);-dN/dV  (mV^{-1})", tag.c_str()));
    grDer->SetLineColor(kRed+1); grDer->SetMarkerColor(kRed+1);
    grDer->SetMarkerStyle(20); grDer->SetMarkerSize(0.4); grDer->SetLineWidth(2);

    TCanvas* cCal = new TCanvas(Form("cCal_%s", tag.c_str()),
        Form("Calibration — %s", tag.c_str()), 950, 1000);
    TPad* pad1 = new TPad(Form("p1_%s",tag.c_str()), "scan",  0.0, 0.5, 1.0, 1.0);
    TPad* pad2 = new TPad(Form("p2_%s",tag.c_str()), "deriv", 0.0, 0.0, 1.0, 0.5);
    for (TPad* p : {pad1, pad2}) {
        p->SetLeftMargin(PAD_LEFT); p->SetRightMargin(PAD_RIGHT);
        p->SetTopMargin(PAD_TOP);   p->SetBottomMargin(PAD_BOTTOM);
        p->SetGrid(); p->Draw();
    }

    pad1->cd();
    grScan->Draw("AL");
    grScan->GetXaxis()->SetLimits(-15.0, 110.0);
    grScan->GetXaxis()->SetRangeUser(-15.0, 110.0);
    double yMax = *std::max_element(counts.begin(), counts.end());
    for (int i = 0; i < nFit; ++i) {
        TLine* lp = new TLine(pkPos[i], 0, pkPos[i], yMax);
        lp->SetLineColor(kOrange+7); lp->SetLineStyle(2); lp->SetLineWidth(2);
        lp->Draw("same");
    }

    pad2->cd();
    grDer->Draw("APL");
    grDer->GetXaxis()->SetLimits(-15.0, 110.0);
    grDer->GetXaxis()->SetRangeUser(-15.0, 110.0);
    double dMax = *std::max_element(deriv.begin(), deriv.end());
    grDer->GetYaxis()->SetRangeUser(-dMax * 0.05, dMax * 1.25);
    for (int i = 0; i < nFit; ++i) {
        TLine* lp = new TLine(pkPos[i], 0, pkPos[i], dMax);
        lp->SetLineColor(kOrange+7); lp->SetLineStyle(2); lp->SetLineWidth(2);
        lp->Draw("same");
    }
    // Gaussian diagnostic fit on the dominant (1 p.e.) peak
    if (!pkPos.empty()) {
        double ctr = pkPos[0];
        double sg  = 1.5;
        for (int i = 0; i < n_der; ++i) {
            if (std::abs(thr_der[i] - ctr) < 0.3) {
                double amp0 = deriv[i];
                for (int j = i+1; j < n_der && thr_der[j] < ctr+20.0; ++j) {
                    if (deriv[j] <= amp0*0.5) {
                        sg = std::max(FIT_SIGMA_MIN,
                             std::min((thr_der[j]-ctr)/1.177, FIT_SIGMA_MAX));
                        break;
                    }
                }
                break;
            }
        }
        TF1* fG = new TF1(Form("fG1_%s",tag.c_str()), "gaus",
                           std::max(5.0, ctr-3*sg), ctr+3*sg);
        fG->SetParameters(dMax*0.5, ctr, sg);
        fG->SetParLimits(1, 5.0, MAX_THR);
        fG->SetParLimits(2, FIT_SIGMA_MIN, FIT_SIGMA_MAX);
        fG->SetLineColor(kGreen+2); fG->SetLineWidth(2);
        grDer->Fit(fG, "RQ");
        fG->Draw("same");
        double mu = fG->GetParameter(1), muE = fG->GetParError(1);
        double si = fG->GetParameter(2), siE = fG->GetParError(2);
        TPaveText* pt = new TPaveText(0.55, 0.68, 0.94, 0.88, "NDC");
        pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
        pt->SetTextFont(42);  pt->SetTextSize(0.038);
        pt->AddText(Form("Gain   = %.3f mV/p.e.", res.m));
        pt->AddText(Form("Offset = %.3f mV",      res.q));
        pt->AddText(Form("#mu_{1} = %.3f #pm %.3f mV",  mu, muE));
        pt->AddText(Form("#sigma_{1} = %.3f #pm %.3f mV", si, siE));
        pt->Draw();
    }

    cCal->Update(); cCal->Modified();
    cCal->SaveAs(ctx.png(Form("calibration_%s.png", tag.c_str())).c_str());

    if (scanOut) {
        scanOut->tag        = tag;
        scanOut->thresholds = thresholds;
        scanOut->counts     = counts;
        scanOut->thr_der    = thr_der;
        scanOut->deriv      = deriv;
        scanOut->pkPos      = pkPos;
        scanOut->gain       = res.m;
        scanOut->offset     = res.q;
    }
    return res;
}

// ════════════════════════════════════════════════════════════
//  drawCalibOverlay()
//  Overlay of threshold scans and normalised derivatives for
//  multiple files. Use after calling calibrateSpectrum() on each.
// ════════════════════════════════════════════════════════════
static void drawCalibOverlay(const std::vector<CalibScanData>& allData, OutCtx& ctx) {
    if (allData.empty()) return;
    int nf = (int)allData.size();

    TCanvas* cScan  = new TCanvas("cCalOverlayScan",
        "Calibration — threshold scan overlay",  1200, 700);
    TCanvas* cDeriv = new TCanvas("cCalOverlayDeriv",
        "Calibration — derivative overlay (normalised)", 1200, 700);

    for (TCanvas* cv : {cScan, cDeriv}) {
        cv->SetGrid(); cv->SetTicks(1,1);
        cv->SetLeftMargin(PAD_LEFT); cv->SetRightMargin(PAD_RIGHT);
        cv->SetTopMargin(PAD_TOP);   cv->SetBottomMargin(0.12f);
    }

    double legH  = std::min(0.07*nf + 0.04, 0.55);
    double legX2 = 1.0 - PAD_RIGHT - 0.01;
    double legX1 = legX2 - 0.32;
    double legY2 = 1.0 - PAD_TOP  - 0.02;
    double legY1 = legY2 - legH;
    double legTxt = std::max(0.024, 0.040 - 0.002*nf);
    TLegend* legS = new TLegend(legX1, legY1, legX2, legY2);
    TLegend* legD = new TLegend(legX1, legY1, legX2, legY2);
    for (auto* leg : {legS, legD}) {
        leg->SetBorderSize(1); leg->SetFillStyle(1001);
        leg->SetFillColor(0);  leg->SetTextSize(legTxt);
    }

    bool firstS = true, firstD = true;
    std::vector<TGraph*> gScans, gDerivs;

    for (int fi = 0; fi < nf; ++fi) {
        const CalibScanData& d = allData[fi];
        int col = calibColor(fi);
        int nS = (int)d.thresholds.size();
        int nD = (int)d.thr_der.size();

        TGraph* gS = new TGraph(nS, d.thresholds.data(), d.counts.data());
        gS->SetTitle(";Threshold (mV);Counts");
        gS->SetLineColor(col); gS->SetLineWidth(2);
        cScan->cd();
        gS->Draw(firstS ? "AL" : "L same");
        if (firstS) {
            gS->GetXaxis()->SetTitleSize(0.055); gS->GetXaxis()->SetLabelSize(0.048);
            gS->GetYaxis()->SetTitleSize(0.055); gS->GetYaxis()->SetLabelSize(0.048);
            gS->GetYaxis()->SetTitleOffset(1.05);
            cScan->SetLogy(1);
            firstS = false;
        }
        legS->AddEntry(gS, d.tag.c_str(), "l");
        gScans.push_back(gS);

        double peak = 0;
        for (int i = 0; i < nD; ++i)
            if (d.thr_der[i] > 5.0 && d.deriv[i] > peak) peak = d.deriv[i];
        std::vector<double> dn(nD);
        for (int i = 0; i < nD; ++i)
            dn[i] = (peak > 0) ? d.deriv[i] / peak : d.deriv[i];

        TGraph* gD = new TGraph(nD, d.thr_der.data(), dn.data());
        gD->SetTitle(";Threshold (mV);-dN/dV  (normalised to peak)");
        gD->SetLineColor(col); gD->SetLineWidth(2);
        cDeriv->cd();
        gD->Draw(firstD ? "AL" : "L same");
        if (firstD) {
            gD->GetXaxis()->SetTitleSize(0.055); gD->GetXaxis()->SetLabelSize(0.048);
            gD->GetYaxis()->SetTitleSize(0.055); gD->GetYaxis()->SetLabelSize(0.048);
            gD->GetYaxis()->SetTitleOffset(1.05);
            firstD = false;
        }
        legD->AddEntry(gD, d.tag.c_str(), "l");
        gDerivs.push_back(gD);

        cDeriv->cd();
        for (double pk : d.pkPos) {
            TLine* lp = new TLine(pk, -0.15, pk, 1.0);
            lp->SetLineColor(col); lp->SetLineStyle(2); lp->SetLineWidth(1);
            lp->Draw("same");
        }
    }

    // Adaptive X range: cut where all files fall below 0.5% of their max
    double xmax = 0;
    for (auto* g : gScans) {
        double ymax = 0;
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y; g->GetPoint(i, x, y);
            if (y > ymax) ymax = y;
        }
        double thr = ymax * 0.005;
        for (int i = g->GetN()-1; i >= 0; --i) {
            double x, y; g->GetPoint(i, x, y);
            if (y > thr) { if (x > xmax) xmax = x; break; }
        }
    }
    double xlo = -2.0, xhi = xmax * 1.03;
    for (auto* g : gScans)  g->GetXaxis()->SetLimits(xlo, xhi);
    for (auto* g : gDerivs) {
        g->GetXaxis()->SetLimits(xlo, xhi);
        g->GetYaxis()->SetRangeUser(-0.20, 1.20);
    }

    cScan ->cd(); legS->Draw();
    cDeriv->cd(); legD->Draw();
    cScan ->Update(); cScan ->Modified();
    cDeriv->Update(); cDeriv->Modified();
    cScan ->SaveAs(ctx.png("calib_all_scans.png").c_str());
    cDeriv->SaveAs(ctx.png("calib_all_derivs.png").c_str());
}
