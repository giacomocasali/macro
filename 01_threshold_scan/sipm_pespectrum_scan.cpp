// final.cpp
// Automatically detects all data.vbias_{V}_L.root files in the current directory.
// Headers:
//   header/Config.h           — parameters, FileInfo, file discovery
//   header/SignalProcessing.h — baseline, filter, peak finder, laser trigger
//   header/Plotting.h         — scan canvas, laser timing canvas
//
// Derivative: pointwise backward difference
//   -dN/dV[i] = -(N[i] - N[i-1]) / STEP_SCAN
// evaluated at threshold = thresholds[i] - STEP_SCAN/2  (bin centre).
// Error (Poisson):
//   err[i] = sqrt(N[i] + N[i-1]) / STEP_SCAN

#include "../header/Config.h"
#include "../header/SignalProcessing.h"
#include "../header/Plotting.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGraphErrors.h>

// ════════════════════════════════════════════════════════════
//  PROCESS A SINGLE FILE
// ════════════════════════════════════════════════════════════
static void processFile(const FileInfo& info,
                        double cutoff_MHz,
                        double t_trig_start,
                        double t_trig_end,
                        bool   doLaser) {

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

    int j_trig_start, j_trig_end;
    triggerWindowIndices(t1, N, t_trig_start, t_trig_end,
                         j_trig_start, j_trig_end, info.tag);

    std::cout << "\n[" << info.tag << "]  "
              << treeCh1->GetEntries() << " events,  fs=" << fs_MHz << " MHz\n"
              << "  Trigger window: [" << t1[j_trig_start] << ", "
              << t1[j_trig_end] << "] ns  (samples "
              << j_trig_start << "-" << j_trig_end << ")\n";

    // ── Threshold grid ───────────────────────────────────────
    int n_pts = (int)std::round((MAX_THR - MIN_THR) / STEP_SCAN) + 1;
    std::vector<double> thresholds(n_pts), counts(n_pts, 0.0);
    for (int k = 0; k < n_pts; ++k) thresholds[k] = MIN_THR + k * STEP_SCAN;

    Long64_t nEntries = treeCh1->GetEntries();

    // ── Event loop ───────────────────────────────────────────
    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 500 == 0) {
            if (gROOT->IsInterrupted()) break;
            std::cout << "\r  Progress: " << i << " / " << nEntries << std::flush;
            gSystem->ProcessEvents();
        }
        treeCh1->GetEntry(i);
        std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
        std::vector<double> af = butterworthLowPass(
            correctBaseline(v_t, v_a, t_base_start, t_base_end),
            cutoff_MHz, fs_MHz);

        for (int k = 0; k < n_pts; ++k) {
            double thr = thresholds[k];
            bool crossed = false;
            if (thr > 0) {
                bool armed = false;
                for (int j = j_trig_start; j <= j_trig_end; ++j) {
                    if (!armed) { if (af[j] < thr * 0.5) armed = true; }
                    else        { if (af[j] > thr) { crossed = true; break; } }
                }
            } else {
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

    // ── Poisson errors on counts ─────────────────────────────
    std::vector<double> errors(n_pts);
    for (int k = 0; k < n_pts; ++k)
        errors[k] = (counts[k] > 0) ? std::sqrt(counts[k]) : 0.0;

    // ── Centred derivative  -dN/dV ───────────────────────────
    // -dN/dV[k] = -(N[k+h] - N[k-h]) / (2*h*STEP_SCAN)
    // where h = STEP_DER / STEP_SCAN (number of bins per half-step).
    // Using STEP_DER = 0.5 mV → h=5 bins: smooths over 10 bins,
    // reducing noise by ~sqrt(10) vs pointwise — much better for fit.
    // Evaluated at thresholds[k] (exact bin centre, no offset needed).
    // Poisson error: sqrt(N[k+h] + N[k-h]) / (2*STEP_DER)
    int h     = (int)std::round(STEP_DER / STEP_SCAN);
    if (h < 1) h = 1;
    int n_der = n_pts - 2*h;
    std::vector<double> thr_der(n_der), deriv(n_der), err_der(n_der);
    for (int k = h; k < n_pts - h; ++k) {
        int idx      = k - h;
        thr_der[idx] = thresholds[k];
        deriv[idx]   = -(counts[k+h] - counts[k-h]) / (2.0 * STEP_DER);
        double eU    = (counts[k+h] > 0) ? std::sqrt(counts[k+h]) : 0.0;
        double eD    = (counts[k-h] > 0) ? std::sqrt(counts[k-h]) : 0.0;
        err_der[idx] = std::sqrt(eU*eU + eD*eD) / (2.0 * STEP_DER);
    }

    // ── Peak finding + sigma-gated candidate selection ───────
    // findFirstPeak uses run-length rise+fall on the smoothed derivative.
    // We iterate over candidates, keeping the first with sigma <= SIGMA_PHYS_MAX.
    const double REFMAX_THR_VAL = 8.0;
    const double SIGMA_PHYS_MAX = 5.0;
    const int    MAX_CAND       = 10;

    // Estimate sigma from right half-max
    auto estimateSigma = [&](int iP) -> double {
        if (iP < 0) return FIT_SIGMA_MAX + 1.0;
        double amp0   = deriv[iP];
        double center = thr_der[iP];
        for (int i = iP + 1; i < n_der; ++i) {
            if (thr_der[i] > center + FIT_SIGMA_MAX * 4.0) break;
            if (deriv[i] <= amp0 * 0.5)
                return std::max(FIT_SIGMA_MIN,
                       std::min((thr_der[i] - center) / 1.177, FIT_SIGMA_MAX));
        }
        return FIT_SIGMA_MAX + 1.0;
    };

    // Build candidate list
    std::vector<int> candidates;
    {
        int start = 0;
        for (int cand = 0; cand < MAX_CAND; ++cand) {
            std::vector<double> thr_slice(thr_der.begin() + start, thr_der.end());
            std::vector<double> der_slice(deriv.begin()   + start, deriv.end());
            int rel = findFirstPeak(thr_slice, der_slice);
            if (rel < 0) break;
            int abs_idx = start + rel;
            if (thr_der[abs_idx] < REFMAX_THR_VAL) { start = abs_idx + 1; continue; }
            candidates.push_back(abs_idx);
            start = abs_idx + 1;
            if (start >= n_der) break;
        }
    }

    // Pick first candidate with sigma <= SIGMA_PHYS_MAX
    int iPeak = -1;
    double bestSigma = 1e30;
    for (int idx : candidates) {
        double sg = estimateSigma(idx);
        std::cout << "  [CAND] peak at " << thr_der[idx]
                  << " mV, sigma_est=" << sg << " mV\n";
        if (sg <= SIGMA_PHYS_MAX) { iPeak = idx; break; }
        if (sg < bestSigma) { bestSigma = sg; iPeak = idx; }
    }
    if (iPeak >= 0 && estimateSigma(iPeak) > SIGMA_PHYS_MAX)
        std::cout << "  [WARN] No narrow peak found; using best candidate at "
                  << thr_der[iPeak] << " mV (sigma_est="
                  << estimateSigma(iPeak) << " mV)\n";

    double sigmaEstVal = estimateSigma(iPeak);
    if (sigmaEstVal > FIT_SIGMA_MAX) sigmaEstVal = FIT_SIGMA_MAX;

    // ── Gaussian fit ─────────────────────────────────────────
    TF1* fitGaus = nullptr;
    double fitMean=0, fitMeanErr=0, fitSigma=0, fitSigmaErr=0;
    double fitAmp=0, fitAmpErr=0, fitChi2=0;
    int    fitNdf=0;

    if (iPeak >= 0) {
        double center   = thr_der[iPeak];
        double amp0     = deriv[iPeak];
        double sigmaEst = sigmaEstVal;
        double fitLo    = std::max(REFMAX_THR_VAL, center - 3.0 * sigmaEst);
        double fitHi    = center + 3.0 * sigmaEst;

        std::vector<double> xe(n_der, 0.0);
        TGraphErrors grFit(n_der,
            thr_der.data(), deriv.data(),
            xe.data(),      err_der.data());
        fitGaus = new TF1(("gaus_" + info.tag).c_str(), "gaus", fitLo, fitHi);
        fitGaus->SetParameters(amp0, center, sigmaEst);
        fitGaus->SetParNames("Amplitude", "Mean", "Sigma");
        fitGaus->SetParLimits(1, REFMAX_THR_VAL, MAX_THR);
        fitGaus->SetParLimits(2, FIT_SIGMA_MIN, FIT_SIGMA_MAX);
        fitGaus->SetLineColor(kGreen+2); fitGaus->SetLineWidth(2);
        grFit.Fit(fitGaus, "RQ");

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
            if (r > 5.0) std::cout << "  *** WARNING: chi2/ndf=" << r << " > 5 ***";
            if (r < 0.3) std::cout << "  *** WARNING: chi2/ndf=" << r << " < 0.3 ***";
        }
        std::cout << "\n";
    } else {
        std::cerr << "  [WARN] No peak found for " << info.tag << "\n";
    }

    // ── Scan canvas ──────────────────────────────────────────
    drawScanCanvas(info, thresholds, counts,
                   thr_der, deriv, err_der,
                   fitGaus, fitMean, fitMeanErr,
                   fitSigma, fitSigmaErr, fitAmp,
                   fitChi2, fitNdf,
                   t_trig_start, t_trig_end);

    // ── Save to ROOT file ────────────────────────────────────
    std::string outname = "scan_" + info.tag + ".root";
    TFile* fout = new TFile(outname.c_str(), "RECREATE");

    double b_vbias=info.vbias, b_filter=info.filter, b_cutoff=cutoff_MHz;
    double b_tbase_start=t_base_start, b_tbase_end=t_base_end;
    double b_ttrig_start=t_trig_start, b_ttrig_end=t_trig_end;
    double b_thr, b_count, b_err_count, b_deriv_val, b_err_deriv;

    TTree* tScan = new TTree("scan", "Threshold scan");
    tScan->Branch("vbias",        &b_vbias,       "vbias/D");
    tScan->Branch("filter",       &b_filter,      "filter/D");
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

    TTree* tDeriv = new TTree("deriv", "Derivative -dN/dV (pointwise)");
    tDeriv->Branch("vbias",        &b_vbias,       "vbias/D");
    tDeriv->Branch("filter",       &b_filter,      "filter/D");
    tDeriv->Branch("cutoff_MHz",   &b_cutoff,      "cutoff_MHz/D");
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
    tFit->Branch("filter",       &b_filter,      "filter/D");
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

    // ── Laser timing canvas ──────────────────────────────────
    if (doLaser) {
        TFile* rawForLaser = TFile::Open(info.path.c_str(), "READ");
        if (rawForLaser && !rawForLaser->IsZombie()) {
            makeLaserCanvas(rawForLaser, info, cutoff_MHz,
                            t_trig_start, t_trig_end, info.tag);
            rawForLaser->Close();
        }
    }
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════
void sipm_pespectrum_scan() {
    gStyle->SetOptStat(0);

    std::vector<FileInfo> files = findInputFiles();
    if (files.empty()) {
        std::cerr << "No data.vbias_{V}_L.root files found.\n"; return;
    }

    double cutoff_MHz;
    std::cout << "\nLow-pass filter cut-off [MHz]: ";
    std::cin >> cutoff_MHz;

    // Read waveform time axis from first file
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

    std::cout << "\nSettings summary:\n"
              << "  LP cut-off     : " << cutoff_MHz    << " MHz\n"
              << "  Trigger window : [" << t_trig_start << ", "
                                        << t_trig_end   << "] ns\n\n";

    // ── Laser canvas prompt ──────────────────────────────────
    bool doLaser = false;
    {
        char ans = 0;
        while (ans != 'y' && ans != 'n') {
            std::cout << "Produce laser timing canvas? [y/n]: ";
            std::cin >> ans;
            if (ans != 'y' && ans != 'n')
                std::cout << "  Please enter 'y' or 'n'.\n";
        }
        doLaser = (ans == 'y');
    }
    std::cout << "\n";

    for (auto& f : files)
        processFile(f, cutoff_MHz, t_trig_start, t_trig_end, doLaser);

    std::cout << "\n=== All files processed. ===\n";
}
