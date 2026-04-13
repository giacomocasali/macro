#pragma once
// CalibIO.h
// Save and load calibration results to/from a ROOT file.
//
// File format: calib_vbias<V>_cut<C>mhz.root
//   TTree "calib" with one entry: gain, offset, laser_thr,
//   cutoff_MHz, trig_start, trig_end, ok, thresholds[], counts[]
//
// Usage:
//   saveCalibration(cal, thresholds, counts, vbias, cutoff_MHz, dataDir);
//   CalibResult cal;
//   if (loadCalibration(cal, vbias, cutoff_MHz, dataDir)) { ... }

#include "Calibration.h"
#include <string>
#include <vector>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

// Build calibration file path
static std::string calibFilePath(int vbias, double cutoff_MHz,
                                  const std::string& dataDir = "../../data") {
    return Form("%s/calib_vbias%d_cut%dmhz.root",
                dataDir.c_str(), vbias, (int)std::round(cutoff_MHz));
}

// Save calibration to ROOT file
static void saveCalibration(const CalibResult&         cal,
                             const std::vector<double>& thresholds,
                             const std::vector<double>& counts,
                             int                        vbias,
                             double                     cutoff_MHz,
                             const std::string&         dataDir = "../../data") {
    std::string path = calibFilePath(vbias, cutoff_MHz, dataDir);
    TFile* f = new TFile(path.c_str(), "RECREATE");
    if (!f || f->IsZombie()) {
        std::cerr << "[CalibIO] Cannot write: " << path << "\n";
        return;
    }

    Double_t gain         = cal.m;
    Double_t offset       = cal.q;
    Double_t laser_thr    = cal.laser_thr;
    Double_t cutoff_save  = cutoff_MHz;
    Double_t trig_start   = cal.t_trig_start;
    Double_t trig_end     = cal.t_trig_end;
    Int_t    ok           = cal.ok ? 1 : 0;
    Int_t    n_pts        = (Int_t)std::min(thresholds.size(), counts.size());

    static const int MAXPTS = 2000;
    Double_t thr_arr[MAXPTS] = {}, cnt_arr[MAXPTS] = {};
    for (int i = 0; i < std::min(n_pts, MAXPTS); ++i) {
        thr_arr[i] = thresholds[i];
        cnt_arr[i] = counts[i];
    }

    TTree* t = new TTree("calib", Form("Calibration vbias%d cut%.0fMHz",
                                        vbias, cutoff_MHz));
    t->Branch("gain",        &gain,        "gain/D");
    t->Branch("offset",      &offset,      "offset/D");
    t->Branch("laser_thr",   &laser_thr,   "laser_thr/D");
    t->Branch("cutoff_MHz",  &cutoff_save, "cutoff_MHz/D");
    t->Branch("trig_start",  &trig_start,  "trig_start/D");
    t->Branch("trig_end",    &trig_end,    "trig_end/D");
    t->Branch("ok",          &ok,          "ok/I");
    t->Branch("n_pts",       &n_pts,       "n_pts/I");
    t->Branch("thresholds",  thr_arr,      "thresholds[n_pts]/D");
    t->Branch("counts",      cnt_arr,      "counts[n_pts]/D");
    t->Fill();

    f->Write(); f->Close(); delete f;
    std::cout << "  [CalibIO] Saved: " << path << "\n";
}

// Loads calibration from ROOT file.
// Returns true if file exists and is valid, false otherwise
static bool loadCalibration(CalibResult&  cal,
                             int           vbias,
                             double        cutoff_MHz,
                             const std::string& dataDir = "../../data") {
    std::string path = calibFilePath(vbias, cutoff_MHz, dataDir);
    if (gSystem->AccessPathName(path.c_str())) return false;  // file not found
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) return false;

    TTree* t = (TTree*)f->Get("calib");
    if (!t || t->GetEntries() < 1) { f->Close(); return false; }

    Double_t gain = 0, offset = 0, laser_thr = 10.0;
    Double_t cutoff_save = 200.0, trig_start = 95.0, trig_end = 125.0;
    Int_t    ok = 0;
    t->SetBranchAddress("gain",   &gain);
    t->SetBranchAddress("offset", &offset);
    t->SetBranchAddress("ok",     &ok);
    // Optional branches — backward compatible with old cache files
    if (t->GetBranch("laser_thr"))  t->SetBranchAddress("laser_thr",  &laser_thr);
    if (t->GetBranch("cutoff_MHz")) t->SetBranchAddress("cutoff_MHz", &cutoff_save);
    if (t->GetBranch("trig_start")) t->SetBranchAddress("trig_start", &trig_start);
    if (t->GetBranch("trig_end"))   t->SetBranchAddress("trig_end",   &trig_end);
    t->GetEntry(0);
    f->Close(); delete f;

    cal.m            = gain;
    cal.q            = offset;
    cal.ok           = (ok != 0);
    cal.laser_thr    = laser_thr;
    cal.cutoff_MHz   = cutoff_save;
    cal.t_trig_start = trig_start;
    cal.t_trig_end   = trig_end;

    std::cout << "  [CalibIO] Loaded from " << path << "\n"
              << "  Gain = " << cal.m << " mV/p.e."
              << "   Offset = " << cal.q << " mV\n";
    return cal.ok;
}
