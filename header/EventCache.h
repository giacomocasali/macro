#pragma once
// EventCache.h
// TOT event cache: saves/loads the TOTEvent vector to/from a ROOT file.
// Avoids re-scanning all waveforms (~170 s) on repeated runs with same parameters.
//
// File name: events_vbias<V>_let<L>pe_cut<C>mhz_lthr<T>mV[_nofilt].root
// Any parameter change → different filename → automatic cache miss.
// laser_thr is stored in the file as a branch (laser_thr_saved) and verified
// on load. Old files without this branch are accepted with a warning.

#include "TOTAnalysis.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTree.h>

static std::string eventCachePath(int vbias, double frac_pe, double cutoff_MHz,
                                   double laser_thr = 10.0,
                                   const std::string& dataDir = "../../data",
                                   bool use_filter = true,
                                   const std::string& suffix = "") {
    const char* fsuf = use_filter ? "" : "_nofilt";
    // laser_thr rounded to 0.1 mV in the filename to avoid floating-point noise.
    // Old files with integer rounding (lthr<N>mV) will be cache misses → safe.
    int lthr_int = (int)std::floor(laser_thr);
    int lthr_dec = (int)std::round((laser_thr - lthr_int) * 10.0);
    if (lthr_dec >= 10) { lthr_int += 1; lthr_dec = 0; }
    return Form("%s/events_vbias%d_let%.2fpe_cut%dmhz_lthr%d.%dmV%s%s.root",
                dataDir.c_str(), vbias, frac_pe,
                (int)std::round(cutoff_MHz),
                lthr_int, lthr_dec, fsuf, suffix.c_str());
}

static void saveEvents(const std::vector<TOTEvent>& events,
                       int vbias, double frac_pe, double cutoff_MHz,
                       double laser_thr,
                       const std::string& dataDir = "../../data",
                       bool use_filter = true,
                       const std::string& suffix = "") {
    if (events.empty()) return;
    std::string path = eventCachePath(vbias, frac_pe, cutoff_MHz, laser_thr, dataDir, use_filter, suffix);
    TFile* f = new TFile(path.c_str(), "RECREATE");
    if (!f || f->IsZombie()) {
        std::cerr << "[EventCache] Cannot write: " << path << "\n";
        delete f; return;
    }

    Double_t tot, delta_t, amp_max;
    Int_t    n_pe;
    // laser_thr_saved: stored as a per-entry branch (constant value).
    // loadEvents() reads it back and rejects the cache if it differs by >0.1 mV.
    Double_t laser_thr_saved = laser_thr;
    TTree* t = new TTree("events",
        Form("TOTEvents vbias%d let%.2fpe cut%.0fMHz lthr%.4fmV",
             vbias, frac_pe, cutoff_MHz, laser_thr));
    t->Branch("tot",             &tot,             "tot/D");
    t->Branch("delta_t",         &delta_t,         "delta_t/D");
    t->Branch("amp_max",         &amp_max,         "amp_max/D");
    t->Branch("n_pe",            &n_pe,            "n_pe/I");
    t->Branch("laser_thr_saved", &laser_thr_saved, "laser_thr_saved/D");

    for (auto& e : events) {
        tot = e.tot; delta_t = e.delta_t;
        amp_max = e.amp_max; n_pe = e.n_pe;
        t->Fill();
    }

    f->Write(); f->Close(); delete f;
    std::cout << "  [EventCache] Saved " << events.size()
              << " events -> " << path << "\n";
}

// Returns true if cache found AND laser_thr matches, false = cache miss.
static bool loadEvents(std::vector<TOTEvent>& events,
                       int vbias, double frac_pe, double cutoff_MHz,
                       double laser_thr,
                       const std::string& dataDir = "../../data",
                       bool use_filter = true,
                       const std::string& suffix = "") {
    std::string path = eventCachePath(vbias, frac_pe, cutoff_MHz, laser_thr, dataDir, use_filter, suffix);
    if (gSystem->AccessPathName(path.c_str())) return false;

    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) return false;

    TTree* t = (TTree*)f->Get("events");
    if (!t || t->GetEntries() == 0) { f->Close(); delete f; return false; }

    Double_t tot, delta_t, amp_max;
    Int_t    n_pe;
    t->SetBranchAddress("tot",     &tot);
    t->SetBranchAddress("delta_t", &delta_t);
    t->SetBranchAddress("amp_max", &amp_max);
    t->SetBranchAddress("n_pe",    &n_pe);

    // Verify laser_thr_saved if the branch exists.
    // Old files without this branch are accepted with a warning.
    Double_t laser_thr_saved = laser_thr;
    bool has_thr_branch = (t->GetBranch("laser_thr_saved") != nullptr);
    if (has_thr_branch)
        t->SetBranchAddress("laser_thr_saved", &laser_thr_saved);

    t->GetEntry(0);
    if (has_thr_branch && std::abs(laser_thr_saved - laser_thr) > 0.1) {
        std::cout << "  [EventCache] laser_thr mismatch: cache="
                  << laser_thr_saved << " mV, requested="
                  << laser_thr << " mV — cache miss.\n";
        f->Close(); delete f; return false;
    }
    if (!has_thr_branch)
        std::cout << "  [EventCache] WARNING: old cache without laser_thr_saved"
                  << " — accepted. Delete cache to regenerate with verification.\n";

    Long64_t N = t->GetEntries();
    events.reserve((size_t)N);
    // First entry already in memory from the verification read above
    events.push_back({tot, delta_t, amp_max, (int)n_pe});
    for (Long64_t i = 1; i < N; ++i) {
        t->GetEntry(i);
        events.push_back({tot, delta_t, amp_max, (int)n_pe});
    }
    f->Close(); delete f;

    std::cout << "  [EventCache] Loaded " << events.size()
              << " events from " << path << "\n";
    return true;
}
