#pragma once
// EventCache.h — v3
// Aggiunge branch delta_t_lo e delta_t_hi:
//   delta_t_lo = t_rise(0.2*LET) - t_laser  [ns]
//   delta_t_hi = t_rise(0.8*LET) - t_laser  [ns]
//   Sentinella TRISE_INVALID (-999) se il segnale non supera quella soglia.
//
// Schema "v3" — le vecchie cache "v2" saranno cache miss automatici
// (nome file diverso). Rigenera con sipm_tot_analysis.
//
// Backward compat: loadEvents accetta cache senza delta_t_lo/hi
// (old v2) con warning, riempie i campi con TRISE_INVALID.

#include "TOTAnalysis.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>

static constexpr const char* EVENT_CACHE_SCHEMA = "v3";

static std::string eventCachePath(int vbias, double frac_pe, double cutoff_MHz,
                                   double laser_thr = 10.0,
                                   const std::string& dataDir = "../../data",
                                   bool use_filter = true,
                                   const std::string& suffix = "") {
    const char* fsuf = use_filter ? "" : "_nofilt";
    extern int g_analysis_mode;
    const char* msuf = (g_analysis_mode == 1) ? "_loose" : "";
    int lthr_int = (int)std::floor(laser_thr);
    int lthr_dec = (int)std::round((laser_thr - lthr_int) * 10.0);
    if (lthr_dec >= 10) { lthr_int += 1; lthr_dec = 0; }
    return Form("%s/events_%s_vbias%d_let%.2fpe_cut%dmhz_lthr%d.%dmV%s%s%s.root",
                dataDir.c_str(), EVENT_CACHE_SCHEMA, vbias, frac_pe,
                (int)std::round(cutoff_MHz),
                lthr_int, lthr_dec, fsuf, msuf, suffix.c_str());
}

// ── saveEvents: mantenuta per compatibilita' legacy ──────────────────────
static void saveEvents(const std::vector<TOTEvent>& events,
                       int vbias, double frac_pe, double cutoff_MHz,
                       double laser_thr,
                       const std::string& dataDir = "../../data",
                       bool use_filter = true,
                       const std::string& suffix = "") {
    if (events.empty()) return;
    std::string path = eventCachePath(vbias, frac_pe, cutoff_MHz,
                                       laser_thr, dataDir, use_filter, suffix);
    if (!gSystem->AccessPathName(path.c_str())) {
        TFile* fCheck = TFile::Open(path.c_str(), "READ");
        if (fCheck && !fCheck->IsZombie()) {
            auto* complete = dynamic_cast<TParameter<int>*>(fCheck->Get("cache_complete"));
            if (complete && complete->GetVal() == 1) {
                fCheck->Close(); delete fCheck;
                std::cout << "  [EventCache] Already valid, skip: " << path << "\n";
                return;
            }
            fCheck->Close(); delete fCheck;
        }
    }
    TFile* f = new TFile(path.c_str(), "RECREATE");
    if (!f || f->IsZombie()) {
        std::cerr << "[EventCache] Cannot write: " << path << "\n";
        delete f; return;
    }

    Double_t tot, delta_t, amp_max, delta_t_lo, delta_t_hi;
    Double_t laser_thr_saved = laser_thr;
    Int_t    n_pe;
    TTree* t = new TTree("events",
        Form("TOTEvents v3 vbias%d let%.2fpe cut%.0fMHz lthr%.4fmV",
             vbias, frac_pe, cutoff_MHz, laser_thr));
    t->Branch("tot",             &tot,             "tot/D");
    t->Branch("delta_t",         &delta_t,         "delta_t/D");
    t->Branch("amp_max",         &amp_max,         "amp_max/D");
    t->Branch("n_pe",            &n_pe,            "n_pe/I");
    t->Branch("laser_thr_saved", &laser_thr_saved, "laser_thr_saved/D");
    t->Branch("delta_t_lo",      &delta_t_lo,      "delta_t_lo/D");
    t->Branch("delta_t_hi",      &delta_t_hi,      "delta_t_hi/D");

    for (auto& e : events) {
        tot        = e.tot;        delta_t    = e.delta_t;
        amp_max    = e.amp_max;    n_pe       = e.n_pe;
        delta_t_lo = e.delta_t_lo; delta_t_hi = e.delta_t_hi;
        t->Fill();
    }
    f->Write(); f->Close(); delete f;
    std::cout << "  [EventCache] Saved " << events.size()
              << " events -> " << path << "\n";
}

// ── loadEvents ───────────────────────────────────────────────────────────
static bool loadEvents(std::vector<TOTEvent>& events,
                       int vbias, double frac_pe, double cutoff_MHz,
                       double laser_thr,
                       const std::string& dataDir = "../../data",
                       bool use_filter = true,
                       const std::string& suffix = "") {
    std::string path = eventCachePath(vbias, frac_pe, cutoff_MHz,
                                       laser_thr, dataDir, use_filter, suffix);
    if (gSystem->AccessPathName(path.c_str())) return false;

    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) return false;

    TTree* t = (TTree*)f->Get("events");
    if (!t || t->GetEntries() == 0) { f->Close(); delete f; return false; }

    Double_t tot, delta_t, amp_max;
    Double_t delta_t_lo = TRISE_INVALID, delta_t_hi = TRISE_INVALID;
    Int_t    n_pe;
    t->SetBranchAddress("tot",     &tot);
    t->SetBranchAddress("delta_t", &delta_t);
    t->SetBranchAddress("amp_max", &amp_max);
    t->SetBranchAddress("n_pe",    &n_pe);

    Double_t laser_thr_saved = laser_thr;
    bool has_thr = (t->GetBranch("laser_thr_saved") != nullptr);
    bool has_lo  = (t->GetBranch("delta_t_lo")      != nullptr);
    bool has_hi  = (t->GetBranch("delta_t_hi")      != nullptr);

    if (has_thr) t->SetBranchAddress("laser_thr_saved", &laser_thr_saved);
    if (has_lo)  t->SetBranchAddress("delta_t_lo",      &delta_t_lo);
    if (has_hi)  t->SetBranchAddress("delta_t_hi",      &delta_t_hi);

    t->GetEntry(0);
    if (has_thr && std::abs(laser_thr_saved - laser_thr) > 0.1) {
        std::cout << "  [EventCache] laser_thr mismatch — cache miss.\n";
        f->Close(); delete f; return false;
    }
    if (!has_lo || !has_hi)
        std::cout << "  [EventCache] WARNING: cache v2 senza delta_t_lo/hi"
                  << " — rigenera con sipm_tot_analysis per l'analisi slope.\n";

    Long64_t N = t->GetEntries();
    events.reserve((size_t)N);
    events.push_back({tot, delta_t, amp_max, (int)n_pe,
                      has_lo ? delta_t_lo : TRISE_INVALID,
                      has_hi ? delta_t_hi : TRISE_INVALID});
    for (Long64_t i = 1; i < N; ++i) {
        t->GetEntry(i);
        events.push_back({tot, delta_t, amp_max, (int)n_pe,
                          has_lo ? delta_t_lo : TRISE_INVALID,
                          has_hi ? delta_t_hi : TRISE_INVALID});
    }
    f->Close(); delete f;
    std::cout << "  [EventCache] Loaded " << events.size()
              << " events from " << path << "\n";
    return true;
}
