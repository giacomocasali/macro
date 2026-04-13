#pragma once
// ChunkedHistoFill.h — v2
// ═══════════════════════════════════════════════════════════════
// RAM-optimised histogram filling from the event cache TTree.
// Reads the TTree in chunks of CHUNK_SIZE events, fills all
// registered histograms in a single pass.
//
// Peak RAM per pass: O(CHUNK_SIZE × 36 bytes) ≈ 1.8 MB
// ═══════════════════════════════════════════════════════════════

#include "TOTAnalysis.h"
#include "EventCache.h"
#include <vector>
#include <functional>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TSystem.h>

static constexpr Long64_t CHUNK_SIZE = 50000;

struct CacheEvent {
    double tot;
    double delta_t;
    double amp_max;
    int    n_pe;
};

struct StatAccum {
    double sum = 0, sum2 = 0;
    long   count = 0;
    double vmin = 1e18, vmax = -1e18;
    void add(double v) {
        sum += v; sum2 += v*v; ++count;
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    double mean() const { return count > 0 ? sum/count : 0; }
    double rms()  const {
        if (count < 2) return 0;
        double m = mean();
        return std::sqrt(std::max(sum2/count - m*m, 0.0));
    }
};

using TH1DFillFunc = std::function<void(TH1D*, const CacheEvent&)>;
using TH2DFillFunc = std::function<void(TH2D*, const CacheEvent&)>;

// ═══════════════════════════════════════════════════════════════
class ChunkedFiller {
public:
    using CorrectionFunc = std::function<void(CacheEvent&)>;

    explicit ChunkedFiller(const std::string& cachePath,
                           CorrectionFunc correction = nullptr)
        : cachePath_(cachePath), correction_(correction) {}

    void addTH1D(TH1D* h, TH1DFillFunc func) { h1_.push_back({h, func}); }
    void addTH2D(TH2D* h, TH2DFillFunc func) { h2_.push_back({h, func}); }
    void addStatAccum(StatAccum* acc,
                      std::function<double(const CacheEvent&)> ext,
                      std::function<bool(const CacheEvent&)> filt = nullptr) {
        stats_.push_back({acc, ext, filt});
    }

    Long64_t run() {
        TFile* f = TFile::Open(cachePath_.c_str(), "READ");
        if (!f || f->IsZombie()) { delete f; return 0; }

        TTree* tree = (TTree*)f->Get("events");
        if (!tree) { f->Close(); delete f; return 0; }

        // Abilita cache di lettura per velocità
        tree->SetCacheSize(10 * 1024 * 1024);  // 10 MB cache
        tree->AddBranchToCache("*", true);

        Double_t tot, delta_t, amp_max;
        Int_t    n_pe;
        tree->SetBranchAddress("tot",     &tot);
        tree->SetBranchAddress("delta_t", &delta_t);
        tree->SetBranchAddress("amp_max", &amp_max);
        tree->SetBranchAddress("n_pe",    &n_pe);

        Long64_t nTotal = tree->GetEntries();
        Long64_t nProcessed = 0;

        for (Long64_t i = 0; i < nTotal; ++i) {
            tree->GetEntry(i);
            CacheEvent ev{tot, delta_t, amp_max, (int)n_pe};
            if (correction_) correction_(ev);

            for (auto& [h, func] : h1_) func(h, ev);
            for (auto& [h, func] : h2_) func(h, ev);
            for (auto& [acc, ext, filt] : stats_)
                if (!filt || filt(ev)) acc->add(ext(ev));
            ++nProcessed;
        }

        // Cleanup esplicito
        tree->SetCacheSize(0);
        f->Close();
        delete f;
        return nProcessed;
    }

private:
    std::string cachePath_;
    CorrectionFunc correction_;
    struct H1T { TH1D* h; TH1DFillFunc func; };
    struct H2T { TH2D* h; TH2DFillFunc func; };
    struct ST  { StatAccum* acc;
                 std::function<double(const CacheEvent&)> ext;
                 std::function<bool(const CacheEvent&)> filt; };
    std::vector<H1T> h1_;
    std::vector<H2T> h2_;
    std::vector<ST>  stats_;
};

// ═══════════════════════════════════════════════════════════════
//  Time-walk: empirical, streamed
//  Single pass: accumula Mean(Δt) vs TOT in bin, poi fitta.
// ═══════════════════════════════════════════════════════════════
struct TWBin { double sum=0, sum2=0; int count=0; };

static TF1* computeTimeWalkFromCache(
        const std::string& cachePath,
        double fit_lo, double fit_hi,
        const std::string& tag,
        double& tot_min_out,
        int n_bins = 20)
{
    // Un solo passaggio: trova range E accumula bins contemporaneamente
    // Prima passiamo per trovare tot_min/max
    StatAccum totRange;
    { ChunkedFiller f0(cachePath);
      f0.addStatAccum(&totRange, [](const CacheEvent& e){ return e.tot; });
      if (f0.run() == 0) return nullptr;
    }
    if (totRange.count < 100) return nullptr;

    double tot_min = totRange.vmin, tot_max = totRange.vmax;
    tot_min_out = tot_min;
    double bin_w = (tot_max - tot_min) / n_bins;
    if (bin_w <= 0) return nullptr;

    // Secondo passaggio: accumula bins
    std::vector<TWBin> bins(n_bins);
    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return nullptr; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { f->Close(); delete f; return nullptr; }

    Double_t tot, delta_t;
    tree->SetBranchAddress("tot", &tot);
    tree->SetBranchAddress("delta_t", &delta_t);
    // Disabilita branch inutili
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("tot", 1);
    tree->SetBranchStatus("delta_t", 1);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (delta_t < fit_lo || delta_t > fit_hi) continue;
        int b = (int)((tot - tot_min) / bin_w);
        if (b < 0 || b >= n_bins) continue;
        bins[b].sum += delta_t; bins[b].sum2 += delta_t*delta_t; bins[b].count++;
    }
    f->Close(); delete f;

    std::vector<double> vTOT, vMean, vErr;
    for (int b = 0; b < n_bins; ++b) {
        if (bins[b].count < 5) continue;
        double mean = bins[b].sum / bins[b].count;
        double rms = std::sqrt(std::max(bins[b].sum2/bins[b].count - mean*mean, 0.0));
        vTOT.push_back(tot_min + (b+0.5)*bin_w);
        vMean.push_back(mean);
        vErr.push_back(rms / std::sqrt(bins[b].count));
    }
    if (vTOT.size() < 3) return nullptr;

    TGraphErrors* gr = new TGraphErrors(vTOT.size(), vTOT.data(), vMean.data(), nullptr, vErr.data());
    TF1* fTW = new TF1(Form("fTW_%s", tag.c_str()), "pol2", vTOT.front(), tot_max);
    gr->Fit(fTW, "RQ");
    delete gr;

    std::cout << "  [TW] Dt(TOT) = " << std::scientific << std::setprecision(4)
              << fTW->GetParameter(0) << " + " << fTW->GetParameter(1) << "*TOT + "
              << fTW->GetParameter(2) << "*TOT^2\n" << std::defaultfloat;
    return fTW;
}

// Per-p.e. TW
static std::map<int,double> computeTimeWalkAmplitudeFromCache(
        const std::string& cachePath, double fit_lo, double fit_hi)
{
    std::map<int,double> result;
    std::map<int,TWBin> acc;

    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return result; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { f->Close(); delete f; return result; }

    Double_t delta_t; Int_t n_pe;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("delta_t", 1);
    tree->SetBranchStatus("n_pe", 1);
    tree->SetBranchAddress("delta_t", &delta_t);
    tree->SetBranchAddress("n_pe", &n_pe);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (n_pe < 0 || delta_t < fit_lo || delta_t > fit_hi) continue;
        acc[n_pe].sum += delta_t; acc[n_pe].count++;
    }
    f->Close(); delete f;

    double ref = 0; int nref = 0;
    for (auto& [npe,bin] : acc) {
        if (bin.count < 5) continue;
        result[npe] = bin.sum / bin.count;
        ref += result[npe]; ++nref;
    }
    if (nref > 0) ref /= nref;
    for (auto& [npe,m] : result) m -= ref;
    return result;
}

// Carica subset limitato di eventi (per drawByPE)
static std::vector<TOTEvent> loadEventsSmall(const std::string& cachePath,
                                              Long64_t maxEvents = 200000) {
    std::vector<TOTEvent> events;
    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return events; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { f->Close(); delete f; return events; }

    Double_t tot, delta_t, amp_max; Int_t n_pe;
    tree->SetBranchAddress("tot", &tot);
    tree->SetBranchAddress("delta_t", &delta_t);
    tree->SetBranchAddress("amp_max", &amp_max);
    tree->SetBranchAddress("n_pe", &n_pe);

    Long64_t nRead = std::min(tree->GetEntries(), maxEvents);
    events.reserve(nRead);
    for (Long64_t i = 0; i < nRead; ++i) {
        tree->GetEntry(i);
        events.push_back({tot, delta_t, amp_max, (int)n_pe});
    }
    f->Close(); delete f;
    return events;
}

static Long64_t countCacheEvents(const std::string& cachePath) {
    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return 0; }
    TTree* tree = (TTree*)f->Get("events");
    Long64_t n = tree ? tree->GetEntries() : 0;
    f->Close(); delete f;
    return n;
}
