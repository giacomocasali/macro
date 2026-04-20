#pragma once
// ChunkedHistoFill.h — v3
// Aggiunge delta_t_lo e delta_t_hi a CacheEvent e ChunkedFiller.
// Il resto invariato rispetto a v2.

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
#include <TParameter.h>

static constexpr Long64_t CHUNK_SIZE = 50000;

struct CacheEvent {
    double tot;
    double delta_t;
    double amp_max;
    int    n_pe;
    double delta_t_lo;  // t_rise(0.2*LET) - t_laser  (TRISE_INVALID se assente)
    double delta_t_hi;  // t_rise(0.8*LET) - t_laser  (TRISE_INVALID se assente)
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

        tree->SetCacheSize(10 * 1024 * 1024);
        tree->AddBranchToCache("*", true);

        Double_t tot, delta_t, amp_max;
        Double_t delta_t_lo = TRISE_INVALID, delta_t_hi = TRISE_INVALID;
        Int_t    n_pe;
        tree->SetBranchAddress("tot",     &tot);
        tree->SetBranchAddress("delta_t", &delta_t);
        tree->SetBranchAddress("amp_max", &amp_max);
        tree->SetBranchAddress("n_pe",    &n_pe);

        // Legge i nuovi branch se presenti (cache v3), altrimenti usa sentinella
        bool has_lo = (tree->GetBranch("delta_t_lo") != nullptr);
        bool has_hi = (tree->GetBranch("delta_t_hi") != nullptr);
        if (has_lo) tree->SetBranchAddress("delta_t_lo", &delta_t_lo);
        if (has_hi) tree->SetBranchAddress("delta_t_hi", &delta_t_hi);

        Long64_t nTotal = tree->GetEntries();
        Long64_t nProcessed = 0;

        for (Long64_t i = 0; i < nTotal; ++i) {
            tree->GetEntry(i);
            CacheEvent ev{tot, delta_t, amp_max, (int)n_pe,
                          has_lo ? delta_t_lo : TRISE_INVALID,
                          has_hi ? delta_t_hi : TRISE_INVALID};
            if (correction_) correction_(ev);

            for (auto& [h, func] : h1_) func(h, ev);
            for (auto& [h, func] : h2_) func(h, ev);
            for (auto& [acc, ext, filt] : stats_)
                if (!filt || filt(ev)) acc->add(ext(ev));
            ++nProcessed;
        }

        tree->SetCacheSize(0);
        f->Close(); delete f;
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
//  Time-walk: FitResiduals con esponenziale negativo, streaming.
//  Invariato rispetto a v2.
// ═══════════════════════════════════════════════════════════════
struct TWBin { double sum=0, sum2=0; int count=0; };

static TF1* computeTimeWalkFromCache(
        const std::string& cachePath,
        double fit_lo, double fit_hi,
        const std::string& tag,
        double& tot_min_out,
        int n_bins = 20)
{
    const double tw_lo = fit_lo - 5.0;
    const double tw_hi = fit_hi + 5.0;

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

    std::vector<TWBin> bins(n_bins);
    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return nullptr; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { f->Close(); delete f; return nullptr; }

    Double_t tot, delta_t;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("tot", 1);
    tree->SetBranchStatus("delta_t", 1);
    tree->SetBranchAddress("tot", &tot);
    tree->SetBranchAddress("delta_t", &delta_t);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (delta_t < tw_lo || delta_t > tw_hi) continue;
        int b = (int)((tot - tot_min) / bin_w);
        if (b < 0 || b >= n_bins) continue;
        bins[b].sum += delta_t; bins[b].sum2 += delta_t*delta_t; bins[b].count++;
    }
    f->Close(); delete f;

    std::vector<double> vTOT, vMean, vErr;
    for (int b = 0; b < n_bins; ++b) {
        if (bins[b].count < 5) continue;
        double mean = bins[b].sum / bins[b].count;
        double rms  = std::sqrt(std::max(bins[b].sum2/bins[b].count - mean*mean, 0.0));
        vTOT .push_back(tot_min + (b+0.5)*bin_w);
        vMean.push_back(mean);
        vErr .push_back(rms / std::sqrt(bins[b].count));
    }
    if (vTOT.size() < 3) return nullptr;

    TGraphErrors* gr = new TGraphErrors(vTOT.size(),
        vTOT.data(), vMean.data(), nullptr, vErr.data());
    TF1* fTW = new TF1(Form("fTW_%s", tag.c_str()),
                        "[0] + [1]*exp(-x/[2])", vTOT.front(), tot_max);
    const double y_tail = vMean.back();
    const double y_head = vMean.front();
    const double tau0   = std::max((vTOT.back() - vTOT.front()) * 0.35, 1e-3);
    fTW->SetParameters(y_tail, y_head - y_tail, tau0);
    fTW->SetParLimits(2, 1e-4, 1e6);
    gr->Fit(fTW, "RQ");
    delete gr;

    std::cout << "  [TW] p0=" << fTW->GetParameter(0)
              << "  p1=" << fTW->GetParameter(1)
              << "  tau=" << fTW->GetParameter(2) << " ns\n";
    return fTW;
}

// Per-p.e. TW — invariato
static std::map<int,double> computeTimeWalkAmplitudeFromCache(
        const std::string& cachePath, double fit_lo, double fit_hi)
{
    std::map<int,double> result;
    std::map<int,TWBin>  acc;

    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return result; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { f->Close(); delete f; return result; }

    Double_t delta_t; Int_t n_pe;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("delta_t", 1);
    tree->SetBranchStatus("n_pe", 1);
    tree->SetBranchAddress("delta_t", &delta_t);
    tree->SetBranchAddress("n_pe",    &n_pe);

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

// loadEventsSmall — invariato
static std::vector<TOTEvent> loadEventsSmall(const std::string& cachePath,
                                              Long64_t maxEvents = 200000) {
    std::vector<TOTEvent> events;
    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return events; }
    TTree* tree = (TTree*)f->Get("events");
    if (!tree) { f->Close(); delete f; return events; }

    Double_t tot, delta_t, amp_max;
    Double_t delta_t_lo = TRISE_INVALID, delta_t_hi = TRISE_INVALID;
    Int_t    n_pe;
    tree->SetBranchAddress("tot",     &tot);
    tree->SetBranchAddress("delta_t", &delta_t);
    tree->SetBranchAddress("amp_max", &amp_max);
    tree->SetBranchAddress("n_pe",    &n_pe);
    bool has_lo = (tree->GetBranch("delta_t_lo") != nullptr);
    bool has_hi = (tree->GetBranch("delta_t_hi") != nullptr);
    if (has_lo) tree->SetBranchAddress("delta_t_lo", &delta_t_lo);
    if (has_hi) tree->SetBranchAddress("delta_t_hi", &delta_t_hi);

    Long64_t nRead = std::min(tree->GetEntries(), maxEvents);
    events.reserve(nRead);
    for (Long64_t i = 0; i < nRead; ++i) {
        tree->GetEntry(i);
        events.push_back({tot, delta_t, amp_max, (int)n_pe,
                          has_lo ? delta_t_lo : TRISE_INVALID,
                          has_hi ? delta_t_hi : TRISE_INVALID});
    }
    f->Close(); delete f;
    return events;
}

static Long64_t countCacheEvents(const std::string& cachePath) {
    TFile* f = TFile::Open(cachePath.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return 0; }
    auto* complete = dynamic_cast<TParameter<int>*>(f->Get("cache_complete"));
    if (!complete || complete->GetVal() != 1) { f->Close(); delete f; return 0; }
    TTree* tree = (TTree*)f->Get("events");
    Long64_t n = tree ? tree->GetEntries() : 0;
    f->Close(); delete f;
    return n;
}
