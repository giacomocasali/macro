#pragma once
// VbiasAnalysis_v2_diagnostics.h
// AGGIUNTO: 
//   - branch filter_status per diagnostica efficienza filtri
//   - branch n_events_accepted per statistica
//   - parametri adattivi per soglie basse

#include "TOTAnalysis_adaptive.h"

// Codici filtro
enum FilterReject {
    FILTER_PASS = 0,
    FILTER_EDGE_LEADING = 1,
    FILTER_EDGE_TRAILING = 2,
    FILTER_BASELINE_DIRTY = 3,
    FILTER_HYST_RISE = 4,
    FILTER_PRE_CHECK = 5,
    FILTER_FALL_CONFIRM = 6,
    FILTER_POST_FALL = 7,
    FILTER_END_TRACE = 8,
    FILTER_TOT_RANGE = 9,
    FILTER_NO_LASER = 10
};

// Versione diagnostica CON PARAMETRI ADATTIVI — restituisce codice rigetto
static std::tuple<double,double,int> computeTOT_diag(
    const std::vector<double>& time, const std::vector<double>& amp,
    double threshold, int j_start, int j_end, bool use_filter)
{
    AdaptiveParams p = getAdaptiveParams(threshold, use_filter);
    
    const double edge_limit = p.edge_thr_frac * threshold;
    if (edgeMedian(amp, 0, 100) >= edge_limit)
        return {-1, -1, FILTER_EDGE_LEADING};

    const double thr_lo = threshold * p.hyst_frac;
    double t_rise = -1, t_fall = -1;
    bool armed_rise = false;

    for (int j = j_start; j <= j_end; ++j) {
        if (!armed_rise) {
            if (amp[j] < thr_lo) armed_rise = true;
        } else {
            if (amp[j] >= threshold) {
                double tr = time[j-1] + (threshold - amp[j-1])
                            * (time[j] - time[j-1]) / (amp[j] - amp[j-1]);
                int pre_start = std::max(j_start, j - 10);
                int pre_n = j - pre_start;
                if (pre_n > 0 &&
                    edgeMedian(amp, pre_start, pre_n) >= p.pre_check_frac * threshold) {
                    armed_rise = false;
                    continue;
                }
                t_rise = tr;
                break;
            }
        }
    }
    if (t_rise < 0) return {-1, -1, FILTER_HYST_RISE};

    for (int j = j_start; j <= j_end; ++j) {
        if (time[j] <= t_rise) continue;
        if (amp[j] >= threshold) continue;

        int j_cand = j;
        double t_cand = (j_cand > 0)
            ? time[j_cand-1] + (threshold - amp[j_cand-1])
              * (time[j_cand] - time[j_cand-1]) / (amp[j_cand] - amp[j_cand-1])
            : time[j_cand];

        int j_win_end = std::min(j_cand + p.confirm_window, j_end);
        int n_sub = 0;
        for (int k = j_cand; k <= j_win_end; ++k)
            if (amp[k] < threshold) ++n_sub;

        if (n_sub >= p.min_below) {
            t_fall = t_cand;
            break;
        }
        j = j_win_end;
    }
    if (t_fall < 0) return {-1, -1, FILTER_FALL_CONFIRM};

    {
        int jf_post = 0;
        for (int j = j_start; j <= j_end; ++j)
            if (time[j] >= t_fall) { jf_post = j; break; }

        int post_limit = std::min(jf_post + 50, j_end + 1);
        int post_n = post_limit - jf_post;
        if (post_n > 0 && edgeMedian(amp, jf_post, post_n) >= threshold)
            return {-1, -1, FILTER_POST_FALL};
    }

    {
        const int tail_n = 50;
        int tail_start = std::max(j_start, j_end + 1 - tail_n);
        int tail_len = j_end + 1 - tail_start;
        if (tail_len > 0 && edgeMedian(amp, tail_start, tail_len) >= thr_lo)
            return {-1, -1, FILTER_END_TRACE};
    }

    return {t_rise, t_fall, FILTER_PASS};
}


// Versione con diagnostica E PARAMETRI ADATTIVI
static void collectTOTEvents_fileByFile(
    const std::map<int,std::string>& runFiles,
    double cutoff_MHz, double fs_MHz,
    int j_trig_start, int j_trig_end,
    double let_thr, const CalibResult& cal,
    const std::string& tag, const std::string& dataDir,
    int vbias, double frac_pe, bool use_filter,
    const std::string& cache_suffix = "",
    TH2D* hPersA = nullptr, TH2D* hPersB = nullptr)
{
    const int N = 1024;

    // STAMPA PARAMETRI ADATTIVI
    printAdaptiveParams(let_thr, use_filter);

    std::string cachePath = eventCachePath(vbias, frac_pe, cutoff_MHz,
                                            cal.laser_thr, dataDir, use_filter, cache_suffix);
    std::string tmpCachePath = cachePath + ".tmp";
    if (!gSystem->AccessPathName(tmpCachePath.c_str())) gSystem->Unlink(tmpCachePath.c_str());
    TFile* fOut = new TFile(tmpCachePath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) { delete fOut; return; }

    Double_t tot_br, delta_t_br, amp_max_br, laser_thr_br = cal.laser_thr;
    Int_t n_pe_br, filter_status_br;
    Long64_t n_events_accepted_br;
    TTree* tOut = new TTree("events", Form("TOT vbias%d let%.2fpe", vbias, frac_pe));
    tOut->Branch("tot", &tot_br, "tot/D");
    tOut->Branch("delta_t", &delta_t_br, "delta_t/D");
    tOut->Branch("amp_max", &amp_max_br, "amp_max/D");
    tOut->Branch("n_pe", &n_pe_br, "n_pe/I");
    tOut->Branch("laser_thr_saved", &laser_thr_br, "laser_thr_saved/D");
    tOut->Branch("filter_status", &filter_status_br, "filter_status/I");
    tOut->SetAutoSave(50000);

    // Contatori diagnostica
    long totalEntries = 0, nAccepted = 0;
    long rejectCount[11] = {0}; // indici 0-10

    AdaptiveParams p = getAdaptiveParams(let_thr, use_filter);

    for (auto& [r, path] : runFiles) {
        TFile* fTmp = TFile::Open(path.c_str(), "READ");
        if (fTmp && !fTmp->IsZombie()) {
            TTree* tr = (TTree*)fTmp->Get("ch1");
            if (tr) totalEntries += tr->GetEntries();
            fTmp->Close();
        }
        delete fTmp;
    }

    ProgressBar bar(totalEntries, "TOT " + tag + (use_filter ? "" : " [NO FILT]"));
    long globalIdx = 0;
    bool interrupted = false;

    for (auto& [runNum, path] : runFiles) {
        TFile* fIn = TFile::Open(path.c_str(), "READ");
        if (!fIn || fIn->IsZombie()) { delete fIn; continue; }

        TTree* treeCh1 = (TTree*)fIn->Get("ch1");
        TTree* treeLaser = (TTree*)fIn->Get("laser");
        if (!treeCh1 || !treeLaser) { fIn->Close(); delete fIn; continue; }

        Double_t t1[N], a1[N], tL[N], aL[N];
        treeCh1->SetBranchAddress("time", t1);
        treeCh1->SetBranchAddress("amplitude", a1);
        treeLaser->SetBranchAddress("time", tL);
        treeLaser->SetBranchAddress("amplitude", aL);

        treeCh1->SetCacheSize(2 * 1024 * 1024);
        treeLaser->SetCacheSize(2 * 1024 * 1024);

        Long64_t nFile = treeCh1->GetEntries();
        Long64_t nLaser = treeLaser->GetEntries();
        if (nLaser < nFile) nFile = nLaser;

        for (Long64_t i = 0; i < nFile; ++i) {
            if (bar.update(globalIdx, rejectCount[FILTER_NO_LASER],
                          rejectCount[FILTER_FALL_CONFIRM] + rejectCount[FILTER_HYST_RISE],
                          nAccepted)) {
                interrupted = true;
                goto done;
            }
            ++globalIdx;

            treeCh1->GetEntry(i);
            treeLaser->GetEntry(i);

            double t_laser = laserTriggerTime(tL, aL, N, cal.laser_thr);
            if (t_laser < -900.0) {
                ++rejectCount[FILTER_NO_LASER];
                continue;
            }

            std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
            bool bl_ok = true;
            std::vector<double> af = correctBaseline(v_t, v_a,
                BASELINE_START, BASELINE_END, p.baseline_max_rms, &bl_ok);
            if (!bl_ok) {
                ++rejectCount[FILTER_BASELINE_DIRTY];
                continue;
            }
            if (use_filter) af = butterworthLowPass(af, cutoff_MHz, fs_MHz);

            double amp_max = *std::max_element(af.begin()+j_trig_start,
                                                af.begin()+j_trig_end+1);

            auto [t_rise, t_fall, filt_code] = computeTOT_diag(v_t, af, let_thr,
                j_trig_start, j_trig_end, use_filter);

            if (filt_code != FILTER_PASS) {
                ++rejectCount[filt_code];
                continue;
            }

            double tot = t_fall - t_rise;
            double delta_t = t_rise - t_laser;
            if (tot <= 0.0 || tot >= 150.0) {
                ++rejectCount[FILTER_TOT_RANGE];
                continue;
            }

            int n_pe = estimatePE(amp_max, cal);

            fOut->cd();
            tot_br = tot;
            delta_t_br = delta_t;
            amp_max_br = amp_max;
            n_pe_br = n_pe;
            filter_status_br = FILTER_PASS;
            tOut->Fill();
            ++nAccepted;

            auto fillPers = [&](TH2D* h, double tot_lo, double tot_hi) {
                if (!h || tot < tot_lo || tot >= tot_hi) return;
                int ns = (int)v_t.size();
                for (int j = 0; j < ns; ++j) h->Fill(v_t[j] - t_laser, af[j]);
            };
            fillPers(hPersA, 40.0, 43.0);
            fillPers(hPersB, 50.0, 60.0);
        }

        treeCh1->SetCacheSize(0);
        treeLaser->SetCacheSize(0);
        fIn->Close();
        delete fIn;
    }

done:
    bar.done();
    fOut->cd();
    
    // SALVA n_events_accepted come parametro
    n_events_accepted_br = nAccepted;
    TParameter<Long64_t> paramN("n_events_accepted", n_events_accepted_br);
    paramN.Write();
    
    if (!interrupted) {
        TParameter<int> complete("cache_complete", 1);
        complete.Write();
    }
    fOut->Write();
    fOut->Close();
    delete fOut;

    // STAMPA DIAGNOSTICA
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║           FILTER EFFICIENCY DIAGNOSTIC                   ║\n";
    std::cout << "╠════════════════════════════════════════════════════════════╣\n";
    std::cout << "║ Total events:       " << std::setw(10) << totalEntries << "                          ║\n";
    std::cout << "║ Accepted:           " << std::setw(10) << nAccepted
              << " (" << std::setw(5) << std::fixed << std::setprecision(1)
              << 100.0*nAccepted/totalEntries << "%)            ║\n";
    std::cout << "╠════════════════════════════════════════════════════════════╣\n";
    std::cout << "║ REJECTION BREAKDOWN:                                      ║\n";
    std::cout << "║ No laser:           " << std::setw(10) << rejectCount[FILTER_NO_LASER]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_NO_LASER]/totalEntries << "%)            ║\n";
    std::cout << "║ Dirty baseline:     " << std::setw(10) << rejectCount[FILTER_BASELINE_DIRTY]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_BASELINE_DIRTY]/totalEntries << "%)            ║\n";
    std::cout << "║ Edge leading:       " << std::setw(10) << rejectCount[FILTER_EDGE_LEADING]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_EDGE_LEADING]/totalEntries << "%)            ║\n";
    std::cout << "║ Edge trailing:      " << std::setw(10) << rejectCount[FILTER_EDGE_TRAILING]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_EDGE_TRAILING]/totalEntries << "%)            ║\n";
    std::cout << "║ Hysteresis rise:    " << std::setw(10) << rejectCount[FILTER_HYST_RISE]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_HYST_RISE]/totalEntries << "%)            ║\n";
    std::cout << "║ Pre-check:          " << std::setw(10) << rejectCount[FILTER_PRE_CHECK]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_PRE_CHECK]/totalEntries << "%)            ║\n";
    std::cout << "║ Fall confirm:       " << std::setw(10) << rejectCount[FILTER_FALL_CONFIRM]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_FALL_CONFIRM]/totalEntries << "%)            ║\n";
    std::cout << "║ Post-fall:          " << std::setw(10) << rejectCount[FILTER_POST_FALL]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_POST_FALL]/totalEntries << "%)            ║\n";
    std::cout << "║ End trace:          " << std::setw(10) << rejectCount[FILTER_END_TRACE]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_END_TRACE]/totalEntries << "%)            ║\n";
    std::cout << "║ TOT range:          " << std::setw(10) << rejectCount[FILTER_TOT_RANGE]
              << " (" << std::setw(5) << 100.0*rejectCount[FILTER_TOT_RANGE]/totalEntries << "%)            ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";

    if (interrupted) {
        gSystem->Unlink(tmpCachePath.c_str());
        std::cout << "  [Cache] interrupted: cache not updated.\n";
        return;
    }
    if (!gSystem->AccessPathName(cachePath.c_str())) gSystem->Unlink(cachePath.c_str());
    if (gSystem->Rename(tmpCachePath.c_str(), cachePath.c_str()) != 0) {
        std::cerr << "  [Cache] rename failed: " << tmpCachePath << " -> " << cachePath << "\n";
        return;
    }
}
