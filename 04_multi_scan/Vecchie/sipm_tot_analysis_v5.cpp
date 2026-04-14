/**
 * sipm_tot_analysis_v5.cpp
 * ========================
 * Multi-run, multi-Vbias SiPM TOT analysis — cached calibration + events.
 *
 * Compile: .L sipm_tot_analysis_v5.cpp+
 * Run:     sipm_tot_analysis_v5()
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/TimingCorrection.h"
#include "../header/VbiasAnalysis.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>
#include <limits>
#include <TApplication.h>

#include <TStyle.h>
#include <TSystem.h>

// ============================================================
//  findCalibratedCutoffs
//  Scansiona dataDir cercando calib_vbias<V>_cut<C>mhz.root.
//  Restituisce i cutoff trovati in ordine CRESCENTE.
// ============================================================
static std::vector<double> findCalibratedCutoffs(int vbias,
                                                  const std::string& dataDir)
{
    std::vector<double> cutoffs;
    const std::string prefix = "calib_vbias" + std::to_string(vbias) + "_cut";
    const std::string suffix = "mhz.root";

    void* dirp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp) return cutoffs;
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        std::string fname(entry);
        if (fname.size() < prefix.size() + suffix.size()) continue;
        if (fname.substr(0, prefix.size()) != prefix) continue;
        if (fname.substr(fname.size() - suffix.size()) != suffix) continue;
        std::string mid = fname.substr(prefix.size(),
            fname.size() - prefix.size() - suffix.size());
        try {
            double co = std::stod(mid);
            if (co > 0) cutoffs.push_back(co);
        } catch (...) {}
    }
    gSystem->FreeDirectory(dirp);
    std::sort(cutoffs.begin(), cutoffs.end());
    return cutoffs;
}

// ============================================================
void sipm_tot_analysis_v5()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();
    const std::string dataDir = DATA_DIR;

    std::cout << "\n+==========================================================+\n"
              << "|    SiPM TOT ANALYSIS v5 -- cached calib + events        |\n"
              << "+==========================================================+\n"
              << "  Data dir:    " << dataDir << "\n"
              << "  Calib cache: calib_vbias<V>_cut<C>mhz.root\n"
              << "  Event cache: events_vbias<V>_let<L>pe_cut<C>mhz_lthr<T>mV.root\n\n";

    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        std::cout << prompt << std::flush;
        while (line.empty()) std::getline(std::cin, line);
        return line;
    };

    // --- 1. Lista Vbias -----------------------------------------
    std::vector<int> vbiasList;
    {
        std::string line = readLine(
            "\nVbias values to analyse (space-separated, e.g.  53 54 55):\n> ");
        std::stringstream ss(line); int v;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) { std::cerr << "[ERROR] No Vbias specified.\n"; return; }

    // --- 2. Threshold LET ---------------------------------------
    std::vector<double> fracs_pe;
    {
        std::string line = readLine("LET thresholds  (p.e., e.g.  0.5 1.0 2.0):\n> ");
        std::stringstream ss(line); double v;
        while (ss >> v) if (v > 0) fracs_pe.push_back(v);
    }
    if (fracs_pe.empty()) fracs_pe.push_back(1.0);

    // --- 3. Selezione cutoff + peak finder ----------------------
    std::cout << "\n--- Available calibrations + Peak finder ---\n";

    std::map<int, CalibResult> calMap;  // calibrazione scelta per ogni Vbias

    for (int vbias : vbiasList) {
        auto cutoffs = findCalibratedCutoffs(vbias, dataDir);

        if (cutoffs.empty()) {
            std::cerr << "  [WARN] Vbias=" << vbias
                      << ": no calibration found in " << dataDir << "\n"
                      << "         Run sipm_calibrate() first.\n";
            continue;
        }

        std::cout << "\n  Vbias=" << vbias << "  Calibrations found: ";
        for (double co : cutoffs) std::cout << (int)co << " MHz  ";
        std::cout << "\n";

        double chosen = cutoffs.back();
        if (cutoffs.size() > 1) {
            std::cout << "  Which cutoff to use? Enter value in MHz [default "
                      << (int)chosen << "]: " << std::flush;
            std::string line;
            std::getline(std::cin, line);
            if (!line.empty()) {
                try {
                    double inp = std::stod(line);
                    double best = chosen, bestDiff = 1e9;
                    for (double co : cutoffs) {
                        if (std::abs(co - inp) < bestDiff) {
                            bestDiff = std::abs(co - inp);
                            best = co;
                        }
                    }
                    chosen = best;
                } catch (...) {}
            }
        }

        CalibResult cal;
        if (!loadCalibration(cal, vbias, chosen, dataDir)) {
            std::cerr << "  [ERROR] Cannot load calibration for Vbias="
                      << vbias << " cutoff=" << chosen << " MHz\n";
            continue;
        }
        if (cal.m <= 0 || cal.m > 100.0) {
            std::cerr << "  [ERROR] Suspicious gain=" << cal.m
                      << " mV/p.e. for Vbias=" << vbias << " — skipping.\n";
            continue;
        }
        calMap[vbias] = cal;

        std::cout << "  --> Using cutoff=" << (int)cal.cutoff_MHz << " MHz"
                  << "  Trigger=[" << std::fixed << std::setprecision(1)
                  << cal.t_trig_start << ", " << cal.t_trig_end << "] ns"
                  << "  Laser thr=" << cal.laser_thr << " mV\n";

        for (double frac : fracs_pe) {
            auto [tpk, slo, shi] = peakFinderFromCache(
                vbias, frac, cal.cutoff_MHz, cal.laser_thr, dataDir);
            if (tpk > 0) {
                std::cout << "  Vbias=" << vbias << "  LET=" << frac << " p.e."
                          << "  --> Delta_t peak at "
                          << std::fixed << std::setprecision(2) << tpk << " ns"
                          << "   suggested fit window: ["
                          << std::setprecision(1) << slo << ", " << shi << "] ns\n";
            } else {
                std::cout << "  Vbias=" << vbias << "  LET=" << frac << " p.e."
                          << "  --> event cache not found."
                          << " Suggested wide window: ["
                          << std::setprecision(1)
                          << cal.t_trig_start << ", " << cal.t_trig_end << "] ns\n";
            }
        }
    }

    if (calMap.empty()) {
        std::cerr << "[ERROR] No calibration available. "
                  << "Run sipm_calibrate() first.\n";
        return;
    }

    // --- 4. Finestra di fit -------------------------------------
    double fit_lo, fit_hi;
    std::cout << "\nFit window (Delta_t)  start [ns]: " << std::flush;
    std::cin >> fit_lo;
    std::cout << "Fit window (Delta_t)  end   [ns]: " << std::flush;
    std::cin >> fit_hi;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (fit_lo >= fit_hi) { std::cerr << "[ERROR] Invalid fit window.\n"; return; }

    TWMethod tw_method = askTimeWalkMethod();

    char ans = 0;
    while (ans != 'y' && ans != 'n') {
        std::cout << "Per-p.e. timing analysis? [y/n]: " << std::flush;
        std::cin >> ans;
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    bool do_pe = (ans == 'y');

    char filt_ans = 0;
    while (filt_ans != 'y' && filt_ans != 'n') {
        std::cout << "Apply low-pass filter? [y/n]: " << std::flush;
        std::cin >> filt_ans;
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    bool use_filter = (filt_ans == 'y');

    if (tw_method == TWMethod::NONE)
        std::cout << "\n  [WARNING] Time-walk correction disabled.\n"
                  << "  Results may be biased for threshold-based timing.\n";

    // --- 5. Riepilogo impostazioni ------------------------------
    std::cout << "\n+==========================================================+\n"
              << "|  Settings summary\n"
              << "|  Vbias list     : ";
    for (int v : vbiasList) std::cout << v << " ";
    std::cout << "V\n"
              << "|  Fit window     : [" << fit_lo << ", " << fit_hi << "] ns\n"
              << "|  LET thresholds : ";
    for (double f : fracs_pe) std::cout << f << " ";
    std::cout << "p.e.\n"
              << "|  Time-walk corr : "
              << (tw_method == TWMethod::EMPIRICAL ? "Empirical" :
                  tw_method == TWMethod::AMPLITUDE ? "Amplitude" : "None") << "\n"
              << "|  Per-p.e. anal. : " << (do_pe ? "yes" : "no") << "\n"
              << "|  Low-pass filter : " << (use_filter ? "yes" : "NO") << "\n";
    for (auto& [vb, cal] : calMap)
        std::cout << "|  Vbias=" << vb
                  << "  cutoff=" << (int)cal.cutoff_MHz << " MHz"
                  << "  trig=[" << cal.t_trig_start
                  << ", " << cal.t_trig_end << "] ns\n";
    std::cout << "+==========================================================+\n\n";

    // --- 6. Loop Vbias ------------------------------------------
    std::map<int, std::map<double, std::pair<double,double>>> sigmaResults;
    auto t0 = std::chrono::steady_clock::now();

    for (int vi = 0; vi < (int)vbiasList.size(); ++vi) {
        int vbias = vbiasList[vi];
        std::cout << "\n[" << (vi+1) << "/" << vbiasList.size()
                  << "]  ===  Vbias = " << vbias << " V  ===\n";

        if (calMap.find(vbias) == calMap.end()) {
            std::cerr << "[ERROR] No calibration for Vbias=" << vbias << "\n";
            continue;
        }
        const CalibResult& cal = calMap[vbias];

        std::cout << "  Cutoff=" << cal.cutoff_MHz << " MHz"
                  << "  Trig=[" << cal.t_trig_start
                  << ", " << cal.t_trig_end << "] ns"
                  << "  Laser thr=" << cal.laser_thr << " mV\n";

        auto tVb = std::chrono::steady_clock::now();

        // NOTA: processOneVbias riceve cal come argomento esplicito.
        // Non ricarica la calibrazione da disco — usa quella già scelta.
        auto res = processOneVbias(vbias, fracs_pe,
                                   cal.cutoff_MHz,
                                   cal.t_trig_start, cal.t_trig_end,
                                   fit_lo, fit_hi,
                                   tw_method, do_pe, use_filter,
                                   ctx,
                                   cal);   // <-- calibrazione passata esplicitamente
        if (!res.empty()) sigmaResults[vbias] = res;

        double dt = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tVb).count();
        std::cout << "  Vbias=" << vbias << " done in "
                  << std::fixed << std::setprecision(1) << dt << " s\n";
    }

    // --- 7. Summary cross-Vbias ---------------------------------
    drawVbiasSummary(sigmaResults, ctx);

    double dtTot = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    std::cout << "\n+==========================================================+\n"
              << "|  DONE   total time : "
              << std::fixed << std::setprecision(1) << dtTot << " s\n"
              << "|  PNG  : " << ctx.pngDir  << "\n"
              << "|  ROOT : " << ctx.rootDir << "\n"
              << "|  All canvases are open for inspection.\n"
              << "|  Type  .q  in the ROOT prompt to exit.\n"
              << "+==========================================================+\n";

    if (gApplication) gApplication->Run(kTRUE);
}
