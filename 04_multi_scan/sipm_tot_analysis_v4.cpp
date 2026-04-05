/**
 * sipm_tot_analysis_v4.cpp
 * ========================
 * Multi-run, multi-Vbias SiPM TOT analysis for ALCOR readout.
 * Files expected: ../../data/data.vbias_<V>_run_<N>.root
 *
 * Compile: .L sipm_tot_analysis_v4.cpp+
 * Run:     sipm_tot_analysis_v4()
 *
 * All physics/plotting logic lives in the headers:
 *   ProgressBar.h    -- live terminal progress bar
 *   VbiasAnalysis.h  -- collectTOTEvents_v4(), processOneVbias()
 *   VbiasSummary.h   -- drawMultiLETOverlay(), drawVbiasSummary()
 *   (+ the existing v3 headers: Calibration, TOTAnalysis, TOTPlotting, ...)
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/Calibration.h"
#include "../header/TimingCorrection.h"
#include "../header/VbiasAnalysis.h"   // pulls in VbiasSummary.h + ProgressBar.h

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <chrono>
#include <limits>

#include <TStyle.h>

// ============================================================
void sipm_tot_analysis_v4()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    // NOTE: ROOT::EnableImplicitMT() is intentionally NOT called.
    //       It is incompatible with the manual SetBranchAddress + GetEntry
    //       event loop: parallel threads would write into the same waveform
    //       buffers, corrupting data silently.

    OutCtx ctx = createOutputDirs();

    std::cout << "\n+==========================================================+\n"
              << "|    SiPM TOT ANALYSIS v4 -- multi-run / multi-Vbias       |\n"
              << "+==========================================================+\n";

    // Helper: reads a full line, skipping blank lines left by previous cin >>
    // This avoids the fragile cin.ignore() + getline pattern.
    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        std::cout << prompt << std::flush;
        // Skip any blank lines left by previous cin >> calls
        while (line.empty()) std::getline(std::cin, line);
        return line;
    };

    // --- 1. Vbias list ------------------------------------------
    std::vector<int> vbiasList;
    {
        std::string line = readLine(
            "\nVbias values to analyse (space-separated, e.g.  53 54 55):\n> ");
        std::stringstream ss(line); int v;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) { std::cerr << "[ERROR] No Vbias specified.\n"; return; }

    // --- 2. Common parameters (asked once for all Vbias) --------
    double cutoff_MHz;
    std::cout << "Low-pass filter cut-off [MHz]: " << std::flush;
    std::cin >> cutoff_MHz;

    double t_trig_start, t_trig_end;
    std::cout << "Trigger window  start [ns]: " << std::flush; std::cin >> t_trig_start;
    std::cout << "Trigger window  end   [ns]: " << std::flush; std::cin >> t_trig_end;
    if (t_trig_start >= t_trig_end) {
        std::cerr << "[ERROR] Invalid trigger window.\n"; return;
    }

    double fit_lo, fit_hi;
    std::cout << "Fit window (Delta_t)  start [ns]: " << std::flush; std::cin >> fit_lo;
    std::cout << "Fit window (Delta_t)  end   [ns]: " << std::flush; std::cin >> fit_hi;
    if (fit_lo >= fit_hi) { std::cerr << "[ERROR] Invalid fit window.\n"; return; }

    std::vector<double> fracs_pe;
    {
        std::string line = readLine("LET thresholds  (p.e., e.g.  0.5 1.0 2.0):\n> ");
        std::stringstream ss(line); double v;
        while (ss >> v) if (v > 0) fracs_pe.push_back(v);
    }
    if (fracs_pe.empty()) fracs_pe.push_back(1.0);

    TWMethod tw_method = askTimeWalkMethod();

    char ans = 0;
    while (ans != 'y' && ans != 'n') {
        std::cout << "Per-p.e. timing analysis? [y/n]: " << std::flush; std::cin >> ans;
    }
    bool do_pe = (ans == 'y');

    // --- 3. Settings summary ------------------------------------
    std::cout << "\n+==========================================================+\n"
              << "|  Settings summary\n"
              << "|  Vbias list     : ";
    for (int v : vbiasList) std::cout << v << " ";
    std::cout << "V\n"
              << "|  LP cut-off     : " << cutoff_MHz << " MHz\n"
              << "|  Trigger window : [" << t_trig_start << ", " << t_trig_end << "] ns\n"
              << "|  Fit window     : [" << fit_lo << ", " << fit_hi << "] ns\n"
              << "|  LET thresholds : ";
    for (double f : fracs_pe) std::cout << f << " ";
    std::cout << "p.e.\n"
              << "|  Time-walk corr : "
              << (tw_method == TWMethod::EMPIRICAL ? "Empirical" :
                  tw_method == TWMethod::AMPLITUDE ? "Amplitude" : "None") << "\n"
              << "|  Per-p.e. anal. : " << (do_pe ? "yes" : "no") << "\n"
              << "+==========================================================+\n\n";

    // --- 4. Vbias loop ------------------------------------------
    // sigmaResults[vbias][frac_pe] = {sigma_ns, sigmaErr_ns}
    std::map<int, std::map<double, std::pair<double,double>>> sigmaResults;
    auto t0 = std::chrono::steady_clock::now();

    for (int vi = 0; vi < (int)vbiasList.size(); ++vi) {
        int vbias = vbiasList[vi];
        std::cout << "\n[" << (vi+1) << "/" << vbiasList.size()
                  << "]  ===  Vbias = " << vbias << " V  ===\n";

        auto tVb = std::chrono::steady_clock::now();
        auto res = processOneVbias(vbias, fracs_pe, cutoff_MHz,
                                   t_trig_start, t_trig_end,
                                   fit_lo, fit_hi,
                                   tw_method, do_pe, ctx);
        if (!res.empty()) sigmaResults[vbias] = res;

        double dt = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tVb).count();
        std::cout << "  Vbias=" << vbias << " done in "
                  << std::fixed << std::setprecision(1) << dt << " s\n";
    }

    // --- 5. Cross-Vbias summary ---------------------------------
    drawVbiasSummary(sigmaResults, ctx);

    double dtTot = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    std::cout << "\n+==========================================================+\n"
              << "|  DONE   total time : "
              << std::fixed << std::setprecision(1) << dtTot << " s\n"
              << "|  PNG  : " << ctx.pngDir  << "\n"
              << "|  ROOT : " << ctx.rootDir << "\n"
              << "+==========================================================+\n";
}
