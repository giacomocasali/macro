// NOTE: You will need the header files you already included below 
// (Config.h, OutputManager.h, etc.) present in the "../header/" directory.
// No other new documents are required.

/**
 * sipm_tot_analysis_v5.cpp
 * ========================
 * Multi-run, multi-Vbias SiPM TOT analysis — versione con cache.
 *
 * Differenze rispetto a v4:
 * - La calibrazione viene letta da  ../../data/calib_vbias<V>_cut<C>MHz.root
 * se disponibile (prodotto da sipm_calibrate.cpp), altrimenti viene
 * eseguita sul primo run e salvata automaticamente per la prossima volta.
 * - Gli eventi TOT vengono salvati in  ../../data/events_vbias<V>_let<L>pe_cut<C>MHz.root
 * dopo il primo loop. Le esecuzioni successive con gli stessi parametri
 * li ricaricano in <1 s invece di riscorrere 1M waveform (~170 s).
 * - Le canvas delle slices sono rimosse (inutili con poca statistica).
 *
 * Flusso di lavoro consigliato:
 * 1. Esegui sipm_calibrate() una volta per ogni Vbias × cutoff.
 * 2. Esegui sipm_tot_analysis_v5() per l'analisi completa.
 * Al primo run crea la cache eventi; dal secondo in poi è veloce.
 * 3. Se cambi fit window o time walk: nessun ricalcolo necessario
 * (gli eventi sono già in cache).
 * 4. Se cambi LET o cutoff: cache miss automatico, loop rieseguito.
 *
 * Compile: .L sipm_tot_analysis_v5.cpp+
 * Run:     sipm_tot_analysis_v5()
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
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
void sipm_tot_analysis_v5()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();
    const std::string dataDir = DATA_DIR;

    std::cout << "\n+==========================================================+\n"
              << "|    SiPM TOT ANALYSIS v5 -- cached calib + events        |\n"
              << "+==========================================================+\n"
              << "  Calibration cache: ../../data/calib_vbias<V>_cut<C>MHz.root\n"
              << "  Event cache:       ../../data/events_vbias<V>_let<L>pe_cut<C>MHz.root\n"
              << "  Run sipm_calibrate() first if no calib cache exists.\n\n";

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

    // NOTE: cutoff_MHz and trigger window are read from the calibration
    // cache file (saved by sipm_calibrate). No need to enter them again.

    // --- 2. LET thresholds (needed for cache peek) ---------------
    std::vector<double> fracs_pe;
    {
        std::string line = readLine("LET thresholds  (p.e., e.g.  0.5 1.0 2.0):\n> ");
        std::stringstream ss(line); double v;
        while (ss >> v) if (v > 0) fracs_pe.push_back(v);
    }
    if (fracs_pe.empty()) fracs_pe.push_back(1.0);

    // --- 3. Peak finder dalla cache (se disponibile) -------------
    // Carica la cache eventi per ogni Vbias × LET e stima la
    // posizione del picco Δt PRIMA di chiedere la fit window.
    // Così l'utente sa già dove sta il picco quando inserisce i limiti.
    std::cout << "\n--- Peak finder (lettura dalla cache, se disponibile) ---\n";
    for (int vbias : vbiasList) {
        CalibResult calPk;
        bool found = false;
        for (double co : {200.0, 150.0, 100.0, 50.0, 300.0, 250.0})
            if (loadCalibration(calPk, vbias, co, dataDir)) { found = true; break; }
        if (!found) continue;

        for (double frac : fracs_pe) {
            auto [tpk, slo, shi] = peakFinderFromCache(
                vbias, frac, calPk.cutoff_MHz, calPk.laser_thr, dataDir);
            if (tpk > 0) {
                std::cout << "  Vbias=" << vbias << "  LET=" << frac << " p.e."
                          << "  →  picco a " << std::fixed << std::setprecision(2)
                          << tpk << " ns"
                          << "   finestra suggerita: ["
                          << slo << ", " << shi << "] ns\n";
            } else {
                std::cout << "  Vbias=" << vbias << "  LET=" << frac
                          << " p.e.  →  cache non trovata (il peak finder"
                          << " sarà eseguito dopo il loop eventi)\n";
            }
        }
    }

    // --- 4. Fit window (ora l'utente ha le informazioni necessarie) ---
    double fit_lo, fit_hi;
    std::cout << "\nFit window (Delta_t)  start [ns]: " << std::flush;
    std::cin >> fit_lo;
    std::cout << "Fit window (Delta_t)  end   [ns]: " << std::flush;
    std::cin >> fit_hi;
    if (fit_lo >= fit_hi) { std::cerr << "[ERROR] Invalid fit window.\n"; return; }

    TWMethod tw_method = askTimeWalkMethod();

    char ans = 0;
    while (ans != 'y' && ans != 'n') {
        std::cout << "Per-p.e. timing analysis? [y/n]: " << std::flush; std::cin >> ans;
    }
    bool do_pe = (ans == 'y');

    // --- 5. Settings summary ------------------------------------
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
              << "|  (cutoff & trigger window: read from calibration cache)\n"
              << "+==========================================================+\n\n";

    // --- 6. Vbias loop ------------------------------------------
    std::map<int, std::map<double, std::pair<double,double>>> sigmaResults;
    auto t0 = std::chrono::steady_clock::now();

    for (int vi = 0; vi < (int)vbiasList.size(); ++vi) {
        int vbias = vbiasList[vi];
        std::cout << "\n[" << (vi+1) << "/" << vbiasList.size()
                  << "]  ===  Vbias = " << vbias << " V  ===\n";

        // Load calibration to get cutoff_MHz, trigger window, laser_thr.
        // processOneVbias will reload it internally too — this peek is
        // just to extract the acquisition parameters before the call.
        CalibResult calPeek;
        bool peekOk = false;
        for (double co : {200.0, 150.0, 100.0, 300.0, 250.0}) {
            if (loadCalibration(calPeek, vbias, co, dataDir)) { peekOk = true; break; }
        }
        if (!peekOk) {
            std::cerr << "[ERROR] No calibration cache found for Vbias=" << vbias
                      << ". Run sipm_calibrate() first.\n";
            continue;
        }

        std::cout << "  [From cache] Cutoff=" << calPeek.cutoff_MHz << " MHz"
                  << "  Trig=[" << calPeek.t_trig_start
                  << ", " << calPeek.t_trig_end << "] ns"
                  << "  Laser thr=" << calPeek.laser_thr << " mV\n";

        auto tVb = std::chrono::steady_clock::now();
        auto res = processOneVbias(vbias, fracs_pe, calPeek.cutoff_MHz,
                                   calPeek.t_trig_start, calPeek.t_trig_end,
                                   fit_lo, fit_hi,
                                   tw_method, do_pe, ctx);
        if (!res.empty()) sigmaResults[vbias] = res;

        double dt = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tVb).count();
        std::cout << "  Vbias=" << vbias << " done in "
                  << std::fixed << std::setprecision(1) << dt << " s\n";
    }

    // --- 7. Cross-Vbias summary ---------------------------------
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