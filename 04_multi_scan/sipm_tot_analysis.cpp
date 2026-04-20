/**
 * sipm_tot_analysis_v5_lowram_interactive.cpp
 *
 * Compile: .L sipm_tot_analysis_v5_lowram_interactive.cpp+
 * Run:     sipm_tot_analysis_v5_interactive()
 *
 * Same analysis flow as v5_lowram, with a safer interactive ending:
 * canvases are reopened and kept editable without forcing a nested
 * gApplication->Run(), which can look like a freeze in some setups.
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/TimingCorrection.h"
#include "../header/VbiasAnalysis_v2.h"

// Definizione variabile globale modalità analisi
// 0 = ORIGINAL (restrittivo), 1 = LOOSE (permissivo)
int g_analysis_mode = 0;

#include <algorithm>
#include <chrono>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <TStyle.h>
#include <TSystem.h>

static std::vector<double> findCalibratedCutoffsInteractive(int vbias, const std::string& dataDir) {
    std::vector<double> cutoffs;
    std::string prefix = "calib_vbias" + std::to_string(vbias) + "_cut";
    std::string suffix = "mhz.root";

    void* dirp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp) return cutoffs;

    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        std::string fname(entry);
        if (fname.size() < prefix.size() + suffix.size()) continue;
        if (fname.substr(0, prefix.size()) != prefix) continue;
        if (fname.substr(fname.size() - suffix.size()) != suffix) continue;

        std::string mid = fname.substr(prefix.size(), fname.size() - prefix.size() - suffix.size());
        try {
            double co = std::stod(mid);
            if (co >= 0) cutoffs.push_back(co);
        } catch (...) {}
    }

    gSystem->FreeDirectory(dirp);
    std::sort(cutoffs.begin(), cutoffs.end());
    return cutoffs;
}

void sipm_tot_analysis() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    const std::string dataDir = DATA_DIR;

    std::cout << "\n+==========================================================+\n"
              << "|  SiPM TOT ANALYSIS v5 (LOW-RAM + INTERACTIVE END)       |\n"
              << "+==========================================================+\n"
              << "  Data dir: " << dataDir << "\n\n";

    // ── MODE SELECTION ──────────────────────────────────────
    std::cout << "\n+--- ANALYSIS MODE ---+\n"
              << "|  0 = ORIGINAL (restrictive, high quality)\n"
              << "|  1 = LOOSE    (permissive, more statistics)\n"
              << "+---------------------+\n";
    {
        std::string line;
        while (true) {
            std::cout << "Select mode [0/1]: " << std::flush;
            if (!std::getline(std::cin, line)) { std::exit(1); }
            if (line == "0") { g_analysis_mode = 0; break; }
            if (line == "1") { g_analysis_mode = 1; break; }
            std::cerr << "  [!] Enter 0 or 1.\n";
        }
        std::cout << "  --> Mode: " 
                  << (g_analysis_mode == 0 ? "ORIGINAL" : "LOOSE") << "\n\n";
    }

    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        while (true) {
            std::cout << prompt << std::flush;
            if (!std::getline(std::cin, line)) {
                std::cerr << "\n[ERROR] stdin closed or unexpected EOF.\n";
                std::exit(1);
            }
            auto b = line.find_first_not_of(" \t\r\n");
            if (b != std::string::npos) {
                line = line.substr(b);
                return line;
            }
        }
    };

    auto readDouble = [&readLine](const std::string& prompt) -> double {
        while (true) {
            std::string line = readLine(prompt);
            try {
                std::size_t pos = 0;
                double val = std::stod(line, &pos);
                while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos]))) ++pos;
                if (pos == line.size()) return val;
            } catch (...) {}
            std::cerr << "  [!] Invalid value, try again.\n";
        }
    };

    auto readYN = [&readLine](const std::string& prompt) -> bool {
        while (true) {
            std::string line = readLine(prompt);
            if (!line.empty()) {
                char c = static_cast<char>(std::tolower(static_cast<unsigned char>(line[0])));
                if (c == 'y') return true;
                if (c == 'n') return false;
            }
            std::cerr << "  [!] Type y or n.\n";
        }
    };

    std::vector<int> vbiasList;
    {
        std::string line = readLine("Vbias (es. 53 54 55):\n> ");
        std::stringstream ss(line);
        int v = 0;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) {
        std::cerr << "[ERROR] No Vbias.\n";
        return;
    }

    std::vector<double> fracs_pe;
    {
        std::string line = readLine("LET thresholds (p.e., es. 0.5 1.0 2.0):\n> ");
        std::stringstream ss(line);
        double v = 0.0;
        while (ss >> v) if (v > 0) fracs_pe.push_back(v);
    }
    if (fracs_pe.empty()) fracs_pe.push_back(1.0);

    std::map<int, CalibResult> calMap;
    for (int vbias : vbiasList) {
        auto cutoffs = findCalibratedCutoffsInteractive(vbias, dataDir);
        if (cutoffs.empty()) {
            std::cerr << "  [WARN] Vbias=" << vbias << ": no cal.\n";
            continue;
        }

        double chosen = cutoffs.back();
        if (cutoffs.size() > 1) {
            std::cout << "  Vbias=" << vbias << " cutoff: ";
            for (double co : cutoffs) std::cout << static_cast<int>(co) << " ";
            std::cout << "MHz\n";

            std::string line = readLine("  Which one? [" + std::to_string(static_cast<int>(chosen)) + "]: ");
            if (!line.empty()) {
                try {
                    double inp = std::stod(line);
                    double best = chosen;
                    double bestD = std::numeric_limits<double>::max();
                    for (double co : cutoffs) {
                        double d = std::abs(co - inp);
                        if (d < bestD) {
                            bestD = d;
                            best = co;
                        }
                    }
                    chosen = best;
                } catch (...) {}
            }
        }

        CalibResult cal;
        if (!loadCalibration(cal, vbias, chosen, dataDir)) continue;
        if (cal.m <= 0 || cal.m > 100) continue;
        calMap[vbias] = cal;

        std::cout << "  Vbias=" << vbias
                  << " gain=" << cal.m << " mV/p.e. trig=["
                  << cal.t_trig_start << "," << cal.t_trig_end << "]\n";
    }
    if (calMap.empty()) {
        std::cerr << "[ERROR] No calibration.\n";
        return;
    }

    double fit_lo = readDouble("\nFit window start [ns]: ");
    double fit_hi = readDouble("Fit window end   [ns]: ");
    if (fit_lo >= fit_hi) {
        std::cerr << "[ERROR] fit_lo must be < fit_hi.\n";
        return;
    }

    TWMethod tw_method = askTimeWalkMethod();
    bool do_pe = readYN("Per-p.e.? [y/n]: ");
    bool use_filter = readYN("LP filter? [y/n]: ");

    // Build run tag: vbias list + LET thresholds + TW + filter
    // e.g. "vbias53_54_55__let0.10_0.50_1.00__tw__filt"
    {
        // nothing here yet — tag built below
    }
    std::string runTag;
    {
        // Vbias
        runTag = "vbias";
        for (int v : vbiasList) runTag += std::to_string(v) + "_";
        if (!runTag.empty() && runTag.back() == '_') runTag.pop_back();
        // LET
        runTag += "__let";
        for (double f : fracs_pe) runTag += Form("%.2f_", f);
        if (!runTag.empty() && runTag.back() == '_') runTag.pop_back();
        // TW
        runTag += (tw_method == TWMethod::FIT_RESIDUALS) ? "__tw" : "__notw";
        // Filter
        runTag += use_filter ? "__filt" : "__nofilt";
        // Per-p.e. analysis
        runTag += do_pe ? "__byPE" : "__noPE";
    }
    OutCtx ctx = createOutputDirs(runTag);

    std::cout << "\n+==========================================================+\n|  Vbias: ";
    for (int v : vbiasList) std::cout << v << " ";
    std::cout << "V\n|  Fit: [" << fit_lo << "," << fit_hi << "] ns\n|  LET: ";
    for (double f : fracs_pe) std::cout << f << " ";
    std::cout << "p.e.\n|  TW: "
              << (tw_method == TWMethod::FIT_RESIDUALS ? "FitResidualsExp" : "None")
              << " Filter: " << (use_filter ? "yes" : "NO")
              << "\n+==========================================================+\n\n";

    std::map<int, std::map<double, std::pair<double, double>>> sigmaResults;
    auto t0 = std::chrono::steady_clock::now();
    for (int vbias : vbiasList) {
        auto it = calMap.find(vbias);
        if (it == calMap.end()) continue;

        const CalibResult& cal = it->second;
        auto tVb = std::chrono::steady_clock::now();
        auto res = processOneVbias_v2(vbias, fracs_pe, cal.cutoff_MHz,
                                      cal.t_trig_start, cal.t_trig_end,
                                      fit_lo, fit_hi, tw_method, do_pe, use_filter,
                                      ctx, cal);
        if (!res.empty()) sigmaResults[vbias] = res;

        std::cout << "  Vbias=" << vbias << " done in "
                  << std::fixed << std::setprecision(1)
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - tVb).count()
                  << " s\n";
    }

    drawVbiasSummary(sigmaResults, ctx);
    std::cout << "\n+==========================================================+\n|  DONE  "
              << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
              << " s\n|  PNG:  " << ctx.pngDir
              << "\n|  ROOT: " << ctx.rootDir
              << "\n+==========================================================+\n";

    // Keep interactive workflow: reopen canvases, then return to ROOT prompt.
    // No nested gApplication->Run() here, to avoid apparent hangs at the end.
    ctx.reopenAllCanvases(40);
    std::cout << "\n  Canvases reopened and editable.\n"
              << "  You can zoom, fit, edit from GUI, and save manually.\n"
              << "  All canvas .root files: " << ctx.canvasDir << "\n";
}
