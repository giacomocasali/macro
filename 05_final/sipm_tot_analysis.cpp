/**
 * sipm_tot_analysis.cpp
 *
 * Compile: .L sipm_tot_analysis.cpp+
 * Run:     sipm_tot_analysis()
 *
 * Analysis flow: interactive directory/mode/vbias/LET selection,
 * calibration load, per-vbias TOT analysis, VbiasSummary plot.
 * Canvases are reopened and kept editable at the end without forcing
 * a nested gApplication->Run().
 *
 * FIX (bug #ODR): g_analysis_mode is now defined ONCE in Config.h
 *   as `inline int g_analysis_mode = 0;`. The old bare definition
 *   `int g_analysis_mode = 0;` on line 21 has been removed. Every
 *   other TU that previously re-defined it with
 *   `#ifndef G_ANALYSIS_MODE_DEFINED / inline int g_analysis_mode = 0; / #endif`
 *   must also drop those blocks and rely on Config.h.
 *
 * FIX (include order): STL and ROOT headers now come BEFORE project
 *   headers so that every project header finds its dependencies
 *   already in scope regardless of include order inside those headers.
 *
 * FIX (g_data_dir_override): g_data_dir_override is already declared
 *   inline in Config.h; no duplicate declaration needed here.
 */

// ── Standard library ────────────────────────────────────────────────────────
#include <algorithm>
#include <chrono>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ── ROOT ─────────────────────────────────────────────────────────────────────
#include <TStyle.h>
#include <TSystem.h>

// ── Project headers ──────────────────────────────────────────────────────────
// Config.h MUST be first: it defines g_analysis_mode (inline) and
// g_data_dir_override (inline), both referenced by the headers below.
#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/TimingCorrection.h"
#include "../header/VbiasAnalysis_v2.h"

// ── findCalibratedCutoffsInteractive ─────────────────────────────────────────
// Scans dataDir for calibration files matching
//   calib_vbias<vbias>_cut<N>mhz.root
// and returns the sorted list of cutoff values found.
static std::vector<double> findCalibratedCutoffsInteractive(
        int vbias, const std::string& dataDir)
{
    std::vector<double> cutoffs;
    const std::string prefix = "calib_vbias" + std::to_string(vbias) + "_cut";
    const std::string suffix = "mhz.root";

    void* dirp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp) return cutoffs;

    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        std::string fname(entry);
        if (fname.size() < prefix.size() + suffix.size()) continue;
        if (fname.substr(0, prefix.size()) != prefix)     continue;
        if (fname.substr(fname.size() - suffix.size()) != suffix) continue;

        const std::string mid = fname.substr(
            prefix.size(),
            fname.size() - prefix.size() - suffix.size());
        try {
            double co = std::stod(mid);
            if (co >= 0) cutoffs.push_back(co);
        } catch (...) {}
    }

    gSystem->FreeDirectory(dirp);
    std::sort(cutoffs.begin(), cutoffs.end());
    return cutoffs;
}

// ── sipm_tot_analysis ─────────────────────────────────────────────────────────
void sipm_tot_analysis()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // ── 0. DATA DIRECTORY SELECTION ─────────────────────────────────────────
    // Scans the parent of DATA_DIR (= "data/") and lets the user pick a
    // subdirectory. The choice is written to g_data_dir_override (Config.h)
    // so that every downstream function that calls effectiveDataDir() picks
    // it up automatically without receiving an explicit argument.
    std::string dataDir;
    {
        // Go one level up from DATA_DIR to find sibling subdirectories.
        std::string root = DATA_DIR;
        const size_t slash = root.find_last_of("/\\");
        if (slash != std::string::npos) root = root.substr(0, slash);

        std::cout << "\n+==========================================================+\n"
                  << "|  SiPM TOT ANALYSIS (LOW-RAM + INTERACTIVE END)           |\n"
                  << "+==========================================================+\n"
                  << "  Dataset root: " << root << "\n\n";

        // List subdirectories.
        std::vector<std::string> subs;
        void* dp = gSystem->OpenDirectory(root.c_str());
        if (dp) {
            const char* ent = nullptr;
            while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
                std::string s(ent);
                if (s == "." || s == "..") continue;
                const std::string full = root + "/" + s;
                FileStat_t st;
                if (gSystem->GetPathInfo(full.c_str(), st) == 0 &&
                    R_ISDIR(st.fMode))
                    subs.push_back(s);
            }
            gSystem->FreeDirectory(dp);
        }
        std::sort(subs.begin(), subs.end());

        if (subs.empty()) {
            std::cerr << "[ERROR] No subdirectories found in " << root << "\n";
            return;
        }

        std::cout << "  Available subdirectories:\n";
        for (size_t i = 0; i < subs.size(); ++i)
            std::cout << "    [" << (i + 1) << "] " << subs[i] << "\n";
        std::cout << "\n  Choose number (or full path): " << std::flush;

        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cerr << "[ERROR] stdin closed.\n";
            return;
        }

        // Trim whitespace.
        const auto b = line.find_first_not_of(" \t\r\n");
        const auto e = line.find_last_not_of(" \t\r\n");
        if (b == std::string::npos) {
            std::cerr << "[ERROR] Empty input.\n";
            return;
        }
        line = line.substr(b, e - b + 1);

        // Numeric index → pick from list; anything else → use as path.
        try {
            size_t pos = 0;
            const int idx = std::stoi(line, &pos);
            if (pos == line.size() && idx >= 1 &&
                idx <= static_cast<int>(subs.size()))
                dataDir = root + "/" + subs[idx - 1];
            else
                dataDir = line;
        } catch (...) {
            dataDir = line;
        }

        if (gSystem->AccessPathName(dataDir.c_str()) != 0) {
            std::cerr << "[ERROR] Directory not accessible: " << dataDir << "\n";
            return;
        }

        // Publish the choice globally so every downstream helper sees it.
        g_data_dir_override = dataDir;
        std::cout << "\n  --> Using: " << dataDir << "\n\n";
    }

    // All interactive prompts below are wrapped in a single try/catch so that
    // a redirected stdin that reaches EOF before all prompts are answered
    // unwinds the stack cleanly instead of calling std::exit(1) while TFiles
    // may be open.
    try {

    // ── 1. ANALYSIS MODE ────────────────────────────────────────────────────
    std::cout << "\n+--- ANALYSIS MODE ---+\n"
              << "|  0 = ORIGINAL (restrictive, high quality)\n"
              << "|  1 = LOOSE    (permissive, more statistics)\n"
              << "+---------------------+\n";
    {
        std::string line;
        while (true) {
            std::cout << "Select mode [0/1]: " << std::flush;
            if (!std::getline(std::cin, line)) {
                std::cerr << "\n[ERROR] stdin closed unexpectedly — aborting.\n";
                throw std::runtime_error("stdin EOF");
            }
            if (line == "0") { g_analysis_mode = 0; break; }
            if (line == "1") { g_analysis_mode = 1; break; }
            std::cerr << "  [!] Enter 0 or 1.\n";
        }
        std::cout << "  --> Mode: "
                  << (g_analysis_mode == 0 ? "ORIGINAL" : "LOOSE") << "\n\n";
    }

    // ── 2. INPUT HELPERS ────────────────────────────────────────────────────

    // readLine: loops until the user types a non-blank line.
    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        while (true) {
            std::cout << prompt << std::flush;
            if (!std::getline(std::cin, line)) {
                std::cerr << "\n[ERROR] stdin closed.\n";
                throw std::runtime_error("stdin EOF");
            }
            const auto b = line.find_first_not_of(" \t\r\n");
            if (b != std::string::npos) {
                return line.substr(b);
            }
        }
    };

    // readDouble: keeps asking until the user enters a valid floating-point number.
    auto readDouble = [&readLine](const std::string& prompt) -> double {
        while (true) {
            const std::string line = readLine(prompt);
            try {
                std::size_t pos = 0;
                const double val = std::stod(line, &pos);
                // Accept only if the entire token is numeric (trailing spaces ok).
                while (pos < line.size() &&
                       std::isspace(static_cast<unsigned char>(line[pos])))
                    ++pos;
                if (pos == line.size()) return val;
            } catch (...) {}
            std::cerr << "  [!] Invalid number, try again.\n";
        }
    };

    // readYN: keeps asking until the user types y or n (case-insensitive).
    auto readYN = [&readLine](const std::string& prompt) -> bool {
        while (true) {
            const std::string line = readLine(prompt);
            if (!line.empty()) {
                const char c = static_cast<char>(
                    std::tolower(static_cast<unsigned char>(line[0])));
                if (c == 'y') return true;
                if (c == 'n') return false;
            }
            std::cerr << "  [!] Type y or n.\n";
        }
    };

    // ── 3. VBIAS LIST ───────────────────────────────────────────────────────
    std::vector<int> vbiasList;
    {
        const std::string line = readLine("Vbias (e.g. 53 54 55):\n> ");
        std::stringstream ss(line);
        int v = 0;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) {
        std::cerr << "[ERROR] No Vbias values.\n";
        return;
    }

    // ── 4. LET THRESHOLDS ───────────────────────────────────────────────────
    std::vector<double> fracs_pe;
    {
        std::cout << "\nLET thresholds — choose mode:\n"
                  << "  [1] Range: start end step  (e.g. 0.4 0.8 0.05)\n"
                  << "  [2] List:  space-separated  (e.g. 0.5 1.0 2.0)\n"
                  << "> " << std::flush;

        std::string line;
        if (!std::getline(std::cin, line)) throw std::runtime_error("stdin EOF");

        // Trim.
        {
            const auto b = line.find_first_not_of(" \t\r\n");
            if (b == std::string::npos) { line = ""; }
            else {
                const auto e = line.find_last_not_of(" \t\r\n");
                line = line.substr(b, e - b + 1);
            }
        }

        std::vector<double> nums;
        {
            std::stringstream ss(line);
            double v = 0;
            while (ss >> v) nums.push_back(v);
        }

        // Treat as a range when exactly 3 positive numbers are given and
        // step fits at least once inside [lo, hi].
        // FIX: old test `nums[2] < (nums[1]-nums[0])*2` was inverted.
        bool isRange = false;
        if (nums.size() == 3) {
            const double lo   = nums[0];
            const double hi   = nums[1];
            const double step = nums[2];
            if (lo > 0 && hi > lo && step > 0 && step <= (hi - lo) + 1e-9)
                isRange = true;
        }

        if (isRange) {
            const double lo   = nums[0];
            const double hi   = nums[1];
            const double step = nums[2];
            // Robust iteration: compute number of steps up front to avoid
            // floating-point accumulation drift.
            const int nSteps = static_cast<int>(std::round((hi - lo) / step));
            for (int i = 0; i <= nSteps; ++i) {
                const double v = lo + i * step;
                if (v > 0 && v <= hi + 1e-9) fracs_pe.push_back(v);
            }
            std::cout << "  Range [" << lo << ", " << hi << "] step " << step
                      << " -> " << fracs_pe.size() << " thresholds\n";
        } else {
            for (double v : nums)
                if (v > 0) fracs_pe.push_back(v);
        }

        if (fracs_pe.empty()) fracs_pe.push_back(1.0);
    }

    // ── 5. CALIBRATION LOAD ─────────────────────────────────────────────────
    std::map<int, CalibResult> calMap;
    for (int vbias : vbiasList) {
        auto cutoffs = findCalibratedCutoffsInteractive(vbias, dataDir);
        if (cutoffs.empty()) {
            std::cerr << "  [WARN] Vbias=" << vbias
                      << ": no calibration file found.\n";
            continue;
        }

        // Default: use the largest cutoff found.
        double chosen = cutoffs.back();

        // If multiple cutoffs exist, let the user pick one.
        if (cutoffs.size() > 1) {
            std::cout << "  Vbias=" << vbias << " available cutoffs: ";
            for (double co : cutoffs)
                std::cout << static_cast<int>(co) << " ";
            std::cout << "MHz\n";

            const std::string line = readLine(
                "  Which one? [" +
                std::to_string(static_cast<int>(chosen)) + "]: ");
            if (!line.empty()) {
                try {
                    const double inp = std::stod(line);
                    double best  = chosen;
                    double bestD = std::numeric_limits<double>::max();
                    for (double co : cutoffs) {
                        const double d = std::abs(co - inp);
                        if (d < bestD) { bestD = d; best = co; }
                    }
                    chosen = best;
                } catch (...) {}
            }
        }

        CalibResult cal;
        if (!loadCalibration(cal, vbias, chosen, dataDir)) continue;

        // FIX (bug #6): use GAIN_MIN/GAIN_MAX from Config.h instead of
        // hard-coded range [0, 100]. Devices with gain > 100 mV/p.e. were
        // silently discarded.
        if (cal.m < GAIN_MIN || cal.m > GAIN_MAX) {
            std::cerr << "  [WARN] Vbias=" << vbias
                      << ": gain=" << cal.m
                      << " mV/p.e. outside valid range ["
                      << GAIN_MIN << ", " << GAIN_MAX << "] — skipped.\n";
            continue;
        }
        calMap[vbias] = cal;

        std::cout << "  Vbias=" << vbias
                  << "  gain=" << cal.m << " mV/p.e."
                  << "  trig=[" << cal.t_trig_start
                  << ", " << cal.t_trig_end << "]\n";
    }
    if (calMap.empty()) {
        std::cerr << "[ERROR] No valid calibration loaded.\n";
        return;
    }

    // ── 6. FIT WINDOW ───────────────────────────────────────────────────────
    const double fit_lo = readDouble("\nFit window start [ns]: ");
    const double fit_hi = readDouble("Fit window end   [ns]: ");
    if (fit_lo >= fit_hi) {
        std::cerr << "[ERROR] fit_lo must be < fit_hi.\n";
        return;
    }

    // ── 7. OPTIONS ──────────────────────────────────────────────────────────
    const TWMethod tw_method  = askTimeWalkMethod();
    const bool     do_pe      = readYN("Per-p.e. analysis? [y/n]: ");
    const bool     use_filter = readYN("LP filter? [y/n]: ");

    // ── 8. RUN TAG ──────────────────────────────────────────────────────────
    // Format: vbias53_54_55__let0.40_0.50__tw__filt__byPE
    std::string runTag;
    {
        runTag = "vbias";
        for (int v : vbiasList) runTag += std::to_string(v) + "_";
        if (!runTag.empty() && runTag.back() == '_') runTag.pop_back();

        runTag += "__let";
        for (double f : fracs_pe) runTag += Form("%.2f_", f);
        if (!runTag.empty() && runTag.back() == '_') runTag.pop_back();

        runTag += (tw_method == TWMethod::FIT_RESIDUALS) ? "__tw" : "__notw";
        runTag += use_filter ? "__filt"  : "__nofilt";
        runTag += do_pe      ? "__byPE"  : "__noPE";
    }
    OutCtx ctx = createOutputDirs(runTag);

    // ── 9. SUMMARY HEADER ───────────────────────────────────────────────────
    std::cout << "\n+==========================================================+\n"
              << "|  Vbias: ";
    for (int v : vbiasList) std::cout << v << " ";
    std::cout << "V\n"
              << "|  Fit:   [" << fit_lo << ", " << fit_hi << "] ns\n"
              << "|  LET:   ";
    for (double f : fracs_pe) std::cout << f << " ";
    std::cout << "p.e.\n"
              << "|  TW:    "
              << (tw_method == TWMethod::FIT_RESIDUALS
                      ? "FitResidualsExp" : "None")
              << "  Filter: " << (use_filter ? "yes" : "NO")
              << "\n+==========================================================+\n\n";

    // ── 10. MAIN LOOP ────────────────────────────────────────────────────────
    std::map<int, std::map<double, std::pair<double, double>>> sigmaResults;
    const auto t0 = std::chrono::steady_clock::now();

    for (int vbias : vbiasList) {
        const auto it = calMap.find(vbias);
        if (it == calMap.end()) continue;

        const CalibResult& cal = it->second;
        const auto tVb = std::chrono::steady_clock::now();

        auto res = processOneVbias_v2(
            vbias, fracs_pe,
            cal.cutoff_MHz, cal.t_trig_start, cal.t_trig_end,
            fit_lo, fit_hi,
            tw_method, do_pe, use_filter,
            ctx, cal);

        if (!res.empty()) sigmaResults[vbias] = res;

        std::cout << "  Vbias=" << vbias << " done in "
                  << std::fixed << std::setprecision(1)
                  << std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - tVb).count()
                  << " s\n";
    }

    // ── 11. VBIAS SUMMARY ────────────────────────────────────────────────────
    drawVbiasSummary(sigmaResults, ctx);

    std::cout << "\n+==========================================================+\n"
              << "|  DONE  "
              << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - t0).count()
              << " s\n"
              << "+==========================================================+\n";

    // ── 12. INTERACTIVE END ──────────────────────────────────────────────────
    // Reopen saved canvases for interactive inspection.
    // No nested gApplication->Run() — avoid apparent hangs.
    ctx.reopenAllCanvases(40);

    } catch (const std::runtime_error& e) {
        std::cerr << "\n[ABORT] " << e.what() << "\n";
    }
}
