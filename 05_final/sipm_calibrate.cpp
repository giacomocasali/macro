/**
 * sipm_calibrate.cpp
 * ==================
 * Standalone program for p.e. calibration via threshold scan.
 * Reads the FIRST run file for a given Vbias, runs the calibration,
 * saves the result to  <dataDir>/calib_vbias<V>_cut<C>MHz.root
 * and produces the calibration PNG.
 *
 * Run ONCE PER VBIAS (or whenever the LP cutoff changes).
 * The result is then read by sipm_tot_analysis without re-running
 * the threshold scan (saves ~5 s and keeps the logic clean).
 *
 * Compile: .L sipm_calibrate.cpp+
 * Run:     sipm_calibrate()
 *
 * FIX (include order): STL and ROOT headers now come before project
 *   headers.
 * FIX (dead include): removed <TApplication.h> — TApplication is
 *   never used here.
 * FIX (g_data_dir_override): g_data_dir_override is now set after
 *   the user selects the data directory, so any downstream call to
 *   effectiveDataDir() returns the correct path rather than the
 *   hardcoded DATA_DIR_LOW.
 * FIX (dataDir shadow): removed the per-loop `const std::string
 *   dataDir = userDataDir;` re-declaration; the outer `userDataDir`
 *   variable is used directly everywhere instead.
 * FIX (g_analysis_mode ODR): no definition here; the symbol lives
 *   in Config.h as `inline int g_analysis_mode = 0;`.
 */

// ── Standard library ────────────────────────────────────────────────────────
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ── ROOT ─────────────────────────────────────────────────────────────────────
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

// ── Project headers ──────────────────────────────────────────────────────────
#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/FilterDiagnostics.h"

void sipm_calibrate()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();

    std::cout << "\n+==========================================================+\n"
              << "|        SiPM CALIBRATION -- threshold scan                 |\n"
              << "+==========================================================+\n";

    // Helper: skip blank lines left by previous std::cin >> calls.
    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        std::cout << prompt << std::flush;
        while (line.empty()) std::getline(std::cin, line);
        return line;
    };

    // ── 0. Data directory ────────────────────────────────────────────────────
    std::string userDataDir;
    {
        const std::string line = readLine(
            "\nData directory (leave blank to use default: "
            + DATA_DIR + "):\n> ");
        if (!line.empty() && line != " ") {
            userDataDir = line;
            // Strip trailing slashes.
            while (!userDataDir.empty() && userDataDir.back() == '/')
                userDataDir.pop_back();
            std::cout << "  Using directory: " << userDataDir << "\n";
        } else {
            userDataDir = DATA_DIR;
            std::cout << "  Using default directory: " << userDataDir << "\n";
        }
    }
    // Publish the choice so every downstream call to effectiveDataDir()
    // returns userDataDir instead of the hardcoded DATA_DIR_LOW.
    g_data_dir_override = userDataDir;

    // ── 1. Input parameters ──────────────────────────────────────────────────
    std::vector<int> vbiasList;
    {
        const std::string line = readLine(
            "\nVbias values to calibrate (e.g. 53 54 55):\n> ");
        std::stringstream ss(line);
        int v = 0;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) {
        std::cerr << "[ERROR] No Vbias specified.\n";
        return;
    }

    // Low-pass filter: default OFF (raw data).
    double cutoff_MHz = 0.0;   // <= 0 means no filter
    {
        char filt = 0;
        while (filt != 'y' && filt != 'n') {
            std::cout << "Apply low-pass filter? [y/n] (default: n): "
                      << std::flush;
            std::cin >> filt;
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (filt == 'y') {
            std::cout << "Low-pass filter cut-off [MHz]: " << std::flush;
            std::cin >> cutoff_MHz;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            if (cutoff_MHz <= 0) {
                std::cerr << "[ERROR] Cut-off must be > 0 MHz.\n";
                return;
            }
            std::cout << "  Filter: ON at " << cutoff_MHz << " MHz\n";
        } else {
            std::cout << "  Filter: OFF (using raw data)\n";
        }
    }

    double t_trig_start = 0.0;
    double t_trig_end   = 0.0;
    std::cout << "Trigger window start [ns]: " << std::flush;
    std::cin >> t_trig_start;
    std::cout << "Trigger window end   [ns]: " << std::flush;
    std::cin >> t_trig_end;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Illumination level.
    bool highLight = false;
    {
        std::string ll;
        while (ll != "low" && ll != "high") {
            std::cout << "Illumination level [low/high]"
                         " (low = few p.e./event, high = many p.e./event): "
                      << std::flush;
            std::getline(std::cin, ll);
            if (ll.empty()) ll = "low";
            if (ll != "low" && ll != "high")
                std::cout << "  [ERROR] Enter 'low' or 'high'.\n";
        }
        highLight = (ll == "high");
        std::cout << "  Illumination: "
                  << (highLight ? "HIGH (median spacing algorithm)\n"
                                : "LOW  (dominant peak algorithm)\n");
    }

    // ── 2. Per-vbias calibration loop ────────────────────────────────────────
    for (int vbias : vbiasList) {
        std::cout << "\n--- Vbias = " << vbias << " V ---\n";

        // Find the first run file — supports TWO naming conventions:
        //   (A) data.vbias_{v}_run_{N}.root          (original)
        //   (B) data_x_{x}_y_{y}_vbias_{v}.root      (2D scan serpentina)
        const std::string patternA = makeRunPattern(vbias);
        const std::string patternB = "_vbias_" + std::to_string(vbias) + ".root";

        std::map<int, std::string> foundRuns;
        void* dirp = gSystem->OpenDirectory(userDataDir.c_str());
        if (!dirp) {
            std::cerr << "[ERROR] Cannot open " << userDataDir << "\n";
            continue;
        }

        int syntheticKey = 0;
        const char* entry = nullptr;
        while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
            std::string fname(entry);
            if (fname.size() < 5 ||
                fname.substr(fname.size() - 5) != ".root") continue;

            // Convention A: data.vbias_<V>_run_<N>.root
            if (fname.find(patternA) != std::string::npos) {
                const size_t pos = fname.find("_run_");
                if (pos == std::string::npos) continue;
                try {
                    std::string sub = fname.substr(pos + 5);
                    const size_t dot = sub.find(".root");
                    if (dot != std::string::npos) sub = sub.substr(0, dot);
                    foundRuns[std::stoi(sub)] = userDataDir + "/" + fname;
                } catch (...) {
                    std::cerr << "[WARN] Cannot parse run number from: "
                              << fname << "\n";
                }
                continue;
            }

            // Convention B: data_x_<X>_y_<Y>_vbias_<V>.root
            // Prefer the canonical x=90 y=0 file if it exists.
            if (fname.find(patternB) != std::string::npos) {
                const std::string forced =
                    userDataDir + "/data_x_90_y_0_vbias_"
                    + std::to_string(vbias) + ".root";
                if (!gSystem->AccessPathName(forced.c_str())) {
                    foundRuns[0] = forced;
                } else {
                    foundRuns[syntheticKey++] = userDataDir + "/" + fname;
                }
            }
        }
        gSystem->FreeDirectory(dirp);

        if (foundRuns.empty()) {
            std::cerr << "[WARN] No files for Vbias=" << vbias << "\n";
            continue;
        }

        const std::string& calFile = foundRuns.begin()->second;
        std::cout << "  Using: " << calFile << "\n";

        // ── Read sampling rate from the first waveform ────────────────────
        double fs_MHz = 0.0;
        {
            TFile* f0 = TFile::Open(calFile.c_str(), "READ");
            if (!f0 || f0->IsZombie()) {
                std::cerr << "[ERROR] Cannot open file.\n";
                delete f0;
                continue;
            }
            TTree* tr = static_cast<TTree*>(f0->Get("ch1"));
            if (tr && tr->GetEntries() > 0 && tr->GetBranch("time")) {
                const int N0 = 1024;
                Double_t tb[N0];
                tr->SetBranchAddress("time", tb);
                tr->GetEntry(0);
                const double dt = tb[1] - tb[0];
                if (dt > 0)
                    fs_MHz = 1000.0 / dt;
                else
                    std::cerr << "[ERROR] Invalid time spacing in waveform.\n";
            }
            f0->Close();
            delete f0;
            f0 = nullptr;
        }
        if (fs_MHz <= 0) {
            std::cerr << "[ERROR] Cannot read sampling rate.\n";
            continue;
        }
        std::cout << "  Sampling rate: " << fs_MHz << " MHz\n";

        // ── Run calibration ───────────────────────────────────────────────
        TFile* fCal = TFile::Open(calFile.c_str(), "READ");
        if (!fCal || fCal->IsZombie()) {
            std::cerr << "[ERROR] Cannot open file for calibration.\n";
            delete fCal;
            continue;
        }
        TTree* tCal = static_cast<TTree*>(fCal->Get("ch1"));
        if (!tCal) {
            fCal->Close();
            delete fCal;
            std::cerr << "[ERROR] No ch1 tree in " << calFile << "\n";
            continue;
        }

        const std::string calTag = "vbias" + std::to_string(vbias);
        CalibScanData scanData;
        CalibResult cal = calibrateSpectrum(
            tCal, cutoff_MHz, fs_MHz, calTag, ctx, &scanData, highLight);
        fCal->Close();
        delete fCal;
        fCal = nullptr;

        if (!cal.ok) {
            std::cerr << "[ERROR] Calibration failed for Vbias="
                      << vbias << "\n";
            continue;
        }
        if (cal.m <= 0 || cal.m > GAIN_MAX) {
            std::cerr << "[ERROR] Suspicious gain=" << cal.m
                      << " mV/p.e. for Vbias=" << vbias
                      << " — calibration may have failed"
                         " (valid range: (0, " << GAIN_MAX << "] mV/p.e.).\n";
            continue;
        }

        // Store acquisition parameters so they get saved to disk.
        cal.cutoff_MHz   = cutoff_MHz;
        cal.t_trig_start = t_trig_start;
        cal.t_trig_end   = t_trig_end;
        cal.highLight    = highLight;

        // ── FilterDiag: diagnostics + initial laser_thr estimate (20%) ───
        const double diag_cutoff =
            (cutoff_MHz > 0) ? cutoff_MHz : 500.0;
        drawFilterDiagnostics(calFile, diag_cutoff, fs_MHz,
                              t_trig_start, t_trig_end,
                              cal, calTag, ctx);

        // ── Override laser_thr to 5% of median peak ───────────────────────
        // drawFilterDiagnostics sets cal.laser_thr to 20% of peak.
        // 20% is too high for precise timing — the laser rise time causes
        // a ~15-20 ns shift in t_laser, making delta_t negative.
        {
            TFile* fLas = TFile::Open(calFile.c_str(), "READ");
            if (fLas && !fLas->IsZombie()) {
                TTree* trL = static_cast<TTree*>(fLas->Get("laser"));
                if (trL) {
                    const int NL = 1024;
                    Double_t tLL[NL], aLL[NL];
                    trL->SetBranchAddress("time",      tLL);
                    trL->SetBranchAddress("amplitude", aLL);

                    std::vector<double> peaks;
                    // FIX (bug #13): 3000 samples for a stable median when
                    // the amplitude distribution has an afterpulse tail.
                    const Long64_t nS =
                        std::min(static_cast<Long64_t>(3000),
                                 trL->GetEntries());
                    peaks.reserve(static_cast<size_t>(nS));

                    for (Long64_t ii = 0; ii < nS; ++ii) {
                        trL->GetEntry(ii);
                        // Baseline: median of pre-signal samples.
                        std::vector<double> pre;
                        for (int j = 0; j < NL; ++j)
                            if (tLL[j] < BASELINE_END) pre.push_back(aLL[j]);
                        double off = 0.0;
                        if (!pre.empty()) {
                            auto tmp = pre;
                            std::nth_element(tmp.begin(),
                                             tmp.begin() + tmp.size() / 2,
                                             tmp.end());
                            off = tmp[tmp.size() / 2];
                        }
                        double pk = -1e9;
                        for (int j = 0; j < NL; ++j)
                            if (aLL[j] - off > pk) pk = aLL[j] - off;
                        if (pk > 0) peaks.push_back(pk);
                    }

                    if (!peaks.empty()) {
                        std::nth_element(peaks.begin(),
                                         peaks.begin() + peaks.size() / 2,
                                         peaks.end());
                        const double medAmp = peaks[peaks.size() / 2];
                        const double thr5pct =
                            std::max(10.0, std::min(medAmp * 0.05, 50.0));
                        std::cout << "  [Laser] Overriding threshold: "
                                  << cal.laser_thr << " → " << thr5pct
                                  << " mV (5% of " << medAmp << " mV)\n";
                        cal.laser_thr = thr5pct;
                    }
                }
                fLas->Close();
                delete fLas;
            }
        }

        // ── Save calibration to disk ──────────────────────────────────────
        saveCalibration(cal, scanData.thresholds, scanData.counts,
                        vbias, cutoff_MHz, userDataDir);

        std::cout << "  Gain       = " << cal.m         << " mV/p.e.\n"
                  << "  Offset     = " << cal.q         << " mV\n"
                  << "  Laser thr  = " << cal.laser_thr << " mV\n"
                  << "  Cutoff     = " << cal.cutoff_MHz << " MHz\n"
                  << "  Trig win   = [" << cal.t_trig_start
                  << ", " << cal.t_trig_end << "] ns\n"
                  << "  Light mode = "
                  << (cal.highLight ? "HIGH" : "LOW") << "\n";
    }

    std::cout << "\n+==========================================================+\n"
              << "|  CALIBRATION DONE\n"
              << "|  PNG: " << ctx.pngDir << "\n"
              << "+==========================================================+\n";
}
