/**
 * sipm_calibrate.cpp
 * ==================
 * Standalone program for p.e. calibration via threshold scan.
 * Reads the FIRST run file for a given Vbias, runs the calibration,
 * saves the result to  ../../data/calib_vbias<V>_cut<C>MHz.root
 * and produces the calibration PNG.
 *
 * Run ONCE PER VBIAS (or whenever the LP cutoff changes).
 * The result is then read by sipm_tot_analysis_v5 without re-running
 * the threshold scan (saves ~5 s and keeps the logic clean).
 *
 * Compile: .L sipm_calibrate.cpp+
 * Run:     sipm_calibrate()
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/FilterDiagnostics.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <sstream>
#include <limits>
#include <TApplication.h>

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>

void sipm_calibrate()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();

    std::cout << "\n+==========================================================+\n"
              << "|        SiPM CALIBRATION -- threshold scan                 |\n"
              << "+==========================================================+\n";

    // Helper: skip blank lines left by previous cin >> calls
    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        std::cout << prompt << std::flush;
        while (line.empty()) std::getline(std::cin, line);
        return line;
    };

    // --- 1. Input parameters ------------------------------------
    std::vector<int> vbiasList;
    {
        std::string line = readLine("\nVbias values to calibrate (e.g.  53 54 55):\n> ");
        std::stringstream ss(line); int v;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) { std::cerr << "[ERROR] No Vbias specified.\n"; return; }

    // Low-pass filter: default OFF (raw data)
    double cutoff_MHz = 0.0;  // <= 0 means no filter
    {
        char filt = 0;
        while (filt != 'y' && filt != 'n') {
            std::cout << "Apply low-pass filter? [y/n] (default: n): " << std::flush;
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

    double t_trig_start, t_trig_end;
    std::cout << "Trigger window start [ns]: " << std::flush;
    std::cin >> t_trig_start;
    std::cout << "Trigger window end   [ns]: " << std::flush;
    std::cin >> t_trig_end;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // --- 2. Loop over each Vbias --------------------------------
    for (int vbias : vbiasList) {
        std::cout << "\n--- Vbias = " << vbias << " V ---\n";

        // Find the first run file
        const std::string dataDir = DATA_DIR;
        const std::string pattern = "data.vbias_" + std::to_string(vbias) + "_run_";
        std::map<int, std::string> foundRuns;
        void* dirp = gSystem->OpenDirectory(dataDir.c_str());
        if (!dirp) { std::cerr << "[ERROR] Cannot open " << dataDir << "\n"; continue; }
        const char* entry;
        while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
            std::string fname(entry);
            if (fname.find(pattern) == std::string::npos) continue;
            if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
            size_t pos = fname.find("_run_");
            if (pos == std::string::npos) continue;
            try {
                std::string sub = fname.substr(pos + 5);
                size_t dot = sub.find(".root");
                if (dot != std::string::npos) sub = sub.substr(0, dot);
                foundRuns[std::stoi(sub)] = dataDir + "/" + fname;
            } catch (...) {
                std::cerr << "[WARN] Cannot parse run number from: " << fname << "\n";
            }
        }
        gSystem->FreeDirectory(dirp);

        if (foundRuns.empty()) {
            std::cerr << "[WARN] No files for Vbias=" << vbias << "\n"; continue;
        }
        const std::string& calFile = foundRuns.begin()->second;
        std::cout << "  Using: " << calFile << "\n";

        // Read sampling rate
        double fs_MHz = 0;
        {
            TFile* f0 = TFile::Open(calFile.c_str(), "READ");
            if (!f0 || f0->IsZombie()) { std::cerr << "[ERROR] Cannot open.\n"; continue; }
            TTree* tr = (TTree*)f0->Get("ch1");
            if (tr && tr->GetEntries() > 0 && tr->GetBranch("time")) {
                const int N0=1024; Double_t tb[N0];
                tr->SetBranchAddress("time", tb); tr->GetEntry(0);
                double dt = tb[1] - tb[0];
                if (dt > 0) fs_MHz = 1000.0 / dt;
                else std::cerr << "[ERROR] Invalid time spacing in waveform.\n";
            }
            f0->Close();
            delete f0; f0 = nullptr;
        }
        if (fs_MHz <= 0) { std::cerr << "[ERROR] Cannot read sampling rate.\n"; continue; }
        std::cout << "  Sampling rate: " << fs_MHz << " MHz\n";

        // Run calibration
        TFile* fCal = TFile::Open(calFile.c_str(), "READ");
        if (!fCal || fCal->IsZombie()) { std::cerr << "[ERROR] Cannot open.\n"; continue; }
        TTree* tCal = (TTree*)fCal->Get("ch1");
        if (!tCal) { fCal->Close(); std::cerr << "[ERROR] No ch1 tree.\n"; continue; }

        std::string calTag = "vbias" + std::to_string(vbias);
        CalibScanData scanData;
        CalibResult cal = calibrateSpectrum(tCal, cutoff_MHz, fs_MHz,
                                             calTag, ctx, &scanData);
        fCal->Close(); delete fCal; fCal = nullptr;

        if (!cal.ok) {
            std::cerr << "[ERROR] Calibration failed for Vbias=" << vbias << "\n";
            continue;
        }
        if (cal.m <= 0 || cal.m > 100.0) {
            std::cerr << "[ERROR] Suspicious gain=" << cal.m
                      << " mV/p.e. for Vbias=" << vbias
                      << " — calibration may have failed.\n";
            continue;
        }

        // Store acquisition parameters in cal so they get saved to disk
        cal.cutoff_MHz   = cutoff_MHz;
        cal.t_trig_start = t_trig_start;
        cal.t_trig_end   = t_trig_end;

        // FilterDiag: runs diagnostics and estimates laser amplitude.
        // It sets cal.laser_thr to 20% of peak — we override to 5% below.
        double diag_cutoff = (cutoff_MHz > 0) ? cutoff_MHz : 500.0;
        drawFilterDiagnostics(calFile, diag_cutoff, fs_MHz,
                              t_trig_start, t_trig_end, cal, calTag, ctx);

        // Override laser_thr: 5% of median peak (not 20%).
        // 20% is too high for precise timing — the laser rise time causes
        // a ~15-20 ns shift in t_laser, making delta_t negative.
        {
            TFile* fLas = TFile::Open(calFile.c_str(), "READ");
            if (fLas && !fLas->IsZombie()) {
                TTree* trL = (TTree*)fLas->Get("laser");
                if (trL) {
                    const int NL = 1024;
                    Double_t tLL[NL], aLL[NL];
                    trL->SetBranchAddress("time", tLL);
                    trL->SetBranchAddress("amplitude", aLL);
                    std::vector<double> peaks;
                    Long64_t nS = std::min((Long64_t)500, trL->GetEntries());
                    for (Long64_t ii = 0; ii < nS; ++ii) {
                        trL->GetEntry(ii);
                        std::vector<double> pre;
                        for (int j = 0; j < NL; ++j)
                            if (tLL[j] < BASELINE_END) pre.push_back(aLL[j]);
                        double off = 0;
                        if (!pre.empty()) {
                            auto tmp = pre;
                            std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end());
                            off = tmp[tmp.size()/2];
                        }
                        double pk = -1e9;
                        for (int j = 0; j < NL; ++j)
                            if (aLL[j] - off > pk) pk = aLL[j] - off;
                        if (pk > 0) peaks.push_back(pk);
                    }
                    if (!peaks.empty()) {
                        std::nth_element(peaks.begin(), peaks.begin()+peaks.size()/2, peaks.end());
                        double medAmp = peaks[peaks.size()/2];
                        double thr5pct = std::max(10.0, std::min(medAmp * 0.05, 50.0));
                        std::cout << "  [Laser] Overriding threshold: "
                                  << cal.laser_thr << " → " << thr5pct
                                  << " mV (5% of " << medAmp << " mV)\n";
                        cal.laser_thr = thr5pct;
                    }
                }
                fLas->Close(); delete fLas;
            }
        }

        // Save everything (gain, offset, laser_thr, cutoff, trig window)
        saveCalibration(cal, scanData.thresholds, scanData.counts,
                        vbias, cutoff_MHz, dataDir);

        std::cout << "  Gain       = " << cal.m << " mV/p.e.\n"
                  << "  Offset     = " << cal.q << " mV\n"
                  << "  Laser thr  = " << cal.laser_thr << " mV\n"
                  << "  Cutoff     = " << cal.cutoff_MHz << " MHz\n"
                  << "  Trig win   = [" << cal.t_trig_start
                  << ", " << cal.t_trig_end << "] ns\n";
    }

    std::cout << "\n+==========================================================+\n"
              << "|  CALIBRATION DONE\n"
              << "|  PNG: " << ctx.pngDir << "\n"
              << "+==========================================================+\n";
}
