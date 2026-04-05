/**
 * sipm_calibrate.cpp
 * ==================
 * Programma standalone per la calibrazione p.e. tramite threshold scan.
 * Legge il PRIMO run file per un dato Vbias, esegue la calibrazione,
 * salva il risultato in  ../../data/calib_vbias<V>_cut<C>MHz.root
 * e produce il PNG di calibrazione.
 *
 * Va eseguito UNA VOLTA PER VBIAS (o quando cambia il cutoff LP).
 * Il risultato viene poi letto da sipm_tot_analysis_v5 senza rieseguire
 * il threshold scan (risparmio di ~5 s, ma soprattutto chiarezza logica).
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

    // Helper: salta righe vuote lasciate da cin >> precedenti
    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        std::cout << prompt << std::flush;
        while (line.empty()) std::getline(std::cin, line);
        return line;
    };

    // --- 1. Input parametri ------------------------------------
    std::vector<int> vbiasList;
    {
        std::string line = readLine("\nVbias values to calibrate (e.g.  53 54 55):\n> ");
        std::stringstream ss(line); int v;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) { std::cerr << "[ERROR] No Vbias specified.\n"; return; }

    double cutoff_MHz;
    std::cout << "Low-pass filter cut-off [MHz]: " << std::flush;
    std::cin >> cutoff_MHz;

    double t_trig_start, t_trig_end;
    std::cout << "Trigger window start [ns]: " << std::flush;
    std::cin >> t_trig_start;
    std::cout << "Trigger window end   [ns]: " << std::flush;
    std::cin >> t_trig_end;

    // --- 2. Loop su ogni Vbias --------------------------------
    for (int vbias : vbiasList) {
        std::cout << "\n--- Vbias = " << vbias << " V ---\n";

        // Cerca il primo run file
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
            } catch (...) {}
        }
        gSystem->FreeDirectory(dirp);

        if (foundRuns.empty()) {
            std::cerr << "[WARN] No files for Vbias=" << vbias << "\n"; continue;
        }
        const std::string& calFile = foundRuns.begin()->second;
        std::cout << "  Using: " << calFile << "\n";

        // Leggi fs_MHz
        double fs_MHz = 0;
        {
            TFile* f0 = TFile::Open(calFile.c_str(), "READ");
            if (!f0 || f0->IsZombie()) { std::cerr << "[ERROR] Cannot open.\n"; continue; }
            TTree* tr = (TTree*)f0->Get("ch1");
            if (tr) {
                const int N0=1024; Double_t tb[N0];
                tr->SetBranchAddress("time", tb); tr->GetEntry(0);
                fs_MHz = 1000.0 / (tb[1] - tb[0]);
            }
            f0->Close();
        }
        if (fs_MHz <= 0) { std::cerr << "[ERROR] Cannot read sampling rate.\n"; continue; }
        std::cout << "  Sampling rate: " << fs_MHz << " MHz\n";

        // Calibrazione
        TFile* fCal = TFile::Open(calFile.c_str(), "READ");
        if (!fCal || fCal->IsZombie()) { std::cerr << "[ERROR] Cannot open.\n"; continue; }
        TTree* tCal = (TTree*)fCal->Get("ch1");
        if (!tCal) { fCal->Close(); std::cerr << "[ERROR] No ch1 tree.\n"; continue; }

        std::string calTag = "vbias" + std::to_string(vbias);
        CalibScanData scanData;
        CalibResult cal = calibrateSpectrum(tCal, cutoff_MHz, fs_MHz,
                                             calTag, ctx, &scanData);
        fCal->Close();

        if (!cal.ok) {
            std::cerr << "[ERROR] Calibration failed for Vbias=" << vbias << "\n";
            continue;
        }

        // Store acquisition parameters in cal so they get saved to disk
        cal.cutoff_MHz   = cutoff_MHz;
        cal.t_trig_start = t_trig_start;
        cal.t_trig_end   = t_trig_end;

        // FilterDiag: estimates laser amplitude and sets cal.laser_thr
        // MUST be called before saveCalibration to persist the correct threshold
        drawFilterDiagnostics(calFile, cutoff_MHz, fs_MHz,
                              t_trig_start, t_trig_end, cal, calTag, ctx);

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
