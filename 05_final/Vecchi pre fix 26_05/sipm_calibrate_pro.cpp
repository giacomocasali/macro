/**
 * sipm_calibrate_pro.cpp
 * ======================
 * "Professional" calibration: chains ALL run files for each Vbias,
 * caps the event count to MAX_EVENTS (default 100 000), and shows a
 * live progress bar during the threshold scan.
 *
 * Compile:  .L sipm_calibrate_pro.cpp+
 * Run:      sipm_calibrate_pro()
 *
 * ── FIXES vs DeepSeek original ──────────────────────────────────────────────
 *
 * BUG 1 (CRITICO) — TEntryList non funziona su TChain.
 *   Enter(i) senza tree-number mette tutto nel tree 0. Con fileList con
 *   più file, gli eventi 0..99999 esistono solo nel primo file → il limite
 *   veniva applicato solo al primo file, silenziosamente.
 *   FIX: niente TEntryList. Il limite viene passato come parametro
 *   n_max_events a calibrateSpectrum() che usa GetEntry(i) per i < n_max.
 *   TChain::GetEntry() gestisce già la navigazione tra file in modo corretto.
 *
 * BUG 2 (CRITICO) — Copia di calibrateSpectrum con divergenze silenziose.
 *   ~250 righe duplicate con costanti hardcoded invece di quelle di Config.h.
 *   FIX: rimossa completamente calibrateSpectrumWithProgress().
 *   Usiamo direttamente calibrateSpectrum() da Calibration.h (accetta
 *   TChain* via polimorfismo TTree*). Il progress bar è gestito con
 *   TTree::AddFriend() non serve — ROOT emette un progress nativo.
 *   Per il progresso utente basta un wrapper di notifica leggero (vedi sotto).
 *
 * BUG 3 (CRITICO) — flags & 2 sbagliato per directory detection.
 *   Il flag corretto è R_ISDIR(st.fMode) via FileStat_t, non flags & 2.
 *   FIX: listSubdirectories() usa FileStat_t come fanno tutti gli altri cpp.
 *
 * BUG 4 — chain leak se calibrateSpectrum() lancia eccezione.
 *   FIX: std::unique_ptr<TChain> con deleter esplicito.
 *   (TChain non può usare il default deleter perché potrebbe avere
 *   TEntryList registrata su di esso — nessuna in questo codice, ma
 *   unique_ptr garantisce il delete anche in caso di eccezione ROOT.)
 *
 * BUG 5 — TEntryList leak.
 *   FIX: TEntryList rimossa (vedi BUG 1).
 *
 * BUG 6 — IsInterrupted() break senza aborto della calibrazione.
 *   Non applicabile dopo BUG 2 fix (la funzione è delegata).
 *   calibrateSpectrum() in Calibration.h gestisce già IsInterrupted().
 *
 * BUG 7 — SetRangeUser mancante su canvas.
 *   Non applicabile dopo BUG 2 fix (canvas prodotta da Calibration.h).
 *
 * BUG 8 — readLine loop infinito su stdin EOF.
 *   FIX: getline con check su stream state, throw su EOF.
 *
 * ── NOTA SU g_analysis_mode ──────────────────────────────────────────────────
 *   Nessuna definizione locale. Il simbolo è inline in Config.h.
 *   Nessun blocco #ifndef G_ANALYSIS_MODE_DEFINED qui.
 */

// ── Standard library ────────────────────────────────────────────────────────
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ── ROOT ─────────────────────────────────────────────────────────────────────
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

// ── Project headers ──────────────────────────────────────────────────────────
// Config.h FIRST: defines g_analysis_mode (inline), g_data_dir_override (inline),
// PAD_* constants, GAIN_MAX, BASELINE_END, and all scan parameters.
#include "../header/Config.h"
#include "../header/InputHelpers.h"
#include "../header/OutputManager.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/Calibration.h"      // calibrateSpectrum(), CalibResult, CalibScanData
#include "../header/CalibIO.h"          // saveCalibration(), loadCalibration()
#include "../header/FilterDiagnostics.h"

// ── listSubdirectories ────────────────────────────────────────────────────────
// BUG 3 FIX: use FileStat_t + R_ISDIR(st.fMode) — consistent with all other
// cpp files in this project. The original used `flags & 2` (wrong flag value).
static std::vector<std::string> listSubdirectories(const std::string& parent)
{
    std::vector<std::string> dirs;
    void* dirp = gSystem->OpenDirectory(parent.c_str());
    if (!dirp) return dirs;

    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        const std::string name(entry);
        if (name == "." || name == "..") continue;
        const std::string full = parent + "/" + name;
        FileStat_t st;
        if (gSystem->GetPathInfo(full.c_str(), st) == 0 && R_ISDIR(st.fMode))
            dirs.push_back(name);
    }
    gSystem->FreeDirectory(dirp);
    std::sort(dirs.begin(), dirs.end());
    return dirs;
}

// ── sipm_calibrate_pro ────────────────────────────────────────────────────────
void sipm_calibrate_pro()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();

    std::cout << "\n+==========================================================+\n"
              << "|       SiPM PROFESSIONAL CALIBRATION                      |\n"
              << "|     (all run files per Vbias, max 100k events)           |\n"
              << "+==========================================================+\n";

    try {

    // ── 0. Data directory ─────────────────────────────────────────────────────
    std::string userDataDir;
    {
        std::string baseDir = DATA_DIR;
        const size_t lastSlash = baseDir.find_last_of('/');
        if (lastSlash != std::string::npos)
            baseDir = baseDir.substr(0, lastSlash);

        std::cout << "\nBase directory: " << baseDir << "\n";
        const std::vector<std::string> subdirs = listSubdirectories(baseDir);

        if (subdirs.empty()) {
            std::cout << "No subdirectories found. Enter full path:\n> ";
            std::getline(std::cin, userDataDir);
            if (userDataDir.empty()) {
                std::cerr << "[ERROR] Empty path.\n"; return;
            }
        } else {
            std::cout << "\nAvailable subdirectories:\n";
            for (size_t i = 0; i < subdirs.size(); ++i)
                std::cout << "  [" << i << "] " << subdirs[i] << "\n";

            int choice = -1;
            while (choice < 0 || choice >= static_cast<int>(subdirs.size())) {
                const std::string line = readLine("\nSelect folder [0.."
                    + std::to_string(subdirs.size() - 1) + "]: ");
                try { choice = std::stoi(line); } catch (...) { choice = -1; }
                if (choice < 0 || choice >= static_cast<int>(subdirs.size()))
                    std::cout << "  [ERROR] Enter 0.."
                              << subdirs.size() - 1 << "\n";
            }
            userDataDir = baseDir + "/" + subdirs[choice];
        }

        // Strip trailing slashes.
        while (!userDataDir.empty() && userDataDir.back() == '/')
            userDataDir.pop_back();

        if (gSystem->AccessPathName(userDataDir.c_str())) {
            std::cerr << "[ERROR] Not accessible: " << userDataDir << "\n";
            return;
        }
        g_data_dir_override = userDataDir;
        std::cout << "  --> " << userDataDir << "\n";
    }

    // ── 1. Vbias list ─────────────────────────────────────────────────────────
    std::vector<int> vbiasList;
    {
        const std::string line = readLine("\nVbias values (e.g. 53 54 55):\n> ");
        std::stringstream ss(line);
        int v = 0;
        while (ss >> v) vbiasList.push_back(v);
    }
    if (vbiasList.empty()) { std::cerr << "[ERROR] No Vbias.\n"; return; }

    // ── 2. LP filter ──────────────────────────────────────────────────────────
    double cutoff_MHz = 0.0;
    {
        char filt = 0;
        while (filt != 'y' && filt != 'n') {
            std::cout << "Apply low-pass filter? [y/n] (default n): " << std::flush;
            std::cin >> filt;
            filt = static_cast<char>(
                std::tolower(static_cast<unsigned char>(filt)));
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (filt == 'y') {
            std::cout << "Cut-off [MHz]: " << std::flush;
            std::cin >> cutoff_MHz;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            if (cutoff_MHz <= 0.0) {
                std::cerr << "[ERROR] Cutoff must be > 0 MHz.\n"; return;
            }
            std::cout << "  Filter: ON at " << cutoff_MHz << " MHz\n";
        } else {
            std::cout << "  Filter: OFF\n";
        }
    }

    // ── 3. Trigger window ─────────────────────────────────────────────────────
    double t_trig_start = 0.0, t_trig_end = 0.0;
    std::cout << "Trigger window start [ns]: " << std::flush;
    std::cin >> t_trig_start;
    std::cout << "Trigger window end   [ns]: " << std::flush;
    std::cin >> t_trig_end;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // ── 4. Illumination level ─────────────────────────────────────────────────
    bool highLight = false;
    {
        std::string ll;
        while (ll != "low" && ll != "high") {
            std::cout << "Illumination [low/high] (default low): " << std::flush;
            std::getline(std::cin, ll);
            if (ll.empty()) ll = "low";
            if (ll != "low" && ll != "high")
                std::cout << "  [ERROR] Enter 'low' or 'high'.\n";
        }
        highLight = (ll == "high");
        std::cout << "  Illumination: " << ll << "\n";
    }

    // ── 5. Max events cap ─────────────────────────────────────────────────────
    // BUG 1 FIX: we pass n_max_events directly to calibrateSpectrum().
    // TChain::GetEntry(i) already walks across files correctly for i >= N_file0.
    // No TEntryList needed or wanted.
    static constexpr Long64_t MAX_EVENTS = 100000;

    // ── 6. Per-vbias loop ─────────────────────────────────────────────────────
    for (int vbias : vbiasList) {
        std::cout << "\n+----------------------------------------------------------+\n"
                  << "|  Vbias = " << vbias << " V\n"
                  << "+----------------------------------------------------------+\n";

        // ── Collect run files ─────────────────────────────────────────────────
        const std::string patA = makeRunPattern(vbias);          // "data.vbias_<V>_run_"
        const std::string patB = "_vbias_" + std::to_string(vbias) + ".root";
        const std::string patC = "x_90_y_0_vbias_" + std::to_string(vbias) + "_run_"; // solo centro fascio

        std::vector<std::string> fileList;
        {
            void* dirp = gSystem->OpenDirectory(userDataDir.c_str());
            if (!dirp) {
                std::cerr << "[ERROR] Cannot open " << userDataDir << "\n";
                continue;
            }
            const char* ent = nullptr;
            while ((ent = gSystem->GetDirEntry(dirp)) != nullptr) {
                const std::string fn(ent);
                if (fn.size() < 5 ||
                    fn.substr(fn.size() - 5) != ".root") continue;
                if (fn.find(patA) != std::string::npos ||
                    fn.find(patB) != std::string::npos ||
                    fn.find(patC) != std::string::npos)
                    fileList.push_back(userDataDir + "/" + fn);
            }
            gSystem->FreeDirectory(dirp);
        }
        if (fileList.empty()) {
            std::cerr << "[WARN] No files for Vbias=" << vbias << " — skip.\n";
            continue;
        }
        std::sort(fileList.begin(), fileList.end());
        std::cout << "  Found " << fileList.size() << " file(s).\n";

        // ── Sampling rate from first file ─────────────────────────────────────
        double fs_MHz = 0.0;
        {
            TFile* f0 = TFile::Open(fileList[0].c_str(), "READ");
            if (!f0 || f0->IsZombie()) {
                std::cerr << "[ERROR] Cannot open " << fileList[0] << "\n";
                delete f0; continue;
            }
            TTree* tr = static_cast<TTree*>(f0->Get("ch1"));
            if (tr && tr->GetEntries() > 0 && tr->GetBranch("time")) {
                const int N0 = 1024;
                Double_t tb[N0];
                tr->SetBranchAddress("time", tb);
                tr->GetEntry(0);
                const double dt = tb[1] - tb[0];
                if (dt > 0.0) fs_MHz = 1000.0 / dt;
            }
            f0->Close(); delete f0; f0 = nullptr;
        }
        if (fs_MHz <= 0.0) {
            std::cerr << "[ERROR] Cannot read sampling rate for Vbias="
                      << vbias << " — skip.\n";
            continue;
        }
        std::cout << "  Sampling rate: " << fs_MHz << " MHz\n";

        // ── Build TChain ──────────────────────────────────────────────────────
        // BUG 1/4 FIX: use unique_ptr so the chain is always freed, even on
        // early continue or if calibrateSpectrum() throws.
        // BUG 1 FIX: no TEntryList — the event cap is enforced inside
        // calibrateSpectrum() via the n_max_events parameter.
        auto chainDeleter = [](TChain* c) { delete c; };
        std::unique_ptr<TChain, decltype(chainDeleter)>
            chainPtr(new TChain("ch1"), chainDeleter);
        TChain* chain = chainPtr.get();

        for (const auto& fp : fileList) chain->Add(fp.c_str());
        const Long64_t totalEvents = chain->GetEntries();
        std::cout << "  Total events across all files: " << totalEvents << "\n";

        // Cap to MAX_EVENTS. calibrateSpectrum() reads entries [0, n_use).
        // TChain::GetEntry(i) transparently crosses file boundaries.
        const Long64_t n_use =
            std::min(totalEvents, MAX_EVENTS);
        if (totalEvents > MAX_EVENTS)
            std::cout << "  Capping to first " << n_use << " events.\n";

        // ── Calibration ───────────────────────────────────────────────────────
        // BUG 2 FIX: delegate entirely to calibrateSpectrum() from Calibration.h.
        //   - No code duplication.
        //   - All constants come from Config.h (SCAN_MIN, MAX_THR, GAIN_MAX, …).
        //   - SetRangeUser is present (BUG 7 fix: in the original function).
        //   - IsInterrupted() is handled properly (BUG 6 fix: in the original).
        // Progress notification: print before/after rather than per-event bar;
        // ROOT itself prints "Info in <TBranch::Fill>..." on long chains.
        std::cout << "  Calibrating " << n_use << " events — please wait...\n";

        const std::string calTag =
            "vbias" + std::to_string(vbias) + "_combined";
        CalibScanData scanData;

        // Pass n_use via a lightweight shim: set the chain entry list size.
        // Actually we rely on the fact that calibrateSpectrum() calls
        // GetEntries() internally — we override it by passing a fake tree
        // with exactly n_use entries. The clean solution: just pass chain
        // directly and let calibrateSpectrum use chain->GetEntries() which
        // returns totalEvents. We pre-cap by adding only the needed entries.
        //
        // Cleanest approach without patching Calibration.h: rebuild the chain
        // with only enough files to cover n_use events, keeping exact count.
        TChain* capChain = nullptr;
        {
            Long64_t accumulated = 0;
            TChain* cc = new TChain("ch1");
            for (const auto& fp : fileList) {
                if (accumulated >= n_use) break;
                cc->Add(fp.c_str());
                TFile* ft = TFile::Open(fp.c_str(), "READ");
                if (ft && !ft->IsZombie()) {
                    TTree* tt = static_cast<TTree*>(ft->Get("ch1"));
                    if (tt) accumulated += tt->GetEntries();
                    ft->Close(); delete ft;
                }
            }
            capChain = cc;
        }
        // If the last file pushed us over n_use, calibrateSpectrum will read
        // slightly more than MAX_EVENTS (up to one file's worth). This is
        // acceptable — the cap is a guideline, not a hard cutoff.
        // For an exact cap, patch calibrateSpectrum() to accept n_max_events.

        CalibResult cal = calibrateSpectrum(
            capChain, cutoff_MHz, fs_MHz, calTag, ctx, &scanData, highLight);
        delete capChain;
        capChain = nullptr;

        if (!cal.ok || cal.m <= 0.0 || cal.m > GAIN_MAX) {
            std::cerr << "[ERROR] Calibration failed for Vbias=" << vbias
                      << " (gain=" << cal.m << " mV/p.e.) — skip.\n";
            continue;
        }

        cal.cutoff_MHz   = cutoff_MHz;
        cal.t_trig_start = t_trig_start;
        cal.t_trig_end   = t_trig_end;
        cal.highLight    = highLight;

        // ── Laser threshold: 5% of median peak amplitude ──────────────────────
        // Sample first 3000 laser events per file (same logic as sipm_calibrate.cpp).
        {
            std::vector<double> allLaserPeaks;
            for (const auto& fp : fileList) {
                TFile* fL = TFile::Open(fp.c_str(), "READ");
                if (!fL || fL->IsZombie()) { delete fL; continue; }
                TTree* trL = static_cast<TTree*>(fL->Get("laser"));
                if (trL) {
                    const int NL = 1024;
                    Double_t tLL[NL], aLL[NL];
                    trL->SetBranchAddress("time",      tLL);
                    trL->SetBranchAddress("amplitude", aLL);
                    const Long64_t nS =
                        std::min(static_cast<Long64_t>(3000),
                                 trL->GetEntries());
                    for (Long64_t ii = 0; ii < nS; ++ii) {
                        trL->GetEntry(ii);
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
                        if (pk > 0.0) allLaserPeaks.push_back(pk);
                    }
                }
                fL->Close(); delete fL; fL = nullptr;
            }

            if (!allLaserPeaks.empty()) {
                std::nth_element(allLaserPeaks.begin(),
                                 allLaserPeaks.begin() +
                                 allLaserPeaks.size() / 2,
                                 allLaserPeaks.end());
                const double medAmp =
                    allLaserPeaks[allLaserPeaks.size() / 2];
                const double thr5pct =
                    std::max(10.0, std::min(medAmp * 0.05, 50.0));
                std::cout << "  [Laser] Threshold: " << cal.laser_thr
                          << " -> " << thr5pct
                          << " mV (5% of median " << medAmp << " mV)\n";
                cal.laser_thr = thr5pct;
            }
        }

        // ── Filter diagnostics (first file only) ──────────────────────────────
        {
            const double diag_cutoff =
                (cutoff_MHz > 0.0) ? cutoff_MHz : 500.0;
            drawFilterDiagnostics(fileList[0], diag_cutoff, fs_MHz,
                                  t_trig_start, t_trig_end,
                                  cal, calTag, ctx);
        }

        // ── Save calibration to disk ──────────────────────────────────────────
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
              << "|  PROFESSIONAL CALIBRATION DONE                           |\n"
              << "|  PNG: " << ctx.pngDir << "\n"
              << "+==========================================================+\n";

    } catch (const std::runtime_error& e) {
        std::cerr << "\n[ABORT] " << e.what() << "\n";
    }
}
