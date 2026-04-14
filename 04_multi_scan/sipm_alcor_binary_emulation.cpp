/**
 * sipm_alcor_binary_emulation.cpp
 *
 * Compile: .L sipm_alcor_binary_emulation.cpp+
 * Run:     sipm_alcor_binary_emulation()
 *
 * Purpose:
 *   Emulate an ALCOR-like binary readout:
 *   - no laser coincidence
 *   - no absolute signal arrival time
 *   - per event output is only hit/no-hit above threshold
 *
 * It keeps the same calibration chain (gain/offset from CalibIO) so thresholds
 * can still be expressed in p.e., but timing-based TOT analysis is not used.
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/SignalProcessing.h"
#include "../header/ProgressBar.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

namespace {

struct BinaryStats {
    long nEvents = 0;
    long nHit = 0;
    long nDirtyBaseline = 0;
    double eff = 0.0;
    double effErr = 0.0;
    std::string ampTmpRoot;
};

static std::vector<double> findCalibratedCutoffs(int vbias, const std::string& dataDir) {
    std::vector<double> cutoffs;
    const std::string prefix = "calib_vbias" + std::to_string(vbias) + "_cut";
    const std::string suffix = "mhz.root";

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

static std::map<int, std::string> findRunsForVbias(int vbias, const std::string& dataDir) {
    std::map<int, std::string> runs;
    const std::string pattern = "data.vbias_" + std::to_string(vbias) + "_run_";
    void* dirp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp) return runs;

    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        std::string fname(entry);
        if (fname.find(pattern) == std::string::npos) continue;
        if (fname.size() < 5 || fname.substr(fname.size() - 5) != ".root") continue;
        size_t pos = fname.find("_run_");
        if (pos == std::string::npos) continue;
        try {
            std::string num = fname.substr(pos + 5);
            size_t dot = num.find(".root");
            if (dot != std::string::npos) num = num.substr(0, dot);
            int runId = std::stoi(num);
            runs[runId] = dataDir + "/" + fname;
        } catch (...) {}
    }
    gSystem->FreeDirectory(dirp);
    return runs;
}

static double estimateFsMHzFromFirstRun(const std::string& path) {
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        delete f;
        return 0.0;
    }
    TTree* tr = (TTree*)f->Get("ch1");
    if (!tr || !tr->GetBranch("time") || tr->GetEntries() < 1) {
        f->Close();
        delete f;
        return 0.0;
    }
    const int N = 1024;
    Double_t t[N] = {};
    tr->SetBranchAddress("time", t);
    tr->GetEntry(0);
    double dt = t[1] - t[0];
    f->Close();
    delete f;
    if (!std::isfinite(dt) || dt <= 0) return 0.0;
    return 1000.0 / dt;
}

static BinaryStats analyseBinaryOneLET(const std::map<int, std::string>& runs,
                                       int jStart, int jEnd,
                                       double letThr, bool useFilter,
                                       double cutoffMHz, double fsMHz,
                                       const std::string& tag,
                                       OutCtx& ctx) {
    BinaryStats st;
    const int N = 1024;
    const double baselineMaxRMS = useFilter ? 2.0 : 4.0;

    TH1D* hAmp = new TH1D(Form("hAmp_%s", tag.c_str()),
                          Form("A_{max} distribution %s;A_{max} (mV);Events", tag.c_str()),
                          500, -5, 200);
    hAmp->SetDirectory(nullptr);

    long totalToProcess = 0;
    for (const auto& kv : runs) {
        TFile* f = TFile::Open(kv.second.c_str(), "READ");
        if (!f || f->IsZombie()) { delete f; continue; }
        TTree* tr = (TTree*)f->Get("ch1");
        if (tr) totalToProcess += tr->GetEntries();
        f->Close();
        delete f;
    }

    ProgressBar bar(totalToProcess, "ALCOR binary " + tag + (useFilter ? "" : " [NO FILT]"));
    long globalIdx = 0;

    for (const auto& kv : runs) {
        TFile* f = TFile::Open(kv.second.c_str(), "READ");
        if (!f || f->IsZombie()) { delete f; continue; }

        TTree* tr = (TTree*)f->Get("ch1");
        if (!tr) { f->Close(); delete f; continue; }

        Double_t t[N], a[N];
        tr->SetBranchAddress("time", t);
        tr->SetBranchAddress("amplitude", a);
        tr->SetCacheSize(2 * 1024 * 1024);

        Long64_t nFile = tr->GetEntries();
        for (Long64_t i = 0; i < nFile; ++i) {
            if (bar.update(globalIdx, 0, 0, st.nHit)) goto done;
            ++globalIdx;

            tr->GetEntry(i);
            std::vector<double> vt(t, t + N), va(a, a + N);
            bool blOk = true;
            std::vector<double> af = correctBaseline(vt, va, BASELINE_START, BASELINE_END, baselineMaxRMS, &blOk);
            if (!blOk) { ++st.nDirtyBaseline; continue; }
            if (useFilter) af = butterworthLowPass(af, cutoffMHz, fsMHz);

            double ampMax = *std::max_element(af.begin() + jStart, af.begin() + jEnd + 1);
            hAmp->Fill(ampMax);
            ++st.nEvents;
            if (ampMax >= letThr) ++st.nHit;
        }

        tr->SetCacheSize(0);
        f->Close();
        delete f;
    }

done:
    bar.done();
    if (st.nEvents > 0) {
        st.eff = (double)st.nHit / (double)st.nEvents;
        st.effErr = std::sqrt(st.eff * (1.0 - st.eff) / st.nEvents);
    }

    TCanvas* c = new TCanvas(Form("cBinary_%s", tag.c_str()), "", 900, 550);
    c->SetGrid();
    c->SetLeftMargin(PAD_LEFT);
    c->SetRightMargin(PAD_RIGHT);
    c->SetBottomMargin(PAD_BOTTOM);
    c->SetTopMargin(PAD_TOP);
    hAmp->SetLineColor(kAzure + 1);
    hAmp->SetLineWidth(2);
    hAmp->Draw("HIST");

    double yMax = std::max(1.0, hAmp->GetMaximum() * 1.15);
    hAmp->SetMaximum(yMax);
    TLine* thr = new TLine(letThr, 0.0, letThr, yMax);
    thr->SetLineColor(kRed + 1);
    thr->SetLineStyle(2);
    thr->SetLineWidth(2);
    thr->Draw("same");

    TPaveText* pt = new TPaveText(0.56, 0.68, 0.94, 0.90, "NDC");
    pt->SetBorderSize(1);
    pt->SetFillColor(0);
    pt->SetTextFont(42);
    pt->SetTextSize(0.032);
    pt->AddText(Form("N events = %ld", st.nEvents));
    pt->AddText(Form("N hit = %ld", st.nHit));
    pt->AddText(Form("P(hit) = %.5f #pm %.5f", st.eff, st.effErr));
    pt->AddText(Form("thr = %.2f mV", letThr));
    pt->Draw();

    c->Update();
    c->Modified();
    ctx.savePNG(c, Form("binary_amp_%s.png", tag.c_str()));

    std::string tmp = Form("/tmp/alcor_amp_%s.root", tag.c_str());
    TFile* fTmp = new TFile(tmp.c_str(), "RECREATE");
    hAmp->Write("hAmp");
    fTmp->Write();
    fTmp->Close();
    delete fTmp;
    delete hAmp;
    st.ampTmpRoot = tmp;
    return st;
}

} // namespace

void sipm_alcor_binary_emulation() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();
    const std::string dataDir = DATA_DIR;

    std::cout << "\n+==========================================================+\n"
              << "|  SiPM ALCOR BINARY EMULATION                            |\n"
              << "+==========================================================+\n"
              << "  Data dir: " << dataDir << "\n"
              << "  Mode: no coincidence, no t0, hit/no-hit only\n\n";

    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        while (true) {
            std::cout << prompt << std::flush;
            if (!std::getline(std::cin, line)) {
                std::cerr << "\n[ERROR] stdin closed or unexpected EOF.\n";
                std::exit(1);
            }
            auto b = line.find_first_not_of(" \t\r\n");
            if (b != std::string::npos) return line.substr(b);
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
        std::string line = readLine("Thresholds in p.e. (es. 0.5 1.0 2.0):\n> ");
        std::stringstream ss(line);
        double v = 0.0;
        while (ss >> v) if (v > 0) fracs_pe.push_back(v);
    }
    if (fracs_pe.empty()) fracs_pe.push_back(1.0);

    bool useFilter = readYN("LP filter? [y/n]: ");

    std::map<int, CalibResult> calMap;
    for (int vbias : vbiasList) {
        auto cutoffs = findCalibratedCutoffs(vbias, dataDir);
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
                    double req = std::stod(line);
                    double best = chosen;
                    double dBest = std::numeric_limits<double>::max();
                    for (double co : cutoffs) {
                        double d = std::abs(co - req);
                        if (d < dBest) { dBest = d; best = co; }
                    }
                    chosen = best;
                } catch (...) {}
            }
        }

        CalibResult cal;
        if (!loadCalibration(cal, vbias, chosen, dataDir)) continue;
        if (cal.m <= 0 || cal.m > 100) continue;
        calMap[vbias] = cal;
        std::cout << "  Vbias=" << vbias << " gain=" << cal.m << " mV/p.e. trig=["
                  << cal.t_trig_start << "," << cal.t_trig_end << "]\n";
    }
    if (calMap.empty()) {
        std::cerr << "[ERROR] No calibration loaded.\n";
        return;
    }

    std::cout << "\nWindow mode:\n"
              << "  [0] Full waveform (recommended for unknown arrival time)\n"
              << "  [1] Use calibrated trigger window\n";
    int winMode = 0;
    while (true) {
        std::string line = readLine("Choice [0/1]: ");
        try {
            winMode = std::stoi(line);
            if (winMode == 0 || winMode == 1) break;
        } catch (...) {}
        std::cerr << "  [!] Enter 0 or 1.\n";
    }

    std::cout << "\n+==========================================================+\n"
              << "|  Binary mode summary                                    |\n"
              << "+==========================================================+\n|  Vbias: ";
    for (int v : vbiasList) std::cout << v << " ";
    std::cout << "V\n|  Thresholds (p.e.): ";
    for (double f : fracs_pe) std::cout << f << " ";
    std::cout << "\n|  LP filter: " << (useFilter ? "yes" : "NO")
              << "\n|  Window: " << (winMode == 0 ? "FULL" : "CALIB TRIG")
              << "\n+==========================================================+\n\n";

    auto t0 = std::chrono::steady_clock::now();
    for (int vbias : vbiasList) {
        auto itCal = calMap.find(vbias);
        if (itCal == calMap.end()) continue;
        const CalibResult& cal = itCal->second;

        auto runs = findRunsForVbias(vbias, dataDir);
        if (runs.empty()) {
            std::cerr << "  [WARN] Vbias=" << vbias << ": no runs found.\n";
            continue;
        }

        double fsMHz = estimateFsMHzFromFirstRun(runs.begin()->second);
        if (fsMHz <= 0) {
            std::cerr << "  [WARN] Vbias=" << vbias << ": invalid fs, skipped.\n";
            continue;
        }

        int jStart = 0, jEnd = 1023;
        if (winMode == 1) {
            TFile* f0 = TFile::Open(runs.begin()->second.c_str(), "READ");
            if (f0 && !f0->IsZombie()) {
                TTree* tr = (TTree*)f0->Get("ch1");
                if (tr) {
                    const int N = 1024;
                    Double_t t[N];
                    tr->SetBranchAddress("time", t);
                    tr->GetEntry(0);
                    triggerWindowIndices(t, N, cal.t_trig_start, cal.t_trig_end, jStart, jEnd, "binary");
                }
                f0->Close();
            }
            delete f0;
        }

        std::cout << "+-- Vbias=" << vbias << " V  " << runs.size() << " run(s)\n"
                  << "  Gain=" << cal.m << " mV/p.e.  Offset=" << cal.q << " mV\n"
                  << "  Window idx: [" << jStart << "," << jEnd << "]\n";

        std::vector<double> xThr, yEff, yErr;
        for (double frac : fracs_pe) {
            const double letThr = cal.q + frac * cal.m;
            std::string tag = Form("vbias%d_let%.2fpe", vbias, frac);

            std::cout << "\n  +-- LET=" << frac << " p.e.  thr=" << std::fixed << std::setprecision(2)
                      << letThr << " mV\n";
            BinaryStats bs = analyseBinaryOneLET(runs, jStart, jEnd, letThr, useFilter,
                                                 cal.cutoff_MHz, fsMHz, tag, ctx);
            std::cout << "  [BIN] events=" << bs.nEvents
                      << " hit=" << bs.nHit
                      << " dirtyBL=" << bs.nDirtyBaseline
                      << " P(hit)=" << std::setprecision(5) << bs.eff
                      << " +/- " << bs.effErr << "\n";

            xThr.push_back(frac);
            yEff.push_back(bs.eff);
            yErr.push_back(bs.effErr);
            if (!bs.ampTmpRoot.empty()) gSystem->Unlink(bs.ampTmpRoot.c_str());
        }

        if (!xThr.empty()) {
            TGraph* gr = new TGraph((int)xThr.size(), xThr.data(), yEff.data());
            TCanvas* cEff = new TCanvas(Form("cEff_vbias%d", vbias), "", 850, 520);
            cEff->SetGrid();
            cEff->SetLeftMargin(PAD_LEFT);
            cEff->SetRightMargin(PAD_RIGHT);
            cEff->SetBottomMargin(PAD_BOTTOM);
            cEff->SetTopMargin(PAD_TOP);
            gr->SetTitle(Form("ALCOR binary firing probability — Vbias %d V;Threshold (p.e.);P(hit)", vbias));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(1.0);
            gr->SetMarkerColor(kAzure + 1);
            gr->SetLineColor(kAzure + 1);
            gr->SetLineWidth(2);
            gr->Draw("APL");
            gr->GetYaxis()->SetRangeUser(0.0, 1.02);
            cEff->Update();
            cEff->Modified();
            ctx.savePNG(cEff, Form("binary_eff_vs_pe_vbias%d.png", vbias));
            delete gr;
        }
    }

    std::cout << "\n+==========================================================+\n|  DONE  "
              << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
              << " s\n|  PNG:  " << ctx.pngDir
              << "\n|  ROOT: " << ctx.rootDir
              << "\n+==========================================================+\n";

    ctx.reopenAllCanvases(40);
    std::cout << "\n  Canvases reopened and editable.\n"
              << "  All canvas .root files: " << ctx.canvasDir << "\n";
}

