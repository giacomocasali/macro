#pragma once
// Config.h
// Global parameters, file metadata, and DATA_DIR.
// Edit DATA_DIR when you move or rename the data folder.
// All other constants control scan range, filter, baseline window, and plot margins.

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TSystem.h>

// ── Data directory — change this when you move the data ─────
// Use an absolute path or a path relative to the macro directory.
inline const std::string DATA_DIR =
    "/home/giacomo/caendgz-sipmanalysis/data/10 run filtro 5 53,54,55";

// ── Threshold scan parameters ────────────────────────────────
constexpr double MIN_THR          = -10.0; // mV — negative end shows pedestal in plot
constexpr double MAX_THR          = 110.0; // mV
constexpr double STEP_SCAN        =   0.1; // mV
constexpr double STEP_DER         =   0.2; // mV — derivative step
constexpr double MIN_POSITIVE_THR =   5.0; // mV — p.e. peaks searched above this
constexpr double FIT_SIGMA_MAX    =  15.0; // mV — Gaussian sigma upper bound
constexpr double FIT_SIGMA_MIN    =   0.1; // mV — Gaussian sigma lower bound
constexpr double MIN_RISING_MV    =   1.0; // mV — min rising flank for peak finder

// ── Baseline window (hardware-determined, must be pre-signal) ─
constexpr double BASELINE_START   =   0.0; // ns
constexpr double BASELINE_END     =  30.0; // ns

// ── Canvas margins (same for all pads) ───────────────────────
constexpr float PAD_LEFT   = 0.13f;
constexpr float PAD_RIGHT  = 0.04f;
constexpr float PAD_TOP    = 0.10f;
constexpr float PAD_BOTTOM = 0.13f;

// ── File metadata — parses data.vbias_{V}_F.root filenames ──
struct FileInfo {
    std::string path;
    double      vbias;
    double      filter;
    std::string tag;   // e.g. "vbias55_filter3"
};

static bool parseFilename(const std::string& path, FileInfo& info) {
    size_t slash = path.find_last_of("/\\");
    std::string fname = (slash == std::string::npos) ? path : path.substr(slash + 1);
    const std::string prefix = "data.vbias_{";
    size_t p0 = fname.find(prefix);
    if (p0 == std::string::npos) return false;
    size_t vStart = p0 + prefix.size();
    size_t vEnd   = fname.find('}', vStart);
    if (vEnd == std::string::npos) return false;
    std::string vStr = fname.substr(vStart, vEnd - vStart);
    size_t lStart = vEnd + 2;
    size_t lEnd   = fname.find(".root", lStart);
    if (lEnd == std::string::npos) return false;
    std::string lStr = fname.substr(lStart, lEnd - lStart);
    try { info.vbias = std::stod(vStr); info.filter = std::stod(lStr); }
    catch (...) { return false; }
    info.path = path;
    info.tag  = "vbias" + vStr + "_filter" + lStr;
    return true;
}

// Scans the current directory and returns all matching FileInfo entries,
// sorted by vbias then filter.
static std::vector<FileInfo> findInputFiles() {
    std::vector<FileInfo> files;
    void* dh = gSystem->OpenDirectory(DATA_DIR.c_str());
    if (!dh) { std::cerr << "Cannot open data directory: " << DATA_DIR << "\n"; return files; }
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dh)) != nullptr) {
        std::string fname(entry);
        if (fname.find("data.vbias_{") == std::string::npos) continue;
        if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
        FileInfo info;
        if (parseFilename(DATA_DIR + "/" + fname, info)) {
            files.push_back(info);
            std::cout << "Found: " << fname
                      << "  →  vbias=" << info.vbias
                      << " V,  filter=" << info.filter << "\n";
        } else {
            std::cerr << "Cannot parse filename: " << fname << "\n";
        }
    }
    gSystem->FreeDirectory(dh);
    std::sort(files.begin(), files.end(), [](const FileInfo& a, const FileInfo& b){
        return (a.vbias != b.vbias) ? a.vbias < b.vbias : a.filter < b.filter;
    });
    return files;
}
