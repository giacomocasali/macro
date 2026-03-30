#pragma once
// header/Config.h
// Global analysis parameters, pad geometry, and file metadata.
// Edit this file to change scan ranges, filter thresholds, or plot margins.

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TSystem.h>

// ════════════════════════════════════════════════════════════
//  SCAN & DERIVATIVE PARAMETERS
// ════════════════════════════════════════════════════════════
constexpr double MIN_THR          = -10.0; // minimum threshold (mV)
constexpr double MAX_THR          = 110.0; // maximum threshold (mV)
constexpr double STEP_SCAN        =   0.1; // threshold scan step (mV)
constexpr double STEP_DER         =   0.2; // derivative step for fit/peak finding (mV)
constexpr double MIN_POSITIVE_THR =   5.0; // spike-rejection: fit/zoom only above this (mV)
constexpr double FIT_SIGMA_MAX    =  15.0; // upper bound on Gaussian sigma in fit (mV)
constexpr double FIT_SIGMA_MIN    =   0.1; // lower bound on Gaussian sigma in fit (mV)
constexpr double MIN_RISING_MV    =   1.0; // min rising flank width for peak finder (mV)

// ════════════════════════════════════════════════════════════
//  BASELINE WINDOW  (fixed, hardware-determined)
//  Samples in [BASELINE_START, BASELINE_END) ns are used to
//  compute the median offset. Must be noise-only (pre-signal).
// ════════════════════════════════════════════════════════════
constexpr double BASELINE_START   =   0.0; // ns
constexpr double BASELINE_END     =  30.0; // ns

// ════════════════════════════════════════════════════════════
//  PAD MARGINS  (identical for all pads → same physical X scale)
// ════════════════════════════════════════════════════════════
constexpr float PAD_LEFT   = 0.13f;
constexpr float PAD_RIGHT  = 0.04f;
constexpr float PAD_TOP    = 0.10f;
constexpr float PAD_BOTTOM = 0.13f;

// ════════════════════════════════════════════════════════════
//  FILE METADATA
//  Parses filenames of the form  data.vbias_{55}_3.root
// ════════════════════════════════════════════════════════════
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
    void* dh = gSystem->OpenDirectory("../../data");
    if (!dh) { std::cerr << "Cannot open data directory.\n"; return files; }
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dh)) != nullptr) {
        std::string fname(entry);
        if (fname.find("data.vbias_{") == std::string::npos) continue;
        if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
        FileInfo info;
        if (parseFilename("../../data/" + fname, info)) {
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
