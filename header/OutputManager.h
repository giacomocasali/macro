#pragma once
// OutputManager.h
// Output directory management and sequential file numbering.
//
// Expected layout (canvas/ and file_root/ created once by the user):
//   macro/
//   ├── 04_multi_scan/          <- program runs here
//   ├── canvas/                 <- PNG output base
//   └── file_root/              <- ROOT output base
//
// The program creates timestamped subdirectories automatically:
//   macro/canvas/YYYYMMDD_HHMMSS/    <- numbered PNGs
//   macro/file_root/YYYYMMDD_HHMMSS/ <- numbered ROOT files
//
// Usage:
//   OutCtx ctx = createOutputDirs();
//   canvas->SaveAs(ctx.png("tot_map_vbias55.png").c_str());
//   // → macro/canvas/20260401_143022/01_tot_map_vbias55.png
//   ctx.saveRoot(h2D, h2D_corr, "tag", 1.0);
//   // → macro/file_root/20260401_143022/14_tot_tag_let1.00pe.root

#include <string>
#include <vector>
#include <ctime>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TKey.h>

// OutCtx: holds the three output directories and a global sequential counter.
// One instance created in main, passed by reference to all functions that save files.
//
// Strategy: every canvas is saved as PNG (for quick viewing) AND as a .root
// file (for interactive editing later).  After saving, the canvas is deleted
// from memory so that ROOT does not accumulate hundreds of MB of pixmaps.
// At the very end of the analysis, reopenAllCanvases() can reload them all
// from the .root files so you get the interactive session you want.
struct OutCtx {
    std::string pngDir;      // macro/canvas/YYYYMMDD_HHMMSS
    std::string rootDir;     // macro/file_root/YYYYMMDD_HHMMSS
    std::string canvasDir;   // macro/file_root/YYYYMMDD_HHMMSS/canvases
    int         n = 0;       // global sequential counter

    // Paths of saved canvas .root files, for reopening at the end.
    std::vector<std::string> savedCanvasPaths;

    // Returns the full path for the next PNG and increments the counter.
    std::string png(const std::string& name) {
        return Form("%s/%02d_%s", pngDir.c_str(), ++n, name.c_str());
    }

    // Save canvas as PNG + ROOT, then delete it from memory.
    // The .root file preserves the full TCanvas with all pads, histograms,
    // fits, and TPaveText — you can reopen it with TFile and modify anything.
    void savePNG(TCanvas* c, const std::string& name) {
        if (!c) return;
        // 1. PNG for quick viewing
        std::string pngPath = png(name);
        c->SaveAs(pngPath.c_str());

        // 2. ROOT file for later interactive editing
        std::string rootName = name;
        // .png → .root
        if (rootName.size() > 4 &&
            rootName.substr(rootName.size()-4) == ".png")
            rootName = rootName.substr(0, rootName.size()-4) + ".root";
        else
            rootName += ".root";
        std::string rootPath = Form("%s/%02d_%s",
                                    canvasDir.c_str(), n, rootName.c_str());
        {
            TFile fOut(rootPath.c_str(), "RECREATE");
            c->Write();
            fOut.Close();
        }
        savedCanvasPaths.push_back(rootPath);

        // 3. Free memory — this is the key to not crashing
        delete c;
    }

    // Close and delete all canvases still in memory (safety net).
    static void closeAllCanvases() {
        gROOT->GetListOfCanvases()->Delete();
    }

    // Reopen all saved canvases from their .root files.
    // Call once at the very end of the analysis for interactive inspection.
    // Only loads the last `maxReopen` canvases to avoid re-saturating RAM.
    void reopenAllCanvases(int maxReopen = 40) const {
        int total = (int)savedCanvasPaths.size();
        int skip  = std::max(0, total - maxReopen);
        if (skip > 0)
            std::cout << "  [OUT] " << total << " canvases saved, "
                      << "reopening last " << maxReopen << ".\n"
                      << "        Earlier ones are in: " << canvasDir << "\n";
        int loaded = 0;
        for (int i = skip; i < total; ++i) {
            TFile* f = TFile::Open(savedCanvasPaths[i].c_str(), "READ");
            if (!f || f->IsZombie()) { delete f; continue; }
            TIter next(f->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)next())) {
                if (TString(key->GetClassName()) == "TCanvas") {
                    TCanvas* c = (TCanvas*)key->ReadObj();
                    if (c) {
                        c->SetBit(TObject::kCanDelete, false);
                        c->Draw();
                        ++loaded;
                    }
                }
            }
            // Don't close f — ROOT needs the file open while the canvas lives.
            // The file will be closed when ROOT exits.
        }
        std::cout << "  [OUT] Reopened " << loaded << " canvases for inspection.\n"
                  << "        All " << total << " are also saved as .root in:\n"
                  << "        " << canvasDir << "\n";
    }

    // Returns the full path for the next ROOT file and increments the counter.
    std::string root(const std::string& name) {
        return Form("%s/%02d_%s", rootDir.c_str(), ++n, name.c_str());
    }

    // Save TOT maps (raw and optionally corrected) to a numbered ROOT file.
    void saveRoot(TH2D* h2D, TH2D* h2D_corr,
                  const std::string& tag, double frac_pe) {
        std::string fname = root(Form("tot_%s_let%.2fpe.root",
                                      tag.c_str(), frac_pe));
        TFile* fout = new TFile(fname.c_str(), "RECREATE");
        if (h2D)      h2D->Write("h2D");
        if (h2D_corr) h2D_corr->Write("h2D_corrected");
        fout->Write(); fout->Close(); delete fout;
        std::cout << "  Root file: " << fname << "\n";
    }
};

// Creates timestamped subdirectories and returns a ready-to-use OutCtx.
// Call once at the start of main.
static OutCtx createOutputDirs() {
    time_t now = time(nullptr);
    struct tm* t = localtime(&now);
    char ts[32];
    strftime(ts, sizeof(ts), "%Y%m%d_%H%M%S", t);

    std::string pngDir  = "../canvas/"    + std::string(ts);
    std::string rootDir = "../file_root/" + std::string(ts);

    // Base directories (canvas/ and file_root/) must already exist
    if (gSystem->mkdir(pngDir .c_str(), true) != 0)
        std::cerr << "[OUT] Cannot create " << pngDir  << "\n";
    if (gSystem->mkdir(rootDir.c_str(), true) != 0)
        std::cerr << "[OUT] Cannot create " << rootDir << "\n";

    OutCtx ctx;
    ctx.pngDir    = pngDir;
    ctx.rootDir   = rootDir;
    ctx.canvasDir = rootDir + "/canvases";

    if (gSystem->mkdir(ctx.canvasDir.c_str(), true) != 0)
        std::cerr << "[OUT] Cannot create " << ctx.canvasDir << "\n";

    std::cout << "  PNG  folder : " << pngDir  << "\n"
              << "  ROOT folder : " << rootDir << "\n"
              << "  Canvas ROOT : " << ctx.canvasDir << "\n";
    return ctx;
}
