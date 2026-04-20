// plot_resolution_vs_threshold_stats.cpp
// Legge N dal tot_map in file_root/ (non dalla cache, che può non esistere)

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <cmath>
#include <regex>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH2D.h>

struct ResPoint {
    double threshold_pe = 0.0;
    double sigma_ns     = 0.0;
    double sigma_err_ns = 0.0;
    long   n_events     = -1;
};

// Cerca tot_map in qualsiasi sottocartella di rootDir
static long readNfromTotMapDir(const std::string& rootDir, int vbias, double thr) {
    // pattern: tot_vbias55_let0.85pe_let...root
    char thrBuf[32];
    snprintf(thrBuf, sizeof(thrBuf), "%.2f", thr);
    std::string want = Form("tot_vbias%d_let%spe_", vbias, thrBuf);

    // Scansiona rootDir e una profondità di sottocartelle
    std::vector<std::string> searchDirs = {rootDir};
    {
        void* d = gSystem->OpenDirectory(rootDir.c_str());
        if (d) {
            const char* e = nullptr;
            while ((e = gSystem->GetDirEntry(d)) != nullptr) {
                std::string s(e);
                if (s == "." || s == "..") continue;
                searchDirs.push_back(rootDir + "/" + s);
            }
            gSystem->FreeDirectory(d);
        }
    }

    for (const auto& dir : searchDirs) {
        void* d = gSystem->OpenDirectory(dir.c_str());
        if (!d) continue;
        const char* e = nullptr;
        while ((e = gSystem->GetDirEntry(d)) != nullptr) {
            std::string fname(e);
            if (fname.find(want) == std::string::npos) continue;
            if (fname.find(".root") == std::string::npos) continue;
            std::string full = dir + "/" + fname;
            TFile* f = TFile::Open(full.c_str(), "READ");
            if (!f || f->IsZombie()) { delete f; continue; }
            TH2D* h2 = (TH2D*)f->Get("h2D_corrected");
            if (!h2) h2 = (TH2D*)f->Get("h2D");
            long n = h2 ? (long)h2->GetEntries() : -1;
            f->Close(); delete f;
            gSystem->FreeDirectory(d);
            return n;
        }
        gSystem->FreeDirectory(d);
    }
    return -1;
}

void plot_resolution_vs_threshold_stats(
    const char* dbPath   = "../file_root/resolution_scan_results.root",
    const char* rootDir  = "../file_root")
{
    gStyle->SetOptStat(0);

    TFile* f = TFile::Open(dbPath, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Cannot open: " << dbPath << "\n";
        delete f; return;
    }
    TTree* t = (TTree*)f->Get("resolution_scan");
    if (!t) {
        std::cerr << "[ERROR] TTree 'resolution_scan' not found\n";
        f->Close(); delete f; return;
    }

    int vbias = 0;
    double threshold_pe = 0.0, sigma_ns = 0.0, sigma_err_ns = 0.0;
    t->SetBranchAddress("vbias",         &vbias);
    t->SetBranchAddress("threshold_pe",  &threshold_pe);
    t->SetBranchAddress("sigma_ns",      &sigma_ns);
    t->SetBranchAddress("sigma_err_ns",  &sigma_err_ns);

    std::map<int, std::vector<ResPoint>> byVbias;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        long n = readNfromTotMapDir(rootDir, vbias, threshold_pe);
        byVbias[vbias].push_back({threshold_pe, sigma_ns, sigma_err_ns, n});
    }
    f->Close(); delete f;

    if (byVbias.empty()) { std::cerr << "[WARN] No rows.\n"; return; }

    // Sort each vbias by threshold
    for (auto& kv : byVbias)
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const ResPoint& a, const ResPoint& b){ return a.threshold_pe < b.threshold_pe; });

    // Check if any N is available
    bool hasN = false;
    for (auto& kv : byVbias)
        for (auto& p : kv.second)
            if (p.n_events > 0) { hasN = true; break; }

    int nPads = hasN ? 2 : 1;
    TCanvas* c = new TCanvas("c_res_vs_thr_stats", "Resolution vs threshold",
                              1000, hasN ? 800 : 500);
    if (hasN) c->Divide(1, 2);

    // ── PAD 1: Resolution ──────────────────────────────────────
    c->cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(hasN ? 0.08 : 0.12);

    TMultiGraph* mg  = new TMultiGraph();
    TLegend*     leg = new TLegend(0.62, 0.68, 0.90, 0.88);
    leg->SetBorderSize(1);

    int color = 2;
    for (auto& kv : byVbias) {
        int vb = kv.first;
        std::vector<double> x, y, ex, ey;
        for (auto& p : kv.second) {
            x.push_back(p.threshold_pe);
            y.push_back(p.sigma_ns);
            ex.push_back(0.0);
            ey.push_back(p.sigma_err_ns);
        }
        auto* gr = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
        gr->SetName(Form("gr_vb%d", vb));
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
        gr->SetMarkerColor(color); gr->SetLineColor(color); gr->SetLineWidth(2);
        mg->Add(gr, "LP");
        leg->AddEntry(gr, Form("Vbias %d V", vb), "lp");
        ++color; if (color == 5) ++color;
    }
    mg->SetTitle("Timing resolution vs p.e. threshold;Threshold (p.e.);#sigma_{#Delta t} (ns)");
    mg->Draw("A");
    mg->GetYaxis()->SetTitleOffset(1.4);
    leg->Draw();

    // ── PAD 2: N events ───────────────────────────────────────
    if (hasN) {
        c->cd(2);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogy();

        TMultiGraph* mgN = new TMultiGraph();
        color = 2;
        for (auto& kv : byVbias) {
            int vb = kv.first;
            std::vector<double> x, y, ex, ey;
            for (auto& p : kv.second) {
                if (p.n_events <= 0) continue;
                x.push_back(p.threshold_pe);
                y.push_back((double)p.n_events);
                ex.push_back(0.0); ey.push_back(0.0);
            }
            if (x.empty()) { ++color; if (color==5) ++color; continue; }
            auto* gr = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
            gr->SetName(Form("grN_vb%d", vb));
            gr->SetMarkerStyle(21); gr->SetMarkerSize(1.0);
            gr->SetMarkerColor(color); gr->SetLineColor(color); gr->SetLineWidth(2);
            mgN->Add(gr, "LP");
            ++color; if (color==5) ++color;
        }
        mgN->SetTitle("Statistics per threshold;Threshold (p.e.);N events accepted");
        mgN->Draw("A");
        mgN->GetYaxis()->SetTitleOffset(1.4);
    } else {
        std::cout << "[WARN] No N found in tot_maps — skipping stats pad\n";
    }

    gSystem->mkdir("../canvas", true);
    const char* out = "../canvas/resolution_vs_threshold_with_stats.png";
    c->SaveAs(out);
    std::cout << "[OK] Plot saved: " << out << "\n";
    std::cout << "[OK] Input DB  : " << dbPath << "\n";
}
