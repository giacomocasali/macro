/**
 * sipm_vbias_summary.cpp
 * ======================
 * Two-in-one summary plot vs V_bias:
 *
 *   Canvas 1 — BREAKDOWN VOLTAGE
 *     Reads calib_vbias<V>_cut<C>mhz.root, fits gain vs Vbias linearly,
 *     extracts V_bd = -b/a (intercept with gain = 0).
 *
 *   Canvas 2 — DETECTION PROBABILITY vs Vbias  (optional)
 *     Reads map_results_vbias<V>_let<F>_pdet_*.root produced by
 *     sipm_pos_scan(), filters for a user-chosen position (x, y),
 *     plots P_det vs Vbias for every frac_pe found.
 *     If no map_results files are found, this canvas is skipped with a warning.
 *
 * Compile:  .L sipm_vbias_summary.cpp+
 * Run:      sipm_vbias_summary()
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TParameter.h>

#include "../header/Config.h"
#include "../header/InputHelpers.h"
#include "../header/CalibIO.h"
#include "../header/OutputManager.h"

// ── colour palette ────────────────────────────────────────────────────────
static const int PALETTE[] = {
    kAzure+1, kRed+1, kGreen+2, kOrange+7,
    kMagenta+1, kCyan+2, kViolet+1, kYellow+2
};
static int paletteColor(int i) { return PALETTE[i % 8]; }

// ════════════════════════════════════════════════════════════════════════════
//  helpers
// ════════════════════════════════════════════════════════════════════════════

struct CalibFile { int vbias; double cutoff_MHz; };

static std::vector<CalibFile> findCalibFiles(const std::string& dir)
{
    static const std::regex re(R"(^calib_vbias(\d+)_cut(\d+)mhz\.root$)");
    std::vector<CalibFile> v;
    void* dp = gSystem->OpenDirectory(dir.c_str());
    if (!dp) return v;
    const char* ent;
    while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
        std::smatch m;
        std::string fn(ent);
        if (!std::regex_match(fn, m, re)) continue;
        try { v.push_back({std::stoi(m[1]), std::stod(m[2])}); } catch (...) {}
    }
    gSystem->FreeDirectory(dp);
    std::sort(v.begin(), v.end(),
        [](const CalibFile& a, const CalibFile& b){ return a.vbias < b.vbias; });
    return v;
}

struct MapFile { int vbias; double frac_pe; std::string mode; std::string path; };

static std::vector<MapFile> findMapFiles(const std::string& dir)
{
    static const std::regex re(
        R"(^map_results_vbias(\d+)_let([\d\.]+?)(?:_loose)?(?:_(pdet_(?:accepted|crossing)))?\.root$)");
    std::vector<MapFile> v;
    void* dp = gSystem->OpenDirectory(dir.c_str());
    if (!dp) return v;
    const char* ent;
    while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
        std::smatch m;
        std::string fn(ent);
        if (!std::regex_match(fn, m, re)) continue;
        try {
            v.push_back({ std::stoi(m[1].str()),
                          std::stod(m[2].str()),
                          m[3].matched ? m[3].str() : std::string("unknown"),
                          dir + "/" + fn });
        } catch (...) {}
    }
    gSystem->FreeDirectory(dp);
    std::sort(v.begin(), v.end(),
        [](const MapFile& a, const MapFile& b){
            return a.vbias != b.vbias ? a.vbias < b.vbias : a.frac_pe < b.frac_pe; });
    return v;
}

static long readLongParam(TFile* f, const char* name)
{
    { auto* p = dynamic_cast<TParameter<long>*>(f->Get(name)); if (p) return p->GetVal(); }
    { auto* p = dynamic_cast<TParameter<int>* >(f->Get(name)); if (p) return (long)p->GetVal(); }
    return -1L;
}

// Read p_det for position (x_t, y_t) from a map_results TTree.
// Returns -1 if not found.
static double readPdet(const std::string& path,
                        double x_t, double y_t, double tol,
                        long& n_laser_out, long& n_cross_out)
{
    n_laser_out = n_cross_out = -1;
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return -1; }
    TTree* t = static_cast<TTree*>(f->Get("map"));
    if (!t) { f->Close(); delete f; return -1; }

    Double_t x = 0, y = 0, pdet = -1;
    Long64_t nlas = 0, ncrs = 0;

    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("p_det", &pdet);
    if (t->GetBranch("n_laser"))    t->SetBranchAddress("n_laser",    &nlas);
    if (t->GetBranch("n_crossing")) t->SetBranchAddress("n_crossing", &ncrs);

    const Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (std::abs(x - x_t) < tol && std::abs(y - y_t) < tol) {
            n_laser_out = static_cast<long>(nlas);
            n_cross_out = static_cast<long>(ncrs);
            f->Close(); delete f;
            return pdet;
        }
    }
    f->Close(); delete f;
    return -1;
}

static double binomErr(double p, long n)
{
    if (n <= 0 || p < 0 || p > 1) return 0;
    return std::sqrt(p * (1.0 - p) / (double)n);
}

// ════════════════════════════════════════════════════════════════════════════
//  Canvas 1 — gain vs Vbias + V_bd fit
// ════════════════════════════════════════════════════════════════════════════
static void plotBreakdown(const std::string& dataDir, OutCtx& ctx)
{
    auto calibFiles = findCalibFiles(dataDir);
    if (calibFiles.empty()) {
        std::cerr << "  [BD] No calib files found — skipping breakdown canvas.\n";
        return;
    }

    // Pick cutoff with most Vbias points
    std::map<double,int> cnt;
    for (auto& cf : calibFiles) cnt[cf.cutoff_MHz]++;
    double cut = calibFiles[0].cutoff_MHz; int best = 0;
    for (auto& [c,n] : cnt) if (n > best) { best = n; cut = c; }

    struct GP { int vb; double g, eg; };
    std::vector<GP> pts;
    for (auto& cf : calibFiles) {
        if (cf.cutoff_MHz != cut) continue;
        CalibResult cal;
        if (!loadCalibration(cal, cf.vbias, cf.cutoff_MHz, dataDir)) continue;
        if (cal.m < GAIN_MIN || cal.m > GAIN_MAX) continue;
        pts.push_back({ cf.vbias, cal.m, cal.m * 0.03 });
        std::cout << "  [BD] Vbias=" << cf.vbias
                  << "  gain=" << std::fixed << std::setprecision(3)
                  << cal.m << " mV/p.e.\n";
    }

    if (pts.size() < 2) {
        std::cerr << "  [BD] Need >= 2 points. Got " << pts.size() << ".\n";
        return;
    }

    const int N = (int)pts.size();
    std::vector<double> vb(N), g(N), evb(N, 0), eg(N);
    for (int i = 0; i < N; ++i) { vb[i]=pts[i].vb; g[i]=pts[i].g; eg[i]=pts[i].eg; }

    TGraphErrors* gr = new TGraphErrors(N, vb.data(), g.data(), evb.data(), eg.data());
    gr->SetName("gr_gain"); gr->SetTitle(";V_{bias} (V);Gain (mV/p.e.)");
    gr->SetMarkerStyle(20); gr->SetMarkerSize(1.3);
    gr->SetMarkerColor(kAzure+1); gr->SetLineColor(kAzure+1); gr->SetLineWidth(2);

    double vlo = vb.front()-2, vhi = vb.back()+2;
    TF1* fl = new TF1("fl_gain","[0]*x+[1]", vlo, vhi);
    fl->SetParameters((g.back()-g.front())/(vb.back()-vb.front()), g.front());
    fl->SetLineColor(kRed+1); fl->SetLineWidth(2);
    TFitResultPtr fr = gr->Fit(fl,"SR");

    const double a=fl->GetParameter(0), b=fl->GetParameter(1);
    const double ea=fl->GetParError(0), eb=fl->GetParError(1);
    if (std::abs(a) < 1e-9) { std::cerr << "  [BD] Fit slope ~ 0.\n"; return; }

    const double Vbd     = -b / a;
    const double Vbd_err = std::abs(Vbd)*std::sqrt((ea/a)*(ea/a)+(eb/b)*(eb/b));
    const double chi2ndf = fr->Ndf() > 0 ? fr->Chi2()/fr->Ndf() : -1;

    std::cout << "\n  [BD] V_bd = " << std::fixed << std::setprecision(3)
              << Vbd << " ± " << Vbd_err << " V"
              << "  (chi2/ndf=" << chi2ndf << ")\n";

    TCanvas* c = new TCanvas("cBD","Breakdown voltage",850,600);
    c->SetGrid();
    c->SetLeftMargin(PAD_LEFT); c->SetRightMargin(PAD_RIGHT+0.02f);
    c->SetBottomMargin(PAD_BOTTOM); c->SetTopMargin(PAD_TOP);

    double yhi = g.back()*1.3;
    gr->GetXaxis()->SetRangeUser(std::min(vlo, Vbd-1), vhi);
    gr->GetYaxis()->SetRangeUser(0, yhi);
    gr->GetXaxis()->SetTitleSize(0.048f); gr->GetYaxis()->SetTitleSize(0.048f);
    gr->Draw("AP");

    TF1* fext = new TF1("fext_gain","[0]*x+[1]",
                         std::min(vlo,Vbd-0.5), vhi);
    fext->SetParameters(a,b);
    fext->SetLineColor(kRed+1); fext->SetLineWidth(2); fext->SetLineStyle(2);
    fext->Draw("same");

    TLine* lv = new TLine(Vbd, 0, Vbd, yhi*0.5);
    lv->SetLineColor(kGreen+2); lv->SetLineWidth(2); lv->SetLineStyle(2);
    lv->Draw();

    TPaveText* pt = new TPaveText(0.14,0.58,0.54,0.88,"NDC");
    pt->SetBorderSize(1); pt->SetFillColorAlpha(kWhite,0.85);
    pt->SetTextFont(42); pt->SetTextSize(0.036); pt->SetTextAlign(12);
    pt->AddText("Fit:  gain = a #times V_{bias} + b");
    pt->AddText(Form("a = %.4f #pm %.4f mV/p.e./V", a, ea));
    pt->AddText(Form("b = %.3f #pm %.3f mV/p.e.", b, eb));
    pt->AddText(Form("V_{bd} = %.3f #pm %.3f V", Vbd, Vbd_err));
    pt->AddText(Form("#chi^{2}/ndf = %.2f   cutoff = %.0f MHz", chi2ndf, cut));
    pt->Draw();

    TLatex* lx = new TLatex(Vbd, -yhi*0.055,
                              Form("V_{bd}=%.2fV",Vbd));
    lx->SetTextColor(kGreen+2); lx->SetTextSize(0.034);
    lx->SetTextAlign(22); lx->Draw();

    c->Update(); c->Modified();
    ctx.savePNG(c, "breakdown_voltage.png");
}

// ════════════════════════════════════════════════════════════════════════════
//  Canvas 2 — P_det vs Vbias for a fixed (x, y) position
// ════════════════════════════════════════════════════════════════════════════
static void plotPdetVsVbias(const std::string& dataDir, OutCtx& ctx)
{
    auto mapFiles = findMapFiles(dataDir);
    if (mapFiles.empty()) {
        std::cout << "  [Pdet] No map_results_*.root found — "
                     "run sipm_pos_scan() first.\n";
        return;
    }

    // Show available positions from the first file
    std::set<double> xs_avail, ys_avail;
    {
        TFile* f0 = TFile::Open(mapFiles[0].path.c_str(), "READ");
        if (f0 && !f0->IsZombie()) {
            TTree* t = static_cast<TTree*>(f0->Get("map"));
            if (t) {
                Double_t x=0, y=0;
                t->SetBranchAddress("x",&x); t->SetBranchAddress("y",&y);
                for (Long64_t i=0; i<t->GetEntries(); ++i) {
                    t->GetEntry(i); xs_avail.insert(x); ys_avail.insert(y);
                }
            }
            f0->Close(); delete f0;
        }
    }

    std::cout << "\n  [Pdet] Available x positions:";
    for (double v : xs_avail) std::cout << " " << v;
    std::cout << "\n  [Pdet] Available y positions:";
    for (double v : ys_avail) std::cout << " " << v;
    std::cout << "\n";

    // Ask position
    double x_t = 90.0, y_t = 0.0;
    {
        const std::string lx = readLineOrEmpty("  x position [ENTER = 90]: ");
        if (!lx.empty()) try { x_t = std::stod(lx); } catch (...) {}
        const std::string ly = readLineOrEmpty("  y position [ENTER = 0]: ");
        if (!ly.empty()) try { y_t = std::stod(ly); } catch (...) {}
    }
    const double tol = 1.0;
    std::cout << "  --> position (" << x_t << ", " << y_t
              << ")  tol=" << tol << " um\n\n";

    // Group by (frac_pe, mode)
    // key = (frac_pe, mode) → list of (vbias, p_det, err)
    std::map<std::pair<double,std::string>,
             std::vector<std::tuple<int,double,double>>> curves;

    for (auto& mf : mapFiles) {
        long nlas = -1, ncrs = -1;
        double p = readPdet(mf.path, x_t, y_t, tol, nlas, ncrs);
        if (p < 0) {
            std::cout << "  [Pdet] (" << x_t << "," << y_t << ") not found in "
                      << "vbias=" << mf.vbias
                      << " let=" << mf.frac_pe << "\n";
            continue;
        }
        long n_den = (mf.mode == "pdet_crossing" && ncrs > 0) ? ncrs : nlas;
        double err = binomErr(p, n_den);
        curves[{mf.frac_pe, mf.mode}].emplace_back(mf.vbias, p*100.0, err*100.0);
        std::cout << "  [Pdet] vbias=" << mf.vbias
                  << "  let=" << mf.frac_pe
                  << "  [" << mf.mode << "]"
                  << "  P_det=" << std::fixed << std::setprecision(1)
                  << p*100.0 << "%\n";
    }

    if (curves.empty()) {
        std::cerr << "  [Pdet] No data found for position ("
                  << x_t << "," << y_t << ").\n";
        return;
    }

    TCanvas* c = new TCanvas("cPdet",
        Form("P_{{det}} vs V_{{bias}}  pos(%.0f,%.0f)", x_t, y_t),
        900, 600);
    c->SetGrid();
    c->SetLeftMargin(PAD_LEFT); c->SetRightMargin(0.05f);
    c->SetBottomMargin(PAD_BOTTOM); c->SetTopMargin(PAD_TOP);

    TLegend* leg = new TLegend(0.14, 0.65, 0.55, 0.88);
    leg->SetBorderSize(1); leg->SetFillColorAlpha(kWhite,0.85);
    leg->SetTextFont(42); leg->SetTextSize(0.030);

    bool first = true;
    int ci = 0;
    // mode linestyle: accepted=solid, crossing=dashed
    auto modeStyle = [](const std::string& m) {
        return m == "pdet_crossing" ? 2 : 1;
    };
    auto modeLabel = [](const std::string& m) -> std::string {
        if (m == "pdet_crossing")  return "crossing";
        if (m == "pdet_accepted")  return "accepted";
        return m;
    };

    for (auto& [key, pts] : curves) {
        auto& [frac, mode] = key;
        std::sort(pts.begin(), pts.end());  // sort by vbias

        const int n = (int)pts.size();
        std::vector<double> vb(n), pd(n), evb(n,0), epd(n);
        for (int i=0; i<n; ++i) {
            vb[i]  = std::get<0>(pts[i]);
            pd[i]  = std::get<1>(pts[i]);
            epd[i] = std::get<2>(pts[i]);
        }

        TGraphErrors* gr = new TGraphErrors(n,
            vb.data(), pd.data(), evb.data(), epd.data());
        int col = paletteColor(ci++);
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
        gr->SetMarkerColor(col); gr->SetLineColor(col);
        gr->SetLineWidth(2); gr->SetLineStyle(modeStyle(mode));
        gr->SetTitle(Form(";V_{bias} (V);P_{det} (%%)"));

        if (first) {
            gr->GetXaxis()->SetTitleSize(0.048f);
            gr->GetYaxis()->SetTitleSize(0.048f);
            gr->GetYaxis()->SetRangeUser(0, 110);
            gr->Draw("ALP");
            first = false;
        } else {
            gr->Draw("LP same");
        }
        leg->AddEntry(gr,
            Form("LET=%.2f p.e. [%s]", frac, modeLabel(mode).c_str()), "lp");
    }

    // Position label
    TPaveText* ptPos = new TPaveText(0.57, 0.14, 0.94, 0.25, "NDC");
    ptPos->SetBorderSize(1); ptPos->SetFillColorAlpha(kWhite,0.85);
    ptPos->SetTextFont(42); ptPos->SetTextSize(0.034);
    ptPos->AddText(Form("Position: x=%.0f, y=%.0f #mum", x_t, y_t));
    ptPos->Draw();

    leg->Draw();
    c->Update(); c->Modified();
    ctx.savePNG(c, Form("pdet_vs_vbias_x%.0f_y%.0f.png", x_t, y_t));
}

// ════════════════════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════════════════════
void sipm_vbias_summary()
{
    g_data_dir_override = "";
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    std::cout << "\n+==========================================================+\n"
              << "|  SiPM VBIAS SUMMARY — breakdown + P_det vs Vbias        |\n"
              << "+==========================================================+\n\n";

    // ── data directory ───────────────────────────────────────────────────────
    std::string dataDir;
    {
        std::string root = DATA_DIR;
        while (!root.empty() && root.back() == '/') root.pop_back();
        const size_t sl = root.find_last_of("/\\");
        if (sl != std::string::npos) root = root.substr(0, sl);

        std::vector<std::string> subs;
        void* dp = gSystem->OpenDirectory(root.c_str());
        if (dp) {
            const char* ent;
            while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
                std::string s(ent);
                if (s == "." || s == "..") continue;
                FileStat_t st;
                if (gSystem->GetPathInfo((root+"/"+s).c_str(), st) == 0 &&
                    R_ISDIR(st.fMode)) subs.push_back(s);
            }
            gSystem->FreeDirectory(dp);
        }
        std::sort(subs.begin(), subs.end());

        int defIdx = -1;
        const std::string preferred = "giacomo_serpentina";
        for (size_t i = 0; i < subs.size(); ++i)
            if (subs[i] == preferred) { defIdx = (int)i+1; break; }
        if (defIdx < 0 && !subs.empty()) defIdx = 1;

        if (!subs.empty()) {
            std::cout << "  Folders in " << root << ":\n";
            for (size_t i = 0; i < subs.size(); ++i)
                std::cout << "    [" << (i+1) << "] " << subs[i]
                          << ((int)i+1 == defIdx ? "   <-- default" : "") << "\n";
            const std::string l = readLineOrEmpty(
                Form("\n  Choose [n] or path  [ENTER = %d]: ", defIdx));
            if (l.empty())
                dataDir = root + "/" + subs[defIdx-1];
            else {
                try {
                    size_t p=0; int idx=std::stoi(l,&p);
                    dataDir = (p==l.size() && idx>=1 && idx<=(int)subs.size())
                              ? root+"/"+subs[idx-1] : l;
                } catch (...) { dataDir = l; }
            }
        } else {
            dataDir = readLine("  Data folder path: ");
        }
        while (!dataDir.empty() && dataDir.back()=='/') dataDir.pop_back();
        if (gSystem->AccessPathName(dataDir.c_str())) {
            std::cerr << "[ERR] Not accessible: " << dataDir << "\n"; return;
        }
        g_data_dir_override = dataDir;
        std::cout << "  --> " << dataDir << "\n\n";
    }

    OutCtx ctx = createOutputDirs("vbias_summary");

    // Canvas 1: breakdown
    std::cout << "--- BREAKDOWN VOLTAGE ---\n";
    plotBreakdown(dataDir, ctx);

    // Canvas 2: P_det vs Vbias
    std::cout << "\n--- DETECTION PROBABILITY vs Vbias ---\n";
    plotPdetVsVbias(dataDir, ctx);

    std::cout << "\n+==========================================================+\n"
              << "|  DONE — canvas in: " << ctx.pngDir << "\n"
              << "+==========================================================+\n";
}
