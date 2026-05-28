/**
 * sipm_dark_noise_map.cpp
 * =======================
 * Analisi noise/DCR su prese dati senza laser (vbias=0 o qualsiasi).
 *
 * MISURE PER POSIZIONE:
 *   - noise_rms   : RMS della baseline (mV)
 *   - peak_to_peak: ampiezza picco-picco nella finestra di osservazione (mV)
 *   - dcr_rate    : N_crossing / (N_waveform * T_window) in Hz
 *   - n_crossing  : conteggio grezzo dei crossing
 *   - n_waveforms : numero di waveform processate
 *
 * INPUT:
 *   <dataDir>/data_x_<X>_y_<Y>_vbias_<V>_run_<N>.root  (multirun automatico)
 *   oppure
 *   <dataDir>/data_x_<X>_y_<Y>_vbias_<V>.root
 *
 * OUTPUT:
 *   <dataDir>/dark_noise_map.root  con TTree "dark_map"
 *   PNG: dark_noise_rms_map.png, dark_dcr_rate_map.png,
 *        dark_peak_to_peak_map.png, dark_n_crossing_map.png
 *
 * Compile:  .L sipm_dark_noise_map.cpp+
 * Run:      sipm_dark_noise_map()
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TROOT.h>

#include "../header/Config.h"
#include "../header/InputHelpers.h"

// ─── Costanti ────────────────────────────────────────────────────────────────
static const int    N_SAMPLES        = 1024;
static const double BASELINE_START_NS = 0.0;
static const double BASELINE_END_NS   = 30.0;
static const float  PAD_L = 0.13f, PAD_R = 0.08f,
                    PAD_T = 0.10f, PAD_B = 0.13f;

// ─── Structs ──────────────────────────────────────────────────────────────────
struct DarkPoint {
    double x = 0.0, y = 0.0;
    int    vbias = 0;
    double noise_rms    = 0.0;
    double peak_to_peak = 0.0;
    double dcr_rate     = 0.0;
    long   n_crossing   = 0;
    long   n_waveforms  = 0;
    double t_window_ns  = 0.0;
    bool   ok = false;
};

struct PositionKey {
    double x; double y; int vbias;
    bool operator<(const PositionKey& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return vbias < o.vbias;
    }
};

// ─── Parse filename ───────────────────────────────────────────────────────────
static bool parsePositionFile(const std::string& fname,
                               double& x, double& y, int& vbias, int& run)
{
    static const std::regex reRun(
        R"(^data_x_([\-\d]+)_y_([\-\d]+)_vbias_(\d+)_run_(\d+)\.root$)");
    static const std::regex reNoRun(
        R"(^data_x_([\-\d]+)_y_([\-\d]+)_vbias_(\d+)\.root$)");
    std::smatch m;
    if (std::regex_match(fname, m, reRun)) {
        try { x=std::stod(m[1]); y=std::stod(m[2]); vbias=std::stoi(m[3]); run=std::stoi(m[4]); }
        catch(...){ return false; }
        return true;
    }
    if (std::regex_match(fname, m, reNoRun)) {
        try { x=std::stod(m[1]); y=std::stod(m[2]); vbias=std::stoi(m[3]); run=0; }
        catch(...){ return false; }
        return true;
    }
    return false;
}

// ─── Mediana ──────────────────────────────────────────────────────────────────
static double medianVec(std::vector<double>& v) {
    if (v.empty()) return 0.0;
    std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
    return v[v.size()/2];
}

// ─── Analisi una posizione ────────────────────────────────────────────────────
static DarkPoint analyzePosition(const PositionKey& key,
                                  const std::map<int,std::string>& runs,
                                  double crossing_thr_mV,
                                  double obs_start_ns,
                                  double obs_end_ns)
{
    DarkPoint dp;
    dp.x = key.x; dp.y = key.y; dp.vbias = key.vbias;
    dp.t_window_ns = obs_end_ns - obs_start_ns;

    long   n_wf=0, n_cross=0;
    double sum_rms=0.0, sum_ptp=0.0;

    for (const auto& [runNum, path] : runs) {
        TFile* fIn = TFile::Open(path.c_str(), "READ");
        if (!fIn || fIn->IsZombie()) { delete fIn; continue; }
        TTree* tr = (TTree*)fIn->Get("ch1");
        if (!tr) { fIn->Close(); delete fIn; continue; }

        Double_t t1[N_SAMPLES], a1[N_SAMPLES];
        tr->SetBranchAddress("time",      t1);
        tr->SetBranchAddress("amplitude", a1);
        tr->SetCacheSize(2*1024*1024);

        Long64_t nEv = tr->GetEntries();
        for (Long64_t i = 0; i < nEv; ++i) {
            tr->GetEntry(i);

            // Baseline mediana
            std::vector<double> pre;
            pre.reserve(64);
            for (int j=0; j<N_SAMPLES; ++j)
                if (t1[j] >= BASELINE_START_NS && t1[j] < BASELINE_END_NS)
                    pre.push_back(a1[j]);
            if (pre.empty()) continue;
            double offset = medianVec(pre);

            // Noise RMS
            double sum2=0.0;
            for (double v : pre) { double d=v-offset; sum2+=d*d; }
            sum_rms += std::sqrt(sum2/pre.size());

            // Peak-to-peak nella finestra
            double vmin=1e9, vmax=-1e9;
            for (int j=0; j<N_SAMPLES; ++j) {
                if (t1[j] < obs_start_ns || t1[j] > obs_end_ns) continue;
                double v = a1[j]-offset;
                if (v<vmin) vmin=v;
                if (v>vmax) vmax=v;
            }
            if (vmax > vmin) sum_ptp += (vmax-vmin);

            // Crossing: primo attraversamento di crossing_thr_mV
            bool crossed=false;
            for (int j=1; j<N_SAMPLES && !crossed; ++j) {
                if (t1[j] < obs_start_ns) continue;
                if (t1[j] > obs_end_ns)   break;
                double v0=a1[j-1]-offset, v1=a1[j]-offset;
                if (v0 < crossing_thr_mV && v1 >= crossing_thr_mV) {
                    ++n_cross; crossed=true;
                }
            }
            ++n_wf;
        }
        fIn->Close(); delete fIn;
    }

    if (n_wf == 0) return dp;

    dp.noise_rms    = sum_rms / n_wf;
    dp.peak_to_peak = sum_ptp / n_wf;
    dp.n_crossing   = n_cross;
    dp.n_waveforms  = n_wf;
    double T_s = dp.t_window_ns * 1e-9;
    dp.dcr_rate = (T_s > 0.0) ? n_cross/(static_cast<double>(n_wf)*T_s) : 0.0;
    dp.ok = true;
    return dp;
}

// ─── Build edges TH2D ────────────────────────────────────────────────────────
static std::vector<double> buildEdges(const std::vector<double>& v) {
    std::vector<double> e;
    if (v.empty()) return e;
    if (v.size()==1) { e.push_back(v[0]-5.0); e.push_back(v[0]+5.0); return e; }
    e.push_back(v[0]-0.5*(v[1]-v[0]));
    for (size_t i=1;i<v.size();++i) e.push_back(0.5*(v[i-1]+v[i]));
    e.push_back(v.back()+0.5*(v.back()-v[v.size()-2]));
    return e;
}

// ─── Mappa 2D ────────────────────────────────────────────────────────────────
static void drawMap(const std::vector<DarkPoint>& pts,
                    const std::string& qty, const std::string& zTitle,
                    const std::string& outPath)
{
    std::set<double> uxs, uys;
    for (const auto& p : pts) if (p.ok) { uxs.insert(p.x); uys.insert(p.y); }
    if (uxs.empty()) return;

    std::vector<double> xv(uxs.begin(),uxs.end()), yv(uys.begin(),uys.end());
    auto xe=buildEdges(xv), ye=buildEdges(yv);

    TH2D* h = new TH2D(("h_"+qty).c_str(),
        Form(";x (mm);y (mm);%s",zTitle.c_str()),
        (int)xv.size(),xe.data(),(int)yv.size(),ye.data());
    h->SetDirectory(nullptr);

    for (const auto& p : pts) {
        if (!p.ok) continue;
        double val=0.0;
        if      (qty=="noise_rms")    val=p.noise_rms;
        else if (qty=="peak_to_peak") val=p.peak_to_peak;
        else if (qty=="dcr_rate")     val=p.dcr_rate;
        else if (qty=="n_crossing")   val=static_cast<double>(p.n_crossing);
        int bx=h->GetXaxis()->FindBin(p.x);
        int by=h->GetYaxis()->FindBin(p.y);
        h->SetBinContent(bx,by,val);
    }

    TCanvas* c=new TCanvas(("c_"+qty).c_str(),(zTitle+" map").c_str(),900,750);
    c->SetLeftMargin(PAD_L); c->SetRightMargin(PAD_R);
    c->SetBottomMargin(PAD_B); c->SetTopMargin(PAD_T);
    gStyle->SetPalette(kBird);
    h->SetContour(64);
    h->GetXaxis()->SetTitleSize(0.048f);
    h->GetYaxis()->SetTitleSize(0.048f);
    h->Draw("COLZ");
    c->Update(); c->Modified();
    c->SaveAs(outPath.c_str());
    delete c; delete h;
}

// ─── MAIN ─────────────────────────────────────────────────────────────────────
void sipm_dark_noise_map()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::cout << "\n+==========================================================+\n"
              << "|  SiPM DARK NOISE MAP — noise RMS, DCR, peak-to-peak      |\n"
              << "+==========================================================+\n\n";

    // ── 0. Selezione cartella (identica a sipm_pos_scan) ─────────────────────
    std::string dataDir;
    {
        const char* envFull = std::getenv("SIPM_DATA_DIR");
        if (envFull && envFull[0] != '\0') {
            dataDir = envFull;
        } else {
            // Deriva root = parent di DATA_DIR
            std::string root = DATA_DIR;
            while (!root.empty() && root.back()=='/') root.pop_back();
            const size_t sl = root.find_last_of("/\\");
            if (sl != std::string::npos) root = root.substr(0,sl);

            std::vector<std::string> subs;
            void* dp = gSystem->OpenDirectory(root.c_str());
            if (dp) {
                const char* ent=nullptr;
                while ((ent=gSystem->GetDirEntry(dp))!=nullptr) {
                    std::string s(ent);
                    if (s=="."||s=="..") continue;
                    FileStat_t st;
                    if (gSystem->GetPathInfo((root+"/"+s).c_str(),st)==0 &&
                        R_ISDIR(st.fMode)) subs.push_back(s);
                }
                gSystem->FreeDirectory(dp);
            }
            std::sort(subs.begin(),subs.end());

            if (!subs.empty()) {
                std::cout << "  Cartelle in " << root << ":\n";
                for (size_t i=0;i<subs.size();++i)
                    std::cout << "    [" << (i+1) << "] " << subs[i] << "\n";
                const std::string l = readLineOrEmpty(
                    Form("\n  Scegli [n] o path [ENTER=1]: "));
                if (l.empty()) {
                    dataDir = root+"/"+subs[0];
                } else {
                    try {
                        size_t pos=0;
                        int idx=std::stoi(l,&pos);
                        if (pos==l.size()&&idx>=1&&idx<=(int)subs.size())
                            dataDir=root+"/"+subs[idx-1];
                        else dataDir=l;
                    } catch(...){ dataDir=l; }
                }
            } else {
                dataDir = readLine("  Cartella dati (path): ");
            }
        }

        while (!dataDir.empty() && dataDir.back()=='/') dataDir.pop_back();
        if (gSystem->AccessPathName(dataDir.c_str())) {
            std::cerr << "[ERR] Non accessibile: " << dataDir << "\n";
            return;
        }
        std::cout << "  --> " << dataDir << "\n\n";
    }

    // ── 1. Parametri analisi ──────────────────────────────────────────────────
    double crossing_thr = readDoubleDefault("Soglia crossing [mV, default=5.0]: ", 5.0);
    double obs_start    = readDoubleDefault("Finestra inizio [ns, default=0.0]: ",  0.0);
    double obs_end      = readDoubleDefault("Finestra fine   [ns, default=204.6]: ", 204.6);

    std::cout << "\nParametri:\n"
              << "  Soglia crossing : " << crossing_thr << " mV\n"
              << "  Finestra        : [" << obs_start << ", " << obs_end << "] ns\n\n";

    // ── 2. Scan cartella ──────────────────────────────────────────────────────
    std::map<PositionKey, std::map<int,std::string>> posMap;
    {
        void* dp = gSystem->OpenDirectory(dataDir.c_str());
        if (!dp) { std::cerr << "[ERR] Cannot open: " << dataDir << "\n"; return; }
        const char* ent=nullptr;
        while ((ent=gSystem->GetDirEntry(dp))!=nullptr) {
            std::string fn(ent);
            if (fn.size()<5||fn.substr(fn.size()-5)!=".root") continue;
            if (fn.find("data_x_")==std::string::npos) continue;
            if (fn.find("dark_noise")!=std::string::npos) continue;
            double x,y; int vbias,run;
            if (!parsePositionFile(fn,x,y,vbias,run)) continue;
            posMap[{x,y,vbias}][run] = dataDir+"/"+fn;
        }
        gSystem->FreeDirectory(dp);
    }

    if (posMap.empty()) {
        std::cerr << "[ERR] Nessun file trovato in: " << dataDir << "\n";
        return;
    }
    std::cout << "Trovate " << posMap.size() << " posizioni.\n\n";

    // ── 3. Analisi per posizione ──────────────────────────────────────────────
    std::vector<DarkPoint> results;
    results.reserve(posMap.size());
    int iPos=0;
    for (const auto& [key,runs] : posMap) {
        ++iPos;
        std::cout << Form("[%3d/%3d] x=%6.1f y=%6.1f vbias=%d runs=%d  ",
                          iPos,(int)posMap.size(),
                          key.x,key.y,key.vbias,(int)runs.size());
        std::cout.flush();

        DarkPoint dp = analyzePosition(key,runs,crossing_thr,obs_start,obs_end);

        if (dp.ok)
            std::cout << Form("rms=%.2fmV  ptp=%.2fmV  dcr=%.1fHz  cross=%ld/%ld\n",
                              dp.noise_rms,dp.peak_to_peak,
                              dp.dcr_rate,dp.n_crossing,dp.n_waveforms);
        else
            std::cout << "[FAIL]\n";
        results.push_back(dp);
    }

    // ── 4. Scrittura ROOT ─────────────────────────────────────────────────────
    std::string outRoot = dataDir+"/dark_noise_map.root";
    {
        TFile* fOut = new TFile(outRoot.c_str(),"RECREATE");
        if (!fOut||fOut->IsZombie()) {
            std::cerr << "[ERR] Cannot create: " << outRoot << "\n";
            delete fOut; return;
        }

        TParameter<double> pThr("crossing_thr_mV",crossing_thr);
        TParameter<double> pO0("obs_start_ns",obs_start);
        TParameter<double> pO1("obs_end_ns",obs_end);
        pThr.Write(); pO0.Write(); pO1.Write();

        Double_t b_x,b_y,b_rms,b_ptp,b_dcr,b_twin;
        Int_t b_vbias; Long64_t b_ncross,b_nwf;

        TTree* tOut=new TTree("dark_map","Dark noise map");
        tOut->Branch("x",           &b_x,    "x/D");
        tOut->Branch("y",           &b_y,    "y/D");
        tOut->Branch("vbias",       &b_vbias,"vbias/I");
        tOut->Branch("noise_rms",   &b_rms,  "noise_rms/D");
        tOut->Branch("peak_to_peak",&b_ptp,  "peak_to_peak/D");
        tOut->Branch("dcr_rate",    &b_dcr,  "dcr_rate/D");
        tOut->Branch("n_crossing",  &b_ncross,"n_crossing/L");
        tOut->Branch("n_waveforms", &b_nwf,  "n_waveforms/L");
        tOut->Branch("t_window_ns", &b_twin, "t_window_ns/D");

        for (const auto& dp : results) {
            if (!dp.ok) continue;
            b_x=dp.x; b_y=dp.y; b_vbias=dp.vbias;
            b_rms=dp.noise_rms; b_ptp=dp.peak_to_peak;
            b_dcr=dp.dcr_rate;
            b_ncross=dp.n_crossing; b_nwf=dp.n_waveforms;
            b_twin=dp.t_window_ns;
            tOut->Fill();
        }
        fOut->Write(); fOut->Close(); delete fOut;
        std::cout << "\nScritto: " << outRoot << "\n";
    }

    // ── 5. PNG ────────────────────────────────────────────────────────────────
    std::string pb = dataDir+"/dark_";
    drawMap(results,"noise_rms",   "Noise RMS (mV)",       pb+"noise_rms_map.png");
    drawMap(results,"peak_to_peak","Peak-to-peak (mV)",    pb+"peak_to_peak_map.png");
    drawMap(results,"dcr_rate",    "DCR rate (Hz)",        pb+"dcr_rate_map.png");
    drawMap(results,"n_crossing",  "N crossing (counts)",  pb+"n_crossing_map.png");
    std::cout << "PNG salvati in: " << dataDir << "/dark_*.png\n";

    // ── 6. Riepilogo ──────────────────────────────────────────────────────────
    double m_rms=0,m_dcr=0; int n_ok=0;
    for (const auto& dp : results) if (dp.ok) {
        m_rms+=dp.noise_rms; m_dcr+=dp.dcr_rate; ++n_ok;
    }
    if (n_ok>0) {
        m_rms/=n_ok; m_dcr/=n_ok;
        std::cout << Form("\n=== RIEPILOGO (%d posizioni) ===\n",n_ok)
                  << Form("  Noise RMS medio : %.2f mV\n",m_rms)
                  << Form("  DCR rate medio  : %.1f Hz\n",m_dcr)
                  << Form("  Finestra        : %.1f ns\n",obs_end-obs_start)
                  << Form("  Soglia          : %.1f mV\n",crossing_thr)
                  << "=================================\n";
    }
}
