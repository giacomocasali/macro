/**
 * sipm_detection_prob.cpp
 *
 * Calcola e plotta la probabilità di rilevazione vs soglia LET per ogni Vbias.
 *
 * TRE definizioni a confronto sullo stesso plot:
 *
 *   P_cross = N_crossing  / N_laser_tot
 *             Edoardo: "quante volte il SiPM ha superato la soglia LET?"
 *             (t_rise valido, prima di qualsiasi check sul TOT/t_fall)
 *
 *   P_tot   = N_accepted  / N_laser_tot
 *             "quanti eventi hanno sia t_rise che t_fall validi e TOT ok?"
 *             (quello che finisce in cache, penalizzato dal TOT killer)
 *
 *   P_sig   = N_sig_netto / N_laser_tot
 *             "segnale vero dopo sottrazione BKG (fit S+B su delta_t)"
 *
 * N_crossing, N_accepted, N_laser_tot sono letti dai TParameter salvati
 * in cache da VbiasAnalysis_v2.h (richiede la patch nCrossing).
 * Cache vecchie senza TParameter: P_cross/P_tot saltate con [WARN],
 * viene calcolata solo P_sig.
 *
 * Compile:  .L sipm_detection_prob.cpp+
 * Run:      sipm_detection_prob()
 */

#include "../header/Config.h"
#include "../header/CalibIO.h"
#include "../header/EventCache.h"
#include "../header/OutputManager.h"
#include "../header/SidebandAnalysis.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TMath.h>


// ══════════════════════════════════════════════════════════════════════════
//  Struct risultato per una singola (vbias, frac_pe)
// ══════════════════════════════════════════════════════════════════════════
struct ProbResult {
    double p_cross  = -1;  double ep_cross = 0;
    double p_tot    = -1;  double ep_tot   = 0;
    double p_sig    = -1;  double ep_sig   = 0;
    long   n_laser  = 0;
    long   n_cross  = -1;
    long   n_acc    = -1;
    double n_sig    = -1;
};

// ── Errore binomiale ──────────────────────────────────────────────────────
static double binomErr(double p, long N) {
    if (N <= 0 || p < 0 || p > 1) return 0;
    return std::sqrt(p * (1.0 - p) / (double)N);
}

// ── Legge TParameter<long> (o <int> per backward compat) ─────────────────
static long readLongParam(TFile* f, const char* name) {
    {auto* p=dynamic_cast<TParameter<long>*>(f->Get(name)); if(p) return p->GetVal();}
    {auto* p=dynamic_cast<TParameter<int> *>(f->Get(name)); if(p) return (long)p->GetVal();}
    return -1L;
}

// ══════════════════════════════════════════════════════════════════════════
//  N_sig via fit S+B (usa SidebandAnalysis.h)
// ══════════════════════════════════════════════════════════════════════════
static std::pair<double,double> computeNsigSB(const std::string& cachePath,
                                               double fit_lo, double fit_hi,
                                               const std::string& tag)
{
    std::vector<TOTEvent> events;
    {
        TFile* f = TFile::Open(cachePath.c_str(), "READ");
        if (!f || f->IsZombie()) { delete f; return {-1,-1}; }
        TTree* tree = (TTree*)f->Get("events");
        if (!tree || tree->GetEntries()==0) { f->Close(); delete f; return {-1,-1}; }

        Double_t tot, delta_t, amp_max; Int_t n_pe;
        tree->SetBranchAddress("tot",     &tot);
        tree->SetBranchAddress("delta_t", &delta_t);
        tree->SetBranchAddress("amp_max", &amp_max);
        tree->SetBranchAddress("n_pe",    &n_pe);

        Long64_t N = tree->GetEntries();
        events.reserve((size_t)N);
        for (Long64_t i=0; i<N; ++i) {
            tree->GetEntry(i);
            events.push_back({tot, delta_t, amp_max, (int)n_pe,
                              TRISE_INVALID, TRISE_INVALID});
        }
        f->Close(); delete f;
    }
    if (events.empty()) return {-1,-1};

    SBResult sb = fitSplusB(events, fit_lo, fit_hi,
                             SB::DEFAULT_SB_WIDTH, SB::DEFAULT_EXT, tag);
    if (!sb.fit_ok) {
        std::cerr << "  [Prob] S+B fit non convergito per " << tag << "\n";
        cleanupSBResult(sb); return {-1,-1};
    }

    double nsig  = sb.N_sig;
    double ensig = std::sqrt(std::max(0.0,nsig) + std::max(0.0,sb.N_bkg));
    std::cout << "  [Prob] S+B → N_sig=" << std::fixed << std::setprecision(0)
              << nsig << " ±" << ensig
              << "  N_bkg=" << sb.N_bkg
              << "  chi2/ndf=" << std::setprecision(2) << sb.chi2ndf << "\n";
    cleanupSBResult(sb);
    return {nsig, ensig};
}

// ══════════════════════════════════════════════════════════════════════════
//  Calcola ProbResult leggendo la cache
// ══════════════════════════════════════════════════════════════════════════
static ProbResult computeProb(const std::string& cachePath,
                               double fit_lo, double fit_hi,
                               const std::string& tag)
{
    ProbResult res;

    // ── Leggi TParameter dalla cache ──────────────────────────────────────
    {
        TFile* f = TFile::Open(cachePath.c_str(), "READ");
        if (!f || f->IsZombie()) { delete f; return res; }
        res.n_laser = readLongParam(f, "n_laser_tot");
        res.n_cross = readLongParam(f, "n_crossing");
        res.n_acc   = readLongParam(f, "n_accepted");

        // Fallback n_acc: se assente usa GetEntries del tree
        if (res.n_acc < 0) {
            TTree* t=(TTree*)f->Get("events");
            if (t) { res.n_acc=(long)t->GetEntries();
                std::cout<<"  [WARN] n_accepted assente (cache vecchia)"
                         <<" — uso GetEntries="<<res.n_acc<<"\n"; }
        }
        if (res.n_cross < 0)
            std::cout<<"  [WARN] n_crossing assente (cache vecchia)"
                     <<" — P_cross non disponibile.\n";
        if (res.n_laser <= 0)
            std::cout<<"  [WARN] n_laser_tot assente (cache vecchia)"
                     <<" — P_cross/P_tot non calcolabili.\n";
        f->Close(); delete f;
    }

    double nl = (res.n_laser > 0) ? (double)res.n_laser : -1.0;

    if (nl > 0) {
        if (res.n_cross >= 0) {
            res.p_cross  = std::min(1.0, (double)res.n_cross / nl);
            res.ep_cross = binomErr(res.p_cross, res.n_laser);
        }
        if (res.n_acc >= 0) {
            res.p_tot  = std::min(1.0, (double)res.n_acc / nl);
            res.ep_tot = binomErr(res.p_tot, res.n_laser);
        }
    }

    // ── N_sig via S+B ────────────────────────────────────────────────────
    auto [nsig, ensig] = computeNsigSB(cachePath, fit_lo, fit_hi, tag);
    if (nsig >= 0 && nl > 0) {
        res.n_sig  = nsig;
        res.p_sig  = std::min(1.0, std::max(0.0, nsig / nl));
        res.ep_sig = std::max(0.0, ensig / nl);
    }

    std::cout << "  [Prob] N_laser=" << res.n_laser
              << "  N_cross=" << res.n_cross
              << "  N_acc=" << res.n_acc
              << "  N_sig=" << std::fixed << std::setprecision(0) << res.n_sig << "\n"
              << "  [Prob] P_cross=" << std::setprecision(1)
              << (res.p_cross>=0?res.p_cross*100:-1) << "%"
              << "  P_tot=" << (res.p_tot>=0?res.p_tot*100:-1) << "%"
              << "  P_sig=" << (res.p_sig>=0?res.p_sig*100:-1) << "%\n";
    return res;
}

// ══════════════════════════════════════════════════════════════════════════
//  Trova cache per (vbias, frac, cal)
// ══════════════════════════════════════════════════════════════════════════
static std::string findCache(int vbias, double frac, const CalibResult& cal,
                              const std::string& dataDir, bool use_filter)
{
    std::string full = eventCachePath(vbias, frac, cal.cutoff_MHz,
                                       cal.laser_thr, dataDir, use_filter);
    std::string cpFile;
    { size_t sl=full.find_last_of("/\\");
      std::string base=(sl==std::string::npos)?full:full.substr(sl+1);
      cpFile = (base.size()>5 && base.substr(base.size()-5)==".root")
               ? base.substr(0,base.size()-5) : base; }

    std::string best_clean, best_any;
    void* dp=gSystem->OpenDirectory(dataDir.c_str());
    if (dp) { const char* ent;
        while((ent=gSystem->GetDirEntry(dp))!=nullptr){
            std::string fn(ent);
            if(fn.find(cpFile)!=0) continue;
            if(fn.size()<5||fn.substr(fn.size()-5)!=".root") continue;
            bool bkg=(fn.find("_bkgsub")!=std::string::npos);
            if(!bkg&&best_clean.empty()) best_clean=fn;
            if(best_any.empty()) best_any=fn;
        }gSystem->FreeDirectory(dp);}
    std::string chosen=best_clean.empty()?best_any:best_clean;
    return chosen.empty()?"":dataDir+"/"+chosen;
}

// ══════════════════════════════════════════════════════════════════════════
//  Costruisce TGraphErrors con stile
// ══════════════════════════════════════════════════════════════════════════
static TGraphErrors* makeGraph(const std::vector<double>& x,
                                const std::vector<double>& y,
                                const std::vector<double>& ey,
                                int color, int marker, int lStyle)
{
    std::vector<double> ex(x.size(),0.0);
    TGraphErrors* gr=new TGraphErrors((int)x.size(),
        x.data(),y.data(),ex.data(),ey.data());
    gr->SetMarkerStyle(marker); gr->SetMarkerSize(1.0);
    gr->SetMarkerColor(color);  gr->SetLineColor(color);
    gr->SetLineWidth(2);        gr->SetLineStyle(lStyle);
    return gr;
}

// ══════════════════════════════════════════════════════════════════════════
//  MAIN
// ══════════════════════════════════════════════════════════════════════════
void sipm_detection_prob() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // ── Selezione directory ───────────────────────────────────────────────
    std::string dataDir;
    {
        std::string root=DATA_DIR;
        size_t slash=root.find_last_of("/\\");
        if(slash!=std::string::npos) root=root.substr(0,slash);
        std::cout<<"\n+==========================================================+\n"
                 <<"|  SiPM DETECTION PROBABILITY vs LET                      |\n"
                 <<"+==========================================================+\n"
                 <<"  Root: "<<root<<"\n\n";

        std::vector<std::string> subs;
        void* dp=gSystem->OpenDirectory(root.c_str());
        if(dp){const char* ent;
            while((ent=gSystem->GetDirEntry(dp))!=nullptr){
                std::string s(ent);
                if(s=="."||s=="..") continue;
                FileStat_t st;
                if(gSystem->GetPathInfo((root+"/"+s).c_str(),st)==0&&R_ISDIR(st.fMode))
                    subs.push_back(s);
            }gSystem->FreeDirectory(dp);}
        std::sort(subs.begin(),subs.end());
        if(subs.empty()){std::cerr<<"[ERR] Nessuna sottocartella.\n";return;}
        for(size_t i=0;i<subs.size();++i) std::cout<<"  ["<<(i+1)<<"] "<<subs[i]<<"\n";
        std::cout<<"\n  Scegli [n] o path: "<<std::flush;
        std::string line; if(!std::getline(std::cin,line)) return;
        auto b=line.find_first_not_of(" \t\r\n"), e=line.find_last_not_of(" \t\r\n");
        if(b==std::string::npos) return;
        line=line.substr(b,e-b+1);
        try{size_t pos=0;int idx=std::stoi(line,&pos);
            if(pos==line.size()&&idx>=1&&idx<=(int)subs.size())
                dataDir=root+"/"+subs[idx-1];
            else dataDir=line;}catch(...){dataDir=line;}
        if(gSystem->AccessPathName(dataDir.c_str())){std::cerr<<"[ERR] Non accessibile.\n";return;}
        g_data_dir_override=dataDir;
        std::cout<<"  --> "<<dataDir<<"\n\n";
    }

    // ── Helper ────────────────────────────────────────────────────────────
    auto readLine=[](const std::string& p)->std::string{
        std::string l;std::cout<<p<<std::flush;std::getline(std::cin,l);
        auto b=l.find_first_not_of(" \t\r\n");
        return(b==std::string::npos)?std::string(""):l.substr(b);};
    auto readDouble=[&readLine](const std::string& p)->double{
        while(true){std::string l=readLine(p);try{return std::stod(l);}catch(...){}}};
    auto buildRange=[](double lo,double hi,double step)->std::vector<double>{
        std::vector<double> o;int n=(int)std::round((hi-lo)/step);
        for(int i=0;i<=n;++i){double v=lo+i*step;if(v<=hi+1e-9)o.push_back(v);}return o;};

    // ── Input ─────────────────────────────────────────────────────────────
    std::vector<int> vbiasList;
    {std::string l=readLine("Vbias lista (es. 53 54 55): ");
     std::stringstream ss(l);int v;while(ss>>v)vbiasList.push_back(v);}
    if(vbiasList.empty()){std::cerr<<"[ERR] Nessun Vbias.\n";return;}

    std::vector<double> fracList;
    {std::cout<<"\nLET frac_pe — range (start end step) o lista:\n";
     std::string l=readLine("> ");
     std::vector<double> nums; std::stringstream ss(l); double x;
     while(ss>>x) nums.push_back(x);
     bool isRange=(nums.size()==3&&nums[0]>0&&nums[1]>nums[0]
                   &&nums[2]>0&&nums[2]<=(nums[1]-nums[0])+1e-9);
     if(isRange) fracList=buildRange(nums[0],nums[1],nums[2]);
     else for(double v:nums) if(v>0) fracList.push_back(v);
     if(fracList.empty()) fracList.push_back(1.0);
     std::cout<<"  → "<<fracList.size()<<" soglie\n";}

    double cutoff_global=0;
    {std::string l=readLine("Cutoff MHz [0=auto]: ");
     if(!l.empty())try{cutoff_global=std::stod(l);}catch(...){};}
    double fit_lo=readDouble("Fit window start [ns]: ");
    double fit_hi=readDouble("Fit window end   [ns]: ");
    bool use_filter=true;
    {std::string l=readLine("LP filter? [y/n]: ");
     if(!l.empty()&&std::tolower((unsigned char)l[0])=='n')use_filter=false;}

    // ── Calibrazioni ──────────────────────────────────────────────────────
    std::map<int,CalibResult> calMap;
    for(int vb:vbiasList){
        double cut=cutoff_global;
        if(cut==0){
            std::vector<double> fc;
            std::string pre="calib_vbias"+std::to_string(vb)+"_cut",suf="mhz.root";
            void* dp=gSystem->OpenDirectory(dataDir.c_str());
            if(dp){const char* ent;
                while((ent=gSystem->GetDirEntry(dp))!=nullptr){
                    std::string fn(ent);
                    if(fn.size()<pre.size()+suf.size())continue;
                    if(fn.substr(0,pre.size())!=pre)continue;
                    if(fn.substr(fn.size()-suf.size())!=suf)continue;
                    try{fc.push_back(std::stod(fn.substr(pre.size(),
                        fn.size()-pre.size()-suf.size())));}catch(...){}
                }gSystem->FreeDirectory(dp);}
            std::sort(fc.begin(),fc.end());
            if(fc.empty()){std::cerr<<"  [WARN] Vbias="<<vb<<": no cal.\n";continue;}
            cut=fc.back();
        }
        CalibResult cal;
        if(!loadCalibration(cal,vb,cut,dataDir)){
            std::cerr<<"  [WARN] Vbias="<<vb<<": cal non caricata.\n";continue;}
        calMap[vb]=cal;
        std::cout<<"  Vbias="<<vb<<"  gain="<<cal.m<<" mV/pe  cut="<<cal.cutoff_MHz<<" MHz\n";
    }
    if(calMap.empty()){std::cerr<<"[ERR] Nessuna calibrazione.\n";return;}

    // ── Loop vbias × frac_pe ─────────────────────────────────────────────
    std::map<int,std::map<double,ProbResult>> probMap;
    for(auto& [vb,cal]:calMap){
        for(double frac:fracList){
            std::cout<<"\n  Vbias="<<vb<<"  frac_pe="<<frac<<"\n";
            std::string cp=findCache(vb,frac,cal,dataDir,use_filter);
            if(cp.empty()){std::cerr<<"  [WARN] Cache non trovata.\n";continue;}
            std::cout<<"  Cache: "<<cp.substr(cp.find_last_of("/\\")+1)<<"\n";
            probMap[vb][frac]=computeProb(cp,fit_lo,fit_hi,
                                          Form("vbias%d_let%.2fpe",vb,frac));
        }
    }
    if(probMap.empty()){std::cerr<<"[ERR] Nessun dato.\n";return;}

    // ══════════════════════════════════════════════════════════════════════
    //  PLOT — 3 curve per Vbias, colori distinti, stile linea per tipo
    // ══════════════════════════════════════════════════════════════════════
    OutCtx ctx=createOutputDirs("det_prob");

    static const int COLS[]={kAzure+1,kRed+1,kGreen+2,kOrange+7,kMagenta+1,kCyan+2};
    const int NC=(int)(sizeof(COLS)/sizeof(COLS[0]));

    TCanvas* c=new TCanvas("cDetProb","P_{det} vs LET",1100,700);
    c->SetGrid();
    c->SetLeftMargin(PAD_LEFT); c->SetRightMargin(PAD_RIGHT);
    c->SetBottomMargin(PAD_BOTTOM); c->SetTopMargin(PAD_TOP);

    // Due legende: una per Vbias (colore), una per tipo P (stile)
    int nVb=(int)probMap.size();
    TLegend* legVb=new TLegend(0.55,0.68,0.77,0.68+0.055*nVb);
    legVb->SetBorderSize(1);legVb->SetFillStyle(1001);legVb->SetFillColor(0);
    legVb->SetTextFont(42);legVb->SetTextSize(0.032);
    legVb->SetHeader("Vbias","C");

    TLegend* legType=new TLegend(0.78,0.68,0.97,0.68+0.055*3);
    legType->SetBorderSize(1);legType->SetFillStyle(1001);legType->SetFillColor(0);
    legType->SetTextFont(42);legType->SetTextSize(0.032);
    legType->SetHeader("Definizione P_{det}","C");

    // Entries dummy per legType
    // FIX: store pointers — these graphs are NOT drawn on canvas so ROOT
    // does not take ownership; delete them explicitly after the plot is saved.
    std::vector<TGraphErrors*> dummyGraphs;
    {auto mk=[&dummyGraphs](int mk,int ls)->TGraphErrors*{
        TGraphErrors*g=new TGraphErrors(1);
        g->SetMarkerStyle(mk);g->SetLineStyle(ls);
        g->SetMarkerColor(kGray+2);g->SetLineColor(kGray+2);
        g->SetLineWidth(2);
        dummyGraphs.push_back(g);
        return g;};
     legType->AddEntry(mk(20,1),"P_{cross} = N_{cross}/N_{laser}  (solo t_{rise})","lp");
     legType->AddEntry(mk(24,2),"P_{tot}   = N_{acc}/N_{laser}    (TOT ok)","lp");
     legType->AddEntry(mk(22,3),"P_{sig}   = N_{sig}/N_{laser}    (fit S+B)","lp");}

    bool first=true;
    double ymax=5;
    int ci=0;

    for(auto& [vb,letmap]:probMap){
        int col=COLS[ci%NC];

        // Raccogli punti per i tre tipi
        struct Pts { std::vector<double> x,y,ey; };
        Pts pc,pt,ps;

        for(auto& [frac,pr]:letmap){
            double xv=frac*100.0;
            if(pr.p_cross>=0){pc.x.push_back(xv);pc.y.push_back(pr.p_cross*100);pc.ey.push_back(pr.ep_cross*100);ymax=std::max(ymax,pr.p_cross*100);}
            if(pr.p_tot  >=0){pt.x.push_back(xv);pt.y.push_back(pr.p_tot  *100);pt.ey.push_back(pr.ep_tot  *100);ymax=std::max(ymax,pr.p_tot  *100);}
            if(pr.p_sig  >=0){ps.x.push_back(xv);ps.y.push_back(pr.p_sig  *100);ps.ey.push_back(pr.ep_sig  *100);ymax=std::max(ymax,pr.p_sig  *100);}
        }

        auto draw=[&](Pts& p,int marker,int ls){
            if(p.x.empty()) return;
            TGraphErrors* gr=makeGraph(p.x,p.y,p.ey,col,marker,ls);
            gr->SetTitle(";Soglia LET (% p.e.);P_{det} (%)");
            if(first){
                gr->Draw("APL");
                gr->GetYaxis()->SetRangeUser(0,std::min(110.0,ymax*1.2+5));
                gr->GetXaxis()->SetTitleSize(0.048f);gr->GetYaxis()->SetTitleSize(0.048f);
                gr->GetXaxis()->SetLabelSize(0.040f);gr->GetYaxis()->SetLabelSize(0.040f);
                first=false;
            } else { gr->Draw("PL SAME"); }
        };
        draw(pc,20,1);
        draw(pt,24,2);
        draw(ps,22,3);

        // Entry Vbias nella legenda colore
        // FIX: gd not drawn on canvas → ROOT doesn't own it → must delete manually.
        TGraphErrors* gd=new TGraphErrors(1);
        gd->SetLineColor(col);gd->SetMarkerColor(col);
        gd->SetMarkerStyle(20);gd->SetLineWidth(2);
        legVb->AddEntry(gd,Form("Vbias = %d V",vb),"lp");
        dummyGraphs.push_back(gd);
        ++ci;
    }

    // Linea tratteggiata al 50%
    if(!fracList.empty()){
        double x0=fracList.front()*100-2, x1=fracList.back()*100+2;
        TLine* l=new TLine(x0,50,x1,50);
        l->SetLineColor(kGray+1);l->SetLineStyle(2);l->SetLineWidth(1);l->Draw();
    }

    legVb->Draw(); legType->Draw();

    // Info box
    {TPaveText* pt=new TPaveText(0.14,0.72,0.40,0.88,"NDC");
     pt->SetFillColorAlpha(kWhite,0.85);pt->SetBorderSize(1);
     pt->SetTextFont(42);pt->SetTextSize(0.029);pt->SetTextAlign(12);
     pt->AddText("Fit #Deltat window:");
     pt->AddText(Form("  [%.1f, %.1f] ns",fit_lo,fit_hi));
     pt->AddText(Form("  LP filter: %s",use_filter?"ON":"OFF"));
     pt->Draw();}

    c->Update(); c->Modified();
    ctx.savePNG(c,"detection_probability_vs_LET.png");

    // FIX: delete dummy graphs (not owned by any canvas)
    for (auto* g : dummyGraphs) delete g;
    dummyGraphs.clear();

    // ── Tabella riassuntiva ───────────────────────────────────────────────
    std::cout<<"\n+-- RIASSUNTO --\n"<<std::fixed<<std::setprecision(1);
    std::cout<<std::setw(7)<<"Vbias"
             <<std::setw(10)<<"LET(%pe)"
             <<std::setw(13)<<"P_cross(%)"
             <<std::setw(13)<<"P_tot(%)"
             <<std::setw(13)<<"P_sig(%)"<<"\n"
             <<"-----------------------------------------------------------\n";
    for(auto& [vb,letmap]:probMap)
        for(auto& [frac,pr]:letmap)
            std::cout<<std::setw(7)<<vb
                     <<std::setw(10)<<frac*100
                     <<std::setw(13)<<(pr.p_cross>=0?pr.p_cross*100:-1)
                     <<std::setw(13)<<(pr.p_tot  >=0?pr.p_tot  *100:-1)
                     <<std::setw(13)<<(pr.p_sig  >=0?pr.p_sig  *100:-1)<<"\n";

    std::cout<<"\n+==========================================================+\n"
             <<"|  DONE — sipm_detection_prob                              |\n"
             <<"+==========================================================+\n";
}
