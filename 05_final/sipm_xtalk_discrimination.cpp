//std::string pat = Form("calib_vbias%d_cut", vbias);
#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/CalibIO.h"
#include "../header/EventCache.h"
#include "../header/SignalProcessing.h"
#include "../header/ButterworthFilter.h"
#include "../header/TOTAnalysis.h"
#include "../header/Calibration.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cctype>
#include <limits>
#include <cstdlib>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>

#ifndef G_ANALYSIS_MODE_DEFINED
inline int g_analysis_mode = 0;
#define G_ANALYSIS_MODE_DEFINED
#endif

static constexpr int    N_SAMPLES = 1024;
static constexpr double AMP_LO_PE = 1.5;
static constexpr double AMP_HI_PE = 2.5;
static constexpr int    SMOOTH_K  = 3;

static constexpr double M_RANGE_MIN = 0.25;
static constexpr double M_RANGE_MAX = 1.25;
static constexpr int    M_BINS      = 200;

static std::string resolveHighDir(const std::string& arg_high)
{
    // FIX (v4 2026-04-29): il fallback hardcoded "../../data/data_high" non
    // corrisponde alla struttura corrente. Config.h definisce DATA_DIR_HIGH
    // (= data_filter_3 = high light, 1 run) come sorgente canonica.
    // Priorità: arg esplicito > env var SIPM_DATA_DIR_HIGH > DATA_DIR_HIGH.
    if (!arg_high.empty()) return arg_high;
    const char* env = std::getenv("SIPM_DATA_DIR_HIGH");
    if (env && env[0] != '\0') return std::string(env);
    return DATA_DIR_HIGH;
}

struct MEvent   { double tot; double m; double amp_pe; };
struct DvdtEvent { double rise_time; double peak2_ratio; double amp_pe; };

static void extractMEvents(const std::string& cachePath,
                            double pe_lo, double pe_hi,
                            double gain_mV_pe,
                            std::vector<MEvent>& out, long max_ev=200000)
{
    if (cachePath.empty()) return;
    
    TFile* f=TFile::Open(cachePath.c_str(),"READ");
    if (!f||f->IsZombie()){delete f;return;}
    
    TTree* tr=(TTree*)f->Get("events");
    if (!tr||!tr->GetBranch("delta_t_lo")||!tr->GetBranch("delta_t_hi")){
        std::cerr<<"[ERROR] Cache senza delta_t_lo/hi. Rigenera con sipm_tot_analysis v3.\n";
        f->Close();delete f;return;
    }

    Double_t tot,amp_max,delta_t_lo,delta_t_hi; Int_t n_pe;
    tr->SetBranchAddress("tot",        &tot);
    tr->SetBranchAddress("amp_max",    &amp_max);
    tr->SetBranchAddress("n_pe",       &n_pe);
    tr->SetBranchAddress("delta_t_lo", &delta_t_lo);
    tr->SetBranchAddress("delta_t_hi", &delta_t_hi);

    if (gain_mV_pe <= 0) {
        auto* pg=dynamic_cast<TParameter<double>*>(f->Get("gain_mV_pe"));
        if (pg) gain_mV_pe=pg->GetVal();
        if (gain_mV_pe <= 0) gain_mV_pe = 30.0;
    }

    // FIX DOCUMENTAZIONE: delta_frac DEVE corrispondere alla differenza delle
    // soglie usate in VbiasAnalysis_v2.h per calcolare delta_t_lo e delta_t_hi.
    // VbiasAnalysis_v2.h linee 140-141:
    //   let_thr_lo = 0.20 * let_thr  → delta_t_lo
    //   let_thr_hi = 0.80 * let_thr  → delta_t_hi
    // Quindi delta_frac = 0.80 - 0.20 = 0.60. Corretto.
    // Se mai cambi quelle costanti in VbiasAnalysis, cambia anche questo.
    const double delta_frac = 0.6;
    Long64_t nev=tr->GetEntries();
    out.reserve(std::min((long)nev,max_ev));

    for (Long64_t i=0; i<nev && (long)out.size()<max_ev; ++i) {
        tr->GetEntry(i);
        if (delta_t_lo < TRISE_INVALID+1.0) continue;
        if (delta_t_hi < TRISE_INVALID+1.0) continue;
        double amp_pe = (gain_mV_pe>0) ? amp_max/gain_mV_pe : (double)n_pe;
        if (amp_pe<AMP_LO_PE||amp_pe>AMP_HI_PE) continue;
        double m = (delta_t_hi - delta_t_lo) / delta_frac;
        out.push_back({tot,m,amp_pe});
    }
    f->Close();delete f;
    std::cout<<"  [m] extracted "<<out.size()<<" events from "<<cachePath<<"\n";
}

static std::vector<double> boxSmooth(const std::vector<double>& v, int k)
{
    int n=(int)v.size();
    std::vector<double> s(n,0.0);
    for (int i=0;i<n;++i){
        double sum=0; int cnt=0;
        for (int d=-k;d<=k;++d){int j=i+d;if(j>=0&&j<n){sum+=v[j];++cnt;}}
        s[i]=cnt>0?sum/cnt:0.0;
    }
    return s;
}

static DvdtEvent computeDvdt(const std::vector<double>& t,
                              const std::vector<double>& a,
                              double t_win_start, double t_win_end,
                              double gain_mV_pe)
{
    DvdtEvent f; f.rise_time=-1;
    int n=(int)t.size();
    if (n<10||gain_mV_pe<=0) return f;

    int j0=0,j1=n-1;
    for (int j=0;j<n;++j)  if(t[j]>=t_win_start-5.0){j0=j;break;}
    for (int j=n-1;j>=0;--j) if(t[j]<=t_win_end+5.0){j1=j;break;}
    if (j1-j0<5) return f;

    int j_peak=j0;
    for (int j=j0;j<=j1;++j) if(a[j]>a[j_peak]) j_peak=j;
    double amp_max=a[j_peak];
    if (amp_max<AMP_LO_PE*gain_mV_pe||amp_max>AMP_HI_PE*gain_mV_pe) return f;
    f.amp_pe=amp_max/gain_mV_pe;

    double thr10=0.10*amp_max, thr90=0.90*amp_max;
    double t10=-1,t90=-1;
    for (int j=j0;j<j_peak&&(t10<0||t90<0);++j){
        if(t10<0&&a[j]<thr10&&a[j+1]>=thr10){double da=a[j+1]-a[j];
            t10=(std::abs(da)>1e-12)?t[j]+(thr10-a[j])*(t[j+1]-t[j])/da:t[j];}
        if(t90<0&&a[j]<thr90&&a[j+1]>=thr90){double da=a[j+1]-a[j];
            t90=(std::abs(da)>1e-12)?t[j]+(thr90-a[j])*(t[j+1]-t[j])/da:t[j];}
    }
    if (t10<0||t90<0||t90<=t10) return f;
    f.rise_time=t90-t10;

    int j_start = j0;
    int j_end = std::min(j_peak + (int)((j_peak-j0)*0.2), j1);
    int nw=j_end-j_start+1;
    if (nw<5) return f;

    std::vector<double> der(nw,0.0);
    for (int i=1;i<nw-1;++i){int jj=j_start+i;double dt=t[jj+1]-t[jj-1];
        if(dt>1e-9)der[i]=(a[jj+1]-a[jj-1])/dt;}
    {double dt0=t[j_start+1]-t[j_start],dte=t[j_end]-t[j_end-1];
     der[0]=(dt0>1e-9)?(a[j_start+1]-a[j_start])/dt0:0.0;
     der[nw-1]=(dte>1e-9)?(a[j_end]-a[j_end-1])/dte:0.0;}

    auto sder=boxSmooth(der,SMOOTH_K);

    int i_pk1=0;
    for (int i=1;i<nw;++i) if(sder[i]>sder[i_pk1]) i_pk1=i;
    double pk1=sder[i_pk1];
    if (pk1<=0) return f;

    bool in_val=false; 
    double pk2_after=0.0;
    for (int i=i_pk1+1;i<nw-1;++i){
        if(!in_val && sder[i]<0.05*pk1) in_val=true;
        if(in_val && sder[i]>sder[i-1] && sder[i]>sder[i+1] && sder[i]>pk2_after) 
            pk2_after=sder[i];
    }
    
    double pk2_before=0.0;
    for (int i=1; i<i_pk1-2; ++i){
        if(sder[i]>sder[i-1] && sder[i]>sder[i+1] && sder[i]>pk2_before && sder[i]>0.15*pk1)
            pk2_before=sder[i];
    }
    
    double slope_mean=0, slope_var=0;
    int n_slope=0;
    for(int i=0; i<nw && i<i_pk1; ++i) {
        if(sder[i] > 0.3*pk1) {
            slope_mean += sder[i];
            n_slope++;
        }
    }
    if(n_slope>2) {
        slope_mean /= n_slope;
        for(int i=0; i<nw && i<i_pk1; ++i) {
            if(sder[i] > 0.3*pk1) {
                double d = sder[i] - slope_mean;
                slope_var += d*d;
            }
        }
        slope_var = std::sqrt(slope_var/n_slope) / slope_mean;
    }
    
    double pk2 = std::max(pk2_after, pk2_before);
    
    if(pk2 < 0.15*pk1 && slope_var > 0) {
        f.peak2_ratio = slope_var;
    } else {
        f.peak2_ratio = pk2/pk1;
    }
    
    return f;
}

static std::map<int,std::string> findRunFiles(int vbias, const std::string& dataDir)
{
    std::map<int,std::string> runs;
    // FIX: prova entrambi i pattern. La cartella high light usa filename
    // SENZA graffe ("data.vbias_55_run_3.root") mentre la cartella low light
    // usa il pattern con graffe ("data.vbias_{55}_run_3.root").
    const std::string patBrace = makeRunPattern(vbias);                          // con graffe
    const std::string patPlain = RUN_PREFIX + std::to_string(vbias) + "_run_";   // senza graffe
    void* dirp=gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp){
        std::cerr<<"[WARN] Cannot open dir: "<<dataDir<<"\n"; return runs;
    }
    const char* entry;
    while((entry=gSystem->GetDirEntry(dirp))!=nullptr){
        std::string fn(entry);
        // Match con graffe o senza graffe.
        bool match = (fn.find(patBrace) != std::string::npos);
        if (!match) {
            // Per il pattern senza graffe, dobbiamo evitare match parziali:
            // "data.vbias_55_run_" matcha anche "data.vbias_555_run_". Verifichiamo
            // che subito dopo RUN_PREFIX ci sia il numero esatto seguito da "_run_".
            size_t p = fn.find(patPlain);
            if (p != std::string::npos) {
                // Verifica boundary: prima di patPlain c'e' inizio stringa o '/'
                // e dopo patPlain c'e' una cifra (numero del run).
                bool start_ok = (p == 0);
                bool end_ok = (p + patPlain.size() < fn.size())
                              && std::isdigit((unsigned char)fn[p+patPlain.size()]);
                if (start_ok && end_ok) match = true;
            }
        }
        if (!match) continue;
        if(fn.size()<5||fn.substr(fn.size()-5)!=".root") continue;
        size_t pos=fn.find("_run_");
        if(pos==std::string::npos) continue;
        try{
            std::string sub=fn.substr(pos+5);
            size_t dot=sub.find(".root");
            if(dot!=std::string::npos) sub=sub.substr(0,dot);
            runs[std::stoi(sub)]=dataDir+"/"+fn;
        }catch(...){}
    }
    gSystem->FreeDirectory(dirp);
    return runs;
}

static void extractDvdt(const std::map<int,std::string>& runs,
                         double gain_mV_pe, double cutoff_MHz, double fs_MHz,
                         double t_win_start, double t_win_end,
                         const std::string& label,
                         std::vector<DvdtEvent>& out, long max_ev=1000000)
{
    const double bl_rms_max=(cutoff_MHz>0)?2.0:4.0;
    long total=0;
    for (auto& [runN,path]:runs){
        if(total>=max_ev) break;
        TFile* fIn=TFile::Open(path.c_str(),"READ");
        if(!fIn||fIn->IsZombie()){delete fIn;continue;}
        TTree* tr=(TTree*)fIn->Get("ch1");
        if(!tr){fIn->Close();delete fIn;continue;}

        Double_t t_arr[N_SAMPLES],a_arr[N_SAMPLES];
        tr->SetBranchAddress("time",     t_arr);
        tr->SetBranchAddress("amplitude",a_arr);
        tr->SetCacheSize(2*1024*1024);

        Long64_t nev=tr->GetEntries();
        std::cout<<"  ["<<label<<" dV/dt] run "<<runN<<"  "<<nev<<" events\n";

        for(Long64_t i=0;i<nev&&total<max_ev;++i){
            tr->GetEntry(i);
            std::vector<double> tv(t_arr,t_arr+N_SAMPLES);
            std::vector<double> av(a_arr,a_arr+N_SAMPLES);

            bool bl_ok=true;
            auto ac=correctBaseline(tv,av,BASELINE_START,BASELINE_END,bl_rms_max,&bl_ok);
            if(!bl_ok) continue;

            std::vector<double> af=ac;
            if(cutoff_MHz>0&&fs_MHz>0)
                af=butterworthLowPass(ac,fs_MHz,cutoff_MHz,4);

            DvdtEvent ev=computeDvdt(tv,af,t_win_start,t_win_end,gain_mV_pe);
            if(ev.rise_time<0) continue;
            out.push_back(ev);
            ++total;
        }
        fIn->Close();delete fIn;
    }
    std::cout<<"  ["<<label<<" dV/dt] selected: "<<out.size()<<"\n";
}

static TH1D* makeH1(const char* nm,const char* tt,
                    int nb,double lo,double hi,int col, int fillStyle=0, int lstyle=1)
{
    auto* h=new TH1D(nm,tt,nb,lo,hi);
    h->SetDirectory(nullptr);
    h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(lstyle);
    if(fillStyle>0){h->SetFillStyle(fillStyle);h->SetFillColorAlpha(col,0.35f);}
    h->GetXaxis()->SetTitleSize(0.045f);h->GetYaxis()->SetTitleSize(0.045f);
    h->GetXaxis()->SetLabelSize(0.040f);h->GetYaxis()->SetLabelSize(0.040f);
    h->GetXaxis()->SetTitleFont(42);h->GetYaxis()->SetTitleFont(42);
    return h;
}
static void normArea(TH1D* h){if(h->Integral()>0)h->Scale(1.0/h->Integral());}

// Risultato del fit doppia gaussiana sulla distribuzione di m
struct DG2Result {
    bool   ok = false;
    double mu1=0, sigma1=0, mu1Err=0, sigma1Err=0;  // peak 1 (lower)
    double mu2=0, sigma2=0, mu2Err=0, sigma2Err=0;  // peak 2 (higher)
    double A1=0, A2=0;
    int    N=0;
    TF1*   f_total = nullptr;  // fit totale (drawn esternamente)
    TF1*   f_g1    = nullptr;  // gaussiana 1 (drawn esternamente)
    TF1*   f_g2    = nullptr;  // gaussiana 2 (drawn esternamente)
};

// Double Gaussian fit della distribuzione m. Funziona sia per low che per high
// light. Restituisce i parametri estratti e i TF1 da disegnare. Il chiamante
// decide cosa disegnare e prende possesso dei TF1 (deve deletare a fine vita).
//
// FIX: prima la funzione era pensata solo per high light (suffisso fisso "f2gauss"
// causava collisioni di nome se chiamata due volte). Ora il suffisso e'
// passato per parametro.
//
// Per low light: ci si aspetta UN picco principale con eventuale spalla di
// crosstalk. Il fit a doppia gaussiana e' comunque utile per quantificare
// la coda. Se peak2 ha ampiezza < 5% di peak1, il fit non e' affidabile.
static DG2Result fitDoubleGaussian(TH1D* h, const std::string& suffix,
                                    int color1=kGreen+2, int color2=kMagenta+2)
{
    DG2Result res;
    if (!h || h->GetEntries() < 100) return res;
    res.N = (int)h->GetEntries();

    // Find two peaks
    int bin1 = h->GetMaximumBin();
    double peak1_pos = h->GetBinCenter(bin1);
    double peak1_val = h->GetBinContent(bin1);

    // Suppress first peak to find second
    TH1D* hTemp = (TH1D*)h->Clone(Form("hTemp_%s", suffix.c_str()));
    hTemp->SetDirectory(nullptr);
    for (int i = bin1-5; i <= bin1+5; ++i) {
        if (i > 0 && i <= h->GetNbinsX()) hTemp->SetBinContent(i, 0);
    }
    int bin2 = hTemp->GetMaximumBin();
    double peak2_pos = hTemp->GetBinCenter(bin2);
    double peak2_val = hTemp->GetBinContent(bin2);
    delete hTemp;

    // Ensure peak1 < peak2 (peak1 e' SEMPRE quello a m piu' basso)
    if (peak1_pos > peak2_pos) {
        std::swap(peak1_pos, peak2_pos);
        std::swap(peak1_val, peak2_val);
    }

    // Define double Gaussian con suffisso unico per evitare collisione nomi
    TF1* f = new TF1(Form("f2gauss_%s", suffix.c_str()),
        "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2)",
        h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

    f->SetParameters(peak1_val, peak1_pos, 0.02, peak2_val, peak2_pos, 0.02);
    f->SetParLimits(1, peak1_pos-0.05, peak1_pos+0.05);
    f->SetParLimits(4, peak2_pos-0.05, peak2_pos+0.05);
    f->SetParLimits(2, 0.005, 0.1);
    f->SetParLimits(5, 0.005, 0.1);
    f->SetParLimits(0, 0, peak1_val*5);
    f->SetParLimits(3, 0, peak2_val*5);

    int status = h->Fit(f, "RQN");
    res.ok = (status == 0 || status == 4000);

    // Singole gaussiane per visualizzazione
    TF1* g1 = new TF1(Form("g1_%s", suffix.c_str()), "gaus",
                       h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    g1->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
    g1->SetLineColor(color1); g1->SetLineStyle(2); g1->SetLineWidth(2);

    TF1* g2 = new TF1(Form("g2_%s", suffix.c_str()), "gaus",
                       h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    g2->SetParameters(f->GetParameter(3), f->GetParameter(4), f->GetParameter(5));
    g2->SetLineColor(color2); g2->SetLineStyle(2); g2->SetLineWidth(2);

    f->SetLineColor(kBlack); f->SetLineWidth(2);

    // Estrai parametri
    res.A1     = f->GetParameter(0);
    res.mu1    = f->GetParameter(1);
    res.mu1Err = f->GetParError(1);
    res.sigma1 = std::abs(f->GetParameter(2));
    res.sigma1Err = f->GetParError(2);
    res.A2     = f->GetParameter(3);
    res.mu2    = f->GetParameter(4);
    res.mu2Err = f->GetParError(4);
    res.sigma2 = std::abs(f->GetParameter(5));
    res.sigma2Err = f->GetParError(5);
    res.f_total = f;
    res.f_g1    = g1;
    res.f_g2    = g2;

    if (res.ok) {
        std::cout << "  [Fit-2G " << suffix << "] Peak 1: mu="
                  << std::fixed << std::setprecision(4) << res.mu1
                  << " +/- " << res.mu1Err
                  << "  sigma=" << res.sigma1 << " +/- " << res.sigma1Err << "\n";
        std::cout << "  [Fit-2G " << suffix << "] Peak 2: mu="
                  << res.mu2 << " +/- " << res.mu2Err
                  << "  sigma=" << res.sigma2 << " +/- " << res.sigma2Err << "\n";
    } else {
        std::cout << "  [Fit-2G " << suffix << "] fit non convergito (status="
                  << status << ")\n";
    }

    return res;
}

static double fisherCut(const std::vector<double>& v1,
                        const std::vector<double>& v2,
                        double& sep_sigma)
{
    sep_sigma=0;
    if(v1.empty()||v2.empty()) return 0;
    auto mean=[](const std::vector<double>& v)->double{
        return std::accumulate(v.begin(),v.end(),0.0)/v.size();};
    auto var=[&](const std::vector<double>& v)->double{
        double m=mean(v),s2=0;for(auto x:v)s2+=(x-m)*(x-m);return s2/v.size();};
    double m1=mean(v1),m2=mean(v2),s1=var(v1),s2=var(v2);
    if(s1+s2<=0) return 0.5*(m1+m2);
    double w1=(double)v1.size()/(v1.size()+v2.size());
    double w2=(double)v2.size()/(v1.size()+v2.size());
    double cut=(w1*m1+w2*m2);
    double pooled=std::sqrt(w1*s1+w2*s2);
    if(pooled>0) sep_sigma=std::abs(m2-m1)/pooled;
    return cut;
}

static std::vector<double> scanAvailableFilters(int vbias, const std::string& dataDir)
{
    std::vector<double> cutoffs;
    std::string pat = Form("calib_vbias%d_cut", vbias);
    
    std::cout << "[DEBUG] Scanning directory: " << dataDir << "\n";
    std::cout << "[DEBUG] Looking for pattern: " << pat << "\n";
    
    void* dirp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dirp) {
        std::cout << "[DEBUG] Failed to open directory!\n";
        return cutoffs;
    }
    
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr) {
        checkEsc(); // Hook to stop process on ESC key
        
        std::string fn(entry);
        if (fn == "." || fn == "..") continue;
        
        if (fn.find(pat) != std::string::npos) {
            std::cout << "[DEBUG] Found potential match: " << fn << "\n";
        } else {
            continue;
        }
        
        if (fn.size() < 5 || fn.substr(fn.size()-5) != ".root") {
            std::cout << "[DEBUG] Rejected (not .root): " << fn << "\n";
            continue;
        }
        
        size_t pos = fn.find("_cut");
        if (pos == std::string::npos) continue;
        
        size_t pos2 = fn.find("mhz", pos);
        if (pos2 == std::string::npos) {
            std::cout << "[DEBUG] Rejected ('mhz' not found): " << fn << "\n";
            continue;
        }
        
        try {
            std::string sub = fn.substr(pos+4, pos2-pos-4);
            std::cout << "[DEBUG] Extracted cutoff string: '" << sub << "'\n";
            
            double cut = std::stod(sub);
            std::cout << "[DEBUG] Parsed cutoff value: " << cut << " MHz\n";
            
            if (std::find(cutoffs.begin(), cutoffs.end(), cut) == cutoffs.end()) {
                cutoffs.push_back(cut);
            }
        } catch(...) {
            std::cout << "[DEBUG] Exception during std::stod parsing!\n";
        }
    }
    gSystem->FreeDirectory(dirp);
    std::sort(cutoffs.begin(), cutoffs.end());
    return cutoffs;
}

void analyzeVbias(int vbias, double thr_opt, double cutoff_MHz,
                  double fit_lo, double fit_hi)
{
    const std::string dataDirLow  = DATA_DIR;
    const std::string dataDirHigh = resolveHighDir("");
    
    std::cout<<"\n═══════════════════════════════════════════════════\n"
             <<"  CROSSTALK vs 2PE DISCRIMINATION\n"
             <<"═══════════════════════════════════════════════════\n"
             <<"  Low light  → "<<dataDirLow<<"\n"
             <<"  High light → "<<dataDirHigh<<"\n\n";

    CalibResult cal;
    if (!loadCalibration(cal,vbias,cutoff_MHz,dataDirLow)||!cal.ok){
        std::cerr<<"[ERROR] No calibration found for vbias="<<vbias
                 <<" cutoff="<<cutoff_MHz<<" MHz\n";
        return;
    }

    const double gain_mV_pe  = cal.m;
    const double t_win_start = cal.t_trig_start;
    const double t_win_end   = cal.t_trig_end;
    const double pe_lo = AMP_LO_PE, pe_hi = AMP_HI_PE;

    double fs_MHz=5000.0;
    {
        auto runs=findRunFiles(vbias,dataDirLow);
        if (!runs.empty()){
            TFile* f0=TFile::Open(runs.begin()->second.c_str(),"READ");
            if(f0&&!f0->IsZombie()){
                TTree* tr=(TTree*)f0->Get("ch1");
                if(tr&&tr->GetEntries()>1){
                    Double_t tb[N_SAMPLES];
                    tr->SetBranchAddress("time",tb);tr->GetEntry(0);
                    double dt=tb[1]-tb[0];if(dt>0)fs_MHz=1000.0/dt;
                }
                f0->Close();delete f0;
            }
        }
    }

    std::cout<<"  Gain = "<<gain_mV_pe<<" mV/p.e.\n"
             <<"  fs   = "<<fs_MHz<<" MHz\n"
             <<"  Trigger window = ["<<t_win_start<<", "<<t_win_end<<"] ns\n"
             <<"  Fit range      = ["<<fit_lo<<", "<<fit_hi<<"] ns\n\n";

    bool highDirOk = (gSystem->AccessPathName(dataDirHigh.c_str()) == 0);
    if (!highDirOk)
        std::cerr<<"[WARN] High light dir not found: "<<dataDirHigh<<"\n\n";

    std::cout<<"--- METODO m (slope) ---\n";
    
    int j_start=0, j_end=N_SAMPLES-1;
    {
        auto runs = findRunFiles(vbias, dataDirLow);
        if (!runs.empty()) {
            TFile* f0 = TFile::Open(runs.begin()->second.c_str(),"READ");
            if (f0 && !f0->IsZombie()) {
                TTree* tr = (TTree*)f0->Get("ch1");
                if (tr) {
                    Double_t tb[N_SAMPLES];
                    tr->SetBranchAddress("time", tb);
                    tr->GetEntry(0);
                    triggerWindowIndices(tb, N_SAMPLES, t_win_start, t_win_end, 
                                         j_start, j_end);
                }
                f0->Close(); delete f0;
            }
        }
    }
    
    int cal_m_tenth = (int)std::lround(gain_mV_pe * 10.0);
    int cal_q_tenth = (int)std::lround(cal.q * 10.0);
    extern int g_analysis_mode;
    std::string cacheSuffix = Form("_twv1_js%d_je%d_m%d_q%d_mode%d",
                                    j_start, j_end, cal_m_tenth, cal_q_tenth, 
                                    g_analysis_mode);
    
    bool use_filter = (cutoff_MHz > 0);
    std::string cacheLo = eventCachePath(vbias, thr_opt, cutoff_MHz, 
                                          cal.laser_thr, dataDirLow, 
                                          use_filter, cacheSuffix);
    std::string cacheHi = "";
    if (highDirOk) {
        cacheHi = eventCachePath(vbias, thr_opt, cutoff_MHz, 
                                  cal.laser_thr, dataDirHigh,
                                  use_filter, cacheSuffix);
        if (gSystem->AccessPathName(cacheHi.c_str())) cacheHi = "";
    }
    
    if (gSystem->AccessPathName(cacheLo.c_str())) {
        std::cout << "  [m] No cache found: " << cacheLo << "\n";
        cacheLo = "";
    } else {
        std::cout << "  [m] Using cache: " << cacheLo << "\n";
    }
    if (!cacheHi.empty()) {
        std::cout << "  [m] High light cache: " << cacheHi << "\n";
    }

    std::vector<MEvent> mEvLo, mEvHi;
    extractMEvents(cacheLo, pe_lo, pe_hi, gain_mV_pe, mEvLo);
    if (!cacheHi.empty())
        extractMEvents(cacheHi, pe_lo, pe_hi, gain_mV_pe, mEvHi);
    
    if (mEvLo.empty())
        std::cout << "  [m] No cache found - skipping m-method analysis\n";

    std::cout<<"\n--- METODO dV/dt ---\n";
    auto runsLo = findRunFiles(vbias, dataDirLow);
    auto runsHi = highDirOk ? findRunFiles(vbias, dataDirHigh)
                            : std::map<int,std::string>{};

    std::vector<DvdtEvent> dvLo, dvHi;
    extractDvdt(runsLo, gain_mV_pe, cutoff_MHz, fs_MHz,
                t_win_start, t_win_end, "low",  dvLo, 1000000);
    if (!runsHi.empty())
        extractDvdt(runsHi, gain_mV_pe, cutoff_MHz, fs_MHz,
                    t_win_start, t_win_end, "high", dvHi, 1000000);

    if (mEvLo.empty() && dvLo.empty()){
        std::cerr<<"[ERROR] No events extracted. Check cache/raw files.\n";
        return;
    }
    
    if (mEvLo.empty()) {
        std::cout << "\n[INFO] Running dV/dt-only analysis (no m-method cache)\n";
    }

    TH1D *hM_lo=nullptr, *hM_hi=nullptr;
    TH2D *h2M_lo=nullptr, *h2M_hi=nullptr;
    
    if (!mEvLo.empty() || !mEvHi.empty()) {
        hM_lo = makeH1("hM_lo", ";m (ns/p.e.);Freq. norm.", M_BINS, M_RANGE_MIN, M_RANGE_MAX, kAzure+1, 3244);
        hM_hi = makeH1("hM_hi", ";m (ns/p.e.);Freq. norm.", M_BINS, M_RANGE_MIN, M_RANGE_MAX, kRed+1,   3245);

        h2M_lo = new TH2D("h2M_lo",
            Form("TOT vs m  Vbias %d V  low light;TOT (ns);m (ns/p.e.)",vbias),
            50,30,80, M_BINS, M_RANGE_MIN, M_RANGE_MAX);
        h2M_hi = new TH2D("h2M_hi",
            Form("TOT vs m  Vbias %d V  high light;TOT (ns);m (ns/p.e.)",vbias),
            50,30,80, M_BINS, M_RANGE_MIN, M_RANGE_MAX);
        h2M_lo->SetDirectory(nullptr); h2M_hi->SetDirectory(nullptr);

        for (auto& e:mEvLo){ hM_lo->Fill(e.m); h2M_lo->Fill(e.tot,e.m); }
        for (auto& e:mEvHi){ hM_hi->Fill(e.m); h2M_hi->Fill(e.tot,e.m); }
    }

    TH1D* hPR_lo = makeH1("hPR_lo",";dV/dt_{2}/dV/dt_{1};Freq. norm.", 100,0,0.6, kAzure+1, 3244);
    TH1D* hPR_hi = makeH1("hPR_hi",";dV/dt_{2}/dV/dt_{1};Freq. norm.", 100,0,0.6, kRed+1,   3245);
    TH1D* hRT_lo = makeH1("hRT_lo",";Rise time 10%-90% (ns);Freq. norm.",100,0.5,2.5,kAzure+1, 3244);
    TH1D* hRT_hi = makeH1("hRT_hi",";Rise time 10%-90% (ns);Freq. norm.",100,0.5,2.5,kRed+1,   3245);

    for (auto& e:dvLo){ hPR_lo->Fill(e.peak2_ratio); hRT_lo->Fill(e.rise_time); }
    for (auto& e:dvHi){ hPR_hi->Fill(e.peak2_ratio); hRT_hi->Fill(e.rise_time); }

    std::vector<TH1D*> all_h = {hPR_lo, hPR_hi, hRT_lo, hRT_hi};
    if (hM_lo) all_h.push_back(hM_lo);
    if (hM_hi) all_h.push_back(hM_hi);
    for(auto* h : all_h) normArea(h);

    double sep_m=0, sep_pr=0, cut_m=0, cut_pr=0;
    std::vector<double> vm_lo,vm_hi,vpr_lo,vpr_hi;
    
    if (!mEvLo.empty() || !mEvHi.empty()) {
        for(auto& e:mEvLo)  vm_lo .push_back(e.m);
        for(auto& e:mEvHi)  vm_hi .push_back(e.m);
        cut_m = fisherCut(vm_lo, vm_hi, sep_m);
        std::cout<<Form("\n[Fisher m]     cut=%.3f ns/p.e.  sep=%.2f sigma\n",cut_m, sep_m);
    }
    
    for(auto& e:dvLo)   vpr_lo.push_back(e.peak2_ratio);
    for(auto& e:dvHi)   vpr_hi.push_back(e.peak2_ratio);
    cut_pr = fisherCut(vpr_lo,vpr_hi,sep_pr);
    std::cout<<Form("[Fisher dV/dt] cut=%.3f           sep=%.2f sigma\n", cut_pr,sep_pr);

    OutCtx ctx = createOutputDirs("xtalk_disc");
    
    if (gSystem->AccessPathName("../canvas")) {
        std::cerr << "[WARN] ../canvas/ non esiste — la creo\n";
        gSystem->mkdir("../canvas", kTRUE);
    }
    if (gSystem->AccessPathName("../file_root")) {
        std::cerr << "[WARN] ../file_root/ non esiste — la creo\n";
        gSystem->mkdir("../file_root", kTRUE);
    }

    auto makeLegSimple = [](TH1D* h, const char* label, int n, double x1, double y1, double w=0.38, double hh=0.12)->TLegend*{
        auto* leg = new TLegend(x1,y1,x1+w,y1+hh);
        leg->SetBorderSize(1); leg->SetFillColor(0);
        leg->SetTextFont(42); leg->SetTextSize(0.032);
        leg->AddEntry(h, Form("%s  N=%d",label,n), "lf");
        return leg;
    };
    auto makePT = [&](const char* line2)->TPaveText*{
        auto* pt=new TPaveText(PAD_LEFT+0.02,0.78,0.50,0.92,"NDC");
        pt->SetBorderSize(1);pt->SetFillColor(0);
        pt->SetTextFont(42); pt->SetTextSize(0.032);
        pt->AddText(Form("Vbias = %d V",vbias));
        pt->AddText(line2);
        return pt;
    };

    // ── Multifit doppia gaussiana su low e high light ─────────────
    // Utile per quantificare la separazione fra picco principale (segnale)
    // e spalla di crosstalk in entrambe le condizioni di luce.
    DG2Result dg_lo, dg_hi;
    if (hM_lo && !mEvLo.empty())
        dg_lo = fitDoubleGaussian(hM_lo, Form("lo_v%d",vbias), kAzure+2, kCyan+1);
    if (hM_hi && !mEvHi.empty())
        dg_hi = fitDoubleGaussian(hM_hi, Form("hi_v%d",vbias), kRed+2, kOrange+7);

    // ── Canvas: m - low light ──────────────────────────────────────
    // Disegna l'histogram e il fit a doppia gaussiana (se convergito).
    if (hM_lo && !mEvLo.empty()) {
        auto* c = new TCanvas(Form("c1_m_low_%d",vbias), "m - low light", 900,650);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = hM_lo->GetMaximum()*1.40;
        hM_lo->GetYaxis()->SetRangeUser(0,ymax);
        hM_lo->Draw("HIST");
        if (dg_lo.ok) {
            dg_lo.f_total->Draw("same");
            dg_lo.f_g1->Draw("same");
            dg_lo.f_g2->Draw("same");
            TLegend* leg = new TLegend(0.50, 0.55, 0.93, 0.88);
            leg->SetBorderSize(1); leg->SetFillColor(0);
            leg->SetTextFont(42); leg->SetTextSize(0.028);
            leg->AddEntry(hM_lo, Form("Low light  N=%d", dg_lo.N), "lf");
            leg->AddEntry(dg_lo.f_total, "Double Gaussian fit", "l");
            leg->AddEntry(dg_lo.f_g1,
                Form("Peak 1: #mu=%.4f #pm %.4f  #sigma=%.4f",
                     dg_lo.mu1, dg_lo.mu1Err, dg_lo.sigma1), "l");
            leg->AddEntry(dg_lo.f_g2,
                Form("Peak 2: #mu=%.4f #pm %.4f  #sigma=%.4f",
                     dg_lo.mu2, dg_lo.mu2Err, dg_lo.sigma2), "l");
            leg->Draw();
        } else {
            makeLegSimple(hM_lo, "Low light", (int)mEvLo.size(), 0.65,0.80)->Draw();
        }
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c1_m_low_vbias%d.png",vbias));
    }

    // ── Canvas: m - high light ─────────────────────────────────────
    if(!mEvHi.empty()){
        auto* c = new TCanvas(Form("c1_m_high_%d",vbias), "m - high light", 900,650);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = hM_hi->GetMaximum()*1.40;
        hM_hi->GetYaxis()->SetRangeUser(0,ymax);
        hM_hi->Draw("HIST");
        if (dg_hi.ok) {
            dg_hi.f_total->Draw("same");
            dg_hi.f_g1->Draw("same");
            dg_hi.f_g2->Draw("same");
            TLegend* leg = new TLegend(0.50, 0.55, 0.93, 0.88);
            leg->SetBorderSize(1); leg->SetFillColor(0);
            leg->SetTextFont(42); leg->SetTextSize(0.028);
            leg->AddEntry(hM_hi, Form("High light  N=%d", dg_hi.N), "lf");
            leg->AddEntry(dg_hi.f_total, "Double Gaussian fit", "l");
            leg->AddEntry(dg_hi.f_g1,
                Form("Peak 1: #mu=%.4f #pm %.4f  #sigma=%.4f",
                     dg_hi.mu1, dg_hi.mu1Err, dg_hi.sigma1), "l");
            leg->AddEntry(dg_hi.f_g2,
                Form("Peak 2: #mu=%.4f #pm %.4f  #sigma=%.4f",
                     dg_hi.mu2, dg_hi.mu2Err, dg_hi.sigma2), "l");
            leg->Draw();
        }
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c1_m_high_vbias%d.png",vbias));
    }

    // ── Canvas: m - overlay con DOPPIO multifit per confronto ────
    // Mostra sovrapposti low e high light + entrambi i fit + valori
    // sigma/mu per verificare la compatibilita' tra le due distribuzioni.
    if (hM_lo || hM_hi) {
        TH1D* hRef_lo = hM_lo ? hM_lo : hM_hi;
        TH1D* hRef_hi = hM_hi ? hM_hi : hM_lo;
        auto* c = new TCanvas(Form("c1_m_overlay_%d",vbias), "m - overlay", 1100, 750);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = std::max(hRef_lo->GetMaximum(), hRef_hi->GetMaximum())*1.45;
        hRef_lo->GetYaxis()->SetRangeUser(0,ymax);
        hRef_lo->Draw("HIST");
        if(hM_hi && !mEvHi.empty()) hRef_hi->Draw("HIST SAME");

        // Disegno fit low (se convergito)
        if (dg_lo.ok && hM_lo) {
            // Cloni indipendenti per non interferire col canvas low gia' salvato
            TF1* fT_lo = (TF1*)dg_lo.f_total->Clone(Form("ovl_fT_lo_v%d",vbias));
            TF1* fG1_lo = (TF1*)dg_lo.f_g1->Clone(Form("ovl_g1_lo_v%d",vbias));
            TF1* fG2_lo = (TF1*)dg_lo.f_g2->Clone(Form("ovl_g2_lo_v%d",vbias));
            fT_lo->SetLineColor(kAzure+1); fT_lo->SetLineStyle(1);
            fG1_lo->SetLineColor(kAzure+2); fG1_lo->SetLineStyle(2);
            fG2_lo->SetLineColor(kCyan+1);  fG2_lo->SetLineStyle(2);
            fT_lo->Draw("same"); fG1_lo->Draw("same"); fG2_lo->Draw("same");
        }
        // Disegno fit high (se convergito)
        if (dg_hi.ok && hM_hi && !mEvHi.empty()) {
            TF1* fT_hi = (TF1*)dg_hi.f_total->Clone(Form("ovl_fT_hi_v%d",vbias));
            TF1* fG1_hi = (TF1*)dg_hi.f_g1->Clone(Form("ovl_g1_hi_v%d",vbias));
            TF1* fG2_hi = (TF1*)dg_hi.f_g2->Clone(Form("ovl_g2_hi_v%d",vbias));
            fT_hi->SetLineColor(kRed+1); fT_hi->SetLineStyle(1);
            fG1_hi->SetLineColor(kRed+2); fG1_hi->SetLineStyle(2);
            fG2_hi->SetLineColor(kOrange+7); fG2_hi->SetLineStyle(2);
            fT_hi->Draw("same"); fG1_hi->Draw("same"); fG2_hi->Draw("same");
        }

        // Legenda completa con tutti i parametri di fit
        TLegend* leg = new TLegend(0.45, 0.45, 0.95, 0.90);
        leg->SetBorderSize(1); leg->SetFillColor(0);
        leg->SetTextFont(42); leg->SetTextSize(0.024);
        if(hM_lo) leg->AddEntry(hM_lo,
            Form("Low light  N=%d",(int)mEvLo.size()), "lf");
        if(hM_hi && !mEvHi.empty())
            leg->AddEntry(hM_hi,
                Form("High light  N=%d",(int)mEvHi.size()), "lf");
        if (dg_lo.ok && hM_lo) {
            leg->AddEntry((TObject*)nullptr, "  --- Low light fit ---", "");
            leg->AddEntry((TObject*)nullptr,
                Form("  Peak1: #mu=%.4f#pm%.4f  #sigma=%.4f",
                     dg_lo.mu1, dg_lo.mu1Err, dg_lo.sigma1), "");
            leg->AddEntry((TObject*)nullptr,
                Form("  Peak2: #mu=%.4f#pm%.4f  #sigma=%.4f",
                     dg_lo.mu2, dg_lo.mu2Err, dg_lo.sigma2), "");
        }
        if (dg_hi.ok && hM_hi && !mEvHi.empty()) {
            leg->AddEntry((TObject*)nullptr, "  --- High light fit ---", "");
            leg->AddEntry((TObject*)nullptr,
                Form("  Peak1: #mu=%.4f#pm%.4f  #sigma=%.4f",
                     dg_hi.mu1, dg_hi.mu1Err, dg_hi.sigma1), "");
            leg->AddEntry((TObject*)nullptr,
                Form("  Peak2: #mu=%.4f#pm%.4f  #sigma=%.4f",
                     dg_hi.mu2, dg_hi.mu2Err, dg_hi.sigma2), "");
        }
        // Compatibilita' tra fit (numero di sigma di separazione tra le mu corrispondenti)
        if (dg_lo.ok && dg_hi.ok) {
            // Confronto picco 1 lo vs picco 1 hi
            double dmu1 = dg_hi.mu1 - dg_lo.mu1;
            double sig1 = std::sqrt(dg_lo.mu1Err*dg_lo.mu1Err
                                  + dg_hi.mu1Err*dg_hi.mu1Err);
            double comp1 = (sig1 > 0) ? std::abs(dmu1)/sig1 : 0.0;
            // Confronto picco 2
            double dmu2 = dg_hi.mu2 - dg_lo.mu2;
            double sig2 = std::sqrt(dg_lo.mu2Err*dg_lo.mu2Err
                                  + dg_hi.mu2Err*dg_hi.mu2Err);
            double comp2 = (sig2 > 0) ? std::abs(dmu2)/sig2 : 0.0;
            leg->AddEntry((TObject*)nullptr,
                Form("  Compat. #mu_{1}: %.1f#sigma   #mu_{2}: %.1f#sigma",
                     comp1, comp2), "");
        }
        leg->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c1_m_overlay_vbias%d.png",vbias));
    }

    if (h2M_lo && (h2M_lo->GetEntries() > 0 || h2M_hi->GetEntries() > 0)) {
        auto* c=new TCanvas(Form("c2_scatter_%d",vbias),"TOT vs m",1200,500);
        c->Divide(2,1);
        for(int side=1;side<=2;++side){
            TVirtualPad* pad = c->cd(side);
            if (!pad) continue; // safety check
            pad->SetLogz(); 
            pad->SetGrid();
            pad->SetLeftMargin(PAD_LEFT); 
            pad->SetBottomMargin(PAD_BOTTOM);
            pad->SetTopMargin(PAD_TOP);   
            pad->SetRightMargin(0.15f);
            auto* h = (side==1) ? h2M_lo : h2M_hi;
            if (!h || h->GetEntries() == 0) continue;
            h->Draw("COLZ");
        }
        ctx.savePNG(c,Form("c2_scatter_m_vbias%d.png",vbias));
    }

    {
        auto* c = new TCanvas(Form("c3_dvdt_low_%d",vbias), "dV/dt ratio - low light", 800,600);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = hPR_lo->GetMaximum()*1.30;
        hPR_lo->GetYaxis()->SetRangeUser(0,ymax);
        hPR_lo->Draw("HIST");
        makeLegSimple(hPR_lo, "Low light", (int)dvLo.size(), 0.65,0.80)->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c3_dvdt_low_vbias%d.png",vbias));
    }
    if(!dvHi.empty()){
        auto* c = new TCanvas(Form("c3_dvdt_high_%d",vbias), "dV/dt ratio - high light", 800,600);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = hPR_hi->GetMaximum()*1.30;
        hPR_hi->GetYaxis()->SetRangeUser(0,ymax);
        hPR_hi->Draw("HIST");
        makeLegSimple(hPR_hi, "High light", (int)dvHi.size(), 0.65,0.80)->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c3_dvdt_high_vbias%d.png",vbias));
    }
    {
        auto* c = new TCanvas(Form("c3_dvdt_overlay_%d",vbias), "dV/dt ratio - overlay", 800,600);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = std::max(hPR_lo->GetMaximum(), hPR_hi->GetMaximum())*1.30;
        hPR_lo->GetYaxis()->SetRangeUser(0,ymax);
        hPR_lo->Draw("HIST");
        if(!dvHi.empty()) hPR_hi->Draw("HIST SAME");
        TLegend* leg = new TLegend(0.65,0.70,0.90,0.88);
        leg->AddEntry(hPR_lo, Form("Low light  N=%d",(int)dvLo.size()), "lf");
        if(!dvHi.empty()) leg->AddEntry(hPR_hi, Form("High light   N=%d",(int)dvHi.size()), "lf");
        leg->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c3_dvdt_overlay_vbias%d.png",vbias));
    }

    {
        auto* c = new TCanvas(Form("c4_rt_low_%d",vbias), "Rise time - low light", 800,600);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = hRT_lo->GetMaximum()*1.30;
        hRT_lo->GetYaxis()->SetRangeUser(0,ymax);
        hRT_lo->Draw("HIST");
        makeLegSimple(hRT_lo, "Low light", (int)dvLo.size(), 0.65,0.80)->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c4_rise_time_low_vbias%d.png",vbias));
    }
    if(!dvHi.empty()){
        auto* c = new TCanvas(Form("c4_rt_high_%d",vbias), "Rise time - high light", 800,600);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = hRT_hi->GetMaximum()*1.30;
        hRT_hi->GetYaxis()->SetRangeUser(0,ymax);
        hRT_hi->Draw("HIST");
        makeLegSimple(hRT_hi, "High light", (int)dvHi.size(), 0.65,0.80)->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c4_rise_time_high_vbias%d.png",vbias));
    }
    {
        auto* c = new TCanvas(Form("c4_rt_overlay_%d",vbias), "Rise time - overlay", 800,600);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        double ymax = std::max(hRT_lo->GetMaximum(), hRT_hi->GetMaximum())*1.30;
        hRT_lo->GetYaxis()->SetRangeUser(0,ymax);
        hRT_lo->Draw("HIST");
        if(!dvHi.empty()) hRT_hi->Draw("HIST SAME");
        TLegend* leg = new TLegend(0.65,0.70,0.90,0.88);
        leg->AddEntry(hRT_lo, Form("Low light  N=%d",(int)dvLo.size()), "lf");
        if(!dvHi.empty()) leg->AddEntry(hRT_hi, Form("High light   N=%d",(int)dvHi.size()), "lf");
        leg->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c4_rise_time_overlay_vbias%d.png",vbias));
    }

    if (!vm_lo.empty() || !vm_hi.empty()) {
        double mMin=vm_lo.empty()?M_RANGE_MIN:*std::min_element(vm_lo.begin(),vm_lo.end());
        double mMax=vm_lo.empty()?M_RANGE_MAX:*std::max_element(vm_lo.begin(),vm_lo.end());
        if(!vm_hi.empty()){
            mMin=std::min(mMin,*std::min_element(vm_hi.begin(),vm_hi.end()));
            mMax=std::max(mMax,*std::max_element(vm_hi.begin(),vm_hi.end()));
        }
        double mRange=mMax-mMin>1e-9?mMax-mMin:1.0;

        TH1D* hMN_lo=makeH1("hMN_lo",";Normalized discriminant;Freq. norm.",100,0,0.3,kAzure+1,0);
        TH1D* hMN_hi=makeH1("hMN_hi",";Normalized discriminant;Freq. norm.",100,0,0.3,kAzure-4,0);
        TH1D* hPN_lo=makeH1("hPN_lo",";Normalized discriminant;Freq. norm.",100,0,0.3,kRed+1,0,2);
        TH1D* hPN_hi=makeH1("hPN_hi",";Normalized discriminant;Freq. norm.",100,0,0.3,kOrange+7,0,2);

        for(double v:vm_lo)  hMN_lo->Fill((v-mMin)/mRange);
        for(double v:vm_hi)  hMN_hi->Fill((v-mMin)/mRange);
        for(double v:vpr_lo) hPN_lo->Fill(v);
        for(double v:vpr_hi) hPN_hi->Fill(v);
        for(auto* h:{hMN_lo,hMN_hi,hPN_lo,hPN_hi}) normArea(h);

        double ymax=0;
        for(auto* h:{hMN_lo,hMN_hi,hPN_lo,hPN_hi}) ymax=std::max(ymax,h->GetMaximum());
        ymax*=1.35;

        auto* c=new TCanvas(Form("c5_confronto_%d",vbias),"Discriminant comparison",900,600);
        c->SetGrid(); c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
        hMN_lo->GetYaxis()->SetRangeUser(0,ymax);
        hMN_lo->Draw("HIST"); hMN_hi->Draw("HIST SAME");
        hPN_lo->Draw("HIST SAME"); hPN_hi->Draw("HIST SAME");

        auto* leg=new TLegend(0.50,0.62,0.88,0.88);
        leg->SetBorderSize(1); leg->SetFillColor(0);
        leg->SetTextFont(42); leg->SetTextSize(0.028);
        leg->AddEntry(hMN_lo,"m low light",  "lf");
        leg->AddEntry(hMN_hi,"m high light",  "lf");
        leg->AddEntry(hPN_lo,"dV/dt low light","l");
        leg->AddEntry(hPN_hi,"dV/dt high light","l");
        leg->Draw();
        makePT(Form("Amp [%.1f,%.1f] p.e.",AMP_LO_PE,AMP_HI_PE))->Draw();
        ctx.savePNG(c,Form("c5_confronto_vbias%d.png",vbias));
        // FIX double-free: hMN/hPN drawn su c → NON deletare dopo savePNG.
        // SetDirectory(nullptr) + drawn → distrutti da ~TPad con la canvas.
        // delete hMN_lo; delete hMN_hi; delete hPN_lo; delete hPN_hi; // RIMOSSO
    }

    // FIX BUG FISICO: il vecchio codice faceva h2c->Fill(mEvLo[i], dvLo[i])
    // accoppiando l'i-esimo evento del metodo-m con l'i-esimo del metodo-dV/dt.
    // I due vettori vengono da sorgenti DIVERSE (cache vs raw) e l'accoppiamento
    // per indice non ha senso fisico → correlazione completamente spuria.
    // Plot rimosso. Se vuoi una correlazione vera serve un loop evento-per-evento
    // sul raw con entrambe le variabili calcolate insieme (come in extractDvdt
    // + slope contestuale).
    if (false && !mEvLo.empty() && !dvLo.empty()) {  // DISABILITATO

        double mMin=*std::min_element(vm_lo.begin(),vm_lo.end());
        double mMax=*std::max_element(vm_lo.begin(),vm_lo.end());
        double mRange = (mMax - mMin);
        if(mRange < 1e-9) mRange = 1.0;

        TH2D* h2c=new TH2D(Form("h2corr_%d",vbias),
            Form("Correlation m vs dV/dt  Vbias %d V  low light;"
                 "m normalized;dV/dt_{2}/dV/dt_{1}",vbias),
            100,0,0.2,100,0,0.2);
        h2c->SetDirectory(nullptr);
        size_t Nc=std::min(mEvLo.size(),dvLo.size());
        for(size_t i=0;i<Nc;++i)
            h2c->Fill((mEvLo[i].m-mMin)/mRange,dvLo[i].peak2_ratio);

        if (h2c->GetEntries() > 0) {
            auto* c=new TCanvas(Form("c6_corr_%d",vbias),"Correlation m vs dV/dt",750,600);
            c->SetLogz(); c->SetGrid();
            c->SetLeftMargin(PAD_LEFT); c->SetBottomMargin(PAD_BOTTOM);
            c->SetTopMargin(PAD_TOP);  c->SetRightMargin(0.15f);
            h2c->Draw("COLZ");
            ctx.savePNG(c,Form("c6_correlation_vbias%d.png",vbias));
        }
        delete h2c;
    }

    // FIX: ogni esecuzione produce un root file CON timestamp nel nome.
    // Prima xtalk_disc_vbias%d.root sovrascriveva sempre lo stesso file
    // (perdita risultati precedenti). Ora il timestamp ISO-like garantisce
    // unicita': xtalk_disc_vbias55_20260429_143021.root
    time_t now_t = time(nullptr);
    struct tm* tm_now = localtime(&now_t);
    char ts_buf[32];
    std::strftime(ts_buf, sizeof(ts_buf), "%Y%m%d_%H%M%S", tm_now);
    std::string rootOut=Form("%s/xtalk_disc_vbias%d_%s.root",
                              ctx.rootDir.c_str(), vbias, ts_buf);
    TFile* fOut=new TFile(rootOut.c_str(),"RECREATE");
    if (hM_lo) hM_lo->Write(); 
    if (hM_hi) hM_hi->Write();
    hPR_lo->Write();hPR_hi->Write();
    hRT_lo->Write();hRT_hi->Write();
    if (h2M_lo) h2M_lo->Write();
    if (h2M_hi) h2M_hi->Write();
    TParameter<double> pCutM ("cut_m",        cut_m);
    TParameter<double> pSepM ("sep_m_sigma",  sep_m);
    TParameter<double> pCutPR("cut_dvdt",     cut_pr);
    TParameter<double> pSepPR("sep_dvdt_sigma",sep_pr);
    TParameter<double> pGain ("gain_mV_pe",   gain_mV_pe);
    pCutM.Write();pSepM.Write();pCutPR.Write();pSepPR.Write();pGain.Write();

    // FIX: salva anche i parametri del fit a doppia gaussiana per low/high.
    if (dg_lo.ok) {
        TParameter<double>("fit_lo_mu1",       dg_lo.mu1).Write();
        TParameter<double>("fit_lo_mu1_err",   dg_lo.mu1Err).Write();
        TParameter<double>("fit_lo_sigma1",    dg_lo.sigma1).Write();
        TParameter<double>("fit_lo_sigma1_err",dg_lo.sigma1Err).Write();
        TParameter<double>("fit_lo_mu2",       dg_lo.mu2).Write();
        TParameter<double>("fit_lo_mu2_err",   dg_lo.mu2Err).Write();
        TParameter<double>("fit_lo_sigma2",    dg_lo.sigma2).Write();
        TParameter<double>("fit_lo_sigma2_err",dg_lo.sigma2Err).Write();
    }
    if (dg_hi.ok) {
        TParameter<double>("fit_hi_mu1",       dg_hi.mu1).Write();
        TParameter<double>("fit_hi_mu1_err",   dg_hi.mu1Err).Write();
        TParameter<double>("fit_hi_sigma1",    dg_hi.sigma1).Write();
        TParameter<double>("fit_hi_sigma1_err",dg_hi.sigma1Err).Write();
        TParameter<double>("fit_hi_mu2",       dg_hi.mu2).Write();
        TParameter<double>("fit_hi_mu2_err",   dg_hi.mu2Err).Write();
        TParameter<double>("fit_hi_sigma2",    dg_hi.sigma2).Write();
        TParameter<double>("fit_hi_sigma2_err",dg_hi.sigma2Err).Write();
    }
    fOut->Write();fOut->Close();delete fOut;

    // FIX double-free: tutti gli histogram qui sotto (hM_lo, hM_hi, hPR_*,
    // hRT_*, h2M_*) sono stati drawn sui canvas tramite ctx.savePNG.
    // I canvas sono gia' stati distrutti da savePNG → ROOT cancella anche
    // gli histogram associati. Deletarli qui = double-free.
    // Vecchio codice:
    //   if (hM_lo) delete hM_lo; ... if (h2M_hi) delete h2M_hi;
    // Rimosso. Memory leak < crash da double-free.
    //
    // Anche dg_lo.f_total, dg_lo.f_g1, dg_lo.f_g2 sono drawn sui canvas
    // low e overlay (e clonati per overlay). Stessa logica: non deletare,
    // ROOT li cancella con i pad.

    std::cout<<"\n[OK] Canvas → "<<ctx.pngDir<<"\n"
             <<"[OK] ROOT   → "<<rootOut   <<"\n"
             <<Form("[OK] Taglio m       = %.3f ns/p.e.  (%.2f sigma)\n",cut_m, sep_m)
             <<Form("[OK] Taglio dV/dt   = %.3f          (%.2f sigma)\n",cut_pr,sep_pr);
}

void sipm_xtalk_discrimination(int vbias = -1, 
                                double thr_opt = -1.0,
                                double fit_lo = -1.0,
                                double fit_hi = -1.0)
{
    std::cout << "\n╔═══════════════════════════════════════════════════╗\n";
    std::cout << "║  SiPM CROSSTALK vs 2PE DISCRIMINATION            ║\n";
    std::cout << "╚═══════════════════════════════════════════════════╝\n\n";
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    
    if (vbias < 0) {
        std::cout << "Vbias [V]: "; 
        std::cin >> vbias;
    }
    if (thr_opt < 0) {
        std::cout << "LET threshold [p.e.]: "; 
        std::cin >> thr_opt;
    }
    if (fit_lo < 0) {
        std::cout << "Fit range [ns] (lo hi): "; 
        std::cin >> fit_lo >> fit_hi;
    }
    
    auto filters = scanAvailableFilters(vbias, DATA_DIR);
    if (filters.empty()) {
        std::cerr << "[ERROR] No calibration files found for vbias=" 
                  << vbias << "\n";
        return;
    }
    
    double cutoff_MHz = 0;
    if (filters.size() == 1) {
        cutoff_MHz = filters[0];
    } else {
        std::cout << "\nAvailable filters [MHz]: ";
        for (size_t i=0; i<filters.size(); ++i) {
            std::cout << filters[i];
            if (i < filters.size()-1) std::cout << ", ";
        }
        std::cout << "\nSelect cutoff [MHz]: ";
        std::cin >> cutoff_MHz;
        
        if (std::find(filters.begin(), filters.end(), cutoff_MHz) == filters.end()) {
            std::cerr << "[ERROR] Invalid cutoff " << cutoff_MHz << " MHz\n";
            return;
        }
    }
    
    std::cout << "\n>>> Running analysis:\n"
              << "    Vbias     = " << vbias << " V\n"
              << "    Threshold = " << thr_opt << " p.e.\n"
              << "    Cutoff    = " << cutoff_MHz << " MHz\n"
              << "    Fit range = [" << fit_lo << ", " << fit_hi << "] ns\n\n";
    
    analyzeVbias(vbias, thr_opt, cutoff_MHz, fit_lo, fit_hi);
    
    std::cout << "\n[DONE] Analysis complete.\n";
}
