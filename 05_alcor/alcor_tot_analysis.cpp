// =============================================================================
// alcor_tot_analysis.cpp
// ALCOR SiPM Calibration — Time over Threshold (ToT) Analysis
//
// Usage:
//   root -l 'alcor_tot_analysis.cpp+'
//   then call: alcor_tot_analysis()
//
// Input format (auto-detected at runtime):
//   ROOT TTree  — branches: t_leading [double], t_trailing [double],
//                            t_laser [double]  (optional laser reference)
//   CSV         — columns:  t_leading, t_trailing, t_laser (header line)
//   SIMULATION  — internal Monte Carlo (for testing without real data)
//
// TDC conversion:
//   ToT_ns = (t_trailing - t_leading) * lsb_ns
//   Default lsb_ns = 1.0  (already in ns). Change at startup prompt
//   if ALCOR delivers raw TDC ticks (typical: 0.05 ns = 50 ps LSB).
//
// Analyses performed (15 canvases, all saved to plots/):
//   01a  ToT spectrum — log scale
//   01b  ToT spectrum — linear scale
//   02   ToT vs p.e. calibration curve  (non-linear fit)
//   03   Delta-t spectrum — raw (before time-walk correction)
//   04   2D map: Delta-t vs ToT  (time walk visible as diagonal band)
//   05   Time-walk correction fit  (t_lead vs ToT)
//   06   Delta-t spectrum — after time-walk correction
//   07   Before vs after overlay
//   08   ToT slices: timing distribution per p.e. band
//   09   Timing resolution vs p.e. — before correction
//   10   Timing resolution vs p.e. — after correction
//   11   SNR: peak separation / sigma per p.e. band
//   12   ToT stability: running mean over event index
//   13   Inter-arrival time distribution  (dark rate proxy)
//   14   Gain linearity: ToT peak position vs n p.e.
//   15   Summary panel
//
// Official code language: English
// =============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <numeric>

#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TLine.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TRandom3.h>

// =============================================================================
// SIMULATION PARAMETERS  (used only in simulation mode)
// =============================================================================
static const int    SIM_N_EVENTS    = 200000;
static const double SIM_TDC_BIN_NS  = 0.050;   // 50 ps — ALCOR nominal LSB
static const double SIM_T0          = 50.0;     // ns
static const double SIM_LASER_JIT   = 0.10;     // ns sigma
static const double SIM_TW_A        = 2.5;      // time-walk amplitude (ns)
static const double SIM_SIGMA_STOC  = 0.20;     // ns
static const double SIM_SIGMA_CONST = 0.05;     // ns
static const double SIM_TOT_BASE    = 10.0;     // ns at 1 p.e.
static const double SIM_TOT_SLOPE   = 8.0;      // ns/p.e.
static const double SIM_TOT_SIGMA   = 1.5;      // ns smearing per p.e.
static const double SIM_MU_LASER    = 2.5;      // mean p.e. per pulse
static const int    SIM_NPE_MAX     = 8;

// =============================================================================
// INTERNAL EVENT STRUCTURE
// =============================================================================
struct AlcorEvent {
    double tot_raw;    // raw ToT (ticks or ns, before lsb conversion)
    double t_lead;     // leading-edge time (same units as tot_raw)
    double t_laser;    // laser reference time (-999 if not available)
};

// =============================================================================
// UTILITY: repeat a string n times (ASCII box drawing)
// =============================================================================
static std::string repeatStr(const std::string& s, size_t n) {
    std::string out; out.reserve(s.size()*n);
    for(size_t i=0;i<n;++i) out+=s;
    return out;
}

static void printBox(const std::vector<std::string>& lines) {
    size_t maxLen = 0;
    for(const auto& l:lines) if(l.size()>maxLen) maxLen=l.size();
    size_t w = maxLen+4;
    std::cout<<"\n  +"<<repeatStr("-",w)<<"+"<<std::endl;
    for(const auto& l:lines){
        size_t pad=w-l.size()-2;
        std::cout<<"  |  "<<l<<std::string(pad,' ')<<"|"<<std::endl;
    }
    std::cout<<"  +"<<repeatStr("-",w)<<"+"<<std::endl;
}

// =============================================================================
// UTILITY: yes/no prompt
// =============================================================================
static char askYesNo(const std::string& prompt) {
    char c=0;
    while(c!='y'&&c!='n'){
        std::cout<<prompt;
        std::string line; std::getline(std::cin,line);
        size_t p=line.find_first_not_of(" \t\r\n");
        if(p!=std::string::npos) c=std::tolower((unsigned char)line[p]);
        if(c!='y'&&c!='n') std::cout<<"  [ERROR] Please enter 'y' or 'n'."<<std::endl;
    }
    return c;
}

// =============================================================================
// UTILITY: prompt a positive double
// =============================================================================
static double promptPositiveDouble(const std::string& msg, double defaultVal) {
    double v = -1.0;
    while(v<=0.0){
        std::cout<<msg<<" [default "<<defaultVal<<"]: ";
        std::string line; std::getline(std::cin,line);
        if(line.empty()||line.find_first_not_of(" \t\r\n")==std::string::npos)
            return defaultVal;
        try{ v=std::stod(line); }
        catch(...){ v=-1.0; }
        if(v<=0.0) std::cout<<"  [ERROR] Please enter a positive number."<<std::endl;
    }
    return v;
}

// =============================================================================
// UTILITY: prompt an integer in range
// =============================================================================
static int promptInt(const std::string& msg, int defaultVal, int lo, int hi) {
    int v = -1;
    while(v<lo||v>hi){
        std::cout<<msg<<" [default "<<defaultVal<<", range "<<lo<<"-"<<hi<<"]: ";
        std::string line; std::getline(std::cin,line);
        if(line.empty()||line.find_first_not_of(" \t\r\n")==std::string::npos)
            return defaultVal;
        try{ v=std::stoi(line); }
        catch(...){ v=-1; }
        if(v<lo||v>hi) std::cout<<"  [ERROR] Please enter a value between "<<lo<<" and "<<hi<<"."<<std::endl;
    }
    return v;
}

// =============================================================================
// FIT FUNCTIONS
// =============================================================================

// Sum of N Gaussians — par[0]=N, then (A,mu,sigma) triplets
Double_t multiGauss(Double_t *x, Double_t *par) {
    int n=(int)par[0]; double s=0;
    for(int i=0;i<n;i++) s+=par[1+3*i]*TMath::Gaus(x[0],par[2+3*i],par[3+3*i]);
    return s;
}

// Asymmetric q-Gaussian (best for SiPM timing slices)
static Double_t qGaussAsym(Double_t* x, Double_t* par) {
    double xx=x[0],A=par[0],mu=par[1],sig=par[2],q1=par[3],q2=par[4];
    if(sig<=0) return 0;
    const double eps=1e-6;
    if(xx<=mu){
        if(TMath::Abs(q1-1.0)<eps) return A*TMath::Exp(-0.5*TMath::Power((xx-mu)/sig,2));
        double arg=1-(1-q1)*(1.0/(3-q1))*TMath::Power((xx-mu)/sig,2);
        return (arg<=0)?0:A*TMath::Power(arg,1.0/(1-q1));
    } else {
        double arg=1-(1-q2)*(1.0/(3-q2))*TMath::Power((xx-mu)/sig,2);
        return (arg<=0)?0:A*TMath::Power(arg,1.0/(1-q2));
    }
}

// ToT vs p.e. non-linear calibration: ToT = p0 + p1*log(1 + p2*n)
// Physically motivated: large signals saturate faster
static Double_t totCalibFunc(Double_t* x, Double_t* par) {
    double n=x[0];
    if(n<=0) return par[0];
    return par[0]+par[1]*TMath::Log(1.0+par[2]*n);
}

// Time-walk model: t_lead = t0 + A / sqrt(ToT)
static Double_t timeWalkFunc(Double_t* x, Double_t* par) {
    return par[0]+par[1]/TMath::Sqrt(x[0]+1e-6);
}

// Resolution model: sigma = sqrt(sigma_stoc^2/n + sigma_const^2)
static TF1* makeSigmaFit(const char* name, double xlo, double xhi) {
    TF1* f=new TF1(name,"sqrt(([0]*[0])/x + [1]*[1])",xlo,xhi);
    f->SetParameters(0.2,0.05);
    f->SetParLimits(0,0.0,5.0);
    f->SetParLimits(1,0.0,5.0);
    f->SetParNames("#sigma_{stoc}","#sigma_{const}");
    return f;
}

// =============================================================================
// ADAPTIVE q-GAUSSIAN FIT ON A SLICE
// Returns {mean, mean_err, sigma, sigma_err, TF1*}
// =============================================================================
static std::tuple<double,double,double,double,TF1*>
fitSlice(TH1D* h, const std::string& fname, double maxWin=2.0) {
    if(!h||h->GetEntries()<20) return {0,0,0,0,nullptr};
    int maxBin=h->GetMaximumBin();
    double center=h->GetBinCenter(maxBin);
    double halfMax=h->GetMaximum()*0.5;
    double bw=h->GetBinWidth(1);
    int lB=maxBin,rB=maxBin;
    while(lB>1&&h->GetBinContent(lB)>halfMax) lB--;
    while(rB<h->GetNbinsX()&&h->GetBinContent(rB)>halfMax) rB++;
    double sigEst=std::max((rB-lB)*bw/2.35,2.0*bw);
    double fitWin=std::min(3.0*sigEst,maxWin);
    TF1* f=new TF1(fname.c_str(),qGaussAsym,center-fitWin,center+fitWin,5);
    f->SetParameters(h->GetMaximum(),center,sigEst,1.0,1.3);
    f->FixParameter(3,1.0);
    f->SetNpx(1000);
    h->GetXaxis()->SetRangeUser(center-5*sigEst,center+5*sigEst);
    h->Fit(f,"IMREQ");
    return {f->GetParameter(1),f->GetParError(1),
            std::abs(f->GetParameter(2)),f->GetParError(2),f};
}

// =============================================================================
// QUANTISE to TDC bin
// =============================================================================
static double quantise(double t, double lsb) {
    return TMath::Floor(t/lsb+0.5)*lsb;
}

// =============================================================================
// DATA LOADERS
// =============================================================================

// --- ROOT loader ---
static bool loadROOT(const std::string& path,
                     double lsb_ns,
                     std::vector<AlcorEvent>& events)
{
    TFile* f=TFile::Open(path.c_str(),"READ");
    if(!f||f->IsZombie()){
        std::cerr<<"[ERROR] Cannot open ROOT file: "<<path<<std::endl;
        return false;
    }

    // Try to find TTree — accept common names
    TTree* tree=nullptr;
    for(const char* nm:{"alcor","data","hits","events","ch1","tree"}){
        tree=(TTree*)f->Get(nm);
        if(tree) break;
    }
    if(!tree){
        std::cerr<<"[ERROR] No recognisable TTree found in "<<path<<std::endl;
        f->Close(); return false;
    }

    // Detect branch names flexibly
    std::string bLead="", bTrail="", bLaser="";
    TObjArray* branches=tree->GetListOfBranches();
    for(int i=0;i<branches->GetEntries();i++){
        std::string bn=branches->At(i)->GetName();
        std::string bnl=bn;
        std::transform(bnl.begin(),bnl.end(),bnl.begin(),::tolower);
        if(bLead.empty()&&(bnl.find("lead")!=std::string::npos||bnl=="t_lead"||bnl=="tlead"||bnl=="t1"))
            bLead=bn;
        if(bTrail.empty()&&(bnl.find("trail")!=std::string::npos||bnl=="t_trail"||bnl=="ttrail"||bnl=="t2"))
            bTrail=bn;
        if(bLaser.empty()&&(bnl.find("laser")!=std::string::npos||bnl.find("ref")!=std::string::npos||bnl=="tref"))
            bLaser=bn;
    }

    if(bLead.empty()||bTrail.empty()){
        std::cerr<<"[ERROR] Could not find leading/trailing edge branches."<<std::endl;
        std::cerr<<"  Available branches: ";
        for(int i=0;i<branches->GetEntries();i++)
            std::cerr<<branches->At(i)->GetName()<<" ";
        std::cerr<<std::endl;
        f->Close(); return false;
    }

    std::cout<<"[INFO] Using branches: lead='"<<bLead
             <<"'  trail='"<<bTrail<<"'";
    if(!bLaser.empty()) std::cout<<"  laser='"<<bLaser<<"'";
    std::cout<<std::endl;

    Double_t vLead=0,vTrail=0,vLaser=-999;
    tree->SetBranchAddress(bLead.c_str(),&vLead);
    tree->SetBranchAddress(bTrail.c_str(),&vTrail);
    if(!bLaser.empty()) tree->SetBranchAddress(bLaser.c_str(),&vLaser);

    Long64_t n=tree->GetEntries();
    events.reserve(n);
    for(Long64_t i=0;i<n;i++){
        tree->GetEntry(i);
        double tot_raw=vTrail-vLead;
        if(tot_raw<=0) continue;
        AlcorEvent ev;
        ev.tot_raw=tot_raw;
        ev.t_lead =vLead;
        ev.t_laser=bLaser.empty()?-999.0:vLaser;
        events.push_back(ev);
    }
    f->Close();
    std::cout<<"[INFO] Loaded "<<events.size()<<" valid events from ROOT."<<std::endl;
    return true;
}

// --- CSV loader ---
static bool loadCSV(const std::string& path,
                    std::vector<AlcorEvent>& events)
{
    std::ifstream fin(path);
    if(!fin.is_open()){
        std::cerr<<"[ERROR] Cannot open CSV: "<<path<<std::endl;
        return false;
    }
    std::string line;
    std::getline(fin,line); // skip header
    int nLoaded=0;
    while(std::getline(fin,line)){
        if(line.empty()) continue;
        std::replace(line.begin(),line.end(),',',' ');
        std::istringstream ss(line);
        double tl,tt,tref=-999;
        ss>>tl>>tt;
        if(ss.fail()) continue;
        ss>>tref;
        double tot=tt-tl;
        if(tot<=0) continue;
        AlcorEvent ev; ev.tot_raw=tot; ev.t_lead=tl; ev.t_laser=tref;
        events.push_back(ev);
        nLoaded++;
    }
    std::cout<<"[INFO] Loaded "<<nLoaded<<" valid events from CSV."<<std::endl;
    return nLoaded>0;
}

// --- Simulation ---
static void simulate(double lsb_ns, std::vector<AlcorEvent>& events) {
    TRandom3 rng(42);
    events.reserve(SIM_N_EVENTS);
    for(int ev=0;ev<SIM_N_EVENTS;ev++){
        int npe=rng.Poisson(SIM_MU_LASER);
        if(npe<1) npe=1;
        if(npe>SIM_NPE_MAX) npe=SIM_NPE_MAX;
        double tot_true=SIM_TOT_BASE+SIM_TOT_SLOPE*(npe-1)
                        +rng.Gaus(0,SIM_TOT_SIGMA*TMath::Sqrt((double)npe));
        if(tot_true<2.0) tot_true=2.0;
        double sig_t=TMath::Sqrt(SIM_SIGMA_STOC*SIM_SIGMA_STOC/npe
                                 +SIM_SIGMA_CONST*SIM_SIGMA_CONST);
        double t_sipm=SIM_T0+rng.Gaus(0,sig_t);
        double t_lead_true=t_sipm+SIM_TW_A/TMath::Sqrt(tot_true);
        double t_laser=SIM_T0+rng.Gaus(0,SIM_LASER_JIT);
        AlcorEvent e;
        e.tot_raw =quantise(tot_true,    SIM_TDC_BIN_NS)/lsb_ns;
        e.t_lead  =quantise(t_lead_true, SIM_TDC_BIN_NS)/lsb_ns;
        e.t_laser =quantise(t_laser,     SIM_TDC_BIN_NS)/lsb_ns;
        events.push_back(e);
    }
    std::cout<<"[INFO] Simulated "<<events.size()<<" events."<<std::endl;
}

// =============================================================================
// CANVAS HELPER: standard margins
// =============================================================================
static TCanvas* makeCanvas(const char* name, const char* title,
                            int w=900, int h=700,
                            bool logY=false, bool logZ=false) {
    TCanvas* c=new TCanvas(name,title,w,h);
    c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12); c->SetTopMargin(0.07);
    if(logY) c->SetLogy();
    if(logZ) c->SetLogz();
    c->SetGrid();
    return c;
}

// =============================================================================
// MAIN
// =============================================================================
void alcor_tot_analysis() {

    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetTitleSize(0.045,"XY");
    gStyle->SetLabelSize(0.040,"XY");
    TVirtualFitter::SetDefaultFitter("Minuit2");

    gSystem->Exec("mkdir -p plots");

    // =========================================================================
    // INTERACTIVE STARTUP
    // =========================================================================
    printBox({"ALCOR SiPM Calibration -- ToT Analysis",
              "Official code language: English"});

    // --- LSB ---
    std::cout<<"\n  [STEP 1/4]  TDC conversion factor"<<std::endl;
    std::cout<<"  If ALCOR delivers raw ticks, enter the LSB in ns."<<std::endl;
    std::cout<<"  Typical ALCOR LSB = 0.05 ns (50 ps)."<<std::endl;
    std::cout<<"  If data are already in ns, enter 1.0."<<std::endl;
    double lsb_ns = promptPositiveDouble("  LSB [ns/tick]", 1.0);

    // --- number of p.e. peaks to fit ---
    std::cout<<"\n  [STEP 2/4]  Expected number of p.e. peaks in spectrum"<<std::endl;
    int nPeaks = promptInt("  Number of peaks", 8, 2, 15);

    // --- laser reference available? ---
    std::cout<<"\n  [STEP 3/4]  Laser / external reference"<<std::endl;
    bool hasLaser = (askYesNo("  Do you have a laser reference channel? [y/n]: ")=='y');

    // --- data source ---
    std::cout<<"\n  [STEP 4/4]  Data source"<<std::endl;
    printBox({"[1] ROOT file  (alcor_data.root)",
              "[2] CSV file   (alcor_data.csv)",
              "[3] Simulation (built-in Monte Carlo)"});
    int srcChoice=0;
    while(srcChoice<1||srcChoice>3){
        std::cout<<"  Your choice [1/2/3]: ";
        std::string line; std::getline(std::cin,line);
        try{ srcChoice=std::stoi(line); }catch(...){ srcChoice=0; }
        if(srcChoice<1||srcChoice>3)
            std::cout<<"  [ERROR] Please enter 1, 2, or 3."<<std::endl;
    }

    // =========================================================================
    // LOAD / SIMULATE DATA
    // =========================================================================
    std::vector<AlcorEvent> rawEvents;

    if(srcChoice==1){
        std::cout<<"  Enter ROOT file path [default: alcor_data.root]: ";
        std::string p; std::getline(std::cin,p);
        if(p.empty()||p.find_first_not_of(" \t")==std::string::npos) p="../../data/alcor_data.root";
        if(!loadROOT(p,lsb_ns,rawEvents)) return;
    } else if(srcChoice==2){
        std::cout<<"  Enter CSV file path [default: alcor_data.csv]: ";
        std::string p; std::getline(std::cin,p);
        if(p.empty()||p.find_first_not_of(" \t")==std::string::npos) p="../../data/alcor_data.csv";
        if(!loadCSV(p,rawEvents)) return;
    } else {
        std::cout<<"[INFO] Running in simulation mode..."<<std::endl;
        simulate(lsb_ns,rawEvents);
        hasLaser=true; // simulation always includes laser reference
    }

    if(rawEvents.empty()){
        std::cerr<<"[ERROR] No events loaded. Exiting."<<std::endl;
        return;
    }

    // =========================================================================
    // CONVERT RAW -> PHYSICAL UNITS
    // =========================================================================
    // tot_ns = tot_raw * lsb_ns
    // dt_ns  = (t_lead - t_laser) * lsb_ns  [only if laser available]
    int N = (int)rawEvents.size();
    std::vector<double> vToT(N), vTlead(N), vDt(N);
    bool anyLaser=false;

    for(int i=0;i<N;i++){
        vToT[i]   = rawEvents[i].tot_raw * lsb_ns;
        vTlead[i] = rawEvents[i].t_lead  * lsb_ns;
        if(hasLaser && rawEvents[i].t_laser>-998){
            vDt[i] = (rawEvents[i].t_lead - rawEvents[i].t_laser)*lsb_ns;
            anyLaser=true;
        } else {
            vDt[i] = -999.0;
        }
    }
    if(!anyLaser) hasLaser=false;

    double totMin=*std::min_element(vToT.begin(),vToT.end());
    double totMax=*std::max_element(vToT.begin(),vToT.end());
    double totRange=totMax-totMin;

    std::cout<<"\n[INFO] ToT range: "<<totMin<<" -- "<<totMax<<" ns"<<std::endl;
    std::cout<<"[INFO] Events with laser ref: "<<(hasLaser?"yes":"no")<<std::endl;

    // =========================================================================
    // PASS 1: ToT SPECTRUM + PEAK FINDING + CALIBRATION
    // =========================================================================
    int nBinsToT=400;
    double totLo=std::max(0.0,totMin-totRange*0.05);
    double totHi=totMax+totRange*0.05;

    TH1D *hToT = new TH1D("hToT",
        "ALCOR ToT Spectrum;ToT (ns);counts",
        nBinsToT, totLo, totHi);
    for(int i=0;i<N;i++) hToT->Fill(vToT[i]);

    // Peak search — try progressively lower thresholds until nPeaks found
    TSpectrum spec(nPeaks+5);
    std::vector<double> peakPos;
    for(double thr : {0.003, 0.005, 0.010, 0.020, 0.050}){
        int nFound=spec.Search(hToT, 2, "goff", thr);
        double *xPk=spec.GetPositionX();
        peakPos.assign(xPk, xPk+nFound);
        std::sort(peakPos.begin(), peakPos.end());
        // Remove peaks too close to each other (< 40% of expected spacing)
        double expSpacing=(totMax-totMin)/(nPeaks+1);
        std::vector<double> filtered;
        for(int k=0;k<(int)peakPos.size();k++){
            if(filtered.empty()||peakPos[k]-filtered.back()>expSpacing*0.4)
                filtered.push_back(peakPos[k]);
        }
        peakPos=filtered;
        if((int)peakPos.size()>=nPeaks) break;
        std::cout<<"[INFO] TSpectrum thr="<<thr<<" found "<<peakPos.size()
                 <<" peaks, retrying..."<<std::endl;
    }
    if((int)peakPos.size()>nPeaks) peakPos.resize(nPeaks);
    int nFit=(int)peakPos.size();
    std::cout<<"[INFO] Using "<<nFit<<" peaks for calibration fit."<<std::endl;

    // Estimate inter-peak spacing for parameter initialisation
    double spacing=(nFit>1)?(peakPos.back()-peakPos.front())/(nFit-1):
                            (totMax-totMin)/(nPeaks+1);

    // Multi-Gaussian fit on spectrum
    TF1 *fSpec=new TF1("fSpec",multiGauss,totLo,totHi,1+3*nFit);
    fSpec->FixParameter(0,nFit);
    fSpec->SetNpx(2000);
    for(int i=0;i<nFit;i++){
        double xest  = peakPos[i];
        double aest  = hToT->GetBinContent(hToT->FindBin(xest));
        double sigest= spacing*0.12*TMath::Sqrt((double)(i+1));
        fSpec->SetParameter(1+3*i, aest);
        fSpec->SetParLimits(1+3*i, 0.0, aest*3.0);
        fSpec->SetParameter(2+3*i, xest);
        fSpec->SetParLimits(2+3*i, xest-spacing*0.45, xest+spacing*0.45);
        fSpec->SetParameter(3+3*i, sigest);
        fSpec->SetParLimits(3+3*i, lsb_ns*2, spacing*0.55);
    }
    TVirtualFitter::SetDefaultFitter("Minuit2");
    hToT->Fit(fSpec,"RQMBL");

    // Collect fitted peak positions and sigmas
    std::vector<double> vToT_mean(nFit),vToT_sigma(nFit),vPE(nFit);
    for(int i=0;i<nFit;i++){
        vToT_mean[i] =fSpec->GetParameter(2+3*i);
        vToT_sigma[i]=std::abs(fSpec->GetParameter(3+3*i));
        vPE[i]=(double)(i+1);
        std::cout<<Form("  Peak %d p.e.: ToT = %.3f +/- %.3f ns",
            i+1, vToT_mean[i], vToT_sigma[i])<<std::endl;
    }

    // =========================================================================
    // CANVAS 01a: ToT spectrum — log scale
    // =========================================================================
    TCanvas *cToTlog=makeCanvas("cToTlog","ToT Spectrum (log)",900,700,true);
    hToT->SetLineColor(kTeal+2); hToT->SetFillColorAlpha(kTeal,0.2);
    hToT->SetMinimum(0.5);
    hToT->Draw("HIST");
    fSpec->SetLineColor(kRed+1); fSpec->SetLineWidth(2); fSpec->Draw("SAME");
    // Mark peaks
    for(int i=0;i<nFit;i++){
        TLine *l=new TLine(vToT_mean[i],hToT->GetMinimum(),
                           vToT_mean[i],hToT->GetBinContent(hToT->FindBin(vToT_mean[i]))*0.8);
        l->SetLineColor(kOrange+7); l->SetLineStyle(2); l->Draw();
    }
    cToTlog->SaveAs("plots/alcor_01a_tot_spectrum_log.png");

    // =========================================================================
    // CANVAS 01b: ToT spectrum — linear scale
    // =========================================================================
    TCanvas *cToTlin=makeCanvas("cToTlin","ToT Spectrum (linear)",900,700);
    hToT->Draw("HIST");
    fSpec->Draw("SAME");
    cToTlin->SaveAs("plots/alcor_01b_tot_spectrum_linear.png");

    // =========================================================================
    // CANVAS 02: ToT vs p.e. calibration curve (non-linear fit)
    // =========================================================================
    std::vector<double> vPE_err(nFit,0.0),vToT_err(nFit);
    for(int i=0;i<nFit;i++) vToT_err[i]=vToT_sigma[i];

    TGraphErrors *grCal=new TGraphErrors(nFit,&vPE[0],&vToT_mean[0],
                                          &vPE_err[0],&vToT_err[0]);
    grCal->SetTitle("ToT Calibration Curve;photoelectrons (n);ToT_{peak} (ns)");
    grCal->SetMarkerStyle(21); grCal->SetMarkerSize(1.3);
    grCal->SetMarkerColor(kAzure+2); grCal->SetLineColor(kAzure+2);

    // Linear fit (primary — physically correct for ALCOR ToT model)
    TF1 *fCalLin=new TF1("fCalLin","pol1",0.5,nFit+0.5);
    fCalLin->SetLineColor(kRed+1); fCalLin->SetLineWidth(2);
    grCal->Fit(fCalLin,"RQ");
    double gain_ns_pe = fCalLin->GetParameter(1);
    double offset_ns  = fCalLin->GetParameter(0);

    // Non-linear fit (kept for comparison — useful with real ALCOR data)
    TF1 *fCal=new TF1("fCal",totCalibFunc,0.5,nFit+0.5,3);
    fCal->SetParameters(offset_ns, gain_ns_pe*(double)nFit, 1.0/(double)nFit);
    fCal->SetParLimits(0, offset_ns-20, offset_ns+20);
    fCal->SetParLimits(1, 0.0, gain_ns_pe*(double)nFit*5.0);
    fCal->SetParLimits(2, 1e-4, 10.0);
    fCal->SetParNames("ToT_{0}","A","#lambda");
    fCal->SetLineColor(kGray+2); fCal->SetLineStyle(2); fCal->SetLineWidth(2);
    grCal->Fit(fCal,"RQB+");

    TCanvas *cCal=makeCanvas("cCal","ToT Calibration",900,700);
    grCal->GetYaxis()->SetTitleOffset(1.5);
    grCal->Draw("AP");
    fCalLin->Draw("SAME"); fCal->Draw("SAME");
    TLegend *legCal=new TLegend(0.18,0.55,0.65,0.88);
    legCal->SetBorderSize(1); legCal->SetTextSize(0.030);
    legCal->AddEntry(grCal,"ToT peak positions","pe");
    legCal->AddEntry(fCalLin,"Linear fit  (primary)","l");
    legCal->AddEntry(fCal,"Log fit  (reference)","l");
    legCal->AddEntry((TObject*)0,
        Form("Gain   = %.3f #pm %.3f ns/p.e.",
             fCalLin->GetParameter(1),fCalLin->GetParError(1)),"");
    legCal->AddEntry((TObject*)0,
        Form("Offset = %.3f #pm %.3f ns",
             fCalLin->GetParameter(0),fCalLin->GetParError(0)),"");
    legCal->Draw();
    cCal->SaveAs("plots/alcor_02_tot_calibration.png");

    // =========================================================================
    // FROM HERE: timing analysis (only if laser reference available)
    // =========================================================================
    if(!hasLaser){
        std::cout<<"\n[INFO] No laser reference — skipping timing plots 03-15."<<std::endl;
        printBox({"Analysis complete (no laser reference).",
                  "Plots 01-02 saved to plots/"});
        gSystem->Exec("find . -maxdepth 1 -name '*.so' -o -name '*.d' -o -name '*.pcm' | xargs rm -f 2>/dev/null");
        return;
    }

    // Filter events with valid laser reference
    std::vector<double> vToT_ev, vDt_ev, vTlead_ev;
    for(int i=0;i<N;i++){
        if(vDt[i]<-998) continue;
        vToT_ev.push_back(vToT[i]);
        vDt_ev.push_back(vDt[i]);
        vTlead_ev.push_back(vTlead[i]);
    }
    int Nev=(int)vToT_ev.size();

    double dtMin=*std::min_element(vDt_ev.begin(),vDt_ev.end());
    double dtMax=*std::max_element(vDt_ev.begin(),vDt_ev.end());
    double dtC=(dtMin+dtMax)*0.5;
    double dtHW=std::max((dtMax-dtMin)*0.6, 5.0);

    // =========================================================================
    // CANVAS 03: raw Delta-t spectrum
    // =========================================================================
    TH1D *hDtRaw=new TH1D("hDtRaw",
        "Raw #Delta-t Spectrum;#Delta-t (ns);counts",
        500, dtC-dtHW, dtC+dtHW);
    for(int i=0;i<Nev;i++) hDtRaw->Fill(vDt_ev[i]);

    TCanvas *cDtRaw=makeCanvas("cDtRaw","Raw Delta-t",900,700,true);
    hDtRaw->SetLineColor(kTeal+3); hDtRaw->SetFillColorAlpha(kTeal,0.2);
    hDtRaw->SetMinimum(0.5);
    hDtRaw->Draw("HIST");
    cDtRaw->SaveAs("plots/alcor_03_dt_raw.png");

    // =========================================================================
    // CANVAS 04: 2D map Delta-t vs ToT (time walk)
    // =========================================================================
    TH2D *h2D=new TH2D("h2D",
        "Time-Walk Map;ToT (ns);#Delta-t (ns)",
        200, totLo, totHi, 400, dtC-dtHW, dtC+dtHW);
    for(int i=0;i<Nev;i++) h2D->Fill(vToT_ev[i],vDt_ev[i]);

    TCanvas *c2D=makeCanvas("c2D","2D Time-Walk Map",1000,700,false,true);
    c2D->SetRightMargin(0.13);
    h2D->SetContour(100); h2D->GetYaxis()->SetTitleOffset(1.5);
    h2D->Draw("COLZ");
    cDtRaw->SaveAs("plots/alcor_03_dt_raw.png");
    c2D->SaveAs("plots/alcor_04_timewalk_map.png");

    // =========================================================================
    // CANVAS 05: time-walk correction fit t_lead vs ToT
    // Profile histogram for robust fit
    // =========================================================================
    TH2D *h2Dwalk=new TH2D("h2Dwalk","",200,totLo,totHi,400,
                            *std::min_element(vTlead_ev.begin(),vTlead_ev.end())-1,
                            *std::max_element(vTlead_ev.begin(),vTlead_ev.end())+1);
    for(int i=0;i<Nev;i++) h2Dwalk->Fill(vToT_ev[i],vTlead_ev[i]);
    TProfile *pWalk=(TProfile*)h2Dwalk->ProfileX("pWalk");

    TF1 *fWalk=new TF1("fWalk",timeWalkFunc,totLo,totHi,2);
    fWalk->SetParameters(pWalk->GetBinContent(pWalk->GetNbinsX()/2), 5.0);
    fWalk->SetParNames("t_{0}","A_{TW}");
    fWalk->SetLineColor(kRed+1); fWalk->SetLineWidth(2);
    pWalk->Fit(fWalk,"RQ");

    TCanvas *cWalk=makeCanvas("cWalk","Time-Walk Correction",900,700);
    pWalk->SetTitle("Time-Walk Profile: t_{lead} vs ToT;ToT (ns);t_{lead} (ns)");
    pWalk->SetMarkerStyle(20); pWalk->SetMarkerSize(0.6);
    pWalk->SetMarkerColor(kAzure+2);
    pWalk->GetYaxis()->SetTitleOffset(1.5);
    pWalk->Draw("PE");
    fWalk->Draw("SAME");
    TLegend *legWalk=new TLegend(0.55,0.75,0.92,0.88);
    legWalk->SetBorderSize(1); legWalk->SetTextSize(0.030);
    legWalk->AddEntry(pWalk,"Profile data","pe");
    legWalk->AddEntry(fWalk,Form("t_{0}+A/#sqrt{ToT},  A=%.3f ns",fWalk->GetParameter(1)),"l");
    legWalk->Draw();
    cWalk->SaveAs("plots/alcor_05_timewalk_correction.png");

    // Apply time-walk correction: dt_corr = dt - (f_walk(ToT) - t0)
    double t0_tw=fWalk->GetParameter(0);
    std::vector<double> vDt_corr(Nev);
    for(int i=0;i<Nev;i++)
        vDt_corr[i]=vDt_ev[i]-(fWalk->Eval(vToT_ev[i])-t0_tw);

    double dtcC=0; for(int i=0;i<Nev;i++) dtcC+=vDt_corr[i]; dtcC/=Nev;
    // Re-estimate range after correction
    double dtcMin=*std::min_element(vDt_corr.begin(),vDt_corr.end());
    double dtcMax=*std::max_element(vDt_corr.begin(),vDt_corr.end());
    double dtcHW=std::max((dtcMax-dtcMin)*0.6,3.0);

    // =========================================================================
    // CANVAS 06: corrected Delta-t spectrum
    // =========================================================================
    TH1D *hDtCorr=new TH1D("hDtCorr",
        "Corrected #Delta-t Spectrum;#Delta-t_{corr} (ns);counts",
        500, dtcC-dtcHW, dtcC+dtcHW);
    for(int i=0;i<Nev;i++) hDtCorr->Fill(vDt_corr[i]);

    TCanvas *cDtCorr=makeCanvas("cDtCorr","Corrected Delta-t",900,700,true);
    hDtCorr->SetLineColor(kOrange+7); hDtCorr->SetFillColorAlpha(kOrange,0.2);
    hDtCorr->SetMinimum(0.5);
    hDtCorr->Draw("HIST");
    cDtCorr->SaveAs("plots/alcor_06_dt_corrected.png");

    // =========================================================================
    // CANVAS 07: before vs after overlay
    // =========================================================================
    TCanvas *cOver=makeCanvas("cOver","Before vs After Correction",900,700,true);
    TH1D *hRawC=(TH1D*)hDtRaw->Clone("hRawC");
    TH1D *hCorrC=(TH1D*)hDtCorr->Clone("hCorrC");
    hRawC->SetLineColor(kTeal+3); hRawC->SetFillStyle(0); hRawC->SetLineWidth(2);
    hCorrC->SetLineColor(kOrange+7); hCorrC->SetFillStyle(0); hCorrC->SetLineWidth(2);
    // Align to common x-range
    double xlo_ov=std::min(hRawC->GetXaxis()->GetXmin(),hCorrC->GetXaxis()->GetXmin());
    double xhi_ov=std::max(hRawC->GetXaxis()->GetXmax(),hCorrC->GetXaxis()->GetXmax());
    hRawC->GetXaxis()->SetRangeUser(xlo_ov,xhi_ov);
    hRawC->SetTitle("Before vs After Time-Walk Correction;#Delta-t (ns);counts");
    hRawC->SetMinimum(0.5);
    hRawC->GetYaxis()->SetTitleOffset(1.5);
    hRawC->Draw("HIST");
    hCorrC->Draw("HIST SAME");
    TLegend *legOver=new TLegend(0.55,0.75,0.90,0.88);
    legOver->SetBorderSize(1); legOver->SetTextSize(0.028);
    legOver->AddEntry(hRawC,"Raw","l");
    legOver->AddEntry(hCorrC,"Time-walk corrected","l");
    legOver->Draw();
    cOver->SaveAs("plots/alcor_07_before_after.png");

    // =========================================================================
    // PASS 2: build 2D maps (ToT vs corrected dt) for slice analysis
    // =========================================================================
    TH2D *h2DCorr=new TH2D("h2DCorr",
        "Corrected Map;ToT (ns);#Delta-t_{corr} (ns)",
        200, totLo, totHi, 400, dtcC-dtcHW, dtcC+dtcHW);
    for(int i=0;i<Nev;i++) h2DCorr->Fill(vToT_ev[i],vDt_corr[i]);

    // =========================================================================
    // CANVAS 08: timing slices per p.e. band
    // =========================================================================
    int nCols=4, nRows=(nFit+nCols-1)/nCols;
    TCanvas *cSlices=new TCanvas("cSlices","Timing Slices per p.e. Band",
                                  1200,nRows*350);
    cSlices->Divide(nCols,nRows);

    std::vector<double> vBand,vBandZ,vSigRaw,vESigRaw,vSigCorr,vESigCorr;

    for(int n=0;n<nFit;n++){
        cSlices->cd(n+1);
        gPad->SetLeftMargin(0.20); gPad->SetLogy();

        double tlo=vToT_mean[n]-1.5*vToT_sigma[n];
        double thi=vToT_mean[n]+1.5*vToT_sigma[n];
        int bLo=h2D->GetXaxis()->FindBin(tlo);
        int bHi=h2D->GetXaxis()->FindBin(thi);

        TH1D *hSraw =h2D->ProjectionY(Form("hSraw_%d",n),bLo,bHi);
        TH1D *hScorr=h2DCorr->ProjectionY(Form("hScorr_%d",n),
                         h2DCorr->GetXaxis()->FindBin(tlo),
                         h2DCorr->GetXaxis()->FindBin(thi));

        hSraw->SetTitle(Form("%d p.e. band;#Delta-t (ns);counts",n+1));
        hSraw->SetLineColor(kTeal+3);
        hSraw->SetFillColorAlpha(kTeal,0.15);
        hSraw->GetYaxis()->SetTitleOffset(2.0);
        hSraw->GetYaxis()->SetTitleSize(0.07);
        hSraw->GetXaxis()->SetTitleSize(0.07);
        hSraw->GetXaxis()->SetLabelSize(0.07);
        hSraw->GetYaxis()->SetLabelSize(0.07);

        double slYmax=TMath::Max(hSraw->GetMaximum(),hScorr->GetMaximum())*2.0;
        if(slYmax<1) slYmax=10;
        hSraw->SetMinimum(0.5); hSraw->SetMaximum(slYmax);
        hSraw->Draw("HIST");

        hScorr->SetLineColor(kOrange+7); hScorr->SetFillStyle(0);
        hScorr->SetLineWidth(2); hScorr->Draw("HIST SAME");

        // Fit corrected slice
        if(hScorr->GetEntries()>20){
            auto res=fitSlice(hScorr,Form("fSlC_%d",n),2.0);
            if(std::get<4>(res)){
                std::get<4>(res)->SetLineColor(kRed+1);
                std::get<4>(res)->SetLineWidth(2);
                std::get<4>(res)->Draw("SAME");
                vBand.push_back((double)(n+1)); vBandZ.push_back(0.0);
                vSigCorr.push_back(std::get<2>(res));
                vESigCorr.push_back(std::get<3>(res));
            }
        }
        // Fit raw slice (invisible, just for resolution extraction)
        if(hSraw->GetEntries()>20){
            auto res=fitSlice(hSraw,Form("fSlR_%d",n),2.0);
            if(std::get<4>(res)){
                vSigRaw.push_back(std::get<2>(res));
                vESigRaw.push_back(std::get<3>(res));
            }
        }
    }
    cSlices->SaveAs("plots/alcor_08_slices.png");

    // =========================================================================
    // CANVAS 09: resolution before correction
    // =========================================================================
    int nRes=(int)std::min({vBand.size(),vSigRaw.size(),vBandZ.size()});
    if(nRes>0){
        TCanvas *cRes1=makeCanvas("cRes1","Resolution (raw)",900,700);
        TGraphErrors *grRes1=new TGraphErrors(nRes,&vBand[0],&vSigRaw[0],
                                               &vBandZ[0],&vESigRaw[0]);
        grRes1->SetTitle("Timing Resolution (raw);"
            "photoelectrons (n);#sigma_{#Delta-t} (ns)");
        grRes1->SetMarkerStyle(21); grRes1->SetMarkerSize(1.2);
        grRes1->SetMarkerColor(kTeal+3);
        grRes1->GetYaxis()->SetTitleOffset(1.5);
        TF1 *fSig1=makeSigmaFit("fSig1",0.5,nRes+0.5);
        fSig1->SetLineColor(kBlack); fSig1->SetLineStyle(2);
        grRes1->Fit(fSig1,"RQB"); grRes1->Draw("AP");
        TLegend *lR1=new TLegend(0.45,0.65,0.90,0.88);
        lR1->SetBorderSize(1); lR1->SetTextSize(0.028);
        lR1->AddEntry(grRes1,"Raw resolution","pe");
        lR1->AddEntry(fSig1,"#sqrt{#sigma_{stoc}^{2}/n + #sigma_{const}^{2}}","l");
        lR1->AddEntry((TObject*)0,Form("#sigma_{stoc}=%.4f #pm %.4f ns",
            fSig1->GetParameter(0),fSig1->GetParError(0)),"");
        lR1->AddEntry((TObject*)0,Form("#sigma_{const}=%.4f #pm %.4f ns",
            fSig1->GetParameter(1),fSig1->GetParError(1)),"");
        lR1->Draw();
        cRes1->SaveAs("plots/alcor_09_resolution_raw.png");
    }

    // =========================================================================
    // CANVAS 10: resolution after correction
    // =========================================================================
    int nResC=(int)std::min({vBand.size(),vSigCorr.size(),vBandZ.size()});
    if(nResC>0){
        TCanvas *cRes2=makeCanvas("cRes2","Resolution (corrected)",900,700);
        TGraphErrors *grRes2=new TGraphErrors(nResC,&vBand[0],&vSigCorr[0],
                                               &vBandZ[0],&vESigCorr[0]);
        grRes2->SetTitle("Timing Resolution (corrected);"
            "photoelectrons (n);#sigma_{#Delta-t} (ns)");
        grRes2->SetMarkerStyle(21); grRes2->SetMarkerSize(1.2);
        grRes2->SetMarkerColor(kOrange+7);
        grRes2->GetYaxis()->SetTitleOffset(1.5);
        TF1 *fSig2=makeSigmaFit("fSig2",0.5,nResC+0.5);
        fSig2->SetLineColor(kBlack); fSig2->SetLineStyle(2);
        grRes2->Fit(fSig2,"RQB"); grRes2->Draw("AP");
        TLegend *lR2=new TLegend(0.45,0.65,0.90,0.88);
        lR2->SetBorderSize(1); lR2->SetTextSize(0.028);
        lR2->AddEntry(grRes2,"Corrected resolution","pe");
        lR2->AddEntry(fSig2,"#sqrt{#sigma_{stoc}^{2}/n + #sigma_{const}^{2}}","l");
        lR2->AddEntry((TObject*)0,Form("#sigma_{stoc}=%.4f #pm %.4f ns",
            fSig2->GetParameter(0),fSig2->GetParError(0)),"");
        lR2->AddEntry((TObject*)0,Form("#sigma_{const}=%.4f #pm %.4f ns",
            fSig2->GetParameter(1),fSig2->GetParError(1)),"");
        lR2->Draw();
        cRes2->SaveAs("plots/alcor_10_resolution_corrected.png");
    }

    // =========================================================================
    // CANVAS 11: SNR — peak separation vs sigma per p.e. band
    // =========================================================================
    if(nFit>=2){
        TCanvas *cSNR=makeCanvas("cSNR","Peak SNR per p.e. Band",900,700);
        std::vector<double> vSNR,vSNRpe;
        for(int i=0;i<nFit-1;i++){
            double sep=vToT_mean[i+1]-vToT_mean[i];
            double sig=0.5*(vToT_sigma[i]+vToT_sigma[i+1]);
            if(sig>0){ vSNR.push_back(sep/sig); vSNRpe.push_back(i+1.5); }
        }
        TGraph *grSNR=new TGraph((int)vSNRpe.size(),&vSNRpe[0],&vSNR[0]);
        grSNR->SetTitle("Peak Separation / Sigma vs p.e.;"
            "p.e. band (midpoint);#Delta_{peak} / #sigma");
        grSNR->SetMarkerStyle(20); grSNR->SetMarkerSize(1.3);
        grSNR->SetMarkerColor(kViolet+2); grSNR->SetLineColor(kViolet+2);
        grSNR->GetYaxis()->SetTitleOffset(1.5);
        grSNR->Draw("APL");
        // Reference line at SNR=2 (minimum resolvable)
        TLine *lSNR2=new TLine(vSNRpe.front()-0.3,2.0,vSNRpe.back()+0.3,2.0);
        lSNR2->SetLineColor(kRed); lSNR2->SetLineStyle(2); lSNR2->Draw();
        TLatex *txtSNR=new TLatex(vSNRpe.back()+0.1,2.05,"SNR = 2");
        txtSNR->SetTextSize(0.030); txtSNR->SetTextColor(kRed); txtSNR->Draw();
        cSNR->SaveAs("plots/alcor_11_snr.png");
    }

    // =========================================================================
    // CANVAS 12: ToT stability — running mean over event index
    // =========================================================================
    {
        int block=std::max(1,Nev/200); // 200 points
        int nPts=Nev/block;
        std::vector<double> vIdx(nPts),vMeanTot(nPts);
        for(int i=0;i<nPts;i++){
            double sum=0;
            for(int j=i*block;j<(i+1)*block&&j<Nev;j++) sum+=vToT_ev[j];
            vIdx[i]=i*block; vMeanTot[i]=sum/block;
        }
        TGraph *grStab=new TGraph(nPts,&vIdx[0],&vMeanTot[0]);
        grStab->SetTitle("ToT Running Mean (Stability);"
            "event index;ToT_{mean} (ns)");
        grStab->SetLineColor(kAzure+2); grStab->SetLineWidth(1);
        grStab->GetYaxis()->SetTitleOffset(1.5);
        TCanvas *cStab=makeCanvas("cStab","ToT Stability",900,700);
        grStab->Draw("AL");
        // Reference: overall mean
        double totMeanAll=std::accumulate(vToT_ev.begin(),vToT_ev.end(),0.0)/Nev;
        TLine *lMean=new TLine(0,totMeanAll,vIdx.back(),totMeanAll);
        lMean->SetLineColor(kRed); lMean->SetLineStyle(2); lMean->Draw();
        cStab->SaveAs("plots/alcor_12_tot_stability.png");
    }

    // =========================================================================
    // CANVAS 13: inter-arrival time (dark rate proxy)
    // Uses t_lead of consecutive events — only meaningful if events are
    // time-ordered and represent continuous acquisition
    // =========================================================================
    {
        std::vector<double> vIAT;
        for(int i=1;i<Nev;i++){
            double dt_iat=(vTlead_ev[i]-vTlead_ev[i-1])*lsb_ns;
            if(dt_iat>0&&dt_iat<1e6) vIAT.push_back(dt_iat);
        }
        if(!vIAT.empty()){
            double iatMax=*std::max_element(vIAT.begin(),vIAT.end());
            TH1D *hIAT=new TH1D("hIAT",
                "Inter-Arrival Time;#Delta t_{arrival} (ns);counts",
                200, 0, std::min(iatMax*1.05,iatMax));
            for(double v:vIAT) hIAT->Fill(v);
            TCanvas *cIAT=makeCanvas("cIAT","Inter-Arrival Time",900,700,true);
            hIAT->SetLineColor(kSpring-1); hIAT->SetFillColorAlpha(kSpring,0.2);
            hIAT->SetMinimum(0.5); hIAT->Draw("HIST");
            // Exponential fit (Poisson process)
            TF1 *fExp=new TF1("fExp","[0]*exp(-x/[1])",
                hIAT->GetBinCenter(1), hIAT->GetBinCenter(hIAT->GetNbinsX()));
            fExp->SetParameters(hIAT->GetMaximum(), hIAT->GetMean());
            fExp->SetLineColor(kRed+1); fExp->SetLineWidth(2);
            hIAT->Fit(fExp,"RQ");
            fExp->Draw("SAME");
            TLegend *lIAT=new TLegend(0.55,0.75,0.90,0.88);
            lIAT->SetBorderSize(1); lIAT->SetTextSize(0.028);
            lIAT->AddEntry(hIAT,"Inter-arrival times","l");
            lIAT->AddEntry(fExp,Form("Exp fit: #tau = %.2f ns",fExp->GetParameter(1)),"l");
            lIAT->Draw();
            cIAT->SaveAs("plots/alcor_13_interarrival.png");
        }
    }

    // =========================================================================
    // CANVAS 14: gain linearity — ToT peak vs n p.e. (same as 02 but focused)
    // =========================================================================
    {
        TCanvas *cLin=makeCanvas("cLin","ToT Peak Linearity",900,700);
        TGraphErrors *grLin=new TGraphErrors(nFit,&vPE[0],&vToT_mean[0],
                                              &vPE_err[0],&vToT_err[0]);
        grLin->SetTitle("ToT Peak Linearity;"
            "photoelectrons (n);ToT_{peak} (ns)");
        grLin->SetMarkerStyle(21); grLin->SetMarkerSize(1.3);
        grLin->SetMarkerColor(kAzure+2);
        grLin->GetYaxis()->SetTitleOffset(1.5);
        // Fit both linear and non-linear, display residuals
        TF1 *fL=new TF1("fL2","pol1",0.5,nFit+0.5);
        fL->SetLineColor(kRed+1); fL->SetLineWidth(2);
        grLin->Fit(fL,"RQ");
        grLin->Draw("AP"); fL->Draw("SAME");
        TLegend *lLin=new TLegend(0.18,0.68,0.60,0.88);
        lLin->SetBorderSize(1); lLin->SetTextSize(0.030);
        lLin->AddEntry(grLin,"ToT peak positions","pe");
        lLin->AddEntry(fL,"Linear fit","l");
        lLin->AddEntry((TObject*)0,
            Form("Gain = %.3f #pm %.3f ns/p.e.",
                 fL->GetParameter(1),fL->GetParError(1)),"");
        lLin->AddEntry((TObject*)0,
            Form("Offset = %.3f #pm %.3f ns",
                 fL->GetParameter(0),fL->GetParError(0)),"");
        lLin->Draw();
        cLin->SaveAs("plots/alcor_14_linearity.png");
    }

    // =========================================================================
    // CANVAS 15: summary panel
    // =========================================================================
    {
        TCanvas *cSum=new TCanvas("cSum","Summary",1200,800);
        cSum->Divide(3,2);

        cSum->cd(1); hToT->Draw("HIST"); fSpec->Draw("SAME");
        cSum->cd(2); h2D->Draw("COLZ");
        cSum->cd(3); hDtRaw->SetLineColor(kTeal+3); hDtRaw->Draw("HIST");
                     hDtCorr->SetLineColor(kOrange+7); hDtCorr->Draw("HIST SAME");
        cSum->cd(4);
        {
            TGraphErrors *grS=new TGraphErrors(nFit,&vPE[0],&vToT_mean[0],
                                                &vPE_err[0],&vToT_err[0]);
            grS->SetTitle(";n p.e.;ToT_{peak} (ns)");
            grS->SetMarkerStyle(20); grS->SetMarkerColor(kAzure+2);
            grS->Draw("AP"); fCal->Draw("SAME");
        }
        if(nResC>0){
            cSum->cd(5);
            TGraphErrors *grRS=new TGraphErrors(nResC,&vBand[0],&vSigCorr[0],
                                                 &vBandZ[0],&vESigCorr[0]);
            grRS->SetTitle(";n p.e.;#sigma_{#Delta-t} (ns)");
            grRS->SetMarkerStyle(21); grRS->SetMarkerColor(kOrange+7);
            grRS->Draw("AP");
        }
        cSum->cd(6);
        TPaveText *pSum=new TPaveText(0.05,0.05,0.95,0.95,"brNDC");
        pSum->SetBorderSize(1); pSum->SetFillColor(kWhite);
        pSum->SetTextFont(42); pSum->SetTextSize(0.055); pSum->SetTextAlign(12);
        pSum->AddText("ALCOR ToT Analysis Summary");
        pSum->AddText(Form("Events analysed: %d",Nev));
        pSum->AddText(Form("LSB: %.4f ns/tick",lsb_ns));
        pSum->AddText(Form("p.e. peaks found: %d",nFit));
        if(nFit>=2)
            pSum->AddText(Form("ToT gain: %.3f ns/p.e.",
                (vToT_mean.back()-vToT_mean[0])/(nFit-1)));
        if(nResC>0&&!vSigCorr.empty())
            pSum->AddText(Form("Best #sigma_{#Delta-t}: %.3f ns @ %d p.e.",
                *std::min_element(vSigCorr.begin(),vSigCorr.end()),
                (int)(std::min_element(vSigCorr.begin(),vSigCorr.end())-vSigCorr.begin())+1));
        pSum->Draw();
        cSum->SaveAs("plots/alcor_15_summary.png");
    }

    // =========================================================================
    // TERMINAL SUMMARY
    // =========================================================================
    printBox({"Analysis complete!",
              Form("Events analysed : %d", Nev),
              Form("LSB             : %.4f ns/tick", lsb_ns),
              Form("p.e. peaks fit  : %d", nFit),
              "All plots saved to plots/"});

    std::cout<<"\n  ToT calibration table:"<<std::endl;
    std::cout<<"  "<<repeatStr("-",40)<<std::endl;
    printf("  %-8s  %-12s  %-12s\n","p.e.","ToT mean (ns)","ToT sigma (ns)");
    std::cout<<"  "<<repeatStr("-",40)<<std::endl;
    for(int i=0;i<nFit;i++)
        printf("  %-8d  %-12.4f  %-12.4f\n",i+1,vToT_mean[i],vToT_sigma[i]);
    std::cout<<"  "<<repeatStr("-",40)<<std::endl;

    // Cleanup ACLiC files (done via shell script — see run_alcor.sh)
    std::cout<<"\n  [INFO] To clean ACLiC files after ROOT exits,"<<std::endl;
    std::cout<<"  run:  ./run_alcor.sh"<<std::endl;
    std::cout<<"\n  Goodbye! All canvases remain open until you close ROOT.\n"<<std::endl;
}
