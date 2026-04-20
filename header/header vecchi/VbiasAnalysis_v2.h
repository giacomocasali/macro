#pragma once
// VbiasAnalysis_v2.h — v4
// FIX RAM: processa un file alla volta (no TChain → da 14 GB a ~2 GB)
// FIX laser: 5% del picco
// FIX edge: adattivo a use_filter

#include "Config.h"
#include "OutputManager.h"
#include "SignalProcessing.h"
#include "Calibration.h"
#include "CalibIO.h"
#include "EventCache.h"
#include "TOTAnalysis.h"
#include "TimingCorrection.h"
#include "TOTPlotting.h"
#include "VbiasSummary.h"
#include "FilterDiagnostics.h"
#include "ProgressBar.h"
#include "SidebandAnalysis.h"
#include "WaveformPlotter.h"
#include "PersistencePlot.h"
#include "ChunkedHistoFill.h"
#include "TOTBandPersistence.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <functional>
#include <limits>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TBox.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TROOT.h>

// Laser threshold: 5% del picco
static double estimateLaserThreshold(const std::string& firstFile, double& median_amp_out) {
    median_amp_out=0;
    TFile*f=TFile::Open(firstFile.c_str(),"READ");if(!f||f->IsZombie()){delete f;return 10;}
    TTree*tr=(TTree*)f->Get("laser");if(!tr){f->Close();delete f;return 10;}
    const int N=1024;Double_t tL[N],aL[N];
    tr->SetBranchAddress("time",tL);tr->SetBranchAddress("amplitude",aL);
    std::vector<double>peaks;Long64_t nS=std::min((Long64_t)500,tr->GetEntries());
    for(Long64_t i=0;i<nS;++i){tr->GetEntry(i);
        std::vector<double>pre;for(int j=0;j<N;++j)if(tL[j]<BASELINE_END)pre.push_back(aL[j]);
        double off=0;if(!pre.empty()){auto tmp=pre;std::nth_element(tmp.begin(),tmp.begin()+tmp.size()/2,tmp.end());off=tmp[tmp.size()/2];}
        double pk=-1e9;for(int j=0;j<N;++j)if(aL[j]-off>pk)pk=aL[j]-off;
        if(pk>0)peaks.push_back(pk);}
    f->Close();delete f;
    if(peaks.empty())return 10;
    std::nth_element(peaks.begin(),peaks.begin()+peaks.size()/2,peaks.end());
    median_amp_out=peaks[peaks.size()/2];
    double thr=std::max(10.,std::min(median_amp_out*0.05,50.));
    std::cout<<"  [Laser] Median amp: "<<std::fixed<<std::setprecision(1)<<median_amp_out<<" mV → thr: "<<thr<<" mV (5%)\n";
    return thr;
}


// ═══════════════════════════════════════════════════════════════
//  collectTOTEvents_fileByFile
//  Processa UN FILE ALLA VOLTA — apre, scorre, chiude.
//  RAM: ~30 MB per file (1 TTree ch1 + 1 TTree laser + buffer).
//  Tutti i risultati finiscono nello stesso TTree di output.
// ═══════════════════════════════════════════════════════════════
static void collectTOTEvents_fileByFile(
        const std::map<int,std::string>& runFiles,
        double cutoff_MHz, double fs_MHz,
        int j_trig_start, int j_trig_end,
        double let_thr, const CalibResult& cal,
        const std::string& tag, const std::string& dataDir,
        int vbias, double frac_pe, bool use_filter,
        TH2D* hPersA = nullptr,   // persistenza banda A (TOT 40-43 ns)
        TH2D* hPersB = nullptr)   // persistenza banda B (TOT 50-60 ns)
{
    const int N = 1024;
    const double edge_thr_frac = use_filter ? 0.5 : 100.0;

    // Apri il file di output cache
    std::string cachePath = eventCachePath(vbias, frac_pe, cutoff_MHz,
                                            cal.laser_thr, dataDir, use_filter);
    TFile* fOut = new TFile(cachePath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) { delete fOut; return; }

    Double_t tot_br, delta_t_br, amp_max_br, laser_thr_br = cal.laser_thr;
    Int_t n_pe_br;
    TTree* tOut = new TTree("events", Form("TOT vbias%d let%.2fpe", vbias, frac_pe));
    tOut->Branch("tot", &tot_br, "tot/D");
    tOut->Branch("delta_t", &delta_t_br, "delta_t/D");
    tOut->Branch("amp_max", &amp_max_br, "amp_max/D");
    tOut->Branch("n_pe", &n_pe_br, "n_pe/I");
    tOut->Branch("laser_thr_saved", &laser_thr_br, "laser_thr_saved/D");
    tOut->SetAutoSave(50000);

    long totalEntries = 0, nNoLaser = 0, nNoTOT = 0, nAccepted = 0;
    long nDirtyBaseline = 0;

    // ┌──────────────────────────────────────────────────────┐
    // │  Baseline RMS threshold:                              │
    // │  Senza filtro: il rumore RF ha RMS ~1-2 mV.           │
    // │  Un dark count nei primi 30 ns porta RMS a 5-10+ mV.  │
    // │  Soglia: 4.0 mV senza filtro, 2.0 mV con filtro.     │
    // │  Abbastanza lasco da non rigettare il rumore normale, │
    // │  abbastanza stretto da beccare un segnale in baseline.│
    // └──────────────────────────────────────────────────────┘
    const double baseline_max_rms = use_filter ? 2.0 : 4.0;

    // Conta totale eventi per progress bar
    for (auto& [r, path] : runFiles) {
        TFile* fTmp = TFile::Open(path.c_str(), "READ");
        if (fTmp && !fTmp->IsZombie()) {
            TTree* tr = (TTree*)fTmp->Get("ch1");
            if (tr) totalEntries += tr->GetEntries();
            fTmp->Close();
        }
        delete fTmp;
    }

    ProgressBar bar(totalEntries, "TOT " + tag + (use_filter ? "" : " [NO FILT]"));
    long globalIdx = 0;

    // ── Processa un file alla volta ──────────────────────────
    for (auto& [runNum, path] : runFiles) {
        TFile* fIn = TFile::Open(path.c_str(), "READ");
        if (!fIn || fIn->IsZombie()) { delete fIn; continue; }

        TTree* treeCh1 = (TTree*)fIn->Get("ch1");
        TTree* treeLaser = (TTree*)fIn->Get("laser");
        if (!treeCh1 || !treeLaser) { fIn->Close(); delete fIn; continue; }

        Double_t t1[N], a1[N], tL[N], aL[N];
        treeCh1->SetBranchAddress("time", t1);
        treeCh1->SetBranchAddress("amplitude", a1);
        treeLaser->SetBranchAddress("time", tL);
        treeLaser->SetBranchAddress("amplitude", aL);

        // Limita la cache di lettura per risparmiare RAM
        treeCh1->SetCacheSize(2 * 1024 * 1024);   // 2 MB
        treeLaser->SetCacheSize(2 * 1024 * 1024);

        Long64_t nFile  = treeCh1->GetEntries();
        Long64_t nLaser = treeLaser->GetEntries();
        if (nLaser < nFile) {
            std::cerr << "  [WARN] run " << runNum
                      << ": ch1=" << nFile << " entries, laser=" << nLaser
                      << " - processing only " << nLaser << " aligned events.\n";
            nFile = nLaser;
        }

        for (Long64_t i = 0; i < nFile; ++i) {
            if (bar.update(globalIdx, nNoLaser, nNoTOT, nAccepted)) goto done;
            ++globalIdx;

            treeCh1->GetEntry(i);
            treeLaser->GetEntry(i);

            double t_laser = laserTriggerTime(tL, aL, N, cal.laser_thr);
            if (t_laser < -900.0) { ++nNoLaser; continue; }

            std::vector<double> v_t(t1, t1+N), v_a(a1, a1+N);
            bool bl_ok = true;
            std::vector<double> af = correctBaseline(v_t, v_a,
                BASELINE_START, BASELINE_END, baseline_max_rms, &bl_ok);
            if (!bl_ok) { ++nDirtyBaseline; continue; }
            if (use_filter) af = butterworthLowPass(af, cutoff_MHz, fs_MHz);

            double amp_max = *std::max_element(af.begin()+j_trig_start,
                                                af.begin()+j_trig_end+1);

            auto [t_rise, t_fall] = computeTOT(v_t, af, let_thr,
                j_trig_start, j_trig_end,
                100, edge_thr_frac, 0.5, 10, 0.3, 50, false, 50, 10);

            if (t_rise < 0 || t_fall < 0) { ++nNoTOT; continue; }

            double tot = t_fall - t_rise;
            double delta_t = t_rise - t_laser;
            if (tot <= 0.0 || tot >= 150.0) continue;

            int n_pe = estimatePE(amp_max, cal);

            // Scrivi nel file di output (fOut è ancora la current directory)
            fOut->cd();
            tot_br = tot; delta_t_br = delta_t;
            amp_max_br = amp_max; n_pe_br = n_pe;
            tOut->Fill();
            ++nAccepted;

            // Persistenza per banda TOT
            auto fillPers = [&](TH2D* h, double tot_lo, double tot_hi) {
                if (!h || tot < tot_lo || tot >= tot_hi) return;
                int ns = (int)v_t.size();
                for (int j = 0; j < ns; ++j) h->Fill(v_t[j] - t_laser, af[j]);
            };
            fillPers(hPersA, 40.0, 43.0);
            fillPers(hPersB, 50.0, 60.0);
        }

        // Chiudi il file di input — libera TUTTA la RAM di questo file
        treeCh1->SetCacheSize(0);
        treeLaser->SetCacheSize(0);
        fIn->Close();
        delete fIn;
    }

done:
    bar.done();
    fOut->cd();
    fOut->Write();
    fOut->Close();
    delete fOut;
    std::cout << "  [Cache] " << nAccepted << " events"
              << "  noLaser=" << nNoLaser
              << "  dirtyBL=" << nDirtyBaseline
              << "  noTOT=" << nNoTOT << "\n";
}


// peakFinderFromCache — invariato
static std::tuple<double,double,double> peakFinderFromCache(
        int vbias, double frac_pe, double cutoff_MHz, double laser_thr,
        const std::string& dataDir = DATA_DIR) {
    std::string path=eventCachePath(vbias,frac_pe,cutoff_MHz,laser_thr,dataDir);
    if(gSystem->AccessPathName(path.c_str()))return{-1,-1,-1};
    TH1D*hTOT=new TH1D("hTp","",1500,0,150);hTOT->SetDirectory(nullptr);
    {ChunkedFiller f(path);f.addTH1D(hTOT,[](TH1D*h,const CacheEvent&e){h->Fill(e.tot);});if(f.run()==0){delete hTOT;return{-1,-1,-1};}}
    double total=hTOT->Integral(),cumul=0,tot_cut=15;
    for(int b=1;b<=hTOT->GetNbinsX();++b){cumul+=hTOT->GetBinContent(b);if(cumul>=total*0.15){tot_cut=hTOT->GetBinCenter(b);break;}}
    delete hTOT;tot_cut=std::max(2.,std::min(tot_cut,15.));
    TH1D*hPk=new TH1D("hPF","",1500,-50,200);hPk->SetDirectory(nullptr);
    {ChunkedFiller f(path);f.addTH1D(hPk,[tot_cut](TH1D*h,const CacheEvent&e){if(e.tot>0&&e.tot<tot_cut)h->Fill(e.delta_t);});f.run();}
    double tpk=hPk->GetBinCenter(hPk->GetMaximumBin());delete hPk;
    return{tpk,tpk-5,tpk+5};
}


// analyseOneLET_chunked — invariato (già usa ChunkedFiller, non TChain)
struct OneLETResult{double sigma=0,sigmaErr=0;bool valid=false;std::string h2d_tmpfile;};

static OneLETResult analyseOneLET_chunked(
        const std::string& cachePath, double frac, double let_thr,
        double fit_lo, double fit_hi, TWMethod tw_method, bool do_pe_analysis,
        const std::string& ltag, int vbias, OutCtx& ctx)
{
    OneLETResult result;
    Long64_t nEvents=countCacheEvents(cachePath);
    if(nEvents==0){std::cout<<"  +-- No events.\n";return result;}
    std::cout<<"  [Chunked] "<<nEvents<<" events\n";

    TH2D*h2D=new TH2D(Form("h2D_%s",ltag.c_str()),Form("TOT vs #Deltat LET=%.2f p.e. %s;TOT (ns);#Deltat (ns)",frac,ltag.c_str()),300,0,150,5000,-50,200);h2D->SetDirectory(nullptr);
    TH1D*hTp=new TH1D(Form("hTp_%s",ltag.c_str()),"",1500,0,150);hTp->SetDirectory(nullptr);
    TH1D*hDt=new TH1D(Form("hDt_%s",ltag.c_str()),Form("#Deltat %s;#Deltat (ns);Counts",ltag.c_str()),1250,-50,200);hDt->SetDirectory(nullptr);
    int nBP=std::max(400,(int)std::round((fit_hi-fit_lo)/0.02));
    TH1D*hPr=new TH1D(Form("hPr_%s",ltag.c_str()),Form("Projection %s;#Deltat (ns);Events",ltag.c_str()),nBP,fit_lo-1,fit_hi+1);hPr->SetDirectory(nullptr);
    const double sbw=SB::DEFAULT_SB_WIDTH,sbx=SB::DEFAULT_EXT,gg=SB::SB_GAP;
    double wlo=fit_lo-sbx,whi=fit_hi+sbx;int nBW=std::max(5,(int)std::round((whi-wlo)/SB::BIN_WIDTH));
    TH1D*hW=new TH1D(Form("hW_%s",ltag.c_str()),"",nBW,wlo,whi);hW->SetDirectory(nullptr);hW->Sumw2();
    TH1D*hSL=new TH1D(Form("hSL_%s",ltag.c_str()),"",std::max(5,(int)std::round(sbw/SB::BIN_WIDTH)),fit_lo-gg-sbw,fit_lo-gg);hSL->SetDirectory(nullptr);hSL->Sumw2();
    TH1D*hSR=new TH1D(Form("hSR_%s",ltag.c_str()),"",std::max(5,(int)std::round(sbw/SB::BIN_WIDTH)),fit_hi+gg,fit_hi+gg+sbw);hSR->SetDirectory(nullptr);hSR->Sumw2();
    TH1D*hSg=new TH1D(Form("hSg_%s",ltag.c_str()),"",std::max(5,(int)std::round((fit_hi-fit_lo)/SB::BIN_WIDTH)),fit_lo,fit_hi);hSg->SetDirectory(nullptr);hSg->Sumw2();
    StatAccum dtSt;
    {ChunkedFiller p(cachePath);
     p.addTH2D(h2D,[](TH2D*h,const CacheEvent&e){h->Fill(e.tot,e.delta_t);});
     p.addTH1D(hTp,[](TH1D*h,const CacheEvent&e){h->Fill(e.tot);});
     p.addTH1D(hDt,[](TH1D*h,const CacheEvent&e){h->Fill(e.delta_t);});
     p.addTH1D(hPr,[](TH1D*h,const CacheEvent&e){h->Fill(e.delta_t);});
     p.addTH1D(hW,[wlo,whi](TH1D*h,const CacheEvent&e){if(e.delta_t>=wlo&&e.delta_t<whi)h->Fill(e.delta_t);});
     p.addTH1D(hSL,[fit_lo,gg,sbw](TH1D*h,const CacheEvent&e){if(e.delta_t>=fit_lo-gg-sbw&&e.delta_t<fit_lo-gg)h->Fill(e.delta_t);});
     p.addTH1D(hSR,[fit_hi,gg,sbw](TH1D*h,const CacheEvent&e){if(e.delta_t>=fit_hi+gg&&e.delta_t<fit_hi+gg+sbw)h->Fill(e.delta_t);});
     p.addTH1D(hSg,[fit_lo,fit_hi](TH1D*h,const CacheEvent&e){if(e.delta_t>=fit_lo&&e.delta_t<fit_hi)h->Fill(e.delta_t);});
     p.addStatAccum(&dtSt,[](const CacheEvent&e){return e.delta_t;});
     p.run();}
    double tot_cut=15;{double total=hTp->Integral(),cumul=0;for(int b=1;b<=hTp->GetNbinsX();++b){cumul+=hTp->GetBinContent(b);if(cumul>=total*0.15){tot_cut=hTp->GetBinCenter(b);break;}}}delete hTp;tot_cut=std::max(2.,std::min(tot_cut,15.));
    TH1D*hPF=new TH1D(Form("hPF_%s",ltag.c_str()),Form("Peak finder (TOT<%.1f) %s",tot_cut,ltag.c_str()),1500,-50,200);hPF->SetDirectory(nullptr);
    {ChunkedFiller pf(cachePath);pf.addTH1D(hPF,[tot_cut](TH1D*h,const CacheEvent&e){if(e.tot>0&&e.tot<tot_cut)h->Fill(e.delta_t);});pf.run();}

    // ── Plots ──
    /*
    {double tpk=hPF->GetBinCenter(hPF->GetMaximumBin());std::cout<<"  [Peak] TOT<"<<tot_cut<<" ns Peak="<<std::fixed<<std::setprecision(1)<<tpk<<" ns\n";
     TCanvas*c=new TCanvas(Form("cPF_%s",ltag.c_str()),"",900,500);c->SetGrid();c->SetLeftMargin(PAD_LEFT);c->SetRightMargin(PAD_RIGHT);c->SetBottomMargin(PAD_BOTTOM);c->SetTopMargin(PAD_TOP);
     hPF->SetLineColor(kAzure+1);hPF->SetLineWidth(2);hPF->GetXaxis()->SetRangeUser(tpk-30,tpk+30);hPF->Draw("HIST");c->Update();c->Modified();ctx.savePNG(c,Form("peak_finder_%s.png",ltag.c_str()));}
    */
    {TCanvas*cFull=new TCanvas(Form("c2Dfull_%s",ltag.c_str()),"",900,700);cFull->SetRightMargin(0.15);cFull->SetLeftMargin(PAD_LEFT);cFull->SetBottomMargin(PAD_BOTTOM);cFull->SetTopMargin(PAD_TOP);cFull->SetLogz();cFull->SetGrid();h2D->Draw("COLZ");cFull->Update();cFull->Modified();ctx.savePNG(cFull,Form("tot_map_%s.png",ltag.c_str()));}
    {int bLo=hPr->FindBin(fit_lo),bHi=hPr->FindBin(fit_hi),bMax=bLo;for(int b=bLo;b<=bHi;++b)if(hPr->GetBinContent(b)>hPr->GetBinContent(bMax))bMax=b;
     double mu=hPr->GetBinCenter(bMax),sigma=0.5;{double sw=0,swx=0,swx2=0;int bl=hPr->FindBin(mu-5),br=hPr->FindBin(mu+5);for(int b=bl;b<=br;++b){double c2=hPr->GetBinContent(b),x=hPr->GetBinCenter(b);sw+=c2;swx+=c2*x;swx2+=c2*x*x;}if(sw>0){mu=swx/sw;sigma=std::max(std::sqrt(std::max(swx2/sw-mu*mu,0.)),0.05);}}
     int yMin=(int)std::lround(mu-1.0*sigma),yMax=(int)std::lround(mu+1.0*sigma);if(yMax<=yMin){yMin=(int)std::floor(mu)-1;yMax=(int)std::ceil(mu)+1;}
     TH2D*h2Zoom=(TH2D*)h2D->Clone(Form("h2Zoom_%s",ltag.c_str()));h2Zoom->SetDirectory(nullptr);
     // Keep native Y resolution for narrow zoom windows (no rebin).
     TCanvas*c=new TCanvas(Form("c2Dzoom_%s",ltag.c_str()),"",900,700);c->SetRightMargin(0.15);c->SetLeftMargin(PAD_LEFT);c->SetBottomMargin(PAD_BOTTOM);c->SetTopMargin(PAD_TOP);c->SetLogz();c->SetGrid();
     h2Zoom->SetTitle(Form("TOT map zoomed around peak %s;TOT (ns);#Deltat (ns)",ltag.c_str()));h2Zoom->GetYaxis()->SetRangeUser(yMin,yMax);h2Zoom->Draw("COLZ");
     TLine*lMu=new TLine(0.0,mu,150.0,mu);lMu->SetLineColor(kRed+1);lMu->SetLineStyle(2);lMu->SetLineWidth(2);lMu->Draw("same");
     c->Update();c->Modified();ctx.savePNG(c,Form("tot_map_zoom_%s.png",ltag.c_str()));delete lMu;delete h2Zoom;}
    {TCanvas*c=new TCanvas(Form("cTOTx_%s",ltag.c_str()),"",900,600);c->SetGrid();c->SetLeftMargin(PAD_LEFT);c->SetRightMargin(PAD_RIGHT);c->SetBottomMargin(PAD_BOTTOM);c->SetTopMargin(PAD_TOP);c->SetLogy();
     TH1D*hTx=h2D->ProjectionX(Form("hTOTx_%s",ltag.c_str()));hTx->SetDirectory(nullptr);
     hTx->SetTitle(Form("TOT distribution %s;TOT (ns);Counts",ltag.c_str()));
     hTx->SetLineColor(kAzure+1);hTx->SetLineWidth(2);hTx->Draw("HIST");
     TPaveText*ptx=new TPaveText(0.55,0.72,0.93,0.88,"NDC");ptx->SetBorderSize(1);ptx->SetFillColor(0);ptx->SetTextFont(42);ptx->SetTextSize(0.036);
     ptx->AddText(Form("Entries: %.0f",hTx->GetEntries()));ptx->AddText(Form("Mean: %.2f ns",hTx->GetMean()));ptx->AddText(Form("RMS:  %.2f ns",hTx->GetRMS()));ptx->Draw();
     c->Update();c->Modified();ctx.savePNG(c,Form("tot_xprojection_%s.png",ltag.c_str()));delete hTx;}
    /*
    {double mn=dtSt.vmin,mx=dtSt.vmax,mg=std::max(0.05*(mx-mn),1.);TCanvas*c=new TCanvas(Form("cDt_%s",ltag.c_str()),"",1100,600);c->SetGrid();c->SetLeftMargin(PAD_LEFT);c->SetRightMargin(PAD_RIGHT);c->SetBottomMargin(PAD_BOTTOM);c->SetTopMargin(PAD_TOP);double yM=hDt->GetMaximum()*1.15;hDt->SetMaximum(yM);hDt->SetLineColor(kAzure+1);hDt->SetLineWidth(2);hDt->GetXaxis()->SetRangeUser(mn-mg,mx+mg);hDt->Draw("HIST");
     TBox*bF=new TBox(fit_lo,0,fit_hi,yM);bF->SetFillColorAlpha(kGreen-9,0.30);bF->SetLineWidth(0);bF->Draw("same");hDt->Draw("HIST same");
     long nNeg=0;int bZ=hDt->FindBin(0);for(int b=1;b<bZ;++b)nNeg+=(long)hDt->GetBinContent(b);
     TPaveText*pt=new TPaveText(0.14,0.72,0.50,0.88,"NDC");pt->SetBorderSize(1);pt->SetFillColor(0);pt->SetTextFont(42);pt->SetTextSize(0.034);pt->AddText(Form("N=%lld",nEvents));pt->AddText(Form("N(#Deltat<0)=%ld (%.1f%%)",nNeg,100.*nNeg/nEvents));pt->AddText(Form("[%.1f, %.1f] ns",mn,mx));pt->Draw();
     c->Update();c->Modified();ctx.savePNG(c,Form("delta_t_full_%s.png",ltag.c_str()));}
    */
    {TCanvas*c=new TCanvas(Form("cPr_%s",ltag.c_str()),"",900,550);c->SetGrid();c->SetLeftMargin(PAD_LEFT);c->SetRightMargin(PAD_RIGHT);c->SetBottomMargin(PAD_BOTTOM);c->SetTopMargin(PAD_TOP);
     int bLo=hPr->FindBin(fit_lo),bHi=hPr->FindBin(fit_hi),bMax=bLo;for(int b=bLo;b<=bHi;++b)if(hPr->GetBinContent(b)>hPr->GetBinContent(bMax))bMax=b;
     double ctr=hPr->GetBinCenter(bMax),sg=0.5;{double sw=0,swx=0,swx2=0;int bl=hPr->FindBin(ctr-5),br=hPr->FindBin(ctr+5);for(int b=bl;b<=br;++b){double c2=hPr->GetBinContent(b),x=hPr->GetBinCenter(b);sw+=c2;swx+=c2*x;swx2+=c2*x*x;}if(sw>0){double m=swx/sw;sg=std::max(std::sqrt(std::max(swx2/sw-m*m,0.)),0.05);ctr=m;}}
     hPr->SetLineColor(kAzure+1);hPr->SetLineWidth(2);hPr->Draw("HIST");
     TF1*fQ=new TF1(Form("fQ_%s",ltag.c_str()),qGaussProj,std::max(fit_lo,ctr-4*sg),std::min(fit_hi,ctr+4*sg),5);
     fQ->SetParNames("A","#mu","#sigma","q_{1}","q_{2}");
     fQ->SetParameters(hPr->GetMaximum(),ctr,sg,1.0,1.2);
     fQ->FixParameter(3,1.0);
     fQ->SetParLimits(4,1.0,2.0);
     hPr->Fit(fQ,"RQ");
     fQ->SetLineColor(kRed+1);fQ->Draw("same");
     TPaveText*ptG=new TPaveText(0.52,0.68,0.93,0.88,"NDC");ptG->SetBorderSize(1);ptG->SetFillColor(0);ptG->SetTextFont(42);ptG->SetTextSize(0.036);
     ptG->AddText("q-Gaussian fit");
     ptG->AddText(Form("#sigma=%.3f#pm%.3f ns",std::abs(fQ->GetParameter(2)),fQ->GetParError(2)));
     ptG->AddText(Form("#mu=%.3f ns",fQ->GetParameter(1)));
     ptG->AddText(Form("q_{2}=%.3f",fQ->GetParameter(4)));
     ptG->Draw();
     c->Update();c->Modified();ctx.savePNG(c,Form("tot_projection_%s.png",ltag.c_str()));}

    // ── Time walk ──
    TH2D*h2Df=h2D;
    if(tw_method==TWMethod::EMPIRICAL){double tm=0;TF1*fTW=computeTimeWalkFromCache(cachePath,fit_lo,fit_hi,ltag,tm);
     if(fTW){double ref=fTW->Eval(tm);TH2D*h2Dc=new TH2D(Form("h2Dc_%s",ltag.c_str()),"",300,0,150,5000,-50,200);h2Dc->SetDirectory(nullptr);
      TH1D*hWc=new TH1D(Form("hWc_%s",ltag.c_str()),"",nBW,wlo,whi);hWc->SetDirectory(nullptr);hWc->Sumw2();
      TH1D*hSLc=new TH1D(Form("hSLc_%s",ltag.c_str()),"",hSL->GetNbinsX(),fit_lo-gg-sbw,fit_lo-gg);hSLc->SetDirectory(nullptr);hSLc->Sumw2();
      TH1D*hSRc=new TH1D(Form("hSRc_%s",ltag.c_str()),"",hSR->GetNbinsX(),fit_hi+gg,fit_hi+gg+sbw);hSRc->SetDirectory(nullptr);hSRc->Sumw2();
      TH1D*hSgc=new TH1D(Form("hSgc_%s",ltag.c_str()),"",hSg->GetNbinsX(),fit_lo,fit_hi);hSgc->SetDirectory(nullptr);hSgc->Sumw2();
      {auto tw=[fTW,ref](CacheEvent&e){e.delta_t-=(fTW->Eval(e.tot)-ref);};ChunkedFiller p2(cachePath,tw);
       p2.addTH2D(h2Dc,[](TH2D*h,const CacheEvent&e){h->Fill(e.tot,e.delta_t);});
       p2.addTH1D(hWc,[wlo,whi](TH1D*h,const CacheEvent&e){if(e.delta_t>=wlo&&e.delta_t<whi)h->Fill(e.delta_t);});
       p2.addTH1D(hSLc,[fit_lo,gg,sbw](TH1D*h,const CacheEvent&e){if(e.delta_t>=fit_lo-gg-sbw&&e.delta_t<fit_lo-gg)h->Fill(e.delta_t);});
       p2.addTH1D(hSRc,[fit_hi,gg,sbw](TH1D*h,const CacheEvent&e){if(e.delta_t>=fit_hi+gg&&e.delta_t<fit_hi+gg+sbw)h->Fill(e.delta_t);});
       p2.addTH1D(hSgc,[fit_lo,fit_hi](TH1D*h,const CacheEvent&e){if(e.delta_t>=fit_lo&&e.delta_t<fit_hi)h->Fill(e.delta_t);});
       p2.run();}
      ctx.saveRoot(h2D,h2Dc,ltag,frac);delete hW;hW=hWc;delete hSL;hSL=hSLc;delete hSR;hSR=hSRc;delete hSg;hSg=hSgc;h2Df=h2Dc;delete fTW;}
    }else if(tw_method==TWMethod::AMPLITUDE){auto mp=computeTimeWalkAmplitudeFromCache(cachePath,fit_lo,fit_hi);
     if(!mp.empty()){TH2D*h2Dc=new TH2D(Form("h2Dc_%s",ltag.c_str()),"",300,0,150,5000,-50,200);h2Dc->SetDirectory(nullptr);
      {auto tw=[&mp](CacheEvent&e){if(e.n_pe<0)return;auto it=mp.find(e.n_pe);if(it!=mp.end())e.delta_t-=it->second;};ChunkedFiller p2(cachePath,tw);p2.addTH2D(h2Dc,[](TH2D*h,const CacheEvent&e){h->Fill(e.tot,e.delta_t);});p2.run();}
      ctx.saveRoot(h2D,h2Dc,ltag,frac);h2Df=h2Dc;}
    }else{ctx.saveRoot(h2D,nullptr,ltag,frac);}

    if(do_pe_analysis){auto ev=loadEventsSmall(cachePath,200000);if(!ev.empty()){if(tw_method==TWMethod::EMPIRICAL){double tm=0;TF1*f=computeTimeWalkFromCache(cachePath,fit_lo,fit_hi,ltag+"_pe",tm);if(f){double r=f->Eval(tm);for(auto&e:ev)e.delta_t-=(f->Eval(e.tot)-r);delete f;}}drawByPE(ev,ltag,fit_lo,fit_hi,ctx);}}

    // ── Sideband ──
    {double nL=hSL->Integral(),nR=hSR->Integral();int bL=hSL->GetNbinsX(),bR=hSR->GetNbinsX();double bpb=0,berr=0;
     if(nL>0&&nR>0){bpb=0.5*(nL/bL+nR/bR);berr=0.5*std::sqrt(nL/(bL*bL)+nR/(bR*bR));}else if(nL>0){bpb=nL/bL;berr=std::sqrt(nL)/bL;}else if(nR>0){bpb=nR/bR;berr=std::sqrt(nR)/bR;}
     std::cout<<"  [SB] Bkg="<<std::fixed<<std::setprecision(3)<<bpb<<" +/- "<<berr<<"\n";
     TH1D*hSub=(TH1D*)hSg->Clone(Form("hSub_%s",ltag.c_str()));hSub->SetDirectory(nullptr);
     for(int b=1;b<=hSub->GetNbinsX();++b){double raw=hSub->GetBinContent(b);hSub->SetBinContent(b,raw-bpb);hSub->SetBinError(b,std::sqrt(raw+bpb+berr*berr));}
     double sigma=0,sigmaErr=0;std::string src="none";TF1*fSB=nullptr;bool fitOk=false;
     if(hW->Integral()>=SB::MIN_EVENTS_FIT){
      int bLo=hW->FindBin(fit_lo+1e-6),bHi=hW->FindBin(fit_hi-1e-6),bPk=bLo;double hMax=0;
      for(int b=bLo;b<=bHi;++b)if(hW->GetBinContent(b)>hMax){hMax=hW->GetBinContent(b);bPk=b;}
      double mu0=hW->GetBinCenter(bPk),wH=(fit_hi-fit_lo)/4,sw=0,swx=0,swx2=0;
      for(int b=hW->FindBin(mu0-wH);b<=hW->FindBin(mu0+wH);++b){double cc=hW->GetBinContent(b),x=hW->GetBinCenter(b);sw+=cc;swx+=cc*x;swx2+=cc*x*x;}
      double s0=2;if(sw>0){double m=swx/sw,v=swx2/sw-m*m;if(v>0)s0=std::max(0.3,std::min(std::sqrt(v),(fit_hi-fit_lo)/2));mu0=m;}
      double C0=0;{double sm=0;int cnt=0;for(int b=1;b<bLo;++b){sm+=hW->GetBinContent(b);++cnt;}for(int b=bHi+1;b<=hW->GetNbinsX();++b){sm+=hW->GetBinContent(b);++cnt;}if(cnt>0)C0=std::max(0.,sm/cnt);}
      fSB=new TF1(Form("fSB_%s",ltag.c_str()),"gaus(0)+pol0(3)",wlo,whi);fSB->SetParameters(std::max(1.,hMax-C0),mu0,s0,C0);fSB->SetParLimits(0,0,1e8);fSB->SetParLimits(2,0.1,fit_hi-fit_lo);if(bpb>0)fSB->SetParLimits(3,0,std::max(bpb*2,1.));
      int st=hW->Fit(fSB,"RQLN");fitOk=(st==0||st==4000);if(fitOk){sigma=std::abs(fSB->GetParameter(2));sigmaErr=fSB->GetParError(2);src="S+B";if(sigma<0.1||sigma>(fit_hi-fit_lo))fitOk=false;}
      if(!fitOk&&hSub->Integral()>10){TF1*fG=new TF1("fg","gaus",std::max(fit_lo,mu0-2*s0),std::min(fit_hi,mu0+2*s0));fG->SetParameters(hSub->GetMaximum(),mu0,s0);fG->SetParLimits(2,0.1,fit_hi-fit_lo);if(hSub->Fit(fG,"RQWN")==0){sigma=std::abs(fG->GetParameter(2));sigmaErr=fG->GetParError(2);src="sub";fitOk=true;}delete fG;}
      if(!fitOk){double ct=hSg->GetBinCenter(hSg->GetMaximumBin()),sg2=std::max(hSg->GetRMS(),0.1);TF1*fG=new TF1("fg2","gaus",std::max(fit_lo,ct-2*sg2),std::min(fit_hi,ct+2*sg2));fG->SetParameters(hSg->GetMaximum(),ct,sg2);hSg->Fit(fG,"RQN");sigma=std::abs(fG->GetParameter(2));sigmaErr=fG->GetParError(2);src="legacy";delete fG;}
      std::cout<<"  [SB] sigma="<<std::fixed<<std::setprecision(3)<<sigma<<" +/- "<<sigmaErr<<" ns ("<<src<<")\n";}
     /*
     {TCanvas*c=new TCanvas(Form("cSB_%s",ltag.c_str()),"",1000,500);c->Divide(2,1);
      c->cd(1);gPad->SetGrid();gPad->SetLeftMargin(PAD_LEFT);gPad->SetRightMargin(PAD_RIGHT);gPad->SetBottomMargin(PAD_BOTTOM);gPad->SetTopMargin(PAD_TOP);hW->SetLineColor(kBlack);hW->SetLineWidth(2);hW->Draw("HIST E");if(fSB&&fitOk){fSB->SetLineColor(kRed);fSB->SetLineWidth(2);fSB->Draw("same");}
      TPaveText*pt=new TPaveText(0.50,0.62,0.93,0.88,"NDC");pt->SetBorderSize(1);pt->SetFillColor(0);pt->SetTextFont(42);pt->SetTextSize(0.036);pt->AddText(Form("#sigma=%.3f#pm%.3f ns",sigma,sigmaErr));pt->AddText(src.c_str());pt->AddText(Form("Bkg=%.2f",bpb));pt->Draw();
      c->cd(2);gPad->SetGrid();gPad->SetLeftMargin(PAD_LEFT);gPad->SetRightMargin(PAD_RIGHT);gPad->SetBottomMargin(PAD_BOTTOM);gPad->SetTopMargin(PAD_TOP);hSub->SetLineColor(kAzure+1);hSub->SetLineWidth(2);hSub->SetMarkerStyle(20);hSub->SetMarkerSize(0.5);hSub->Draw("E");
      c->Update();c->Modified();ctx.savePNG(c,Form("sideband_%s.png",ltag.c_str()));}
     */
     delete fSB;
     if(sigma>0&&sigma<(fit_hi-fit_lo))result={sigma,sigmaErr,true,""};
     delete hSub;}
    // savePNG() now writes both PNG and .root, then deletes the canvas.
    // Histograms are embedded in the saved .root canvas — safe to free here.
    delete hW; delete hSL; delete hSR; delete hSg;
    delete hPF; delete hDt; delete hPr;
    {std::string tmp=Form("/tmp/h2d_%s.root",ltag.c_str());TFile*fT=new TFile(tmp.c_str(),"RECREATE");if(h2Df)h2Df->Write("h2D");fT->Close();delete fT;result.h2d_tmpfile=tmp;}
    if(h2Df!=h2D)delete h2D;delete h2Df;
    return result;
}


// ═══════════════════════════════════════════════════════════════
//  processOneVbias_v2 — usa collectTOTEvents_fileByFile
// ═══════════════════════════════════════════════════════════════
static std::map<double,std::pair<double,double>> processOneVbias_v2(
        int vbias, const std::vector<double>& fracs_pe,
        double cutoff_MHz, double t_trig_start, double t_trig_end,
        double fit_lo, double fit_hi,
        TWMethod tw_method, bool do_pe_analysis, bool use_filter,
        OutCtx& ctx, const CalibResult& cal_in)
{
    std::map<double,std::pair<double,double>> sigmaOut;
    const std::string dataDir=DATA_DIR;
    const std::string pattern="data.vbias_"+std::to_string(vbias)+"_run_";
    std::map<int,std::string> foundRuns;
    {void*dirp=gSystem->OpenDirectory(dataDir.c_str());if(!dirp)return sigmaOut;const char*entry;
     while((entry=gSystem->GetDirEntry(dirp))!=nullptr){std::string fname(entry);if(fname.find(pattern)==std::string::npos)continue;if(fname.size()<5||fname.substr(fname.size()-5)!=".root")continue;size_t pos=fname.find("_run_");if(pos==std::string::npos)continue;try{std::string sub=fname.substr(pos+5);size_t dot=sub.find(".root");if(dot!=std::string::npos)sub=sub.substr(0,dot);foundRuns[std::stoi(sub)]=dataDir+"/"+fname;}catch(...){}}
     gSystem->FreeDirectory(dirp);}
    if(foundRuns.empty())return sigmaOut;
    std::cout<<"\n+-- Vbias="<<vbias<<" V  "<<foundRuns.size()<<" run(s)\n";

    double fs_MHz = 0;
    {
        TFile* f0 = TFile::Open(foundRuns.begin()->second.c_str(), "READ");
        if (!f0 || f0->IsZombie()) { delete f0; return sigmaOut; }
        TTree* tr = (TTree*)f0->Get("ch1");
        if (tr && tr->GetEntries() >= 2 && tr->GetBranch("time")) {
            const int N0 = 1024; Double_t tb[N0] = {};
            tr->SetBranchAddress("time", tb);
            if (tr->GetEntry(0) > 0) {
                double dt = tb[1] - tb[0];
                if (std::isfinite(dt) && dt > 0.0)
                    fs_MHz = 1000.0 / dt;
                else
                    std::cerr << "[ERROR] Vbias=" << vbias
                              << ": invalid time axis (dt=" << dt << " ns).\n";
            }
        }
        f0->Close(); delete f0;
    }
    if (fs_MHz <= 0 || !std::isfinite(fs_MHz)) {
        std::cerr << "[ERROR] Vbias=" << vbias
                  << ": unable to estimate fs_MHz - file skipped.\n";
        return sigmaOut;
    }

    CalibResult cal=cal_in;
    const std::string calTag="vbias"+std::to_string(vbias);

    if(!cal.ok)return sigmaOut;
    std::cout<<"  Gain="<<cal.m<<" mV/p.e.  Offset="<<cal.q<<" mV  laser_thr="<<cal.laser_thr<<" mV\n";

    int j_start=0,j_end=1023;
    {TFile*f0=TFile::Open(foundRuns.begin()->second.c_str(),"READ");if(f0&&!f0->IsZombie()){TTree*tr=(TTree*)f0->Get("ch1");if(tr){const int N0=1024;Double_t tb[N0];tr->SetBranchAddress("time",tb);tr->GetEntry(0);triggerWindowIndices(tb,N0,t_trig_start,t_trig_end,j_start,j_end,calTag);}f0->Close();delete f0;}}
    std::cout<<"  Trigger: ["<<t_trig_start<<","<<t_trig_end<<"] → ["<<j_start<<","<<j_end<<"]\n";

    std::vector<std::string> tmps;
    for(double frac:fracs_pe){
        if(gROOT->IsInterrupted())break;
        double let_thr=cal.q+frac*cal.m;
        std::cout<<"\n  +-- LET="<<frac<<" p.e.  thr="<<std::fixed<<std::setprecision(2)<<let_thr<<" mV\n";
        const std::string ltag=Form("vbias%d_let%.2fpe",vbias,frac);
        std::string cp=eventCachePath(vbias,frac,cutoff_MHz,cal.laser_thr,dataDir,use_filter);
        bool cached=!gSystem->AccessPathName(cp.c_str())&&countCacheEvents(cp)>0;

        /*
        // Persistenze per banda TOT (DISABLED):
        // commentate per usare la cache eventi senza rieseguire il loop raw.
        TH2D* hPersA = new TH2D(Form("hPersA_%s",ltag.c_str()),
            Form("Persistence TOT #in [40,43) ns  %s;t#minust_{laser} (ns);Amplitude (mV)",ltag.c_str()),
            1024,-50,200, 440,-20,200);
        hPersA->SetDirectory(nullptr);
        TH2D* hPersB = new TH2D(Form("hPersB_%s",ltag.c_str()),
            Form("Persistence TOT #in [50,60) ns  %s;t#minust_{laser} (ns);Amplitude (mV)",ltag.c_str()),
            1024,-50,200, 440,-20,200);
        hPersB->SetDirectory(nullptr);
        */

        if(!cached){
            // ┌─────────────────────────────────────────────┐
            // │  FILE-BY-FILE: apre un file, processa,      │
            // │  chiude, poi il prossimo. RAM: ~30 MB.      │
            // │  (vs TChain: ~14 GB con 10 file)            │
            // └─────────────────────────────────────────────┘
            collectTOTEvents_fileByFile(foundRuns, cutoff_MHz, fs_MHz,
                j_start, j_end, let_thr, cal, ltag, dataDir,
                vbias, frac, use_filter, nullptr, nullptr);
        } else {
            // Cache già presente: non rieseguire loop raw.
            std::cout << "  [Cache] hit: " << cp << "\n";
        }

        auto res=analyseOneLET_chunked(cp,frac,let_thr,fit_lo,fit_hi,tw_method,do_pe_analysis,ltag,vbias,ctx);
        if(res.valid)sigmaOut[frac]={res.sigma,res.sigmaErr};
        tmps.push_back(res.h2d_tmpfile);

        /*
        // ── Disegna persistenze (DISABLED) ────────────────────
        auto drawPers = [&](TH2D* h, const char* bandLabel, const char* fname) {
            if (!h || h->GetEntries() == 0) {
                std::cout << "  [Pers] No events in " << bandLabel << "\n";
                delete h; return;
            }
            std::cout << "  [Pers] " << bandLabel << ": " << (long)h->GetEntries() << " samples\n";
            TCanvas* c = new TCanvas(Form("cPers%s_%s", bandLabel, ltag.c_str()),"",1100,600);
            c->SetLeftMargin(PAD_LEFT); c->SetRightMargin(0.13f);
            c->SetBottomMargin(PAD_BOTTOM); c->SetTopMargin(PAD_TOP);
            c->SetGrid(); c->SetLogz(1);
            gStyle->SetPalette(kBird); h->SetContour(100);
            h->GetXaxis()->SetTitleSize(0.048f); h->GetYaxis()->SetTitleSize(0.048f);
            h->GetXaxis()->SetLabelSize(0.040f); h->GetYaxis()->SetLabelSize(0.040f);
            h->GetZaxis()->SetTitle("Counts/bin"); h->GetZaxis()->SetTitleSize(0.040f);
            h->Draw("COLZ");
            TLine* lL=new TLine(0,-20,0,200); lL->SetLineColor(kGreen+2); lL->SetLineStyle(2); lL->SetLineWidth(2); lL->Draw("same");
            TLine* lT=new TLine(-50,let_thr,200,let_thr); lT->SetLineColor(kAzure+1); lT->SetLineStyle(2); lT->SetLineWidth(2); lT->Draw("same");
            c->Update(); c->Modified();
            ctx.savePNG(c, Form("%s_%s.png", fname, ltag.c_str()));
            delete lL; delete lT; delete h;
        };
        drawPers(hPersA, "A_tot40_43", "persistence_tot40_43");
        drawPers(hPersB, "B_tot50_60", "persistence_tot50_60");
        */
    }

    if(fracs_pe.size()>1){std::vector<TH2D*>maps;for(auto&tmp:tmps){if(tmp.empty()){maps.push_back(nullptr);continue;}TFile*f=TFile::Open(tmp.c_str(),"READ");TH2D*h=nullptr;if(f&&!f->IsZombie()){h=(TH2D*)f->Get("h2D");if(h)h->SetDirectory(nullptr);f->Close();}delete f;maps.push_back(h);gSystem->Unlink(tmp.c_str());}
     drawMultiLETOverlay(fracs_pe,maps,calTag,fit_lo,fit_hi,ctx);for(auto*h:maps)delete h;}
    else{for(auto&tmp:tmps)if(!tmp.empty())gSystem->Unlink(tmp.c_str());}
    return sigmaOut;
}
