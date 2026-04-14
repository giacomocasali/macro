/**
 * sipm_tot_analysis_v5_lowram.cpp
 *
 * Compile: .L sipm_tot_analysis_v5_lowram.cpp+
 * Run:     sipm_tot_analysis_v5()
 */

#include "../header/Config.h"
#include "../header/OutputManager.h"
#include "../header/Calibration.h"
#include "../header/CalibIO.h"
#include "../header/TimingCorrection.h"
#include "../header/VbiasAnalysis_v2.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>
#include <limits>
#include <TStyle.h>
#include <TSystem.h>
#include <TApplication.h>

static std::vector<double> findCalibratedCutoffs(int vbias, const std::string& dataDir) {
    std::vector<double> cutoffs;
    std::string prefix="calib_vbias"+std::to_string(vbias)+"_cut", suffix="mhz.root";
    void*dirp=gSystem->OpenDirectory(dataDir.c_str());if(!dirp)return cutoffs;
    const char*entry;
    while((entry=gSystem->GetDirEntry(dirp))!=nullptr){std::string fname(entry);
        if(fname.size()<prefix.size()+suffix.size())continue;
        if(fname.substr(0,prefix.size())!=prefix)continue;
        if(fname.substr(fname.size()-suffix.size())!=suffix)continue;
        std::string mid=fname.substr(prefix.size(),fname.size()-prefix.size()-suffix.size());
        try{double co=std::stod(mid);if(co>=0)cutoffs.push_back(co);}catch(...){}}
    gSystem->FreeDirectory(dirp);std::sort(cutoffs.begin(),cutoffs.end());return cutoffs;
}

void sipm_tot_analysis_v5()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    OutCtx ctx = createOutputDirs();
    const std::string dataDir = DATA_DIR;

    std::cout << "\n+==========================================================+\n"
              << "|    SiPM TOT ANALYSIS v5 (LOW-RAM + FIX)                 |\n"
              << "+==========================================================+\n"
              << "  Data dir: " << dataDir << "\n\n";

    // readLine: legge una riga non vuota; gestisce EOF e errori stdin.
    auto readLine = [](const std::string& prompt) -> std::string {
        std::string line;
        while (true) {
            std::cout << prompt << std::flush;
            if (!std::getline(std::cin, line)) {
                // EOF o errore: esci invece di girare per sempre
                std::cerr << "\n[ERROR] stdin closed or unexpected EOF.\n";
                std::exit(1);
            }
            // Rimuovi whitespace iniziale/finale
            auto b = line.find_first_not_of(" \t\r\n");
            if (b != std::string::npos) { line = line.substr(b); break; }
            // Riga vuota → ripete il prompt
        }
        return line;
    };

    // 1. Vbias
    std::vector<int> vbiasList;
    {std::string line=readLine("Vbias (es. 53 54 55):\n> ");std::stringstream ss(line);int v;while(ss>>v)vbiasList.push_back(v);}
    if(vbiasList.empty()){std::cerr<<"[ERROR] No Vbias.\n";return;}

    // 2. LET
    std::vector<double> fracs_pe;
    {std::string line=readLine("LET thresholds (p.e., es. 0.5 1.0 2.0):\n> ");std::stringstream ss(line);double v;while(ss>>v)if(v>0)fracs_pe.push_back(v);}
    if(fracs_pe.empty())fracs_pe.push_back(1.0);

    // 3. Calibration
    std::map<int,CalibResult> calMap;
    for(int vbias:vbiasList){
        auto cutoffs=findCalibratedCutoffs(vbias,dataDir);
        if(cutoffs.empty()){std::cerr<<"  [WARN] Vbias="<<vbias<<": no cal.\n";continue;}
        double chosen=cutoffs.back();
        if(cutoffs.size()>1){std::cout<<"  Vbias="<<vbias<<" cutoff: ";for(double co:cutoffs)std::cout<<(int)co<<" ";std::cout<<"MHz\n  Which one? ["<<(int)chosen<<"]: "<<std::flush;std::string line;std::getline(std::cin,line);if(!line.empty()){try{double inp=std::stod(line);double best=chosen,bestD=1e9;for(double co:cutoffs)if(std::abs(co-inp)<bestD){bestD=std::abs(co-inp);best=co;}chosen=best;}catch(...){}}}
        CalibResult cal;if(!loadCalibration(cal,vbias,chosen,dataDir))continue;if(cal.m<=0||cal.m>100)continue;
        calMap[vbias]=cal;
        std::cout<<"  Vbias="<<vbias<<" gain="<<cal.m<<" mV/p.e. trig=["<<cal.t_trig_start<<","<<cal.t_trig_end<<"]\n";
    }
    if(calMap.empty()){std::cerr<<"[ERROR] No cal.\n";return;}

    // 4. Fit window — con validazione e recovery da failbit
    // Helper: legge un double dalla riga già letta (riusa readLine per gestire EOF)
    auto readDouble = [&readLine](const std::string& prompt) -> double {
        while (true) {
            std::string line = readLine(prompt);
            try {
                std::size_t pos;
                double val = std::stod(line, &pos);
                // Verifica che tutta la stringa sia stata consumata
                while (pos < line.size() && std::isspace((unsigned char)line[pos])) ++pos;
                if (pos == line.size()) return val;
            } catch (...) {}
            std::cerr << "  [!] Invalid value, try again.\n";
        }
    };

    double fit_lo = readDouble("\nFit window start [ns]: ");
    double fit_hi = readDouble("Fit window end   [ns]: ");
    if (fit_lo >= fit_hi) { std::cerr << "[ERROR] fit_lo must be < fit_hi.\n"; return; }

    // Helper: legge un char y/n con recovery da failbit e EOF
    auto readYN = [&readLine](const std::string& prompt) -> bool {
        while (true) {
            std::string line = readLine(prompt);
            if (!line.empty()) {
                char c = (char)std::tolower((unsigned char)line[0]);
                if (c == 'y') return true;
                if (c == 'n') return false;
            }
            std::cerr << "  [!] Type y or n.\n";
        }
    };

    TWMethod tw_method = askTimeWalkMethod();
    bool do_pe     = readYN("Per-p.e.? [y/n]: ");
    bool use_filter = readYN("LP filter? [y/n]: ");

    // 5. Summary
    std::cout<<"\n+==========================================================+\n|  Vbias: ";for(int v:vbiasList)std::cout<<v<<" ";
    std::cout<<"V\n|  Fit: ["<<fit_lo<<","<<fit_hi<<"] ns\n|  LET: ";for(double f:fracs_pe)std::cout<<f<<" ";
    std::cout<<"p.e.\n|  TW: "<<(tw_method==TWMethod::EMPIRICAL?"Emp":tw_method==TWMethod::AMPLITUDE?"Amp":"None")
             <<" Filter: "<<(use_filter?"yes":"NO")<<"\n+==========================================================+\n\n";

    // 6. Analysis
    std::map<int,std::map<double,std::pair<double,double>>> sigmaResults;
    auto t0=std::chrono::steady_clock::now();
    for(int vi=0;vi<(int)vbiasList.size();++vi){
        int vbias=vbiasList[vi];if(calMap.find(vbias)==calMap.end())continue;
        const CalibResult&cal=calMap[vbias];auto tVb=std::chrono::steady_clock::now();
        auto res=processOneVbias_v2(vbias,fracs_pe,cal.cutoff_MHz,cal.t_trig_start,cal.t_trig_end,fit_lo,fit_hi,tw_method,do_pe,use_filter,ctx,cal);
        if(!res.empty())sigmaResults[vbias]=res;
        std::cout<<"  Vbias="<<vbias<<" done in "<<std::fixed<<std::setprecision(1)<<std::chrono::duration<double>(std::chrono::steady_clock::now()-tVb).count()<<" s\n";
    }
    drawVbiasSummary(sigmaResults,ctx);
    std::cout<<"\n+==========================================================+\n|  DONE  "<<std::fixed<<std::setprecision(1)<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count()<<" s\n|  PNG:  "<<ctx.pngDir<<"\n|  ROOT: "<<ctx.rootDir<<"\n+==========================================================+\n";

    // Reopen the last N canvases from their .root files for interactive editing.
    // All canvases (including earlier ones) are permanently saved as .root files
    // in ctx.canvasDir — you can always open them later with TBrowser or:
    //   TFile f("canvases/XX_name.root"); TCanvas* c = (TCanvas*)f.Get("cName");
    ctx.reopenAllCanvases(40);

    std::cout << "\n  Canvases reopened. Type .q to exit.\n"
              << "  All canvas .root files: " << ctx.canvasDir << "\n";
    if (gApplication) gApplication->Run(kTRUE);
}
