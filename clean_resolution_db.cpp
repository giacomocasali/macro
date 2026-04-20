// clean_resolution_db.cpp
// Rimuove duplicati dal DB, mantiene solo l'ultima entry per (vbias, threshold_pe)
// Utile dopo re-run con parametri diversi (es. cambio fit window)

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>

struct Key {
    int vbias;
    long long thr_key;  // threshold_pe * 1e6 (per evitare float compare)
    bool operator<(const Key& o) const {
        if (vbias != o.vbias) return vbias < o.vbias;
        return thr_key < o.thr_key;
    }
};

struct Entry {
    double threshold_pe;
    double sigma_ns;
    double sigma_err_ns;
};

void clean_resolution_db(const char* dbPath = "../file_root/resolution_scan_results.root") {
    TFile* f = TFile::Open(dbPath, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Cannot open: " << dbPath << "\n";
        delete f;
        return;
    }
    
    TTree* t = (TTree*)f->Get("resolution_scan");
    if (!t) {
        std::cerr << "[ERROR] Tree not found\n";
        f->Close(); delete f;
        return;
    }
    
    Int_t vbias;
    Double_t threshold_pe, sigma_ns, sigma_err_ns;
    t->SetBranchAddress("vbias", &vbias);
    t->SetBranchAddress("threshold_pe", &threshold_pe);
    t->SetBranchAddress("sigma_ns", &sigma_ns);
    t->SetBranchAddress("sigma_err_ns", &sigma_err_ns);
    
    // Mantiene ULTIMA entry per ogni (vbias, threshold)
    std::map<Key, Entry> unique;
    Long64_t nTot = t->GetEntries();
    
    for (Long64_t i = 0; i < nTot; ++i) {
        t->GetEntry(i);
        Key k{vbias, (long long)(threshold_pe * 1e6 + 0.5)};
        unique[k] = {threshold_pe, sigma_ns, sigma_err_ns};
    }
    
    std::cout << "[INFO] Total entries: " << nTot << "\n";
    std::cout << "[INFO] Unique entries: " << unique.size() << "\n";
    std::cout << "[INFO] Duplicates removed: " << nTot - unique.size() << "\n";
    
    f->Close(); delete f;
    
    // Riscrivi file pulito
    TFile* fOut = TFile::Open(dbPath, "RECREATE");
    Int_t db_vbias;
    Double_t db_thr, db_sigma, db_err;
    TTree* tOut = new TTree("resolution_scan", "sigma vs LET threshold");
    tOut->Branch("vbias", &db_vbias, "vbias/I");
    tOut->Branch("threshold_pe", &db_thr, "threshold_pe/D");
    tOut->Branch("sigma_ns", &db_sigma, "sigma_ns/D");
    tOut->Branch("sigma_err_ns", &db_err, "sigma_err_ns/D");
    
    for (auto& [k, e] : unique) {
        db_vbias = k.vbias;
        db_thr = e.threshold_pe;
        db_sigma = e.sigma_ns;
        db_err = e.sigma_err_ns;
        tOut->Fill();
        std::cout << "  Vbias=" << k.vbias 
                  << "  thr=" << e.threshold_pe << " pe"
                  << "  sigma=" << e.sigma_ns << " ns\n";
    }
    
    tOut->Write();
    fOut->Close();
    delete fOut;
    
    std::cout << "[OK] Cleaned DB: " << dbPath << "\n";
}
