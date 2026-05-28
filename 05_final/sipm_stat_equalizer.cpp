/**
 * sipm_stat_equalizer.cpp
 * =======================
 * Legge le cache degli eventi prodotte da sipm_pos_scan() e calcola,
 * per ogni posizione della griglia, quante waveform (run da 100k) servono
 * per raggiungere lo stesso numero di eventi ACCETTATI del punto di
 * riferimento (90, 0) con N_target specificato dall'utente.
 *
 * FISICA:
 *   P_det_accepted(x,y) = n_accepted / n_laser_found  (dalla cache TParameter)
 *   N_waveform_needed(x,y) = ceil( N_target / P_det_accepted(x,y) )
 *   N_run(x,y) = ceil( N_waveform_needed(x,y) / WAVEFORMS_PER_RUN )
 *
 * OUTPUT:
 *   - Tabella a schermo con posizione, P_det%, N_waveform, N_run
 *   - Totale run e tempo stimato
 *
 * Compile:  .L sipm_stat_equalizer.cpp+
 * Run:      sipm_stat_equalizer()
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <TFile.h>
#include <TParameter.h>
#include <TSystem.h>

#include "../header/Config.h"
#include "../header/InputHelpers.h"

// ── Costanti ──────────────────────────────────────────────────────────────────
static constexpr long   WAVEFORMS_PER_RUN = 100000;
static constexpr double MINUTES_PER_RUN   = 7.5;    // minuti per run da 100k

// ── Struct posizione ──────────────────────────────────────────────────────────
struct PosStat {
    double x, y;
    int    vbias;
    double frac_pe;
    long   n_accepted;
    long   n_laser_found;
    long   n_crossing;
    double p_det_acc;    // n_accepted / n_laser_found
    std::string cache_path;
};

// ════════════════════════════════════════════════════════════════════════════
//  readCacheStats: legge TParameter dalla cache ROOT
// ════════════════════════════════════════════════════════════════════════════
static bool readCacheStats(const std::string& path, PosStat& ps)
{
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; return false; }

    auto* complete = dynamic_cast<TParameter<int>*>(f->Get("cache_complete"));
    if (!complete || complete->GetVal() != 1) {
        f->Close(); delete f; return false;
    }

    auto* pacc  = dynamic_cast<TParameter<long>*>(f->Get("n_accepted"));
    auto* plas  = dynamic_cast<TParameter<long>*>(f->Get("n_laser_found"));
    auto* pcros = dynamic_cast<TParameter<long>*>(f->Get("n_crossing"));

    if (!pacc || !plas) {
        f->Close(); delete f; return false;
    }

    ps.n_accepted    = pacc->GetVal();
    ps.n_laser_found = plas->GetVal();
    ps.n_crossing    = pcros ? pcros->GetVal() : 0;
    ps.p_det_acc     = (ps.n_laser_found > 0)
                       ? static_cast<double>(ps.n_accepted) / ps.n_laser_found
                       : 0.0;

    f->Close(); delete f;
    return true;
}

// ════════════════════════════════════════════════════════════════════════════
//  scanCacheFiles: trova tutti i file cache nella directory dati
// ════════════════════════════════════════════════════════════════════════════
static std::vector<PosStat> scanCacheFiles(const std::string& dataDir,
                                            int vbias_filter,
                                            double frac_filter)
{
    // Pattern: events_v3_vbias<V>_let<F>pe_cut<C>mhz_lthr<T>mV[_nofilt][_loose]_x<X>_y<Y>.root
    // oppure il nome che usa pos_scan internamente come suffix = _x<X>_y<Y>
    // Vediamo tutti i file events_v3_* e filtriamo per vbias e let

    static const std::regex re(
        R"(^events_v3_vbias(\d+)_let([\d\.]+)pe_cut(\d+)mhz_lthr([\d\.]+)mV.*_x([\-\d]+)_y([\-\d]+)\.root$)");

    std::vector<PosStat> result;

    void* dp = gSystem->OpenDirectory(dataDir.c_str());
    if (!dp) {
        std::cerr << "[ERR] Cannot open: " << dataDir << "\n";
        return result;
    }
    const char* ent = nullptr;
    while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
        std::string fn(ent);
        std::smatch m;
        if (!std::regex_match(fn, m, re)) continue;
        try {
            int    vb  = std::stoi(m[1].str());
            double fr  = std::stod(m[2].str());
            double x   = std::stod(m[5].str());
            double y   = std::stod(m[6].str());

            if (vbias_filter >= 0 && vb != vbias_filter) continue;
            if (frac_filter   >= 0 && std::abs(fr - frac_filter) > 0.001) continue;

            PosStat ps;
            ps.x          = x;
            ps.y          = y;
            ps.vbias      = vb;
            ps.frac_pe    = fr;
            ps.cache_path = dataDir + "/" + fn;

            if (readCacheStats(ps.cache_path, ps))
                result.push_back(ps);
        } catch (...) {}
    }
    gSystem->FreeDirectory(dp);

    // Ordina per (x, y)
    std::sort(result.begin(), result.end(),
        [](const PosStat& a, const PosStat& b){
            return a.x != b.x ? a.x < b.x : a.y < b.y;
        });
    return result;
}

// ════════════════════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════════════════════
void sipm_stat_equalizer()
{
    std::cout << "\n+==========================================================+\n"
              << "|  SiPM STATISTICA EQUALIZZATA — equalizzatore run griglia |\n"
              << "+==========================================================+\n\n";

    // ── 0. Data directory ────────────────────────────────────────────────────
    std::string dataDir;
    {
        std::string root = DATA_DIR;
        while (!root.empty() && root.back()=='/') root.pop_back();
        const size_t sl = root.find_last_of("/\\");
        if (sl != std::string::npos) root = root.substr(0, sl);

        std::vector<std::string> subs;
        void* dp = gSystem->OpenDirectory(root.c_str());
        if (dp) {
            const char* ent = nullptr;
            while ((ent = gSystem->GetDirEntry(dp)) != nullptr) {
                std::string s(ent);
                if (s=="."||s=="..") continue;
                FileStat_t st;
                if (gSystem->GetPathInfo((root+"/"+s).c_str(),st)==0
                    && R_ISDIR(st.fMode)) subs.push_back(s);
            }
            gSystem->FreeDirectory(dp);
        }
        std::sort(subs.begin(), subs.end());

        if (!subs.empty()) {
            std::cout << "  Cartelle disponibili in " << root << ":\n";
            for (size_t i=0; i<subs.size(); ++i)
                std::cout << "    [" << (i+1) << "] " << subs[i]
                          << (i==0 ? "   <-- default" : "") << "\n";
            const std::string l = readLineOrEmpty(
                "\n  Scegli [n] o percorso  [ENTER = 1]: ");
            if (l.empty()) dataDir = root+"/"+subs[0];
            else {
                try {
                    size_t pos=0; int idx=std::stoi(l,&pos);
                    if (pos==l.size() && idx>=1 && idx<=(int)subs.size())
                        dataDir = root+"/"+subs[idx-1];
                    else dataDir = l;
                } catch(...) { dataDir = l; }
            }
        } else {
            dataDir = readLine("  Percorso cartella dati: ");
        }
        while (!dataDir.empty() && dataDir.back()=='/') dataDir.pop_back();
        std::cout << "  --> " << dataDir << "\n\n";
    }

    // ── 1. Vbias e LET ───────────────────────────────────────────────────────
    int    vbias_sel = -1;
    double frac_sel  = -1.0;
    {
        const std::string sv = readLine("  Vbias [es. 57, ENTER = tutti]: ");
        if (!sv.empty()) try { vbias_sel = std::stoi(sv); } catch(...) {}
        const std::string sf = readLine("  LET frac_pe [es. 0.6, ENTER = tutti]: ");
        if (!sf.empty()) try { frac_sel = std::stod(sf); } catch(...) {}
    }

    // ── 2. Scansiona cache ───────────────────────────────────────────────────
    std::cout << "\n  Scansiono cache in " << dataDir << " ...\n";
    auto stats = scanCacheFiles(dataDir, vbias_sel, frac_sel);

    if (stats.empty()) {
        std::cerr << "[ERR] Nessuna cache trovata. Hai già girato sipm_pos_scan?\n";
        std::cerr << "      Verifica che i file si chiamino: events_v3_..._x<X>_y<Y>.root\n";
        return;
    }
    std::cout << "  Trovate " << stats.size() << " posizioni.\n\n";

    // ── 3. Punto di riferimento (90, 0) ──────────────────────────────────────
    // Se non trovato esattamente, usa il punto più vicino
    double ref_x = 90.0, ref_y = 0.0;
    {
        const std::string sr = readLine(
            "  Punto di riferimento (x y) [default 90 0]: ");
        std::istringstream ss(sr);
        double a=0, b=0;
        if (ss>>a>>b) { ref_x=a; ref_y=b; }
    }

    const PosStat* ref_pt = nullptr;
    double minDist = 1e18;
    for (const auto& ps : stats) {
        double d = std::hypot(ps.x - ref_x, ps.y - ref_y);
        if (d < minDist) { minDist = d; ref_pt = &ps; }
    }
    if (!ref_pt) { std::cerr << "[ERR] Nessun punto trovato.\n"; return; }

    std::cout << "\n  Punto di riferimento trovato:\n"
              << std::fixed << std::setprecision(1)
              << "    (x=" << ref_pt->x << ", y=" << ref_pt->y << ")"
              << "  Vbias=" << ref_pt->vbias << " V"
              << "  LET=" << std::setprecision(2) << ref_pt->frac_pe << " p.e.\n"
              << "    n_accepted=" << ref_pt->n_accepted
              << "  n_laser=" << ref_pt->n_laser_found
              << "  P_det_acc=" << std::setprecision(2)
              << ref_pt->p_det_acc * 100.0 << "%\n";

    // ── 4. N_target ──────────────────────────────────────────────────────────
    long N_target = 0;
    {
        const std::string sn = readLine(
            Form("\n  N eventi accettati TARGET per posizione [default = %ld (=quello del ref)]: ",
                 ref_pt->n_accepted));
        if (sn.empty()) N_target = ref_pt->n_accepted;
        else try { N_target = std::stol(sn); }
             catch(...) { N_target = ref_pt->n_accepted; }
    }
    if (N_target <= 0) {
        std::cerr << "[ERR] N_target deve essere > 0.\n"; return;
    }
    std::cout << "  --> N_target = " << N_target << " eventi accettati\n\n";

    // ── 5. Calcolo e tabella ──────────────────────────────────────────────────
    const int W = 78;
    std::cout << std::string(W, '=') << "\n";
    std::cout << std::left
              << std::setw(8)  << "x(mm)"
              << std::setw(8)  << "y(mm)"
              << std::setw(10) << "P_det%"
              << std::setw(12) << "N_acc_ora"
              << std::setw(14) << "N_wave_need"
              << std::setw(8)  << "N_run"
              << std::setw(10) << "t_run(m)"
              << "Note\n";
    std::cout << std::string(W, '-') << "\n";

    long total_runs  = 0;
    int  n_warn_low  = 0;

    for (const auto& ps : stats) {
        double p = ps.p_det_acc;

        long n_waveforms_needed = 0;
        long n_runs_needed      = 0;
        std::string note        = "";

        if (p <= 0.0) {
            // P_det = 0: impossibile stimare
            n_waveforms_needed = -1;
            n_runs_needed      = -1;
            note               = "P_det=0 — skip";
        } else {
            // N_waveform = ceil(N_target / P_det)
            n_waveforms_needed = static_cast<long>(
                std::ceil(static_cast<double>(N_target) / p));
            // N_run = ceil(N_waveform / 100k)
            n_runs_needed = static_cast<long>(
                std::ceil(static_cast<double>(n_waveforms_needed) / WAVEFORMS_PER_RUN));

            // Segnala se la stima è basata su poca statistica
            if (ps.n_laser_found < 10000) {
                note = "! stat bassa (n_laser=" + std::to_string(ps.n_laser_found) + ")";
                ++n_warn_low;
            } else if (p < 0.05) {
                note = "! P_det<5% — stima incerta";
            }
        }

        total_runs += std::max(0L, n_runs_needed);
        double t_min = n_runs_needed > 0
                       ? n_runs_needed * MINUTES_PER_RUN : 0.0;

        // Evidenzia il punto di riferimento
        bool is_ref = (&ps == ref_pt);

        std::cout << std::fixed << std::left
                  << std::setw(8)  << std::setprecision(1) << ps.x
                  << std::setw(8)  << std::setprecision(1) << ps.y
                  << std::setw(10) << std::setprecision(2) << p*100.0
                  << std::setw(12) << ps.n_accepted
                  << std::setw(14) << (n_waveforms_needed>0
                                       ? std::to_string(n_waveforms_needed) : "N/A")
                  << std::setw(8)  << (n_runs_needed>0
                                       ? std::to_string(n_runs_needed) : "N/A")
                  << std::setw(10) << std::setprecision(1) << t_min
                  << (is_ref ? "[REF] " : "")
                  << note << "\n";
    }

    // ── 6. Sommario ──────────────────────────────────────────────────────────
    std::cout << std::string(W, '=') << "\n\n";
    double total_hours = total_runs * MINUTES_PER_RUN / 60.0;
    std::cout << "  N posizioni:    " << stats.size() << "\n"
              << "  N run totali:   " << total_runs << "\n"
              << "  Tempo totale:   " << std::fixed << std::setprecision(1)
              << total_runs * MINUTES_PER_RUN << " min  ("
              << total_hours << " ore)\n"
              << "  Waveform/run:   " << WAVEFORMS_PER_RUN << "\n"
              << "  N_target/pos:   " << N_target << " eventi accettati\n";

    if (n_warn_low > 0)
        std::cout << "\n  ATTENZIONE: " << n_warn_low
                  << " posizioni hanno statistica bassa — stime incerte.\n"
                  << "  Suggerimento: gira prima una run esplorativa da 100k\n"
                  << "  su quelle posizioni per migliorare la stima di P_det.\n";

    // ── 7. Consiglio sulla notte ──────────────────────────────────────────────
    std::cout << "\n  Buio nautico a Bologna (fine maggio): ~22:00 - 04:20 (~6.3 ore)\n"
              << "  Buio astronomico completo:             ~23:00 - 03:35 (~4.5 ore)\n";
    if (total_hours <= 4.5)
        std::cout << "  --> Griglia FATTIBILE nel buio astronomico completo.\n";
    else if (total_hours <= 6.3)
        std::cout << "  --> Griglia FATTIBILE nel buio nautico (22:00-04:20).\n";
    else
        std::cout << "  --> Griglia NON fattibile in una notte ("
                  << std::setprecision(1) << total_hours
                  << " ore > 6.3 ore disponibili).\n"
                  << "     Considera di ridurre N_target o il numero di posizioni.\n";

    std::cout << "\n+==========================================================+\n"
              << "|  DONE\n"
              << "+==========================================================+\n\n";
}
