#pragma once
// ProgressBar.h
// Live ASCII progress bar for long event loops.
// Usage:
//   ProgressBar bar(nEntries, "my loop");
//   for (Long64_t i = 0; i < nEntries; ++i) {
//       if (bar.update(i, nNoLaser, nNoTOT, nAccepted)) break;
//       ...
//   }
//   bar.done();

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <TSystem.h>
#include <TROOT.h>

struct ProgressBar {
    Long64_t total;
    std::string label;
    std::chrono::steady_clock::time_point t0;
    long nNoLaser = 0, nNoTOT = 0, nAccepted = 0;

    static const int PRINT_EVERY = 1000;
    static const int BAR_W       = 40;

    ProgressBar(Long64_t n, const std::string& lbl)
        : total(n), label(lbl), t0(std::chrono::steady_clock::now())
    {
        std::cout << "\n  +- " << label << "  (" << total << " events)\n"
                  << "  |" << std::flush;
    }

    // Call inside the event loop (checks interrupt every PRINT_EVERY events).
    // Returns true if the loop should be interrupted.
    bool update(Long64_t i, long noLaser, long noTOT, long accepted) {
        nNoLaser = noLaser; nNoTOT = noTOT; nAccepted = accepted;

        if (i % PRINT_EVERY != 0) return false;

        // ROOT interrupt check
        if (gROOT->IsInterrupted()) {
            std::cout << "\n  [INTERRUPT] Stopping at event " << i << "\n";
            return true;
        }
        gSystem->ProcessEvents();

        double frac    = (total > 0) ? (double)i / (double)total : 1.0;
        int    filled  = (int)(frac * BAR_W);
        auto   now     = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - t0).count();
        double eta_s   = (frac > 0.01 && frac < 0.999)
                         ? elapsed / frac * (1.0 - frac) : 0.0;

        std::cout << "\r  | [";
        for (int b = 0; b < BAR_W; ++b) std::cout << (b < filled ? '#' : '-');
        std::cout << "] "
                  << std::setw(5) << std::fixed << std::setprecision(1)
                  << frac * 100.0 << "%"
                  << "  acc="     << accepted
                  << "  noLas="   << noLaser
                  << "  noTOT="   << noTOT;
        if (eta_s > 1.0)
            std::cout << "  ETA " << (int)eta_s << "s";
        std::cout << "   " << std::flush;
        return false;
    }

    void done() {
        auto   now = std::chrono::steady_clock::now();
        double dt  = std::chrono::duration<double>(now - t0).count();
        std::cout << "\r  | [";
        for (int b = 0; b < BAR_W; ++b) std::cout << '#';
        std::cout << "] 100.0%"
                  << "  accepted=" << nAccepted
                  << "  noLaser="  << nNoLaser
                  << "  noTOT="    << nNoTOT
                  << "  time="     << std::fixed << std::setprecision(1) << dt << "s\n"
                  << "  +--------------------------------------\n" << std::flush;
    }
};
