#ifndef UTILS_H
#define UTILS_H

#include <TString.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <TH2D.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TKey.h>

class Utils {
public:
    static Double_t* setupBranch_dgz(TTree* tree, const char* branchName) {
        Double_t* arr = new Double_t[1024];
        tree->SetBranchAddress(branchName, arr);
        return arr;
    }

    static int FindMaxBin(TH1D* hist, double xMin, double xMax) {
        int binMin = hist->GetXaxis()->FindBin(xMin);
        int binMax = hist->GetXaxis()->FindBin(xMax);
        int maxBin = binMin;
        double maxVal = hist->GetBinContent(binMin);
        for (int i = binMin; i <= binMax; ++i) {
            double val = hist->GetBinContent(i);
            if (val > maxVal) {
                maxVal = val;
                maxBin = i;
            }
        }
        return maxBin;
    }

    static double CountEvents(TH1D* hist, int bin_start, int bin_end) {
        int nBins = hist->GetNbinsX();
        bin_start = std::max(1, std::min(bin_start, nBins));
        bin_end = std::max(1, std::min(bin_end, nBins));
        return hist->Integral(bin_start, bin_end);
    }

    // Reso STATIC e aggiornato per compatibilità C++11/14
    static std::vector<std::pair<std::vector<double>, std::vector<double>>> correctWaveforms(
        const std::vector<std::pair<std::vector<double>, std::vector<double>>>& waveforms,
        double pre_signal_end = 30.0, bool nocorr = false
    ) {
        std::vector<std::pair<std::vector<double>, std::vector<double>>> corrected_waveforms;
        corrected_waveforms.reserve(waveforms.size());

        for (const auto& entry : waveforms) {
            const std::vector<double>& time = entry.first;
            const std::vector<double>& amp  = entry.second;

            std::vector<double> pre_amp;
            pre_amp.reserve(amp.size());

            for (size_t i = 0; i < time.size(); ++i) {
                if (time[i] <= pre_signal_end)
                    pre_amp.push_back(amp[i]);
            }

            double offset = 0.0;
            if (!pre_amp.empty()) {
                std::vector<double> tmp = pre_amp;
                size_t mid = tmp.size() / 2;
                std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
                offset = tmp[mid]; 
            }

            std::vector<double> corrected_amp = amp;
            if(nocorr == true) {for (auto& a : corrected_amp) a -= offset;}

            corrected_waveforms.emplace_back(time, corrected_amp);
        }
        return corrected_waveforms;
    }

    static double median_error(const std::vector<double>& vec) {
        if (vec.empty()) return 0.0;
        std::vector<double> sorted = vec;
        std::sort(sorted.begin(), sorted.end());
        size_t n = sorted.size();
        double med = (n % 2 == 0) ? 0.5 * (sorted[n/2 - 1] + sorted[n/2]) : sorted[n/2];
        double sum_sq = 0.0;
        for (double v : vec) sum_sq += (v - med) * (v - med);
        double sigma = std::sqrt(sum_sq / (n - 1));
        return 1.253 * sigma / std::sqrt(n);
    }

    static TH1D* projectionY(TH2D* h2, double lower_x, double upper_x, const char* name = "_proj") {
        if (!h2) return nullptr;
        int binx1 = h2->GetXaxis()->FindBin(lower_x);
        int binx2 = h2->GetXaxis()->FindBin(upper_x);
        return h2->ProjectionY(name, binx1, binx2);
    }

    static bool readDataFile(const TString &fullPath, std::vector<double> &col1, std::vector<double> &col2) {
        std::ifstream infile(fullPath.Data());
        if (!infile) return false;
        std::string line;
        std::getline(infile, line);
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            double v1, v2;
            if (iss >> v1 >> v2) {
                col1.push_back(v1);
                col2.push_back(v2);
            }
        }
        return true;
    }

    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> computeDerivative(const std::vector<double> &x, const std::vector<double> &y) {
        std::vector<double> xmid, dydx, err;
        size_t N = std::min(x.size(), y.size());
        for (size_t i = 0; i + 1 < N; ++i) {
            double dx = x[i+1] - x[i];
            double dy = y[i+1] - y[i];
            xmid.push_back((x[i] + x[i+1]) / 2.0);
            dydx.push_back(-dy / dx);
            err.push_back(std::sqrt(std::max(0.0, -dy / dx)));
        }
        return std::make_tuple(xmid, dydx, err);
    }
};

#endif