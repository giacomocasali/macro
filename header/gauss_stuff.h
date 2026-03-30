#include <TSpectrum.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <iostream>

class GaussStuff {
public:

    struct PeakParams {
        double mean;
        double sigma;
        double amplitude;
    };


    static TF1* multi_gauss(TH1F* hist, int n_peaks_requested, const TString& fit_name = "fit") {
        if (!hist)
            return nullptr;
        /*double integral = hist->Integral("width"); // considera il bin width
        if (integral > 0) {
            hist->Scale(1.0 / integral);
            std::cout << "[Info] Histogram normalized. Integral = " << hist->Integral("width") << std::endl;
        } else {
            std::cerr << "[Warning] Histogram integral = 0. Skipping normalization." << std::endl;
        }*/

        // Trova picchi con TSpectrum
        TSpectrum spectrum(n_peaks_requested);
        int found = spectrum.Search(hist, 1, "", 0.0001);
        if (found < 1) {
            std::cerr << "No peaks found.\n";
            return nullptr;
        }

        double* xpeaks = spectrum.GetPositionX();
        double* ypeaks = spectrum.GetPositionY();
        std::vector<PeakParams> peaks;

        double max_mV = 400;  // valore massimo ammissibile in mV

        for (int i = 0; i < found; ++i) {
            if (xpeaks[i] > max_mV) continue;  // ⬅️ taglio sulla posizione, non sull’altezza
            peaks.push_back({xpeaks[i], hist->GetBinWidth(1), ypeaks[i]});
        }

        if (peaks.empty()) {
            std::cerr << "No peaks found below max mV.\n";
            return nullptr;
        }



        // Ordina per posizione
        std::sort(peaks.begin(), peaks.end(), [](const PeakParams& a, const PeakParams& b) {
            return a.mean < b.mean;
        });

        // Fit gaussiano individuale per stimare parametri
        for (auto& pk : peaks) {
            double range_min = pk.mean - 1;
            double range_max = pk.mean + 1;
            TF1 gausFit("tmpGaus", "gaus", range_min, range_max);
            gausFit.SetParameters(pk.amplitude, pk.mean, 0.5);
            hist->Fit(&gausFit, "RQ0"); // fit silenzioso
            pk.amplitude = gausFit.GetParameter(0);
            pk.mean = gausFit.GetParameter(1);
            pk.sigma = std::abs(gausFit.GetParameter(2));
        }

        // Calcola range di fit complessivo
        double fit_min = findLeftFitBoundary(hist, peaks.front().mean);
        double fit_max = peaks.back().mean + 3*peaks.back().sigma;

        // Multi-gaussiana: crea la funzione
        int n_gaus = peaks.size();
        TF1* fitFunc = new TF1(fit_name, [peaks, n_gaus](double *x, double *par) {
            double sum = 0;
            for (int i = 0; i < n_gaus; ++i) {
                double amp = par[i*3 + 0];
                double mean = par[i*3 + 1];
                double sigma = par[i*3 + 2];
                sum += amp * TMath::Gaus(x[0], mean, sigma, true);
            }
            return sum;
        }, fit_min, fit_max, n_gaus * 3);

        fitFunc->SetNpx(1000);

        for (int i = 0; i < n_gaus; ++i) {
            fitFunc->SetParameter(i*3 + 0, peaks[i].amplitude);
            fitFunc->SetParName(i*3 + 0, Form("A_{%d}", i+1));
            fitFunc->SetParameter(i*3 + 1, peaks[i].mean);
            fitFunc->SetParName(i*3 + 1, Form("#mu_{%d}", i+1));
            fitFunc->SetParameter(i*3 + 2, peaks[i].sigma);
            fitFunc->SetParName(i*3 + 2, Form("#sigma_{%d}", i+1));
            
        }
        
        // Indici della seconda gaussiana
        int idx_mean  = 1*3 + 1;   // 4
        int idx_sigma = 1*3 + 2;   // 5

       

        hist->Fit(fitFunc, "IMREQ");
        double mean  = fitFunc->GetParameter(idx_mean);
        double emean = fitFunc->GetParError(idx_mean);

        double sigma  = fitFunc->GetParameter(idx_sigma);
        double esigma = fitFunc->GetParError(idx_sigma);

        std::cout << "Seconda gaussiana:" << std::endl;
        std::cout << "Media  = " << mean  << " ± " << emean  << std::endl;
        std::cout << "Sigma  = " << sigma << " ± " << esigma << std::endl;

        return fitFunc;
    }



    
    static TF1* FitQGaussPlusConst(TH1* proj, const char* fitName) {
        // Creo la funzione di fit con nome passato come argomento
        TF1* qGaussPlusConst = new TF1(fitName,
            [](Double_t* x, Double_t* par) {
                Double_t par_q[5] = { par[0], par[1], par[2], par[3], par[4] };
                return qGaussAsym(x, par_q) + par[5];
            }, 30, 50, 6);

        qGaussPlusConst->SetNpx(1000);
        qGaussPlusConst->SetParNames("A", "mean", "sigma", "q1", "q2", "offset");

        // Stima iniziale parametri
        int maxbin = proj->GetMaximumBin();
        double A0 = proj->GetMaximum();
        double sum = 0, m = 0, s2 = 0;
        for (int b = maxbin - 5; b <= maxbin + 5; ++b) {
            double x = proj->GetBinCenter(b);
            double y = proj->GetBinContent(b);
            sum += y; m += x * y;
        }
        m /= sum;
        for (int b = maxbin - 5; b <= maxbin + 5; ++b) {
            double x = proj->GetBinCenter(b);
            double y = proj->GetBinContent(b);
            s2 += y * TMath::Power(x - m, 2);
        }
        double sigma0 = TMath::Sqrt(s2 / sum);

        double offset0 = 0;
        int bmin = proj->FindBin(30);
        int bmax = proj->FindBin(33);
        for (int b = bmin; b <= bmax; ++b)
            offset0 += proj->GetBinContent(b);
        offset0 /= (bmax - bmin + 1);

        qGaussPlusConst->SetParameters(A0, m, sigma0, 1.0, 1.3, offset0);
        qGaussPlusConst->FixParameter(3, 1.0);
        qGaussPlusConst->SetParLimits(4, 1.0, 1.5);
        const double offset_threshold = 1e-3;
        if (offset0 < offset_threshold) {
            qGaussPlusConst->FixParameter(5, 0.0);
        } else {
            qGaussPlusConst->SetParLimits(5, 0.0, 1000);
        }

        proj->Fit(qGaussPlusConst, "IMREQ");

        return qGaussPlusConst;
    }
    
    static TF1* FitQGauss(TH1* proj, const char* fitName) {
        // Creo la funzione di fit con nome passato come argomento
        TF1* qGauss = new TF1(fitName,
            [](Double_t* x, Double_t* par) {
                // par[0] = A, par[1] = mean, par[2] = sigma, par[3] = q1, par[4] = q2
                return qGaussAsym(x, par);
            }, 30, 50, 5);

        qGauss->SetNpx(1000);
        qGauss->SetParNames("A", "mean", "sigma", "q1", "q2");

        // Stima iniziale parametri
        int maxbin = proj->GetMaximumBin();
        double A0 = proj->GetMaximum();
        double sum = 0, m = 0, s2 = 0;
        for (int b = maxbin - 5; b <= maxbin + 5; ++b) {
            double x = proj->GetBinCenter(b);
            double y = proj->GetBinContent(b);
            sum += y;
            m += x * y;
        }
        m /= sum;
        for (int b = maxbin - 5; b <= maxbin + 5; ++b) {
            double x = proj->GetBinCenter(b);
            double y = proj->GetBinContent(b);
            s2 += y * TMath::Power(x - m, 2);
        }
        double sigma0 = TMath::Sqrt(s2 / sum);

        // Setta parametri iniziali
        qGauss->SetParameters(A0, m, sigma0, 1.0, 1.3);

        // puoi decidere se fissare q1=1 come nel tuo Fit precedente
        qGauss->FixParameter(3, 1.0);

        proj->Fit(qGauss, "IMREQ");

        return qGauss;
    }


    static std::tuple<double,double,double,double,TF1*> q_gauss_projection(TH1D* hProj){
        if (!hProj || hProj->GetEntries() == 0) return {0,0,0,0,nullptr};

        double mean_guess  = hProj->GetMean();
        double sigma_guess = hProj->GetRMS();
        double A_guess     = hProj->GetMaximum();

        TF1* fFit = new TF1("qGaussFitQ1Fixed", qGaussAsym,
                            hProj->GetXaxis()->GetXmin(),
                            hProj->GetXaxis()->GetXmax(), 5);

        fFit->SetParameters(A_guess, mean_guess, sigma_guess, 1.0, 1.3);

        // 🔹 Assegna nomi ai parametri
        fFit->SetParName(0, "A");
        fFit->SetParName(1, "#mu");
        fFit->SetParName(2, "#sigma");
        fFit->SetParName(3, "q_{1}");
        fFit->SetParName(4, "q_{2}");

        // 🔹 Fissa q1 = 1
        fFit->FixParameter(3, 1.0);
        fFit->SetNpx(1000);
        hProj->Fit(fFit,"IMRE");

        double mean      = fFit->GetParameter(1);
        double mean_err  = fFit->GetParError(1);
        double sigma     = fFit->GetParameter(2);
        double sigma_err = fFit->GetParError(2);

        return {mean, mean_err, sigma, sigma_err, fFit};
    }


    static std::tuple<double,double,double,double,TF1*> gauss_fit_projection(TH1D* hProj) {
        if (!hProj || hProj->GetEntries() == 0) return {0,0,0,0,nullptr};

        double mean_guess  = hProj->GetMean();
        double sigma_guess = hProj->GetRMS();
        double A_guess     = hProj->GetMaximum();

        TF1* fFit = new TF1("gaussFit", "gaus",
                            hProj->GetXaxis()->GetXmin(),
                            hProj->GetXaxis()->GetXmax());
        fFit->SetParameters(A_guess, mean_guess, sigma_guess);

        // 🔹 Assegna nomi ai parametri
        fFit->SetParName(0, "A");
        fFit->SetParName(1, "#mu");
        fFit->SetParName(2, "#sigma");

        fFit->SetNpx(1000);
        hProj->Fit(fFit, "IMREQ");

        double mean      = fFit->GetParameter(1);
        double mean_err  = fFit->GetParError(1);
        double sigma     = fFit->GetParameter(2);
        double sigma_err = fFit->GetParError(2);

        return {mean, mean_err, sigma, sigma_err, fFit};
    }

    
private:

    static double findLeftFitBoundary(TH1F* hist, double mean) {
        int meanBin = hist->FindBin(mean);
        double prevContent = hist->GetBinContent(meanBin);
        int counter_up = 0;

        for (int bin = meanBin - 1; bin > 1; --bin) {
            double content = hist->GetBinContent(bin);
            if (content > prevContent) {
                counter_up++;
                if (counter_up > 1) {
                    return hist->GetBinCenter(bin + 1);
                }
            } else {
                counter_up = 0;
            }
            prevContent = content;
        }
        return hist->GetBinLowEdge(1);
    }
    
    static Double_t qGaussAsym(Double_t* x, Double_t* par) {
        double xx    = x[0];
        double A     = par[0];
        double mean  = par[1];
        double sigma = par[2];
        double q1    = par[3];
        double q2    = par[4];

        if (sigma <= 0) return 0;
        const double eps = 1e-6;

        if (xx <= mean) {
            if (TMath::Abs(q1 - 1.0) < eps) {
                double arg_sq = TMath::Power((xx - mean) / sigma, 2);
                return A * TMath::Exp(-0.5 * arg_sq);
            } else {
                double arg = 1 - (1 - q1) * (1. / (3 - q1)) * TMath::Power((xx - mean) / sigma, 2);
                if (arg <= 0) return 0;
                return A * TMath::Power(arg, 1. / (1 - q1));
            }
        } else {
            double arg = 1 - (1 - q2) * (1. / (3 - q2)) * TMath::Power((xx - mean) / sigma, 2);
            if (arg <= 0) return 0;
            return A * TMath::Power(arg, 1. / (1 - q2));
        }
    }

    
};
