// Sostituisce drawDeltaTHistogram in TOTPlotting.h
// Mostra TUTTO il range reale dei delta_t (inclusi i negativi),
// non solo la finestra di fit.
// Il range dell'asse X va dal minimo reale al massimo reale dei dati,
// con un margine del 5% su entrambi i lati.
// La finestra di fit [fit_lo, fit_hi] è mostrata come banda colorata.
static void drawDeltaTHistogram(const std::vector<TOTEvent>& events,
                                const std::string& tag,
                                double fit_lo, double fit_hi,
                                OutCtx& ctx,
                                bool corrected = false) {
    if (events.empty()) return;
    std::string suffix = corrected ? "_corr" : "";

    // Trova il range reale dei dati
    double dt_min =  1e9, dt_max = -1e9;
    for (auto& e : events) {
        if (e.delta_t < dt_min) dt_min = e.delta_t;
        if (e.delta_t > dt_max) dt_max = e.delta_t;
    }
    // Aggiungi margine del 5% su entrambi i lati
    double margin = std::max(0.05 * (dt_max - dt_min), 1.0);
    double x_lo = dt_min - margin;
    double x_hi = dt_max + margin;

    // Binning: 0.2 ns/bin sul range completo
    // (0.02 ns/bin sarebbe troppo fine per range molto estesi con delta_t negativi)
    int nBins = std::max(200, (int)std::round((x_hi - x_lo) / 0.2));

    TH1D* h = new TH1D(
        Form("hDeltaTFull%s_%s", suffix.c_str(), tag.c_str()),
        Form("#Deltat histogram%s   %s;#Deltat (ns);Counts",
             corrected ? " (corrected)" : "", tag.c_str()),
        nBins, x_lo, x_hi);
    h->SetDirectory(nullptr);
    for (auto& e : events)
        h->Fill(e.delta_t);

    TCanvas* c = new TCanvas(
        Form("cDeltaTFull%s_%s", suffix.c_str(), tag.c_str()),
        Form("#Deltat full range%s — %s", corrected ? " (corr)" : "", tag.c_str()),
        1100, 600);
    c->SetGrid();
    c->SetLeftMargin(PAD_LEFT);
    c->SetRightMargin(PAD_RIGHT);
    c->SetBottomMargin(PAD_BOTTOM);
    c->SetTopMargin(PAD_TOP);

    // Banda colorata per la finestra di fit
    double yMax = h->GetMaximum() * 1.15;
    h->SetMaximum(yMax);
    h->SetLineColor(kAzure+1); h->SetLineWidth(2);
    h->Draw("HIST");

    // Banda verde chiaro: finestra di fit
    TBox* bFit = new TBox(fit_lo, 0.0, fit_hi, yMax);
    bFit->SetFillColorAlpha(kGreen-9, 0.30);
    bFit->SetLineWidth(0);
    bFit->Draw("same");
    h->Draw("HIST same");  // ridisegna sopra la banda

    // Linee verticali ai bordi della finestra di fit
    TLine* lLo = new TLine(fit_lo, 0.0, fit_lo, yMax * 0.9);
    TLine* lHi = new TLine(fit_hi, 0.0, fit_hi, yMax * 0.9);
    lLo->SetLineColor(kGreen+2); lLo->SetLineStyle(2); lLo->SetLineWidth(2);
    lHi->SetLineColor(kGreen+2); lHi->SetLineStyle(2); lHi->SetLineWidth(2);
    lLo->Draw("same"); lHi->Draw("same");

    // Linea verticale a Δt = 0
    if (x_lo < 0.0 && x_hi > 0.0) {
        TLine* lZero = new TLine(0.0, 0.0, 0.0, yMax * 0.85);
        lZero->SetLineColor(kRed+1); lZero->SetLineStyle(3); lZero->SetLineWidth(1);
        lZero->Draw("same");
        delete lZero;
    }

    // Statistiche nel box
    long nNeg = 0, nPos = 0;
    for (auto& e : events) {
        if (e.delta_t < 0) ++nNeg;
        else ++nPos;
    }
    TPaveText* pt = new TPaveText(0.14, 0.72, 0.50, 0.88, "NDC");
    pt->SetBorderSize(1); pt->SetFillColor(0); pt->SetFillStyle(1001);
    pt->SetTextFont(42); pt->SetTextSize(0.034);
    pt->AddText(Form("N total = %d", (int)events.size()));
    pt->AddText(Form("N(#Deltat < 0) = %ld  (%.1f%%)",
                     nNeg, 100.0*nNeg/events.size()));
    pt->AddText(Form("Range: [%.1f, %.1f] ns", dt_min, dt_max));
    pt->Draw();

    TLegend* leg = new TLegend(0.55, 0.78, 0.93, 0.88);
    leg->SetBorderSize(1); leg->SetFillColor(0); leg->SetTextSize(0.033);
    leg->AddEntry(h,    "#Deltat (full range)", "l");
    leg->AddEntry(bFit, "Fit window", "f");
    if (x_lo < 0.0) leg->AddEntry(lLo, "#Deltat = 0", "l");  // lLo è il primo verde
    leg->Draw();

    c->Update(); c->Modified();
    ctx.savePNG(c, Form("delta_t_full%s_%s.png", suffix.c_str(), tag.c_str()));
    delete h;
    delete bFit;
    delete lLo;
    delete lHi;
    delete pt;
    delete leg;
    // delete c;  // FIX: canvas resta aperta
}
