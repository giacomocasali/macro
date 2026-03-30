# caendgz-sipmanalysis — Macro ROOT per analisi SiPM

Raccolta di macro C++ per ROOT per l'analisi di segnali SiPM (Silicon Photomultiplier)
acquisiti con oscilloscopio digitale a 5 GS/s e chip TDC ALCOR.

---

## Struttura del progetto

```
caendgz-sipmanalysis/analysis/
├── header/                        ← header condivisi (Utils.h, Config.h, ecc.)
├── 01_threshold_scan/             ← scan della soglia, calibrazione guadagno
├── 02_timing/                     ← analisi temporale LET e CFD
├── 03_waveform/                   ← analisi forme d'onda
├── 04_multi_scan/                 ← overlay di scan multipli
└── 05_alcor/                      ← analisi chip TDC ALCOR (ToT)
```

---

## Come eseguire le macro

**Regola fondamentale:** lanciare ROOT sempre dalla cartella `analysis/`,
non dall'interno delle sottocartelle. Gli `#include` degli header
sono scritti relativamente a questa posizione.

```bash
cd caendgz-sipmanalysis/analysis
root -l '01_threshold_scan/sipm_threshold_scan_full.cpp+'
```

Per le macro interattive (che chiedono parametri a runtime):

```bash
root -l
.L 02_timing/sipm_timing_unfiltered_v4_adaptive.cpp+
sipm_timing_unfiltered_v4_adaptive()
```

Per uscire da ROOT in qualsiasi momento: `Ctrl+C` o tasto `Esc` nelle macro
che supportano il loop di interruzione.

---

## 01 — Threshold scan

### `sipm_threshold_scan_full.cpp`
**Funzione:** `sipm_threshold_scan_full()`
**Input:** tutti i file `data.vbias_{V}_L.root` nella directory corrente (auto-detect)
**Output:** `scan_vbiasV_luxL.root` + `scan_vbiasV_luxL.png` + `laser_timing_*.png`

La macro principale per la calibrazione del guadagno SiPM. Per ogni file trovato:
- chiede la frequenza di taglio del filtro passa-basso e la finestra di trigger
- esegue lo scan della soglia da -10 a 110 mV a passi di 0.1 mV con logica ad arming (50%)
- calcola la derivata -dN/dV e fittia il primo picco p.e. con una gaussiana
- produce canvas con scan, derivata, linee di riferimento a 0.5/1/1.5 p.e.
- salva i risultati (scan, derivata, fit) in un file ROOT di output
- opzionalmente produce il canvas di timing laser

```bash
root -l '01_threshold_scan/sipm_threshold_scan_full.cpp+'
# > LP cutoff [MHz]: 500
# > Trigger start [ns]: 45
# > Trigger end [ns]: 55
# > Laser canvas? [y/n]: y
```

---

### `sipm_threshold_scan.cpp`
**Funzione:** `sipm_threshold_scan()`
**Input:** tutti i file `data.vbias_{V}_L.root` nella directory corrente (auto-detect)
**Output:** `scan_vbiasV_filterL.root` + PNG del canvas

Versione modulare con logica separata in `header/Config.h`, `header/SignalProcessing.h`,
`header/Plotting.h`. Funzionalmente equivalente a `sipm_threshold_scan_full.cpp`
ma con un algoritmo di peak-finding migliorato (selezione del candidato per sigma fisica).

```bash
root -l '01_threshold_scan/sipm_threshold_scan.cpp+'
```

---

### `threshold_scan2.cpp`
**Funzione:** `threshold_scan2()`
**Input:** `data.vbias_{55}.root`
**Output:** `threshold_scan.png`, `threshold_scan.root`

Versione compatta e autonoma (nessun header esterno). Scan da -100 a +100 mV
con logica ad arming. Utile per test rapidi su un singolo file.

```bash
root -l '01_threshold_scan/threshold_scan2.cpp+'
# > Frequenza di taglio LP [MHz]: 500
```

---

### `sipm_LET_threshold_scan.cpp`
**Funzione:** `sipm_LET_threshold_scan()`
**Input:** `data/data.vbias_{40}.root`
**Output:** cartelle `scan_f0.2/`, `scan_f0.4/`, `scan_f0.6/`, `scan_f0.8/`
con PNG `01_overlay.png`, `02_slices.png`, `03_final_trend.png`

Scan automatico del LET su frazioni del guadagno (0.2, 0.4, 0.6, 0.8 × gain).
Per ogni frazione: calibra il gain con TSpectrum, costruisce la mappa 2D
Δt vs ampiezza, fittia le slice con q-Gaussiana asimmetrica, estrae
σ_stocastica e σ_costante della risoluzione temporale.

```bash
root -l '01_threshold_scan/sipm_LET_threshold_scan.cpp+'
```

---

### `sipm_LET_scan.cpp`
**Funzione:** `sipm_LET_scan()`
**Input:** `data/data.vbias_{40}.root`
**Output:** stesse cartelle di `sipm_LET_threshold_scan.cpp`

Versione rinominata di `sipm_LET_threshold_scan.cpp` con nome funzione aggiornato.
Contenuto identico.

---

## 02 — Timing analysis

Tutte le macro di questa famiglia:
1. costruiscono lo spettro di ampiezza e calibrano il guadagno (multi-gaussiana + linearità)
2. misurano il tempo di arrivo con discriminatore LET (Leading Edge Threshold)
3. misurano il tempo di arrivo con discriminatore CFD (Constant Fraction, 50% del picco, scan backwards)
4. producono mappe 2D Δt vs ampiezza, slice per p.e., trend di media e risoluzione

---

### `sipm_timing_windowed.cpp`
**Funzione:** `sipm_timing_windowed()`
**Input:** `data.vbias_{40}.root`
**Output:** `plot_01_spectrum.png` … `plot_14_resolution_CFD.png`

Versione originale con finestra di segnale fissa [45, 55] ns.
14 canvas salvati in PNG.

```bash
root -l '02_timing/sipm_timing_windowed.cpp+'
```

---

### `sipm_timing_unfiltered_v1.cpp`
**Funzione:** `sipm_timing_unfiltered_v1()`
**Input:** `data.vbias_{40}.root`
**Output:** `uf_plot_01_spectrum.png` … `uf_plot_14_resolution_CFD.png`

Prima versione unfiltered: nessuna finestra temporale, waveform completa.
Finestra di fit sul timing fissa a ±2 ns.

```bash
root -l '02_timing/sipm_timing_unfiltered_v1.cpp+'
```

---

### `sipm_timing_unfiltered_v2.cpp`
**Funzione:** `sipm_timing_unfiltered_v2()`
**Input:** `data.vbias_{40}.root`
**Output:** `plots/uf_plot_01_spectrum_linear.png` … `plots/uf_plot_15_resolution_CFD.png`

Aggiunge rispetto a v1: spettro in scala lineare, legenda monospaziata,
banda σ sulla linearità, vincoli positivi su σ_stocastica e σ_costante nel fit.
15 canvas salvati in `plots/`.

```bash
root -l '02_timing/sipm_timing_unfiltered_v2.cpp+'
```

---

### `sipm_timing_unfiltered_v3_adaptive.cpp`
**Funzione:** `sipm_timing_unfiltered_v3_adaptive()`
**Input:** `data.vbias_{54}.root`
**Output:** `f7_plot_01_spectrum.png` … `f7_plot_14_resolution_CFD.png`

Aggiunge rispetto a v2: finestra di fit adattiva sul timing (stima FWHM → ±3σ),
zoom display a ±5σ. 14 canvas.

```bash
root -l '02_timing/sipm_timing_unfiltered_v3_adaptive.cpp+'
```

---

### `sipm_timing_unfiltered_v4_adaptive.cpp`
**Funzione:** `sipm_timing_unfiltered_v4_adaptive()`
**Input:** `data.vbias_{55}.root`
**Output:** `f7_plot_01_spectrum.png` … `f7_plot_14_resolution_CFD.png`

Come v3 ma senza vincoli sui parametri del fit multi-gaussiana (constraints rimossi).
Versione più recente e raccomandato per dati a Vbias = 55 V.

```bash
root -l '02_timing/sipm_timing_unfiltered_v4_adaptive.cpp+'
```

---

### `sipm_timing_LET_vs_CFD.cpp`
**Funzione:** `sipm_timing_LET_vs_CFD()`
**Input:** `data.vbias_{40}.root`
**Dipendenze:** `header/Utils.h`, `header/gauss_stuff.h`
**Output:** canvas interattivi (non salvati automaticamente)

Confronto diretto LET vs CFD sulla stessa acquisizione. Produce:
- spettro di ampiezza con fit multi-gaussiana
- canvas di linearità
- mappe 2D Δt vs ampiezza per LET e CFD
- canvas comparativo finale con trend di media e risoluzione sovrapposti

```bash
cd caendgz-sipmanalysis/analysis
root -l '02_timing/sipm_timing_LET_vs_CFD.cpp+'
```

---

## 03 — Waveform analysis

### `sipm_risetime_analysis.cpp`
**Funzione:** `sipm_risetime_analysis()`
**Input:** `data.vbias_{40}.root`
**Dipendenze:** `header/Utils.h`
**Output:** canvas interattivi (non salvati automaticamente)

Misura il rise time (dal crossing LET al picco locale) per ogni evento.
Produce mappa 2D rise time vs ampiezza con linee di separazione per p.e.,
e proiezioni 1D con palette arcobaleno per 0–7 p.e.

```bash
cd caendgz-sipmanalysis/analysis
root -l '03_waveform/sipm_risetime_analysis.cpp+'
```

---

### `waveform_2D_alignment.C`
**Funzione:** `waveform_2D_alignment()`
**Input:** `data.vbias_{54}.root`
**Output:** canvas TH2D interattivo (non salvato)

Allinea tutte le waveform al tempo di trigger laser (primo crossing a 50 mV)
e le accumula in un istogramma 2D: asse X = tempo relativo al laser [ns],
asse Y = ampiezza corretta [mV]. Scala Z logaritmica. Utile per
visualizzare la forma media del segnale e identificare la finestra ottimale.

```bash
root -l '03_waveform/waveform_2D_alignment.C+'
```

---

### `waveform_lowpass.cpp`
**Funzione:** `waveform_lowpass()`
**Input:** `data.vbias_{55}.root`
**Output:** `waveform_entry{N}_run{K}_fc{F}MHz.png`

Strumento interattivo in due fasi:
1. **Preview**: estrae waveform casuali, mostra corretta in baseline, chiede conferma
2. **Analisi**: applica filtro Butterworth 4° ordine alla frequenza scelta,
   sovrappone raw e filtrata, salva PNG. Ripetibile con cutoff diversi.

```bash
root -l '03_waveform/waveform_lowpass.cpp+'
# Fase 1: y/n per accettare la waveform mostrata
# Fase 2: inserire cutoff in MHz (es. 500)
```

---

### `waveform_analysis.cpp`
**Funzione:** `waveform_analysis()`
**Input:** `data.vbias_{53}.root`
**Output:** `waveform_entry{N}_run{K}_fc{F}MHz_thr{T}mV.png` + `waveform_results.csv`

Analisi completa a due canali (ch1 + laser). Per ogni waveform selezionata estrae:
ampiezza di picco, integrale del impulso, rise time (10%–90%), fall time (90%–10%),
tempo di crossing alla soglia su entrambi i canali.
Esporta i risultati evento per evento in CSV.

```bash
root -l '03_waveform/waveform_analysis.cpp+'
# > LP cutoff [MHz]: 500
# > Threshold [mV]: 5
```

---

### `signal_classifier.cpp`
**Funzione:** `signal_classifier()`
**Input:** `data.vbias_{53}.root`
**Output:** `sc_01_snr_distribution.png` … `sc_06_signal_grid.png`

Classifica ogni waveform come segnale o rumore usando due criteri:
SNR = A_max / σ_noise ≥ soglia e t_picco in finestra temporale [t_lo, t_hi].
Produce distribuzioni di SNR e ampiezza per le due classi, mappa 2D A_max vs SNR,
confronto waveform migliore vs rumore, e griglia di waveform segnale filtrate.

```bash
root -l '03_waveform/signal_classifier.cpp+'
# > SNR threshold: 5
# > t_sig_lo [ns]: 45
# > t_sig_hi [ns]: 55
# > Cutoff [MHz]: 500
```

---

### `rising_edge_heatmap.cpp`
**Funzione:** `rising_edge_heatmap()`
**Input:** `data.vbias_{55}.root`
**Output:** `reh_01_raw_heatmap.png`, `reh_02_filtered_heatmap.png`, `reh_03_amplitude_spectrum.png`

Costruisce heatmap 2D del fronte di salita allineato al tempo CFD (50% del picco,
scan backwards dal picco). Solo eventi segnale (SNR ≥ soglia e picco in finestra).
Produce heatmap raw e filtrata sovrapposto, spettro di ampiezza degli eventi selezionati.

```bash
root -l '03_waveform/rising_edge_heatmap.cpp+'
# > SNR threshold: 5
# > t_sig_lo [ns]: 45
# > t_sig_hi [ns]: 55
# > Window pre-CFD [ns]: 10
# > Window post-CFD [ns]: 30
# > LP cutoff [MHz]: 500
```

---

## 04 — Multi-scan overlay

### `read_all_scans.cpp`
**Funzione:** `read_all_scans()`
**Input:** tutti i file `scan_vbias{V}_filter{F}.root` nella directory corrente
**Output:** `all_scans.png`, `all_derivs.png`

Legge i file ROOT prodotti da `sipm_threshold_scan_full.cpp` o `sipm_threshold_scan.cpp`
e sovrappone tutti i scan su due canvas:
- canvas 1: N(soglia) in scala logaritmica, colori distinti per file
- canvas 2: derivata -dN/dV normalizzata al proprio massimo

Range X adattivo (taglia le code a zero), leggenda automatica con Vbias e filtro.

```bash
root -l '04_multi_scan/read_all_scans.cpp+'
```

---

### `pulse_interarrival_time.cpp`
**Funzione:** `pulse_interarrival_time()`
**Input:** `calibration.root`
**Dipendenze:** `header/Utils.h`
**Output:** canvas TH1D interattivo (non salvato)

Misura il tempo di inter-arrivo tra i primi due impulsi di ogni waveform
(proxy del dark count rate). Usa una soglia fissa a 10 mV con dead time di 12 ns.
Istogramma con scala Y logaritmica.

```bash
cd caendgz-sipmanalysis/analysis
root -l '04_multi_scan/pulse_interarrival_time.cpp+'
```

---

## 05 — ALCOR TDC analysis

### `alcor_tot_analysis.cpp`
**Funzione:** `alcor_tot_analysis()`
**Input:** ROOT file con TTree `alcor` (rami `t_leading`, `t_trailing`) +
opzionale TTree `laser` (ramo `t_ref`), oppure CSV, oppure simulazione interna
**Output:** `plots/alcor_01a_tot_log.png` … `plots/alcor_15_summary.png`

Analisi completa del chip TDC ALCOR basata sul Time over Threshold (ToT).
Chiede interattivamente: LSB del TDC in ns (tipicamente 0.05 ns = 50 ps),
numero di picchi p.e. da fittare, presenza del canale laser.
Produce 15 canvas:
- spettro ToT (log e lineare), calibrazione ToT vs p.e.
- Δt raw, mappa 2D Δt vs ToT (time walk visibile)
- fit e correzione del time walk
- risoluzione temporale prima e dopo la correzione
- SNR, stabilità ToT, inter-arrival time, linearità del guadagno, pannello riassuntivo

```bash
root -l '05_alcor/alcor_tot_analysis.cpp+'
# > TDC LSB [ns]: 0.05
# > Number of p.e. peaks: 8
# > Laser reference? [y/n]: y
# > Input source [1=ROOT, 2=CSV, 3=simulation]: 1
# > File path: alcor_data.root
```

---

## Note tecniche

**Filtro passa-basso:** tutte le macro usano un filtro Butterworth IIR del 4° ordine
implementato come cascata di sezioni biquad. Frequenza di taglio tipica: 500 MHz
su acquisizioni a 5 GS/s (Nyquist = 2500 MHz).

**Correzione baseline:** calcolata come mediana dei campioni nella finestra pre-segnale
[0, 30] ns. La mediana è robusta agli spike di rumore impulsivo.

**Discriminatore CFD:** il crossing al 50% del picco viene cercato scansionando
all'indietro dal picco — questo garantisce di trovare il fronte di salita
e non un ricrossing spurio sul fronte di discesa.

**File di dati attesi:** i file root hanno la struttura:
- TTree `ch1`: rami `time[1024]` (ns) e `amplitude[1024]` (mV)
- TTree `laser`: stessa struttura, canale di riferimento laser
