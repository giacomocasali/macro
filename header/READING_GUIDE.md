# Guida alla lettura del codice — SiPM TOT Analysis

## Come leggere questo codice senza impazzire

Il codice è diviso in due macro (`sipm_calibrate.cpp` e `sipm_tot_analysis_v5.cpp`)
e una serie di header nella cartella `header/`. Ogni header fa una cosa sola.
L'ordine qui sotto è l'ordine causale: ogni file dipende solo da quelli sopra di lui.

---

## Ordine di lettura

```
Config.h              → costanti globali e DATA_DIR
ButterworthFilter.h   → filtro IIR passa-basso
SignalProcessing.h    → correzione baseline, crossing laser, finestra trigger
Calibration.h         → scan soglia → derivata → picchi → gain [mV/p.e.]
CalibIO.h             → salva/carica calibrazione su file ROOT
FilterDiagnostics.h   → diagnostica visiva del filtro + stima soglia laser
TOTAnalysis.h         → computeTOT(), estimatePE(), collectTOTEvents()
EventCache.h          → salva/carica vettore eventi su file ROOT
TimingCorrection.h    → correzione time walk (metodo empirico o per p.e.)
SidebandAnalysis.h    → stima fondo + fit S+B + sottrazione sideband
TOTPlotting.h         → tutte le canvas (mappa 2D, proiezione, slice, per p.e.)
VbiasSummary.h        → overlay multi-LET e grafico σ vs Vbias
VbiasAnalysis.h       → pipeline completa per un Vbias (usa tutto quanto sopra)
```

I due file di classe (`WaveformPlotter.h`, `PersistencePlotter.h`) e
`OutputManager.h` sono indipendenti e possono essere letti in qualsiasi momento.
`ProgressBar.h` è solo estetica per il terminale.

---

## Flusso di esecuzione

### sipm_calibrate.cpp
```
main
 └─ per ogni Vbias:
     ├─ legge il primo run file → campionamento fs_MHz
     ├─ calibrateSpectrum()         [Calibration.h]
     │    ├─ stima noise_rms baseline (ARM_LEVEL)
     │    ├─ event loop: discriminatore con ARM_LEVEL fisso
     │    ├─ derivata -dN/dV
     │    ├─ trova picchi, stima gain_raw dal gap 1→2 p.e.
     │    ├─ assegna n_pe = round((pos - offset_guess)/gain_raw)
     │    └─ fit lineare posizione_picco vs n_pe → gain, offset
     ├─ drawFilterDiagnostics()     [FilterDiagnostics.h]
     │    └─ stima soglia laser = 20% del picco laser mediano
     └─ saveCalibration()           [CalibIO.h]
```

### sipm_tot_analysis_v5.cpp
```
main
 ├─ per ogni Vbias: loadCalibration() o calibrateSpectrum()
 └─ processOneVbias()               [VbiasAnalysis.h]
      ├─ legge tutti i run → TChain ch1 + laser
      ├─ per ogni LET threshold:
      │    ├─ loadEvents() (cache) o collectTOTEvents_v4()
      │    │    └─ per ogni evento:
      │    │         ├─ laserTriggerTime()    [SignalProcessing.h]
      │    │         ├─ correctBaseline()     [SignalProcessing.h]
      │    │         ├─ butterworthLowPass()  [ButterworthFilter.h]
      │    │         ├─ amp_max nella finestra trigger
      │    │         ├─ computeTOT()          [TOTAnalysis.h]
      │    │         └─ estimatePE()          [TOTAnalysis.h]
      │    ├─ applyTimeWalkCorrection()       [TimingCorrection.h]
      │    ├─ fillTH2D(), drawTOTMap()        [TOTPlotting.h]
      │    ├─ drawGlobalProjection()          [TOTPlotting.h]
      │    ├─ drawSlicesAndTrend()            [TOTPlotting.h]
      │    └─ fitSplusB()                     [SidebandAnalysis.h]
      └─ drawVbiasSummary()                  [VbiasSummary.h]
```

---

## Parametri critici in Config.h

| Costante         | Valore | Significato                                      |
|------------------|--------|--------------------------------------------------|
| `DATA_DIR`       | path   | **Cambia questo** se sposti i dati               |
| `BASELINE_START` | 0 ns   | Inizio finestra baseline (deve essere pre-segnale)|
| `BASELINE_END`   | 30 ns  | Fine finestra baseline                           |
| `MIN_POSITIVE_THR` | 5 mV | I picchi della calibrazione sono cercati sopra questa soglia |

---

## La funzione più importante: computeTOT()

Sta in `TOTAnalysis.h`. Prende un waveform filtrato e restituisce `{t_rise, t_fall}`
in nanosecondi, oppure `{-1, -1}` se l'evento è scartato.

**Condizioni di scarto (nell'ordine in cui sono applicate):**

1. **Edge-stability**: la mediana dei primi/ultimi 100 campioni deve essere
   sotto `threshold * 0.5`. Scarta waveform con baseline sollevata.

2. **Rising-edge hysteresis**: il discriminatore si arma quando il segnale scende
   sotto `threshold * 0.5` (50% della soglia), poi triggera quando sale sopra
   `threshold`. Scarta code discendenti di impulsi precedenti.

3. **Pre-crossing check**: la mediana dei 10 campioni prima del crossing deve
   essere sotto `threshold * 0.3`. Secondo controllo che il segnale venga
   dalla baseline, non da una coda.

4. **Anti-ringing**: `t_fall` è accettato solo dopo 3 campioni consecutivi sotto
   la soglia. Evita di triggerare sul ringing appena sotto soglia.

5. **Post-fall baseline check** ← *il più importante per qualità dei dati*:
   la mediana di tutti i campioni da `t_fall` alla fine della finestra deve
   essere sotto `threshold * 0.5`. **Scarta pile-up** (secondo impulso prima
   che il primo sia decaduto), **fronti non allineati** (breve dip sotto soglia
   seguito da risalita), e qualunque evento in cui il segnale non torna
   a baseline nella finestra.

---

## Struttura dei file cache

I file cache evitano di rileggere tutti i waveform ad ogni analisi.

**Calibrazione** (scritta da `sipm_calibrate`, letta da `sipm_tot_analysis_v5`):
```
calib_vbias<V>_cut<C>mhz.root
  └─ TTree "calib": gain, offset, laser_thr, cutoff_MHz, trig_start, trig_end
```

**Eventi** (scritta e letta da `sipm_tot_analysis_v5`):
```
events_vbias<V>_let<L>pe_cut<C>mhz_lthr<T>mV[_nofilt].root
  └─ TTree "events": tot, delta_t, amp_max, n_pe  (uno per evento accettato)
```
Se cambi qualsiasi parametro (Vbias, LET, cutoff, soglia laser, filtro on/off)
il nome del file cambia e si ottiene automaticamente un cache miss → ricalcolo.

---

## Cosa fare se qualcosa non torna

### Il gain della calibrazione è sbagliato
Guarda il PNG `calibration_vbias<V>.png`. Pannello superiore: curva N(V).
Pannello inferiore: derivata -dN/dV con i picchi marcati. Controlla:
- I picchi sono equidistanti?
- Il print `[CAL] gain_raw estimate` è ragionevole (5–50 mV/p.e.)?
- I picchi usati nel fit (`[CAL] Peaks used for fit`) hanno i numeri `n=1,2,3,...` giusti?

### Il fit di Δt è larghissimo o non converge
Guarda prima `peak_finder_vbias<V>_let<L>pe.png`. Indica dove sta il picco.
Se il picco non è centrato nella tua finestra di fit → aggiusta `fit_lo`/`fit_hi`.

### Vedo eventi con TOT enormi o Δt fuori range
Controlla `persistence_vbias<V>_let<L>pe.png` e `waveforms_smallTOT_*.png`.
I waveform con TOT grande sono pile-up o segnali ad alta ampiezza. Il post-fall
check in `computeTOT()` dovrebbe averli già rimossi — se ne vedi ancora,
la finestra di trigger è troppo larga o la soglia LET è troppo bassa.

### La soglia laser sembra sbagliata
Guarda `laser_diagnostics_vbias<V>.png`. La linea verde è `recThr` (20% del picco
mediano), quella arancione è 10 mV (il vecchio default). Se la distribuzione
`t_laser` è rumorosa o bimodale con la linea verde, prova ad alzare `ARM_SIGMA`
in `FilterDiagnostics.h` (riga con `recThr = medLaserAmp * 0.20`).
