# Reading Guide — SiPM TOT Analysis

## How to read this codebase

The active workflow is split into two main macros (`sipm_calibrate.cpp` and
`sipm_tot_analysis.cpp`) plus a set of headers in `header/`.
Each header has a narrow responsibility. The order below is a practical causal
order: each file mostly depends on files listed above it.

---

## Suggested Reading Order

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
ChunkedHistoFill.h    → riempimento istogrammi a chunk (risparmio RAM)
VbiasSummary.h        → overlay multi-LET e grafico σ vs Vbias
VbiasAnalysis_v2.h    → pipeline completa per un Vbias (low-RAM, no TChain)
```

`WaveformPlotter.h`, `PersistencePlot.h`, and `OutputManager.h` can be read at
any time. `ProgressBar.h` is terminal UX only.

> **Version note**: `VbiasAnalysis_v2.h` replaces `VbiasAnalysis.h` for the
> current low-RAM pipeline. The active entrypoint is `sipm_tot_analysis.cpp`.

---

## Runtime Flow

### `sipm_calibrate.cpp`
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

### `sipm_tot_analysis.cpp`
```
main
 ├─ per ogni Vbias: loadCalibration() [CalibIO.h]
 │    └─ se più cutoff disponibili → scelta interattiva
 └─ processOneVbias_v2()            [VbiasAnalysis_v2.h]
      ├─ estimateLaserThreshold()    [internal to VbiasAnalysis_v2.h]
      │    └─ 5% of median laser peak (20% is still used in FilterDiagnostics)
      ├─ per ogni LET threshold:
      │    ├─ cache hit => reuse events ROOT directly (fast path)
      │    ├─ cache miss => collectTOTEvents_fileByFile()
      │    │    └─ open/scan/close one ROOT file at a time (no TChain)
      │    │         ├─ laserTriggerTime()    [SignalProcessing.h]
      │    │         ├─ correctBaseline()     [SignalProcessing.h]
      │    │         │    └─ scarta se RMS > 4 mV (2 mV con filtro)
      │    │         ├─ butterworthLowPass()  [ButterworthFilter.h]
      │    │         ├─ amp_max nella finestra trigger
      │    │         ├─ computeTOT()          [TOTAnalysis.h]
      │    │         └─ estimatePE()          [TOTAnalysis.h]
      │    ├─ applyTimeWalkCorrection()       [TimingCorrection.h]
      │    └─ analyseOneLET_chunked()
      │         ├─ build TH2/TH1 from cache
      │         ├─ save 4 active canvases:
      │         │    1) `tot_map`
      │         │    2) `tot_map_zoom` (Y fixed to 45..47 ns)
      │         │    3) `tot_xprojection`
      │         │    4) `tot_projection` (q-Gaussian fit)
      │         └─ run sideband sigma extraction (numerical, plot currently off)
      └─ drawVbiasSummary()                  [VbiasSummary.h]
```

---

## Main Differences v1 -> v2 (low-RAM)

| Aspetto | v1 (VbiasAnalysis.h) | v2 (VbiasAnalysis_v2.h) |
|---|---|---|
| Lettura dati | `TChain` (tutti i file in RAM) | Un file alla volta, apri/chiudi |
| RAM tipica | ~14 GB | ~2 GB |
| Laser threshold | 20% median peak (FilterDiagnostics) | 5% median peak (local estimate) |
| `edge_thr_frac` | fisso | adattivo: 0.5 con filtro, 100 senza |
| Baseline RMS cut | assente | 4 mV (raw) / 2 mV (filtrato) |
| Histogram fill | all in memory | chunked fill (`ChunkedHistoFill.h`) |
| Cache format | serialized vector | ROOT TTree with `laser_thr_saved` |

---

## Critical Parameters in `Config.h`

| Costante           | Valore | Significato                                        |
|--------------------|--------|----------------------------------------------------|
| `DATA_DIR`         | path   | Update this if data path changes                    |
| `BASELINE_START`   | 0 ns   | Baseline window start (must be pre-signal)          |
| `BASELINE_END`     | 30 ns  | Baseline window end                                 |
| `MIN_POSITIVE_THR` | 5 mV   | Calibration peaks searched above this threshold     |

---

## Core Function: `computeTOT()`

Defined in `TOTAnalysis.h`. It takes a waveform and returns `{t_rise, t_fall}`
in ns, or `{-1, -1}` if event quality checks fail.

**Rejection logic (in order):**

1. **Edge stability**: first/last-edge medians must be below `edge_thr_frac*threshold`.

2. **Rising-edge hysteresis**: arm below `threshold*hyst_frac`, then fire above `threshold`.

3. **Pre-crossing check**: local median before `t_rise` must stay below `pre_check_frac*threshold`.

4. **Falling confirmation**:
   - candidate `t_fall` is the first threshold crossing below threshold,
   - accept only if at least `min_below` samples are below threshold in the next
     `confirm_window` samples (current defaults in pipeline: 50, 10).

5. **Post-fall and tail checks**:
   post-fall median and end-of-trace median must be compatible with recovery.

---

## Cache Files

Cache files avoid rescanning all waveforms on repeated runs.

**Calibration** (written by `sipm_calibrate`, read by `sipm_tot_analysis`):
```
calib_vbias<V>_cut<C>mhz.root
  └─ TTree "calib": gain, offset, laser_thr, cutoff_MHz, trig_start, trig_end
```

**Events** (written/read by `sipm_tot_analysis`):
```
events_vbias<V>_let<L>pe_cut<C>mhz_lthr<T>mV[_nofilt].root
  └─ TTree "events": tot, delta_t, amp_max, n_pe, laser_thr_saved
```
If key cache parameters change (Vbias, LET, cutoff, laser threshold, filter on/off),
filename changes and cache miss triggers recomputation.

> `laser_thr_saved` is stored in event cache for threshold consistency checks.

---

## Current Interactive Behavior Notes

The input handlers in `sipm_tot_analysis.cpp` are robust to malformed input and EOF.

| Problema | Sintomo | Fix applicato |
|---|---|---|
| Loop infinito su EOF/stdin chiuso | `getline` fallisce, `line` resta vuota, loop non termina | `readLine` controlla il valore di ritorno di `getline`; su EOF chiama `std::exit(1)` |
| `fit_lo`/`fit_hi` non inizializzati | Input non numerico lascia failbit e valori garbage | `readDouble` usa `std::stod` + verifica consumo stringa; ripete il prompt finché valido |
| Loop y/n bloccato su failbit | `cin >> ans` fallisce, failbit persiste, il loop non si sblocca | `readYN` legge una riga intera con `readLine`, poi controlla il primo carattere; nessun `cin >>` diretto |

`readLine`, `readDouble`, and `readYN` are local lambdas and do not modify the
physics/event-selection logic.

---

## Troubleshooting

### Calibration gain looks wrong
Check `calibration_vbias<V>.png`:
- Are derivative peaks roughly equally spaced?
- Is `[CAL] gain_raw estimate` in a reasonable range?

### Delta-t fit is broad or unstable
Use `tot_projection_*` and adjust `fit_lo/fit_hi`.
Remember: fit is q-Gaussian, seeded by a local Gaussian around the peak.

### Why no persistence plots?
The persistence plotting path is currently disabled in `VbiasAnalysis_v2.h` to
keep cache-hit runs fast and avoid rescanning raw files.

### Zoom map window changes unexpectedly
The active zoom map is hardcoded to Y in `45..47 ns`.
If plot does not reflect code edits, restart ROOT and rebuild the macro (`++`)
to avoid stale compiled symbols.

### RAM usage still high
v2 scans one file at a time (~tens of MB per input file). If needed, tune
chunk sizes in `ChunkedHistoFill.h` and ROOT cache size (`SetCacheSize`).
