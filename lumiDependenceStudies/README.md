# Luminosity-dependence studies

Studies of how the ZDC and HF observables change with instantaneous luminosity (`instLumi`, in units of 10^33 cm^-2 s^-1) for the Dzero PbPb UPC forests. Run the steps below **in order** from this folder (`lumiDependenceStudies/`).

Overview:

```
forest(s) --runlumi_list--> run/lumi list --(omstools instlumi.py)--> instlumi.csv
forest(s) + instlumi.csv --addlumi--> forest + instLumi branches --analyze--> ntuple
ntuple --> ZDCStudies/plotAndFit, plotMeanSigmaVsLumi   and   HFStudies/plotHF
```

## 0. Build

```bash
make                      # runlumi_list.exe, addlumi.exe, analyze.exe
make -C ZDCStudies        # plotAndFit.exe, plotMeanSigmaVsLumi.exe
make -C HFStudies         # plotHF.exe
```

The Makefiles need ROOT (`root-config`) and the shared headers in `../include/` (`xjjanauti.h`, ...).

## 1. List the run/lumisection pairs in the input forest

```bash
./runlumi_list.exe <forest.root[,forest2.root,...]> runLumi_dir/runLumi_inp_<TAG>.txt
```

Scans the `Tree` of the forest(s) and writes one `<run> <lumi>` pair per line (unique, sorted). Existing lists are kept in `runLumi_dir/`, e.g. `runLumi_inp_Dzero_260426-yrefmva_PbPbUPC_HIForward0_Dpt-2.txt`.

## 2. Download the instantaneous luminosity per lumisection (OMS)

Uses the `filling_scheme_mod` branch of omstools: <https://github.com/kovbalcsab/omstools/tree/filling_scheme_mod>

```bash
git clone -b filling_scheme_mod https://github.com/kovbalcsab/omstools.git
cd omstools
pip3 install -r requirements.txt
# put your OMS API credentials in env.py:
#   CLIENT_ID = '...'
#   CLIENT_SECRET = '...'
python3 instlumi.py --inputtxt <path>/runLumi_dir/runLumi_inp_<TAG>.txt [--outcsv <out.csv>]
```

- Input: the run/lumi list from step 1 (`<run> <ls>` per line).
- Default output: `outcsv/instlumi.csv` (sorted by run, then ls). `init_lumi`, `end_lumi`, `avg_lumi` (= (init+end)/2) are in 10^33 cm^-2 s^-1.
- `addlumi.cc` reads columns `0..4` (run, ls, init, end, avg) and column `8` (the `lumi_trend` column; `"leveled"` marks a lumi-leveled lumisection, `null` values are read as 0). Check that your csv has this column, otherwise everything is treated as not leveled.
- Needs network access to OMS; run where you have credentials (e.g. lxplus).

## 3. Attach the luminosity to the forest

```bash
./addlumi.exe <forest.root> <instlumi.csv> rootfiles/<name>_wLumi.root [overwrite=0]
```

Copies the forest (the input is never modified) and adds the branches `instLumiInit`, `instLumiEnd`, `instLumi` (avg) and `instLumiLeveled` to `Tree`. Lumisections missing from the csv get `instLumi = -1`. It refuses to overwrite an existing output unless `overwrite = 1`.


## 4. Make the analysis ntuple

```bash
./analyze.exe <..._wLumi.root> <ntuple.root> [triggerSelection=0] [goodlumi.json]
./analyze.exe --help
```

Reads the `_wLumi` file through `lumimessenger.h`, keeps only events passing `triggerSelection` and writes a TNtuple with `instLumi, lumiLeveled, ZDCsumPlus, ZDCsumMinus, HFEMaxPlus_eta5, HFEMaxMinus_eta5, nVtx, nTrackInAcceptanceHP, Run, Lumi`. Per-forest ntuples (HIForward0-23) are in `rootfiles/tuple_*_wLumi.root`; merge with `hadd` if you use several.

`triggerSelection` (optional integer, default 0; anything else is rejected):

| value | selection |
|---|---|
| `0` | any zero-bias UPC path (`isZeroBias*`) |
| `-1` | `isNotBptxOR` or `isUnpairedBunchBptxMinus` or `isUnpairedBunchBptxPlus` |
| `-2` | `isUnpairedBunchBptxMinus` |
| `-3` | `isUnpairedBunchBptxPlus` |
| `-4` | `isNotBptxOR` |

It also refuses identical input/output paths and warns if no event passes.

## 5. Define the luminosity bins (only when the dataset changes)

`params.h` holds 10 equal-population `instLumi` bins (`params::lumibins`). Recompute the quantiles and edit `params.h` if you change dataset or selection.

## 6. Plots and fits

Both `plotAndFit.exe` and `plotHF.exe` take the ntuple, a `TAG`, and an optional `leveled` flag (3rd argument):

- `-1` (default): no cut
- `0`: non-leveled lumisections only (output tag gets `_notleveled`)
- `1`: leveled lumisections only (output tag gets `_leveled`)

An optional 4th argument `noZDCcut` (default 0): `1` applies no ZDC value cut (every event fills both the An0n and 0nAn histograms, so they are identical; in `plotHF` all HF histograms use the wide range) and appends `_noZDCcut` to the output tag. Default is 0.

Output goes to `outputs/<TAG>[_leveled|_notleveled][_noZDCcut]/` (`Main.root` + PDFs).

### ZDC

```bash
cd ZDCStudies
./plotAndFit.exe <ntuple.root> <TAG> [leveled] [noZDCcut]     # ZDC spectra per lumi bin (0nAn / An0n), fits
./plotMeanSigmaVsLumi.exe <TAG>                    # reads outputs/<TAG>/Main.root -> mean_vs_lumi.pdf, sigma_vs_lumi.pdf
```

Run `plotAndFit.exe` first. Give `plotMeanSigmaVsLumi.exe` the full output tag, including any `_leveled`/`_notleveled`/`_noZDCcut` suffix.

### HF

```bash
cd HFStudies
./plotHF.exe <ntuple.root> <TAG> [leveled] [noZDCcut]         # HF energy spectra per lumi bin, fraction above 16 vs lumi
```

## Example (HIForward0)

```bash
TAG=20260928_Dzero_260426-yrefmva_PbPbUPC_HIForward0_Dpt-2_wLumi
./analyze.exe <forest>_wLumi.root ntuple.root   # or reuse rootfiles/tuple_*_wLumi.root
cd ZDCStudies && ./plotAndFit.exe ../ntuple.root $TAG 1 && ./plotMeanSigmaVsLumi.exe ${TAG}_leveled
cd ../HFStudies && ./plotHF.exe ../ntuple.root $TAG 1
```

## 7. Lumi-dependent ZDC recalibration (prototype)

Corrects ZDCsum Plus/Minus as a function of instLumi so that the noise core (q90/q95 of the 0n side), the 1n peak (mean and width) and the 2n peak are the same as in the lowest-lumi bin. Monotone cubic (PCHIP, C1) map per side with knots `0, n90, n95, m1-s1, m1, m1+s1, m2`, each knot position a pol2 in instLumi (`ZDCStudies/zdcLumiCorrection.h`).

```bash
cd ZDCStudies
make calibZDC.exe plotAndFit.exe
./calibZDC.exe <ntuple.root> <TAG>               # UNCORRECTED data, non-leveled -> ./zdcLumiCorrParams.h + outputs/calib_<TAG>/
make plotAndFit.exe                              # rebuild so the new parameters are compiled in
./plotAndFit.exe <ntuple.root> <TAG> 0 0 1       # applyCorr=1 (5th arg, smooth map) -> outputs/<TAG>_notleveled_corr/ (map_*, slope_*, *_ratio_to_bin0 pdfs)
./plotAndFit.exe <ntuple.root> <TAG> 0 0 2       # applyCorr=2: legacy piecewise-linear map -> ..._corrlin/ (for comparison)
./plotMeanSigmaVsLumi.exe <TAG>_notleveled_corr
```

Always derive the parameters from the uncorrected ntuple (`calibZDC` never applies the correction). `analyze.cc` / the ntuple are untouched; apply with `zdccorr::correct(zdc, instLumi, zdccorr::kPlus|kMinus)`.
