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
./analyze.exe <..._wLumi.root> <ntuple.root>
```

Reads the `_wLumi` file through `lumimessenger.h`, keeps only zero-bias-triggered events (`isZeroBias*`) and writes a TNtuple with `instLumi, lumiLeveled, ZDCsumPlus, ZDCsumMinus, HFEMaxPlus_eta5, HFEMaxMinus_eta5, nVtx, nTrackInAcceptanceHP`. Per-forest ntuples (HIForward0-23) are in `rootfiles/tuple_*_wLumi.root`; merge with `hadd` if you use several.

## 5. Define the luminosity bins (only when the dataset changes)

`params.h` holds 10 equal-population `instLumi` bins (`params::lumibins`). Recompute the quantiles and edit `params.h` if you change dataset or selection.

## 6. Plots and fits

Both `plotAndFit.exe` and `plotHF.exe` take the ntuple, a `TAG`, and an optional `leveled` flag:

- `-1` (default): no cut
- `0`: non-leveled lumisections only (output tag gets `_notleveled`)
- `1`: leveled lumisections only (output tag gets `_leveled`)

Output goes to `outputs/<TAG>[_leveled|_notleveled]/` (`Main.root` + PDFs).

### ZDC

```bash
cd ZDCStudies
./plotAndFit.exe <ntuple.root> <TAG> [leveled]     # ZDC spectra per lumi bin (0nAn / An0n), fits
./plotMeanSigmaVsLumi.exe <TAG>                    # reads outputs/<TAG>/Main.root -> mean_vs_lumi.pdf, sigma_vs_lumi.pdf
```

Run `plotAndFit.exe` first. Give `plotMeanSigmaVsLumi.exe` the full output tag, including any `_leveled`/`_notleveled` suffix.

### HF

```bash
cd HFStudies
./plotHF.exe <ntuple.root> <TAG> [leveled]         # HF energy spectra per lumi bin, fraction above 16 vs lumi
```

## Example (HIForward0)

```bash
TAG=20260928_Dzero_260426-yrefmva_PbPbUPC_HIForward0_Dpt-2_wLumi
./analyze.exe <forest>_wLumi.root ntuple.root   # or reuse rootfiles/tuple_*_wLumi.root
cd ZDCStudies && ./plotAndFit.exe ../ntuple.root $TAG 1 && ./plotMeanSigmaVsLumi.exe ${TAG}_leveled
cd ../HFStudies && ./plotHF.exe ../ntuple.root $TAG 1
```
