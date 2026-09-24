# Cavity coupling calculator

This repository contains Python code for computing transverse-traceless (TT) gravitational-wave coupling from COMSOL cavity-field CSV data (where it is assumed the magnetic field points in the z-direction). In addition, it can compute couplings to axion and dark photon particles from the same data.

Two full mode files are included in `data/`:

- `TE011`: `FLASH_LF_TE011_218.3MHz.csv`
- `TM010`: `FLASH_LF_TM010_129.0MHz.csv`

## Repository structure

```text
.
├── src/
│   ├── GW_coupling_TT_E_D_fields.py
│   ├── coupling.py
│   └── gw_coupling_tt.py
├── data/
│   ├── FLASH_LF_TE011_218.3MHz.csv
│   └── FLASH_LF_TM010_129.0MHz.csv
├── results/                 # created locally; ignored by Git
├── requirements.txt
├── .gitignore
└── README.md
```

`GW_coupling_TT_E_D_fields.py` and `gw_coupling_tt.py` are old versions, where the same functionality is now implemented in `coupling.py`

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Run one selected mode

Run the TE011 mode:

```bash
python src/gw_coupling_tt.py --mode TE011 --n-angles 201 --n-processes 6 --output-dir results
```

Run the TM010 mode:

```bash
python src/gw_coupling_tt.py --mode TM010 --n-angles 201 --n-processes 6 --output-dir results
```

For a quick test, reduce the angular grid:

```bash
python src/gw_coupling_tt.py --mode TE011 --n-angles 5 --n-processes 1 --output-dir results
```

`--n-angles 201` means a `201 x 201` scan over beta and phi. This is the production-style setting. `--n-angles 5` is only a quick sanity check.

### Other couplings

Run coupling calculation for axions to the TM010 mode:

```bash
python src/coupling.py --mode TM010 --source axion --ouput-dir results
```

Run coupling calculation for dark photons to the TM010 mode:

```bash
python src/coupling.py --mode TM010 --source dp --ouput-dir results
```

The TT-gauge calculation is also possible in the same script:

```bash
python src/coupling.py --mode TM010 --source gw --ouput-dir results
```

## What `--mode` means

`--mode` selects one of the predefined CSV files in `data/`:

```python
MODE_FILES = {
    "TE011": "FLASH_LF_TE011_218.3MHz.csv",
    "TM010": "FLASH_LF_TM010_129.0MHz.csv",
}
```

So `--mode TE011` is equivalent to explicitly running:

```bash
python src/gw_coupling_tt.py --csv data/FLASH_LF_TE011_218.3MHz.csv
```

This is different from the older `--filename-filter` option. A filename filter only searches for matching filenames inside a folder. For this repository, `--mode` is clearer because the two supported modes are already known.

## Advanced explicit CSV usage

You can still provide a CSV file directly:

```bash
python src/gw_coupling_tt.py --csv data/FLASH_LF_TM010_129.0MHz.csv --n-angles 201 --n-processes 6
```

The folder/filter mode is kept only for batch processing many CSV files:

```bash
python src/gw_coupling_tt.py --folder data --filename-filter TE011 --n-angles 201 --n-processes 6
```

## Input CSV format

The CSV file should be exported from COMSOL with eight header rows. After the header, the code expects these columns:

```text
x, y, z, Ex, Ey, Ez, normE, normD
```

Coordinates are assumed to be in millimeters. Complex field values may use `i` as the imaginary unit.

The filename should contain the mode name and frequency, for example:

```text
FLASH_LF_TE011_218.3MHz.csv
FLASH_LF_TM010_129.0MHz.csv
```

The script extracts `TE011` or `TM010` as the mode name and `218.3 MHz` or `129.0 MHz` as the frequency.

## Output

The script saves a pickle file in the output directory:

```text
results/TT_gauge_TE011_218.3000MHz.pkl
results/TT_gauge_TM010_129.0000MHz.pkl
```

Each output file is a dictionary with two keys, `"results"` and `"summary"`:

- `"results"` is a pandas DataFrame with columns:

  ```text
  beta, phi, coupling_parallel, coupling_cross
  ```

- `"summary"` is a dictionary of scan-level statistics:

  ```text
  mean_p, max_p, mean_c, max_c, mean, max
  ```

  where `mean_p`/`max_p` and `mean_c`/`max_c` are the mean and max of the parallel and cross couplings respectively, and `mean`/`max` are the mean and max of the combined (parallel + cross) coupling.

### Output from `coupling.py`

For an axion:

```text
results/axion_TE011_218.3000MHz.pkl
results/axion_TM010_129.0000MHz.pkl
```

Each output file is a dictionary with two keys, `"results"` and `"summary"`:

- `"results"` is a pandas DataFrame with a single column:

  ```text
  coupling
  ```

- `"summary"` is a dictionary:

  ```text
  mean, max
  ```

  Since the axion coupling is a single value rather than a scan, `mean` and `max` are both equal to that value.

For a dark photon:

```text
results/DP_TE011_218.3000MHz.pkl
results/DP_TM010_129.0000MHz.pkl
```

Each output file is a dictionary with two keys, `"results"` and `"summary"`:

- `"results"` is a pandas DataFrame with columns:

  ```text
  beta, phi, coupling
  ```

- `"summary"` is a dictionary:

  ```text
  mean, max
  ```

  the mean and max coupling over the beta/phi scan.

### Loading the output

```python
import pickle

with open("results/TT_gauge_TM010_129.0000MHz.pkl", "rb") as f:
    data = pickle.load(f)

df = data["results"]        # pandas DataFrame with the full scan
summary = data["summary"]   # dict of mean/max coupling values
```
