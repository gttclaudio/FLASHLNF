# Cavity coupling calculator from analytical modes

This repository contains Python code for computing transverse-traceless (TT) gravitational-wave, axion, scalar and dark photon coupling for ideal cylindrical, spherical, and rectangular cavities, using analytical expressions for cavity modes.

The implementation evaluates the directional coupling response of cavity eigenmodes to an incident gravitational wave by scanning over propagation angles $\beta$ and $\phi$, computing both GW polarization states (plus (parallel) and cross), or the coupling averaged over $\beta$ and $\phi$ for isotropically distributed scalar momenta or dark photon directions.

The code supports automatic mode normalization, parallel evaluation, and export of coupling maps in a standard format.

## Repository structure

```text
.
├── geometry/
│   ├── base.py
│   ├── cylindrical.py
│   ├── spherical.py
│   ├── rectangular.py
│   └── integration.py
├── modes/
│   ├── base.py
│   ├── cylindrical.py
│   ├── spherical.py
│   └── rectangular.py
├── coupling/
│   ├── coupling.py
│   ├── sources.py
│   └── utils.py
├── coupling_calculator.py  
├── results/                  # created locally; ignored by Git
├── requirements.txt
└── README.md
```

## Installation

Create a virtual environment and install dependencies:

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

## Run examples

### Cylindrical cavity

Compute GW coupling for the TM010 mode:

```bash
python coupling_calculator.py \
    --geometry cylindrical \
    --mode TM010 \
    --R 0.05 \
    --L 0.05 \
    --source "gw" \
    --n-beta 151 \
    --n-phi 51 \
    --n-processes 8 \
    --output-dir results
```

### Spherical cavity

```bash
python coupling_calculator.py \
    --geometry spherical \
    --mode TM010 \
    --R 0.05 \
    --source "gw" \
    --n-beta 151 \
    --n-phi 51 \
    --n-processes 8
```

### Rectangular cavity

```bash
python coupling_calculator.py \
    --geometry rectangular \
    --mode TE101 \
    --a 0.10 \
    --b 0.05 \
    --c 0.08 \
    --source "gw" \
    --n-beta 151 \
    --n-phi 51 \
    --n-processes 8
```

## Parameters

### Geometry selection

```text
--geometry
```

Supported values:

* `cylindrical`
* `spherical`
* `rectangular`

### Source selection

```text
--source
```

Supported values:

* `gw`
* `dp`
* `axion`
* `scalar`

**Axion** couplings only compute a single value. **GW** coupling averages over propagation direction. **DP** coupling (dark photon) averages over dark photon direction. **Scalar** coupling allows for an extra parameter --wavenumber, which takes in the wavenumber of the scalar in $\mathrm{m}^{-1}$, allowing for DM or non-DM scalar coupling calculations.

### Mode selection

```text
--mode TM010
```

Mode names follow standard cavity notation:

```text
TM010
TM011
TE101
...
```

Internally, cylindrical and spherical cavities evaluate both degenerate mode families automatically:

```text
TMa
TMb
```

while rectangular cavities evaluate a single mode family.

### Geometry parameters

Cylindrical:

```text
--R radius [m]
--L cavity length [m]
```

Spherical:

```text
--R radius [m]
```

Rectangular:

```text
--a x-dimension [m]
--b y-dimension [m]
--c z-dimension [m]
```

### Simulation parameters

```text
--n-beta        Number of β samples
--n-phi         Number of φ samples
--n-processes   Number of parallel processes
```

## Output

Results are stored automatically into geometry-specific folders:

```text
results/
├── cylindrical/
│   └── TT_gauge_TM010_129.0000MHz.pkl
├── spherical/
│   └── TT_gauge_TM011_340.0000MHz.pkl
└── rectangular/
    └── TT_gauge_TE101_218.3000MHz.pkl
```

Each output file contains a pickled pandas DataFrame.

Columns:

```text
beta
phi
coupling_parallel
coupling_cross
```

Example:

| beta | phi  | coupling_parallel | coupling_cross |
| ---- | ---- | ----------------- | -------------- |
| 0.00 | 0.00 | 0.031             | 0.017          |
| 0.31 | 0.00 | 0.042             | 0.024          |

Additional metadata are stored in `DataFrame.attrs`.

Common metadata:

```text
geometry
mode
frequency_mhz
n_beta
n_phi
```

For rectangular cavities:

```text
a
b
c
```

For spherical cavities:
```text
R
```

For cylindrical cavities:
```text
R
L
```
