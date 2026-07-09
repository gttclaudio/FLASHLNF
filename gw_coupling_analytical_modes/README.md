# Cavity coupling calculator from analytical modes

This repository contains Python code for computing transverse-traceless (TT) gravitational-wave, axion, scalar and dark photon coupling for ideal cylindrical, spherical, and rectangular cavities, using analytical expressions for cavity modes.

The implementation evaluates the directional coupling response of cavity eigenmodes to an incident gravitational wave by scanning over propagation angles $\beta$ and $\phi$, computing both GW polarization states (plus (parallel) and cross), or the coupling averaged over $\beta$ and $\phi$ for isotropically distributed scalar momenta or dark photon directions.

The code supports automatic mode normalization, parallel evaluation, and export of coupling maps in a standard format.

## Definitions of coupling coefficients
In the following expressions we define $\vec{E}_n(\vec{x})$ and $\vec{B}_n(\vec{x})$ the electric and magnetic field, respectively, for the cavity mode $n$. The cavity is taken to have a volume $V$ and have a homogeneous magnetic field in the $z$-direction $\vec{B} = B_0 \hat{z}$. Angles $\beta$ and $\phi$ are taken to represent the polar angles **towards which** a particular vector is pointing.

### Axions
$$
    C_\mathrm{axion}= \frac{1}{B_0^2 V} \frac{ \left| \int \enspace \vec{E}_n(\vec{x}) \cdot \vec{B} \enspace  dV\right|^2}{\int |\vec{E}_n(\vec{x})| ^2 \enspace dV},
$$

where for the axion case $\beta$ and $\phi$ are not parameters, and only a single value of $C_\mathrm{axion}$ is computed.

### Gravitational waves in TT gauge

$$
    C_\mathrm{gw}^{+, \times} (\beta, \phi) = \frac{1}{B_0^2 V} \frac{ \left| \int \enspace \vec{E}_n(\vec{x}) \cdot \hat{j}^{\mathrm{+, \times}}_\mathrm{eff} (\vec{x}; \beta, \phi) \enspace dV \right|^2}{ \int |\vec{E}_n(\vec{x})|^2 \enspace dV},
$$

where 

$$
\hat{j}^\mathrm{+,\times}_\mathrm{eff} (\vec{x}; \beta, \phi) = \frac{\vec{j}^\mathrm{+,\times}_\mathrm{eff} (\vec{x}; \beta, \phi)}{\omega_0 h_{+, \times}} = e^{-i \vec{k}(\beta, \phi) \cdot \vec{x}} \vec{B}_{+, \times}(\beta, \phi),
$$

for $\vec{k}(\beta, \phi)$ describing the wave vector describing a gravitational wave propagating in the $\beta, \phi$ direction, and $\vec{B}_{+, \times}(\beta, \phi)$ defined as

$$
\mathbf{B}_{+} = (\hat{B}_{\perp} \cdot \hat{e}_2) \mathbf{e}_1 + (\vec{B}_{\perp} \cdot \hat{e}_1) \hat{e}_2,
$$

$$
\mathbf{B}_{\times} = -(\vec{B}_{\perp} \cdot \hat{e}_1) \hat{e}_1 + (\vec{B}_{\perp} \cdot \hat{e}_2) \hat{e}_2,
$$

$$
B_{\parallel} = (\vec{B} \cdot \hat{k}) \hat{k}, \quad B_{\perp} = \mathbf{B} - (\mathbf{B} \cdot \hat{k}) \hat{k}.
$$

Where $\hat{k}(\beta, \phi)$ is the unit vector of the direction of gravitational wave propagation, and $\hat{e}_1(\beta, \phi)$ and $\hat{e}_2(\beta, \phi)$ are perpendicular unit vectors, which define the plus and cross polarisations.

### Dark photons

$$
    C_\mathrm{dp} (\beta, \phi) = \frac{1}{V} \frac{ \left| \int \enspace \vec{E}_n(\vec{x}) \cdot \hat{e}_A(\beta, \phi) \enspace  dV\right| ^2}{\int |\vec{E}_n(\vec{x})|^2 \enspace dV},
$$

where $\hat{e}_A(\beta, \phi)$ is the unit vector describing the polarisation direction of the vector field $A$, parametrised by polar angles $\beta$ and $\phi$. The script also computes the average over angles $\beta, \phi$ for an isotropic dark photon field.


### Scalars

$$
    C_\mathrm{scalar} (\beta, \phi) = \frac{1}{B_0^2 V} \frac{ \left| \int \enspace \vec{B}_n(\vec{x}) \cdot \vec{B} \enspace e^{i \vec{k} (\beta, \phi)\cdot \vec{x}} \enspace  dV\right|^2}{\int |\vec{B}_n(\vec{x})|^2 \enspace dV},
$$

where $\vec{k} (\beta, \phi)$ is the wave vector of the scalar particle. Then $\beta, \phi$ are averaged over for an isotropic momentum distribution. The phase term is important to include for scalars, because even for scalars constituting dark matter, the constant contribution from $e^{i \vec{k} (\beta, \phi) \cdot \vec{x}} \approx  1 + \vec{k} (\beta, \phi) \cdot \vec{x}$ leads to a 0 coupling.


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
