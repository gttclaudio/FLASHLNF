"""
Compute transverse-traceless (TT) gravitational-wave coupling from the
``E_D_fields_FLASH_cavity`` CSV data set.

This is an adaptation of ``GW_coupling_TT.py`` (which reads the COMSOL exports in
``OLD_data``) to the newer ``E_D_fields_FLASH_cavity`` exports. The integration
physics is identical to the original script; only the CSV loader and the
frequency handling differ.

Differences handled here
------------------------
1. File format. The new CSVs have a single header row (no 8 metadata rows) and
   the columns

       mode_label, alfa, x_m, y_m, z_m,
       Re_Ex_Vpm, Im_Ex_Vpm, Re_Ey_Vpm, Im_Ey_Vpm, Re_Ez_Vpm, Im_Ez_Vpm,
       Eabs_Vpm, Re_Dx_Cpm2, Im_Dx_Cpm2, Re_Dy_Cpm2, Im_Dy_Cpm2,
       Re_Dz_Cpm2, Im_Dz_Cpm2, Dabs_Cpm2

   Coordinates are in METERS (not mm), and the electric / displacement fields are
   given as separate real and imaginary columns instead of COMSOL "i" strings.
   The loader converts these into the same schema the original integration code
   expects: ``x, y, z`` in millimetres and complex ``Ex, Ey, Ez, normE, normD``.

2. Regular-grid check. The new files are sampled on a regular grid restricted to
   the inside of the cavity (every retained (x, y) column carries the full set of
   z nodes). ``check_regular_grid`` verifies the uniform-spacing condition and is
   run for every file before integrating.

3. Frequency. The new filenames carry only the mode and tuner angle (e.g.
   ``TM010_tuners_0deg``) and NO frequency, but the TT phase factor needs
   ``k = 2*pi*f/c``. Frequencies are read from the spreadsheet
   ``Modes_freq_ FLASH.xlsx`` (a table of eigenfrequency [MHz] vs. mode and tuner
   position 0/45/90/135/180 deg). Each CSV's frequency is looked up by its
   (mode, tuner-angle) pair. Files with no matching entry in the spreadsheet are
   SKIPPED, and a summary listing the missing (mode, tuner) pairs is printed at
   the end.
"""

import argparse
import itertools
import multiprocessing as mp
import os
import re

import numpy as np
import pandas as pd
from tqdm import tqdm


EPSILON_0 = 8.854e-12
SPEED_OF_LIGHT = 2.998e8

# Spreadsheet of mode eigenfrequencies vs. tuner position, used to look up the
# frequency for each CSV file. See load_freq_table().
DEFAULT_FREQ_XLSX = "Modes_freq_ FLASH.xlsx"


def load_freq_table(xlsx_path):
    """
    Read the mode/tuner frequency spreadsheet into a lookup dict.

    The spreadsheet is laid out as a table whose rows are mode names (e.g.
    ``TM010``) and whose columns are tuner positions labelled ``0°``, ``45°``,
    ``90°``, ``135°``, ``180°``. Cells that are not numeric (e.g. ``-`` for a
    configuration that was not simulated) are skipped.

    Returns
    -------
    dict[(str, int), float]
        Mapping ``(mode, tuner_angle_deg) -> frequency [MHz]``.
    """
    raw = pd.read_excel(xlsx_path, header=None)

    # Locate the header row holding the tuner-position labels and which column
    # each tuner angle lives in.
    angle_pat = re.compile(r"^(\d+)\s*°$")
    header_row, tuner_cols = None, {}
    for i in range(len(raw)):
        cols = {
            col: int(angle_pat.match(str(raw.iat[i, col]).strip()).group(1))
            for col in range(raw.shape[1])
            if angle_pat.match(str(raw.iat[i, col]).strip())
        }
        if cols:
            header_row, tuner_cols = i, cols
            break
    if header_row is None:
        raise ValueError(f"Could not find a tuner-position header row in '{xlsx_path}'.")

    mode_pat = re.compile(r"^(TM|TE)\d{3}$")
    table = {}
    for i in range(header_row + 1, len(raw)):
        for col in range(raw.shape[1]):
            cell = str(raw.iat[i, col]).strip()
            if not mode_pat.match(cell):
                continue
            for ccol, angle in tuner_cols.items():
                try:
                    freq = float(raw.iat[i, ccol])
                except (ValueError, TypeError):
                    continue
                if np.isfinite(freq):
                    table[(cell, angle)] = freq
            break  # one mode per row
    return table


def parse_mode_tuner(mode_label):
    """
    Split a ``mode_label`` such as ``TM010_tuners_45deg`` into ``("TM010", 45)``.

    Returns ``(mode, tuner_angle_deg)`` or ``(mode, None)`` if the tuner angle
    cannot be parsed.
    """
    parts = mode_label.split("_")
    mode = parts[0]
    tuner = None
    for token in parts[1:]:
        m = re.match(r"^(\d+)deg$", token)
        if m:
            tuner = int(m.group(1))
            break
    return mode, tuner


# ---------------------------------------------------------------------------
# Grid weighting (unchanged from GW_coupling_TT.py)
# ---------------------------------------------------------------------------
def _axis_node_weights(unique_coords):
    """Return sorted grid nodes and trapezoidal node weights along one axis."""
    nodes = np.sort(np.unique(unique_coords))
    if len(nodes) < 2:
        return nodes, np.ones_like(nodes, dtype=float)

    spacing = np.diff(nodes)
    weights = np.empty_like(nodes, dtype=float)
    weights[0] = spacing[0] / 2.0
    weights[-1] = spacing[-1] / 2.0
    if len(nodes) > 2:
        weights[1:-1] = (spacing[:-1] + spacing[1:]) / 2.0
    return nodes, weights


def _local_z_node_weights_per_xy(df, z_nodes_global, z_weights_global):
    """
    Assign z-direction integration weights for each local (x, y) column.

    Some (x, y) positions may have only one sampled z point. In that case, the
    nearest global z-node weight is used as a fallback.
    """
    dz_out = np.empty(len(df), dtype=float)

    for (_, _), block in df.groupby(["x", "y"], sort=False):
        z_values = block["z"].to_numpy()
        row_index = block.index.to_numpy()

        z_sorted = np.sort(z_values)
        if len(z_sorted) >= 2:
            _, z_weights_local = _axis_node_weights(z_sorted)
            sort_order = np.argsort(z_values)
            dz_out[row_index[sort_order]] = z_weights_local
            continue

        nearest_index = np.searchsorted(z_nodes_global, z_values[0])
        if nearest_index == len(z_nodes_global) or (
            nearest_index > 0
            and abs(z_values[0] - z_nodes_global[nearest_index - 1])
            < abs(z_values[0] - z_nodes_global[nearest_index])
        ):
            nearest_index = max(nearest_index - 1, 0)
        dz_out[row_index] = z_weights_global[nearest_index]

    return dz_out


def _xy_node_weights(df):
    """Assign per-point trapezoidal integration weights in the x and y directions."""
    x_nodes, x_weights = _axis_node_weights(df["x"])
    y_nodes, y_weights = _axis_node_weights(df["y"])

    x_index = np.searchsorted(x_nodes, df["x"].to_numpy())
    y_index = np.searchsorted(y_nodes, df["y"].to_numpy())
    return x_weights[x_index], y_weights[y_index]


def f_tt(z_m, wave_number):
    """TT phase factor for a propagation coordinate z_m and wave number k."""
    return np.exp(1j * z_m * wave_number)


# ---------------------------------------------------------------------------
# TT-coupling integral (unchanged from GW_coupling_TT.py)
# ---------------------------------------------------------------------------
def compute_tt_coupling_from_data(df_raw, beta, freq_hz, phi=0.0):
    """
    Compute the two TT-coupling contributions for one pair of scan angles.

    Parameters
    ----------
    df_raw : pandas.DataFrame
        Field-map data with coordinates [mm] and complex electric/displacement
        fields (columns x, y, z, Ex, Ey, Ez, normE, normD).
    beta : float
        Polar angle of the GW propagation direction in radians.
    freq_hz : float
        Mode frequency in Hz.
    phi : float
        Azimuthal rotation angle in radians.

    Returns
    -------
    tuple[complex, complex]
        Normalized parallel-like and cross-like TT-coupling terms.
    """
    df = df_raw.copy()
    df = df.groupby(["x", "y", "z"], as_index=False).mean()

    dx_pt, dy_pt = _xy_node_weights(df[["x", "y", "z"]])
    z_nodes_global, z_weights_global = _axis_node_weights(df["z"])
    dz_pt = _local_z_node_weights_per_xy(df[["x", "y", "z"]], z_nodes_global, z_weights_global)

    volume_weights = dx_pt * dy_pt * dz_pt * 1e-9  # mm^3 -> m^3

    x_mm = df["x"].to_numpy()
    y_mm = df["y"].to_numpy()
    z_mm = df["z"].to_numpy()

    x = x_mm * 1e-3
    y = y_mm * 1e-3
    z = z_mm * 1e-3
    z_center_m = 0.5 * (z_mm.min() + z_mm.max()) * 1e-3
    z_centered = z - z_center_m

    ex = df["Ex"].to_numpy(complex)
    ey = df["Ey"].to_numpy(complex)
    ez = df["Ez"].to_numpy(complex)
    norm_e = df["normE"].to_numpy(complex)
    norm_d = df["normD"].to_numpy(complex)

    cavity_volume = float(np.sum(volume_weights))
    normalization_integral = np.sum(volume_weights * norm_e * np.conjugate(norm_d) / EPSILON_0)

    y_phi = y * np.cos(phi) - x * np.sin(phi)
    ex_phi = ex * np.cos(phi) + ey * np.sin(phi)
    ey_phi = ey * np.cos(phi) - ex * np.sin(phi)

    wave_number = 2.0 * np.pi * freq_hz / SPEED_OF_LIGHT
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)

    propagation_coordinate = sin_beta * y_phi + cos_beta * z_centered
    tt_phase = f_tt(propagation_coordinate, wave_number)

    j_eff_x = sin_beta * tt_phase
    j_eff_y = sin_beta * cos_beta * tt_phase
    j_eff_z = sin_beta**2 * tt_phase

    parallel_integral = np.sum(volume_weights * ex_phi * j_eff_x)
    cross_integral = np.sum(volume_weights * (ey_phi * j_eff_y - ez * j_eff_z))

    denominator = cavity_volume * normalization_integral
    coupling_parallel = np.abs(parallel_integral) ** 2 / denominator
    coupling_cross = np.abs(cross_integral) ** 2 / denominator

    return coupling_parallel, coupling_cross


def worker_angle(args):
    """Multiprocessing wrapper for one beta/phi grid point."""
    beta, df, freq_hz, phi = args
    return compute_tt_coupling_from_data(df, beta, freq_hz, phi)


# ---------------------------------------------------------------------------
# New-format CSV loading
# ---------------------------------------------------------------------------
def mode_label_from_filename(filename):
    """Recover the COMSOL ``mode_label`` (e.g. ``TM010_tuners_0deg``) from a filename."""
    base = os.path.basename(filename)
    return base.split("_E_D_field")[0]


def load_field_csv(csv_path):
    """
    Load one ``E_D_fields_FLASH_cavity`` CSV and return a DataFrame in the schema
    expected by ``compute_tt_coupling_from_data``: coordinates in millimetres and
    complex ``Ex, Ey, Ez, normE, normD``.
    """
    raw = pd.read_csv(csv_path)

    out = pd.DataFrame(
        {
            "x": raw["x_m"].to_numpy(float) * 1e3,  # m -> mm
            "y": raw["y_m"].to_numpy(float) * 1e3,
            "z": raw["z_m"].to_numpy(float) * 1e3,
            "Ex": raw["Re_Ex_Vpm"].to_numpy(float) + 1j * raw["Im_Ex_Vpm"].to_numpy(float),
            "Ey": raw["Re_Ey_Vpm"].to_numpy(float) + 1j * raw["Im_Ey_Vpm"].to_numpy(float),
            "Ez": raw["Re_Ez_Vpm"].to_numpy(float) + 1j * raw["Im_Ez_Vpm"].to_numpy(float),
            # Eabs / Dabs are real magnitudes; cast to complex to match the OLD
            # loader, which produced complex normE / normD.
            "normE": raw["Eabs_Vpm"].to_numpy(float).astype(complex),
            "normD": raw["Dabs_Cpm2"].to_numpy(float).astype(complex),
        }
    )
    return out.dropna()


def check_regular_grid(df, atol_mm=1e-3):
    """
    Verify that the (x, y, z) sampling is a regular grid.

    Returns a dict with, per axis, the number of unique nodes, the spacing, and
    whether the spacing is uniform within ``atol_mm``. Also reports whether every
    (x, y) column carries the same number of z nodes (the "inside cavity" regular
    grid condition).
    """
    report = {}
    uniform = True
    for axis in ("x", "y", "z"):
        nodes = np.sort(df[axis].unique())
        if len(nodes) < 2:
            report[axis] = {"n": len(nodes), "spacing": None, "uniform": False}
            uniform = False
            continue
        spacing = np.diff(nodes)
        axis_uniform = bool(np.allclose(spacing, spacing[0], atol=atol_mm))
        uniform = uniform and axis_uniform
        report[axis] = {
            "n": len(nodes),
            "spacing": float(spacing[0]),
            "uniform": axis_uniform,
            "range": (float(nodes[0]), float(nodes[-1])),
        }

    z_per_xy = df.groupby(["x", "y"])["z"].nunique()
    report["z_per_xy_min"] = int(z_per_xy.min())
    report["z_per_xy_max"] = int(z_per_xy.max())
    report["uniform_columns"] = bool(z_per_xy.min() == z_per_xy.max())
    report["is_regular"] = bool(uniform and report["uniform_columns"])
    return report


def print_grid_report(report):
    for axis in ("x", "y", "z"):
        info = report[axis]
        if info["spacing"] is None:
            print(f"  {axis}: n={info['n']} (single node)")
        else:
            lo, hi = info["range"]
            flag = "OK" if info["uniform"] else "NON-UNIFORM"
            print(
                f"  {axis}: n={info['n']:3d}  spacing={info['spacing']:.4f} mm  "
                f"range=[{lo:.1f}, {hi:.1f}] mm  [{flag}]"
            )
    print(
        f"  z per (x,y) column: min={report['z_per_xy_min']} max={report['z_per_xy_max']}  "
        f"-> regular grid: {report['is_regular']}"
    )


# ---------------------------------------------------------------------------
# Scan driver
# ---------------------------------------------------------------------------
def run_scan(
    folder_base="E_D_fields_FLASH_cavity",
    mode_filter="",
    z_min_mm=0.0,
    z_max_mm=1500.0,
    n_angles=201,
    n_processes=6,
    freq_xlsx=DEFAULT_FREQ_XLSX,
):
    """Run the TT-coupling scan for all matching CSV files in the input folder."""
    freq_table = load_freq_table(freq_xlsx)
    print(f"Loaded {len(freq_table)} (mode, tuner) frequencies from '{freq_xlsx}'.")

    input_files = sorted(
        filename
        for filename in os.listdir(folder_base)
        if filename.endswith(".csv")
        and (mode_filter is None or mode_filter == "" or mode_filter in filename)
    )

    print(f"Found {len(input_files)} CSV file(s) in '{folder_base}'.")
    skipped_modes = []

    for filename in input_files:
        mode_label = mode_label_from_filename(filename)
        mode, tuner = parse_mode_tuner(mode_label)

        if (mode, tuner) not in freq_table:
            print(f"[SKIP] {filename}: no frequency for (mode={mode}, tuner={tuner} deg) in '{freq_xlsx}'.")
            skipped_modes.append((mode, tuner))
            continue

        freq_hz = freq_table[(mode, tuner)] * 1e6
        print(f"\n[RUN ] {mode_label}: {mode} mode, tuner {tuner} deg, freq = {freq_hz / 1e6:.3f} MHz")

        df = load_field_csv(os.path.join(folder_base, filename))

        df = df[(df["z"] >= z_min_mm) & (df["z"] <= z_max_mm)].reset_index(drop=True)

        report = check_regular_grid(df)
        print_grid_report(report)
        if not report["is_regular"]:
            print(
                "  WARNING: grid is not strictly regular; trapezoidal weights are "
                "still applied but verify the sampling."
            )

        phis = np.linspace(0.0, 2.0 * np.pi, n_angles)
        betas = np.linspace(0.0, 2.0 * np.pi, n_angles)
        angle_pairs = list(itertools.product(betas, phis))
        tasks = [(beta, df, freq_hz, phi) for beta, phi in angle_pairs]

        with mp.Pool(processes=n_processes) as pool:
            results = list(tqdm(pool.imap(worker_angle, tasks), total=len(tasks)))

        records = [
            {"b": beta, "phi": phi, "num_p": coupling_p, "num_c": coupling_c}
            for (beta, phi), (coupling_p, coupling_c) in zip(angle_pairs, results)
        ]
        df_results = pd.DataFrame(records)

        output_path = os.path.join(
            folder_base,
            f"TT_gauge_{mode_label}_{freq_hz / 1e6:.4f}MHz.pkl",
        )
        df_results.to_pickle(output_path)
        print(f"  Saved: {output_path}")

    if skipped_modes:
        unique_skipped = sorted(set(skipped_modes), key=lambda mt: (mt[0], -1 if mt[1] is None else mt[1]))
        print("\n" + "=" * 70)
        print("Frequency information missing for the following (mode, tuner) pairs:")
        for mode, tuner in unique_skipped:
            print(f"  - {mode}, tuner {tuner} deg")
        print(
            f"\nAdd the eigenfrequency [MHz] for each pair to '{freq_xlsx}'\n"
            "(matching the mode row and tuner-position column), then re-run."
        )
        print("=" * 70)


def parse_args():
    """Parse command-line options for batch processing."""
    parser = argparse.ArgumentParser(
        description="Compute TT gravitational-wave coupling from E_D_fields_FLASH_cavity CSV files."
    )
    parser.add_argument(
        "--folder",
        default="E_D_fields_FLASH_cavity",
        help="Folder containing input CSV files.",
    )
    parser.add_argument(
        "--mode-filter",
        default="",
        help="Only process CSV files containing this string. Empty string processes all CSV files.",
    )
    parser.add_argument("--z-min-mm", type=float, default=0.0, help="Minimum z value kept in the scan [mm].")
    parser.add_argument("--z-max-mm", type=float, default=1500.0, help="Maximum z value kept in the scan [mm].")
    parser.add_argument("--n-angles", type=int, default=201, help="Number of beta and phi grid points.")
    parser.add_argument("--n-processes", type=int, default=6, help="Number of multiprocessing workers.")
    parser.add_argument(
        "--freq-xlsx",
        default=DEFAULT_FREQ_XLSX,
        help="Spreadsheet of mode/tuner eigenfrequencies [MHz] used to look up each file's frequency.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    mode_filter = args.mode_filter if args.mode_filter else None
    run_scan(
        folder_base=args.folder,
        mode_filter=mode_filter,
        z_min_mm=args.z_min_mm,
        z_max_mm=args.z_max_mm,
        n_angles=args.n_angles,
        n_processes=args.n_processes,
        freq_xlsx=args.freq_xlsx,
    )
