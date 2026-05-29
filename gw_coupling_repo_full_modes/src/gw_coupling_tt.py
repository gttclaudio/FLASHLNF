"""
Compute transverse-traceless (TT) gravitational-wave coupling from cavity-field CSV data.

This module reads electromagnetic field maps exported from COMSOL, assigns finite-volume
weights to sampled grid points, and evaluates two TT-coupling terms over a scan of
incident angles. Each input mode file is saved as a pickle file containing a DataFrame
with the scan angles and coupling values.

Main modules
------------
numpy
    Numerical array operations, complex arithmetic, trigonometric functions, and sums.
pandas
    CSV input, table cleanup, grouping repeated grid points, and saving results.
pathlib
    Safer, platform-independent path handling for input and output files.
itertools
    Cartesian product of the beta and phi angle grids.
multiprocessing
    Parallel evaluation of independent angle points.
tqdm
    Progress bar for the scan.
argparse
    Command-line options for selecting a mode, CSV files, angle grid, z range,
    and worker processes.

Expected CSV format
-------------------
The CSV file is expected to contain eight columns after skipping the first eight header
rows:
    x, y, z, Ex, Ey, Ez, normE, normD
Coordinates are assumed to be in mm. Field values may be real or complex strings using
"i" as the imaginary unit.
"""

import argparse
import itertools
import multiprocessing as mp
import re
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


EPSILON_0 = 8.854e-12
SPEED_OF_LIGHT = 2.998e8

MODE_FILES = {
    "TE011": "FLASH_LF_TE011_218.3MHz.csv",
    "TM010": "FLASH_LF_TM010_129.0MHz.csv",
}


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

    Some (x, y) positions may have only one sampled z point. In that case, the nearest
    global z-node weight is used as a fallback.
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
    """Return the TT phase factor for propagation coordinate z_m and wave number k."""
    return np.exp(1j * z_m * wave_number)


def _convert_complex_columns(df):
    """Convert COMSOL-style complex strings to Python complex values."""
    for column in ["Ex", "Ey", "Ez", "normE", "normD"]:
        if df[column].dtype == object:
            df[column] = df[column].apply(lambda value: complex(str(value).replace("i", "j")))
    return df


def compute_tt_coupling_from_data(df_raw, beta, freq_hz, phi=0.0):
    """
    Compute the two TT-coupling contributions for one pair of scan angles.

    Parameters
    ----------
    df_raw : pandas.DataFrame
        Field-map data with coordinates and complex electric/displacement fields.
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
    df = _convert_complex_columns(df_raw.copy())
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
    cross_integral = np.sum(volume_weights * (ey_phi * j_eff_y + ez * j_eff_z))

    denominator = cavity_volume * normalization_integral
    coupling_parallel = np.abs(parallel_integral) ** 2 / denominator
    coupling_cross = np.abs(cross_integral) ** 2 / denominator

    return coupling_parallel, coupling_cross


def _worker_angle(args):
    """Multiprocessing wrapper for one beta/phi grid point."""
    beta, df, freq_hz, phi = args
    return compute_tt_coupling_from_data(df, beta, freq_hz, phi)


def load_field_csv(csv_path):
    """Load one field-map CSV file and assign the expected column names."""
    df = pd.read_csv(csv_path, skiprows=8)
    df.columns = ["x", "y", "z", "Ex", "Ey", "Ez", "normE", "normD"]
    return df.dropna()


def parse_mode_and_frequency(csv_path):
    """
    Infer the cavity mode name and frequency from a filename.

    Expected filename examples:
        FLASH_LF_TE011_218.3MHz.csv
        FLASH_LF_TE011_218.3MHz_sample.csv
    """
    name = Path(csv_path).name
    mode_match = re.search(r"_(TE\d+|TM\d+)_", name)
    freq_match = re.search(r"_([0-9]+(?:\.[0-9]+)?)MHz", name)

    if mode_match is None or freq_match is None:
        raise ValueError(
            f"Could not infer mode/frequency from filename: {name}. "
            "Expected a name like FLASH_LF_TE011_218.3MHz.csv."
        )

    mode = mode_match.group(1)
    freq_mhz = float(freq_match.group(1))
    return mode, freq_mhz * 1e6


def run_scan_for_csv(csv_path, z_min_mm=0.0, z_max_mm=1500.0, n_angles=201, n_processes=6, output_dir=None):
    """Run the TT-coupling scan for one CSV file."""
    csv_path = Path(csv_path)
    output_dir = Path(output_dir) if output_dir is not None else csv_path.parent

    mode, freq_hz = parse_mode_and_frequency(csv_path)
    print(f"Input: {csv_path}")
    print(f"Mode: {mode}, frequency: {freq_hz / 1e6:.3f} MHz")

    df = load_field_csv(csv_path)
    print("z range before cut [mm]:", np.max(df["z"]) - np.min(df["z"]))

    df = df[(df["z"] >= z_min_mm) & (df["z"] <= z_max_mm)]
    print(
        "x unique:", len(df["x"].unique()),
        "y unique:", len(df["y"].unique()),
        "z unique:", len(df["z"].unique()),
    )

    single_z_groups = sum(1 for _, group in df.groupby(["x", "y"]) if group["z"].nunique() == 1)
    total_groups = df.groupby(["x", "y"]).ngroups
    print(f"(x,y) groups: {total_groups}, single-z groups: {single_z_groups}")

    phis = np.linspace(0.0, 2.0 * np.pi, n_angles)
    betas = np.linspace(0.0, 2.0 * np.pi, n_angles)
    angle_pairs = list(itertools.product(betas, phis))
    tasks = [(beta, df, freq_hz, phi) for beta, phi in angle_pairs]

    if n_processes == 1:
        results = [_worker_angle(task) for task in tqdm(tasks, total=len(tasks))]
    else:
        with mp.Pool(processes=n_processes) as pool:
            results = list(tqdm(pool.imap(_worker_angle, tasks), total=len(tasks)))

    records = [
        {"beta": beta, "phi": phi, "coupling_parallel": coupling_p, "coupling_cross": coupling_c}
        for (beta, phi), (coupling_p, coupling_c) in zip(angle_pairs, results)
    ]
    df_results = pd.DataFrame(records)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"TT_gauge_{mode}_{freq_hz / 1e6:.4f}MHz.pkl"
    df_results.to_pickle(output_path)
    print(f"Saved: {output_path}")
    return output_path


def discover_csv_files(folder, filename_filter=""):
    """Return CSV files in a folder, optionally restricted by a substring filter."""
    folder = Path(folder)
    return sorted(
        path for path in folder.glob("*.csv")
        if "MHz" in path.name and (not filename_filter or filename_filter in path.name)
    )



def resolve_mode_csv(mode, data_dir="data"):
    """Return the CSV path corresponding to a named mode stored in data_dir."""
    if mode not in MODE_FILES:
        valid_modes = ", ".join(MODE_FILES)
        raise ValueError(f"Unknown mode {mode!r}. Available modes: {valid_modes}")
    return Path(data_dir) / MODE_FILES[mode]

def run_scan(csv_files=None, folder=None, filename_filter="", **kwargs):
    """
    Run the TT-coupling scan for explicit CSV files or for matching files in a folder.

    Use csv_files for clear examples and production jobs with known inputs.
    Use folder plus filename_filter for batch processing many modes in one directory.
    """
    if csv_files is None:
        if folder is None:
            raise ValueError("Provide either csv_files or folder.")
        csv_files = discover_csv_files(folder, filename_filter=filename_filter)

    if not csv_files:
        raise FileNotFoundError("No matching CSV files were found.")

    outputs = []
    for csv_file in csv_files:
        outputs.append(run_scan_for_csv(csv_file, **kwargs))
    return outputs


def parse_args():
    """Parse command-line options for mode-based or explicit CSV processing."""
    parser = argparse.ArgumentParser(
        description="Compute TT gravitational-wave coupling from cavity-field CSV files."
    )
    parser.add_argument(
        "--mode",
        choices=sorted(MODE_FILES),
        help="Mode to process using the built-in data/MODE filename map, e.g. TE011 or TM010.",
    )
    parser.add_argument(
        "--csv",
        nargs="+",
        help="Explicit CSV file(s) to process. Overrides --mode if provided.",
    )
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Folder containing the full CSV files used by --mode. Default: data.",
    )
    parser.add_argument(
        "--folder",
        help="Optional folder scan for advanced batch use. Used only when neither --csv nor --mode is provided.",
    )
    parser.add_argument(
        "--filename-filter",
        default="",
        help="Optional substring filter for --folder mode. Not needed for --mode or --csv.",
    )
    parser.add_argument("--z-min-mm", type=float, default=0.0, help="Minimum z value kept in the scan [mm].")
    parser.add_argument("--z-max-mm", type=float, default=1500.0, help="Maximum z value kept in the scan [mm].")
    parser.add_argument("--n-angles", type=int, default=201, help="Number of beta and phi grid points.")
    parser.add_argument("--n-processes", type=int, default=6, help="Number of multiprocessing workers.")
    parser.add_argument("--output-dir", default="results", help="Folder for output pickle files. Default: results.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    common_kwargs = dict(
        z_min_mm=args.z_min_mm,
        z_max_mm=args.z_max_mm,
        n_angles=args.n_angles,
        n_processes=args.n_processes,
        output_dir=args.output_dir,
    )

    if args.csv:
        run_scan(csv_files=args.csv, **common_kwargs)
    elif args.mode:
        run_scan(csv_files=[resolve_mode_csv(args.mode, args.data_dir)], **common_kwargs)
    else:
        run_scan(
            folder=args.folder or args.data_dir,
            filename_filter=args.filename_filter,
            **common_kwargs,
        )
