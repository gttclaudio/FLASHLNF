import argparse, os
import numpy as np
import pandas as pd

from geometry import CylindricalCavity, SphericalCavity, RectangularCavity
from modes import CylindricalMode, SphericalMode, RectangularMode

from coupling.coupling import CouplingStrength
from coupling.utils import mean_calc


def parse_args():
    parser = argparse.ArgumentParser(description="Compute cavity coupling strength C(β, φ) to various source terms, geometries and modes.")
    
    # Geometry selection
    parser.add_argument("--geometry", choices=["rectangular", "cylindrical", "spherical"], default="cylindrical", help="Cavity geometry type")
    
    # Source term selection
    parser.add_argument("--source", choices=["gw", "axion", "scalar", "dp"], default="gw", help="Source term type")
    parser.add_argument("--wavenumber", type=float, default=1, help="Particle wavenumber in units of [m^-1]")

    # Mode selection
    parser.add_argument("--mode", default="TM010", help="Mode name in format TM010.")
    
    # Cylindrical/Spherical cavity
    parser.add_argument("--R", type=float, default=0.05, help="Cavity radius [m]")
    parser.add_argument("--L", type=float, default=0.05,  help="Cylindrical cavity height [m]")

    # Rectangular cavity, it is assumed the magnetic field is in the z-direction
    parser.add_argument("--a", type=float, default=None, help="Rectangular cavity x-dimension length [m]")
    parser.add_argument("--b", type=float, default=None, help="Rectangular cavity y-dimension length [m]")
    parser.add_argument("--c", type=float, default=None, help="Rectangular cavity z-dimension length [m]")
    
    # Simulation parameters
    parser.add_argument("--n-beta", type=int, default=11, help="Number of beta angles to sample")
    parser.add_argument("--n-phi", type=int, default=1, help="Number of phi angles to sample")
    parser.add_argument("--n-processes", type=int, default=1, help="Number of processors for parallel execution")
    
    # Results directory
    parser.add_argument("--output-dir", type=str, default="results")
    
    return parser.parse_args()

def make_filename(args, freq_mhz):
    if args.source == "gw":
        return f"TT_gauge_{args.mode}_{freq_mhz:.4f}MHz.pkl"
    if args.source == "dp":
        return f"DP_{args.mode}_{freq_mhz:.4f}MHz.pkl"
    if args.source == "axion":
        return f"axion_{args.mode}_{freq_mhz:.4f}MHz.pkl"
    if args.source == "scalar":
        return f"scalar_{args.mode}_{freq_mhz:.4f}MHz.pkl"
    
def compute_mode_sum(cavity, mode_class, mode_names, mode_ind, source, k_scale=1.0, beta_vals=None, phi_vals=None, pol=None, nproc=1):
    results = []

    for mode_name in mode_names:

        mode = mode_class(indices=mode_ind, mode_name=mode_name, cavity=cavity)

        mode.normalize()

        solver = CouplingStrength(cavity=cavity, mode=mode, source=source, 
                                  beta_vals=beta_vals, phi_vals=phi_vals, k_scale=k_scale,
                                  pol=pol, nproc=nproc)

        results.append(solver.run())

    C_total = results[0]
    for r in results[1:]:
        C_total = C_total + r

    return C_total


def main():
    args = parse_args()

    mode_fam = args.mode[:2]
    mode_ind = [int(i) for i in args.mode[2:]]

    # --- Create cavity and mode ---
    if args.geometry == "cylindrical":
        cavity = CylindricalCavity(R=args.R, L=args.L)
        mode_class = CylindricalMode
        mode_name_arr = [mode_fam + "a", mode_fam + "b"]
        
    elif args.geometry ==  "spherical":
        cavity = SphericalCavity(R=args.R)
        mode_class = SphericalMode
        mode_name_arr = [mode_fam + "a", mode_fam + "b"]
        
    elif args.geometry ==  "rectangular":
        cavity = RectangularCavity(a=args.a, b=args.b, c=args.c)
        mode_class = RectangularMode
        mode_name_arr = [mode_fam]

    mode_name = mode_name_arr[0]
    mode = mode_class(indices=mode_ind, mode_name=mode_name, cavity=cavity)
    freq_mhz = mode.omega() / (2 * np.pi * 1e6)
    print(f"[INFO] Mode {args.mode} frequency f = {freq_mhz:.4f} MHz.")

    beta_vals = np.linspace(0.0, np.pi, args.n_beta)
    phi_vals = np.linspace(0.0, 2.0 * np.pi, args.n_phi) 

    if args.source == "gw":

        C = {"plus": None, "cross": None}

        for pol in ["plus", "cross"]:
            C[pol] = compute_mode_sum(cavity=cavity, mode_class=mode_class, mode_names=mode_name_arr, mode_ind=mode_ind, source=args.source, 
                                      beta_vals=beta_vals, phi_vals=phi_vals, pol=pol, nproc=args.n_processes)

            mean_C = mean_calc(C[pol], beta_vals)
            max_C = np.max(C[pol])

            print(f"[INFO] Results for coupling strength in the {pol} polarisation:")
            print(f"⟨C(β, φ)⟩ = {mean_C:.4f}, Cₘₐₓ = {max_C:.4f}")

        # Build dataframe
        records = []

        for i_phi, phi in enumerate(phi_vals):
            for i_beta, beta in enumerate(beta_vals):

                records.append({
                    "beta": beta,
                    "phi": phi,
                    "coupling_parallel": C["plus"][i_phi, i_beta],
                    "coupling_cross": C["cross"][i_phi, i_beta],
                })

    elif args.source == "dp":

        C = compute_mode_sum(cavity=cavity, mode_class=mode_class, mode_names=mode_name_arr, mode_ind=mode_ind, source=args.source, 
                             beta_vals=beta_vals, phi_vals=phi_vals, pol=None, nproc=args.n_processes)
        
        mean_C = mean_calc(C, beta_vals)
        max_C = np.max(C)

        print(f"[INFO] Results for coupling strength to dark photon:")
        print(f"⟨C(β, φ)⟩ = {mean_C:.4f}, Cₘₐₓ = {max_C:.4f}")

                # Build dataframe
        records = []

        for i_phi, phi in enumerate(phi_vals):
            for i_beta, beta in enumerate(beta_vals):

                records.append({
                    "beta": beta,
                    "phi": phi,
                    "coupling": C[i_phi, i_beta],
                })
    
    elif args.source == "axion":

        C = compute_mode_sum(cavity=cavity, mode_class=mode_class, mode_names=mode_name_arr, mode_ind=mode_ind, source=args.source, 
                             beta_vals=None, phi_vals=None, pol=None, nproc=args.n_processes)

        print(f"[INFO] Results for coupling strength to axion:")
        print(f"C = {C:.4f}")

        records = []
        records.append({"coupling": C})


    elif args.source == "scalar":

        C = compute_mode_sum(cavity=cavity, mode_class=mode_class, mode_names=mode_name_arr, mode_ind=mode_ind, source=args.source, 
                             beta_vals=beta_vals, phi_vals=phi_vals, k_scale=args.wavenumber, pol=None, nproc=args.n_processes)
        
        mean_C = mean_calc(C, beta_vals)
        max_C = np.max(C)

        print(f"[INFO] Results for coupling strength to scalar:")
        print(f"⟨C(β, φ)⟩ = {mean_C:.4f}, Cₘₐₓ = {max_C:.4f}")

                # Build dataframe
        records = []

        for i_phi, phi in enumerate(phi_vals):
            for i_beta, beta in enumerate(beta_vals):

                records.append({
                    "beta": beta,
                    "phi": phi,
                    "coupling": C[i_phi, i_beta],
                })

    filename = make_filename(args, freq_mhz)

    df = pd.DataFrame(records)

    df.attrs = {
        "geometry": args.geometry, "mode": args.mode,
        "frequency_mhz": freq_mhz,
    }

    if args.source in ["gw", "dp"]:
        df.attrs.update({
            "n_beta": args.n_beta, "n_phi": args.n_phi,
        })

    if args.geometry == "rectangular":
        df.attrs.update({
            "a": args.a, "b": args.b, "c": args.c,
        })

    if args.geometry == "cylindrical":
        df.attrs.update({
            "R": args.R, "L": args.L,
        })
        
    if args.geometry == "spherical":
        df.attrs.update({
            "R": args.R,
        })

    save_dir = os.path.join(args.output_dir, args.geometry)
    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, filename)
    df.to_pickle(filepath)

    print(f"[INFO] Results saved to {filepath}")

if __name__ == "__main__":
    main()
