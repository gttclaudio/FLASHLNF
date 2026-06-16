import numpy as np
from tqdm import tqdm
from coupling.utils import h_monochromatic, compute_k_pol, make_jeff
from multiprocessing import Pool

def compute_coupling(args):
    cavity, mode, B, pol, omega, k, e1, e2, t = args
    def hplus(tau):  return h_monochromatic(amplitude=1.0, tau=tau, omega=omega)
    def hcross(tau): return h_monochromatic(amplitude=1.0, tau=tau, omega=omega, phase=np.pi/2)

    jeff = make_jeff(B=B, cavity=cavity, hplus=hplus, hcross=hcross, k=k, e1=e1, e2=e2)[pol]
    V = cavity.volume()
    def E1(Y): return mode.E(Y)
    def E2(Y): return jeff(Y, t)

    coupling = cavity.overlap_integral(E1, E2, method="nquad", epsabs=1e-8, epsrel=1e-6, limit=80, complex_value=True) / np.sqrt(V)

    return coupling


class CouplingStrength:
    def __init__(self, cavity, mode, theta_vals, phi_vals: int = 0, B=(0.0, 0.0, 1.0), pol: str = "cross", nproc: int = 1):
        self.cavity = cavity
        self.mode = mode
        self.B = np.asarray(B, dtype=float)
        self.pol = str(pol)
        self.nproc = int(nproc)
        self.theta_vals = theta_vals
        self.phi_vals = phi_vals

        self.omega = mode.omega()
        
    def run(self, t=0.0, k=None):
        directions = []
        for phi in self.phi_vals:
            for theta in self.theta_vals:
                k, e1, e2 = compute_k_pol(theta, phi)
                directions.append((k, e1, e2))
        
        eta_k = []
        mode_numbers = ''.join(map(str, self.mode.indices))
        desc = f"[INFO] Computing coupling for mode {self.mode.mode_name}_{mode_numbers} with {self.pol} polarization"

        args = [(self.cavity, self.mode, self.B, self.pol, self.omega, k, e1, e2, t)
            for (k, e1, e2) in directions]

        with Pool(processes=self.nproc) as pool:
            it = pool.imap(compute_coupling, args)
            eta_k = list(tqdm(it, total=len(args), desc=desc, unit="k"))

        eta_k = np.array(eta_k)
        eta_reshaped = eta_k.reshape(len(self.phi_vals), len(self.theta_vals))
        return eta_reshaped