import numpy as np
from tqdm import tqdm
from coupling.utils import compute_k_pol
from coupling.sources import gw_coupling, axion_coupling, dp_coupling, scalar_coupling
from multiprocessing import Pool

SOURCE_REGISTRY = {
    "gw": gw_coupling,
    "axion": axion_coupling,
    "dp": dp_coupling,
    "scalar": scalar_coupling,
}

class CouplingStrength:
    def __init__(self, cavity, mode, source, beta_vals, phi_vals, B=(0.0, 0.0, 1.0), pol: str = "cross", nproc: int = 1):
        self.cavity = cavity
        self.mode = mode

        self.source = source
        self.kernel = SOURCE_REGISTRY[source]

        self.B = np.asarray(B, dtype=float)
        self.pol = str(pol)
        self.nproc = int(nproc)

        self.beta_vals = beta_vals
        self.phi_vals = phi_vals

        self.omega = mode.omega()

    def requires_direction_scan(self):

        return self.source in {"gw", "dp"}

    def _run_directional(self):

        directions = []

        for phi in self.phi_vals:
            for beta in self.beta_vals:

                k, e1, e2 = compute_k_pol(beta, phi)

                directions.append((k, e1, e2))

        args = []

        if self.source == "gw":

            args = [
                (self.cavity, self.mode, self.B, self.pol, self.omega, k, e1, e2)
                for k, e1, e2 in directions
            ]

        elif self.source == "dp":

            args = [
                (self.cavity, self.mode, k)
                for k, _, _ in directions
            ]

        with Pool(self.nproc) as pool:

            C = list(tqdm(pool.imap(self.kernel, args), total=len(args)))

        return np.array(C).reshape(len(self.phi_vals), len(self.beta_vals))

    def _run_single(self):

        if self.source == "axion":

            args = (self.cavity, self.mode,self.B)

        elif self.source == "scalar":

            args = (self.cavity, self.mode, self.B,)

        return self.kernel(args)
    
    def run(self):
    
        if self.requires_direction_scan():
            return self._run_directional()

        return self._run_single()