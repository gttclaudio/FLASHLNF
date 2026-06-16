import numpy as np
from scipy.special import jv, jvp, jn_zeros, jnp_zeros
from .base import CavityMode
from scipy.constants import c as c_cnst

class CylindricalMode(CavityMode):
    """
    Cylindrical cavity modes:
    TMa, TMb, TEa, TEb
    Indices: n,p,q
    """
    def __init__(self, indices, mode_name, cavity):
        self.n, self.p, self.q = indices
        super().__init__(indices, mode_name, cavity)
        self.root = self._find_root()
        
    # ---------------- mode index validation ----------------
    def _is_zero_mode(self):
        n, p, q = self.n, self.p, self.q
    
        if self.mode_name == "TEa":
            if n == 0 and q == 0:
                return True
    
        if self.mode_name == "TEb":
            if q == 0:
                return True
    
        if self.mode_name == "TMa":
            if n == 0:
                return True
    
        return False

    def _validate(self):
        n, p, q = self.n, self.p, self.q

        if self.mode_name not in ("TMa", "TMb", "TEa", "TEb"):
            raise ValueError(f"Invalid mode_name {self.mode_name} for cylindrical cavity")
        
        if any(i < 0 for i in (n, p, q)):
            raise ValueError("Mode indices must be non-negative")
            
        if p == 0:
            raise ValueError("Cylindrical modes require p >= 1")
            
        if self.mode_name.startswith("TE") and q == 0:
            raise ValueError("TE cylindrical modes require q >= 1")
            
    def _find_root(self):
        if self.mode_name in ("TMa", "TMb"):
            return jn_zeros(self.n, self.p)[self.p-1]
        elif self.mode_name in ("TEa", "TEb"):
            return jnp_zeros(self.n, self.p)[self.p-1]
            
    # ---------------- wavevector ----------------
    def k(self):
        L, R = self.cavity.L, self.cavity.R
        return np.sqrt((self.q * np.pi/L)**2 + (self.root/R)**2)
        
    def omega(self):
        return c_cnst * self.k()
        
    # ---------------- prenormalized E field ----------------
    def E_prenorm(self, Y):

        if self._is_zero_mode():
            return np.zeros(3, dtype=complex)
            
        k = self.k()
        n, p, q = self.n, self.p, self.q
        L, R = self.cavity.L, self.cavity.R
        root_np = self.root

        Y = np.array(Y, copy=True)  # defensive copy
        Y[0] = max(Y[0], 1e-30)  # avoid division by zero at rho=0
        
        if self.mode_name == 'TEa':
                
            Er   = 1 / (k**2 - (q * np.pi / L)**2 ) * (n / Y[0]) * jv(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.cos(n * Y[1])
            Ephi = - 1 / (k**2 - (q * np.pi / L)**2 ) * (root_np / R) * jvp(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.sin(n * Y[1])
            Ez   = 0.0
            
        elif self.mode_name == 'TEb':
                
            Er   = - 1 / (k**2 - (q * np.pi / L)**2 ) * (n / Y[0]) * jv(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.sin(n * Y[1])
            Ephi = - 1 / (k**2 - (q * np.pi / L)**2 ) * (root_np / R) * jvp(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.cos(n * Y[1])
            Ez   = 0.0
            
        elif self.mode_name == 'TMa':
                
            Er   = - 1 / (k**2 - (q * np.pi / L)**2 ) * (q * np.pi / L) * (root_np / R) * jvp(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.sin(n * Y[1])
            Ephi = - 1 / (k**2 - (q * np.pi / L)**2 ) *  (q * np.pi / L) * (n / Y[0]) * jv(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.cos(n * Y[1])
            Ez   =  jv(n, root_np * Y[0] / R) * np.cos(q * np.pi / L * Y[2]) * np.sin(n * Y[1])
            
        elif self.mode_name == 'TMb':
                
            Er   = - 1 / (k**2 - (q * np.pi / L)**2 ) * (q * np.pi / L) * (root_np / R) * jvp(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.cos(n * Y[1])
            Ephi =  1 / (k**2 - (q * np.pi / L)**2 ) * (q * np.pi / L) * (n / Y[0]) * jv(n, root_np * Y[0] / R) * np.sin(q * np.pi / L * Y[2]) * np.sin(n * Y[1])
            Ez   =  jv(n, root_np * Y[0] / R) * np.cos(q * np.pi / L * Y[2]) * np.cos(n * Y[1])
            
        return np.array([Er, Ephi, Ez])
        
    # ---------------- normalized E field ----------------
    def E(self, Y):
        if self.norm is None:
            raise RuntimeError("Mode not normalized")
            
        return self.E_prenorm(Y)/np.sqrt(self.norm)
        
    # ---------------- normalization ----------------
    def normalize(self):
        if self._is_zero_mode():
            self.norm = 1
        else:    
            def E1(Y): return self.E_prenorm(Y)
            self.norm = self.cavity.overlap_integral(E1, E1)
