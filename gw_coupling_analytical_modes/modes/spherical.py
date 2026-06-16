import numpy as np
from scipy.special import spherical_jn, lpmv
from .base import CavityMode
from scipy.constants import c as c_cnst
from scipy.optimize import brentq
from scipy.special import lpmv, spherical_jn

class SphericalMode(CavityMode):
    """
    Spherical cavity modes:
    TMa, TMb, TEa, TEb
    Indices: m,n,p
    """
    def __init__(self, indices, mode_name, cavity):
        self.m, self.n, self.p = indices
        super().__init__(indices, mode_name, cavity)
        self.root = self._find_root()
        
    # ---------------- mode index validation ----------------
    def _is_zero_mode(self):
        m, n, p = self.m, self.n, self.p

        if self.mode_name == "TMb":
            if m == 0:
                return True

        if self.mode_name == "TEb":
            if m == 0:
                return True

        return False
    
    def _validate(self):
        m, n, p = self.m, self.n, self.p

        if self.mode_name not in ("TEa", "TEb", "TMa", "TMb"):
            raise ValueError(f"Invalid mode_name {self.mode_name} for spherical cavity")

        if abs(m) > n:
            raise ValueError(f"Spherical modes require |m| =< n")
            
        if n <= 0 or p <= 0:
            raise ValueError(f"Spherical modes require n >= 1 and p >= 1")

    # ---------------- root finding ----------------
    def _find_root(self):
        if self.mode_name in ("TMa", "TMb"):
            func = jn_hat_der
        elif self.mode_name in ("TEa", "TEb"):
            func = spherical_jn
        return find_zero(self.n, self.p, func)

    # ---------------- wavevector ----------------
    def k(self):
        return self.root / self.cavity.R

    def omega(self):
        return c_cnst * self.k()

    # ---------------- prenormalized E field ----------------
    def E_prenorm(self, Y):   
        if self._is_zero_mode():
            return np.zeros(3, dtype=complex)
        
        k = self.k()
        m, n, p = self.m, self.n, self.p
        R = self.cavity.R
        root_np = self.root
        
        Y[0] = max(Y[0], 1e-12)
        Y[1] = max(Y[1], 1e-12)
        
        sin_theta_safe = np.sin(Y[1])
        
        if self.mode_name == 'TEa':

            Er = 0.0
            Etheta = m / (k * Y[0] * sin_theta_safe) * jn_hat(n, root_np * Y[0] / R) * lpmv(m, n, np.cos(Y[1])) * np.sin(m * Y[2])
            Ephi = 1 / (k * Y[0]) * jn_hat(n, root_np * Y[0] / R) * dPnm_dtheta(n, m, Y[1]) * np.cos(m * Y[2])

        elif self.mode_name == 'TEb':

            Er = 0.0
            Etheta = - m / (k * Y[0] * sin_theta_safe) * jn_hat(n, root_np * Y[0] / R) * lpmv(m, n, np.cos(Y[1])) * np.cos(m * Y[2])
            Ephi = 1 / (k * Y[0]) * jn_hat(n, root_np * Y[0] / R) * dPnm_dtheta(n, m, Y[1]) * np.sin(m * Y[2])
            
        elif self.mode_name == 'TMa':

            Er = - n * (n + 1) / (k * Y[0]**2) * jn_hat(n, root_np * Y[0] / R) * lpmv(m, n, np.cos(Y[1])) * np.cos(m * Y[2])
            Etheta = - 1 / (Y[0]) * jn_hat_der(n, root_np * Y[0] / R) * dPnm_dtheta(n, m, Y[1]) * np.cos(m * Y[2])
            Ephi = m / (Y[0] * sin_theta_safe) * jn_hat_der(n, root_np * Y[0] / R) * lpmv(m, n, np.cos(Y[1])) * np.sin(m * Y[2])

        elif self.mode_name == 'TMb':

            Er = - n * (n + 1) / (k * Y[0]**2) * jn_hat(n, root_np * Y[0] / R) * lpmv(m, n, np.cos(Y[1])) * np.sin(m * Y[2])
            Etheta = - 1 / (Y[0]) * jn_hat_der(n, root_np * Y[0] / R) * dPnm_dtheta(n, m, Y[1]) * np.sin(m * Y[2])
            Ephi = - m / (Y[0] * sin_theta_safe) * jn_hat_der(n, root_np * Y[0] / R) * lpmv(m, n, np.cos(Y[1])) * np.cos(m * Y[2])
            
        else:
            raise ValueError(f"Invalid mode {self.mode_name} {self.mode_ind}")

        return np.array([Er, Etheta, Ephi])
        
    # ---------------- normalized E field ----------------
    def E(self, Y):
        if self.norm is None:
            raise RuntimeError("Mode not normalized")
        return self.E_prenorm(Y) / np.sqrt(self.norm)
        
    # ---------------- normalization ----------------
    def normalize(self):
        if self._is_zero_mode():
            self.norm = 1
        else:
            def E1(Y): return self.E_prenorm(Y)
            self.norm = self.cavity.overlap_integral(E1, E1)

# ---------------- helper functions ----------------
def jn_hat(n, x):
    return x * spherical_jn(n, x)

def jn_hat_der(n, x):
    return spherical_jn(n, x) + x * spherical_jn(n, x, derivative=True)

def lpmv_ext(m, n, x):
    if abs(m) > n:
        return 0.0
    return lpmv(m, n, x)

def dPnm_dtheta(n, m, theta):
    #avoid division by zero at poles
    sin_t = np.sin(theta)
    sin_safe = np.where(np.abs(sin_t) > 1e-12, sin_t,
                        np.sign(sin_t + 1e-300) * 1e-12)
    return (n*np.cos(theta) * lpmv_ext(m, n, np.cos(theta))
          - (n+m) * lpmv_ext(m, n-1, np.cos(theta))) / sin_safe

def derivative_spherical_jn(n, x):
    if n==0:
        return -spherical_jn(1, x)
    return spherical_jn(n-1, x) - (n+1)/ x * spherical_jn(n, x)

def find_zero(n, p, function, x_max=100, num_points=10000):
    xs = np.linspace(0.01, x_max, num_points)
    ys = [function(n, x) for x in xs]
    zeros = []
    for i in range(len(xs)-1):
        if ys[i]*ys[i+1] < 0:
            zero = brentq(lambda x: function(n, x), xs[i], xs[i+1])
            zeros.append(zero)
            if len(zeros) > p:
                break
    if not zeros:
        raise ValueError(f"No zero found for n={n}, p={p}")
    return zeros[-1]