import numpy as np
from .base import CavityMode
from scipy.constants import c as c_cnst

class RectangularMode(CavityMode):
    """
    Rectangular cavity modes:
    TM, TE
    Indices: m,n,p
    """
    def __init__(self, indices, mode_name, cavity):
        self.m, self.n, self.p = indices
        super().__init__(indices, mode_name, cavity)
        
    # ---------------- mode index validation ----------------
    def _is_zero_mode(self):
        return False


    def _validate(self):
        m, n, p = self.m, self.n, self.p

        if self.mode_name not in ("TE", "TM"):
            raise ValueError(f"Invalid mode_name {self.mode_name} for rectangular cavity")

        if any(i < 0 for i in (m, n, p)):
            raise ValueError("Mode indices must be non-negative")
            
        if m == 0 or n == 0:
                raise ValueError("Rectangular TM modes require m >= 1 and n >= 1")

        if self.mode_name == "TE":
            if p == 0:
                raise ValueError("Rectangular TE modes require p >= 1")
            

    # ---------------- wavevector ----------------
    def k(self):
        a, b, c = self.cavity.a, self.cavity.b, self.cavity.c
        return np.sqrt((self.m*np.pi/a)**2 + (self.n*np.pi/b)**2 + (self.p*np.pi/c)**2)

    def omega(self):
        return c_cnst * self.k()

    # ---------------- prenormalized E field ----------------
    def E_prenorm(self, Y):

        if self._is_zero_mode():
            return np.zeros(3, dtype=complex)
        
        k = self.k()
        a, b, c = self.cavity.a, self.cavity.b, self.cavity.c
        m, n, p = self.m, self.n, self.p
        
        if self.mode_name == 'TE':

            Ex = - 1j / (k**2 - (p * np.pi / c)**2 ) * (n * np.pi / b) * np.cos(m * np.pi / a * Y[0]) * np.sin(n * np.pi / b * Y[1]) * np.sin(p * np.pi / c * Y[2])
            Ey = 1j / (k**2 - (p * np.pi / c)**2 )  * (m * np.pi / a) * np.sin(m * np.pi / a * Y[0]) * np.cos(n * np.pi / b * Y[1]) * np.sin(p * np.pi / c * Y[2])
            Ez =  0.0

        elif self.mode_name == 'TM':

            Ex = -1 / (k**2 - (p * np.pi / c)**2 ) * (m * np.pi / a) * (p * np.pi / c)  * np.cos(m * np.pi / a * Y[0]) * np.sin(n * np.pi / b * Y[1]) * np.sin(p * np.pi / c * Y[2])
            Ey = 1 / (k**2 - (p * np.pi / c)**2 ) * (n * np.pi / b) * (p * np.pi / c)  * np.sin(m * np.pi / a * Y[0]) * np.cos(n * np.pi / b * Y[1]) * np.sin(p * np.pi / c * Y[2])
            Ez = np.sin(m * np.pi / a * Y[0]) * np.sin(n * np.pi / b * Y[1]) * np.cos(p * np.pi / c * Y[2])

        return np.array([Ex, Ey, Ez])

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
