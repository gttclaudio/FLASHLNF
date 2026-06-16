import numpy as np
from scipy import integrate
from .base import Cavity

class CylindricalCavity(Cavity):
    def __init__(self, R, L, rho_min=1e-30, **params):
        super().__init__(**params)
        self.R = R
        self.L = L
        self.rho_min = rho_min

    # ---------------- integration geometry ----------------
    @property
    def bounds(self):
        return [
            (self.rho_min, self.R),
            (0.0, 2 * np.pi),
            (0.0, self.L),
        ]

    def jacobian(self, Y):
        return Y[0]

    # ---------------- coordinate conversion ----------------
    def cart_to_native(self, X):
        x, y, z = X
        rho = np.sqrt(x*x + y*y)
        phi = np.arctan2(y, x)
        return np.array([rho, phi, z])
    
    def native_to_cart(self, Y):
        rho, phi, z = Y
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return np.array([x, y, z])

    def cart_vec_to_native(self, V, coords):
        rho, phi, z = coords
        c, s = np.cos(phi), np.sin(phi)
        return np.array([
            c*V[0] + s*V[1],
            -s*V[0] + c*V[1],
            V[2]
        ])
        
    # ---------------- center ----------------
    def center(self):
        return np.array([0.0, 0.0, self.L/2])

    # ---------------- domain ----------------
    def inside(self, Y):
        rho, _, z = Y
        return (rho <= self.R) and (0.0 <= z <= self.L)

    # ---------------- slice ----------------
    def slice_limits(self, k, n_phi=500, n_z=500):
        z = np.linspace(0, self.L, n_z)
        phi = np.linspace(0.0, 2*np.pi, n_phi)

        Phi, Z = np.meshgrid(phi, z)
        X = self.R * np.cos(Phi)
        Y = self.R * np.sin(Phi)

        points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        x_par = [np.dot(k, p) - k[2]*self.L/2 for p in points]
        return min(x_par), max(x_par)
    
    def perp_lim(self):
        return np.sqrt(self.L**2 + 4*self.R**2)/2

    # ---------------- volume ----------------
    def volume(self):
        return np.pi * self.R**2 * self.L

