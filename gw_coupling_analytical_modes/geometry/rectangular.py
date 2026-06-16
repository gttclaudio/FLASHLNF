import numpy as np
from .base import Cavity
from scipy import integrate

class RectangularCavity(Cavity):
    def __init__(self, a, b, c, **params):
        super().__init__(**params)
        self.a = a
        self.b = b
        self.c = c

    # ---------------- integration geometry ----------------
    @property
    def bounds(self):
        return [
            (0.0, self.a),
            (0.0, self.b),
            (0.0, self.c),
        ]

    def jacobian(self, Y):
        return 1.0
    
    # ---------------- coordinate conversion ----------------
    def cart_to_native(self, X):
        # For a rectangle, Cartesian coords are already native
        return np.array(X)
    
    def native_to_cart(self, Y):
        # For a rectangle, Cartesian coords are already native
        return np.array(Y)

    def cart_vec_to_native(self, V, *args):
        # Cartesian basis is the natural basis
        return np.array(V)
        
    # ---------------- center ----------------
    def center(self):
        return np.array([self.a/2, self.b/2, self.c/2])

    # ---------------- domain ----------------
    def inside(self, Y):
        x, y, z = Y
        return (0 <= x <= self.a) and (0 <= y <= self.b) and (0 <= z <= self.c)

    # ---------------- slice ----------------
    def slice_limits(self, k):
        # Compute min/max along k using the 8 vertices
        p0 = np.array([self.a/2, self.b/2, self.c/2])
        vertices = np.array([[x,y,z] for x,y,z in
                             [(0,0,0),(0,0,self.c),(0,self.b,0),(0,self.b,self.c),
                              (self.a,0,0),(self.a,0,self.c),(self.a,self.b,0),(self.a,self.b,self.c)]])
        x_par = [np.dot(k,v)-np.dot(k,p0) for v in vertices]
        return min(x_par), max(x_par)

    def perp_lim(self):
        return np.sqrt(self.a**2 + self.b**2 + self.c**2)/2
        
    # ---------------- volume ----------------
    def volume(self):
        return self.a * self.b * self.c