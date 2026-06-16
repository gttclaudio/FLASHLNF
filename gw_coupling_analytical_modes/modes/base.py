import numpy as np
from abc import ABC, abstractmethod
from scipy.constants import c as c_cnst

class CavityMode(ABC):
    def __init__(self, indices, mode_name, cavity):
        """
        indices: tuple of mode numbers (m,n,p) or (n,p,q) depending on cavity
        mode_name: string like 'TE', 'TM', 'TMa', 'TMb', etc.
        cavity: instance of cavity class
        """
        self.indices = indices
        self.mode_name = mode_name
        self.cavity = cavity
        self.norm = None
        
        self._validate()
        
    def _is_zero_mode(self):
        """
        Check whether the mode is a zero mode.
        """
        pass

    @abstractmethod
    def _validate(self):
        """
        Validate indices and mode_name for this geometry.
        Must raise ValueError on invalid combinations.
        """
        pass

    @abstractmethod
    def k(self):
        """Return magnitude of wavevector."""
        pass

    def omega(self):
        """Angular frequency."""
        return c_cnst * self.k()

    @abstractmethod
    def E_prenorm(self, Y):
        """Prenormalized E field at point Y."""
        pass

    def E(self, Y):
        """Normalized E field."""
        if self.norm is None:
            raise RuntimeError("Mode not normalized")
        return self.E_prenorm(Y) / np.sqrt(self.norm)

    def normalize(self):
        """Compute normalization factor from cavity overlap integral."""
        def E1(Y): return self.E_prenorm(Y)
        self.norm = self.cavity.overlap_integral(E1, E1)
