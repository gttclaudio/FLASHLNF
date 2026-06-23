import numpy as np
from coupling.utils import h_monochromatic, make_jeff

c_cnst = 3e8

def gw_coupling(args):
    cavity, mode, B_field, pol, omega, k, e1, e2 = args
    def hplus(tau):  return h_monochromatic(amplitude=1.0, tau=tau, omega=omega)
    def hcross(tau): return h_monochromatic(amplitude=1.0, tau=tau, omega=omega, phase=np.pi/2)

    jeff = make_jeff(B=B_field, cavity=cavity, hplus=hplus, hcross=hcross, k=k, e1=e1, e2=e2)[pol]
    V = cavity.volume()
    def E1(Y): return mode.E(Y)
    def E2(Y): return jeff(Y, t=0.0)

    coupling = np.abs(cavity.overlap_integral(E1, E2,
                                              method="nquad", epsabs=1e-8, epsrel=1e-6,
                                              limit=80, complex_value=True))**2 / V

    return coupling

def axion_coupling(args):
    cavity, mode, B_field = args

    V = cavity.volume()
    def E1(Y): return mode.E(Y)
    def E2(Y): return cavity.cart_vec_to_native(B_field, Y)

    coupling = np.abs(cavity.overlap_integral(E1, E2, 
                                       method="nquad", epsabs=1e-8, 
                                       epsrel=1e-6, limit=80, complex_value=True))**2 / V
    


    return coupling

def dp_coupling(args):
    cavity, mode, k = args
  
    V = cavity.volume()
    def E1(Y): return mode.E(Y)
    def E2(Y): return cavity.cart_vec_to_native(k, Y)

    coupling = np.abs(cavity.overlap_integral(E1, E2,
                                            method="nquad", epsabs=1e-8, 
                                            epsrel=1e-6, limit=80, complex_value=True))**2 / V
    


    return coupling

def scalar_coupling(args):
    cavity, mode, B_field, k, k_scale = args

    k = k * k_scale
  
    V = cavity.volume()
    def E1(Y): return mode.B(Y)
    def E2(Y): return cavity.cart_vec_to_native(B_field, Y) * np.exp(1j * (np.dot(k, cavity.native_to_cart(Y)))) # And possibly a phase term 

    coupling = np.abs(cavity.overlap_integral(E1, E2,
                                              method="nquad", epsabs=1e-8, epsrel=1e-6, 
                                              limit=80, complex_value=True))**2 / V


    return coupling