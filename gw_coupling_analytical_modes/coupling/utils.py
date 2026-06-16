import numpy as np

c_cnst = 2.99792458e8

def h_monochromatic(amplitude, tau, omega, phase=0):
    return amplitude * np.exp(1j * (omega * tau + phase))

def mean_calc(eta, theta):
    eta_sum, sin_sum = 0.0, 0.0
    for row in eta:          
        for i, element in enumerate(row):
            sin_theta = np.sin(theta[i])
            eta_sum += element * sin_theta
            sin_sum += sin_theta
    
    result = eta_sum / sin_sum if sin_sum > 0 else 0.0

    return result

def compute_k_pol(theta, phi):

    k = -np.array([np.sin(theta)*np.cos(phi), 
                  np.sin(theta)*np.sin(phi), 
                  np.cos(theta)
                 ])

    e1 = np.array([-np.sin(phi), 
                   np.cos(phi), 
                   0.0
                 ])

    e2 = np.array([np.cos(theta)*np.cos(phi), 
                   np.cos(theta)*np.sin(phi), 
                   -np.sin(theta)
                 ])

    return k, e1, e2

def decompose_B(B, k, e1, e2):
    """
    Compute B_plus and B_cross polarizations perpendicular to k
    """
    Bperp = B - np.dot(B, k) * k
    B_plus  = np.dot(Bperp, e2) * e1 + np.dot(Bperp, e1) * e2
    B_cross = -np.dot(Bperp, e1) * e1 + np.dot(Bperp, e2) * e2
    
    return B_plus, B_cross

def make_jeff(B, cavity, hplus, hcross, k, e1, e2):

    B_plus, B_cross = decompose_B(B, k, e1, e2)
    center = cavity.center()

    def tau(Y, t):
        return t - np.vdot(k, cavity.native_to_cart(Y) - center) / c_cnst

    def jeff_from_B(B, h):
        def jeff(Y, t):
            return cavity.cart_vec_to_native(h(tau(Y, t)) * B, Y)
        return jeff

    jeff_plus  = jeff_from_B(B_plus,  hplus)
    jeff_cross = jeff_from_B(B_cross, hcross)

    def jeff_full(Y, t):
        return jeff_plus(Y, t) + jeff_cross(Y, t)

    return {
        "plus":  jeff_plus,
        "cross": jeff_cross,
        "mix":   jeff_full,
    }
