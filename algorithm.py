import numpy as np
import warnings

np.seterr(invalid='raise')

N_A = 6.02214076e23
k_B = 1.380649e-23

class InputError(Exception):
    pass

def golden_search(f, a, b, tol=1e-8):
    # ChatGPT's golden-search algorithm
    golden_ratio = (np.sqrt(5) - 1) / 2
    c = b - golden_ratio * (b - a)
    d = a + golden_ratio * (b - a)
    
    while abs(b - a) > tol:
        if abs(f(c)) < abs(f(d)):
            b = d
        else:
            a = c
        c = b - golden_ratio * (b - a)
        d = a + golden_ratio * (b - a)

    root = (a + b) / 2
    return root

def inches_to_meters(length_inches):
    return 0.0254 * length_inches