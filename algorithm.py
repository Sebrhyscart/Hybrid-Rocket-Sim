import numpy as np
from CoolProp.CoolProp import PropsSI

class InputError(Exception):
    pass

def golden_search(f, a, b, tol=1e-6):
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
