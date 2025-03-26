import numpy as np
import matplotlib.pyplot as plt

def pot(x,params):
    Q, e, a, x0, eps0 = params
    return  -((Q*e)/4*np.pi*eps0)/ np.sqrt((x-x0)**2 + a**2)


