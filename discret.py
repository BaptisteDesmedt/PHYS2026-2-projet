import numpy as np
import matplotlib.pyplot as plt

def discretiser(x_min, x_max, num_intervalles, pot, params):
    x_bins = np.linspace(x_min, x_max, num_intervalles + 1)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2

    V_bins = pot(x_centers, params)
    V_steps = np.repeat(V_bins, 2)
    x_steps = np.repeat(x_bins, 2)[1:-1]
    
    return x_steps, V_steps, x_bins, V_bins