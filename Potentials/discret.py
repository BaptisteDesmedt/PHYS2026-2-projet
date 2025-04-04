import numpy as np
import matplotlib.pyplot as plt


def discretiser(x_min, x_max, num_intervalles, pot, params):
    x_bins = np.linspace(x_min, x_max, num_intervalles + 1)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2

    V_bins = pot(x_centers, params)
    V_bins[0] = 0.0
    V_bins[-1] = 0.0
    return x_bins, V_bins, x_centers

def discretiser_period(x_min, x_max, num_intervalles, pot,L,N, params):
    x_bins = np.linspace(x_min, x_max, num_intervalles + 1)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2

    V_bins = pot(x_centers,params,N,L)
    V_bins[0] = 0.0
    V_bins[-1] = 0.0
    
    return x_bins, V_bins, x_centers