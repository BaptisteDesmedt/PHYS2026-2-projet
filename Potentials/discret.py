import numpy as np
import matplotlib.pyplot as plt


def discretiser(x_min, x_max, num_intervalles, pot, params):
    x_bins = np.linspace(x_min, x_max, num_intervalles + 1)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2

    V_bins = pot(x_centers, params)
    V_bins[0] = 0.0
    V_bins[-1] = 0.0
    x_centers[0] = x_min
    x_centers[-1] = x_max
    return x_bins, V_bins, x_centers

def discretiser_period(x_min, x_max, num_intervalles, pot,L,N, params):
    num_intervallesC = num_intervalles -2
    delta_x = (x_max - x_min) / (num_intervallesC + 1)

    x_bins = np.zeros(num_intervalles + 1)
    x_bins[0] = x_min
    x_bins[-1] = x_max
    x_bins[1] = x_min + delta_x / 2
    x_bins[-2] = x_max - delta_x / 2

    for i in range(2, num_intervallesC + 1):
        x_bins[i] = x_min + delta_x / 2 + (i - 1) * delta_x
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2
    V_bins = pot(x_centers, params, N, L)
    V_bins[0] = 0.0
    V_bins[-1] = 0.0
    x_centers[0] = x_min
    x_centers[-1] = x_max
    
    return x_bins, V_bins, x_centers