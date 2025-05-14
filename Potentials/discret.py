import numpy as np

def discretiser(x_min, x_max, num_intervalles, pot, params):
    """
    Discretize a potential function over a given interval.

    Parameters:
    -----------
    x_min : float
        Minimum x value of the interval
    x_max : float
        Maximum x value of the interval
    num_intervalles : int
        Number of intervals for discretization
    pot : callable
        Potential function to discretize
    params : tuple
        Parameters for the potential function

    Returns:
    --------
    tuple
        x_bins : ndarray
            Boundaries of the intervals
        V_bins : ndarray
            Potential values at the center of each interval
        x_centers : ndarray
            Centers of the intervals
    """
    x_bins = np.linspace(x_min, x_max, num_intervalles + 1)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2

    V_bins = pot(x_centers, params)
    V_bins[0] = 0.0
    V_bins[-1] = 0.0
    x_centers[0] = x_min
    x_centers[-1] = x_max
    return x_bins, V_bins, x_centers

def discretiser_period(x_min, x_max, num_intervalles, pot, L, N, params):
    """
    Discretize a periodic potential function over a given interval.

    Parameters:
    -----------
    x_min : float
        Minimum x value of the interval
    x_max : float
        Maximum x value of the interval
    num_intervalles : int
        Number of intervals for discretization
    pot : callable
        Potential function to discretize
    L : float
        Period length of the potential
    N : int
        Number of periods
    params : tuple
        Parameters for the potential function

    Returns:
    --------
    tuple
        x_bins : ndarray
            Boundaries of the intervals
        V_bins : ndarray
            Potential values at the center of each interval
        x_centers : ndarray
            Centers of the intervals with boundary conditions applied
    """
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
    # Initialize V_bins to match x_centers size
    V_bins = np.zeros(len(x_centers))
    
    # Loop through x_centers indices only
    for i in range(len(x_centers)):
        V_bins[i] = pot(x_centers[i], params, N, L)
    
    V_bins[0] = 0.0
    V_bins[-1] = 0.0
    x_centers[0] = x_min
    x_centers[-1] = x_max
    
    return x_bins, V_bins, x_centers
