import numpy as np  # Importing the numpy library for numerical operations
import matplotlib.pyplot as plt  # Importing matplotlib for plotting (not used in this code)

def pot(x, params):
    """
    Calculate the potential at a point x due to a charge distribution.

    Parameters:
    x : float
        The position where the potential is calculated.
    params : tuple
        A tuple containing the parameters (Q, e, a, x0, eps0):
        Q : float
            Charge magnitude.
        e : float
            Elementary charge.
        a : float
            A constant related to the charge distribution.
        x0 : float
            The position of the charge distribution.
        eps0 : float
            Permittivity of free space.

    Returns:
    float
        The potential at position x.
    """
    Q, e, a, x0, eps0 = params  # Unpacking the parameters
    # Calculate the potential using the given formula
    return -((Q * e) / (4 * np.pi * eps0)) / np.sqrt((x - x0) ** 2 + a ** 2)

def pot_period(x,N,L,params):
    Q, e, a, x0, eps0 = params
    x_wells = np.array([n * L for n in range(N)])
    V_wells = np.zeros_like(x)
    for x_i in x_wells:
        V_wells += -((Q * e) / (4 * np.pi * eps0)) / np.sqrt((x - x0) ** 2 + a ** 2)

    V_interaction = np.zeros_like(x)

    for i in range(len(x_wells) - 1):
        x_i, x_j = x_wells[i], x_wells[i+1]
        r_ij = np.abs(x_i - x_j)
        V_interaction += -((Q * e) / (4 * np.pi * eps0)) / np.sqrt((x - x0) ** 2 + a ** 2)
    return V_interaction + V_wells


