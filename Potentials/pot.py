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


def pot_period(x, params, N, L):
    """
    Potentiel périodique de N puits Coulombiens équidistants.
    
    Paramètres:
    x : array
        Positions où calculer le potentiel
    params : tuple
        (Q, e, a, x0, eps0)
    N : int
        Nombre de puits
    L : float
        Distance entre les puits (en mètres)
    
    Retourne:
    V : array
        Potentiel total
    """
    Q, e, a, x0, eps0 = params
    
    # Positions des puits [x0, x0+L, x0+2L, ..., x0+(N-1)L]
    x_wells = x0 + np.arange(N) * L
    
    # Calcul du potentiel total
    V = np.zeros_like(x)
    for xi in x_wells:
        V += - (Q * e) / (4 * np.pi * eps0) / np.sqrt((x - xi)**2 + a**2)
    
    return V