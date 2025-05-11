import numpy as np 

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
    Q, e, a, x0, eps0 = params  
   
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
    
    
    x_wells = x0 + np.arange(N) * L
    
    
    V = np.zeros_like(x)
    for xi in x_wells:
        V += - (Q * e) / (4 * np.pi * eps0)/np.sqrt((x - xi)**2 + a**2)
    return V

def potential_voisin(x, params, N, L):
    """
    Calcule le potentiel périodique 1D pour un cristal avec N puits coulombiens adoucis,
    en considérant uniquement le puits le plus proche de chaque position x.

    Paramètres :
    x : array_like      → Positions (en mètres) où calculer le potentiel.
    Q : float           → Charge des noyaux (en multiples de la charge élémentaire e).
    e : float           → Charge élémentaire (1.602e-19 C).
    epsilon_0 : float   → Permittivité du vide (8.854e-12 F/m).
    a : float           → Paramètre d'adoucissement (en mètres).
    L : float           → Distance entre les puits (en mètres).
    N : int             → Nombre de puits.
    
    Retourne :
    V : ndarray         → Potentiel aux positions x (en volts).
    """
    
    Q, e, a, x0, eps0 = params
    x = np.asarray(x)
    V = np.zeros_like(x)
    x_wells = x0 + np.arange(N)* L  
    

    k_floor = np.floor((x - x0) / L).astype(int)
    left_well = np.clip(k_floor, 0, N-1)

    right_well = np.clip(k_floor + 1, 0, N-1)
    x_left = x_wells[left_well]
    dl = np.sqrt((x-x_left)**2 + a**2)
    x_right = x_wells[right_well]
    dr = np.sqrt((x-x_right)**2 + a**2)
    dist_l = np.sqrt(x**2 +x_left**2)
    dist_r = np.sqrt(x**2 + x_right**2)
    if dr < dl : 
        V +=  - (Q * e) / (4 * np.pi * eps0 )/dr
    else:
        V += - (Q * e) / (4 * np.pi * eps0 )/dl
    return V

