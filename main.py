import pot 
import discret
import numpy as np
import matplotlib.pyplot as plt

def main(): 
    a = 0.01e-9
    e = 1.6e-19
    Q = e
    eps0 = 8.85e-12
    x0 = 0
    params = (Q, e, a, x0, eps0)
    L = 0.05e-9
    x_min = -L
    x_max = L
    num_intervalles = 100
    x_steps, V_steps = discret.discretiser(x_min, x_max, num_intervalles, pot.pot, params)
    x = np.linspace(x_min, x_max, 1000)
    plt.plot(x_steps, V_steps, label="Discrete Steps")
    plt.plot(x, pot.pot(x, params), label="Continuous Function", linestyle='-')
    plt.legend()
    plt.show()  

if __name__ == "__main__":
    main()

