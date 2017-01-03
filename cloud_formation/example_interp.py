# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on Tampere University of Technology course FYS-1320 ESIMERKKI_interp.m by Sampo Saari 2014-10-14

import numpy as np
import matplotlib.pyplot as plt
import sympy
import sympy.solvers

x = np.arange(1, 10, 0.01)


def y(x):
    return (x**2) * np.log((np.sin(x)**2 + 1)*x + 1/x + np.log(x) + 1e-4*(x**5))

plt.plot(x, y(x))
plt.show()

y_val = 200

x_var = sympy.Symbol("x_var")
x_val = sympy.solvers.solve( (x_var**2) * sympy.log((sympy.sin(x_var)**2 + 1)*x_var + 1/x_var + sympy.log(x_var) + 1e-4*(x_var**5) ) - y_val, x_var )

print(x_val)

# TODO
