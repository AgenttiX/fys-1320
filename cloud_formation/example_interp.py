# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on Tampere University of Technology course FYS-1320 ESIMERKKI_interp.m by Sampo Saari 2014-10-14

# Our work is licensed with Creative Commons Attribution 4.0 International
# https://creativecommons.org/licenses/by/4.0/
# However, the license of the example code our work is based on is unclear and thereby so is the license
# for those parts of this code that are based on it

import numpy as np
import matplotlib.pyplot as plt
# import sympy
# import sympy.solvers

# Create x-vector (a Numpy array)
x_vec = np.arange(1, 10, 0.01)


# Function that cannot be analytically solved
def y(x):
    return (x**2) * np.log((np.sin(x)**2 + 1)*x) + 1/x + np.log(x) + 1e-4*(x**5)

# Compute y for the values of x-vector
y_vec = y(x_vec)

# Draw a plot of the function
plt.plot(x_vec, y_vec)

# Interpolate x at y=200 using the vectors
y_val = 200
x_val = np.interp(y_val, y_vec, x_vec)

# If you want to get rid of the PyCharm "Expected type 'Union[ndarray, iterable]', got 'int' instead" warning, you
# can use this version instead, but then the np.interp function returns the float within a list
# x_val = np.interp([y_val], y_vec, x_vec)

print(x_val)

# Show the plot
plt.show()

# Try to find an analytical solution (you'll see that sympy will produce an error)
# x_var = sympy.Symbol("x_var")
# x_val = sympy.solvers.solve( (x_var**2) * sympy.log((sympy.sin(x_var)**2 + 1)*x_var + 1/x_var + sympy.log(x_var) + 1e-4*(x_var**5) ) - y_val, x_var )
