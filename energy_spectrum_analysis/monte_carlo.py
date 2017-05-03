# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# This program creates a 2D graph of a Monte-Carlo integral

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import toolbox
import numpy as np
import numpy.random

import pyqtgraph as pg

stepcount = int(1e4)
pointcount = int(1e4)  # MAX 8
y_max = 1600
cauchy_params = [8.32514938e-15, 1.88310185e-16, 9.32332574e-13]
print("Analytical:", cauchy_params[2])

# x axis limits
a = 6.08859964548e-15
b = 9.3720685487e-15

x_diff = b-a


x = np.arange(start=a, stop=b, step=(b-a)/stepcount)
cauchy_vec = toolbox.cauchy(x, cauchy_params[0], cauchy_params[1], cauchy_params[2])


x_points = a + numpy.random.rand(pointcount)*(b-a)
y_points = numpy.random.rand(pointcount)*y_max

cauchy_points = toolbox.cauchy(x_points, cauchy_params[0], cauchy_params[1], cauchy_params[2])

undervec = np.greater(cauchy_points, y_points)
under = np.sum(undervec)

result = (under/pointcount)*(b-a)*y_max

print("Steps", stepcount)
print("Points", pointcount)
print("Result", result)

points = np.array([x_points, y_points])

x_points_under = x_points[undervec == True]
x_points_over = x_points[undervec == False]

y_points_under = y_points[undervec == True]
y_points_over = y_points[undervec == False]

qapp = pg.mkQApp()

pg.setConfigOptions(antialias=True, background="w", foreground="k")

win = pg.GraphicsWindow(title="Monte Carlo")
plot = win.addPlot(title="Monte Carlo")

plot.plot(x_points_over/1.60218e-16, y_points_over, pen=None, symbol="o", symbolPen=(180, 220, 0))
plot.plot(x_points_under/1.60218e-16, y_points_under, pen=None, symbol="o", symbolPen=(70, 200, 255))
plot.plot(x/1.60218e-16, cauchy_vec, pen=pg.mkPen((255, 0, 0), width=2))

qapp.exec_()
