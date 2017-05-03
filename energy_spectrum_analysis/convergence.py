# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# This program creates a 3D graph of Monte-Carlo integral error as a function of both the point and step counts

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
import matplotlib.pyplot as plt

import plotly.offline
import plotly.graph_objs as go

stepcounts = numpy.logspace(2, 6, num=6)
pointcounts = numpy.logspace(2, 8)  # MAX 8
print(pointcounts)
y_max = 1600
cauchy_params = [8.32514938e-15, 1.88310185e-16, 9.32332574e-13]
print("Analytical:", cauchy_params[2])

# x axis limits
a = 6.08859964548e-15
b = 9.3720685487e-15

x_diff = b-a

results = numpy.zeros((stepcounts.size, pointcounts.size))

for i_stepcount, stepcount in enumerate(stepcounts):
    x = np.arange(start=a, stop=b, step=(b-a)/stepcount)
    cauchy_vec = toolbox.cauchy(x, cauchy_params[0], cauchy_params[1], cauchy_params[2])

    for i_pointcount, pointcount in enumerate(pointcounts):
        x_points = a + numpy.random.rand(pointcount)*(b-a)
        y_points = numpy.random.rand(pointcount)*y_max

        cauchy_points = toolbox.cauchy(x_points, cauchy_params[0], cauchy_params[1], cauchy_params[2])

        undervec = np.greater(cauchy_points, y_points)
        under = np.sum(undervec)

        result = (under/pointcount)*(b-a)*y_max

        results[i_stepcount, i_pointcount] = result

        print("Steps", stepcount)
        print("Points", pointcount)
        print("Result", result)

"""
under1 = 0
for i in range(x.size):
    if y_points[i] < toolbox.cauchy(x_points[i], cauchy_params[0], cauchy_params[1], cauchy_params[2]):
        under1 += 1
        #plt.plot(x_points[i], y_points[i], "ro")
    else:
        pass
        #plt.plot(x_points[i], y_points[i], "ko")
print(under1)
"""

plt.plot(cauchy_vec)
plt.show()

plotly.offline.init_notebook_mode()

# surf = go.Surface(z=results, x=pointcounts, y=stepcounts)
surf = go.Surface(
    x=np.log10(pointcounts),
    y=np.log10(stepcounts),
    z=results*1e13
)

layout = go.Layout(
    title="Monte Carlo",
    xaxis=dict(title="foo")
    # xaxis=dict(type="log", autorange=True),
    # yaxis=dict(type="log", autorange=True)
)

figure = go.Figure(data=[surf], layout=layout)
plotly.offline.plot(figure, filename="test.html")

print(results.max())
print(results.min())
