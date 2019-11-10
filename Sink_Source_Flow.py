# Source and Sink flow

# Import Dependencies
import math
import numpy as np
import matplotlib.pyplot as plt


# Setup the grid used for the calculatios

N = 50
x_start, x_end = -2.0, 2.0
y_start, y_end = -1.0, 1.0
x = np.linspace(x_start, x_end, N)
y = np.linspace(y_start, y_end, N)
X, Y = np.meshgrid(x, y)    # Generates a meshgrid

# Calculate Source Flow
strength_source = 5.0
x_source, y_source = -1.0, 0   # Location of the Source

# Calculate the flow on the meshgrid

u_source = (strength_source / (2 * math.pi) *
            (X - x_source) / ((X - x_source)**2 + (Y - y_source)**2))
v_source = (strength_source / (2 * math.pi) *
            (Y - y_source) / ((X - x_source)**2 + (Y - y_source)**2))

# plot the streamlines
width = 10.0
height = (y_end - y_start) / (x_end - x_start) * width
plt.figure(figsize=(width, height))
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.xlim(x_start, x_end)
plt.ylim(y_start, y_end)
plt.streamplot(X, Y, u_source, v_source,
               density=2, linewidth=1, arrowsize=2, arrowstyle='->')
plt.scatter(x_source, y_source,
            color='#CD2305', s=80, marker='o')
plt.show()
print(v_source)
# Sink Calculation
strength_sink = -5.0
x_sink, y_sink = 1.0, 0

u_sink = (strength_sink / (2 * math.pi) *
           (X-x_sink) / ((X-x_sink)**2 + (Y-y_sink)**2))
v_sink = (strength_sink / (2 * math.pi) *
           (Y-y_sink) / ((X-x_sink)**2 + (Y-y_sink)**2))

# plot the streamlines

width = 10.0
height = (y_end - y_start) / (x_end - x_start) * width
plt.figure(figsize=(width, height))
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.xlim(x_start, x_end)
plt.ylim(y_start, y_end)
plt.streamplot(X, Y, u_sink, v_sink,
               density=2, linewidth=1, arrowsize=2, arrowstyle='->')
plt.scatter(x_sink, y_sink,
            color='#CD2305', s=80, marker='o')
plt.show()

u_pair = u_source+u_sink
v_pair = v_source+v_sink

width = 10.0
height = (y_end - y_start) / (x_end - x_start) * width
plt.figure(figsize=(width, height))
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.xlim(x_start, x_end)
plt.ylim(y_start, y_end)
plt.streamplot(X, Y, u_pair, v_pair,
               density=2, linewidth=1, arrowsize=2, arrowstyle='->')
plt.scatter(x_sink, y_sink,
            color='#CD2305', s=80, marker='o')
plt.scatter(x_source, y_source,
            color='#CD2305', s=80, marker='o')
plt.show()
