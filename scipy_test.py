
# Test file for control library

# import dependencies
import numpy as np
import scipy
import matplotlib.pyplot as plt
from math import e
from control.matlab import *


# Create System

sys1 = tf(6,[1,6,11,6])

response = step(sys1)

plt.figure('Step Response')
plt.plot(response[1], response[0])
plt.grid()

plt.figure('Bode Plot')
bode(sys1)
plt.show()
