# live plot tracking practice

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import PID
import numpy as np
from scipy.interpolate import BSpline, make_interp_spline


def live_pid(P=0.1, I=0.1, D=0.1)

    pid = PID.PID(P, I, D)
    pid.Setpoint =
