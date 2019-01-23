import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from pylab import *
import json
import pandas as pd
import numpy as np


df_ch3 = pd.read_csv('sapt_force_ch3.csv')
df_ccl3 = pd.read_csv('sapt_force_ch3.csv')

x_ch3 = df_ch3['Coordinate']
x_ccl3 = df_ccl3['Coordinate']

y_ch3 = df_ch3['F_elst']
y_ccl3 = df_ccl3['F_elst']

y_fit_ch3 = UnivariateSpline(x_ch3, y_ch3, s=0)

y_fit_ccl3 = UnivariateSpline(x_ccl3, y_ccl3, s=0)

xs_ch3 = np.linspace(x_min, x_max, 1000)

xs_ccl3 = np.linspace(x_min, x_max, 1000)


ys_ch3 = y_fit_ch3(xs_ch3)

ys_ccl3 = y_fit_ccl3(xs_ccl3)

plt.plot(xs_ch3, ys_ch3)
plt.plot(xs_ccl3, ys_ccl3)
plt.show()
