import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from pylab import *
import json
import pandas as pd
import numpy as np


df_ch3 = pd.read_csv('sapt_force_ch3.csv')
df_ccl3 = pd.read_csv('sapt_force_ccl3.csv')

x_ch3 = df_ch3['Coordinate']
x_ccl3 = df_ccl3['Coordinate']

y_ch3 = df_ch3['F_elst']
y_ccl3 = df_ccl3['F_elst']

y_fit_ch3 = UnivariateSpline(x_ch3, y_ch3, s=0)

y_fit_ccl3 = UnivariateSpline(x_ccl3, y_ccl3, s=0)

xs_ch3 = np.linspace(-6.0, 5.05, 1000)

xs_ccl3 = np.linspace(-6.45, 3.90, 1000)


ys_ch3 = y_fit_ch3(xs_ch3)

ys_ccl3 = y_fit_ccl3(xs_ccl3)

plt.plot(xs_ch3, ys_ch3, label='X=H')
plt.plot(xs_ccl3, ys_ccl3, label='X=Cl')

plt.legend(prop={'size': 60})
plt.legend(frameon=False)

plt.show()
