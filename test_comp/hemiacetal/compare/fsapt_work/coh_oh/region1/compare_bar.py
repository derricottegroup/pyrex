import numpy as np
from pylab import *
import matplotlib
import matplotlib.pyplot as plt

#######################
fontsize = 5.0 
ticklabel_size = 8.0
linewidth = 1.0
axiswidth = 1.0
########################

N = 6
elst = (-71.43, -82.84, -88.59, -87.58)
exch = (138.64, 167.27, 176.17, 173.29)
indAB = (-16.97, -22.61, -22.60, -22.02)
indBA = (-30.00, -43.05, -44.84, -44.04)
disp = (-11.76, -13.59, -14.38, -14.23)
total = (8.48, 5.18, 5.75, 5.42)

F = (-11.41, 28.63, -5.64, -13.05, -1.83, -3.3)
Cl = (-17.16, 37.53, -5.63, -14.84, -2.62, -2.73)
Br = (-16.15, 34.65, -5.05, -14.04, -2.47, -3.06)

men_means = (20, 35, 30, 35, 27)
women_means = (25, 32, 34, 20, 25)

ind = np.arange(N) 
width = 0.2       
plt.bar(ind, F, width, label='X = F (R2)', color="red")
plt.bar(ind + width, Cl, width,
    label='X = Cl (R3)', color="green")
plt.bar(ind + 2.0*width, Br, width,
    label='X = Br (R4)', color="blue")


plt.axhline(y=0,color='black',linestyle="--",linewidth=1.5)

plt.ylabel('$\\Delta w_1$ (kcal mol$^{-1}$)', fontsize=10)
plt.title('COH --- OH (Region 1)', fontsize=10)

plt.xticks(ind + 3.0*width / 2, ('elst', 'exch', 'indAB', 'indBA', 'disp', 'total'))
plt.legend(loc='best')
plt.ylim([-20,40])
rc('axes', linewidth=axiswidth)
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(ticklabel_size)
    #tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(ticklabel_size)
    #tick.label1.set_fontweight('bold')
for axis in ['bottom','left']:
    ax.spines[axis].set_linewidth(axiswidth)
ax.xaxis.set_tick_params(width=axiswidth)
ax.yaxis.set_tick_params(width=axiswidth)
plt.legend()
plt.legend(frameon=False, fontsize=fontsize)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(3.5, 2.0)
plt.savefig("bar_plot.pdf")
