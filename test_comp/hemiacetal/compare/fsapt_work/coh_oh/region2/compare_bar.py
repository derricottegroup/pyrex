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
elst = (-31.07, -26.52, -26.98, -27.07)
exch = (65.33, 57.65, 57.91, 57.93)
indAB = (-24.60, -21.41, -21.48, -21.50)
indBA = (-22.72, -20.67, -21.00, -21.09)
disp = (-6.67, -5.85, -5.87, -5.90)
total = (-19.72, -16.81, -17.42, -17.63)

F = (4.55, -7.68, 3.19, 2.05, 0.82, 2.91)
Cl = (4.09, -7.42, 3.12, 1.72, 0.8, 2.3)
Br = (4, -7.40, 3.10, 1.63, 0.77, 2.09)

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

plt.ylabel('$\\Delta w_2$ (kcal mol$^{-1}$)', fontsize=10)
plt.title('COH --- OH (Region 2)', fontsize=10)

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

