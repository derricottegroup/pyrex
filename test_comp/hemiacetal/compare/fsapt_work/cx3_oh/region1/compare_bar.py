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
elst = (-5.51, -9.48, -10.45, -10.13)
exch = (10.50, 11.59, 12.56, 12.84)
indAB = (-0.80, -0.73, -0.90, -1.07)
indBA = (-3.98, -9.80, -9.77, -9.11)
disp = (-1.59, -1.73, -2.41, -2.64)
total = (-1.37, -10.15, -10.96, -10.12)

F = (-3.97, 1.09, 0.07, -5.82, -0.14, -8.78)
Cl = (-4.94, 2.06, -0.1, -5.79, -0.82, -9.59)
Br = (-4.62, 2.34, -0.27, -5.13, -1.05, -8.75)

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
plt.title('CX$_3$ --- OH (Region 1)', fontsize=10)

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

