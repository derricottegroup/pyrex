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
elst = (-0.00, 0.18, 0.21, 0.23)
exch = (2.25, 2.07, 2.07, 2.09)
indAB = (-0.35, -0.29, -0.30, -0.37)
indBA = (-2.11, -3.46, -3.22, -3.01)
disp = (-0.40, -0.35, -0.45, -0.50)
total = (-0.61, -1.84, -1.69, -1.56)

F = (0.18, -0.18, 0.06, -1.35, 0.05, -1.23)
Cl = (0.21, -0.18, 0.05, -1.11, -0.05, -1.08)
Br = (0.23, -0.16, -0.02, -0.9, -0.1, -0.95)

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

