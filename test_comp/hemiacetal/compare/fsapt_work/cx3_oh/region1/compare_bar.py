import numpy as np
from pylab import *
import matplotlib
import matplotlib.pyplot as plt

#######################
fontsize = 10 
ticklabel_size = 16
linewidth = 1.0
axiswidth = 1.0
########################

N = 4
elst = (-5.51, -9.48, -10.45, -10.13)
exch = (10.50, 11.59, 12.56, 12.84)
indAB = (-0.80, -0.73, -0.90, -1.07)
indBA = (-3.98, -9.80, -9.77, -9.11)
disp = (-1.59, -1.73, -2.41, -2.64)
total = (-1.37, -10.15, -10.96, -10.12)
men_means = (20, 35, 30, 35, 27)
women_means = (25, 32, 34, 20, 25)

ind = np.arange(N) 
width = 0.1       
plt.bar(ind, elst, width, label='elst', color="red")
plt.bar(ind + width, exch, width,
    label='exch', color="green")
plt.bar(ind + 2.0*width, indAB, width,
    label='indAB', color="blue")
plt.bar(ind + 3.0*width, indBA, width,
    label='indBA', color="orange")
plt.bar(ind + 4.0*width, disp, width,
    label='disp', color="magenta")
plt.bar(ind + 5.0*width, total, width,
    label='total', color="black")


plt.axhline(y=0,color='black',linestyle="--",linewidth=1.5)

plt.ylabel('$w_1$ (kcal mol$^{-1}$)', fontsize=14)
plt.title('CX$_3$ --- OH (Region 1)', fontsize=14)

plt.xticks(ind + 5.0*width / 2, ('X=H', 'X=F', 'X=Cl', 'X=Br'))
plt.legend(loc='best')
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
fig.set_size_inches(9.5, 4.5)
plt.savefig("bar_plot.pdf")
