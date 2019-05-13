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
elst = (-71.43, -82.84, -88.59, -87.58)
exch = (138.64, 167.27, 176.17, 173.29)
indAB = (-16.97, -22.61, -22.60, -22.02)
indBA = (-30.00, -43.05, -44.84, -44.04)
disp = (-11.76, -13.59, -14.38, -14.23)
total = (8.48, 5.18, 5.75, 5.42)
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
plt.title('COH --- OH (Region 1)', fontsize=14)

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
