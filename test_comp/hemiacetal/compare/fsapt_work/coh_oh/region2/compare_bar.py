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
elst = (-31.07, -26.52, -26.98, -27.07)
exch = (65.33, 57.65, 57.91, 57.93)
indAB = (-24.60, -21.41, -21.48, -21.50)
indBA = (-22.72, -20.67, -21.00, -21.09)
disp = (-6.67, -5.85, -5.87, -5.90)
total = (-19.72, -16.81, -17.42, -17.63)
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

plt.ylabel('$w_2$ (kcal mol$^{-1}$)', fontsize=14)
plt.title('COH --- OH (Region 2)', fontsize=14)

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
