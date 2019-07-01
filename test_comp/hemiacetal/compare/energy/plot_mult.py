import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import pandas as pd
import os
import numpy as np

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

########################################################################################
fontsize = 5.0
ticklabel_size = 8.0
linewidth = 1.0
axiswidth = 1.0
fig_dims = [3.25, 1.75]

csv_files = ['energy_ch3.csv', 'energy_cf3.csv', 'energy_ccl3.csv', 'energy_cbr3.csv']
legend_labels = ['R1', 'R2', 'R3', 'R4']

x_label = "$\\xi$ (au amu$^{1/2}$)"
y_label = "$\\Delta E$ (kcal mol$^{-1}$)"

figure_output_filename = "energy.pdf"
#########################################################################################

gold = '#ffbd39'
turquoise = '#08d9d6'
black = '#252a34'
red = '#ff2e63'

legend_properties = {'weight':'bold'}
colors = [gold,turquoise,red,black]


x = []
y_kcal = []
hartree2kcal = 627.509

for filename in csv_files:
    df_ch3 = pd.read_csv(filename)
    x.append(df_ch3["Coordinate"])
    y = df_ch3["Energy"]
    y_kcal.append((y - y[0])*hartree2kcal)


for i in range(len(legend_labels)):
    plt.plot(x[i],y_kcal[i],linewidth=linewidth,label=legend_labels[i], color=colors[i])

rc('axes', linewidth=axiswidth)
ax = gca()
simpleaxis(ax)
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
plt.legend(frameon=False, fontsize=fontsize, prop=legend_properties)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(fig_dims[0], fig_dims[1])

plt.ylabel(y_label, fontsize=14)
plt.xlabel(x_label, fontsize=14)


plt.savefig(figure_output_filename)
