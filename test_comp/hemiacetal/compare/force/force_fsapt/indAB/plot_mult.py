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

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

fontsize = 8
ticklabel_size = 7 
linewidth = 1.0
axiswidth = 1.0
fig_dims = [3.5, 2.2]

hartree2kcal = 627.509 

csv_files = ['force_indab_ch3.csv', 'force_indab_cf3.csv', 'force_indab_ccl3.csv', 'force_indab_cbr3.csv']
legend_labels = ['X=H', 'X=F', 'X=Cl', 'X=Br']

x_label = "$\\xi$ (au amu$^{1/2}$)"
y_label = "$F_{\\rm indAB}$ (kcal mol$^{-1}$ $\\xi^{-1}$)"

figure_output_filename = "force.svg"

x = []
y_kcal = []

names = ["CH3-OH", "CF3-OH", "CCl3-OH", "CBr3-OH"]
count = 0 
for filename in csv_files:
    df_ch3 = pd.read_csv(filename)
    x.append(df_ch3["Coordinate"])
    y = df_ch3[names[count]]
    y_kcal.append((y - y[0]))
    count = count + 1

for i in range(len(legend_labels)):
    plt.plot(x[i],y_kcal[i],linewidth=linewidth,label=legend_labels[i])

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
ax.xaxis.set_tick_params(width=1)
ax.yaxis.set_tick_params(width=1)
plt.legend()
plt.legend(frameon=False, fontsize=fontsize)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(fig_dims[0], fig_dims[1])

plt.ylabel(y_label, fontsize=fontsize)
plt.xlabel(x_label, fontsize=fontsize)


plt.savefig(figure_output_filename)
