import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from pylab import *
import json
import pandas as pd
import numpy as np

class Params():
    def __init__(self):
        self.do_saveplot = False
        self.fontsize = 16
        self.ticklabel_size = 10
        self.linewidth = 5.0
        self.axiswidth = 2.0
        self.fig_dims = [9.0, 5.0]
        self.coordinate = 'Coordinate'
        json_input = sys.argv[1]
        self.read_input(json_input)
    def read_input(self, json_input):
        json_data=open(json_input).read()
        input_params = json.loads(json_data)
        if 'rexplot' in input_params:
            if 'file' in input_params['rexplot']:
                self.filename = input_params['rexplot']['file']
            if 'properties' in input_params['rexplot']:
                self.properties = input_params['rexplot']['properties']
            if 'coordinate' in input_params['rexplot']:
                self.coordinate = input_params['rexplot']['coordinate']
            if 'x_label' in input_params['rexplot']:    
                self.x_label = input_params['rexplot']['x_label']
            if 'y_label' in input_params['rexplot']:
                self.y_label = input_params['rexplot']['y_label']
            if 'fontsize' in input_params['rexplot']:
                self.fontsize = input_params['rexplot']['fontsize']
            if 'linewidth' in input_params['rexplot']:
                self.linewidth = input_params['rexplot']['linewidth']
            if 'axiswidth' in input_params['rexplot']:
                self.axiswidth = input_params['rexplot']['axiswidth']
            if 'ticklabel_size' in input_params['rexplot']:
                self.ticklabel_size = input_params['rexplot']['ticklabel_size']
            if 'fig_dims' in input_params['rexplot']:
                self.fig_dims = input_params['rexplot']['fig_dims']
            if 'plot_file' in input_params['rexplot']:
                self.do_saveplot = True
                self.plot_file = input_params['rexplot']['plot_file']


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def plot():
    params = Params()
    input_file = params.filename
    df = pd.read_csv(input_file)

    x = df[params.coordinate]
    x_min = x[0]
    x_max = x[len(x)-1]
    y = []
    for i in range(len(params.properties)):
        y.append(df[params.properties[i]])
    
    y_spline = []
    for i in range(len(params.properties)):
        y_fit = UnivariateSpline(x, y[i], s=0)
        y_spline.append(y_fit)

    xs = np.linspace(x_min, x_max, 1000)
    
    y_spline_fit = []
    for i in range(len(params.properties)):
        ys_ = y_spline[i](xs)
        y_spline_fit.append(ys_)
    
    fontsize = params.fontsize
    plt.ylabel(params.y_label, fontsize=fontsize)
    plt.xlabel(params.x_label, fontsize=fontsize)

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    for i in range(len(params.properties)):
        plt.plot(xs, y_spline_fit[i], linewidth=params.linewidth, label=params.properties[i])

    rc('axes', linewidth=params.axiswidth)
    ax = gca()
    simpleaxis(ax) 
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(params.ticklabel_size)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(params.ticklabel_size)
        tick.label1.set_fontweight('bold')
    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(params.axiswidth)
    ax.xaxis.set_tick_params(width=3)
    ax.yaxis.set_tick_params(width=3)
    plt.legend(prop={'size': 60})
    plt.legend(frameon=False)
    fig_dims = params.fig_dims
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(fig_dims[0], fig_dims[1])
    
    if(params.do_saveplot):
        plt.savefig(params.plot_file)
    else:
        plt.show()
