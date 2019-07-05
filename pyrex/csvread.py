import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from pylab import *
import json
import pandas as pd
import numpy as np

class Params():
    def __init__(self):
        self.do_saveplot = False
        self.rel_energy = False
        self.force_min = 0.0
        self.force_max = 0.0
        self.scale = 1.0
        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None
        self.user_color = False
        self.color = []
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
            if 'color' in input_params['rexplot']:
                self.user_color = True
                self.color = input_params['rexplot']['color']
            if 'coordinate' in input_params['rexplot']:
                self.coordinate = input_params['rexplot']['coordinate']
            if 'force_min' in input_params['rexplot']:
                self.force_min = input_params['rexplot']['force_min']
            if 'force_max' in input_params['rexplot']:
                self.force_max = input_params['rexplot']['force_max']
            if 'rel_energy' in input_params['rexplot']:
                self.rel_energy = bool(input_params['rexplot']['rel_energy'])
            if 'scale' in input_params['rexplot']:
                self.scale = input_params['rexplot']['scale'] 
            if 'x_min' in input_params['rexplot']:
                self.x_min = input_params['rexplot']['x_min']
            if 'x_max' in input_params['rexplot']:
                self.x_max = input_params['rexplot']['x_max']
            if 'y_min' in input_params['rexplot']:
                self.y_min = input_params['rexplot']['y_min']
            if 'y_max' in input_params['rexplot']:
                self.y_max = input_params['rexplot']['y_max']
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
    y = []
    for i in range(len(params.properties)):
        array_temp = []
        current_array = df[params.properties[i]]
        if(params.rel_energy):
            for j in range(len(current_array)):
                array_temp.append((current_array[j] - current_array[0])*params.scale)
            y.append(array_temp)
        else:
            y.append(df[params.properties[i]]*params.scale)
    
    y_spline = []
    for i in range(len(params.properties)):
        y_fit = UnivariateSpline(x, y[i], s=0)
        y_spline.append(y_fit)

    xs = np.linspace(x[0], x[len(x)-1], 1000)
    
    y_spline_fit = []
    for i in range(len(params.properties)):
        ys_ = y_spline[i](xs)
        y_spline_fit.append(ys_)
    
    fontsize = params.fontsize
    plt.ylabel(params.y_label, fontsize=fontsize)
    plt.xlabel(params.x_label, fontsize=fontsize)
    if(params.force_min!=0):
        axvline(x=params.force_min, linewidth=3.0, color='k', linestyle=':')
    if(params.force_max!=0):
        axvline(x=params.force_max, linewidth=3.0, color='k', linestyle=':')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    for i in range(len(params.properties)):
        if(params.user_color==True):
            plt.plot(xs, y_spline_fit[i], linewidth=params.linewidth, label=params.properties[i], color=params.color[i])
        else:
            plt.plot(xs, y_spline_fit[i], linewidth=params.linewidth, label=params.properties[i])
    rc('axes', linewidth=params.axiswidth)
    ax = gca()
    simpleaxis(ax) 
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(params.ticklabel_size)
        #tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(params.ticklabel_size)
        #tick.label1.set_fontweight('bold')
    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(params.axiswidth)
    ax.xaxis.set_tick_params(width=3)
    ax.yaxis.set_tick_params(width=3)
    plt.legend(prop={'size': 60})
    plt.legend(frameon=False)
    if(params.y_max and params.y_min):
        plt.ylim(params.y_min, params.y_max)
    else:
        pass
    if(params.x_max!=None and params.x_min!=None):
        plt.xlim(params.x_min, params.x_max)
    else:
        pass
    fig_dims = params.fig_dims
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(fig_dims[0], fig_dims[1])
    
    if(params.do_saveplot):
        plt.savefig(params.plot_file)
    else:
        plt.show()
