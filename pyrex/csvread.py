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
            if 'plot_file' in input_params['rexplot']:
                self.do_saveplot = True
                self.plot_file = input_params['rexplot']['plot_file']


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

    plt.ylabel(params.y_label, fontsize=30)
    plt.xlabel(params.x_label, fontsize=30)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    for i in range(len(params.properties)):
        plt.plot(xs, y_spline_fit[i], linewidth=4.0, label=params.properties[i])

    rc('axes', linewidth=3)
    fontsize = 20
    ax = gca()
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    ax.xaxis.set_tick_params(width=3)
    ax.yaxis.set_tick_params(width=3)
    plt.legend(prop={'size': 15})
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(9.0, 5.0)
    
    if(params.do_saveplot):
        plt.savefig(params.plot_file)
    else:
        plt.show()
    







































