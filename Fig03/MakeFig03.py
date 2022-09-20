#! /usr/bin/env python
import sys
sys.path.append("../util")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os
import errno

# Local imports
from UtilityFcts import *

def Main():
    LegendFSize=24
    AxisFSize=27
    TitleFSize=31
    YLim=[0,6e12]
    YTickLabelSize=24
    fig,ax = CreateSuplotAxesSimple(1,1,10,15)
    PlotFiles=CreateTSInputAnom("reference","Hudson")
    PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax,20)
    PlotFiles=CreateTSInputAnom("reference","Kenzie")
    PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax,20)
    # ax.fill_between([20,88],[-100,-100],[9e12,9e12],
                # color="grey", hatch = '//', label='Analysis Period',
                # facecolor = 'none',alpha=0.6)
    ax.set_ylabel("Ice flux [m$^3$/yr]",fontsize=AxisFSize)
    # ax.set_xlim([-65,-7.5])
    ax.set_xlim([0,65.0])
    ax.set_ylim(YLim)
    ax.legend(loc=1,prop={'size':LegendFSize})
    ax.text(0.21, 0.95, "Control simulation (Ctrl)",color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax.transAxes,fontsize=LegendFSize)
    ax.set_title("Surge frequency",fontweight="bold",fontsize=TitleFSize)
    ax.set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax.tick_params(axis='both',labelsize=YTickLabelSize)
    t = ax.yaxis.get_offset_text()
    t.set_size(LegendFSize)
    plt.tight_layout()
    SavePlot("Fig03","TSSurges_ReferenceSim")

if __name__ == '__main__':
    Main()
