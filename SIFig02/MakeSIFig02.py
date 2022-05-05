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
    LegendFSize=18
    AxisFSize=22
    TitleFSize=27
    YLim=[0,6e8]
    YTickLabelSize=20
    XLim=[0,57.6]
    fig,ax = CreateSuplotAxesSimple(1,1,10,15)
    PlotFiles=CreateTSInputAnom("oce","Hudson")
    PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax)
    ax.set_ylabel("Ice flux [m$^2$/yr]",fontsize=AxisFSize)
    ax.set_xlim(XLim)
    ax.set_ylim(YLim)
    ax.legend(loc=1,prop={'size':LegendFSize})
    ax.set_title("Hudson ice stream",fontweight="bold",fontsize=TitleFSize)
    ax.set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax.tick_params(axis='both',labelsize=YTickLabelSize)
    plt.tight_layout()
    SavePlot("SIFig02","TSSurges_OceHudson")

if __name__ == '__main__':
    Main()
