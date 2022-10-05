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
    XLim=[0,85]
    YLim=[0,9e12]
    YTickLabelSize=20
    LabelPos = [0.03, 0.975]
    PanelLabel=[["(a)","(b)"],["(c)","(d)"],["(e)","(f)"]]
    fig,ax = CreateSuplotAxesSimple(1,2,10,20)
    PlotFiles=CreateTSInputAnom("synchron","Dummy")
    ax[0].fill_between([20,88],[-1000,-1000],[9e15,9e15],
            color="grey", facecolor = 'grey',alpha=0.2)
    PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[0])
    ax[0].set_title("Equal surge interval",fontweight="bold",fontsize=TitleFSize)
    ax[0].legend(loc=1,prop={'size':LegendFSize})
    ax[0].set_ylabel("Ice flux [m$^3$/yr]",fontsize=AxisFSize)
    ax[0].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax[0].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[0].set_xlim(XLim)
    ax[0].set_ylim(YLim)
    ax[0].text(LabelPos[0],LabelPos[1], PanelLabel[0][0],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[0].transAxes,fontsize=AxisFSize)
    PlotFiles[0][0]="Data/HudsonHovShifted_HE12.nc"
    ax[1].fill_between([20,88],[-1000,-1000],[9e15,9e15],
            color="grey", facecolor = 'grey',alpha=0.2)
    PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[1],-3.9)
    ax[1].set_title("Synchronous and shifted surges",fontweight="bold",fontsize=TitleFSize)
    ax[1].legend(loc=1,prop={'size':LegendFSize})
    ax[1].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax[1].set_xlim(XLim)
    ax[1].set_yticklabels([])
    ax[1].set_ylim(YLim)
    ax[1].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[1].text(LabelPos[0],LabelPos[1], PanelLabel[0][1],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[1].transAxes,fontsize=AxisFSize)
    plt.tight_layout()
    # plt.show()
    # sys.exit()
    SavePlot("Fig12","SurgeSynchronicity")

if __name__ == '__main__':
    Main()
