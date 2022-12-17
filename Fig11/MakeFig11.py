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
    LegendFSize=16
    AxisFSize=18
    TitleFSize=27
    XLim=[0,85]
    YLim=[0,9e12]
    YLimDiv=[-575,575]
    YTickLabelSize=17
    PanelLabel=[["(a)","(b)"],["(c)","(d)"],["(e)","(f)"]]
    LabelPos = [0.03, 0.92]
    fig,ax = CreateSuplotAxesSimple(2,2,10,20)
    # Subplot (0,0)
    PlotFiles =CreateTSInputAnom("nosurge","Hudson")
    ax[0][0].fill_between([20,88],[-1000,-1000],[9e15,9e15],
            color="grey", facecolor = 'grey',alpha=0.2)
    PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[0][0])
    ax[0][0].set_title("Hudson ice stream",fontweight="bold",fontsize=TitleFSize)
    ax[0][0].set_ylabel("Ice flux [m$^3$/yr]",fontsize=AxisFSize)
    ax[0][0].legend(loc=1,prop={'size':LegendFSize})
    ax[0][0].set_xlim(XLim)
    ax[0][0].set_ylim(YLim)
    ax[0][0].set_xticklabels([])
    ax[0][0].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[0][0].text(LabelPos[0],LabelPos[1], PanelLabel[0][0],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[0][0].transAxes,fontsize=AxisFSize)
    # Subplot (0,1)
    PlotFiles =CreateTSInputAnom("nosurge","Kenzie")
    print(PlotFiles)
    ax[0][1].fill_between([20,88],[-100,-100],[9e15,9e15],
            color="grey", facecolor = 'grey',alpha=0.2)
    PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[0][1])
    ax[0][1].set_title("Mackenzie ice stream",fontweight="bold",fontsize=TitleFSize)
    ax[0][1].legend(loc=1,prop={'size':LegendFSize})
    ax[0][1].set_xlim(XLim)
    ax[0][1].set_ylim(YLim)
    ax[0][1].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[0][1].set_yticklabels([])
    ax[0][1].set_xticklabels([])
    ax[0][1].text(LabelPos[0],LabelPos[1], PanelLabel[0][1],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[0][1].transAxes,fontsize=AxisFSize)
    # Subplot (1,0)
    ax[1][0].fill_between([20,88],[-1000,-1000],[9e15,9e15],
            color="grey", facecolor = 'grey',alpha=0.2)
    PlotFiles =CreateTSInputDivide("nosurge","Hudson")
    PlotPISMTimeSeries(PlotFiles,"divide_migration",ax[1][0])
    ax[1][0].set_ylabel("Divide migration [km]",fontsize=AxisFSize)
    ax[1][0].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax[1][0].legend(loc=1,prop={'size':LegendFSize})
    ax[1][0].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[1][0].set_xlim(XLim)
    ax[1][0].set_ylim(YLimDiv)
    AddArrowsDS(ax[1][0],loc="upper")
    # Subplot (1,1)
    ax[1][0].text(-0.06,0.95, "West",color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[1][0].transAxes,fontsize=AxisFSize)
    ax[1][0].text(-0.06,0.05, "East",color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[1][0].transAxes,fontsize=AxisFSize)
    ax[1][0].text(LabelPos[0],LabelPos[1], PanelLabel[1][0],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[1][0].transAxes,fontsize=AxisFSize)
    PlotFiles =CreateTSInputDivide("nosurge","Kenzie")
    ax[1][1].fill_between([20,88],[-1000,-1000],[9e15,9e15],
            color="grey", facecolor = 'grey',alpha=0.2)
    PlotPISMTimeSeries(PlotFiles,"divide_migration",ax[1][1])
    ax[1][1].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax[1][1].legend(loc=1,prop={'size':LegendFSize})
    ax[1][1].set_xlim(XLim)
    ax[1][1].set_ylim(YLimDiv)
    ax[1][1].set_yticklabels([])
    ax[1][1].tick_params(axis='both',labelsize=YTickLabelSize)
    AddArrowsDS(ax[1][1])
    ax[1][1].text(LabelPos[0],LabelPos[1], PanelLabel[1][1],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[1][1].transAxes,fontsize=AxisFSize)
    plt.tight_layout()
    # plt.show()
    # sys.exit()

    SavePlot("Fig11","SurgesNoSurgesAndDivide")

if __name__ == '__main__':
    Main()
