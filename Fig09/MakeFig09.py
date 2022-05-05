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
    LegendFSize=20
    AxisFSize=24
    TitleFSize=27
    XLim=[0,57.6]
    LabelPos = [0.03, 0.975]
    T = np.linspace(0,58,100)
    # YLim=[0,6e8]
    YTickLabelSize=20
    PanelLabel=[["(a)","(b)"],["(c)","(d)"],["(e)","(f)"]]
    fig,ax = CreateSuplotAxesSimple(1,2,10,20)
    # Subplot (0,0)
    ax[0].set_title("Hudson ice stream",fontweight="bold",fontsize=TitleFSize)
    PlotFiles =CreateTSInputVolume("Hudson")
    PlotPISMTimeSeries(PlotFiles,"thk",ax[0])
    ax[0].plot(T,np.ones(len(T))*5.07e15,color="grey", linestyle="dashed", \
        linewidth=3, label="Mean")
    ax[0].set_ylabel("Ice volume [m$^3$]",fontsize=AxisFSize)
    ax[0].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax[0].legend(loc=1,prop={'size':LegendFSize})
    ax[0].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[0].set_xlim(XLim)
    ax[0].text(LabelPos[0],LabelPos[1], PanelLabel[0][0],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[0].transAxes,fontsize=AxisFSize)
    AddArrowsVol(ax[0],"Hudson")

    ax[1].set_title("Mackenzie ice stream",fontweight="bold",fontsize=TitleFSize)
    PlotFiles =CreateTSInputVolume("Kenzie")
    PlotPISMTimeSeries(PlotFiles,"thk",ax[1])
    ax[1].plot(T,np.ones(len(T))*4.35e15,color="grey", linestyle="dashed", linewidth=3,label="Mean")
    ax[1].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    ax[1].legend(loc=1,prop={'size':LegendFSize})
    ax[1].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[1].set_xlim(XLim)
    ax[1].text(LabelPos[0],LabelPos[1], PanelLabel[0][1],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[1].transAxes,fontsize=AxisFSize)
    AddArrowsVol(ax[1],"Kenzie")
    plt.tight_layout()
    SavePlot("Fig09","IceVolumeLoss_IceStream")

if __name__ == '__main__':
    Main()
