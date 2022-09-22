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
    YTickLabelSize=17
    LabelFS=20
    AxisFSize=18
    TitleFSize=27
    XLim=[0,85]
    YLim=[0,9e12]
    PanelLabel=[["(a)","(b)"],["(c)","(d)"],["(e)","(f)"]]
    LabelPos = [0.04, 0.92]
    LabPosX = [0.105, 0.120, 0.210]
    fig,ax = CreateSuplotAxesSimple(3,2,15,20)
    Vars=["Freq","climatic_mass_balance","dbdt"]
    PlotLabels=["", "SMB", "Geo. heatflux"]
    ForcingFreq=[0,6500],[3250,6500,13000],[0,0,0]
    for i in range(len(Vars)):
        [PlotFiles,ForcingFreq]=CreateTSInputFreq(Vars[i],"Hudson")
        PlotForcingHEFreq(ForcingFreq,ax[i][0],2.5)
        PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[i][0])
        ax[i][0].fill_between([20,88],[-100,-100],[9e15,9e15],
                color="grey", facecolor = 'grey',alpha=0.2)
        ax[i][0].set_ylabel("Ice flux [m$^3$/yr]",fontsize=AxisFSize)
        ax[i][0].set_xlim(XLim)
        ax[i][0].set_ylim(YLim)
        ax[i][0].text(LabPosX[i], 0.92, PlotLabels[i],color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax[i][0].transAxes,fontsize=LabelFS,backgroundcolor='1.0')
        ax[i][0].text(LabelPos[0],LabelPos[1], PanelLabel[i][0],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[i][0].transAxes,fontsize=LabelFS,backgroundcolor='1.0')
        [PlotFiles, _]=CreateTSInputFreq(Vars[i],"Kenzie")
        PlotForcingHEFreq(ForcingFreq,ax[i][1],2.5)
        PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[i][1])
        ax[i][1].fill_between([20,88],[-100,-100],[9e15,9e15],
                color="grey", facecolor = 'grey',alpha=0.2)
        ax[i][1].legend(loc=1,prop={'size':LegendFSize})
        ax[i][1].set_xlim(XLim)
        ax[i][1].set_ylim(YLim)
        ax[i][1].text(LabPosX[i], 0.92, PlotLabels[i],color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax[i][1].transAxes,fontsize=LabelFS,backgroundcolor='1.0')
        ax[i][1].text(LabelPos[0],LabelPos[1], PanelLabel[i][1],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[i][1].transAxes,fontsize=LabelFS,backgroundcolor='1.0')
        if i == 0:
            ax[i][0].set_title("Hudson ice stream",fontweight="bold",fontsize=TitleFSize)
            ax[i][1].set_title("Mackenzie ice stream",fontweight="bold",fontsize=TitleFSize)
        if i < len(Vars)-1:
            ax[i][1].set_yticklabels([])
            ax[i][1].set_xticklabels([])
            ax[i][0].set_xticklabels([])
            ax[i][0].tick_params(axis='both',labelsize=YTickLabelSize)
            ax[i][1].tick_params(axis='both',labelsize=YTickLabelSize)
        else:
            ax[i][1].set_yticklabels([])
            ax[i][0].tick_params(axis='both',labelsize=YTickLabelSize)
            ax[i][1].tick_params(axis='both',labelsize=YTickLabelSize)
            ax[i][0].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
            ax[i][1].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
    plt.tight_layout()
    SavePlot("Fig10","TSSurges_Freq")

if __name__ == '__main__':
    Main()
