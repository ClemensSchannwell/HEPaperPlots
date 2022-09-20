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
    XLim=[0,85.0]
    YLim=[0,9e12]
    PanelLabel=[["(a)","(b)"],["(c)","(d)"],["(e)","(f)"], ["(g)","(h)"]]
    LabelPos = [0.04, 0.92]
    LabPosX = [0.105, 0.175, 0.15, 0.15]
    fig,ax = CreateSuplotAxesSimple(4,2,15,20)
    Vars=["climatic_mass_balance","bheatflx","ice_surface_temp","sea_level"]
    PlotLabels=["surface mass balance","geothermal heatflux", "ice surface temperature", "sea level"]
    PlotLabels=["SMB","Geo. heatflux", "Surf. temp.", "Sea level"]
    for i in range(len(Vars)):
        PlotFiles=CreateTSInputAnom(Vars[i],"Hudson")
        # if i <= 1:
        print(PlotFiles)
        PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[i][0])
        ax[i][0].set_ylabel("Ice flux [m$^3$/yr]",fontsize=AxisFSize)
        ax[i][0].set_xlim(XLim)
        ax[i][0].set_ylim(YLim)
        ax[i][0].text(LabPosX[i], 0.92, PlotLabels[i],color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax[i][0].transAxes,fontsize=LabelFS)
        ax[i][0].fill_between([20,88],[-100,-100],[9e12,9e12],
                color="grey", hatch = '//', label='Analysis Period',
                facecolor = 'none',alpha=0.6)
        ax[i][0].tick_params(axis='both',labelsize=YTickLabelSize)
        ax[i][0].text(LabelPos[0],LabelPos[1], PanelLabel[i][0],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[i][0].transAxes,fontsize=LabelFS)
        PlotFiles=CreateTSInputAnom(Vars[i],"Kenzie")
        # if i <= 1:
        PlotPISMTimeSeries(PlotFiles,"flux_crossSection",ax[i][1])
        ax[i][1].fill_between([20,88],[-100,-100],[9e12,9e12],
                color="grey", hatch = '//',
                facecolor = 'none',alpha=0.6)
        ax[i][1].legend(loc=1,prop={'size':LegendFSize})
        ax[i][1].set_xlim(XLim)
        ax[i][1].set_ylim(YLim)
        ax[i][1].text(LabPosX[i], 0.92, PlotLabels[i],color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax[i][1].transAxes,fontsize=LabelFS)
        ax[i][1].tick_params(axis='both',labelsize=YTickLabelSize)
        ax[i][1].text(LabelPos[0],LabelPos[1], PanelLabel[i][1],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[i][1].transAxes,fontsize=LabelFS)
        if i == 0:
            ax[i][0].set_title("Hudson ice stream",fontweight="bold",fontsize=TitleFSize)
            ax[i][1].set_title("Mackenzie ice stream",fontweight="bold",fontsize=TitleFSize)
        if i < len(Vars)-1:
            ax[i][1].set_yticklabels([])
            ax[i][1].set_xticklabels([])
            ax[i][0].set_xticklabels([])
        else:
            ax[i][1].set_yticklabels([])
            ax[i][0].set_xlabel("Time [kyrs]",fontsize=AxisFSize)
            ax[i][1].set_xlabel("Time [kyrs]",fontsize=AxisFSize)

    plt.tight_layout()
    # plt.show()
    SavePlot("Fig07","TSSurges_Anomaly")

if __name__ == '__main__':
    Main()
