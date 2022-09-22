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
    Text = ["Ctrl", "SMB-", "Geo-","ST-", "SMB+", "Geo+", "ST+"]
    TextPos = [[6900, 1.06e15],[10300, 1.06e15], [8300, 1.20e15], [8100, 1.10e15], \
            [5300, 0.72e15], [4900, 0.75e15], [6000, 1.02e15]]
    TextPosKen = [[4500, 0.59e15],[7500, 0.38e15], [6000, 0.62e15], [5700, 0.84e15]]
    Colors = ["Black", "Blue", "Red", "Purple", "Blue", "Red", "Purple" ]
    Labels = ["Constant forcing (Ctrl)", "Surface mass balance (SMB)", "Geothermal heatflux (Geo)", "Ice surface temperature (ST)", "SMB", "Geothermal heatflux", "Ice surface temperature"]

    [DeltaVol_Hudson, HEFreq_Hudson] = ComputeSurgeTiming4IceVol("Hudson"\
            ,2.0e12, "anomaly")
    print("Now doing Mackenzie ...")
    [DeltaVol_Kenzie, HEFreq_Kenzie] = ComputeSurgeTiming4IceVol("Kenzie"\
            ,1.0e12, "anomaly")
    HEFreqMean_Hudson = ComputeMean(HEFreq_Hudson)
    HEFreqMean_Kenzie = ComputeMean(HEFreq_Kenzie)
    DeltaVolMean_Hudson = ComputeMean(DeltaVol_Hudson)
    DeltaVolMean_Kenzie = ComputeMean(DeltaVol_Kenzie)
    fig,ax = plt.subplots(1,2,figsize=(15,10))
    TitleFS=25
    AxisFS=22
    LabelFS=19
    TextFS=15

    for i,txt in enumerate(Text):
        if i < 4:
            ax[0].scatter(HEFreqMean_Hudson[i],DeltaVolMean_Hudson[i],s=450,marker="o",c=Colors[i], edgecolors=Colors[i], label=Labels[i])
            ax[1].scatter(HEFreqMean_Kenzie[i],DeltaVolMean_Kenzie[i],s=450,marker="o",c=Colors[i], edgecolors=Colors[i], label=Labels[i])
            ConfidenceEllipse(HEFreq_Kenzie[i],DeltaVol_Kenzie[i],ax[1],n_std=1.0,alpha=0.25,facecolor=Colors[i], edgecolor=Colors[i])
            ax[1].annotate(txt,(TextPosKen[i]),fontsize=TextFS)
        else:
            ax[0].scatter(HEFreqMean_Hudson[i],DeltaVolMean_Hudson[i],s=450,marker="o",c=Colors[i], edgecolors=Colors[i])
        ax[0].annotate(txt,(TextPos[i]),fontsize=TextFS)
        if len(HEFreq_Hudson[i]) > 4:
            ConfidenceEllipse(HEFreq_Hudson[i],DeltaVol_Hudson[i],ax[0],n_std=1.0,alpha=0.25,facecolor=Colors[i], edgecolor=Colors[i])
    ax[0].set_xlabel("Surge Period [yrs]",fontsize=AxisFS)
    ax[0].set_ylabel("$\Delta$ice volume [m$^3$]",fontsize=AxisFS)
    ax[0].tick_params(axis='both',labelsize=LabelFS)
    ax[0].set_title("Hudson ice stream",fontsize=TitleFS,
            fontweight='bold')
    ax[1].set_xlabel("Surge Period [yrs]",fontsize=AxisFS)
    ax[1].tick_params(axis='both',labelsize=LabelFS)
    ax[1].set_yticks([])
    ax[1].set_title("Mackenzie ice stream",fontsize=TitleFS,
            fontweight='bold')
    XLim = [3500,14200]
    YLim = [0.3e15, 1.45e15]
    ax[0].set_xlim(XLim)
    ax[0].set_ylim(YLim)
    ax[0].legend(scatterpoints=1,labelspacing=1.4, frameon=False,loc=3)
    ax[1].set_xlim(XLim)
    ax[1].set_ylim(YLim)
    ax[1].legend(scatterpoints=1,labelspacing=1.4, frameon=False,loc="best")
    ax[0].text(0.045, 0.965, '(a)',color="black", horizontalalignment='center', \
            verticalalignment='center', transform=ax[0].transAxes,fontsize=AxisFS)
    ax[1].text(0.045, 0.965, '(b)',color="black", horizontalalignment='center', \
            verticalalignment='center', transform=ax[1].transAxes,fontsize=AxisFS)
    plt.tight_layout()
    # plt.show()
    SavePlot("Fig08","IceVolume_PeriodPlotAnomaly",600)

if __name__ == '__main__':
    Main()
