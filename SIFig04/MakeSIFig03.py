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
    Text = ["Ctrl","SMB+", "Geo+", "ST+", "SMB+/2", "Geo+/2", "ST+/2"]
    TextPos = [[7100, 1.06e15],[5100, 0.70e15],[4700, 0.75e15], [6500, 1.01e15], \
            [5300, 0.76e15], [6100, 0.77e15], [6900, 0.98e15]]
    Colors = ["Black", "Blue", "Red", "Purple", "Blue", "Red", "Purple" ]
    Labels = ["Constant forcing (Ctrl)", "Surface mass balance (SMB)", "Geothermal heatflux (Geo)", "Ice surface temperature (ST)", "SMB", "Geothermal heatflux", "Ice surface temperature"]

    [DeltaVol_Hudson, HEFreq_Hudson] = ComputeSurgeTiming4IceVol("Hudson"\
            ,2.0e12, "anomalyint")
    HEFreqMean_Hudson = ComputeMean(HEFreq_Hudson)
    DeltaVolMean_Hudson = ComputeMean(DeltaVol_Hudson)
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    TitleFS=25
    AxisFS=22
    LabelFS=19
    TextFS=15

    for i,txt in enumerate(Text):
        if i < 4: 
            ax.scatter(HEFreqMean_Hudson[i],DeltaVolMean_Hudson[i],s=450,marker="o",c=Colors[i], edgecolors=Colors[i], label=Labels[i])
        else:
            ax.scatter(HEFreqMean_Hudson[i],DeltaVolMean_Hudson[i],s=450,marker="o",c=Colors[i], edgecolors=Colors[i])
        ax.annotate(txt,(TextPos[i]),fontsize=TextFS)
        if len(HEFreq_Hudson[i]) > 4 and i!=5:
                ConfidenceEllipse(HEFreq_Hudson[i],DeltaVol_Hudson[i],ax,n_std=1.0,alpha=0.25,facecolor=Colors[i], edgecolor=Colors[i])
    ax.set_xlabel("Surge Period [yrs]",fontsize=AxisFS)
    ax.set_ylabel("$\Delta$ice volume [m$^3$]",fontsize=AxisFS)
    ax.tick_params(axis='both',labelsize=LabelFS)
    ax.set_title("Hudson ice stream",fontsize=TitleFS,
            fontweight='bold')
    XLim = [3500,8100]
    YLim = [0.3e15, 1.25e15]
    ax.set_xlim(XLim)
    ax.set_ylim(YLim)
    ax.legend(scatterpoints=1,labelspacing=1.4, frameon=False,loc=3)
    plt.tight_layout()
    SavePlot("SIFig04","IceVolume_HudsonPosFine",600)

if __name__ == '__main__':
    Main()
