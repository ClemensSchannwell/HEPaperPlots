#! /usr/bin/env python
import sys
sys.path.append("../util")
import matplotlib.pyplot as plt
from matplotlib import colors
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
    YLim=[-1,6e12]
    YTickLabelSize=20
    LabelPos = [0.03, 0.975]
    PanelLabel=[["(a)","(b)"],["(c)","(d)"],["(e)","(f)"]]
    fig,ax = CreateSuplotAxesSimple(1,1,10,20)
    PlotFiles=CreateTSInputAnom("synchron","Dummy")
    FluxH, TimeH = ReadData("Hudson","HE12")
    FluxK, TimeK = ReadData("Kenzie","HE07")
    Start=200
    End=-1
    cmap = plt.cm.cool
    norm = colors.BoundaryNorm(np.arange(-45, 25, 5), cmap.N)
    sc = plt.scatter(FluxH[Start:End],
            FluxK[Start:End],marker='o',c=TimeH[Start:End]/1000,
            cmap=cmap,norm=norm,s=100)
    cbar=fig.colorbar(sc,ticks=np.arange(-45,25,5))
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Time [kyrs BP]',fontsize=18)
    ax.set_ylabel("Ice flux Mackenzie [m$^3$/yr]",fontsize=AxisFSize)
    ax.set_xlabel("Ice flux Hudson [m$^3$/yr]",fontsize=AxisFSize)
    ax.tick_params(axis='both',labelsize=YTickLabelSize)
    plt.savefig("PhasePlotIceFlux.png",bbox_inches='tight',dpi=300)

if __name__ == '__main__':
    Main()
