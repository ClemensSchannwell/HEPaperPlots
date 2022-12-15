#! /usr/bin/env python
import sys
sys.path.append("../util")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os
import errno
import copy

# Local imports
from UtilityFcts import *

def LoadPickle(filename):
    # Get global dictionary
    glob = globals()
    with open(filename, 'rb') as f:
        for k, v in pickle.load(f).items():
        # Set each global variable to the # value from the file
            glob[k] = v

def Main():
    fig,ax = CreateSuplotAxesSimple(1,2,8,20)
    LW=4
    LegendFSize=23
    AxisFSize=26
    TitleFSize=28
    YTickLabelSize=22
    File="Data/FlowLineData_Hudsonfine.nc"
    GradMat = ComputeGradient(File)
    t=np.arange(-65,-7.4,0.01)
    Start=1320
    End = 1410
    IceMask = copy.deepcopy(np.rot90(GradMat[:,Start:End]))
    IceMask[~np.isnan(IceMask)]=0
    IceMask[np.isnan(IceMask)]=1
    GradMat[np.isnan(GradMat)]=0 
    # left panel here
    cmapD = plt.cm.get_cmap('seismic', 15)
    im=ax[0].imshow(np.rot90(GradMat[:,Start:End]),vmin=-6e-3,vmax=6e-3,cmap=cmapD,aspect='auto',extent=[0,2500,0,900])
    t = ax[0].text(
                1100, 145, "Activation wave", ha="center", va="center", rotation=-10, size=25,
                    bbox=dict(boxstyle="larrow,pad=0.6", fc="red", ec="k", lw=2))
    ax[0].contour(np.flipud(IceMask),colors="gray",extent=[0,2500,0,900],linewidth=1)
    ax[0].set_xlabel("Distance [km]",fontsize=AxisFSize)
    ax[0].set_ylabel("Time [yrs]",fontsize=AxisFSize)
    cbar = fig.colorbar(im, fraction=0.066, pad=0.04, ax=ax[0])
    cbar.set_label('surface gradient [-]',fontsize=22)
    cbar.ax.tick_params(labelsize=20)
    ax[0].set_title("Ice surface gradient evolution",fontweight="bold",fontsize=TitleFSize)
    ax[0].text(0.15, 0.03, "upstream",color="red", horizontalalignment='center', verticalalignment='center',
                    transform=ax[0].transAxes,fontsize=LegendFSize, fontweight='bold')
    ax[0].text(0.77, 0.03, "downstream",color="red", horizontalalignment='center', verticalalignment='center',
                    transform=ax[0].transAxes,fontsize=LegendFSize, fontweight='bold')
    ax[0].text(2200, 220, "Front",color="magenta", horizontalalignment='center', verticalalignment='center',
                    fontsize=LegendFSize, fontweight='bold',backgroundcolor='0.95')
    ax[0].text(2200, 280, "B1",color="blue", horizontalalignment='center', verticalalignment='center',
                    fontsize=LegendFSize, fontweight='bold',backgroundcolor='0.95')
    ax[0].text(2200, 450, "B2",color="cyan", horizontalalignment='center', verticalalignment='center',
                    fontsize=LegendFSize, fontweight='bold',backgroundcolor='0.95')
    ax[0].text(2200, 600, "B3",color="black", horizontalalignment='center', verticalalignment='center',
                    fontsize=LegendFSize, fontweight='bold',backgroundcolor='0.95')
    ax[0].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[0].text(0.07,0.95, '(a)',color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[0].transAxes,fontsize=AxisFSize,backgroundcolor='1.00')

    # right panel here
    LoadPickle("Data/VelTS")
    ax[1].plot(Time,VelTrunk,color='magenta',linewidth=LW,label='Front')
    ax[1].plot(Time,VelB1, color='blue',linewidth=LW,label='B1')
    ax[1].plot(Time,VelB2, color='black',linewidth=LW,label='B2')
    ax[1].plot(Time,VelB3, color='cyan',linewidth=LW,label='B3')
    ax[1].set_ylabel("Ice velocity [m/yr]",fontsize=AxisFSize)
    ax[1].set_xlabel("Time [yrs]",fontsize=AxisFSize)
    ax[1].set_ylim([0,44000])
    ax[1].set_xlim([0,900])
    ax[1].tick_params(axis='both',labelsize=YTickLabelSize)
    ax[1].legend(loc="upper right",prop={'size':LegendFSize})
    ax[1].text(0.07,.95, '(b)',color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[1].transAxes,fontsize=AxisFSize)
    plt.tight_layout()
    # plt.show()
    SavePlot("Fig05","ActivationWaveHudson")

if __name__ == '__main__':
    Main()
