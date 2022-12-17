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
    LegendFSize=24
    fig = plt.figure(figsize=(9.3,14))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    File1="Data/SMBMax_IcePerYearWithMask.nc"
    File2="Data/DeltaSMB_IcePerYearWithMask.nc"
    Var = "climatic_mass_balance"
    Extent=[0,10000,0,10000]
    PlotDir = [[File1, "seismic_r", "linear"]]
    PlotDir1 = [[File2, "seismic_r", "linear"]]
    CoastOpts=[["black","2"]]
    CBarOpts=[["-1000","1000","$\Delta$SMB [kg m$^{-2}$ yr$^{-1}$]","12","10"]]
    HudsonLines=CreateRectangleRegion("hudson")
    KenzieLines=CreateRectangleRegion("mackenzie")
    [_,_,im1]=Plot2DPISMField(Var,PlotDir,False,CoastOpts,False,False,False,ax1,Extent,30)
    [_,_,im2]=Plot2DPISMField(Var,PlotDir1,False,CoastOpts,False,False,False,ax2,Extent,15)
    im1.set_clim(-500,500)
    im2.set_clim(-50,50)
    cbar1 = fig.colorbar(im1, fraction=0.046, pad=0.04, ax=ax1)
    cbar1.set_label("SMB [kg m$^{-2}$ yr$^{-1}$]", fontsize=18)
    cbar1.ax.tick_params(labelsize=15)
    cbar2 = fig.colorbar(im2, fraction=0.046,pad=0.04, ax=ax2)
    cbar2.set_label("$\Delta$SMB [kg m$^{-2}$ yr$^{-1}$]", fontsize=18)
    cbar2.ax.tick_params(labelsize=15)
    ax1.set_xticks([])
    PlotRectangle(HudsonLines,ax1)
    PlotRectangle(KenzieLines,ax1)
    MakeAxesLabels([["","Distance [km]",
        "18","16"]],ax1)
    MakeAxesLabels([["Distance [km]","Distance [km]",
        "18","16"]],ax2)
    ax1.text(0.075, 0.95, '(a)',color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax1.transAxes,fontsize=LegendFSize)
    ax2.text(0.075, 0.95, '(b)',color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax2.transAxes,fontsize=LegendFSize)
    plt.tight_layout()
    SavePlot("SIFig01","SMB_SurgeCycle")

if __name__ == '__main__':
    Main()
