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
    PanelLabel=["(a)","(b)"]
    LabelPos = [0.04, 1.03]
    fig1,ax1 = CreateSuplotAxesSimple(2,1,10,11)
    File1="Data/DeltaTasIceMask.nc"
    Var="tas"
    File2="Data/DeltaPrIceMask.nc"
    Var1="pr"
    CoastOpts=[["black","2"]]
    IceMaskOpts=[["gray","6","Ice extent"]]
    HudsonLines=CreateRectangleRegion("hudson")
    MackenzieLines=CreateRectangleRegion("mackenzie")
    Extent=[0,10000,0,10000]
    PlotDir = [[File1, "Reds", "linear"]]
    PlotDir1 = [[File2, "Blues", "linear"]]
    [_,_,im1]= Plot2DPISMField(Var,PlotDir,False,CoastOpts,False,IceMaskOpts,False,ax1[0],Extent)
    [_,_,im2]= Plot2DPISMField(Var1,PlotDir1,False,CoastOpts,False,IceMaskOpts,False,ax1[1],Extent)
    PlotRectangle(HudsonLines,ax1[0])
    PlotRectangle(HudsonLines,ax1[0])
    PlotRectangle(MackenzieLines,ax1[0])
    PlotRectangle(HudsonLines,ax1[1])
    PlotRectangle(HudsonLines,ax1[1])
    PlotRectangle(MackenzieLines,ax1[1])
    im1.set_clim([0,40])
    im2.set_clim([0,750])
    fig1.subplots_adjust(left=0.025)
    cbar_ax = fig1.add_axes([0.637, 0.53, 0.025, 0.35])
    cbar=fig1.colorbar(im1, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('$\Delta$T[K]',fontsize=18)
    fig1.subplots_adjust(left=0.025)
    cbar_ax1 = fig1.add_axes([0.637, 0.11, 0.025, 0.35])
    cbar1=fig1.colorbar(im2, cax=cbar_ax1)
    cbar1.ax.tick_params(labelsize=16)
    cbar1.set_label('$\Delta$Precipitation[K]',fontsize=18)
    MakeAxesLabels([["","Distance [km]",
        "18","16"]],ax1[0])
    ax1[0].tick_params(labelbottom=False)
    MakeAxesLabels([["Distance","Distance [km]",
        "18","16"]],ax1[1])
    ax1[0].text(LabelPos[0],LabelPos[1], PanelLabel[0],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax1[0].transAxes,fontsize=16)
    ax1[1].text(LabelPos[0],LabelPos[1], PanelLabel[1],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax1[1].transAxes,fontsize=16)
    SavePlot("SIFig05","PMIP4_EnsembleRange",600)

if __name__ == '__main__':
    Main()
