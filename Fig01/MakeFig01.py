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
    fig1,ax1 = CreateSuplotAxesSimple(1,1,10,11)
    File1="Data/HELikelyHoodAndIceMask.nc"
    Var="velbar_mag"
    CoastOpts=[["black","2"]]
    IceMaskOpts=[["gray","6","Ice extent"]]
    HudsonLines=CreateRectangleRegion("hudson")
    MackenzieLines=CreateRectangleRegion("mackenzie")
    DivideLine=[[["2670", "3090","4550","3560","blue", "3", "solid"]]]
    LegendFluxGate=[[["200", "700","9700","9700","red",
        "3","solid"]]]
    LegendFlowLine=[[["200", "700","9300","9300","plum",
        "3","solid"]]]
    LegendDivideLine=[[["200", "700","8900","8900","blue",
        "3","solid"]]]
    FluxGateLineKenzie=[[["2810", "3650","5500","5020","red",
        "3","solid"]]]
    FluxGateLineHudson=[[["4460", "4430","2780","2260","red",
        "3","solid"]]]
    MackenzieText=[["Mackenzie", "0.43", "0.61","18", "25", "normal",
        "black", "off"]]
    HudsonText=[["Hudson", "0.53", "0.20","18", "-35", "normal",
        "black","off"]]
    Branch1Text=[["B1", "0.28", "0.22","18", "35", "normal",
        "black","off"]]
    Branch2Text=[["B2", "0.27", "0.33","18", "0", "normal",
        "black","off"]]
    Branch3Text=[["B3", "0.34", "0.37","18", "0", "normal",
        "black","off"]]
    FluxGateText=[["Ice flux gate", "0.175", "0.97","18", "0", "normal",
        "red","off"]]
    FlowLineText=[["Flowline", "0.143", "0.93","18", "0", "normal",
        "plum","off"]]
    DivideLineText=[["Divide line", "0.16", "0.89","18", "0", "normal",
        "blue","off"]]
    CBarOpts=[["0","100","Surge likelyhood [%]","12","10"]]
    Extent=[0,10000,0,10000]
    PlotDir = [[File1, "inferno", "linear"]]
    [_,_,im1]= Plot2DPISMField(Var,PlotDir,False,CoastOpts,False,IceMaskOpts,False,ax1,Extent)
    AddColourBar(CBarOpts,im1,False)
    PlotRectangle(HudsonLines,ax1)
    PlotRectangle(MackenzieLines,ax1)
    AddText(MackenzieText,ax1)
    AddText(HudsonText,ax1)
    AddText(Branch1Text,ax1)
    AddText(Branch2Text,ax1)
    AddText(Branch3Text,ax1)
    AddText(FluxGateText,ax1)
    AddText(FlowLineText,ax1)
    AddText(DivideLineText,ax1)
    PlotRectangle(DivideLine,ax1)
    PlotRectangle(FluxGateLineHudson,ax1)
    PlotRectangle(FluxGateLineKenzie,ax1)
    PlotRectangle(LegendFluxGate,ax1)
    PlotRectangle(LegendFlowLine,ax1)
    PlotRectangle(LegendDivideLine,ax1)
    PlotFlowLine("Hudson","plum",3,ax1)
    PlotFlowLine("Kenzie","plum",3,ax1)
    fig1.subplots_adjust(left=0.025)
    cbar_ax = fig1.add_axes([0.877, 0.100, 0.025, 0.80])
    MakeAxesLabels([["Distance [km]","Distance [km]",
        "18","16"]],ax1)
    cbar=fig1.colorbar(im1, cax=cbar_ax)
    cbar.set_label('Surge likelyhood [%]',fontsize=18)
    cbar.ax.tick_params(labelsize=16)
    im1.set_clim([0,40])
    fig1.suptitle('Surge locations for the Laurentide ice sheet',fontsize=27,fontweight='bold') 
    SavePlot("Fig01","OverviewMap_SurgeLikelyhood",600)

if __name__ == '__main__':
    Main()
