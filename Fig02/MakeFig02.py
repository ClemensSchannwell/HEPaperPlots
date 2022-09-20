#!/usr/bin/env python

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import math
from matplotlib.patches import Rectangle


def AllocateVariables(WidthInInds,TotalNoOfEvents,weights,Var=None):
    MeanVar=np.zeros((int(TotalNoOfEvents),int(WidthInInds)))
    Out=np.zeros((len(weights),1))
    OutFinal=np.zeros((int(WidthInInds/2),1))
    return MeanVar,Out,OutFinal


def ComputeForcingAndMakePlot(Events,WindowWidth,dt,Time):

    TotalNoOfEvents = sum([len(element) for element in Events])
    counter=0
    Colors=CreateColors(TotalNoOfEvents)
    [WinCenter,WidthInInds,HalfWidthInInds,WindowInds] = ComputeWinPar(WindowWidth,dt)
    weights=np.cos(math.pi*WindowInds/(WidthInInds))**2
    [MeanVar,Out,OutFinal] = AllocateVariables(WidthInInds,TotalNoOfEvents,weights)

    Var= ReadNCVar()
    for iEvents in range(len(Events[0])):
        [Center,LC,RC]=ComputeWinBounds(Events[0][iEvents],Time,HalfWidthInInds)
        MeanVar[counter,:]=Var[LC:RC]
        counter=counter+1
    MeanVarInit=MeanVar
    MeanVar=np.mean(MeanVar,axis=0)
    [Out,OutFinal]=ComputeWeightedMean(WinCenter,HalfWidthInInds,MeanVar,Out,OutFinal,weights)
    CreatePlot(MeanVarInit,weights,WinCenter,HalfWidthInInds,OutFinal)

def ComputeWinBounds(Event,Time,HalfWidthInInds):
    Center=int(np.where(abs(Time-Event)<0.05)[0])
    LC = int(Center-HalfWidthInInds)
    RC = int(Center+HalfWidthInInds)
    return Center,LC,RC

def ComputeWeightedMean(WinCenter,HalfWidthInInds,MeanVar,Out,OutFinal,\
        weights,MeanVar1D=False,Out1D=False,OutFinal1D=False):
    IterInv=int(WinCenter+HalfWidthInInds)-1
    Iter = int(WinCenter)
    while Iter <= int(WinCenter + (HalfWidthInInds/2)):
        w = weights[Iter]
        wInv = 1 - weights[Iter]
        Out[Iter] = w * MeanVar[Iter] + wInv * MeanVar[IterInv]
        Iter = Iter + 1
        IterInv= IterInv - 1

    IterInv=int(WinCenter-HalfWidthInInds+1)
    Iter = int(WinCenter)

    while Iter >= int(WinCenter - (HalfWidthInInds/2)):
        w = weights[Iter]
        wInv = 1 - weights[Iter]
        Out[Iter] = w * MeanVar[Iter] + wInv * MeanVar[IterInv]
        Iter = Iter - 1
        IterInv= IterInv + 1
    Start=int(WinCenter-HalfWidthInInds/2)
    End=int(WinCenter+HalfWidthInInds/2)
    OutFinal=Out[Start:End]
    return Out,OutFinal

def ComputeWinPar(WindowWidth,dt):
    WidthInInds=int(WindowWidth/dt)*2
    WinCenter=int(WindowWidth/dt)
    HalfWidthInInds=int(WidthInInds/2)
    WindowInds=np.linspace(-HalfWidthInInds,HalfWidthInInds,WidthInInds+1)
    return WinCenter,WidthInInds,HalfWidthInInds,WindowInds

def CreateColors(TestPointH):
    if isinstance(TestPointH,int):
        Colors = plt.cm.gnuplot(np.linspace(0,1,TestPointH))
    else:
        Colors = plt.cm.gnuplot(np.linspace(0,1,len(TestPointH)))
    return Colors

def CreatePlot(MeanVar,weights,WinCenter,HalfWidthInInds,OutFinal):
    LineStyle=["dashed", "dotted", "dashdot"]
    Labels = ["Event 1", "Event 2", "Event3"]
    Color=["lightgray", "gray", "dimgrey"]
    LineW=4
    LegendFSize=22
    AxisFSize=22
    YTickLabelSize=18
    TitleFSize=26
    MeanForc=np.mean(MeanVar.T,axis=1)

    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(121)
    ax1b=ax1.twinx()
    ax2 = fig.add_subplot(122)
    Time=np.arange(0,0+len(MeanVar.T)*100,100)
    rec1 =Rectangle((3250,4.61e15),6500,1.45e15,linewidth=3,\
            edgecolor="purple", facecolor="none")
    rec2 =Rectangle((0,4.61e15),6500,1.45e15,linewidth=3,\
            edgecolor="purple", facecolor="none")
    rec3 =Rectangle((6500,4.61e15),6400,1.45e15,linewidth=3,\
            edgecolor="purple", facecolor="none")
    for iRun in range(len(MeanVar)):
        ax1.plot(Time,MeanVar[iRun].T,color=Color[iRun],linewidth=LineW,\
                label=Labels[iRun],linestyle=LineStyle[iRun])
    LabelPos = [0.05, 0.965]
    PanelLabel=[["(a)","(b)"],["(c)","(d)"],["(e)","(f)"]]
    ax1.text(LabelPos[0],LabelPos[1], "(a)",color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax1.transAxes,fontsize=28)
    ax2.text(LabelPos[0],LabelPos[1], "(b)",color="black",
            horizontalalignment='center', verticalalignment='center',
            transform=ax2.transAxes,fontsize=28)
    ax1.plot(Time,MeanForc,color="black",linewidth=4,label="Mean")
    p1, =ax1b.plot(Time,weights[0:-1],color="blue",linewidth=3,label="Weight function")
    SetColourOfSecondYAxis(ax1b,p1.get_color())
    ax1.scatter(Time[WinCenter],MeanForc[WinCenter],color='red',s=200)
    ax1.scatter(Time[WinCenter-HalfWidthInInds+1],MeanForc[WinCenter-HalfWidthInInds+1],color='red',s=200)
    ax1.legend(loc=1,prop={'size':LegendFSize})
    ax1.set_title("Heinrich events from ESM",fontsize=TitleFSize)
    ax1.set_ylabel("Ice volume [m$^3$]",fontsize=AxisFSize)
    ax1.set_xlabel("Time [yrs]",fontsize=AxisFSize)
    ax1b.set_ylabel("normalised weight",fontsize=AxisFSize)
    ax1.tick_params(axis='both',labelsize=YTickLabelSize)
    ax1b.tick_params(axis='both',labelsize=YTickLabelSize)
    ax1b.legend(loc=4,prop={'size':LegendFSize})
    ax1.add_patch(rec1)
    ax1.set_xlim([-150,13100])
    ax1.set_ylim([4.6e15,6.2e15])
    ax1b.set_ylim([-0.05,1.05])
    ax2.plot(Time[0:65],OutFinal,color="blue",linewidth=LineW, \
            label="Forcing cycle 1")
    ax2.plot(Time[65::],OutFinal,color="red",linewidth=LineW, \
            label="Forcing cycle 2")
    ax2.legend(loc=9,prop={'size':LegendFSize})
    ax2.set_xlim([-150,13100])
    ax2.set_ylim([4.6e15,6.2e15])
    ax2.set_xlabel("Time [yrs]",fontsize=AxisFSize)
    ax2.tick_params(axis='both',labelsize=YTickLabelSize)
    ax2.set_title("Composite cycle for ice volume",fontsize=TitleFSize)
    ax2.set_yticklabels([])
    ax2.add_patch(rec2)
    ax2.add_patch(rec3)
    plt.tight_layout()
    plt.savefig("Fig02.png",dpi=300)

def ReadNCVar():
    FileName="Data/volumeHudson.nc"
    NCFile=Dataset(FileName)
    Var = NCFile.variables["thk"][:,:]
    Var=Var[:,0,0]
    return Var

def SetColourOfSecondYAxis(ax,color):
    ax.yaxis.label.set_color(color)
    ax.spines["right"].set_edgecolor(color)
    ax.tick_params(axis='y', colors=color,which='both')

def main():
    ###############################################################################
    ### Parameters for script ####################################################
    ##############################################################################
    dt = 100
    Events=[[-28.3,-35.5,-41.6]]
    T = 6500
    Time = np.arange(-64.9,-22.6,0.1)

    ComputeForcingAndMakePlot(Events,T,dt,Time)

if __name__ == "__main__":
    main()
