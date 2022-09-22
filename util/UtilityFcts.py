#!/usr/bin/env python
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import numpy as np
import os
import errno
import pickle
# local import

def AddArrowsDS(ax):
    ax.text(0.14, 0.92, 'Hudson basin',color="blue",
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=16)
    ax.text(0.14, 0.65, 'Mackenzie basin',color="blue",
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=16)
    ax.annotate(s='', xy=(0.14,0.9), xytext=(0.14,0.67),
            arrowprops=dict(arrowstyle='<|-|>',mutation_scale=30,color="blue",
                lw=3),xycoords=ax.transAxes)

def AddArrowsVol(ax,region):
    if region == "Hudson": 
        ax.annotate(text='', xy=(0.08,0.63), xytext=(0.08,0.04),
                arrowprops=dict(arrowstyle='<|-|>',mutation_scale=30,color="blue",
                    lw=5),xycoords=ax.transAxes)
        ax.text(0.17, 0.02, "$\sim$3.8 m sea level",color="blue",
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes,fontsize=19)
    elif region == "Kenzie":
        ax.annotate(text='', xy=(0.14,0.54), xytext=(0.14,0.05),
                arrowprops=dict(arrowstyle='<|-|>',mutation_scale=30,color="blue",
                    lw=5),xycoords=ax.transAxes)
        ax.text(0.18, 0.03, "$\sim$1.8 m sea level",color="blue",
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes,fontsize=19)

def AddColourBar(argCbar,im,cb=True):
        # fmt='%1.0f'
        # cbar=plt.colorbar(im,format=fmt)
        im.set_clim(float(argCbar[0][0]),float(argCbar[0][1]))
        if cb:
            cbar=plt.colorbar(im)
            cbar.set_label(argCbar[0][2], fontsize=int(argCbar[0][3]))
            cbar.ax.tick_params(labelsize=int(argCbar[0][4]))

def AddLine(arg,ax):
    for f in arg:
        xval=np.array([float(f[0]),float(f[1])])
        yval=np.array([float(f[2]),float(f[3])])
        ax.plot(xval,yval,color=f[4],linewidth=f[5],linestyle=f[6])

def AddText(arg,ax):
    for f in arg:
        if f[7] == "on":
            props=dict(facecolor='white',edgecolor='none',pad=0.1)
            ax.text(f[1],f[2],f[0],ha='center', rotation=float(f[4]),fontweight=f[5],color=f[6],\
                va='center',transform=ax.transAxes,fontsize=int(f[3]),bbox=props)
        else:
            ax.text(float(f[1]),float(f[2]),f[0],ha='center', rotation=float(f[4]), \
                va='center',transform=ax.transAxes,fontsize=int(f[3]),fontweight=f[5],color=f[6])

def ComputeDeltaIVol(Ind, Region, Runs):
    [iVol,Time] = ReadData(Region,Runs,"ivol")
    DeltaIVol=[]
    for iEvent in range(len(Ind)-1):
        iVolMax = iVol[Ind[iEvent]]
        iVolMin = np.amin(iVol[Ind[iEvent]:Ind[iEvent+1]]) 
        DeltaIVol.append(iVolMax-iVolMin)
    return np.asarray(DeltaIVol)

def ComputeFlowLine(Region):
    if Region == "Hudson":
        InputFile="Data/DataFlowLineHudson.nc"
        xstart=317
        ystart=258
    elif Region == "Kenzie":
        InputFile="Data/DataFlowLineKenzie.nc"
        xstart=241
        ystart=433
    else:
        print("No valid region option. Choose from Hudson or Kenzie!")
        print("Exiting....")
        sys.exit()
    NCFile=Dataset(InputFile)
    X = NCFile.variables["x"][:]
    Y = NCFile.variables["y"][:]
    [xMat,yMat] = CreateXYMatFromVec(X,Y,len(X),len(Y))
    xVel = Netcdf2Numpy(InputFile,"uvelsurf","no")[0,:,:]
    xVel[xVel==-2e9]=np.nan
    yVel = Netcdf2Numpy(InputFile,"vvelsurf","no")[0,:,:]
    yVel[yVel==-2e9]=np.nan
    [coord,ind] = GenerateFlowLine(xMat[ystart,xstart],yMat[ystart,xstart],xMat,yMat,xVel,yVel,5000)
    return coord,ind,xVel,yVel,xMat,yMat

def ComputeGradient(File):
    NCFile=Dataset(File)
    Distance=NCFile.variables['y'][:]*1000
    Usurf=NCFile.variables['usurf'][:]
    Usurf[Usurf==0]=np.nan
    GradMat=np.zeros((231,len(Usurf[0,:])))
    print(Distance)
    for i in range(len(Usurf[0,:])):
        GradMat[:,i]=np.gradient(Usurf[:,i],Distance)
    return GradMat

    NCFile=Dataset(InputFile)
    X = NCFile.variables["x"][:]
    Y = NCFile.variables["y"][:]
    # Z = NCFile.variables['z'][:]
    [xMat,yMat] = CreateXYMatFromVec(X,Y,len(X),len(Y))
    xVel = Netcdf2Numpy(InputFile,"uvelsurf","no")[0,:,:]
    xVel[xVel==-2e9]=np.nan
    yVel = Netcdf2Numpy(InputFile,"vvelsurf","no")[0,:,:]
    yVel[yVel==-2e9]=np.nan
    [coord,ind] = GenerateFlowLine(xMat[ystart,xstart],yMat[ystart,xstart],xMat,yMat,xVel,yVel,5000)
    return coord,ind,xVel,yVel,xMat,yMat

def ComputeMean(InputTS):
    MeanVec = np.zeros(len(InputTS))
    for ii in range(0,len(InputTS)):
        MeanVec[ii] = np.mean(InputTS[ii])
    return MeanVec

def ComputeSurgeTiming4IceVol(Region,Thresh,RunType="anomaly"):
    HEAll = []
    FluxAll = []
    DeltaIVolAll = []
    Runs = GetRuns(Region,RunType)
    print(Runs)
    for iRun in range(len(Runs)):
        [Flux,Time] = ReadData(Region,Runs[iRun])
        IndStart = int(np.where(Time == -45000)[0][0])
        IndEnd = int(np.where(Time == 20000)[0][0])
        HEEvents=[]
        FluxHE = []
        for i in range(IndStart,IndEnd,21):
            FluxWind= Flux[i:i+21]
            TimeWind = Time[i:i+21]
            if any(FluxWind > Thresh):
                HEEvents.append(TimeWind[np.argmax(FluxWind)])
                FluxHE.append(np.max(FluxWind))
        [HEEvents,FluxHE] = DataCleaning(HEEvents,FluxHE)
        print(HEEvents)
        HEInds = GetIndices4HEEvent(HEEvents,Time)
        DeltaVol=ComputeDeltaIVol(HEInds,Region,Runs[iRun])
        HEFreq=np.diff(HEEvents)
        print(Runs[iRun])
        print(DeltaVol)
        print(HEFreq)
        HEAll.append(HEFreq)
        DeltaIVolAll.append(DeltaVol)
    return DeltaIVolAll, HEAll

def ConfidenceEllipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    if x.size != y.size:
        raise ValueError("x and y must be the same size")
    # sys.exit()
    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
        # Using a special case to obtain the eigenvalues of this
            # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor, **kwargs)
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

def CreateXYMatFromVec(xVec,yVec,MatDimX=400,MatDimY=400):
    xMat=np.tile(xVec, [MatDimY,1])
    yMat=np.tile(yVec, [MatDimX,1])
    yMatT=np.transpose(yMat)
    return xMat,yMatT


def CreateRectangleRegion(region):
    LineThick = "3.5"
    LineStyle="dashed"
    if region == "hudson":
        colour="purple"
        HudsonLine1=[["3120", "4690","2240","2240",colour, LineThick,
            LineStyle]]
        HudsonLine2=[["3120", "4690","3520","3520",colour, LineThick,
            LineStyle]]
        HudsonLine3=[["3120", "3120","2240","3520",colour, LineThick,
            LineStyle]]
        HudsonLine4=[["4690", "4690","2240","3520",colour, LineThick,
            LineStyle]]
        RecLine=[HudsonLine1,HudsonLine2,HudsonLine3,HudsonLine4]
    elif region == "mackenzie":
        colour ="cyan"
        MackenzieLine1=[["3500", "3500","3960","6000",colour, LineThick,
            LineStyle]]
        MackenzieLine2=[["1670", "1670","3960","6000",colour, LineThick,
            LineStyle]]
        MackenzieLine3=[["1670", "3500","3960","3960",colour, LineThick,
            LineStyle]]
        MackenzieLine4=[["1670", "3500","6000","6000",colour, LineThick,
            LineStyle ]]
        RecLine=[MackenzieLine1,MackenzieLine2,MackenzieLine3,MackenzieLine4]
    else:
        print("Unknown region specified. Choose between hudson and mackenzie")
        print("Exiting...")
        sys.exit()
    return RecLine

def CreateSuplotAxes(rows,cols,width=18,height=18):
    fig = plt.figure(figsize=(height,width))
    from mpl_toolkits.axes_grid1 import AxesGrid
    grid = AxesGrid(fig,111,nrows_ncols=(rows,cols),axes_pad=0.35,
            cbar_mode='single',cbar_location='right',cbar_pad=0.1)
    return fig,grid

def CreateSuplotAxesSimple(rows,cols,width=18,height=18):
    fig,ax = plt.subplots(rows,cols,figsize=(height,width))
    return fig,ax

def CreateTSInputFreq(VarName,Region):
    BasePath="/work/ba0989/m300792/HE_Runs/ExperimentsComposite/"
    Colours=["red","black","blue"]
    PlotFiles=[]
    if VarName == "Freq":
        Anomaly=["Ctrl (constant)","Cf (6500 years)"]
        Runs=["HE02","HE01"]
        ForcingFreq=[0,6500]
    elif VarName == "climatic_mass_balance":
        Anomaly=["CfSmb+ (3250 years)", "Cf (6500 years)","CfSmb- (13000 years)"] 
        Runs=["HE08","HE01","HE14"]
        ForcingFreq=[3250,6500,13000]
    elif VarName == "dbdt":
        Anomaly=["CfGeo+ (3250 years)", "Cf (6500 years)","CfGeo- (13000 years)"] 
        Runs=["HE23","HE01","HE24"]
        ForcingFreq=[3250,6500,13000]
    for i in range(len(Runs)):
        List=["Data/"+Region+"Hov_"+Runs[i]+".nc",Colours[i],Anomaly[i],"4.0"]
        PlotFiles.append(List)
    return PlotFiles,ForcingFreq

def CreateTSInputAnomExt(VarName,Region):
    Colours=["red","purple", "black"]
    PlotFiles=[]
    if VarName == "bheatflx":
        Anomaly=["Geo+ (+42 mW m$^{-2}$)", "Geo+ (+21 mW m$^{-2}$)","Ctrl"]
        if Region == "Hudson":
            Runs=["HE06", "HE96", "HE02"]
        elif Region == "Kenzie":
            Runs=["HE17","HE101", "HE02"]
    elif VarName == "climatic_mass_balance":
        Anomaly=["Smb+ (+100 kg m$^{-2}$yr$^{-1}$)","Smb+ (+50 kg m$^{-2}$yr$^{-1}$)", "Ctrl"] 
        if Region == "Hudson":
            Runs=["HE18", "HE94", "HE02"]
        elif Region == "Kenzie":
            Runs=["HE10", "HE99", "HE02"]
    elif VarName == "ice_surface_temp":
        Anomaly=["St+ (+5°C)","St+ (+2.5°C)", "Ctrl"] 
        if Region == "Hudson":
            Runs=["HE12","HE95", "HE02"]
        elif Region == "Kenzie":
            Runs=["HE22", "HE100", "HE02"]
    elif VarName == "reference":
        Runs=["HE02"]
        if Region == "Hudson":
            Anomaly=["Hudson"]
            Colours=["red"]
        elif Region == "Kenzie":
            Anomaly=["Mackenzie"]
            Colours=["black"]
    for i in range(len(Runs)):
        if VarName=="synchron":
            List=["Data/"+Region[i]+"Hov_"+Runs[i]+".nc",Colours[i],Anomaly[i],"4.0"]
        else:
            List=["Data/"+Region+"Hov_"+Runs[i]+".nc",Colours[i],Anomaly[i],"4.0"]
        PlotFiles.append(List)
    return PlotFiles

def CreateTSInputAnom(VarName,Region):
    Colours=["red","black","blue"]
    PlotFiles=[]
    if VarName == "bheatflx":
        Anomaly=["Geo+ (+42 mW m$^{-2}$)","Ctrl","Geo- (-42 mW m$^{-2}$)"]
        if Region == "Hudson":
            Runs=["HE06","HE02","HE13"]
        elif Region == "Kenzie":
            Runs=["HE17","HE02","HE07"]
    elif VarName == "climatic_mass_balance":
        Anomaly=["Smb+ (+100 kg m$^{-2}$yr$^{-1}$)", "Ctrl","Smb- (-100 kg m$^{-2}$yr$^{-1}$)"] 
        if Region == "Hudson":
            Runs=["HE18","HE02","HE19"]
        elif Region == "Kenzie":
            Runs=["HE10","HE02","HE11"]
    elif VarName == "ice_surface_temp":
        Anomaly=["St+ (+5°C)", "Ctrl","St- (-5°C)"] 
        if Region == "Hudson":
            Runs=["HE12","HE02","HE20"]
        elif Region == "Kenzie":
            Runs=["HE22","HE02","HE21"]
    elif VarName == "sea_level":
        Anomaly=["Sl+ (+10 m)", "Ctrl","Sl- (-10 m)"] 
        if Region == "Hudson":
            Runs=["HE45","HE02","HE44"]
        elif Region == "Kenzie":
            Runs=["HE45","HE02","HE44"]
    elif VarName == "oce":
        Anomaly=["Oc+ (+2°C)", "Ctrl"]
        if Region == "Hudson":
            Runs=["HE15","HE02"]
    elif VarName == "nosurge":
        if Region == "Hudson":
            Runs=["HE04","HE02"]
            Anomaly=["CtrlFrozMac", "Ctrl"]
        elif Region == "Kenzie":
            Anomaly=["CtrlFrozHud", "Ctrl"]
            Runs=["HE03","HE02"]
    elif VarName == "synchron":
            Runs=["HE12","HE07"]
            Region=["Hudson","Kenzie"]
            Anomaly=["Hudson (Smb+)","Mackenzie (Geo-)"]
    elif VarName == "reference":
        Runs=["HE02"]
        if Region == "Hudson":
            Anomaly=["Hudson"]
            Colours=["red"]
        elif Region == "Kenzie":
            Anomaly=["Mackenzie"]
            Colours=["black"]
    for i in range(len(Runs)):
        if VarName=="synchron":
            List=["Data/"+Region[i]+"Hov_"+Runs[i]+".nc",Colours[i],Anomaly[i],"4.0"]
        else:
            List=["Data/"+Region+"Hov_"+Runs[i]+".nc",Colours[i],Anomaly[i],"4.0"]
        PlotFiles.append(List)
    return PlotFiles

def CreateTSInputDivide(RunID,Region):
    BasePath="/work/ba0989/m300792/HE_Runs/ExperimentsComposite/"
    Colours=["red","black","blue"]
    PlotFiles=[]
    if RunID == "nosurge":
        if Region == "Hudson":
            Runs=["HE04","HE02"]
            Anomaly=["CtrlFrozMac", "Ctrl"]
        elif Region == "Kenzie":
            Anomaly=["CtrlFrozHud", "Ctrl"]
            Runs=["HE03","HE02"]
    for i in range(len(Runs)):
        List=["Data/DivideHudsonWest_"+Runs[i]+".nc",Colours[i],Anomaly[i],"4.0"]
        PlotFiles.append(List)
    return PlotFiles

def CreateTSInputVolume(Region):
    BasePath="/work/ba0989/m300792/HE_Runs/ExperimentsComposite/"
    Colours=["red","black","blue"]
    PlotFiles=[]
    if Region == "Hudson":
        Runs=["HE06","HE05"]
        Anomaly=["Geo+", "Ctrl"]
    elif Region == "Kenzie":
        Anomaly=["St+", "Ctrl"]
        Runs=["HE22","HE02"]
    for i in range(len(Runs)):
        List=[BasePath+Runs[i]+"/Postprocessing/IceVolume"+Region+".nc",Colours[i],Anomaly[i],"4.0"]
        List=["Data/IceVolume"+Region+"_"+Runs[i]+".nc",Colours[i],Anomaly[i],"4.0"]
        PlotFiles.append(List)
    return PlotFiles



def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

def DataCleaning(EventTS, FluxTS):
    Del = []
    for i in range(0,len(EventTS)-1):
        Delta=EventTS[i+1]-EventTS[i]
        if Delta < 700:
            Del.append(i+1)
    EventTS = np.delete(EventTS,Del)
    FluxTS  = np.delete(FluxTS,Del)
    return EventTS, FluxTS

def ExtrapolateFlowLine(Ind,Coord,xMat,yMat,Distance=30):
    XindStart = Ind[-1,0]
    XindEnd = Ind[-1,0]+round((Ind[-1,0]-Ind[-10,0])/10*Distance)
    YindStart = Ind[-1,1]
    YindEnd = Ind[-1,1]+round((Ind[-1,1]-Ind[-10,1])/10*Distance)
    indExtr = GenerateCrossSectionProfile(YindStart,XindStart,YindEnd,XindEnd)
    IndComb=np.concatenate((Ind,indExtr[1:-1,:]),axis=0)
    CoordComb=np.concatenate((Coord,Ind2XYCoordinates(indExtr[1:-1,:],yMat,xMat)),axis=0)
    return IndComb,CoordComb

def GenerateCrossSectionProfile(xstart,ystart,xend,yend):
    xPos = []
    yPos = []
    SamplingRate=1
    XPosGrid = np.zeros((15002,1))
    YPosGrid = np.zeros((15002,1))
    VecIndAll = np.zeros((15002,2))
    counter=0;
    xPos.append(xstart)
    yPos.append(ystart)
    DirectionVec = np.array([xend-xstart,yend-ystart])
    NormalDirection=DirectionVec/np.linalg.norm(DirectionVec)
    while (np.round(xPos[counter]) != xend or np.round(yPos[counter]) != yend):
        NewPosVec = [xPos[counter],yPos[counter]]+NormalDirection*SamplingRate
        counter=counter+1;
        xPos.append(NewPosVec[0])
        yPos.append(NewPosVec[1])
    XFinal=np.round(xPos)
    YFinal=np.round(yPos)
    IndArray=np.column_stack((YFinal,XFinal))
    IndexArrayTemp = RemoveDuplicates(IndArray)
    IndexArrayFinal=IndexArrayTemp [np.all(IndexArrayTemp !=0,axis=1)]
    return IndexArrayFinal.astype(int)

def GenerateFlowLine(xstart,ystart,XGridIn,YGridIn,vx,vy,SamplingRate=1000,Direction='forward'):
    xPos = []
    yPos = []
    XPosGrid = np.zeros((15002,1))
    YPosGrid = np.zeros((15002,1))
    VecIndAll = np.zeros((15002,2))
    counter=0;
    xPos.append(xstart)
    yPos.append(ystart)
    while counter < 15000:
        tmp=np.sqrt((yPos[counter]-YGridIn)**2+(xPos[counter]-XGridIn)**2)
        VecInd=np.transpose(np.asarray(np.where(tmp==tmp.min())))
        VecIndAll[counter,0]=VecInd[0,1];
        VecIndAll[counter,1]=VecInd[0,0];
        XPosGrid[counter]=(XGridIn[VecInd[0,0],VecInd[0,1]])
        YPosGrid[counter]=(YGridIn[VecInd[0,0],VecInd[0,1]])
        VxPoint=vx[VecInd[0,0],VecInd[0,1]]
        VyPoint=vy[VecInd[0,0],VecInd[0,1]]
        if np.isnan(VxPoint)== True or  np.isnan(VyPoint)== True:
            break;
        DirectionVec = [vx[VecInd[0,0],VecInd[0,1]],vy[VecInd[0,0],VecInd[0,1]]]
        Tmp = np.asarray(DirectionVec)
        if Direction.lower()=="Backward".lower():
            DirectionVec = Tmp*-1.0
        else:
            DirectionVec = Tmp
        NormalDirection=DirectionVec/np.linalg.norm(DirectionVec)
        NewPosVec = [xPos[counter],yPos[counter]]+NormalDirection*SamplingRate
        counter=counter+1;
        xPos.append(NewPosVec[0])
        yPos.append(NewPosVec[1])
    CoordArray=np.hstack((XPosGrid,YPosGrid))
    CoordArrayTemp = RemoveDuplicates(CoordArray)
    IndexArrayTemp = RemoveDuplicates(VecIndAll)
    IndexArrayFinal=IndexArrayTemp [np.all(IndexArrayTemp !=0,axis=1)]
    CoordArrayFinal=CoordArrayTemp [np.all(CoordArrayTemp !=0,axis=1)]
    return CoordArrayFinal,IndexArrayFinal.astype(int)

def GetFirstZeroInAllCols(Var):
    IndRow=(Var==0).argmax(axis=0)
    IndCol=np.arange(0,len(Var[0,:]),1)
    return IndRow,IndCol

def GetIndices4HEEvent(Events,Time):
    Inds = []
    for i in range(len(Events)):
        Inds.append(np.where(Time==Events[i])[0][0])
    return Inds

def GetRuns(Region,RunType="anomaly"):
    if Region == "Hudson":
        if RunType == "anomaly":
            Runs = ["HE02","HE19","HE13","HE20","HE18","HE06","HE12"]
        elif RunType == "anomalyint":
            Runs = ["HE02","HE18","HE06","HE12", "HE94", "HE96", "HE95"]
        else:
            Runs = ["HE01", "HE02","HE08", "HE14", "HE23", "HE24" ]
    else:
        if RunType == "anomaly":
            Runs = ["HE02","HE11","HE07","HE21","HE10","HE17","HE22"]
        else:
            Runs = ["HE01", "HE02","HE08", "HE14", "HE23", "HE24" ]
    return Runs

def Ind2XYCoordinates(IndArray,xMat,yMat):
    XCoord=xMat[IndArray[:,0],IndArray[:,1]]
    YCoord=yMat[IndArray[:,0],IndArray[:,1]]
    XYCoordArray=np.column_stack((XCoord,YCoord))
    return XYCoordArray


def MakeAxesLabels(arg,ax):
        ax.set_xlabel(arg[0][0],fontsize=int(arg[0][2]))
        ax.set_ylabel(arg[0][1],fontsize=int(arg[0][2]))
        ax.tick_params(axis='both',labelsize=int(arg[0][3]))

def MakeFigSurgeBehaviourComb(InputFile,Snapshots,Region,fig,ax3, GL=True):
    LW=2
    YLabelSize=17
    YTickLabelSize=14
    TitleSize=18
    if Region == "Mackenzie":
        counter=0
    else:
        counter=0
    ISName= Region + " ice stream"
    UPanelLabel=["(a)","(c)","(e)","(j)","(m)","(p)"]
    MPanelLabel=["(b)","(d)","(f)","(k)", "(n)","(q)"]
    LabelPos = [0.03, 1.09]
    # GL=True
    DummyY=np.array([-1e9,1e9])
    [Thk,GradThk,Zs,GradZs,Bed,Zb,TempBottom,VelBase,TillW,Bmelt,Dbdt
            ,Time,TauD,TauB,Distance]= ReadAllVars(InputFile)
    TitleStr=["Quiescent phase (" + str((Snapshots[0])/10) + " kyrs)",
            "Pre-surge phase (" + str((Snapshots[1])/10) + " kyrs)", 
            "Surge phase ("+ str((Snapshots[2])/10) + " kyrs)" ]
    for iSlice in range(counter,len(Snapshots)+counter):
        it = Snapshots[iSlice-counter]
        OceanStartInd=np.argmax(Zb[:,it]-Bed[:,it]>2)-1
        ax3[iSlice][0].set_title(TitleStr[iSlice-counter],fontsize=TitleSize,fontweight='bold')
        ax3[iSlice][0].fill_between(Distance,Zb[:,it],Zs[:,it], color='dimgrey')
        ax3[iSlice][0].fill_between(Distance[OceanStartInd:],Zb[OceanStartInd:,it],Bed[OceanStartInd:,it], color='blue')
        ax3[iSlice][0].fill_between(Distance,Bed[:,it],-5000*np.ones(len(Distance)), \
        color='saddlebrown')
        ax3[iSlice][0].set_xlim([0,2450])
        ax3[iSlice][0].set_ylim([-1000,3000])
        ax3[iSlice][0].text(0.40, 0.85, ISName,color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax3[iSlice][0].transAxes,fontsize=YTickLabelSize)
        ax3[iSlice][0].text(LabelPos[0],LabelPos[1], UPanelLabel[iSlice],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax3[iSlice][0].transAxes,fontsize=YLabelSize)
        ax3[iSlice][0].set_ylabel("Elevation [m]",fontsize=YLabelSize)
        ax3[iSlice][0].tick_params(axis='y',labelsize=YTickLabelSize)
        # # Second subplot
        ax2=ax3[iSlice][1].twinx()
        ax21=ax3[iSlice][1].twinx()
        # Set position of secondary y-axis
        ax21.spines['right'].set_position(('outward', 60))
        ax3[iSlice][1].plot(Distance,np.zeros(len(Distance)), color='grey',linestyle='dotted')
        ax3[iSlice][1].plot(Distance,TempBottom[:,it], color='black',linewidth=LW)
        p1, =ax2.plot(Distance,VelBase[:,it], color='blue',linewidth=LW)
        ax2.set_yscale('log')

        p2, =ax21.plot(Distance,TillW[:,it], color='purple',linewidth=LW)
        SetColourOfSecondYAxis(ax21,p2.get_color())
        ax21.set_ylim([0,2.05])
        ax21.set_ylabel("Till water [m]",fontsize=YLabelSize)
        ax21.tick_params(axis='y',labelsize=YTickLabelSize)

        ax3[iSlice][1].set_xlim([0,2450])
        ax3[iSlice][1].set_ylim([-6,0.5])
        ax2.set_ylim([1,65000])
        SetColourOfSecondYAxis(ax2,p1.get_color())
        ax3[iSlice][1].text(LabelPos[0], LabelPos[1], MPanelLabel[iSlice],color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax3[iSlice][1].transAxes,fontsize=YLabelSize)
        ax3[iSlice][1].set_ylabel(r"Ice bottom""\n"r"temperature [°C]",fontsize=YLabelSize)
        ax3[iSlice][1].tick_params(axis='y',labelsize=YTickLabelSize)
        ax2.set_ylabel("Ice Speed [m/yr]",fontsize=YLabelSize)
        ax2.tick_params(axis='y',labelsize=YTickLabelSize)
        if iSlice == 2:
            ax3[iSlice][0].set_xlabel("Distance [km]",fontsize=YLabelSize)
            ax3[iSlice][0].tick_params(axis='x',labelsize=YTickLabelSize)
            ax3[iSlice][1].set_xlabel("Distance [km]",fontsize=YLabelSize)
            ax3[iSlice][1].tick_params(axis='x',labelsize=YTickLabelSize)
        else:
            ax3[iSlice][0].set_xticklabels([])
            ax3[iSlice][1].set_xticklabels([])
        if GL:
            ax3[iSlice][0].plot([Distance[OceanStartInd-1],Distance[OceanStartInd-1]],DummyY,color="red",linewidth=LW)
            ax3[iSlice][1].plot([Distance[OceanStartInd-1],Distance[OceanStartInd-1]],DummyY,color="red",linewidth=LW)

    plt.tight_layout()
    return fig

def Make_Sure_Path_Exists(path):
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

def MaskByValue(Var2Mask,Value=0):
    Var2Mask[Var2Mask==Value] = np.nan
    return Var2Mask

def MaskByVar(Var2Mask,MaskVar,Value=0):
    Var2Mask[MaskVar==Value] = np.nan
    return Var2Mask

def Netcdf2Numpy(NCFile,VarName,latlon='no'):
    File=Dataset(NCFile)
    if latlon.lower() == 'no':
        Var=np.array(File.variables[VarName][:])
        return Var
    else:
        Var=np.array(File.variables[VarName][:])
        Lat=np.array(File.variables['lat'][:])
        Lon=np.array(File.variables['lon'][:])
        return Var,Lat,Lon

def PlotFlowLine(Region,Colour,LineWidth,ax,Res=10):
    [Coord,Ind,_,_,xMat,yMat] = ComputeFlowLine(Region)
    [IndAll,CoordAll] = ExtrapolateFlowLine(Ind,Coord,xMat,yMat)
    ax.plot(IndAll[:,0]*Res,IndAll[:,1]*Res,color=Colour,linewidth=LineWidth,\
            linestyle="dashed", dashes=(10, 3))

def PlotForcingHEFreq(ForcingFreq,ax,LineW=1.0,LineSty="dotted"):
    Colours=["red","black","blue"]
    for iRuns in range(len(ForcingFreq)):
        Fr = ForcingFreq[iRuns]/1000
        if Fr > 0:
            X=np.arange(0+Fr/2,57.6,Fr)
            XComb=np.column_stack((X,X))
            Y=np.ones((XComb.shape))
            Y[:,0]=Y[:,0]*0
            Y[:,1]=Y[:,1]*1e12
            ax.plot(XComb.T,Y.T,linewidth=LineW,color=Colours[iRuns],linestyle=LineSty)

def Plot2DPISMField(argVar,argPlot,argGL,argCoast,argGLExt,argIceMask,\
        argIceMaskExt,ax,Extent,level=15):
    NCFile=Dataset(argPlot[0][0])
    PlotVar = NCFile.variables[argVar][0,:,:]
    CS=[]
    labels=[]
    cmapD = plt.cm.get_cmap(argPlot[0][1], level)
    if argVar != "velsurf_mag":
        PlotVar[PlotVar==0]=np.nan
    if argPlot[0][2] == "log":
        im=ax.imshow(PlotVar,cmap=cmapD,origin="lower",extent=Extent,norm=LogNorm())
    else:
        # im=ax.imshow(PlotVar,cmap=argPlot[0][1],origin="lower",extent=Extent)
        im=ax.imshow(PlotVar,cmap=cmapD,origin="lower",extent=Extent)
    if argCoast:
        Coast = NCFile.variables['mask']
        Coast=Coast[0,:,:]
        Coast[Coast!=4]=0
        Coast[Coast==4]=1
        ax.contour(Coast,colors=argCoast[0][0],extent=Extent,linewidth=float(argCoast[0][1]))
    if argIceMask:
        IceMask = NCFile.variables['mask']
        IceMask=IceMask[0,:,:]
        IceMask[IceMask==4]=0
        IceMask[IceMask!=0]=1
        CS1=ax.contour(IceMask,colors=argIceMask[0][0],extent=Extent,linewidth=float(argIceMask[0][1]))
        if isinstance(CS,list):
            CS.append(CS1.collections[0])
            labels.append(argIceMask[0][2])
        else:
            CS=[CS, CS1.collections[0]]
            labels=[labels,argIceMask[0][2]]
    if argGL:
        GL = NCFile.variables['mask']
        GL=GL[0,:,:]
        GL[GL!=2]=0
        GL[GL==2]=1
        CS1=ax.contour(GL,colors=argGL[0][0],extent=Extent,linewidth=float(argGL[0][1]))
        CS=CS1.collections[0]
        labels=argGL[0][2]
    if argGLExt:
        NCFile=Dataset(argGLExt[0][0])
        GL = NCFile.variables['mask']
        GL=GL[0,:,:]
        GL[GL!=2]=0
        GL[GL==2]=1
        CS3=ax.contour(GL,colors=argGLExt[0][1],extent=Extent,linewidth=float(argGLExt[0][2]))
        CS=[CS, CS3.collections[0]]
        labels=[labels, argGLExt[0][3]]

    if argIceMaskExt:
        NCFile=Dataset(argIceMaskExt[0][0])
        IceMask = NCFile.variables['mask']
        IceMask=IceMask[0,:,:]
        IceMask[IceMask==4]=0
        IceMask[IceMask!=0]=1
        CS1=ax.contour(IceMask,colors=argIceMaskExt[0][1],extent=Extent,linewidth=float(argIceMaskExt[0][2]))
        if isinstance(CS,list):
            CS.append(CS1.collections[0])
            labels.append(argIceMaskExt[0][3])
        else:
            CS=[CS, CS1.collections[0]]
            labels=[labels,argIceMaskExt[0][3]]

    return CS,labels,im

def PlotPISMTimeSeries(arg,var,ax,shift=0,subset=None):
        counter=0
        for f in arg:
            File=f[0]
            Colour=f[1]
            Label=f[2]
            LineWidth=f[3]
            NCFile=Dataset(File)
            Var2Plot = np.asarray(NCFile.variables[var])
            if len(Var2Plot.shape) > 1:
                Var2Plot = np.asarray(NCFile.variables[var][:,0,0])
            Time = NCFile.variables['time'][:]/(365.0*86400*1000)
            if counter == 0:
                Time = Time - Time[0] - shift
                # counter=counter+1
            else:
                Time = Time - Time[0]
            if subset:
                start=int(subset[0][0])
                finish=int(subset[0][1])
                ax.plot(Time[start:finish],Var2Plot[start:finish],Colour,label=Label, linewidth=float(LineWidth));
            else:
                ax.plot(Time,Var2Plot, Colour,label=Label, linewidth=float(LineWidth));

def PlotRectangle(RectObj,ax):
    for iLine in range(len(RectObj)):
        AddLine(RectObj[iLine],ax)

def ProcessVars4Plotting(Zs,Zb,TempBottom,GradThk,GradZs,VelBase,Thk):
    IndRow,IndCol = GetFirstZeroInAllCols(Zs)
    Zs = MaskByValue(Zs)
    Zs = SetInd2Zero(Zs,IndRow,IndCol)
    Zb = SetInd2Zero(Zb,IndRow,IndCol)
    TempBottom= MaskByVar(TempBottom,Thk)
    GradThk= MaskByVar(GradThk,Thk)
    GradZs= MaskByVar(GradZs,Thk)
    VelBase = MaskByValue(VelBase,-2e9)
    return Zs,Zb,TempBottom,GradThk,GradZs,VelBase

def ReadAllVars(FilePath):
    NCFile=Dataset(FilePath)
    Distance = NCFile.variables['y'][:]
    Thk = NCFile.variables['thk'][:]
    GradThk = np.array(np.gradient(Thk,10000,axis=0))
    Zs= NCFile.variables['usurf'][:] 
    GradZs = np.array(np.gradient(Zs,10000,axis=0))
    Bed = NCFile.variables['topg'][:]
    Zb = Zs-Thk
    TempBottom = NCFile.variables['temp_paBott'][:]
    VelBase = NCFile.variables['velbase_mag'][:]
    TillW = NCFile.variables['tillwat'][:]
    Bmelt = NCFile.variables['bmelt'][:]
    Dbdt = NCFile.variables['dbdt'][:]
    Time = NCFile.variables['time'][:]/86400/365
    TauD = NCFile.variables['taud_mag'][:]
    TauB = NCFile.variables['taub_mag'][:]
    [Zs,Zb,TempBottom,GradThk,GradZs,VelBase]=ProcessVars4Plotting(Zs,Zb
            ,TempBottom,GradThk,GradZs,VelBase,Thk)
    return (Thk,GradThk,Zs,GradZs,Bed,Zb,TempBottom,VelBase,TillW,Bmelt,Dbdt
            ,Time,TauD,TauB,Distance)

def ReadData(Region,RunName, var="flux"):
    if var == "flux":
        File="Data/" + Region + "Hov_" + RunName + ".nc"
        NCFile=Dataset(File)
        Flux=NCFile.variables['flux_crossSection'][:]
        Time=NCFile.variables['time'][:]/(86400*365)
        return Flux,Time
    elif var =="ivol":
        File="Data/IceVolume" + Region + "_" + RunName + ".nc"
        NCFile=Dataset(File)
        IVol=NCFile.variables['thk'][:,0,0]
        Time=NCFile.variables['time'][:]/(86400*365)
        return IVol,Time

def RemoveDuplicates(XYInputArray):
    counter=0
    IndUnique = np.zeros((len(XYInputArray),len(XYInputArray[0,:])))
    for iFile in range(0,len(XYInputArray)):
        IndTemp = XYInputArray[iFile,:]
        if iFile==0:
           IndUnique[counter,:] = IndTemp
           counter=counter+1
        elif IndTemp[0]!= IndUnique[counter-1,0] or IndTemp[1]!=IndUnique[counter-1,1]:
           IndUnique[counter,:] = IndTemp
           counter=counter+1
    return IndUnique

def SavePlot(PlotNo,name,res=300):
    SaveStr=PlotNo + "_" + name + ".png"
    plt.savefig(SaveStr,bbox_inches='tight', dpi=int(res))

def SetColourOfSecondYAxis(ax,color):
    ax.yaxis.label.set_color(color)
    ax.spines["right"].set_edgecolor(color)
    ax.tick_params(axis='y', colors=color,which='both')

def SetInd2Zero(Var,IndRow,IndCol):
    Var[IndRow,IndCol]=0
    return Var

if __name__ == '__main__':
    Main()
