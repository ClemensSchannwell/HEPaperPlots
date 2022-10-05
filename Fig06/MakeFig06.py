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
    InputFile="Data/FlowLineData_Kenzie.nc"
    Snapshots=[265,277,292]
    fig,ax3 = CreateSuplotAxesSimple(3,2,10,17)
    figOut=MakeFigSurgeBehaviourComb(InputFile,Snapshots,"Mackenzie",fig,ax3,False)
    SavePlot("Fig06","SurgeBehaviour_Mackenzie")

if __name__ == '__main__':
    Main()
