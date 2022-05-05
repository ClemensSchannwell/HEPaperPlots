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
    InputFile="Data/FlowLineData_Hudson.nc"
    Snapshots=[189,206,212]
    fig,ax3 = CreateSuplotAxesSimple(3,3,10,17)
    figOut=MakeFigSurgeBehaviourComb(InputFile,Snapshots,"Hudson",fig,ax3,False)
    SavePlot("Fig04","SurgeBehaviour_Hudson")

if __name__ == '__main__':
    Main()
