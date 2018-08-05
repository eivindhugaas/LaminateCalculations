
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import csv
import sys
from LaminateFunctions.LaminateFunctions import LaminateTheory as lamt
from LaminateFunctions.LaminateFunctions import StressConcentrationFunctions as scfs
from LaminateFunctions.LaminateFunctions import PressureFunctions as pres
from LaminateFunctions.LaminateFunctions import FailureCriteria as fcrit
from LaminateFunctions.LaminateFunctions import GeometricFunctions as geom
from LaminateFunctions.LaminateFunctions import FatigueFunctions as fatg
from LaminateFunctions.LaminateFunctions import CalculationFunctions as calc
pres=pres()
lamt=lamt()
scfs=scfs()
fcrit=fcrit()
geom=geom()
fatg=fatg()
calc=calc()


#------------- Definition of weakness ------------------
a=50.0 #mm in radial direction or y
d=1.0 #mm in chardist
b=20.0 #mm in axial direction
angleofweakness=90.0 #degrees, defined from a to the axial direction
DamagedLayers=0

#------------- Definiton of pressure vessel ------------

Pressure=62.826571620 #N/mm^2
Diameter=500. #mm

#------------ No of increments --------

increments=100
fatigueanalysis=False
Residualburstvariable=a

#------------- Layup ----------------

# Number, Angle, Thickness, Material. Number 1 is bottom.
LayerData = np.array([
        [28,0.,1.,1],
        [27,90.0,1.,1],
        [26,55.0,1.,1],  
        [25,-55.0,1.,1],
        [24,55.0,1.,1],
        [23,-55.0,1.,1],
        [22,55.0,1.,1],  
        [21,-55.0,1.,1],
        [20,55.0,1.,1],
        [19,-55.0,1.,1],
        [18,55.0,1.,1],  
        [17,-55.0,1.,1],
        [16,55.0,1.,1],
        [15,-55.0,1.,1],
        [14,55.0,1.,1],  
        [13,-55.0,1.,1],
        [12,55.0,1.,1],
        [11,-55.0,1.,1],
        [10,55.0,1.,1],  
        [9,-55.0,1.,1],
        [8,55.0,1.,1],
        [7,-55.0,1.,1],
        [6,55.0,1.,1],  
        [5,-55.0,1.,1],
        [4,55.0,1.,1],
        [3,-55.0,1.,1],        
        [2,90.0,1.,1],  
        [1,0.0,1.,1]
        ])   

Materials = np.array([
     [37537.0,11000.0,5000.0,0.3,0.0879,1000.,500.,200.,100.,1000.,-0.5,10.,30.,0.,1.,10.,23.0103,0.,1.], #material 1 E1, E2, G12, v12,v21 Xt, Xc, Yt, Yc, S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
     [37.5370,11000.0,5.0,0.3,0.0879,1000000.,5000000.0,200.0,100.0,1000000.0,-0.5,10.,30000.,0.,1.,10.,23.0103,0.,1.], #material 2 Fiber Failure
     [37537.0,11.0,5.0,0.3,0.0879,1000.,500.,200000.0,1000000.,1000.,-0.5,10.,30.,0.,1.,10.,23000.103,0.,1.],   #material 3 matrix failure
     [37.5370,11.0,5.0,0.3,0.0879,0.0001,0.0001,0.0001,0.0001,0.0001,-0.5,10.,30000.,0.,1.,10.,23000.103,0.,1.] #Full Failure
     ]) 

#-------------- Sorting Layers ---------------

Layers=[]
for i in range(len(LayerData)):
    for f in range(len(LayerData)):
        if i+1==int(LayerData[f][0]):
            Layers.append(list(LayerData[f]))
Layers=np.array(Layers)      

#------- Layers sorted for correct input to G --------

IntactLayers=Layers    
Miner=MinerSum=np.matrix(np.zeros((len(Layers), 2)))
TotalCycles=0.
ResidualBurstVector=[]

for i in range(0,increments):
    if i>0:
        IntactLayers=DegreadedLayup
        Miner=[]
     
    #---------------------------- Initiated iterative process ------------------------------  

    StrengthVector=calc.LayupStrength(layup=IntactLayers,materials=Materials,Pressure=Pressure,Diameter=Diameter,DamagedLayers=DamagedLayers,a=a,b=b,angleofweakness=angleofweakness,d=d)
        
    FeTsaiWuVector = StrengthVector['FeTsaiWuVector'] 
    FeMaxStressVector = StrengthVector['FeMaxStressVector'] 
    #FeMaxStressFailModeVector = StrengthVector['FeMaxStressFailModeVector']
    ResidualBurstVector=StrengthVector['ResidualBurstVector']
    FatigueLayers=StrengthVector['CyclesUntilFailureVector']
    
    #---------------------------- Established Strength and Fatigue for one layup ------------------------------    
    
    MinCycles=np.matrix(FatigueLayers).min()
    InvertedMiner=np.divide(FatigueLayers,MinCycles)
    Miner=np.divide(1.,InvertedMiner)
    MinerSum=MinerSum+Miner
    TotalCycles=TotalCycles+MinCycles
    ResidualBurst=min(ResidualBurstVector)
    
    #------------ Established Miner sum --------------------
    DegreadedLayup=calc.DegradedLayup(Layup=IntactLayers,MinerSum=MinerSum,FeMaxStressVector=FeMaxStressVector,fatigueanalysis=fatigueanalysis)
    
    resultmatrix=[]
    l=0
    failcount=0
    for l in range(len(DegreadedLayup)):
        if DegreadedLayup[l,3]==4:
            failcount=failcount+1
 
    if failcount==len(DegreadedLayup) or np.array_equal(DegreadedLayup,IntactLayers) or TotalCycles>10**6:
        if not fatigueanalysis:
            print("---------- Final increment is increment number",i,"-------------")
            print("---------- No more layers will fail, final layup and Max Stress exposure factor: -------------")
            for l in range(len(DegreadedLayup)):
                resultmatrix.append((DegreadedLayup[l],FeMaxStressVector[l]))
            print(np.matrix(resultmatrix))
            print(MinerSum)
            print("---------- The layup can take",TotalCycles,"number of cycles---------")
            print("---------- The residual burst pressure of the layup is",ResidualBurst,"Mpa---")
        if fatigueanalysis:
            print("---------- The layup can take",TotalCycles,"number of cycles before failure ---------")
        break
    print("----------- Output increment",i," --------------")

    for l in range(len(DegreadedLayup)):
        resultmatrix.append((DegreadedLayup[l],FeMaxStressVector[l],ResidualBurstVector[l]))
    print(np.matrix(resultmatrix))
    print(MinerSum)
    print(TotalCycles)
          
