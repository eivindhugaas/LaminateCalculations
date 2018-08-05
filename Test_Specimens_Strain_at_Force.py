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

PuckConstants=[4.,2.,0.5]

#------------- Definition of weakness ------------------
a=20.0 #mm in radial direction is the full width of the crack or diameter of ellipse.
d=1.0 #mm in chardist
b=20.0 #mm in axial direction is the full width of the crack or diameter of ellipse.
angleofweakness=0. #degrees, defined from b to the axial direction
DamagedLayers=8 #34 if mongolayup

'''
Schematic shows a crack with fibers in axial, 0 deg, direction, large a, small b and 0 deg angleofweakness.

Helical layers are in global 0 deg, that is in the axial. hoop is in global 90 deg.

\-------------------------------------------------------------\
 \-------------------------------------------------------------\
  \-------------------------------------------------------------\
   )-------------------(<----a---->)|b---------------------------) <----------- axis of 0 deg of weakness
  /-------------------------------------------------------------/
 /-------------------------------------------------------------/
/-------------------------------------------------------------/
'''

increments=100

fatigueanalysis=False #
residualburstanalysis=False #Iterates until it finds the highest pressure the laminate can whithstand. 
damageanalysis=False #Runs one residual burst per damage increment

ProgressiveFailureCriteria='Instantaneous'#'Puck'#'Instantaneous'. This applies for all. Progressive failure
changevector=[0.]#[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40] #For the damageanalysis, need specifying further down which variable the change applies to.
F=100.
Force=[[F],[0.],[0.],[0.],[0.],[0.]]
Width=50 #mm

#------------- Layup ----------------

# Number, Angle, Thickness, Material. Number 1 is bottom.

LayerData = np.array([ 
        [1,0.,1.35,1],
        [2,90.,1.35,1],
        [3,90.,1.35,1],
        [4,0.,1.35,1],
        ]) 

#FilamentLayers=30
#TotalThicknessFilalayer=0.732*2
#ThicknessFWM=TotalThicknessFilalayer/FilamentLayers

#LayerData = np.array([ 
        #[1,1.,((0.678+0.827)/2.)/2.,2],
        #[2,-1.,((0.678+0.827)/2.)/2.,2],
        #[3,(12.7-90.),ThicknessFWM,2],
        #[4,-(12.7-90.),ThicknessFWM,2],
        #[5,(12.7-90.),ThicknessFWM,2],
        #[6,-(12.7-90.),ThicknessFWM,2],
        #[7,(12.7-90.),ThicknessFWM,2],
        #[8,-(12.7-90.),ThicknessFWM,2],        
        #[9,(12.7-90.),ThicknessFWM,2],
        #[10,-(12.7-90.),ThicknessFWM,2],
        #[11,(12.7-90.),ThicknessFWM,2],
        #[12,-(12.7-90.),ThicknessFWM,2],
        #[13,(12.7-90.),ThicknessFWM,2],
        #[14,-(12.7-90.),ThicknessFWM,2],   
        #[15,(12.7-90.),ThicknessFWM,2],
        #[16,-(12.7-90.),ThicknessFWM,2],
        #[17,(12.7-90.),ThicknessFWM,2],
        #[18,-(12.7-90.),ThicknessFWM,2],
        #[19,(12.7-90.),ThicknessFWM,2],
        #[20,-(12.7-90.),ThicknessFWM,2],   
        #[21,(12.7-90.),ThicknessFWM,2],
        #[22,-(12.7-90.),ThicknessFWM,2],
        #[23,(12.7-90.),ThicknessFWM,2],
        #[24,-(12.7-90.),ThicknessFWM,2],
        #[25,(12.7-90.),ThicknessFWM,2],
        #[26,-(12.7-90.),ThicknessFWM,2],   
        #[27,(12.7-90.),ThicknessFWM,2],
        #[28,-(12.7-90.),ThicknessFWM,2],
        #[29,(12.7-90.),ThicknessFWM,2],
        #[30,-(12.7-90.),ThicknessFWM,2],
        #[31,(12.7-90.),ThicknessFWM,2],
        #[32,-(12.7-90.),ThicknessFWM,2],           
        #[33,1.,((0.678+0.827)/2.)/2.,2],
        #[34,-1.,((0.678+0.827)/2.)/2.,2],
#]) 

#LayerData = np.array([ 
        #[1,1.,((0.678+0.827)/2.)/2.,2],
        #[2,-1.,((0.678+0.827)/2.)/2.,2],
        #[3,(12.7-90.),(0.732*2)/4,2],
        #[4,-(12.7-90.),(0.732*2)/4,2],
        #[5,-(12.7-90.),(0.732*2)/4,2],
        #[6,(12.7-90.),(0.732*2)/4,2],
        #[7,-1.,((0.678+0.827)/2.)/2.,2],
        #[8,1.,((0.678+0.827)/2.)/2.,2],
#]) 

    #[E1     , E2    , G12  ,v12,v21    ,Xt   ,Xc  ,Yt ,Yc ,S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
MaterialProperties = np.array([
     [44800.0,12100.0,3400.0,0.3,-0.0879,1006.,487.,46.,132.,49.5,-0.5,           10.,   30.,0.,1.               ,10.,23.0103,0.,1.        ], #material 1 E1, E2, G12, v12,v21 Xt, Xc, Yt, Yc, S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
     [34092.0,11000.0,3070.0,0.3,-0.0879,755.15,365.65,39.,112.,42.,-0.5,           10.,   30.,0.,1.               ,10.,23.0103,0.,1.        ], #material 1 E1, E2, G12, v12,v21 Xt, Xc, Yt, Yc, S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
     ]) 

#-------------- Sorting Layers ---------------


Layers=[]
for i in range(len(LayerData)):
    for f in range(len(LayerData)):
        if i+1==int(LayerData[f][0]):
            Layers.append(list(LayerData[f]))
Layers=np.array(Layers)      

#------- Layers sorted for correct input to G --------

TotalForces=[]
Totaltes=[]
TotalStrain=[]
Break=False

for i in range(0,10000):
    if Break:
        break
    F=F+1
    Force=[[F],[0.],[0.],[0.],[0.],[0.]]    
    StrengthVector=calc.LayupStrengthForceInput(layup=Layers,materials=MaterialProperties,Force=Force,Diameter=0.,DamagedLayers=DamagedLayers,a=a,b=b,angleofweakness=angleofweakness,d=d,ThickWalled=False)
    if float(StrengthVector['StrainVector'][0])>0.01:
        Break=True
        break
    
TotalForce=(Width-a)*F
te=(TotalForce/9.81)/1000.

TotalForces.append(TotalForce)
Totaltes.append(te)

Strain=float(StrengthVector['StrainVector'][0])
TotalStrain.append(Strain)

print(Strain)
print(TotalForce)
print(StrengthVector['LocalStressGlobalCoordVector'])
print(StrengthVector['LocalStressVector'])
#StrengthVector = {'FeTsaiWuVector': FeTsaiWuVector, 'FeMaxStressVector': FeMaxStressVector,'FeHashinVector': FeHashinVector, 'FeMaxStressFailModeVector': FeMaxStressFailModeVector,'CyclesUntilFailureVector': CyclesUntilFailureVector,'ForceVector': Force,'StrainVector': Strain}


    




