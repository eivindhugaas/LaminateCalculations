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
a=10.0 #mm in radial direction is the full width of the crack or diameter of ellipse.
d=1.0 #mm in chardist
b=10.0 #mm in axial direction is the full width of the crack or diameter of ellipse.
angleofweakness=0. #degrees, defined from b to the axial direction
DamagedLayers=0

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
Finital=10.
Force=[[Finital],[0.],[0.],[0.],[0.],[0.]]
Width=50 #mm

#------------- Layup ----------------

# Number, Angle, Thickness, Material. Number 1 is bottom.

LayerData = np.array([ 
        [1,0.,0.32,1],
        [2,90.,0.32,1],
        [3,90.,0.32,1],
        [4,0.,0.32,1],
        ]) 

    #[E1     , E2    , G12  ,v12,v21    ,Xt   ,Xc  ,Yt ,Yc ,S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
MaterialProperties = np.array([
     [44800.0,12100.0,3400.0,0.3,-0.0879,1006.,487.,46.,132.,49.5,-0.5,           10.,   30.,0.,1.               ,10.,23.0103,0.,1.        ], #material 1 E1, E2, G12, v12,v21 Xt, Xc, Yt, Yc, S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
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

for change in changevector:
    a=b=change    
    Break=False
    F=Finital
    for i in range(0,10000):
        if Break:
            break
        F=F+1
        Force=[[F],[0.],[0.],[0.],[0.],[0.]]    
        StrengthVector=calc.LayupStrengthForceInput(layup=Layers,materials=MaterialProperties,Force=Force,Diameter=0.,DamagedLayers=DamagedLayers,a=a,b=b,angleofweakness=angleofweakness,d=d,ThickWalled=False)
        for l in StrengthVector['FeMaxStressVector']:
            if l[0]>1.:
                Break=True
                break
        
    TotalForce=(Width-change)*F
    te=(TotalForce/9.81)/1000.
    
    TotalForces.append(TotalForce)
    Totaltes.append(te)
    
    Strain=float(StrengthVector['StrainVector'][0])
    TotalStrain.append(Strain)
    
print(Totaltes)
plt.plot(changevector, Totaltes)
plt.show()


with open('text2.txt', 'w') as f:
    
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(changevector, Totaltes, TotalStrain, TotalForces))
quit()





