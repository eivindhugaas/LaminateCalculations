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
DamagedLayers=5

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
#------------- Definiton of pressure vessel ------------
pressure=1.0 #N/mm^2
Diameter=500. #mm # Inner Diameter
Thickwalled=False #If true, then multiplies local stress in hoop direction with stress concentration factor to take into account higher stress towards center. Does not affect anything apart from the local hoop direction stress, as such it does not affect the compliance.
#------------ Definition of analysistype --------

increments=100

fatigueanalysis=False #
residualburstanalysis=True #Iterates until it finds the highest pressure the laminate can whithstand. 
damageanalysis=False #Runs one residual burst per damage increment

ProgressiveFailureCriteria='Instantaneous'#'Puck'#'Instantaneous'. This applies for all. Progressive failure
changevector=[0,1,2,3,4,5,6,7,8,9,10] #For the damageanalysis, need specifying further down which variable the change applies to.


#------------- Layup ----------------

# Number, Angle, Thickness, Material. Number 1 is bottom.

LayerData = np.array([ 
        [1,90.,1.,1],
        [2,90.,1.,1],
        [3,45.,1.,1],
        [4,-45.,1.,1],
        [5,0.,1.,1],
        [6,0.,1.,1],
        [7,-45.,1.,1],
        [8,45.,1.,1],
        [9,90.,1.,1],
        [10,90.,1.,1]
        ]) 

LayerData = np.array([ 
        [1,55.,1.,1],
        [2,-55.,1.,1],
        [3,55.,1.,1],
        [4,-55.,1.,1],
        [5,55.,1.,1],
        [6,-55.,1.,1],
        [7,55.,1.,1],
        [8,-55,1.,1],
        [9,55.,1.,1],
        [10,-55.,1.,1]
])   

    #[E1     , E2    , G12  ,v12,v21    ,Xt   ,Xc  ,Yt ,Yc ,S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
MaterialProperties = np.array([
     [40000.0,11000.0,5000.0,0.3,-0.0879,1000.,500.,60.,30.,40.,-0.5,           10.,   30.,0.,1.               ,10.,23.0103,0.,1.        ], #material 1 E1, E2, G12, v12,v21 Xt, Xc, Yt, Yc, S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
     ]) 

#-------------- Sorting Layers ---------------

Layers=[]
for i in range(len(LayerData)):
    for f in range(len(LayerData)):
        if i+1==int(LayerData[f][0]):
            Layers.append(list(LayerData[f]))
Layers=np.array(Layers)      

#------- Layers sorted for correct input to G --------


ResidualBurstHistory=[]

Break=False
#------- Defined Vectors ---------------

for change in changevector: #For every entry in the changevector, go  through this:
    Pressure=pressure #At start of every increment of damage propagation set pressure to originally applied pressure.
    residualconverged=False #Set that the residual has not converged at the start of every increment.
    if damageanalysis: #Set what the changevector should refer to.
        a=b=change #happy
    hitthewall=0 #Reset the hit the wall counter.
    StrainHistory=[]
    ForceHistory=[]
    PostStrain=0.
    ResidualBurstHistory=[]
    for kappakappkappa in range(0,10000): #For every stable iteration of layup.        
        if Break:
            break        
        if PostStrain>0.03:
            break        
        Layup=Layers #Set the layup to the damaged layup found through iterations.  
        Materials=MaterialProperties #Set the material to the material found through iterations.
        Miner=MinerSum=np.matrix(np.zeros((len(Layers), 2))) #Reset the miner sum.
        TotalCycles=0. #Reset cycle counter, will 
        Pressure=Pressure+0.1
        
        for i in range(0,10000): #i is iterations to find a converged solution.
            if i>0: #After first iteration, set the layup to the degreaded layup and reset the Miner for that layup.
                Miner=[]
            
            #---------------------------- Initiated iterative process ------------------------------  
            
            if i==0:
                StrengthVector=calc.LayupStrength(layup=Layup,materials=Materials,Pressure=Pressure,Diameter=Diameter,DamagedLayers=DamagedLayers,a=a,b=b,angleofweakness=angleofweakness,d=d,ThickWalled=Thickwalled)
            else:
                StrengthVector=StrengthVectorpost
              
            #--------------------------- Established Strength which outputs  ---------------------------------
            #{'FeTsaiWuVector','FeMaxStressVector','FeMaxStressFailModeVector','ResidualBurstVector','CyclesUntilFailureVector','ForceVector','StrainVector','FeHashinVector','ResidualBurstVectorHashin'}
            #--------------------------- The strength function outputs all the strength parameters of the layup------------
            
            if ProgressiveFailureCriteria=='Puck':    
                FailureCriteriaVector=StrengthVector['FeHashinVector']
                ResidualBurstVector=StrengthVector['ResidualBurstVectorHashin']
                
            if ProgressiveFailureCriteria=='Instantaneous' or ProgressiveFailureCriteria=='Gradual': 
                FailureCriteriaVector=StrengthVector['FeMaxStressVector'] 
                ResidualBurstVector=StrengthVector['ResidualBurstVector']
            
            #---------------------------- Established Strength and Fatigue Vectors from output of Strength function ------------------------------    
            
            FatigueLayers=StrengthVector['CyclesUntilFailureVector'] #Cycles left in each layer.
            MinCycles=np.matrix(FatigueLayers).min() #The layer with the least cycles left
            InvertedMiner=np.divide(FatigueLayers,MinCycles) 
            Miner=np.divide(1.,InvertedMiner)
            MinerSum=MinerSum+Miner
            TotalCycles=TotalCycles+MinCycles

            #------------ Established Miner sum --------------------
            
            Materials=calc.DegradedMaterials(Layup=Layup,materials=Materials,FailureCriteriaVector=FailureCriteriaVector,
                                                     FatigueCriteriaVector=MinerSum,FailureCriteriaConstants=PuckConstants,
                                                     ProgressiveFailureCriteria=ProgressiveFailureCriteria)
           
            #---------- Established the degraded material for each layer according to the strength functions output -----------
           
            Layup=calc.GraduallyDegradedLayup(Layup=Layup)
            
            #------ Spread correct degraded material out on correct layer, indeed very simple function ------------------
            
            StrengthVectorpost=calc.LayupStrength(layup=Layup,materials=Materials,Pressure=Pressure,Diameter=Diameter,DamagedLayers=DamagedLayers,a=a,b=b,angleofweakness=angleofweakness,d=d,ThickWalled=Thickwalled)
            
            #------------ Established the degraded materials strength -------------
            
            StrainVectorPost=StrengthVectorpost['StrainVector']
            ForceVectorPost=StrengthVectorpost['ForceVector']
            StrainVector=StrengthVector['StrainVector']   
                        
            #------------ Established vector containing the strain in global direction over whole layup. ----------------
            
            PostStrain=StrainVectorPost[1]
            Strain=StrainVector[1]
            PostForce=ForceVectorPost[1]
            
            #------------- Found the strain in Hoop direction for past iteration and this iteration. ----------
            
            Convergence=abs(abs(PostStrain)-abs(Strain))/abs(PostStrain)
            
            #------- Sets the absolute convergence ---------
            
            if Convergence<0.00001 or PostStrain:
                ResidualBurstHistory.append(min(ResidualBurstVector))
                StrainHistory.append(float(PostStrain))
                ForceHistory.append(float(PostForce))
                break
            
            if not residualburstanalysis and not fatigueanalysis and not damageanalysis:
                Break=True
                
    print(ResidualBurstHistory)
#print(ResidualBurstHistory[-1])
plt.plot(StrainHistory, ForceHistory)
plt.show()
