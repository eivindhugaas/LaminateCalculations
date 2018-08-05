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
DamagedLayers=2

'''
Schematic shows a crack with fibers in axial, 0 deg, direction, large a, small b and 0 deg angleofweakness.
\-------------------------------------------------------------\
 \-------------------------------------------------------------\
  \-------------------------------------------------------------\
   )-------------------(<----a---->)|b---------------------------) <----------- axis of 0 deg of weakness
  /-------------------------------------------------------------/
 /-------------------------------------------------------------/
/-------------------------------------------------------------/
'''
#------------- Definiton of pressure vessel ------------

pressure=0.1 #N/mm^2
Diameter=500. #mm # Inner Diameter
Thickwalled=False #If true, then multiplies local stress in hoop direction with stress concentration factor to take into account higher stress towards center. Does not affect anything apart from the local hoop direction stress, as such it does not affect the compliance.
#------------ Definition of analysistype --------

increments=100
fatigueanalysis=False
residualburstanalysis=False
damageanalysis=True
ProgressiveFailureCriteria='Instantaneous'#'Puck'#'Instantaneous'
changevector=[0,1,2,3,4,5]


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

#LayerData = np.array([ 
        #[1,0.,1.,1],
        #[2,0.,1.,1],
        #[3,0.,1.,1],
        #[4,0.,1.,1],
        #[5,0.,1.,1],
        #[6,0.,1.,1],
        #[7,0.,1.,1],
        #[8,0.,1.,1],
        #[9,0.,1.,1],
        #[10,0.,1.,1]
        #])   
    #[E1     , E2    , G12  ,v12,v21    ,Xt   ,Xc  ,Yt ,Yc ,S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
MaterialProperties = np.array([
     [40000.0,11000.0,5000.0,0.3,-0.0879,1000.,500.,60.,30.,40.,-0.5,10.,30.,0.,1.,10.,23.0103,0.,1.], #material 1 E1, E2, G12, v12,v21 Xt, Xc, Yt, Yc, S12, fij, glassfiber: m, log a, k, tref (mm), matrix: m, log a, k, tref (mm).
     #[37.5370,11000.0,5.0,0.3,-0.0879,10000000000.,500000000000.0,200.0,100.0,1000000000000.0,-0.5,10.,30000.,0.,1.,10.,23.0103,0.,1.], #material 2 Fiber Failure
     #[37537.0,11.0,5.0,0.3,-0.0879,1000.,500.,20000000000.0,10000000000.,1000000000000.,-0.5,10.,30.,0.,1.,10.,23000.103,0.,1.],   #material 3 matrix failure
     #[37.5370,11.0,5.0,0.3,-0.0879,10000000000.,5000000000000.,2000000000000.0,100000000000.,100000000000.,-0.5,10.,30000.,0.,1.,10.,23000.103,0.,1.] #Full Failure
     ]) 

#-------------- Sorting Layers ---------------

Layers=[]
for i in range(len(LayerData)):
    for f in range(len(LayerData)):
        if i+1==int(LayerData[f][0]):
            Layers.append(list(LayerData[f]))
Layers=np.array(Layers)      

#------- Layers sorted for correct input to G --------
ResidualBurstVectorI=[]
ResidualBurstHistory=[]
HistoricForce=[]
HistoricStrain=[]
HistoricPressure=[]
HistoricResidualStrain=[]
FeMaxStressVectorFlatHistory=[]
PastPressure=1000000000000000000000000000.

for change in changevector:
    Pressure=pressure
    residualconverged=False
    if damageanalysis:
        a=b=change
    hitthewall=0
    for kappakappkappa in range(0,10000): 
        if residualconverged:
            break
        IntactLayers=Layers  
        Materials=MaterialProperties
        Miner=MinerSum=np.matrix(np.zeros((len(Layers), 2)))
        TotalCycles=0.
        ResidualBurstVector=[]
        convergencecount=0
        InstantaneousConverged=False
        toomuchstrain=False
        for i in range(0,10000):
            if i>0:
                IntactLayers=GraduallyDegradedLayup
                Miner=[]
            
            #---------------------------- Initiated iterative process ------------------------------  
            
            StrengthVector=calc.LayupStrength(layup=IntactLayers,materials=Materials,Pressure=Pressure,Diameter=Diameter,DamagedLayers=DamagedLayers,a=a,b=b,angleofweakness=angleofweakness,d=d,ThickWalled=Thickwalled)

            #--------------------------- Ran the Strength function ---------------------------------
            
            if ProgressiveFailureCriteria=='Puck':    
                FailureCriteriaVector=StrengthVector['FeHashinVector']
                ResidualBurstVector=StrengthVector['ResidualBurstVectorHashin']
                
            if ProgressiveFailureCriteria=='Instantaneous' or ProgressiveFailureCriteria=='Gradual': 
                FailureCriteriaVector=StrengthVector['FeMaxStressVector'] 
                ResidualBurstVector=StrengthVector['ResidualBurstVector']

            
            FeTsaiWuVector = StrengthVector['FeTsaiWuVector'] 
            FeMaxStressVector = StrengthVector['FeMaxStressVector'] 
            FeHashinVector=StrengthVector['FeHashinVector']
          
            #---------------------------- Established Strength and Fatigue Vectors from output of Strength function ------------------------------    
            
            FatigueLayers=StrengthVector['CyclesUntilFailureVector']
            MinCycles=np.matrix(FatigueLayers).min()
            InvertedMiner=np.divide(FatigueLayers,MinCycles)
            Miner=np.divide(1.,InvertedMiner)
            MinerSum=MinerSum+Miner
            TotalCycles=TotalCycles+MinCycles

            #------------ Established Miner sum --------------------
            
            DegradedMaterials=calc.DegradedMaterials(Layup=IntactLayers,materials=Materials,FailureCriteriaVector=FailureCriteriaVector,FatigueCriteriaVector=MinerSum,FailureCriteriaConstants=PuckConstants,ProgressiveFailureCriteria=ProgressiveFailureCriteria)
            
            GraduallyDegradedLayup=calc.GraduallyDegradedLayup(Layup=IntactLayers)
            
            for r in range(len(Materials)):
                if np.array_equal(Materials[r],DegradedMaterials[r]):
                    InstantaneousConverged=True   
                    
            Materials=DegradedMaterials
            
            StrengthVectorpost=calc.LayupStrength(layup=GraduallyDegradedLayup,materials=Materials,Pressure=Pressure,Diameter=Diameter,DamagedLayers=DamagedLayers,a=a,b=b,angleofweakness=angleofweakness,d=d,ThickWalled=Thickwalled)
            
            #------------ Established Degraded layup and its residual burst pressure ---------------------
            
            StrainVectorPost=StrengthVectorpost['StrainVector'] 
            StrainVector=StrengthVector['StrainVector']   
            
            PostStrain=StrainVectorPost[1]
            Strain=StrainVector[1]
            
            convergencecount=convergencecount+1
            if Strain>0.03:
                toomuchstrain=True
            
            if convergencecount>10:
                Pressure=Pressure-0.01
                hitthewall=hitthewall+1
                print("hit wall",hitthewall)

            Convergence=abs(PostStrain)-abs(Strain)
            print("Convergence",Convergence)
            if Convergence<0.0001 or hitthewall>5 or InstantaneousConverged or toomuchstrain:
                if min(ResidualBurstVector)<1000.:
                    ResidualBurst=min(ResidualBurstVector)
                if ProgressiveFailureCriteria=='Puck':
                    Materials=DegradedMaterials       
                       
                if ProgressiveFailureCriteria=='Instantaneous':
                    Materials=MaterialProperties
               
                if ProgressiveFailureCriteria=='Gradual':
                    Materials=DegradedMaterials       
                
                if not fatigueanalysis and not residualburstanalysis and not damageanalysis:
                    print("---------- Final increment is increment number",i,"-------------")
                    print("---------- No more layers will fail at a pressure of",Pressure,". Final layup: -------------")
                    print(DegreadedLayup)
                    print("---------- The residual burst pressure of the layup is",ResidualBurst,"Mpa---")
                    print("---------- The layup can take",TotalCycles,"number of cycles before any plies will begin to fail---------")
                if fatigueanalysis:
                    print("---------- The layup can take",TotalCycles,"number of cycles before failure of all plies at a pressure of",Pressure,"Mpa---------")
                if residualburstanalysis or damageanalysis:
                    
                    ForceVector=StrengthVectorpost['ForceVector']
                    StrainVector=StrengthVectorpost['StrainVector'] 
                    
                    HistoricForce.append(ForceVector[1,0])
                    HistoricStrain.append(StrainVector[1,0])
                    HistoricPressure.append(ResidualBurst)               
                    
                    ResidualBurstHistory.append(ResidualBurst)
                    Pressure=Pressure+(Pressure/100.)
                    
                    #print(Strain)

                    if Strain>0.025 or hitthewall>5:
                        residualconverged=True
                        if damageanalysis:
                            print("-----The residual burst pressure of the layup is",max(ResidualBurstHistory),"Mpa, at variable value",change,"-------")
                            print(ResidualBurstHistory)
                            
                            ResidualBurstVectorI.append(max(ResidualBurstHistory))
                            ResidualBurstHistory=[]
                            
                        else:
                            print("-----The residual burst pressure of the layup is",max(ResidualBurstHistory),"Mpa-------")
                            print(ResidualBurstHistory)
                            

                        
                break            


    
if residualburstanalysis:
    fig=plt.figure()
    fig.suptitle('Residualburst analysis graph',fontsize=14,fontweight='bold')
    
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    
    ax.set_title('Hoop Force vs Hoop Strain')
    
    ax.set_xlabel('Strain in Hoop direction')
    ax.set_ylabel('Force in Hoop Direction (N)')

    #ax.plot(HistoricStrain,HistoricForce,marker='o',color='b',linestyle='-')
    ax.plot(ResidualBurstVectorI,changevector,marker='o',color='b',linestyle='-')
    plt.show()
        
if damageanalysis:
    fig=plt.figure()
    fig.suptitle('Residualburst analysis graph',fontsize=14,fontweight='bold')
      
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
        
    ax.set_title('Hoop Force vs Hoop Strain')
        
    ax.set_xlabel('Change variable')
    ax.set_ylabel('Residual burst')
    
    #ax.plot(HistoricStrain,HistoricForce,marker='o',color='b',linestyle='-')
    ax.plot(changevector,ResidualBurstVectorI,marker='o',color='b',linestyle='-')
    plt.show()
