
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
pres=pres()
lamt=lamt()
scfs=scfs()
fcrit=fcrit()
geom=geom()
fatg=fatg()


#------------- Definition of weakness ------------------
a=40.0 #mm in radial direction or y
d=1.0 #mm in chardist
b=15. #mm in axial direction
angleofweakness=90.0 #degrees, defined from a to the axial direction
DamagedLayers=16

#------------- Definiton of pressure vessel ------------

Pressure=1. #N/mm^2
Diameter=500. #mm

#------------ No of increments --------

increments=100
fatigueanalysis=False
Residualburstanalysis=True
Residualburstvariable=a
VariableVector=np.arange(0.1,20.,0.5)

#------------- Layup ----------------

# Number, Angle, Thickness, Material. Number 1 is bottom.
LayerData = np.array([
        [28,55.0,1.,1],
        [27,-55.0,1.,1],
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
        [2,55.0,1.,1],  
        [1,-55.0,1.,1]
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
ResidualBurstVariableVector=[]
for variable in VariableVector:
    a=variable

    IntactLayers=[]
    IntactLayersInt=[]
    IntactLayers=Layers    
    FeMaxStress=0.
    FatigueLayers=[]
    Miner=MinerSum=np.matrix(np.zeros((len(Layers), 2)))
    FeMaxStressVector=[]
    FeMaxStressFailModeVector=[]
    feTsaiWuVector=[]
    TotalCycles=0.
    ResidualBurstVector=[]
    RelevantResidualBurstVector=[]
    for i in range(0,increments):
        if i>0:
            IntactLayers=IntactLayersInt
            IntactLayersInt=[]
            FatigueLayers=[]
            Miner=[]
            FeMaxStressVector=[]
            FeMaxStressFailModeVector=[]
            feTsaiWuVector=[]
            ResidualBurstVector=[]
            RelevantResidualBurstVector=[]
        for l in range(len(IntactLayers)):
            MaterialNo=int(IntactLayers[l][3])
            
            #------------------------------ Found material number ----------------------------------------------#
            
            G=lamt.Global_Stiffness_Matrix(LayerData=IntactLayers,Materials=Materials)
            
            #------------------------------ Established global stiffness matrix -----------------------------------#
            
            HydrostaticForce=[[pres.AxialLaminateForce(Pressure,Diameter)],[pres.HoopLaminateForce(Pressure,Diameter)]]
            
            #------------------------------ Found Hydrostatic force vector --------------------------------#
            
            G_hsplit=np.hsplit(G,3)
            
            G2by2=np.vsplit(G_hsplit[0],3)[0]
            
            Strain2by2=np.linalg.solve(G2by2,HydrostaticForce)
            
            Strain=np.vstack((Strain2by2,[[0.],[0.],[0.],[0.]]))
            
            Force=G*Strain
                
            #------------------------------ Found full force vector --------------------------------#
            
            MidPlane=lamt.LayerMidPlane(LayerData=IntactLayers,LayerNumber=l+1)
            
            WT=lamt.WallThickness(LayerData=IntactLayers)
            
            #-------------------------------- Found midplane position in layer and wall thickness ------------------------------#
            
            LayerStrain=lamt.LayerStrainfromGlobalStrain(GlobalStrain=Strain,Position=MidPlane)
            
            LocalStrain=lamt.Local_Strain_Transformation_Matrix(LayerData=IntactLayers[l])*LayerStrain
            
            LocalStress=lamt.Local_Stiffness_Matrix(MaterialProperties=Materials[MaterialNo-1])*LocalStrain
            
            #------------------- Established Local stress and strain ---------------------------------#        
            
            TWHSCF=pres.ThickWalledHoopSCF(Pressure=Pressure,Diameter=Diameter,WallThickness=WT,z=MidPlane)
        
            #---------------- Established thick walled cylinder stress amplification in Hoop direction ---------------# 
            
            LocalStressGlobalCoord=lamt.Local_Stress_Transformation_Matrix(LayerData=IntactLayers[l]).getI()*LocalStress
            
            LocalStressGlobalCoord=np.matrix([[float(LocalStressGlobalCoord[0])],[float(LocalStressGlobalCoord[1])*TWHSCF],[float(LocalStressGlobalCoord[2])]])
            
            LocalStress=lamt.Local_Stress_Transformation_Matrix(LayerData=IntactLayers[l])*LocalStressGlobalCoord
            
            #-------------- Established local amplified stress from thick walled distribution. -------------------#
            
            if l>=len(IntactLayers)-DamagedLayers:

                width=geom.WidthNormalToFiber(a=a,b=b,plyangle=IntactLayers[l][1],weaknessangle=angleofweakness)
                SCF=scfs.SCFAtXWhitneyNuismer(distancefromcenterofcrack=(width/2)+d,halfcracklength=width/2)
                LocalStress=np.matrix([[float(LocalStress[0])*SCF],[float(LocalStress[1])],[float(LocalStress[2])]])
        
            #---------------- Established Local Stress Concentration due to crack ------------------#
            
            feTsaiWu=fcrit.TsaiWuExpFact(Material=Materials[MaterialNo-1],LocalStressVector=LocalStress)
            
            FeMaxStress=fcrit.MaxStressExpFact(Material=Materials[MaterialNo-1],LocalStressVector=LocalStress)
            
            FeMaxStressFailMode=fcrit.MaxStressFailMode(Material=Materials[MaterialNo-1],LocalStressVector=LocalStress)
            
            FeMaxStressVector.append(FeMaxStress)
            
            ResidualBurstVector.append(((1./FeMaxStress)*Pressure))
            
            FeMaxStressFailModeVector.append(FeMaxStressFailMode)
            
            feTsaiWuVector.append(feTsaiWu)
            
            #------------------ Established exposure factors ---------------------------# 
            
            CyclesUntilFailure=fatg.CyclesUntilFailure(Material=Materials[MaterialNo-1],Stress=LocalStress,Thickness=IntactLayers[l][2])
                
            FatigueLayers.append(CyclesUntilFailure)
    
            #---------------- Established cycles until failure --------------------------#
        
        MinCycles=np.matrix(FatigueLayers).min()
        
        InvertedMiner=np.divide(FatigueLayers,MinCycles)
    
        Miner=np.divide(1.,InvertedMiner)
        
        MinerSum=MinerSum+Miner
        
        TotalCycles=TotalCycles+MinCycles
        
        #------------ Established Miner sum --------------------
        l=0
        failcount=0.
        RelevantResidualBurstVector=[]
        for l in range(len(Layers)):
            if float(IntactLayers[l][3])==4:
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(4)])
                
            #----------------- Checked if already destroyed, if not then: ---------------#       
            
            elif MinerSum[l,0]>=1.0 and MinerSum[l,1]>=1.0 and fatigueanalysis:
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(4)])
    
            #----------------- Check if destroyed by fatigue in both directions, if not then: ----------------#             
                
            elif float(IntactLayers[l][3])==2 and FeMaxStressVector[l]>1.0:
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(4)])
            elif float(IntactLayers[l][3])==3 and FeMaxStressVector[l]>1.0:
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(4)])
    
            #--------------- Check if already either destroyed by max stress in matrix or fiber, if not: ----------------#            
                     
            elif FeMaxStressVector[l]>1.0 and FeMaxStressFailModeVector[l]=='Xt' or FeMaxStressFailModeVector[l]=='Xc':
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(4)])
                
            #-------------- If fiber failure, then destroyed ------------------------------------------------------------#            
                
            elif FeMaxStressVector[l]>1.0 and FeMaxStressFailModeVector[l]=='Yt' or FeMaxStressFailModeVector[l]=='Yc' or FeMaxStressFailModeVector[l]=='S12':
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(3)])
                
            #-------------- If matrix failure, then compliant matrix. -------------------------------------------------#
                
            #--------------- Check if already either destroyed by max stress in matrix or fiber, if not: ----------------#  
            
            elif MinerSum[l,0]>=1.0 and fatigueanalysis:    
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(4)])         
            
            #-------------- If fiber failure, then destroyed ------------------------------------------------------------#
            
            elif MinerSum[l,1]>=1.0 and fatigueanalysis:
                IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(3)])
           
            #--------------- Check if destroyed by fatigue, if not: ----------------#  
            
            else:
                IntactLayersInt.append(list(IntactLayers[l]))
                
            if float(IntactLayersInt[l][3])!=4:
                RelevantResidualBurstVector.append(ResidualBurstVector[l])
                
        if RelevantResidualBurstVector!=[]:        
            ResidualBurst=min(RelevantResidualBurstVector)
        else:
            ResidualBurst=0.
        IntactLayersInt=np.array(IntactLayersInt) 
        resultmatrix=[]
        l=0
        failcount=0
        for l in range(len(IntactLayersInt)):
            if IntactLayersInt[l,3]==4:
                failcount=failcount+1
        #if not Residualburstanalysis:
        if failcount==len(IntactLayers) or np.array_equal(IntactLayersInt,IntactLayers) or TotalCycles>10**6:
            print("---------- Final increment is increment number",i,"-------------")
            print("---------- No more layers will fail , final layup: -------------")
            for l in range(len(IntactLayersInt)):
                resultmatrix.append((IntactLayersInt[l],FeMaxStressFailModeVector[l],FeMaxStressVector[l]))
            print(np.matrix(resultmatrix))
            print(MinerSum)
            print("---------- This layup can take",TotalCycles,"number of cycles---------")
            print("---------- at a residual burst pressure of",ResidualBurst,"Mpa---")
            break
        if not Residualburstanalysis:
            print("----------- Output increment",i," --------------")
        
            for l in range(len(IntactLayersInt)):
                resultmatrix.append((IntactLayersInt[l],FeMaxStressFailModeVector[l],FeMaxStressVector[l],ResidualBurstVector[l]))
            print(np.matrix(resultmatrix))
            print(MinerSum)
            print(TotalCycles)
    ResidualBurstVariableVector.append(ResidualBurst)        
print(VariableVector,ResidualBurstVariableVector)   
