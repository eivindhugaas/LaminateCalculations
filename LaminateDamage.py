
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
pres=pres()
lamt=lamt()
scfs=scfs()
fcrit=fcrit()
geom=geom()

#------------- Definition of weakness ------------------
a=40.0 #mm in radial direction or y
d=1.0 #mm in chardist
b=20.0 #mm in axial direction
angleofweakness=90.0 #degrees, defined from a to the axial direction
DamagedLayers=6

#------------- Definiton of pressure vessel ------------

Pressure=15.4 #N/mm^2
Diameter=500 #mm

#------------- Layup ----------------

# Number, Angle, Thickness, Material. Number 1 is bottom.
LayerData = np.array([
        [4,90.0,0.3,1],
        [3,0.0,0.3,1],
        [2,0.0,0.3,1],  
        [1,90.0,0.3,1]
        ])   

Materials = np.array([
     [37537.0,11000.0,5000.0,0.3,0.0879,1000.,500.,200.,100.,1000.,-0.5,30.,10.,0.,1.], #material 1 E1, E2, G12, v12,v21 Xt, Xc, Yt, Yc, S12, fij, m, log a, k, tref (mm)
     [37.5370,11000.0,5.0,0.3,0.0879,1000000.,5000000.0,200.0,100.0,1000000.0,-0.5,30.,10.,0.,1.], #material 2 Fiber Failure
     [37537.0,11.0,5.0,0.3,0.0879,1000.,500.,200000.0,1000000.,1000.,-0.5,30.,10.,0.,1.]   #material 3 matrix failure
     ]) 

#-------------- Sorting Layers ---------------

Layers=[]
for i in range(len(LayerData)):
    for f in range(len(LayerData)):
        if i+1==int(LayerData[f][0]):
            Layers.append(list(LayerData[f]))
Layers=np.array(Layers)      

#------- Layers sorted for correct input to G --------

IntactLayers=[]
IntactLayersInt=[]
IntactLayers=Layers    
FeMaxStress=0.

for i in range(0,8):
    print("----------- increment:",i,"---------------") 
    
    if i>0:
        IntactLayers=IntactLayersInt
        IntactLayersInt=[]
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
        
        print (FeMaxStress)
        
        #------------------ Established exposure factors ---------------------------# 
        
        if FeMaxStress>1.0 and FeMaxStressFailMode=='Xt' or FeMaxStressFailMode=='Xc':
            IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(2)])
        elif FeMaxStress>1.0 and FeMaxStressFailMode=='Yt' or FeMaxStressFailMode=='Yc' or FeMaxStressFailMode=='S12':
            IntactLayersInt.append([int(IntactLayers[l][0]),float(IntactLayers[l][1]),float(IntactLayers[l][2]),int(3)])            
        else:
            IntactLayersInt.append(list(IntactLayers[l]))

    IntactLayersInt=np.array(IntactLayersInt) 
    print(IntactLayers)
    if np.array_equal(IntactLayersInt,IntactLayers):
        print("---------- No more layers will fail, final layup: -------------")
        print(IntactLayersInt)
        print("---------- Final increment is increment number",i,"-------------")
        break
