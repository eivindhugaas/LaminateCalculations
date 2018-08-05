import numpy as np
import pandas as pd
import math
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import csv
import sys

class PressureFunctions:
    def __init__(self):
        Check="OK"
    
    def HoopLaminateForce(self,Pressure,Diameter):
        Force=Pressure*(Diameter*0.5)
        return Force  
    

    def AxialLaminateForce(self,Pressure,Diameter):
        Force=Pressure*(Diameter*0.5)*0.5
        return Force
    
    def ThickWalledHoopSCF(self,Pressure=0.,Diameter=0.,WallThickness=0.,z=0.):
        r1=Diameter/2.
        r2=r1+WallThickness
        r=(Diameter/2.)+(WallThickness/2.)+z        
        SCFThickWalled=(((r1**2)*((r**2)+(r2**2)))/(((r**2)*((r2**2)-(r1**2)))))
        SCFHoop=(Diameter/(2*WallThickness))
        SCF=SCFThickWalled/SCFHoop
        
        #This function uses the factor to multiply with pressure to get stress in hoop direction. It divides theese two factors and returns an SCF to multiply into the plys local stress. In other words it does not affect the global stiffness, only the local.
        
        return SCF
        
class StressConcentrationFunctions:
    def __init__(self):
        Check="OK"
        
    def SCFAtXWhitneyNuismer(self,distancefromcenterofcrack=0.,halfcracklength=0.):
        x=distancefromcenterofcrack
        a=halfcracklength
        SCFatx=(x/(((x**2)-(a**2))**0.5))
        return SCFatx
            
class LaminateTheory:
    def __init__(self):
        Check="OK"
        
    def WallThickness(self,LayerData=[[0.]]):
        H = 0
        H_Vector = []
        H_Sum = 0
        Htot = 0
        Htot_Vector = []
        Hk_Vector = []
        for i in range(1, (len(LayerData)+1)):
            H = LayerData[i-1,2]
            H_Vector.append(H)
        H_Sum = np.sum(H_Vector)
        
        #Checked 19.09.2017     
        
        return H_Sum
    
    def LayerStrainfromGlobalStrain(self,GlobalStrain=[0.],Position=0.):
        
        GlobalStrainTension=np.split(GlobalStrain, 2)[0]
        GlobalStrainBending=np.split(GlobalStrain, 2)[1]

        LayerStrain = GlobalStrainTension + (GlobalStrainBending*Position)  
        
        #Checked 19.09.2017    
        
        return LayerStrain            
    
    def Local_Stiffness_Matrix(self,MaterialProperties=[0.]):
        
        V21=-((MaterialProperties[3]*MaterialProperties[1])/MaterialProperties[0])
        
        Q11 = (MaterialProperties[0])/(1-(MaterialProperties[3]*V21))
        Q22 = (MaterialProperties[1])/(1-(MaterialProperties[3]*V21))
        Q12 = ((MaterialProperties[1])*(MaterialProperties[3]))/(1-(MaterialProperties[3]*V21))
        Q66 = MaterialProperties[2]
        LocalStiffnessMatrix = np.matrix(((Q11, Q12, 0),(Q12, Q22, 0),(0, 0, Q66)))           
            
        #Checked 19.09.2017    
        
        return LocalStiffnessMatrix

    def Local_Stress_Transformation_Matrix(self,LayerData=[0.]):       
        
        T_S_11_22 = (math.cos(math.radians(LayerData[1])))**2 
        T_S_12 = (math.sin(math.radians(LayerData[1])))**2 
        T_S_31 = 2*(math.cos(math.radians(LayerData[1])))*(math.sin(math.radians(LayerData[1])))
        T_S_32 = -T_S_31
        T_S_13 = T_S_32/2
        T_S_23 = T_S_31/2
        T_S_33 = T_S_11_22 - T_S_12

        TStressMatrix = np.matrix(((T_S_11_22, T_S_12, T_S_31),(T_S_12, T_S_11_22, T_S_32),(T_S_13, T_S_23, T_S_33)))

        #Checked 19.09.2017

        return TStressMatrix
        
    def Local_Strain_Transformation_Matrix(self,LayerData=[0.]):
                
        T_S_11_22 = (math.cos(math.radians(LayerData[1])))**2 
        T_S_12 = (math.sin(math.radians(LayerData[1])))**2 
        T_S_31 = 2*(math.cos(math.radians(LayerData[1])))*(math.sin(math.radians(LayerData[1])))
        T_S_32 = -T_S_31
        T_S_13 = T_S_32/2
        T_S_23 = T_S_31/2
        T_S_33 = T_S_11_22 - T_S_12

        TStrainMatrix = np.matrix(((T_S_11_22, T_S_12, T_S_31/2),(T_S_12, T_S_11_22, T_S_32/2),(T_S_13*2, T_S_23*2, T_S_33)))
        
        #Checked 19.09.2017
        
        return TStrainMatrix      
    
    def Global_Stiffness_Matrix(self,LayerData=[[0.]],Materials=[[0.]]):
        H = 0
        H_Vector = []
        H_Sum = 0
        Htot = 0
        Htot_Vector = []
        Hk_Vector = []
        
        for i in range(1, (len(LayerData)+1)):
            H = LayerData[i-1,2]
            H_Vector.append(H)
        H_Sum = np.sum(H_Vector)
        
        #------------------- Established total thickness and vector containing all thicknesses -------------
        
        for i in range(1, (len(LayerData)+1)):
            Htot = Htot + LayerData[i-1,2]
            Htot_Vector.append(Htot)
            
        #-------- Established vector containing the top position of each layer relative to bottom layer -------------
        
        for i in range(1, (len(LayerData)+1)):
            Hk = (-H_Sum/2) + Htot_Vector[i-1]
            Hk_Vector.append(Hk)
            
        #-------- Established vector containing the top position of each layer relative to midplane ------------------
        
        Hkmin1_Vector = np.subtract(Hk_Vector,H_Vector)
        
        #------ Established vector containing the bottom position of each layer relative to midplane ---------------
        
        Hkmin1_minHk_Vector = np.subtract(Hk_Vector,Hkmin1_Vector) 
        Hkmin1_minHk_scnd_Vector = np.subtract(np.power(Hk_Vector,2),np.power(Hkmin1_Vector,2))
        Hkmin1_minHk_third_Vector = np.subtract(np.power(Hk_Vector,3),np.power(Hkmin1_Vector,3))
        
        #-----------Established vectors containing the integration factors to multiply with the local stiffness matrices defined in the global stiffness matrix --------------
        
        T_Stress_I_mult_Local_Stiffness_Matrixes=[]
        T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes=[]
        
        A_i=B_i=D_i=[]
        A=B=D=[[0,0,0],[0,0,0],[0,0,0]]
        
        for i in range(len(LayerData)):
            MatNum=math.floor(LayerData[i][3]) # Finds the Layer's material
            Material=np.array(Materials[MatNum-1]) # Chooses correct material in material matrix
        
            Local_Stress_Transformation_Matrix=LaminateTheory().Local_Stress_Transformation_Matrix(LayerData=LayerData[i])
            
            #--------------- Established local stress transformation matrix -----------------
            
            Inverted_Local_Stress_Transformation_Matrix=LaminateTheory().Local_Stress_Transformation_Matrix(LayerData[i]).getI()
            
            #---------- Established local inverted stress transformation matrix -------------
            
            Local_Strain_Transformation_Matrix=LaminateTheory().Local_Strain_Transformation_Matrix(LayerData[i])
            
            #------------- Established local strain transformation matrix -------------------
            
            Inverted_Local_Strain_Transformation_Matrix=LaminateTheory().Local_Strain_Transformation_Matrix(LayerData[i]).getI()

            #------- Established local inverted strain transformation matrix ----------------

            T_Stress_I_mult_Local_Stiffness = Inverted_Local_Stress_Transformation_Matrix*LaminateTheory().Local_Stiffness_Matrix(MaterialProperties=Material)            
            T_Stress_I_mult_Local_Stiffness_mult_T_Strain = T_Stress_I_mult_Local_Stiffness*Local_Strain_Transformation_Matrix
            T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes.append(T_Stress_I_mult_Local_Stiffness_mult_T_Strain)
            
            #-------- Transformed Local stiffness matrix into the global coordinate system to get the local stiffness matrix in the global coordinate system --------
            
            
            A_i = (T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes[i])*((Hkmin1_minHk_Vector[i]))
            
            B_i = 0.5*(T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes[i])*((Hkmin1_minHk_scnd_Vector[i]))
            D_i = (1/3)*(T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes[i])*((Hkmin1_minHk_third_Vector[i]))
            A=A_i + A 
            B=B_i + B 
            D=D_i + D
            
            #--------------------- Established A, B and D matrices --------------------------
            
        E_L=np.concatenate((A, B), axis=0)
        
        #--------------- Established left side of stiffness matrix -------------------
        
        E_R=np.concatenate((B, D), axis=0)
        
        #--------------- Established right side of stiffness matrix -------------------
        
        E = np.hstack((E_L, E_R))
        
        #--------------- Established global stiffness matrix -------------------
        
        #Checked 19.09.2017 and 05.10.2017 benchmarked against Matrixes.xlsx.        
        
        return E  
    
    def LayerMidPlane(self,LayerData=[[0.]],LayerNumber=1000):
        LaminateThickness=0.
        for i in range(len(LayerData)):
            LaminateThickness=LaminateThickness+LayerData[i][2]
        BottomPlane=(-LaminateThickness/2)+(LayerData[0][2]/2)
        LayerMidPlane=BottomPlane
        LayerMidPlaneVector=[BottomPlane]
        for i in range((len(LayerData)-1)):
            LayerMidPlane=LayerMidPlane+(LayerData[i][2]/2)+(LayerData[i+1][2]/2)
            LayerMidPlaneVector.append(LayerMidPlane)
        
        # Checked 05.10.2017    

        return(LayerMidPlaneVector[LayerNumber-1])
    

class FailureCriteria:
    def __init__(self):
        Check="OK"
        
    def TsaiWuExpFact(self,Material=[0.],LocalStressVector=[[0.]]):
        s1=float(LocalStressVector[0])
        s2=float(LocalStressVector[1])
        t12=float(LocalStressVector[2])        
        Xt=Material[5]
        Xc=Material[6]
        Yt=Material[7]
        Yc=Material[8]
        S12=Material[9]
        fij=Material[10]
        F1=(1/Xt)-(1/Xc)
        F11=1/(Xt*Xc)
        F2=(1/Yt)-(1/Yc)
        F22=1/(Yt*Yc)
        #F3=(1/Zt)-(1/Zc)
        #F33=1/(Zt*Zc)
        #F44=1/(S23**2)
        #F55=1/(S13**2)
        F66=1/(S12**2)
        F12=fij*((F11*F22)**2)
        #F13=fij*((F11*F33)**2)
        #F23=fij*((F22*F33)**2)
        a=((F11*(s1**2))+(F22*(s2**2))+(F66*(t12**2)))+(2*F12*s1*s2)
        b=(F1*s1)+(F2*s2)
        c=-1.0
        Rup1=-1.*b
        Rup2a=(b**2.)-(4.0*a*c)
        Rup2=Rup2a**0.5
        Rup=Rup1+Rup2
        Rdown=2.*a
        R=Rup/Rdown
        fe=1./R
        return fe
    
    def MaxStressExpFact(self,Material=[0.],LocalStressVector=[[0.]]):
        s1=float(LocalStressVector[0])
        s2=float(LocalStressVector[1])
        t12=float(LocalStressVector[2]) 
        Xt=Material[5]
        Xc=Material[6]
        Yt=Material[7]
        Yc=Material[8]
        S12=Material[9]
        fij=Material[10]
        if s1>=0.:
            feX=s1/Xt
        elif s1<0.:
            feX=-s1/Xc
        if s2>=0.:
            feY=s2/Yt
        elif s2<0.:
            feY=-s2/Yc
        feS=abs(t12)/S12
        if feX<0. or feY<0. or feS<0.:
            print("Error in failure criteria")
        fe=max([feX,feY,feS])
        return fe
    
    def MaxStressExpVector(self,Material=[0.],LocalStressVector=[[0.]]):
        s1=float(LocalStressVector[0])
        s2=float(LocalStressVector[1])
        t12=float(LocalStressVector[2]) 
        Xt=Material[5]
        Xc=Material[6]
        Yt=Material[7]
        Yc=Material[8]
        S12=Material[9]
        fij=Material[10]

        if s1>=0.:
            feX=s1/Xt
        elif s1<0.:
            feX=-s1/Xc
        if s2>=0.:
            feY=s2/Yt
        elif s2<0.:
            feY=-s2/Yc
        feS=abs(t12)/S12
        if feX<0. or feY<0. or feS<0.:
            print("Error in failure criteria")
        fe=max([feX,feY,feS])
        return np.array([feX,feY,feS])    
        
    def MaxStressFailMode(self,Material=[0.],LocalStressVector=[[0.]]):
        s1=float(LocalStressVector[0])
        s2=float(LocalStressVector[1])
        t12=float(LocalStressVector[2]) 
        Xt=Material[5]
        Xc=Material[6]
        Yt=Material[7]
        Yc=Material[8]
        S12=Material[9]
        fij=Material[10]
        if s1>=0.:
            feX=s1/Xt
            failmodeX='Xt'
        elif s1<0.:
            feX=-s1/Xc
            failmodeX='Xc'
        if s2>=0.:
            feY=s2/Yt
            failmodeY='Yt'
        elif s2<0.:
            feY=-s2/Yc
            failmodeY='Yc'
        feS=abs(t12)/S12
        if feX<0. or feY<0. or feS<0.:
            print("Error in failure criteria")
        fe=max([feX,feY,feS])
        if fe==feX:
            failmode=failmodeX
        if fe==feY:
            failmode=failmodeY
        if fe==feS:
            failmode='S12'
        return failmode        
    def PuckExpVector(self,Material=[0.],LocalStressVector=[[0.]]):
        p1tmin=0.25 #GFRP=0.25 CFRP=0.3
        p1tpluss=0.3 #GFRP=0.3 CFRP=0.35
        s1=float(LocalStressVector[0])
        s2=float(LocalStressVector[1])
        t12=float(LocalStressVector[2]) 
        Xt=Material[5]
        Xc=Material[6]
        Yt=Material[7]
        Yc=Material[8]
        S12=Material[9]
        fij=Material[10]
        
        RAtt=(S12/(2*p1tmin))*(((1+(2*p1tmin*(Yc/S12)))**0.5)-1.)
        pttmin=RAtt*(p1tpluss/S12)
        if s1>=0.:
            feX=(((s1/Xt)**2)+((t12/S12)**2))**0.5
        elif s1<0.:
            feX=((s1/Xc)**2)**0.5
        if s2>=0.:
            feY=(((s2/Yt)**2)+((t12/S12)**2))**0.5
        elif s2<0.:
            a=((s2/(2*S12))**2)+((t12/S12)*2)
            b=((s2/Yc)*(((Yc/(2*S12))**2)-1))
            c=-1
            feY=((-b)+(((b**2)-(4*a*c))**0.5))/(2*a)
        if feX<0. or feY<0.:
            print("Error in failure criteria")
        fe=max([feX,feY])
        return np.array([feX,feY])
    
    def HashinExpVector(self,Material=[0.],LocalStressVector=[[0.]]):
        s1=float(LocalStressVector[0])
        s2=float(LocalStressVector[1])
        t12=float(LocalStressVector[2]) 
        Xt=Material[5]
        Xc=Material[6]
        Yt=Material[7]
        Yc=Material[8]
        S12=Material[9]
        fij=Material[10]
        if s1>=0.:
            feX=(((s1/Xt)**2)+((t12/S12)**2))**0.5
        elif s1<0.:
            feX=((s1/Xc)**2)**0.5
        if s2>=0.:
            feY=(((s2/Yt)**2)+((t12/S12)**2))**0.5
        elif s2<0.:
            a=((s2/(2*S12))**2)+((t12/S12)*2)
            b=((s2/Yc)*(((Yc/(2*S12))**2)-1))
            c=-1
            feY=((-b)+(((b**2)-(4*a*c))**0.5))/(2*a)
        if feX<0. or feY<0.:
            print("Error in failure criteria")
        fe=max([feX,feY])
        return np.array([feX,feY])

class GeometricFunctions:
    def __init__(self):
        Check="OK"
        
    def WidthNormalToFiber(self,a=0.,b=0.,plyangle=0.,weaknessangle=0.):
        plyangle=plyangle+90
        weaknessangle=weaknessangle+90
        x=b*math.cos(math.radians(plyangle+weaknessangle))
        y=a*math.sin(math.radians(plyangle+weaknessangle))
        width=((x**2+y**2)**0.5)

        
        #Checked 05.10.2017 90 is added to avoid negative degrees in case this messes with python.
        return width

class FatigueFunctions:
    def __init__(self):
        Check="OK"
    
    def CyclesUntilFailure(self,Material=[0.],Stress=[0.],Thickness=0.):
        mf=Material[11]
        logaf=Material[12]
        kf=Material[13]
        treff=Material[14]
        mm=Material[15]
        logam=Material[16]
        km=Material[17]
        trefm=Material[18]
        
        deltaf=float(Stress[0])
        deltam=float(Stress[1])
        
        
        t=Thickness
        
        if Material[12]>100. or Material[16]>100.:
            Nf=Nm=10**30
        else:
            Nf=(10.**(logaf))/((deltaf*((t/treff)**kf))**mf)
        
            Nm=(10.**(logam))/((deltam*((t/trefm)**km))**mm)
                
        return [Nf,Nm]
    
class ProgressiveFailureFunctions:
       
    def __init__(self):
        Check="OK"           
        
    def PuckStiffnessReductionFactor(self,FailureCriteria=[],FatigueCriteria=[],PuckConstants=[]):
        fe=FailureCriteria[1]
        c=PuckConstants[0]
        e=PuckConstants[1]
        nr=PuckConstants[2]
        fe=FailureCriteria[1]
        FatigueCriteria=[0.1,0.1]
        
        n=1.0 #Matrix E
        s=1.0 #Fiber E
        g=1.0 #Shear E
        x=1.0 #Strength Fiber
        y=1.0 #Strength Matrix
        xy=1.0 #Strength Shea        
        
        if FailureCriteria[1]>=1.0: #Matrix Failure
            n=((1-nr)/(1+(c*((fe-1)**e))))+nr 
            #print(n)
        if FailureCriteria[0]>1.0: #Fiber Failure
            s=0.1
            x=100000.
        #print(FailureCriteria[1])
        #print(s,n)
        nvector=[s,n,1.,1.,1.,x,x,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]   
        return nvector
    
    def ConstantStiffnessReductionFactor(self,FailureCriteria=[],FatigueCriteria=[],ConstantConstants=[]):
        fe=FailureCriteria[1]
        c=ConstantConstants[0]
        e=ConstantConstants[1]
        nr=ConstantConstants[2]
        fe=FailureCriteria[1]
        FatigueCriteria=[0.1,0.1]
        n=1.0 #Matrix E
        s=1.0 #Fiber E
        g=1.0 #Shear E
        xt=1.0
        xc=1.0 #Strength Fiber
        yt=1.0
        yc=1.0 #Strength Matrix
        xy=1.0 #Strength Shear
        if FailureCriteria[1]>=1.0: #Matrix Failure
            n=0.1
            g=0.1
            yt=10000000000.
            yc=10000000000.
            xy=10000000000.
            #n=((1-nr)/(1+(c*((fe-1)**e))))+nr                   
        if FailureCriteria[0]>=1.0: #Fiber failure
            s=0.1
            n=0.1
            g=0.1
            xt=10000000000.
            xc=10000000000.
            yt=10000000000.
            yc=10000000000.
            xy=10000000000.
        nvector=[s,n,g,1.0,1.0,xt,xc,yt,yc,xy,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]   
    
        #([[40000.0,11000.0,5000.0,0.3,-0.0879,1000.,500.,60.,30.,40.,-0.5,10.,30.,0.,1.,10.,23.0103,0.,1.]]) 
        
        return nvector    
    
    #Checked 04.01.2018, found error in nvector, this is corrected
    
class CalculationFunctions:
       
    def __init__(self):
        Check="OK"       
        
    
    def LayupStrength(self,layup=np.array([[]]),materials=[],Pressure=0.,Diameter=0.,DamagedLayers=0,a=0.,b=0.,angleofweakness=0.,d=0.,ThickWalled=True):
        lamt=LaminateTheory() 
        scfs=StressConcentrationFunctions()
        pres=PressureFunctions()
        fcrit=FailureCriteria()
        geom=GeometricFunctions()
        fatg=FatigueFunctions()  
        progf=ProgressiveFailureFunctions()
        calc=CalculationFunctions()
        
        FeTsaiWuVector=[]
        FeMaxStressVector=[]
        FeHashinVector=[]
        FeMaxStressFailModeVector=[]
        ResidualBurstVector=[]    
        ResidualBurstVectorHashin=[]
        CyclesUntilFailureVector=[]

        for l in range(len(layup)):
    
            MaterialNo=int(layup[l][3])
        
            #------------------------------ Found material number ----------------------------------------------#
            
            G=lamt.Global_Stiffness_Matrix(LayerData=layup,Materials=materials)

            #------------------------------ Established global stiffness matrix -----------------------------------#
            
            HydrostaticForce=[[pres.AxialLaminateForce(Pressure,Diameter+(lamt.WallThickness(LayerData=layup)))],[pres.HoopLaminateForce(Pressure,Diameter+(lamt.WallThickness(LayerData=layup)))]]
            
            # Define Diameter to half into the wall for thickwalled tubes, seeing as the SCF for thick walled benchmarks at the midplane of the laminate.
            
            #------------------------------ Found Hydrostatic force vector --------------------------------#
            
            G_hsplit=np.hsplit(G,3)
            
            G2by2=np.vsplit(G_hsplit[0],3)[0]
            
            Strain2by2=np.linalg.solve(G2by2,HydrostaticForce)
            
            Strain=np.vstack((Strain2by2,[[0.],[0.],[0.],[0.]]))
            
            Force=G*Strain
            
            #Down to here function is checked by benchmarking against Matrixes.xlsx, 05.10.2017
            
            #------------------------------ Found full force vector --------------------------------#
            
            MidPlane=lamt.LayerMidPlane(LayerData=layup,LayerNumber=l+1)
            
            WT=lamt.WallThickness(LayerData=layup)
            
            #-------------------------------- Found midplane position in layer and wall thickness ------------------------------#
            
            LayerStrain=lamt.LayerStrainfromGlobalStrain(GlobalStrain=Strain,Position=MidPlane)
            
            LocalStrain=lamt.Local_Strain_Transformation_Matrix(LayerData=layup[l])*LayerStrain
            
            LocalStress=lamt.Local_Stiffness_Matrix(MaterialProperties=materials[MaterialNo-1])*LocalStrain
            
            
            #All functions checked to here 05.10.2017 or before.


            #------------------- Established Local stress and strain ---------------------------------#        
            
            TWHSCF=pres.ThickWalledHoopSCF(Pressure=Pressure,Diameter=Diameter,WallThickness=WT,z=MidPlane)
            
            #---------------- Established thick walled cylinder stress amplification in Hoop direction ---------------# 
            
            if ThickWalled:
            
                LocalStressGlobalCoord=lamt.Local_Stress_Transformation_Matrix(LayerData=layup[l]).getI()*LocalStress
            
                LocalStressGlobalCoord=np.matrix([[float(LocalStressGlobalCoord[0])],[float(LocalStressGlobalCoord[1])*TWHSCF],[float(LocalStressGlobalCoord[2])]])
            
                LocalStress=lamt.Local_Stress_Transformation_Matrix(LayerData=layup[l])*LocalStressGlobalCoord
            
            #-------------- Established local amplified stress from thick walled distribution. -------------------#
            
            if l>=len(layup)-DamagedLayers:
        
                width=geom.WidthNormalToFiber(a=a,b=b,plyangle=layup[l][1],weaknessangle=angleofweakness)
                SCF=scfs.SCFAtXWhitneyNuismer(distancefromcenterofcrack=(width/2)+d,halfcracklength=width/2)
                LocalStress=np.matrix([[float(LocalStress[0])*SCF],[float(LocalStress[1])],[float(LocalStress[2])]])
        
            #Checked down to here 05.10.2017
        
            #---------------- Established Local Stress Concentration due to crack ------------------#
            
            FeTsaiWu=fcrit.TsaiWuExpFact(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeTsaiWuVector.append(FeTsaiWu)
            
            FeMaxStress=fcrit.MaxStressExpVector(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeMaxStressVector.append(FeMaxStress)
            
            FeMaxStressFailMode=fcrit.MaxStressFailMode(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeMaxStressFailModeVector.append(FeMaxStressFailMode)
            
            FeHashin=fcrit.HashinExpVector(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeHashinVector.append(FeHashin)
            
            ResidualBurstVector.append((1./(FeMaxStress[0]))*Pressure)
            
            ResidualBurstVectorHashin.append((1./(FeHashin[0]))*Pressure)

            #------------------ Established exposure factors ---------------------------# 
            
            CyclesUntilFailure=fatg.CyclesUntilFailure(Material=materials[MaterialNo-1],Stress=LocalStress,Thickness=layup[l][2])    
            CyclesUntilFailureVector.append(CyclesUntilFailure)
        
            #---------------- Established cycles until failure --------------------------#
        
        StrengthVector = {'FeTsaiWuVector': FeTsaiWuVector, 'FeMaxStressVector': FeMaxStressVector, 'FeMaxStressFailModeVector': FeMaxStressFailModeVector,'ResidualBurstVector': ResidualBurstVector, 'CyclesUntilFailureVector': CyclesUntilFailureVector,'ForceVector': Force,'StrainVector': Strain, 'FeHashinVector': FeHashinVector,'ResidualBurstVectorHashin': ResidualBurstVectorHashin}
        
        StrengthVector = pd.Series(StrengthVector)    
        
        return StrengthVector
    
    def LayupStrengthForceInput(self,layup=np.array([[]]),materials=[],Force=[[0.],[0.],[0.],[0.],[0.],[0.]],Diameter=0.,DamagedLayers=0,a=0.,b=0.,angleofweakness=0.,d=0.,ThickWalled=True):
        lamt=LaminateTheory() 
        scfs=StressConcentrationFunctions()
        pres=PressureFunctions()
        fcrit=FailureCriteria()
        geom=GeometricFunctions()
        fatg=FatigueFunctions()  
        progf=ProgressiveFailureFunctions()
        calc=CalculationFunctions()
        
        FeTsaiWuVector=[]
        FeMaxStressVector=[]
        FeHashinVector=[]
        FeMaxStressFailModeVector=[]
        ResidualBurstVector=[]    
        ResidualBurstVectorHashin=[]
        CyclesUntilFailureVector=[]
        LocalStressVector=[]
        LocalStressGlobalCoordVector=[]

        for l in range(len(layup)):
    
            MaterialNo=int(layup[l][3])
        
            #------------------------------ Found material number ----------------------------------------------#
            
            G=lamt.Global_Stiffness_Matrix(LayerData=layup,Materials=materials)

            #------------------------------ Established global stiffness matrix -----------------------------------#
            
            Strain=(np.linalg.inv(G))*Force
            
            #Down to here function is checked by benchmarking against Matrixes.xlsx, 05.10.2017
            
            #------------------------------ Found full force vector --------------------------------#
            
            MidPlane=lamt.LayerMidPlane(LayerData=layup,LayerNumber=l+1)
            
            WT=lamt.WallThickness(LayerData=layup)
            
            #-------------------------------- Found midplane position in layer and wall thickness ------------------------------#
            
            LayerStrain=lamt.LayerStrainfromGlobalStrain(GlobalStrain=Strain,Position=MidPlane)
            
            LocalStrain=lamt.Local_Strain_Transformation_Matrix(LayerData=layup[l])*LayerStrain
            
            LocalStress=lamt.Local_Stiffness_Matrix(MaterialProperties=materials[MaterialNo-1])*LocalStrain
            
            
            #All functions checked to here 05.10.2017 or before.


            #------------------- Established Local stress and strain ---------------------------------#        
            
            #TWHSCF=pres.ThickWalledHoopSCF(Pressure=Pressure,Diameter=Diameter,WallThickness=WT,z=MidPlane)
            
            #---------------- Established thick walled cylinder stress amplification in Hoop direction ---------------# 
            
            if ThickWalled:
            
                LocalStressGlobalCoord=lamt.Local_Stress_Transformation_Matrix(LayerData=layup[l]).getI()*LocalStress
            
                LocalStressGlobalCoord=np.matrix([[float(LocalStressGlobalCoord[0])],[float(LocalStressGlobalCoord[1])*TWHSCF],[float(LocalStressGlobalCoord[2])]])
            
                LocalStress=lamt.Local_Stress_Transformation_Matrix(LayerData=layup[l])*LocalStressGlobalCoord
            
            #-------------- Established local amplified stress from thick walled distribution. -------------------#
            
            if l>=len(layup)-DamagedLayers:
        
                width=geom.WidthNormalToFiber(a=a,b=b,plyangle=layup[l][1],weaknessangle=angleofweakness)
                SCF=scfs.SCFAtXWhitneyNuismer(distancefromcenterofcrack=(width/2)+d,halfcracklength=width/2)
             
                LocalStress=np.matrix([[float(LocalStress[0])*SCF],[float(LocalStress[1])],[float(LocalStress[2])]])
        
            #Checked down to here 05.10.2017
        
            #---------------- Established Local Stress Concentration due to crack ------------------#
            
            LocalStressVector.append(LocalStress)
            
            LocalStressGlobalCoord=lamt.Local_Stress_Transformation_Matrix(LayerData=layup[l]).getI()*LocalStress
            
            LocalStressGlobalCoordVector.append(LocalStressGlobalCoord)
            
            
            #--------------- Established Local Stress Vectors for printout -------------------------#
            
            FeTsaiWu=fcrit.TsaiWuExpFact(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeTsaiWuVector.append(FeTsaiWu)
            
            FeMaxStress=fcrit.MaxStressExpVector(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeMaxStressVector.append(FeMaxStress)
            
            FeMaxStressFailMode=fcrit.MaxStressFailMode(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeMaxStressFailModeVector.append(FeMaxStressFailMode)
            
            FeHashin=fcrit.HashinExpVector(Material=materials[MaterialNo-1],LocalStressVector=LocalStress)
            FeHashinVector.append(FeHashin)
            
            #ResidualBurstVector.append((1./(FeMaxStress[0]))*Pressure)
            
            #ResidualBurstVectorHashin.append((1./(FeHashin[0]))*Pressure)

            #------------------ Established exposure factors ---------------------------# 
            
            CyclesUntilFailure=fatg.CyclesUntilFailure(Material=materials[MaterialNo-1],Stress=LocalStress,Thickness=layup[l][2])    
            CyclesUntilFailureVector.append(CyclesUntilFailure)
        
            #---------------- Established cycles until failure --------------------------#
        
        StrengthVector = {'FeTsaiWuVector': FeTsaiWuVector, 'FeMaxStressVector': FeMaxStressVector,'FeHashinVector': FeHashinVector, 'FeMaxStressFailModeVector': FeMaxStressFailModeVector,'CyclesUntilFailureVector': CyclesUntilFailureVector,'ForceVector': Force,'StrainVector': Strain,'LocalStressVector':LocalStressVector,'LocalStressGlobalCoordVector':LocalStressGlobalCoordVector}
        
        StrengthVector = pd.Series(StrengthVector)    
        
        return StrengthVector
    
#-----------------------------Degraded Layup function is old----------------------------------------------#
    
    def DegradedLayup(self,Layup=[],MinerSum=[0.,0.],FeMaxStressVector=[],fatigueanalysis=True):
        DegradedLayup=[]
        for l in range(len(Layup)):
            if float(Layup[l][3])==4:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(4)])
                
            elif FeMaxStressVector[l][0]>=1.0 and FeMaxStressVector[l][1]>=1.0:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(4)])
                
            #----------------- Checked if already destroyed, if not then: ---------------#       
            
            elif float(Layup[l][3])==2 and FeMaxStressVector[l][1]>=1.0:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(4)])
                
            elif float(Layup[l][3])==3 and FeMaxStressVector[l][0]>=1.0:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(4)])
    
            #--------------- Check if already either destroyed by max stress in matrix or fiber, if not: ----------------#            
                     
            elif FeMaxStressVector[l][0]>=1.0:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(4)])
                
            #-------------- If fiber failure, then destroyed ------------------------------------------------------------#            
                
            elif FeMaxStressVector[l][1]>=1.0:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(3)])
                
            #-------------- If matrix failure, then compliant matrix. -------------------------------------------------#
            
            elif MinerSum[l,0]>=1.0 and MinerSum[l,1]>=1.0 and fatigueanalysis:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(4)])            
            
            
            elif MinerSum[l,0]>=1.0 and fatigueanalysis:    
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(4)])         
    
            
            elif MinerSum[l,1]>=1.0 and fatigueanalysis:
                DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),int(3)])
           
            #--------------- Check if destroyed by fatigue, if not: ----------------#  
            
            else:
                DegradedLayup.append(list(Layup[l]))
        DegradedLayup=np.array(DegradedLayup)     
        
        return DegradedLayup
    
#-----------------------------Degraded Layup function is old----------------------------------------------#
    
    def GraduallyDegradedLayup(self,Layup=[]):
        DegradedLayup=[]
        for l in range(len(Layup)):
            DegradedLayup.append([int(Layup[l][0]),float(Layup[l][1]),float(Layup[l][2]),l+1])
        DegradedLayup=np.array(DegradedLayup)     
        
        return DegradedLayup
        
    def DegradedMaterial(self,Material=[],DegradationConstant=[],ProgressiveFailureCriteria=''):
        
        DegradedMaterial=np.array(Material)*np.array(DegradationConstant)
                
        return DegradedMaterial
    
    def DegradedMaterials(self,Layup=[],materials=[],FailureCriteriaVector=[],FatigueCriteriaVector=[],FailureCriteriaConstants=[],ProgressiveFailureCriteria=''):
        progf=ProgressiveFailureFunctions()
        calc=CalculationFunctions()
        DegradedMaterials=[]
    
        for l in range(len(Layup)):

            if ProgressiveFailureCriteria=='Instantaneous' or ProgressiveFailureCriteria=='Gradual':
                DegradationConstant=progf.ConstantStiffnessReductionFactor(FailureCriteria=FailureCriteriaVector[l],FatigueCriteria=FatigueCriteriaVector[l],ConstantConstants=FailureCriteriaConstants)
            
            if ProgressiveFailureCriteria=='Puck':
                DegradationConstant=progf.PuckStiffnessReductionFactor(FailureCriteria=FailureCriteriaVector[l],FatigueCriteria=FatigueCriteriaVector[l],PuckConstants=FailureCriteriaConstants)
                
            MaterialNo=int(Layup[l][3])
            DegradedMaterial=calc.DegradedMaterial(Material=materials[MaterialNo-1],DegradationConstant=DegradationConstant,ProgressiveFailureCriteria='')
            DegradedMaterials.append(DegradedMaterial)
        
        return DegradedMaterials
    
    
    
    
    
    