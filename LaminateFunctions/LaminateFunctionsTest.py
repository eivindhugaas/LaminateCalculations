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
import unittest



    lamt.WallThickness(self,LayerData=[[0.]])

        
    
    def Local_Strain_from_Global_Strain(self,GlobalStrain=[0.],Position=0.,LocalInverseStrainTransformationMatrix=[[0.]]):
        
        GlobalStrainTension=np.split(GlobalStrain, 2)[0]
        GlobalStrainBending=np.split(GlobalStrain, 2)[1]

        LocalStrainGlobalCoorSystem = GlobalStrainTension + (GlobalStrainBending*Position)  

        LocalStrainLocalCoorSystem=LocalInverseStrainTransformationMatrix*LocalStrainGlobalCoorSystem
        return LocalStrainLocalCoorSystem            
    
    def Local_Stiffness_Matrix(self,MaterialProperties=[0.]):
        Q11 = (MaterialProperties[0])/(1-(MaterialProperties[3]*MaterialProperties[4]))
        Q22 = (MaterialProperties[1])/(1-(MaterialProperties[3]*MaterialProperties[4]))
        Q12 = ((MaterialProperties[1])*(MaterialProperties[3]))/(1-(MaterialProperties[3]*MaterialProperties[4]))
        Q66 = MaterialProperties[2]
        LocalStiffnessMatrix = np.matrix(((Q11, Q12, 0),(Q12, Q22, 0),(0, 0, Q66)))           
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
        for i in range(1, (len(LayerData)+1)):
            Htot = Htot + LayerData[i-1,2]
            Htot_Vector.append(Htot)
        for i in range(1, (len(LayerData)+1)):
            Hk = (-H_Sum/2) + Htot_Vector[i-1]
            Hk_Vector.append(Hk)
        Hkmin1_Vector = np.subtract(Hk_Vector,H_Vector)
        Hkmin1_minHk_Vector = np.subtract(Hk_Vector,Hkmin1_Vector)
        Hkmin1_minHk_scnd_Vector = np.subtract(np.power(Hk_Vector,2),np.power(Hkmin1_Vector,2))
        Hkmin1_minHk_third_Vector = np.subtract(np.power(Hk_Vector,3),np.power(Hkmin1_Vector,3))
        
        T_Stress_I_mult_Local_Stiffness_Matrixes=[]
        T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes=[]
        
        A_i=B_i=D_i=[]
        A=B=D=[[0,0,0],[0,0,0],[0,0,0]]
        
        for i in range(len(LayerData)):
            MatNum=math.floor(LayerData[i][3]) #Finds the Layer's material
            Material=np.array(Materials[MatNum-1]) #Chooses correct material in material matrix
        
            Local_Stress_Transformation_Matrix=LaminateTheory().Local_Stress_Transformation_Matrix(LayerData=LayerData[i])
            Inverted_Local_Stress_Transformation_Matrix=LaminateTheory().Local_Stress_Transformation_Matrix(LayerData[i]).getI()
            Local_Strain_Transformation_Matrix=LaminateTheory().Local_Strain_Transformation_Matrix(LayerData[i])
            Inverted_Local_Strain_Transformation_Matrix=LaminateTheory().Local_Strain_Transformation_Matrix(LayerData[i]).getI()

            T_Stress_I_mult_Local_Stiffness = Inverted_Local_Stress_Transformation_Matrix*LaminateTheory().Local_Stiffness_Matrix(MaterialProperties=Material)
            T_Stress_I_mult_Local_Stiffness_mult_T_Strain = T_Stress_I_mult_Local_Stiffness*Local_Strain_Transformation_Matrix
            T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes.append(T_Stress_I_mult_Local_Stiffness_mult_T_Strain)
            
            A_i = (T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes[i-1])*((Hkmin1_minHk_Vector[i-1]))
            B_i = 0.5*(T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes[i-1])*((Hkmin1_minHk_scnd_Vector[i-1]))
            D_i = (1/3)*(T_Stress_I_mult_Local_Stiffness_mult_T_Strain_Matrixes[i-1])*((Hkmin1_minHk_third_Vector[i-1]))
            A=A_i + A 
            B=B_i + B 
            D=D_i + D
        E_L=np.concatenate((A, B), axis=0)
        E_R=np.concatenate((B, D), axis=0)
        E = np.hstack((E_L, E_R))
        return E  
    
    def LayerMidPlane(self,LayerData=[[0.]],Layer=1000):
        LaminateThickness=0.
        for i in range(len(LayerData)):
            LaminateThickness=LaminateThickness+LayerData[i][2]
        BottomPlane=(-LaminateThickness/2)+(LayerData[0][2]/2)
        LayerMidPlane=BottomPlane
        LayerMidPlaneVector=[BottomPlane]
        for i in range((len(LayerData)-1)):
            LayerMidPlane=LayerMidPlane+LayerData[i][2]/2+LayerData[i+1][2]/2
            LayerMidPlaneVector.append(LayerMidPlane)
        return(LayerMidPlaneVector[Layer-1])

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
        
class GeometricFunctions:
    def __init__(self):
        Check="OK"
        
    def WidthNormalToFiber(self,a=0.,b=0.,plyangle=0.,weaknessangle=0.):
        plyangle=plyangle+90
        weaknessangle=weaknessangle+90
        x=b*math.cos(math.radians(plyangle+weaknessangle))
        y=a*math.sin(math.radians(plyangle+weaknessangle))
        width=((x**2+y**2)**0.5)