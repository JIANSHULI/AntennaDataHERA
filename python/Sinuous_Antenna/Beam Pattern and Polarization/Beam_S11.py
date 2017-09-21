
# coding: utf-8

# In[1]:

#get_ipython().magic(u'matplotlib inline')
#import numpy as n
#import matplotlib.pyplot as p
#import healpy as hp
import numpy as n
import matplotlib
import matplotlib.pyplot as p
import healpy as hp
import copy
import scipy.optimize as op
from matplotlib import pylab as plt
import gainData as gainData
pi=n.pi
c=299792458.

#rotate
def rotateBeam(inputMap,rot=[90,0,0]):
    rotator=hp.Rotator(rot=rot)
    npix=len(inputMap)
    nside=hp.npix2nside(npix)
    theta,phi=hp.pix2ang(nside,range(npix))
    newtheta,newphi=rotator(theta,phi)
    output=hp.get_interp_val(inputMap,newtheta,newphi)
    return output
    
def rotateBeam_to_Z(inputMap,rot=[0,0,90]):
    rotator=hp.Rotator(rot=rot)
    npix=len(inputMap)
    nside=hp.npix2nside(npix)
    theta,phi=hp.pix2ang(nside,range(npix))
    newtheta,newphi=rotator(theta,phi)
    output=hp.get_interp_val(inputMap,newtheta,newphi)
    return output
    
def rotateBeam_to_Z_Yneg(inputMap,rot=[0,0,-90]):
    rotator=hp.Rotator(rot=rot)
    npix=len(inputMap)
    nside=hp.npix2nside(npix)
    theta,phi=hp.pix2ang(nside,range(npix))
    newtheta,newphi=rotator(theta,phi)
    output=hp.get_interp_val(inputMap,newtheta,newphi)
    return output
    
def hpCut(phi,nPix,data):
    nSide=hp.npix2nside(len(data))
    output=n.zeros(nPix)
    thetaVals=n.arange(nPix/2)/(nPix/2.)*pi/2.
    thetaVals=n.hstack([n.flipud(thetaVals),thetaVals,]).T
    phiVals=n.ones(len(thetaVals))
    phi1=phi+pi
    phiVals[:nPix/2]=phi1
    phiVals[nPix/2:]=phi
    output=data[hp.ang2pix(nSide,thetaVals,phiVals)]
    return output


class BeamSinous_Y:
    def __init__(self,dirName,fList,nside,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number,pols=['XX','YY'],rotateY=False):
        self.nf=len(fList)
        self.fAxis=n.zeros(self.nf)
        self.npolsOriginal=len(pols)
        self.npols=max(len(pols),2)
        self.solidAngles=n.zeros((self.npols,self.nf))
        self.effArea=n.zeros_like(self.solidAngles)
        self.ellipticity=n.zeros(self.nf)
        self.nPix=hp.nside2npix(nside)
        self.nSide=nside
        self.pixArea=hp.nside2pixarea(self.nSide)
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
        theta=n.round(n.degrees(theta)).astype(int)
        phi=n.round(n.degrees(phi)).astype(int)
        self.data=n.zeros((self.npols,self.nf,self.nPix))
        self.data_N=n.zeros((self.npols,self.nf,self.nPix))
        if(rotateY):
            pols.append('YY')
        self.pols=pols
        for m in range(self.nf):            
            print m
            tempf=fList[m].split('p')
            self.fAxis[m]=float(tempf[0])*1e6
            if(len(tempf)>1):
                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
            for np in range(self.npolsOriginal):
                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%s-%s_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
                #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
                self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
                self.data_N[np,m,:]=self.data[np,m,:]/self.data[np,m,:].flatten().max(); 
                self.data[np,m,theta>90.]=0.
                self.data_N[np,m,theta>90.]=0.
                self.solidAngles[np,m]=self.pixArea*n.sum(self.data_N[np,m,:])
                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
                                            
                                            
            if(self.npolsOriginal==1):
                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
                self.data_N[1,m,:]=rotateBeam(self.data_N[0,m,:].flatten())
                self.solidAngles[1,m]=self.pixArea*n.sum(self.data_N[1,m,:])
                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)                

class BeamSinous_BackPlane_Y:
    def __init__(self,dirName,fList,nside,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number,pols=['XX','YY'],rotateY=False):
        self.nf=len(fList)
        self.fAxis=n.zeros(self.nf)
        self.npolsOriginal=len(pols)
        self.npols=max(len(pols),2)
        self.solidAngles=n.zeros((self.npols,self.nf))
        self.effArea=n.zeros_like(self.solidAngles)
        self.ellipticity=n.zeros(self.nf)
        self.nPix=hp.nside2npix(nside)
        self.nSide=nside
        self.pixArea=hp.nside2pixarea(self.nSide)
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
        theta=n.round(n.degrees(theta)).astype(int)
        phi=n.round(n.degrees(phi)).astype(int)
        self.data=n.zeros((self.npols,self.nf,self.nPix))
        self.data_N=n.zeros((self.npols,self.nf,self.nPix))
        if(rotateY):
            pols.append('YY')
        self.pols=pols
        for m in range(self.nf):            
            print m
            tempf=fList[m].split('p')
            self.fAxis[m]=float(tempf[0])*1e6
            if(len(tempf)>1):
                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
            for np in range(self.npolsOriginal):
                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,fList[m],Port_Number),skiprows=2);
                #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
                self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
                self.data_N[np,m,:]=self.data[np,m,:]/self.data[np,m,:].flatten().max(); 
                self.data[np,m,theta>90.]=0.
                self.data_N[np,m,theta>90.]=0.
                self.solidAngles[np,m]=self.pixArea*n.sum(self.data_N[np,m,:])
                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
                                            
                                            
            if(self.npolsOriginal==1):
                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
                self.data_N[1,m,:]=rotateBeam(self.data_N[0,m,:].flatten())
                self.solidAngles[1,m]=self.pixArea*n.sum(self.data_N[1,m,:])
                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)                

class BeamSinous_NoDish_Y:
    def __init__(self,dirName,fList,nside,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number,pols=['XX','YY'],rotateY=False):
        self.nf=len(fList)
        self.fAxis=n.zeros(self.nf)
        self.npolsOriginal=len(pols)
        self.npols=max(len(pols),2)
        self.solidAngles=n.zeros((self.npols,self.nf))
        self.effArea=n.zeros_like(self.solidAngles)
        self.ellipticity=n.zeros(self.nf)
        self.nPix=hp.nside2npix(nside)
        self.nSide=nside
        self.pixArea=hp.nside2pixarea(self.nSide)
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
        theta=n.round(n.degrees(theta)).astype(int)
        phi=n.round(n.degrees(phi)).astype(int)
        self.data=n.zeros((self.npols,self.nf,self.nPix))
        self.data_N=n.zeros((self.npols,self.nf,self.nPix))
        if(rotateY):
            pols.append('YY')
        self.pols=pols
        for m in range(self.nf):            
            print m
            tempf=fList[m].split('p')
            self.fAxis[m]=float(tempf[0])*1e6
            if(len(tempf)>1):
                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
            for np in range(self.npolsOriginal):
                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%s-%s_band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
                #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
                self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
                self.data_N[np,m,:]=self.data[np,m,:]/self.data[np,m,:].flatten().max(); 
                self.data[np,m,theta>90.]=0.
                self.data_N[np,m,theta>90.]=0.
                self.solidAngles[np,m]=self.pixArea*n.sum(self.data_N[np,m,:])
                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
                                            
                                            
            if(self.npolsOriginal==1):
                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
                self.data_N[1,m,:]=rotateBeam(self.data_N[0,m,:].flatten())
                self.solidAngles[1,m]=self.pixArea*n.sum(self.data_N[1,m,:])
                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)                

class BeamSinous_NoDish_BackPlane_Y:
    def __init__(self,dirName,fList,nside,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number,pols=['XX','YY'],rotateY=False):
        self.nf=len(fList)
        self.fAxis=n.zeros(self.nf)
        self.npolsOriginal=len(pols)
        self.npols=max(len(pols),2)
        self.solidAngles=n.zeros((self.npols,self.nf))
        self.effArea=n.zeros_like(self.solidAngles)
        self.ellipticity=n.zeros(self.nf)
        self.nPix=hp.nside2npix(nside)
        self.nSide=nside
        self.pixArea=hp.nside2pixarea(self.nSide)
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
        theta=n.round(n.degrees(theta)).astype(int)
        phi=n.round(n.degrees(phi)).astype(int)
        self.data=n.zeros((self.npols,self.nf,self.nPix))
        self.data_N=n.zeros((self.npols,self.nf,self.nPix))
        if(rotateY):
            pols.append('YY')
        self.pols=pols
        for m in range(self.nf):            
            print m
            tempf=fList[m].split('p')
            self.fAxis[m]=float(tempf[0])*1e6
            if(len(tempf)>1):
                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
            for np in range(self.npolsOriginal):
                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,fList[m],Port_Number),skiprows=2);
                #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
                self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
                self.data_N[np,m,:]=self.data[np,m,:]/self.data[np,m,:].flatten().max(); 
                self.data[np,m,theta>90.]=0.
                self.data_N[np,m,theta>90.]=0.
                self.solidAngles[np,m]=self.pixArea*n.sum(self.data_N[np,m,:])
                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
                                            
                                            
            if(self.npolsOriginal==1):
                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
                self.data_N[1,m,:]=rotateBeam(self.data_N[0,m,:].flatten())
                self.solidAngles[1,m]=self.pixArea*n.sum(self.data_N[1,m,:])
                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)                

class BeamSinousPro_BackPlane_Y:
    def __init__(self,dirName,fList,nside,Skirt_Status,BackPlane_Status,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Status,Port_Number,Impedance_Port_Status,Dish_Status,Direction,pols=['XX','YY'],rotateY=False):
        self.nf=len(fList)
        self.fAxis=n.zeros(self.nf)
        self.npolsOriginal=len(pols)
        self.npols=max(len(pols),2)
        self.solidAngles=n.zeros((self.npols,self.nf))
        self.effArea=n.zeros_like(self.solidAngles)
        self.ellipticity=n.zeros(self.nf)
        self.nPix=hp.nside2npix(nside)
        self.nSide=nside
        self.pixArea=hp.nside2pixarea(self.nSide)
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
        theta=n.round(n.degrees(theta)).astype(int)
        phi=n.round(n.degrees(phi)).astype(int)
        self.data=n.zeros((self.npols,self.nf,self.nPix))
        self.data_N=n.zeros((self.npols,self.nf,self.nPix))
        if(rotateY):
            pols.append('YY')
        self.pols=pols
        for m in range(self.nf):            
            print m
            tempf=fList[m].split('p')
            self.fAxis[m]=float(tempf[0])*1e6
            if(len(tempf)>1):
                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
            for np in range(self.npolsOriginal):
                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%s-%s_%s%s%s%s%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,fList[m],Port_Number),skiprows=2);
                #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
                if Direction == 1:
                    self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
                elif Direction == -1:
                    self.data[np,m,:]=rotateBeam_to_Z_Yneg(self.data[np,m,:].flatten())
                self.data_N[np,m,:]=self.data[np,m,:]/self.data[np,m,:].flatten().max();
                #self.data[np,m,theta>90.]=0.
                #self.data_N[np,m,theta>90.]=0.
                self.solidAngles[np,m]=self.pixArea*n.sum(self.data_N[np,m,:])
                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
                                            
                                            
            if(self.npolsOriginal==1):
                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
                self.data_N[1,m,:]=rotateBeam(self.data_N[0,m,:].flatten())
                self.solidAngles[1,m]=self.pixArea*n.sum(self.data_N[1,m,:])
                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)



def calc_Tsky(nu, fq0=.180, T0=180., alpha=-2.5):
    '''Tsky = T0 * (nu/fq0)^alpha. nu, fq0 in GHz, T0 in K.''' 
    return T0 * (nu/fq0)**alpha
    
nu = n.linspace(.050,.250, 1024) 
Trx1 = 100 # K
Trx2 = 75 # K
Tsky = calc_Tsky(nu)

f_obs = 1.5
e_dish = 0.74
Trx_lotarg = n.around(e_dish * (500 - Tsky))
Tsys = n.sqrt(f_obs) * Tsky
Trx_hitarg = e_dish * (Tsys - Tsky)
Trx_targ = n.where(Trx_lotarg > Trx_hitarg, Trx_lotarg, Trx_hitarg).clip(0)
Trx_targ = n.where(Trx1 * n.ones_like(nu) > Trx_targ, Trx1 * n.ones_like(nu), Trx_targ).clip(0)

S11_1 = (1 - Trx1 / Trx_targ).clip(0,1)
S11_2 = (1 - Trx2 / Trx_targ).clip(0,1)
S11_3 = (1 - 85 / Trx_targ).clip(0,1)


N = 0

Growth_Rate_List = ['80']
Outer_Diameter_List = ['175']
Inner_Diameter_List = ['30']
Band_Resistance_List = ['15']
Skirt_Diameter_List = ['1.2']
Skirt_Height_List = ['0.3']
BackPlane_Height_List = ['20']
BackPlane_Diameter_List = ['0.99']

#Frequency_List = ['50','75','100','125','150','175','200','225','250','275','300']
Frequency_List = ['40','50','67','94','121','148','150','175','202','229','250','256','283']    #[40,50,67,94,121,148,150,175,202,229,250,256,283]
Port_Number_List = ['1']
Phi_List = [0,pi/2.0]
PhiDeg_List = [0,90]

Band_Status_List = ['band-no-']
Skirt_Status_List = ['skirt-1.2-0.3-']
#BackPlane_Status_List = ['backcone-15-135-10-0.99','backcone-15-80-5-0.99']
BackPlane_Status_List = ['backplane-50-0.99']
Impedance_Port_Status_List = ['imp-100-']
#Dish_Status_List = ['no-']
Dish_Status_List = ['dish-4.5-']


Simplify = 0
BackPlane = 1
BeamPattern = 1
Dish = 2
S11_Format = 1 # 0:old,1:new
Direction = 1

if Dish == 2:                                                
    for Growth_Rate in Growth_Rate_List:
        for Outer_Diameter in Outer_Diameter_List:
            for Inner_Diameter in Inner_Diameter_List:
                for Port_Number in Port_Number_List:   
            
                    for Impedance_Port_Status in Impedance_Port_Status_List:
                        for Dish_Status in Dish_Status_List:
                            for Band_Status in Band_Status_List:
                                for Skirt_Status in Skirt_Status_List:
                                    for BackPlane_Status in BackPlane_Status_List:
                                        
                                        if BeamPattern == 1:
                                            
                                            BeamSinuousDishBandSkirt = BeamSinousPro_BackPlane_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Status,BackPlane_Status,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Status,Port_Number,Impedance_Port_Status,Dish_Status,Direction, ['XX'],rotateY=True)
                                            
                                            p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,'o',label='Sinuous_Dish-Band-Skirt-BackPlane')
                                            p.xlabel('f (MHz)',fontsize=20)
                                            p.ylabel('$\\xi$',fontsize=20)
                                            p.legend(loc='best',fontsize=10,ncol=1)
                                            p.yscale('log')
                                            p.title('FarEllip_Y_0.%s-%s-%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number))
                                            p.gca().tick_params('x',labelsize=16)
                                            p.gca().tick_params('y',labelsize=16)
                                            p.gcf().set_size_inches([8,6])
                                            p.gca().yaxis.grid(which='minor')
                                            p.grid()   
                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Y_0.%s-%s-%s_%s%s%s%s%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),bbox_inches='tight')
                                            p.close()
                                            
                                            n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Data_Y_0.%s-%s-%s_%s%s%s%s%s-%s.txt' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity],fmt=['%.2f','%.4f']) 
                                            
                                            p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),'o',label='Sinuous_Dish-Band-Skirt-BackPlane')
                                            p.xlabel('f (MHz)',fontsize=20)
                                            p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
                                            p.grid()
                                            #p.ylim(.15,.85)
                                            p.gca().tick_params('y',labelsize=16)
                                            p.gcf().set_size_inches([8,6])
                                            #p.gca().yaxis.grid(which='minor')
                                            p.legend(loc='best',fontsize=10,ncol=1) 
                                            p.title('FarEffecArea_Y_0.%s-%s-%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number))                           
                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Y_0.%s-%s-%s_%s%s%s%s%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),bbox_inches='tight')
                                            p.close()
                                            
                                            n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Data_Y_0.%s-%s-%s_%s%s%s%s%s-%s.txt' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7)],fmt=['%.2f','%.4f']) 
                                            
                                            for m in range(len(Frequency_List)):                                                                    
                                                nth=16000
                                                tha=n.degrees(n.arange(-nth/2,nth/2)*2*pi/nth)
                                                l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
                                                l1=p.plot(tha,hpCut(Phi_List[1],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k')[0]
                                                p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-BackPlane-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-BackPlane-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
                                                p.grid()
                                                p.xlabel('$\\theta$')
                                                p.xlim(-180,180)
                                                p.ylim(-30,30)
                                                p.ylabel('Directivity (dB)')
                                                p.gcf().set_size_inches([10,7])
                                                p.title('FarCut_Y_0.%s-%s-%s_%s%s%s%s%s-%s-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]))
                                                p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_%s%s%s%s%s-%s-%s-%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]),bbox_inches='tight')
                                                p.close()
                                                
                                            phi = m = 0
                                                
                                            for phi in range(len(Phi_List)):
                                                for m in range(len(Frequency_List)):
                                                    if Simplify == 1:                                                        
                                                        if m!=0 and m!=5 and m!=11 and m!=12:                                                                                                                
                                                            nth=16000
                                                            tha=n.degrees(n.arange(-nth/2,nth/2)*2*pi/nth)
                                                    #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
                                                            p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousDishBandSkirtBackPlane-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
                                                    #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
                                                    if Simplify == 0:                                                        
                                                        nth=16000
                                                        tha=n.degrees(n.arange(-nth/2,nth/2)*2*pi/nth)
                                                #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
                                                        p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousDishBandSkirtBackPlane-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
                                                #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
                                                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                                                            ncol=1, mode="expand", borderaxespad=1)
                                                p.grid()
                                                p.xlabel('$\\theta$')
                                                p.xlim(-180,180)
                                                p.ylim(-30,30)
                                                p.ylabel('Directivity (dB)')
                                                p.gcf().set_size_inches([10,7])
                                                p.title('FarCut_Y_0.%s-%s-%s_%s%s%s%s%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number,PhiDeg_List[phi]))
                                                p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_%s%s%s%s%s-%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number,PhiDeg_List[phi]),bbox_inches='tight')
                                                p.close()
                                        
                                                                                
                                        for N in range(2):
                                            PW = S11_Power = N+1      
                                            if S11_Format == 1:                                                
                                                fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%s-%s-%s_%s%s%s%s%s-%s.txt' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number)
                                                fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%s-%s-%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number)
                                                fileNameImpedanceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Impedance_0.%s-%s-%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number)
                                            
                                            elif S11_Format == 0:
                                                fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
                                                fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
                                                fileNameImpedanceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Impedance_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
                                                

                                            FLOW=0.05
                                            FHIGH=0.25


                                            gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
                                                                                fileType='CST_TimeTrace',
                                                                                fMin=FLOW,fMax=FHIGH,
                                                                                comment='S11 derived from cst time domain data')
                                            gainData_cst=gainData.GainData(fileNameS11CST,
                                                                            fileType='CST_S11',
                                                                            fMin=FLOW,fMax=FHIGH,
                                                                            comment='S11 obtained directly from cst')
                                            gainData_impedance=gainData.GainData(fileNameImpedanceCST,
                                                                            fileType='CST_Impedance',
                                                                            fMin=FLOW,fMax=FHIGH,
                                                                            comment='Impedance obtained directly from cst')
                                            #gainData_far=
                                            #gainData_vna=gainData.GainData(fileNameS11VNA,
                                            #							   fileType='VNAHP_S11',
                                            #							   fMin=FLOW,fMax=FHIGH,
                                            #							   comment='s11 obtained from richs vna measurement')

                                            print gainData_cst.gainFrequency.shape
                                            print gainData_impedance.gainFrequency.shape

                                            #first make original plot comparing s11 of time trace and s11 of vna

                                            #p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
                                            p.plot(gainData_timeTrace.tAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
                                            p.plot(gainData_cst.tAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
                                            p.xlim(-30,400)
                                            p.ylim(-70*S11_Power,0)
                                            p.ylabel('|$\widetilde{S}_{11}$|(dB)')
                                            p.xlabel('delay (ns)')
                                            p.legend(loc='best')
                                            p.title('S11_CST_Delay_0.%s-%s-%sPW%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,S11_Power,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number))
                                            #p.show()
                                            p.grid()
                                            #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%s-%s-%sPW%sCr_%s%s%s%s%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,S11_Power,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),bbox_inches='tight')
                                            p.close()

                                            #p.plot(gainData_vna.fAxis,10.*n.log10(n.abs(gainData_vna.gainFrequency)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
                                            p.plot(gainData_timeTrace.fAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainFrequency)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
                                            p.plot(gainData_cst.fAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainFrequency)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
                                            p.plot(nu, 5*PW*n.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
                                            p.xlim(.045,.255)
                                            p.ylim(-25*S11_Power,0)
                                            p.ylabel('|S$_{11}$|(dB)')
                                            p.xlabel('f (GHz)')
                                            p.legend(loc='best')
                                            p.title('S11_CST_Frequency_0.%s-%s-%sPW%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,S11_Power,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number)) 
                                            #p.show()
                                            p.grid()
                                            #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%s-%s-%sPW%sCr_%s%s%s%s%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,S11_Power,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),bbox_inches='tight')
                                            p.close()
                                            
                                            p.plot(gainData_impedance.fAxis,n.abs(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Abs',markersize=4,markeredgecolor='none')
                                            p.xlim(.045,.255)
                                            p.ylabel('|Impedance(Abs)/Ohm')
                                            p.xlabel('f (GHz)')
                                            p.legend(loc='best')
                                            p.title('ImpedanceAbs_CST_Frequency_0.%s-%s-%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number)) 
                                            #p.show()
                                            p.grid()
                                            #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedanceAbs_CST_Frequency_0.%s-%s-%s_%s%s%s%s%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),bbox_inches='tight')
                                            p.close()
                                            
                                            p.plot(gainData_impedance.fAxis,n.angle(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Pha',markersize=4,markeredgecolor='none')
                                            p.xlim(.045,.255)
                                            p.ylabel('|Impedance(Pha)/deg')
                                            p.xlabel('f (GHz)')
                                            p.legend(loc='best')
                                            p.title('ImpedancePha_CST_Frequency_0.%s-%s-%s_%s%s%s%s%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number)) 
                                            #p.show()
                                            p.grid()
                                            #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedancePha_CST_Frequency_0.%s-%s-%s_%s%s%s%s%s-%s.pdf' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Dish_Status,Impedance_Port_Status,Band_Status,Skirt_Status,BackPlane_Status,Port_Number),bbox_inches='tight')
                                            p.close()
    
#else:
#    for N in range(2):
#        PW = S11_Power = N+1                             
#        for Growth_Rate in Growth_Rate_List:
#            for Outer_Diameter in Outer_Diameter_List:
#                for Inner_Diameter in Inner_Diameter_List:                    
#                    for Band_Resistance in Band_Resistance_List:
#                        for Skirt_Diameter in Skirt_Diameter_List:
#                            for Skirt_Height in Skirt_Height_List:
#                                for BackPlane_Height in BackPlane_Height_List:
#                                    for BackPlane_Diameter in BackPlane_Diameter_List:                                    
#                                        for Port_Number in Port_Number_List:                                                
#
#                                            if Dish == 1:
#                                                if BackPlane == 0:
#                                                    if BeamPattern ==1:
#                                                        
#                                                        BeamSinuousDishBandSkirt = BeamSinous_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number, ['XX'],rotateY=True)
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,'o',label='Sinuous_Dish-Band-Skirt')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$\\xi$',fontsize=20)
#                                                        p.legend(loc='best',fontsize=10,ncol=1)
#                                                        p.yscale('log')
#                                                        p.title('FarEllip_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))
#                                                        p.gca().tick_params('x',labelsize=16)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        p.gca().yaxis.grid(which='minor')
#                                                        p.grid()   
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Data_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),'o',label='Sinuous_Dish-Band-Skirt')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
#                                                        p.grid()
#                                                        #p.ylim(.15,.85)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        #p.gca().yaxis.grid(which='minor')
#                                                        p.legend(loc='best',fontsize=10,ncol=1) 
#                                                        p.title('FarEffecArea_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))                           
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Data_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7)],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        for m in range(len(Frequency_List)):                                                                    
#                                                            nth=8000
#                                                            tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                            l1=p.plot(tha,hpCut(Phi_List[1],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k')[0]
#                                                            p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]),bbox_inches='tight')
#                                                            p.close()
#                                                            
#                                                        phi = m = 0
#                                                            
#                                                        for phi in range(len(Phi_List)):
#                                                            for m in range(len(Frequency_List)):
#                                                                if Simplify == 1:                                                        
#                                                                    if m!=0 and m!=5 and m!=11 and m!=12:                                                                                                                
#                                                                        nth=8000
#                                                                        tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                                #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                        p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousDishBandSkirt-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                                #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                                if Simplify == 0:
#                                                                    nth=8000
#                                                                    tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                    p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousDishBandSkirt-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                            #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                                    
#                                                            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                                                                        ncol=2, mode="expand", borderaxespad=1)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,PhiDeg_List[phi]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,PhiDeg_List[phi]),bbox_inches='tight')
#                                                            p.close()
#                                                    
#                                                        
#                                                    fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%s-%s-%s_dish-band-%s-skirt-%s-%s.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
#                                                    fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%s-%s-%s_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
#                                                    fileNameImpedanceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Impedance_0.%s-%s-%s_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
#
#                                                    FLOW=0.05
#                                                    FHIGH=0.25
#
#
#                                                    gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
#                                                                                        fileType='CST_TimeTrace',
#                                                                                        fMin=FLOW,fMax=FHIGH,
#                                                                                        comment='S11 derived from cst time domain data')
#                                                    gainData_cst=gainData.GainData(fileNameS11CST,
#                                                                                    fileType='CST_S11',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='S11 obtained directly from cst')
#                                                    gainData_impedance=gainData.GainData(fileNameImpedanceCST,
#                                                                                    fileType='CST_Impedance',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='Impedance obtained directly from cst')
#
#                                                    #gainData_far=
#                                                    #gainData_vna=gainData.GainData(fileNameS11VNA,
#                                                    #							   fileType='VNAHP_S11',
#                                                    #							   fMin=FLOW,fMax=FHIGH,
#                                                    #							   comment='s11 obtained from richs vna measurement')
#
#                                                    print gainData_cst.gainFrequency.shape
#                                                    print gainData_impedance.gainFrequency.shape
#
#                                                    #first make original plot comparing s11 of time trace and s11 of vna
#
#                                                    #p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.tAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.tAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.xlim(-30,400)
#                                                    p.ylim(-70*S11_Power,0)
#                                                    p.ylabel('|$\widetilde{S}_{11}$|(dB)')
#                                                    p.xlabel('delay (ns)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Delay_0.%s-%s-%s_PW%s_dish-band_%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height))
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%s-%s-%s_PW%s_Cr_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#
#                                                    #p.plot(gainData_vna.fAxis,10.*n.log10(n.abs(gainData_vna.gainFrequency)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.fAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainFrequency)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.fAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainFrequency)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.plot(nu, 5*PW*n.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
#                                                    p.xlim(.045,.255)
#                                                    p.ylim(-25*S11_Power,0)
#                                                    p.ylabel('|S$_{11}$|(dB)')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Frequency_0.%s-%s-%s_PW%s_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%s-%s-%s_PW%s_Cr_dish-band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.abs(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Abs',markersize=4,markeredgecolor='none')
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('Impedance(Abs)/Ohm')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedanceAbs_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedanceAbs_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.angle(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Pha',markersize=4,markeredgecolor='none')                                        
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('|Impedance(Pha)/deg')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedancePha_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedancePha_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#
#                                                
#                                                if BackPlane == 1:
#                                                    if BeamPattern ==1:
#                                                        
#                                                        BeamSinuousDishBandSkirt = BeamSinous_BackPlane_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number, ['XX'],rotateY=True)
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,'o',label='Sinuous_Dish-Band-Skirt-BackPlane')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$\\xi$',fontsize=20)
#                                                        p.legend(loc='best',fontsize=10,ncol=1)
#                                                        p.yscale('log')
#                                                        p.title('FarEllip_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter))
#                                                        p.gca().tick_params('x',labelsize=16)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        p.gca().yaxis.grid(which='minor')
#                                                        p.grid()   
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Data_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),'o',label='Sinuous_Dish-Band-Skirt-BackPlane')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
#                                                        p.grid()
#                                                        #p.ylim(.15,.85)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        #p.gca().yaxis.grid(which='minor')
#                                                        p.legend(loc='best',fontsize=10,ncol=1) 
#                                                        p.title('FarEffecArea_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter))                           
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Data_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7)],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        for m in range(len(Frequency_List)):                                                                    
#                                                            nth=8000
#                                                            tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                            l1=p.plot(tha,hpCut(Phi_List[1],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k')[0]
#                                                            p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-BackPlane-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-BackPlane-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]),bbox_inches='tight')
#                                                            p.close()
#                                                            
#                                                        phi = m = 0
#                                                            
#                                                        for phi in range(len(Phi_List)):
#                                                            for m in range(len(Frequency_List)):
#                                                                if Simplify == 1:                                                        
#                                                                    if m!=0 and m!=5 and m!=11 and m!=12:                                                                                                                
#                                                                        nth=8000
#                                                                        tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                                #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                        p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousDishBandSkirtBackPlane-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                                #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                                if Simplify == 0:                                                        
#                                                                    nth=8000
#                                                                    tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                    p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousDishBandSkirtBackPlane-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                            #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                                                                        ncol=1, mode="expand", borderaxespad=1)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,PhiDeg_List[phi]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,PhiDeg_List[phi]),bbox_inches='tight')
#                                                            p.close()
#                                                    
#                                                    
#                                                    fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
#                                                    fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
#                                                    fileNameImpedanceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Impedance_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
#
#                                                    FLOW=0.05
#                                                    FHIGH=0.25
#
#
#                                                    gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
#                                                                                        fileType='CST_TimeTrace',
#                                                                                        fMin=FLOW,fMax=FHIGH,
#                                                                                        comment='S11 derived from cst time domain data')
#                                                    gainData_cst=gainData.GainData(fileNameS11CST,
#                                                                                    fileType='CST_S11',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='S11 obtained directly from cst')
#                                                    gainData_impedance=gainData.GainData(fileNameImpedanceCST,
#                                                                                    fileType='CST_Impedance',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='Impedance obtained directly from cst')
#                                                    #gainData_far=
#                                                    #gainData_vna=gainData.GainData(fileNameS11VNA,
#                                                    #							   fileType='VNAHP_S11',
#                                                    #							   fMin=FLOW,fMax=FHIGH,
#                                                    #							   comment='s11 obtained from richs vna measurement')
#
#                                                    print gainData_cst.gainFrequency.shape
#                                                    print gainData_impedance.gainFrequency.shape
#
#                                                    #first make original plot comparing s11 of time trace and s11 of vna
#
#                                                    #p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.tAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.tAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.xlim(-30,400)
#                                                    p.ylim(-70*S11_Power,0)
#                                                    p.ylabel('|$\widetilde{S}_{11}$|(dB)')
#                                                    p.xlabel('delay (ns)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Delay_0.%s-%s-%s_PW%s_dish-band_%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter))
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%s-%s-%s_PW%s_Cr_dish-band_%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
#
#                                                    #p.plot(gainData_vna.fAxis,10.*n.log10(n.abs(gainData_vna.gainFrequency)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.fAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainFrequency)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.fAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainFrequency)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.plot(nu, 5*PW*n.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
#                                                    p.xlim(.045,.255)
#                                                    p.ylim(-25*S11_Power,0)
#                                                    p.ylabel('|S$_{11}$|(dB)')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Frequency_0.%s-%s-%s_PW%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%s-%s-%s_PW%s_Cr_dish-band-%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.abs(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Abs',markersize=4,markeredgecolor='none')
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('|Impedance(Abs)/Ohm')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedanceAbs_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedanceAbs_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.angle(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Pha',markersize=4,markeredgecolor='none')
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('|Impedance(Pha)/deg')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedancePha_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedancePha_CST_Frequency_0.%s-%s-%s_dish-band-%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
#
#                                            if Dish == 0:
#                                                if BackPlane == 0:
#                                                    if BeamPattern ==1:
#                                                        
#                                                        BeamSinuousDishBandSkirt = BeamSinous_NoDish_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number, ['XX'],rotateY=True)
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,'o',label='Sinuous_Band-Skirt')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$\\xi$',fontsize=20)
#                                                        p.legend(loc='best',fontsize=10,ncol=1)
#                                                        p.yscale('log')
#                                                        p.title('FarEllip_Y_0.%s-%s-%s_band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))
#                                                        p.gca().tick_params('x',labelsize=16)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        p.gca().yaxis.grid(which='minor')
#                                                        p.grid()   
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Y_0.%s-%s-%s_band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Data_Y_0.%s-%s-%s_band_%s-skirt-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),'o',label='Sinuous_Band-Skirt')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
#                                                        p.grid()
#                                                        #p.ylim(.15,.85)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        #p.gca().yaxis.grid(which='minor')
#                                                        p.legend(loc='best',fontsize=10,ncol=1) 
#                                                        p.title('FarEffecArea_Y_0.%s-%s-%s_band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))                           
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Y_0.%s-%s-%s_band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Data_Y_0.%s-%s-%s_band_%s-skirt-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7)],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        for m in range(len(Frequency_List)):                                                                    
#                                                            nth=8000
#                                                            tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                            l1=p.plot(tha,hpCut(Phi_List[1],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k')[0]
#                                                            p.gcf().legend((l,l1),('Sinuous_Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]),bbox_inches='tight')
#                                                            p.close()
#                                                            
#                                                        phi = m = 0
#                                                            
#                                                        for phi in range(len(Phi_List)):
#                                                            for m in range(len(Frequency_List)):
#                                                                if Simplify == 1:                                                        
#                                                                    if m!=0 and m!=5 and m!=11 and m!=12:                                                                                                                
#                                                                        nth=8000
#                                                                        tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                                #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                        p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousBandSkirt-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                                #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                                if Simplify == 0:
#                                                                    nth=8000
#                                                                    tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                    p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousBandSkirt-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                            #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                                    
#                                                            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                                                                        ncol=2, mode="expand", borderaxespad=1)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,PhiDeg_List[phi]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,PhiDeg_List[phi]),bbox_inches='tight')
#                                                            p.close()
#                                                    
#                                                        
#                                                    fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%s-%s-%s_band-%s-skirt-%s-%s.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
#                                                    fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%s-%s-%s_band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
#                                                    fileNameImpedanceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Impedance_0.%s-%s-%s_band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
#
#                                                    FLOW=0.05
#                                                    FHIGH=0.25
#
#
#                                                    gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
#                                                                                        fileType='CST_TimeTrace',
#                                                                                        fMin=FLOW,fMax=FHIGH,
#                                                                                        comment='S11 derived from cst time domain data')
#                                                    gainData_cst=gainData.GainData(fileNameS11CST,
#                                                                                    fileType='CST_S11',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='S11 obtained directly from cst')
#                                                    gainData_impedance=gainData.GainData(fileNameImpedanceCST,
#                                                                                    fileType='CST_Impedance',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='Impedance obtained directly from cst')
#
#                                                    #gainData_far=
#                                                    #gainData_vna=gainData.GainData(fileNameS11VNA,
#                                                    #							   fileType='VNAHP_S11',
#                                                    #							   fMin=FLOW,fMax=FHIGH,
#                                                    #							   comment='s11 obtained from richs vna measurement')
#
#                                                    print gainData_cst.gainFrequency.shape
#                                                    print gainData_impedance.gainFrequency.shape
#
#                                                    #first make original plot comparing s11 of time trace and s11 of vna
#
#                                                    #p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.tAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.tAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.xlim(-30,400)
#                                                    p.ylim(-70*S11_Power,0)
#                                                    p.ylabel('|$\widetilde{S}_{11}$|(dB)')
#                                                    p.xlabel('delay (ns)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Delay_0.%s-%s-%s_PW%s_band_%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height))
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%s-%s-%s_PW%s_Cr_band_%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#
#                                                    #p.plot(gainData_vna.fAxis,10.*n.log10(n.abs(gainData_vna.gainFrequency)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.fAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainFrequency)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.fAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainFrequency)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.plot(nu, 5*PW*n.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
#                                                    p.xlim(.045,.255)
#                                                    p.ylim(-25*S11_Power,0)
#                                                    p.ylabel('|S$_{11}$|(dB)')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Frequency_0.%s-%s-%s_PW%s_band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%s-%s-%s_PW%s_Cr_band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.abs(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Abs',markersize=4,markeredgecolor='none')
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('Impedance(Abs)/Ohm')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedanceAbs_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedanceAbs_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.angle(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Pha',markersize=4,markeredgecolor='none')                                        
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('|Impedance(Pha)/deg')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedancePha_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedancePha_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                                                    p.close()
#
#                                                
#                                                if BackPlane == 1:
#                                                    if BeamPattern ==1:
#                                                        
#                                                        BeamSinuousDishBandSkirt = BeamSinous_NoDish_BackPlane_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number, ['XX'],rotateY=True)
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,'o',label='Sinuous_Band-Skirt-BackPlane')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$\\xi$',fontsize=20)
#                                                        p.legend(loc='best',fontsize=10,ncol=1)
#                                                        p.yscale('log')
#                                                        p.title('FarEllip_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter))
#                                                        p.gca().tick_params('x',labelsize=16)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        p.gca().yaxis.grid(which='minor')
#                                                        p.grid()   
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Data_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),'o',label='Sinuous_Band-Skirt-BackPlane')
#                                                        p.xlabel('f (MHz)',fontsize=20)
#                                                        p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
#                                                        p.grid()
#                                                        #p.ylim(.15,.85)
#                                                        p.gca().tick_params('y',labelsize=16)
#                                                        p.gcf().set_size_inches([8,6])
#                                                        #p.gca().yaxis.grid(which='minor')
#                                                        p.legend(loc='best',fontsize=10,ncol=1) 
#                                                        p.title('FarEffecArea_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter))                           
#                                                        p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                        p.close()
#                                                        
#                                                        n.savetxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEffecArea_Data_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s.txt'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),n.c_[BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7)],fmt=['%.2f','%.4f']) 
#                                                        
#                                                        for m in range(len(Frequency_List)):                                                                    
#                                                            nth=8000
#                                                            tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                            l1=p.plot(tha,hpCut(Phi_List[1],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k')[0]
#                                                            p.gcf().legend((l,l1),('Sinuous_Band-Skirt-BackPlane-%s'%PhiDeg_List[0],'Sinuous_Band-Skirt-BackPlane-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]),bbox_inches='tight')
#                                                            p.close()
#                                                            
#                                                        phi = m = 0
#                                                            
#                                                        for phi in range(len(Phi_List)):
#                                                            for m in range(len(Frequency_List)):
#                                                                if Simplify == 1:                                                        
#                                                                    if m!=0 and m!=5 and m!=11 and m!=12:                                                                                                                
#                                                                        nth=8000
#                                                                        tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                                #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                        p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousBandSkirtBackPlane-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                                #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                                if Simplify == 0:                                                        
#                                                                    nth=8000
#                                                                    tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                                            #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                                                    p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousBandSkirtBackPlane-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
#                                                            #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
#                                                            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                                                                        ncol=1, mode="expand", borderaxespad=1)
#                                                            p.grid()
#                                                            p.xlabel('$\\theta$')
#                                                            p.xlim(-90,90)
#                                                            p.ylim(-30,30)
#                                                            p.ylabel('Directivity (dB)')
#                                                            p.gcf().set_size_inches([10,7])
#                                                            p.title('FarCut_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,PhiDeg_List[phi]))
#                                                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarField_Y_0.%s-%s-%s_band_%s-skirt-%s-%s-backplane-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter,PhiDeg_List[phi]),bbox_inches='tight')
#                                                            p.close()
#                                                    
#                                                    
#                                                    fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%s-%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
#                                                    fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%s-%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
#                                                    fileNameImpedanceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Impedance_0.%s-%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)
#
#                                                    FLOW=0.05
#                                                    FHIGH=0.25
#
#
#                                                    gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
#                                                                                        fileType='CST_TimeTrace',
#                                                                                        fMin=FLOW,fMax=FHIGH,
#                                                                                        comment='S11 derived from cst time domain data')
#                                                    gainData_cst=gainData.GainData(fileNameS11CST,
#                                                                                    fileType='CST_S11',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='S11 obtained directly from cst')
#                                                    gainData_impedance=gainData.GainData(fileNameImpedanceCST,
#                                                                                    fileType='CST_Impedance',
#                                                                                    fMin=FLOW,fMax=FHIGH,
#                                                                                    comment='Impedance obtained directly from cst')
#                                                    #gainData_far=
#                                                    #gainData_vna=gainData.GainData(fileNameS11VNA,
#                                                    #							   fileType='VNAHP_S11',
#                                                    #							   fMin=FLOW,fMax=FHIGH,
#                                                    #							   comment='s11 obtained from richs vna measurement')
#
#                                                    print gainData_cst.gainFrequency.shape
#                                                    print gainData_impedance.gainFrequency.shape
#
#                                                    #first make original plot comparing s11 of time trace and s11 of vna
#
#                                                    #p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.tAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.tAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.xlim(-30,400)
#                                                    p.ylim(-70*S11_Power,0)
#                                                    p.ylabel('|$\widetilde{S}_{11}$|(dB)')
#                                                    p.xlabel('delay (ns)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Delay_0.%s-%s-%s_PW%s_band_%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter))
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%s-%s-%s_PW%s_Cr_band_%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
#
#                                                    #p.plot(gainData_vna.fAxis,10.*n.log10(n.abs(gainData_vna.gainFrequency)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_timeTrace.fAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainFrequency)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                                                    p.plot(gainData_cst.fAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainFrequency)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                                                    p.plot(nu, 5*PW*n.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
#                                                    p.xlim(.045,.255)
#                                                    p.ylim(-25*S11_Power,0)
#                                                    p.ylabel('|S$_{11}$|(dB)')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('S11_CST_Frequency_0.%s-%s-%s_PW%s_band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%s-%s-%s_PW%s_Cr_band-%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.abs(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Abs',markersize=4,markeredgecolor='none')
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('|Impedance(Abs)/Ohm')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedanceAbs_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedanceAbs_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
#                                                    
#                                                    p.plot(gainData_impedance.fAxis,n.angle(gainData_impedance.gainFrequency),color='k',ls='-',marker='o',label='CST Impedance Pha',markersize=4,markeredgecolor='none')
#                                                    p.xlim(.045,.255)
#                                                    p.ylabel('|Impedance(Pha)/deg')
#                                                    p.xlabel('f (GHz)')
#                                                    p.legend(loc='best')
#                                                    p.title('ImpedancePha_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter)) 
#                                                    #p.show()
#                                                    p.grid()
#                                                    #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/ImpedancePha_CST_Frequency_0.%s-%s-%s_band-%s-skirt-%s-%s-backplane-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,BackPlane_Height,BackPlane_Diameter),bbox_inches='tight')
#                                                    p.close()
