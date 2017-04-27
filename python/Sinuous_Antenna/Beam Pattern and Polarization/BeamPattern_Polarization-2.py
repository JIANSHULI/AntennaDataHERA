
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


# In[2]:

#fList=n.arange(50,210,10)
#Growth_Rate_List = [50]
#Outer_Diameter_List = [275]
#Inner_Diameter_List = []
#Band_Resistance_List = [50]
#Skirt_Diameter_List = [1.2]
#Skirt_Height_List = [0.3]
#Frequency_List = [40,50,67,94,121,148,150,175,202,229,250,256,283]
#Port_Number_List = [1]


# In[3]:

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
#        for m in range(self.nf):            
#            print m
#            tempf=fList[m].split('p')
#            self.fAxis[m]=float(tempf[0])*1e6
#            if(len(tempf)>1):
#                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
#            for np in range(self.npolsOriginal):
#                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far_%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
#                self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
#                self.data[np,m,theta>90.]=0.
#                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
#                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
        for m in range(self.nf):            
            print m
            tempf=fList[m].split('p')
            self.fAxis[m]=float(tempf[0])*1e6
            if(len(tempf)>1):
                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
            for np in range(self.npolsOriginal):
                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                                            data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
                #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
                self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
                self.data_N[np,m,:]=self.data[np,m,:]/self.data[np,m,:].flatten().max(); 
                self.data[np,m,theta>90.]=0.
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


# In[4]:

#beamCylinder=Beam('beamCylinder',['100'],64)
#PW = S11_Power = 2
N = 0

Growth_Rate_List = [49]
Outer_Diameter_List = [225]
Inner_Diameter_List = [30]
Band_Resistance_List = [15]
Skirt_Diameter_List = [1.2]
Skirt_Height_List = [0.3]
Frequency_List = ['40','50','67','94','121','148','150','175','202','229','250','256','283']    #[40,50,67,94,121,148,150,175,202,229,250,256,283]
Port_Number_List = [1]
Phi_List = [0,pi/2.0]
PhiDeg_List = [0,90]

# Y-Direction
for N in range(2):
    PW = S11_Power = N+1                             
    for Growth_Rate in Growth_Rate_List:
        for Outer_Diameter in Outer_Diameter_List:
            for Inner_Diameter in Inner_Diameter_List:                    
                for Band_Resistance in Band_Resistance_List:
                    for Skirt_Diameter in Skirt_Diameter_List:
                        for Skirt_Height in Skirt_Height_List:
                            for Port_Number in Port_Number_List:
                                BeamSinuousDishBandSkirt = BeamSinous_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number, ['XX'],rotateY=True)
                                
                                p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,label='Sinuous_Dish-Band-Skirt',color='k',ls='-',lw=2)
                                p.xlabel('f (MHz)',fontsize=20)
                                p.ylabel('$\\xi$',fontsize=20)
                                p.legend(loc='best',fontsize=10,ncol=1)
                                p.yscale('log')
                                p.title('FarEllip_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))
                                p.gca().tick_params('x',labelsize=16)
                                p.gca().tick_params('y',labelsize=16)
                                p.gcf().set_size_inches([8,6])
                                p.gca().yaxis.grid(which='minor')
                                p.grid()   
                                p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
                                p.close() 
                                
                                p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),label='Sinuous_Dish-Band-Skirt',color='k',ls='-',lw=2)
                                p.xlabel('f (MHz)',fontsize=20)
                                p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
                                p.grid()
                                #p.ylim(.15,.85)
                                p.gca().tick_params('y',labelsize=16)
                                p.gcf().set_size_inches([8,6])
                                #p.gca().yaxis.grid(which='minor')
                                p.legend(loc='best',fontsize=10,ncol=1) 
                                p.title('FarArea_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))                           
                                p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarArea_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
                                p.close()
                                
                                for m in range(len(Frequency_List)):                                                                    
                                    nth=8000
                                    tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
                                    l1=p.plot(tha,hpCut(Phi_List[1],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k')[0]
                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
                                    p.grid()
                                    p.xlabel('$\\theta$')
                                    p.xlim(-90,90)
                                    p.ylim(-30,30)
                                    p.ylabel('Directivity (dB)')
                                    p.gcf().set_size_inches([10,7])
                                    p.title('FarCut_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]))
                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarArea_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]),bbox_inches='tight')
                                    p.close()
                                    
                                phi = m = 0
                                    
                                for phi in range(len(Phi_List)):
                                    for m in range(len(Frequency_List)):
                                        if m!=0 and m!=5 and m!=11 and m!=12:                                                                                                                
                                            nth=8000
                                            tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
                                    #                                    l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
                                            p.plot(tha,hpCut(Phi_List[phi],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),label='SinuousDishBandSkirt-%s-%s' %(PhiDeg_List[phi],Frequency_List[m]))
                                    #                                    p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[0],'Sinuous_Dish-Band-Skirt-%s'%PhiDeg_List[1]),loc='upper center',ncol=2)
                                    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                                               ncol=2, mode="expand", borderaxespad=1)
                                    p.grid()
                                    p.xlabel('$\\theta$')
                                    p.xlim(-90,90)
                                    p.ylim(-30,30)
                                    p.ylabel('Directivity (dB)')
                                    p.gcf().set_size_inches([10,7])
                                    p.title('FarCut_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,PhiDeg_List[phi]))
                                    p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarArea_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,PhiDeg_List[phi]),bbox_inches='tight')
                                    p.close()
                                    
                                fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%i-%i-%i_dish-band-%s-skirt-%s-%s.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
                                fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%i-%i-%i_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
    #                            fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Far_0.%i-%i_dish-band-%s-Skirt-%s-%s-%s-%s.txt' %(Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequecy,Port_Number)
                                #fileNameS11VNA='../reflectometry/RichBradley_GreenBank/TallCylinderGapOverDish_S11_Greenbank_RichBradley.d1'

                                FLOW=0.05
                                FHIGH=0.25


                                gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
                                                                    fileType='CST_TimeTrace',
                                                                    fMin=FLOW,fMax=FHIGH,
                                                                    comment='s11 derived from cst time domain data')
                                gainData_cst=gainData.GainData(fileNameS11CST,
                                                                fileType='CST_S11',
                                                                fMin=FLOW,fMax=FHIGH,
                                                                comment='s11 obtained directly from cst')
                                #gainData_far=
                                #gainData_vna=gainData.GainData(fileNameS11VNA,
                                #							   fileType='VNAHP_S11',
                                #							   fMin=FLOW,fMax=FHIGH,
                                #							   comment='s11 obtained from richs vna measurement')

                                print gainData_cst.gainFrequency.shape

                                #first make original plot comparing s11 of time trace and s11 of vna

                                #p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
                                p.plot(gainData_timeTrace.tAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
                                p.plot(gainData_cst.tAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
                                p.xlim(-30,400)
                                p.ylim(-70*S11_Power,0)
                                p.ylabel('|$\widetilde{S}_{11}$|(dB)')
                                p.xlabel('delay (ns)')
                                p.legend(loc='best')
                                p.title('S11_CST_Delay_0.%i-%i-%i_PW%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height))
                                #p.show()
                                p.grid()
                                #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
                                p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%i-%i-%i_PW%i_Cr_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
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
                                p.title('S11_CST_Frequency_0.%i-%i-%i_PW%i_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
                                #p.show()
                                p.grid()
                                #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
                                p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%i-%i-%i_PW%i_Cr_dish-band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
                                p.close()                                                                                            

#for m in range(7,8):
#    hp.mollview(10*n.log10(BeamSinuousDishBandSkirt.data_N[0,m,:]),rot=(0,90,0),min=-40,max=-0.0001)
#     
#for m in range(7,8):
#    hp.mollview(10*n.log10(BeamSinuousDishBandSkirt.data_N[0,m,:]),rot=(0,90,0))


#class Beam:
#    def __init__(self,dirName,fList,nside,pols=['XX','YY'],rotateY=False):
#        self.nf=len(fList)
#        self.fAxis=n.zeros(self.nf)
#        self.npolsOriginal=len(pols)
#        self.npols=max(len(pols),2)
#        self.solidAngles=n.zeros((self.npols,self.nf))
#        self.effArea=n.zeros_like(self.solidAngles)
#        self.ellipticity=n.zeros(self.nf)
#        self.nPix=hp.nside2npix(nside)
#        self.nSide=nside
#        self.pixArea=hp.nside2pixarea(self.nSide)
#        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
#        theta=n.round(n.degrees(theta)).astype(int)
#        phi=n.round(n.degrees(phi)).astype(int)
#        self.data=n.zeros((self.npols,self.nf,self.nPix))
#        if(rotateY):
#            pols.append('YY')
#        self.pols=pols
#        for m in range(self.nf):            
#            print m
#            tempf=fList[m].split('p')
#            self.fAxis[m]=float(tempf[0])*1e6
#            if(len(tempf)>1):
#                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
#            for np in range(self.npolsOriginal):
#                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far_%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
#                self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
#                self.data[np,m,theta>90.]=0.
#                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
#                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
##        for Growth_Rate in Growth_Rate_List:
##            for Outer_Diameter in Outer_Diameter_List:
##                for Band_Resistance in Band_Resistance_List:
##                    for Skirt_Diameter in Skirt_Diameter_List:
##                        for Skirt_Height in Skirt_Height_List:
##                            for m in range(self.nf):
##                                for Port_Number in Port_Number_List:
##                                    print m
##                                    tempf=fList[m].split('p')
##                                    self.fAxis[m]=float(tempf[0])*1e6
##                                    if(len(tempf)>1):
##                                        self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
##                                    for np in range(self.npolsOriginal):
##                                        #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                        data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequecy,Port_Number),skiprows=2);
##                                        #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                        self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
##                                        self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
##                                        self.data[np,m,theta>90.]=0.
##                                        self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
##                                        self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
#                                        
#            if(self.npolsOriginal==1):
#                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
#                self.solidAngles[1,m]=self.pixArea*n.sum(self.data[1,m,:])
#                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
#            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
#                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)


#class Beam:
#    def __init__(self,dirName,fList,nside,Skirt_Diameter_List,Skirt_Height_List,Growth_Rate_List,Outer_Diameter_List,Inner_Diameter_List,Band_Resistance_List,Port_Number_List,pols=['XX','YY'],rotateY=False):
#        self.nf=len(fList)
#        self.fAxis=n.zeros(self.nf)
#        self.npolsOriginal=len(pols)
#        self.npols=max(len(pols),2)
#        self.solidAngles=n.zeros((self.npols,self.nf))
#        self.effArea=n.zeros_like(self.solidAngles)
#        self.ellipticity=n.zeros(self.nf)
#        self.nPix=hp.nside2npix(nside)
#        self.nSide=nside
#        self.pixArea=hp.nside2pixarea(self.nSide)
#        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
#        theta=n.round(n.degrees(theta)).astype(int)
#        phi=n.round(n.degrees(phi)).astype(int)
#        self.data=n.zeros((self.npols,self.nf,self.nPix))
#        if(rotateY):
#            pols.append('YY')
#        self.pols=pols
#        for m in range(self.nf):            
#            print m
#            tempf=fList[m].split('p')
#            self.fAxis[m]=float(tempf[0])*1e6
#            if(len(tempf)>1):
#                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
#            for np in range(self.npolsOriginal):
#                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far_%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
#                #self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
#                self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
#                self.data[np,m,theta>90.]=0.
#                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
#                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
##        for Growth_Rate in Growth_Rate_List:
##            for Outer_Diameter in Outer_Diameter_List:
##                for Band_Resistance in Band_Resistance_List:
##                    for Skirt_Diameter in Skirt_Diameter_List:
##                        for Skirt_Height in Skirt_Height_List:
##                            for m in range(self.nf):
##                                for Port_Number in Port_Number_List:
##                                    print m
##                                    tempf=fList[m].split('p')
##                                    self.fAxis[m]=float(tempf[0])*1e6
##                                    if(len(tempf)>1):
##                                        self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
##                                    for np in range(self.npolsOriginal):
##                                        #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                        data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequecy,Port_Number),skiprows=2);
##                                        #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                        self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
##                                        self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
##                                        self.data[np,m,theta>90.]=0.
##                                        self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
##                                        self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
#                                        
#            if(self.npolsOriginal==1):
#                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
#                self.solidAngles[1,m]=self.pixArea*n.sum(self.data[1,m,:])
#                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
#            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
#                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)
#                
#class Beam_Y:
#    def __init__(self,dirName,fList,nside,Skirt_Diameter_List,Skirt_Height_List,Growth_Rate_List,Outer_Diameter_List,Inner_Diameter_List,Band_Resistance_List,Port_Number_List,pols=['XX','YY'],rotateY=False):
#        self.nf=len(fList)
#        self.fAxis=n.zeros(self.nf)
#        self.npolsOriginal=len(pols)
#        self.npols=max(len(pols),2)
#        self.solidAngles=n.zeros((self.npols,self.nf))
#        self.effArea=n.zeros_like(self.solidAngles)
#        self.ellipticity=n.zeros(self.nf)
#        self.nPix=hp.nside2npix(nside)
#        self.nSide=nside
#        self.pixArea=hp.nside2pixarea(self.nSide)
#        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
#        theta=n.round(n.degrees(theta)).astype(int)
#        phi=n.round(n.degrees(phi)).astype(int)
#        self.data=n.zeros((self.npols,self.nf,self.nPix))
#        if(rotateY):
#            pols.append('YY')
#        self.pols=pols
#        for m in range(self.nf):            
#            print m
#            tempf=fList[m].split('p')
#            self.fAxis[m]=float(tempf[0])*1e6
#            if(len(tempf)>1):
#                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
#            for np in range(self.npolsOriginal):
#                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far_%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
#                self.data[np,m,:]=rotateBeam_to_Z(self.data[np,m,:].flatten())
#                self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
#                self.data[np,m,theta>90.]=0.
#                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
#                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
##        for Growth_Rate in Growth_Rate_List:
##            for Outer_Diameter in Outer_Diameter_List:
##                for Band_Resistance in Band_Resistance_List:
##                    for Skirt_Diameter in Skirt_Diameter_List:
##                        for Skirt_Height in Skirt_Height_List:
##                            for m in range(self.nf):
##                                for Port_Number in Port_Number_List:
##                                    print m
##                                    tempf=fList[m].split('p')
##                                    self.fAxis[m]=float(tempf[0])*1e6
##                                    if(len(tempf)>1):
##                                        self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
##                                    for np in range(self.npolsOriginal):
##                                        #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                        data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequecy,Port_Number),skiprows=2);
##                                        #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                        self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
##                                        self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
##                                        self.data[np,m,theta>90.]=0.
##                                        self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
##                                        self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
#                                        
#            if(self.npolsOriginal==1):
#                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
#                self.solidAngles[1,m]=self.pixArea*n.sum(self.data[1,m,:])
#                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
#            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
#                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)
#                
#class BeamSinous:
#    def __init__(self,dirName,fList,nside,Skirt_Diameter_List,Skirt_Height_List,Growth_Rate_List,Outer_Diameter_List,Inner_Diameter_List,Band_Resistance_List,Port_Number_List,pols=['XX','YY'],rotateY=False):
#        self.nf=len(fList)
#        self.fAxis=n.zeros(self.nf)
#        self.npolsOriginal=len(pols)
#        self.npols=max(len(pols),2)
#        self.solidAngles=n.zeros((self.npols,self.nf))
#        self.effArea=n.zeros_like(self.solidAngles)
#        self.ellipticity=n.zeros(self.nf)
#        self.nPix=hp.nside2npix(nside)
#        self.nSide=nside
#        self.pixArea=hp.nside2pixarea(self.nSide)
#        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
#        theta=n.round(n.degrees(theta)).astype(int)
#        phi=n.round(n.degrees(phi)).astype(int)
#        self.data=n.zeros((self.npols,self.nf,self.nPix))
#        if(rotateY):
#            pols.append('YY')
#        self.pols=pols
##        for m in range(self.nf):            
##            print m
##            tempf=fList[m].split('p')
##            self.fAxis[m]=float(tempf[0])*1e6
##            if(len(tempf)>1):
##                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
##            for np in range(self.npolsOriginal):
##                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far_%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
##                self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
##                self.data[np,m,theta>90.]=0.
##                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
##                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
#        for Growth_Rate in Growth_Rate_List:
#            for Outer_Diameter in Outer_Diameter_List:
#                for Inner_Diameter in Inner_Diameter_List:                    
#                    for Band_Resistance in Band_Resistance_List:
#                        for Skirt_Diameter in Skirt_Diameter_List:
#                            for Skirt_Height in Skirt_Height_List:
#                                for m in range(self.nf):
#                                    for Port_Number in Port_Number_List:
#                                        print m
#                                        tempf=fList[m].split('p')
#                                        self.fAxis[m]=float(tempf[0])*1e6
#                                        if(len(tempf)>1):
#                                            self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
#                                        for np in range(self.npolsOriginal):
#                                            #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                            data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
#                                            data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
#                                            #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                                            self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
#                                            self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
#                                            self.data[np,m,theta>90.]=0.
#                                            self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
#                                            self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
#                                            
#                                            
#                if(self.npolsOriginal==1):
#                    self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
#                    self.solidAngles[1,m]=self.pixArea*n.sum(self.data[1,m,:])
#                    self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
#                if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
#                    self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)

#class BeamSinous:
#    def __init__(self,dirName,fList,nside,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number,pols=['XX','YY'],rotateY=False):
#        self.nf=len(fList)
#        self.fAxis=n.zeros(self.nf)
#        self.npolsOriginal=len(pols)
#        self.npols=max(len(pols),2)
#        self.solidAngles=n.zeros((self.npols,self.nf))
#        self.effArea=n.zeros_like(self.solidAngles)
#        self.ellipticity=n.zeros(self.nf)
#        self.nPix=hp.nside2npix(nside)
#        self.nSide=nside
#        self.pixArea=hp.nside2pixarea(self.nSide)
#        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
#        theta=n.round(n.degrees(theta)).astype(int)
#        phi=n.round(n.degrees(phi)).astype(int)
#        self.data=n.zeros((self.npols,self.nf,self.nPix))
#        if(rotateY):
#            pols.append('YY')
#        self.pols=pols
##        for m in range(self.nf):            
##            print m
##            tempf=fList[m].split('p')
##            self.fAxis[m]=float(tempf[0])*1e6
##            if(len(tempf)>1):
##                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
##            for np in range(self.npolsOriginal):
##                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far_%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
##                self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
##                self.data[np,m,theta>90.]=0.
##                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
##                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
#        for m in range(self.nf):            
#            print m
#            tempf=fList[m].split('p')
#            self.fAxis[m]=float(tempf[0])*1e6
#            if(len(tempf)>1):
#                self.fAxis[m]+=float(tempf[1])/10.**(len(tempf[1]))*1e6
#            for np in range(self.npolsOriginal):
#                #data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
##                                            data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
#                data=n.loadtxt('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/%s/Far-0.%i-%i_dish-band-%s-skirt-%s-%s-%s-%s.txt' %(dirName,Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,fList[m],Port_Number),skiprows=2);
#                #%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
#                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
#                #self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
#                self.data[np,m,theta>90.]=0.
#                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
#                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
#                                            
#                                            
#            if(self.npolsOriginal==1):
#                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
#                self.solidAngles[1,m]=self.pixArea*n.sum(self.data[1,m,:])
#                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
#            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
#                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)                
#

#BeamSinuousDishBandSkirt = BeamSinous_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter_List,Skirt_Height_List,Growth_Rate_List,Outer_Diameter_List,Inner_Diameter_List,Band_Resistance_List,Port_Number_List, ['XX'],rotateY=True)     

#fstrList=['050','060','070','080','090','100','110','120','130','140','150']
#beamCylinder=Beam('beamCylinderBackPlane',fstrList,64,['XX'],rotateY=True)
#beamNoCylinder=Beam('beamBackPlane',fstrList,64,['XX'],rotateY=True)
#beamBareDipole=Beam('beamBareDipole',fstrList,64,['XX'],rotateY=True)
#beamPaneledCylinder=Beam('beamPaneledCylinder',fstrList,64,['XX'],rotateY=True)
#beamCylinderV2=Beam('beamCylinderBackPlane_v2',fstrList,64,['XX'],rotateY=True)

# Z-Direction
#for Growth_Rate in Growth_Rate_List:
#    for Outer_Diameter in Outer_Diameter_List:
#        for Inner_Diameter in Inner_Diameter_List:                    
#            for Band_Resistance in Band_Resistance_List:
#                for Skirt_Diameter in Skirt_Diameter_List:
#                    for Skirt_Height in Skirt_Height_List:
#                        for Port_Number in Port_Number_List:
#                            BeamSinuousDishBandSkirt = BeamSinous('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number, ['XX'],rotateY=True)
#                            
#                            p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,label='Sinuous_Dish-Band-Skirt',color='k',ls='-',lw=2)
#                            p.xlabel('f (MHz)',fontsize=20)
#                            p.ylabel('$\\xi$',fontsize=20)
#                            p.legend(loc='best',fontsize=10,ncol=1)
#                            p.yscale('log')
#                            p.title('FarEllip_Z_0.%i-%i-%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))
#                            p.gca().tick_params('x',labelsize=16)
#                            p.gca().tick_params('y',labelsize=16)
#                            p.gcf().set_size_inches([8,6])
#                            p.gca().yaxis.grid(which='minor')
#                            p.grid()   
#                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarEllip_Z_0.%i-%i-%i_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                            p.close() 
#                            
#                            p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),label='Sinuous_Dish-Band-Skirt',color='k',ls='-',lw=2)
#                            p.xlabel('f (MHz)',fontsize=20)
#                            p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
#                            p.grid()
#                            #p.ylim(.15,.85)
#                            p.gca().tick_params('y',labelsize=16)
#                            p.gcf().set_size_inches([8,6])
#                            #p.gca().yaxis.grid(which='minor')
#                            p.legend(loc='best',fontsize=10,ncol=1) 
#                            p.title('FarArea_Z_0.%i-%i-%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))                             
#                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarArea_Z_0.%i-%i-%i_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                            p.close()
#                            
#                            for m in range(len(Frequency_List)):                                                                    
#                                nth=800
#                                tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#                                l=p.plot(tha,hpCut(Phi_List[0],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k',ls='--')[0]
#                                l1=p.plot(tha,hpCut(Phi_List[1],nth,10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:])),color='k')[0]
#                                p.gcf().legend((l,l1),('Sinuous_Dish-Band-Skirt-%s'%(PhiDeg_List[0]),'Sinuous_Dish-Band-Skirt-%s'%(PhiDeg_List[1])),loc='upper center',ncol=2)
#                                p.grid()
#                                p.xlabel('$\\theta$')
#                                p.xlim(-90,90)
#                                p.ylim(-30,30)
#                                p.ylabel('Directivity (dB)')
#                                p.gcf().set_size_inches([10,7])
#                                p.title('FarCut_Z_0.%i-%i-%i_dish-band_%s-skirt-%s-%s-%s-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]))
#                                p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarArea_Z_0.%i-%i-%i_dish-band_%s-skirt-%s-%s-%s-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequency_List[m],PhiDeg_List[0],PhiDeg_List[1]),bbox_inches='tight') 
#                                p.close()    
#                                
#                            fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%i-%i-%i_dish-band-%s-skirt-%s-%s.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
#                            fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%i-%i-%i_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height)
##                            fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Far_0.%i-%i_dish-band-%s-Skirt-%s-%s-%s-%s.txt' %(Growth_Rate,Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height,Frequecy,Port_Number)
#                            #fileNameS11VNA='../reflectometry/RichBradley_GreenBank/TallCylinderGapOverDish_S11_Greenbank_RichBradley.d1'
#
#                            FLOW=0.05
#                            FHIGH=0.25
#
#
#                            gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
#                                                                fileType='CST_TimeTrace',
#                                                                fMin=FLOW,fMax=FHIGH,
#                                                                comment='s11 derived from cst time domain data')
#                            gainData_cst=gainData.GainData(fileNameS11CST,
#                                                            fileType='CST_S11',
#                                                            fMin=FLOW,fMax=FHIGH,
#                                                            comment='s11 obtained directly from cst')
#                            #gainData_far=
#                            #gainData_vna=gainData.GainData(fileNameS11VNA,
#                            #							   fileType='VNAHP_S11',
#                            #							   fMin=FLOW,fMax=FHIGH,
#                            #							   comment='s11 obtained from richs vna measurement')
#
#                            print gainData_cst.gainFrequency.shape
#
#                            #first make original plot comparing s11 of time trace and s11 of vna
#
#                            #p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                            p.plot(gainData_timeTrace.tAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                            p.plot(gainData_cst.tAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                            p.xlim(-30,400)
#                            p.ylim(-70*S11_Power,0)
#                            p.ylabel('|$\widetilde{S}_{11}$|(dB)')
#                            p.xlabel('delay (ns)')
#                            p.legend(loc='best')
#                            p.title('S11_CST_Delay_0.%i-%i-%i_PW%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height))
#                            #p.show()
#                            p.grid()
#                            #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
#                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%i-%i-%i_PW%i_Cr_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                            p.close()
#
#                            #p.plot(gainData_vna.fAxis,10.*n.log10(n.abs(gainData_vna.gainFrequency)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
#                            p.plot(gainData_timeTrace.fAxis,10.*S11_Power*n.log10(n.abs(gainData_timeTrace.gainFrequency)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
#                            p.plot(gainData_cst.fAxis,10.*S11_Power*n.log10(n.abs(gainData_cst.gainFrequency)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
#                            p.plot(nu, 5*PW*n.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
#                            p.xlim(.045,.255)
#                            p.ylim(-25*S11_Power,0)
#                            p.ylabel('|S$_{11}$|(dB)')
#                            p.xlabel('f (GHz)')
#                            p.legend(loc='best')
#                            p.title('S11_CST_Frequency_0.%i-%i-%i_PW%i_dish-band-%s-skirt-%s-%s' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height)) 
#                            #p.show()
#                            p.grid()
#                            #p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
#                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%i-%i-%i_PW%i_Cr_dish-band-%s-skirt-%s-%s.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                            p.close()                                                           
                            




# In[14]:

#p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.ellipticity,label='Sinuous_Dish-Band-Skirt',color='k',ls='-',lw=2)
#
##p.plot(beamCylinder.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish',color='k',ls='-',lw=2)
###print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
###print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
###p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
###p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
##p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
##p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)
#
#
#
#p.xlabel('f (MHz)',fontsize=20)
#p.ylabel('$\\xi$',fontsize=20)
#p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.title('Far_0.%i-%i-%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))
##p.ylim(1e-4,1e0)
#
#p.gca().tick_params('x',labelsize=16)
#p.gca().tick_params('y',labelsize=16)
#p.gcf().set_size_inches([8,6])
#p.gca().yaxis.grid(which='minor')
#p.grid()
#p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/Far_0.%i-%i-%i_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#p.close()
##p.savefig('../analysis/compareEllipticity.pdf',bbox_inches='tight')


# In[16]:

#p.plot(beamCylinder.fAxis/1e6,beamCylinder.effArea[0]/(pi*7*7),label='HERA Dish',color='k',ls='-',lw=2)
#print beamCylinder.effArea[0][beamCylinder.fAxis==150e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.effArea[0]/(pi*7*7),label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#print beamNoCylinder.effArea[0][beamNoCylinder.fAxis==150e6]
##p.plot(beamNoBackPlane.fAxis/1e6,beamNoBackPlane.effAreaXX/(pi*7*7),label='No BackPlane')
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.effArea[0]/(pi*7*7),label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.effArea[0]/(pi*7*7),label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.effArea[0]/(pi*7*7),label='HERA Dish v2',color='k',ls=':',lw=2)
#
#
#
##print beamBareDipole.effArea[0][beamBareDipole.fAxis==150e6]
#p.xlabel('f (MHz)',fontsize=20)
#p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)
#p.grid()
#p.ylim(.15,.85)
#p.gca().tick_params('x',labelsize=16)
#p.gca().tick_params('y',labelsize=16)
#p.gcf().set_size_inches([8,6])
##p.gca().yaxis.grid(which='minor')
#p.legend(loc='best',fontsize=10,ncol=4)
##p.savefig('../analysis/Percentage_illumination.pdf')

# Y-Direction
#for Growth_Rate in Growth_Rate_List:
#    for Outer_Diameter in Outer_Diameter_List:
#        for Inner_Diameter in Inner_Diameter_List:                    
#            for Band_Resistance in Band_Resistance_List:
#                for Skirt_Diameter in Skirt_Diameter_List:
#                    for Skirt_Height in Skirt_Height_List:
#                        for Port_Number in Port_Number_List:
#                            BeamSinuousDishBandSkirt = BeamSinous_Y('Sinuous_Antenna', Frequency_List,64,Skirt_Diameter,Skirt_Height,Growth_Rate,Outer_Diameter,Inner_Diameter,Band_Resistance,Port_Number, ['XX'],rotateY=True)
#                            p.plot(BeamSinuousDishBandSkirt.fAxis/1e6,BeamSinuousDishBandSkirt.effArea[0]/(pi*7*7),label='Sinuous_Dish-Band-Skirt',color='k',ls='-',lw=2)
#                            p.xlabel('f (MHz)',fontsize=20)
#                            p.ylabel('$A_{eff}/\\pi r^2$',fontsize=20)                            
#                            p.ylim(.15,.85)
#                            p.title('FarArea_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s' %(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height))                              
#                            p.gca().tick_params('y',labelsize=16)
#                            p.gcf().set_size_inches([8,6])
#                            #p.gca().yaxis.grid(which='minor')
#                            p.legend(loc='best',fontsize=10,ncol=1) 
#                            p.grid()
#                            p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/FarArea_Y_0.%i-%i-%i_dish-band_%s-skirt-%s-%s.pdf'%(Growth_Rate,Inner_Diameter, Outer_Diameter,Band_Resistance,Skirt_Diameter,Skirt_Height),bbox_inches='tight')
#                            p.close() 

# In[7]:

#take a cut through a beam
#def hpCut(phi,nPix,data):
#    nSide=hp.npix2nside(len(data))
#    output=n.zeros(nPix)
#    thetaVals=n.arange(nPix/2)/(nPix/2.)*pi/2.
#    thetaVals=n.hstack([n.flipud(thetaVals),thetaVals,]).T
#    phiVals=n.ones(len(thetaVals))
#    phi1=phi+pi
#    phiVals[:nPix/2]=phi1
#    phiVals[nPix/2:]=phi
#    output=data[hp.ang2pix(nSide,thetaVals,phiVals)]
#    return output
#    


# In[8]:

#nth=800
#tha=n.degrees(n.arange(-nth/2,nth/2)*pi/nth)
#p.plot(tha,hpCut(0,nth,10*n.log10(beamCylinder.data[0,5,:])),color='k',ls='--')
#l1=p.plot(tha,hpCut(pi/2.,nth,10*n.log10(beamCylinder.data[0,5,:])),color='k')[0]
#
#
#p.plot(tha,hpCut(0,nth,10*n.log10(beamNoCylinder.data[0,5,:])),color='r',ls='--')
#l2=p.plot(tha,hpCut(pi/2.,nth,10*n.log10(beamNoCylinder.data[0,5,:])),color='r')[0]
#
##p.plot(tha,hpCut(0,nth,10*n.log10(beamNoBackPlane.dataXX[5,:])),color='b',ls='--')
##l3=p.plot(tha,hpCut(pi/2.,nth,10*n.log10(beamNoBackPlane.dataXX[5,:])),color='b')[0]
#
#
#p.plot(tha,hpCut(0,nth,10*n.log10(beamBareDipole.data[0,5,:])),color='orange',ls='--')
#l3=p.plot(tha,hpCut(pi/2.,nth,10*n.log10(beamBareDipole.data[0,5,:])),color='orange')[0]
#
#p.gcf().legend((l1,l2,l3),('Cylinder and Backplane','Backplane Only','Dipole Only'),loc='upper center',ncol=3)
#p.grid()
#p.xlabel('$\\theta$')
#p.xlim(-90,90)
#p.ylim(-50,0)
#p.ylabel('Directivity (dB)')
#p.gcf().set_size_inches([10,7])


# In[9]:

#for m in range(7,8):
#    hp.mollview(10*n.log10(BeamSinuousDishBandSkirt.data[0,m,:]),rot=(0,90,0),min=-40,max=0)
#
#for m in range(7,8):
#    hp.mollview(10*n.log10(beamNoCylinder.data[0,m,:]),rot=(0,90,0),min=-40,max=0)
#
#
## In[10]:
#
#for m in range(7,8):
#    hp.mollview(10*n.log10(beamBareDipole.data[0,m,:]),rot=(0,90,0),min=-40,max=0)


# In[ ]:



