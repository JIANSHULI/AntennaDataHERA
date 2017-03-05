
# coding: utf-8

# In[2]:

#get_ipython().magic(u'matplotlib inline')
import numpy as n
import matplotlib
import matplotlib.pyplot as p
import healpy as hp
import copy
import scipy.optimize as op

pi=n.pi
c=299792458.


# In[3]:

#take a cut through a beam
def hpCut(phi,nPix,data):
    nSide=hp.npix2nside(len(data))
    output=n.zeros(nPix)
    thetaVals=n.arange(nPix/2)/(nPix/2.)*pi/2.
    thetaVals=n.hstack([n.flipud(thetaVals),thetaVals,]).T
    phiVals=n.ones(len(thetaVals))
    phi1=phi+pi
    phiVals[:nPix/2]=phi1
    phiVals[nPix/2:]=phi
    output=hp.get_interp_val(data,thetaVals,phiVals)
    return output


#rotate
def rotateBeam(inputMap,rot=[90,0,0]):
    rotator=hp.Rotator(rot=rot)
    npix=len(inputMap)
    nside=hp.npix2nside(npix)
    theta,phi=hp.pix2ang(nside,range(npix))
    newtheta,newphi=rotator(theta,phi)
    output=hp.get_interp_val(inputMap,newtheta,newphi)
    return output
    


class Beam:
    def __init__(self,dirName,fList,nside,pols=['XX','YY'],rotateY=False,invert=False):
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
                data=n.loadtxt('../data/beams/%s/%s_%s_%s.txt'%(dirName,dirName,fList[m],self.pols[np]),skiprows=2);
                self.data[np,m,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
                if(invert):
                    self.data[np,m,:]=rotateBeam(self.data[np,m,:].flatten(),rot=[0,180,0])
                self.data[np,m,:]/=self.data[np,m,:].flatten().max(); 
                self.data[np,m,theta>90.]=0.
                self.solidAngles[np,m]=self.pixArea*n.sum(self.data[np,m,:])
                self.effArea[np,m]=(c/(self.fAxis[m]))**2./self.solidAngles[np,m]
            if(self.npolsOriginal==1):
                self.data[1,m,:]=rotateBeam(self.data[0,m,:].flatten())
                self.solidAngles[1,m]=self.pixArea*n.sum(self.data[1,m,:])
                self.effArea[1,m]=(c/(self.fAxis[m]))**2./self.solidAngles[1,m]
            if(len(self.pols)>1 and self.pols[0]=='XX' and self.pols[1]=='YY'):
                self.ellipticity[m]=n.sum((self.data[0,m]-self.data[1,m])**2.)/n.sum((self.data[0,m]+self.data[1,m])**2.)

                
                
class FeedData:
    def __init__(self,dirName,dirNameFeedOnly,fList,nside,pols=['XX','YY'],rotateY=False,invertFeedOnly=True,invertDish=False,dDish=14.,dFocus=4.5):
        self.nf=len(fList)
        self.nSide=nside
        self.nPix=hp.nside2npix(nside)
        self.daveFeedEfficiency=n.zeros(self.nf)
        self.daveTaper=n.zeros(self.nf)
        self.davePol=n.zeros(self.nf)
        self.daveBeamEfficiency=n.zeros(self.nf)
        self.dDish=dDish
        self.dFocus=dFocus
        self.thetaEdge=2.*n.arctan(self.dDish/(4.*self.dFocus))
        self.beamFeedAndDish=Beam(dirName,fList,nside,copy.deepcopy(pols),rotateY,invertDish)
        self.beamFeedOnly=Beam(dirNameFeedOnly,fList,nside,pols,rotateY,invertFeedOnly)
        self.fAxis=self.beamFeedOnly.fAxis
        self.beamFWHM=n.zeros(len(self.fAxis))
        self.gain=4.*pi/self.beamFeedAndDish.solidAngles[0,:]
        theta,phi=hp.pix2ang(self.nSide,range(self.nPix))


        self.nSide=nside
        self.nPix=hp.nside2npix(self.nSide)
        self.pols=self.beamFeedAndDish.pols
        for mm in range(self.nf):
            dataFeed=self.beamFeedOnly.data[0,mm,:]
            dataDish=self.beamFeedAndDish.data[0,mm,:]
            #compute feed efficiency by integrating feed beam intercepted by dish
            selectionEdge=theta<self.thetaEdge
            self.daveFeedEfficiency[mm]=n.sum(dataFeed[selectionEdge])/n.sum(dataFeed)
            #compute taper by averaging gain on ring at dish edge
            phiRange=n.radians(n.arange(0,360))
            thetaRange=n.ones(len(phiRange))*self.thetaEdge
            self.daveTaper[mm]=1./n.mean(hp.get_interp_val(dataFeed,thetaRange,phiRange))
            #compute polarization mismatch by integrating phi=0 and phi=90
            thetaRange=n.radians(n.arange(0,180))
            phiRange=n.ones(len(thetaRange))*0.
            xArc=hp.get_interp_val(dataDish,thetaRange,phiRange)
            yArc=hp.get_interp_val(dataDish,thetaRange,phiRange+pi/2.)
            self.davePol[mm]=n.sum((n.sqrt(xArc)-n.sqrt(yArc))**2.)/n.sum((n.sqrt(xArc)+n.sqrt(yArc))**2.)
            #compute beam efficiency. First determine FWHM by taking average in 360 degree ring. 
            stddevs=[]
            nth=1000
            tha=n.arange(-nth/2,nth/2)*pi/nth
            for nn in range(0,180,10):
                cut=hpCut(n.radians(nn),nth,dataDish)
                try:
                    stddevs.append(op.curve_fit(lambda x,p: n.exp(-x**2./(2.*p**2.)),tha,cut,p0=[n.radians(15.)])[0])
                except Exception as e:
                    print e
                    print 'error'
                    p.plot(tha,cut)
                    p.yscale('log')
                    p.show()
            stddeves=n.array(stddevs)
            stddev=n.max(stddevs)
            fwhm=2.*n.sqrt(2.*n.log(2.))*stddev
            print n.degrees(fwhm)
            self.beamFWHM[mm]=fwhm
            selectionFWHM=theta<fwhm/2.
            self.daveBeamEfficiency[mm]=n.sum(dataDish[selectionFWHM])/n.sum(dataDish)


# In[5]:

#beamCylinder=Beam('beamCylinder',['100'],64)
fstrList=['050','060','070','080','090','100','110','120','130','140','150']
fstrListHigh=['100','110','120','130','140','150','160','170','180','190','200']


'''
FeedHighBand=FeedData('beamHighBand','beamHighBandFeedOnly',fstrListHigh,64,['XX'],rotateY=True,invertFeedOnly=False)
FeedPickets=FeedData('beamPaneledCylinder','beamPicketCylinderFeedOnly_v2',fstrList,64,['XX'],rotateY=True)
FeedCylinder=FeedData('beamCylinderBackPlane_v2','beamCylinderFeedOnly_v2',fstrList,64,['XX'],rotateY=True)
FeedBackPlane=FeedData('beamBackPlane_v2','beamBackPlaneFeedOnly_v2',fstrList,64,['XX'],rotateY=True)
FeedCylinderAdaptive=FeedData('beamCylinderBackPlaneAdaptive_v2','beamCylinderFeedOnly_v2',fstrList,64,['XX'],rotateY=True)
FeedCylinderShortNoGap=FeedData('beamCylinderBackPlaneShortNoGap_v2','beamCylinderBackPlaneShortNoGapFeedOnly_v2',fstrList,64,['XX'],rotateY=True)
FeedCylinderShortGap=FeedData('beamCylinderBackPlaneShortGap_v2','beamCylinderBackPlaneShortGapFeedOnly',fstrList,64,['XX'],rotateY=True)
FeedCylinderNoGap=FeedData('beamCylinderBackPlaneNoGap_v2','beamCylinderBackPlaneNoGapFeedOnly_v2',fstrList,64,['XX'],rotateY=True)
'''
#FeedCTP=FeedData('beamCTPFeedOverDish','beamCTPFeedOnly',fstrList,64,['XX'],rotateY=True,invertFeedOnly=False,invertDish=True)
#FeedCTPShort=FeedData('beamCTPShortFeedOverDish','beamCTPShortFeedOnly',fstrList,64,['XX'],rotateY=True,invertFeedOnly=False,invertDish=True)
FeedCTPShortNoGap=FeedData('beamCTPShortNoGapFeedOverDish','beamCTPShortNoGapFeedOnly',fstrList,64,['XX'],rotateY=True,invertFeedOnly=False,invertDish=True)


# In[8]:

#print FeedPickets.daveBeamEfficiency.shape
p.plot(FeedCTP.fAxis/1e6,10.*n.log10(FeedCTP.gain),label='CTP',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTPShort.fAxis/1e6,10.*n.log10(FeedCTPShort.gain),label='CTP Short',color='r',ls='-',lw=2)
p.plot(FeedCTPShortNoGap.fAxis/1e6,10.*n.log10(FeedCTPShortNoGap.gain),label='CTP Short No Gap',color='grey',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Directivity (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Directivity_CTP_Only.pdf',bbox_inches='tight')


# In[18]:

#print FeedPickets.daveBeamEfficiency.shape
p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.gain),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTP.fAxis/1e6,10.*n.log10(FeedCTP.gain),label='CTP',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Directivity (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Directivity.pdf',bbox_inches='tight')


# In[18]:

#print FeedPickets.daveBeamEfficiency.shape
p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.gain),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Directivity (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Directivity.pdf',bbox_inches='tight')


# In[26]:

print FeedPickets.daveBeamEfficiency.shape
l0=p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.gain),label='Cylinder Feed',color='k',lw=2,ls='-')[0]
p.plot(FeedCylinderNoGap.fAxis/1e6,10.*n.log10(FeedCylinderNoGap.gain),color='k',ls='--',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
l1=p.plot(FeedCylinderShortGap.fAxis/1e6,10.*n.log10(FeedCylinderShortGap.gain),color='r',lw=2,ls='-')[0]
p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Short Cylinder Feed\n No Gap',color='r',ls='--',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
#p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Cylinder Feed',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Directivity (dB)',fontsize=20)
p.gca().legend((l0,l1),('Tall Feed','Short Feed'),loc='best',fontsize=10,ncol=2)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/DirectivityCompareGap.pdf',bbox_inches='tight')


# In[5]:

print FeedPickets.daveBeamEfficiency.shape
p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.gain),label='Cylinder Feed',color='k',ls='-',lw=2)
p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)
p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
#p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Cylinder Feed',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Directivity (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=2)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/Directivity.pdf',bbox_inches='tight')


# In[13]:

print FeedCTP.daveBeamEfficiency



# In[7]:


p.plot(FeedCTP.fAxis/1e6,(FeedCTP.daveBeamEfficiency),label='CTP',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.daveBeamEfficiency),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveBeamEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTPShort.fAxis/1e6,(FeedCTPShort.daveBeamEfficiency),label='CTP Short',color='r',ls='-',lw=2)

p.plot(FeedCTPShortNoGap.fAxis/1e6,(FeedCTPShortNoGap.daveBeamEfficiency),label='CTP Short\n No Gap',color='grey',ls='-',lw=2)


#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Main-lobe efficiency: $\\epsilon_b$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Efficiency.pdf',bbox_inches='tight')


# In[19]:


p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.daveBeamEfficiency),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.daveBeamEfficiency),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveBeamEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTP.fAxis/1e6,(FeedCTP.daveBeamEfficiency),label='CTP',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Main-lobe efficiency: $\\epsilon_b$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Efficiency.pdf',bbox_inches='tight')


# In[9]:

print FeedPickets.daveBeamEfficiency.shape

p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.daveBeamEfficiency),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.daveBeamEfficiency),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveBeamEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.daveBeamEfficiency),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Main-lobe efficiency: $\\epsilon_b$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Efficiency.pdf',bbox_inches='tight')


# In[28]:

print FeedPickets.daveBeamEfficiency.shape
l0=p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.daveBeamEfficiency),label='Cylinder Feed',color='k',lw=2,ls='-')[0]
p.plot(FeedCylinderNoGap.fAxis/1e6,(FeedCylinderNoGap.daveBeamEfficiency),color='k',ls='--',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
l1=p.plot(FeedCylinderShortGap.fAxis/1e6,(FeedCylinderShortGap.daveBeamEfficiency),color='r',lw=2,ls='-')[0]
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.daveBeamEfficiency),label='Short Cylinder Feed\n No Gap',color='r',ls='--',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
#p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Cylinder Feed',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Main Lobe Efficiency, $\\epsilon_b$',fontsize=20)
p.gca().legend((l0,l1),('Tall Feed','Short Feed'),loc='best',fontsize=10,ncol=2)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/EfficiencyCompareGap.pdf',bbox_inches='tight')


# In[6]:

print FeedPickets.daveBeamEfficiency.shape

p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.daveBeamEfficiency),label='Cylinder Feed',color='k',ls='-',lw=2)
p.plot(FeedPickets.fAxis/1e6,(FeedPickets.daveBeamEfficiency),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.daveBeamEfficiency),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)

p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveBeamEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)


#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Main-lobe efficiency: $\\epsilon_b$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=2)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/Efficiency.pdf',bbox_inches='tight')


# In[9]:

#print FeedPickets.daveFeedEfficiency.shape
p.plot(FeedCTP.fAxis/1e6,FeedCTP.daveFeedEfficiency,label='CTP',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,FeedPickets.daveFeedEfficiency,label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveFeedEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTPShort.fAxis/1e6,(FeedCTPShort.daveFeedEfficiency),label='CTP Short',color='r',ls='-',lw=2)
p.plot(FeedCTPShortNoGap.fAxis/1e6,(FeedCTPShortNoGap.daveFeedEfficiency),label='CTP Short No Gap',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Efficiency: $\\epsilon_f$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/feedEfficiency.pdf',bbox_inches='tight')


# In[22]:

#print FeedPickets.daveFeedEfficiency.shape
p.plot(FeedCylinder.fAxis/1e6,FeedCylinder.daveFeedEfficiency,label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,FeedPickets.daveFeedEfficiency,label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveFeedEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTP.fAxis/1e6,(FeedCTP.daveFeedEfficiency),label='CTP',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Efficiency: $\\epsilon_f$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/feedEfficiency.pdf',bbox_inches='tight')


# In[10]:

#print FeedPickets.daveFeedEfficiency.shape
p.plot(FeedCylinder.fAxis/1e6,FeedCylinder.daveFeedEfficiency,label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,FeedPickets.daveFeedEfficiency,label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveFeedEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.daveFeedEfficiency),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Efficiency: $\\epsilon_f$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/feedEfficiency.pdf',bbox_inches='tight')


# In[27]:

print FeedPickets.daveBeamEfficiency.shape
l0=p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.daveFeedEfficiency),label='Cylinder Feed',color='k',lw=2,ls='-')[0]
p.plot(FeedCylinderNoGap.fAxis/1e6,(FeedCylinderNoGap.daveFeedEfficiency),color='k',ls='--',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
l1=p.plot(FeedCylinderShortGap.fAxis/1e6,(FeedCylinderShortGap.daveFeedEfficiency),color='r',lw=2,ls='-')[0]
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.daveFeedEfficiency),label='Short Cylinder Feed\n No Gap',color='r',ls='--',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
#p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Cylinder Feed',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Efficiency,$\\epsilon_f$',fontsize=20)
p.gca().legend((l0,l1),('Tall Feed','Short Feed'),loc='best',fontsize=10,ncol=2)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/feedEfficiencyCompareGap.pdf',bbox_inches='tight')


# In[7]:

print FeedPickets.daveFeedEfficiency.shape
p.plot(FeedCylinder.fAxis/1e6,FeedCylinder.daveFeedEfficiency,label='Cylinder Feed',color='k',ls='-',lw=2)
p.plot(FeedPickets.fAxis/1e6,FeedPickets.daveFeedEfficiency,label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.daveFeedEfficiency),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)
p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.daveFeedEfficiency),label='High Band (100-200MHz)',color='k',ls='--',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Efficiency: $\\epsilon_f$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=2)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/feedEfficiency.pdf',bbox_inches='tight')


# In[10]:

p.plot(FeedCTP.fAxis/1e6,10.*n.log10(FeedCTP.daveTaper),label='CTP',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.daveTaper),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.daveTaper),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTPShort.fAxis/1e6,10.*n.log10(FeedCTPShort.daveTaper),label='CTP Short',color='r',ls='-',lw=2)
p.plot(FeedCTPShortNoGap.fAxis/1e6,10.*n.log10(FeedCTPShortNoGap.daveTaper),label='CTP Short No Gap',color='grey',ls='-',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Taper: $T_f$ (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.ylim(0,22)
#p.savefig('../analysis/feedTaper.pdf',bbox_inches='tight')




# In[23]:

p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.daveTaper),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.daveTaper),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.daveTaper),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTP.fAxis/1e6,10.*n.log10(FeedCTP.daveTaper),label='CTP',color='r',ls='-',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Taper: $T_f$ (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.ylim(0,22)
#p.savefig('../analysis/feedTaper.pdf',bbox_inches='tight')




# In[11]:

p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.daveTaper),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.daveTaper),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.daveTaper),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.daveTaper),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Taper: $T_f$ (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.ylim(0,22)
#p.savefig('../analysis/feedTaper.pdf',bbox_inches='tight')



# In[30]:

print FeedPickets.daveBeamEfficiency.shape
l0=p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.daveTaper),label='Cylinder Feed',color='k',lw=2,ls='-')[0]
p.plot(FeedCylinderNoGap.fAxis/1e6,10.*n.log10(FeedCylinderNoGap.daveTaper),color='k',ls='--',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
l1=p.plot(FeedCylinderShortGap.fAxis/1e6,10.*n.log10(FeedCylinderShortGap.daveTaper),color='r',lw=2,ls='-')[0]
p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.daveTaper),label='Short Cylinder Feed\n No Gap',color='r',ls='--',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
#p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Cylinder Feed',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Taper,$T_f$ (dB)',fontsize=20)
p.gca().legend((l0,l1),('Tall Feed','Short Feed'),loc='best',fontsize=10,ncol=2)
#p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/feedTaperCompareGap.pdf',bbox_inches='tight')


# In[8]:

p.plot(FeedCylinder.fAxis/1e6,10.*n.log10(FeedCylinder.daveTaper),label='Cylinder Feed',color='k',ls='-',lw=2)
p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.daveTaper),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.daveTaper),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)

p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.daveTaper),label='High Band (100-200MHz)',color='k',ls='--',lw=2)


p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Feed Taper: $T_f$ (dB)',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=2)
p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.ylim(0,22)
p.savefig('../analysis/feedTaper.pdf',bbox_inches='tight')


# In[11]:


p.plot(FeedCTP.fAxis/1e6,(FeedCTP.davePol),label='CTP',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.davePol),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.davePol),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTPShort.fAxis/1e6,(FeedCTPShort.davePol),label='CTP Short',color='r',ls='-',lw=2)
p.plot(FeedCTPShortNoGap.fAxis/1e6,(FeedCTPShortNoGap.davePol),label='CTP Short No Gap',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Polarization Mismatch: $\\xi_s$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/PolarizationSimple.pdf',bbox_inches='tight')


# In[24]:


p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.davePol),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.davePol),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.davePol),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCTP.fAxis/1e6,(FeedCTP.davePol),label='CTP',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Polarization Mismatch: $\\xi_s$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/PolarizationSimple.pdf',bbox_inches='tight')


# In[12]:


p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.davePol),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.davePol),label='Picket Feed',color='grey',ls='-',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.davePol),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.davePol),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Polarization Mismatch: $\\xi_s$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/PolarizationSimple.pdf',bbox_inches='tight')


# In[31]:

print FeedPickets.daveBeamEfficiency.shape
l0=p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.davePol),label='Cylinder Feed',color='k',lw=2,ls='-')[0]
p.plot(FeedCylinderNoGap.fAxis/1e6,(FeedCylinderNoGap.davePol),color='k',ls='--',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
l1=p.plot(FeedCylinderShortGap.fAxis/1e6,(FeedCylinderShortGap.davePol),color='r',lw=2,ls='-')[0]
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.davePol),label='Short Cylinder Feed\n No Gap',color='r',ls='--',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
#p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Cylinder Feed',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Polarization Mismatch: $\\xi_s$',fontsize=20)
p.gca().legend((l0,l1),('Tall Feed','Short Feed'),loc='best',fontsize=10,ncol=2)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/PolarizationSimpleCompareGap.pdf',bbox_inches='tight')


# In[9]:


p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.davePol),label='Cylinder Feed',color='k',ls='-',lw=2)
p.plot(FeedPickets.fAxis/1e6,(FeedPickets.davePol),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.davePol),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)

p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.davePol),label='High Band (100-200MHz)',color='k',ls='--',lw=2)


#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Polarization Mismatch: $\\xi_s$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=2)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/PolarizationSimple.pdf',bbox_inches='tight')


# In[10]:

print FeedPickets.daveBeamEfficiency.shape

p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.davePol),label='Cylinder Feed',color='k',ls='-',lw=2)
p.plot(FeedCylinderAdaptive.fAxis/1e6,(FeedCylinderAdaptive.davePol),label='Cylinder Feed Adaptive',color='grey',ls='--',lw=2)


#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Polarization Mismatch: $\\xi_s$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/PolarizationSimple.pdf',bbox_inches='tight')


# In[12]:


p.plot(FeedCTP.fAxis/1e6,(FeedCTP.beamFeedAndDish.ellipticity),label='CTP',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.beamFeedAndDish.ellipticity),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCTPShort.fAxis/1e6,(FeedCTPShort.beamFeedAndDish.ellipticity),label='CTP Short',color='r',ls='-',lw=2)
p.plot(FeedCTPShortNoGap.fAxis/1e6,(FeedCTPShortNoGap.beamFeedAndDish.ellipticity),label='CTP Short No Gap',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Beam Ellipticity: $\\xi$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Polarization.pdf',bbox_inches='tight')


# In[25]:

print FeedPickets.daveBeamEfficiency.shape

p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.beamFeedAndDish.ellipticity),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.beamFeedAndDish.ellipticity),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCTP.fAxis/1e6,(FeedCTP.beamFeedAndDish.ellipticity),label='CTP',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Beam Ellipticity: $\\xi$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Polarization.pdf',bbox_inches='tight')


# In[14]:

print FeedPickets.daveBeamEfficiency.shape

p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.beamFeedAndDish.ellipticity),label='Cylinder Feed',color='k',ls='-',lw=2)
#p.plot(FeedPickets.fAxis/1e6,(FeedPickets.beamFeedAndDish.ellipticity),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.beamFeedAndDish.ellipticity),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Beam Ellipticity: $\\xi$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=4)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
#p.savefig('../analysis/Polarization.pdf',bbox_inches='tight')


# In[32]:

print FeedPickets.daveBeamEfficiency.shape
l0=p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.beamFeedAndDish.ellipticity),label='Cylinder Feed',color='k',lw=2,ls='-')[0]
p.plot(FeedCylinderNoGap.fAxis/1e6,(FeedCylinderNoGap.beamFeedAndDish.ellipticity),color='k',ls='--',lw=2)
#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedPickets.gain),label='Picket Feed',color='grey',ls='-',lw=2)
l1=p.plot(FeedCylinderShortGap.fAxis/1e6,(FeedCylinderShortGap.beamFeedAndDish.ellipticity),color='r',lw=2,ls='-')[0]
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.beamFeedAndDish.ellipticity),label='Short Cylinder Feed\n No Gap',color='r',ls='--',lw=2)
#p.plot(FeedHighBand.fAxis/1e6-50,10.*n.log10(FeedHighBand.gain),label='High Band (100-200MHz)',color='k',ls='--',lw=2)
#p.plot(FeedCylinderShortNoGap.fAxis/1e6,10.*n.log10(FeedCylinderShortNoGap.gain),label='Cylinder Feed',color='r',ls='-',lw=2)


#p.plot(FeedPickets.fAxis/1e6,10.*n.log10(FeedBackPlane.gain),label='Feed with Backplane only',color='grey',ls='-',lw=2)



#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Beam Ellipticity: $\\xi$',fontsize=20)
p.gca().legend((l0,l1),('Tall Feed','Short Feed'),loc='best',fontsize=10,ncol=2)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/PolarizationCompareGap.pdf',bbox_inches='tight')


# In[10]:

print FeedPickets.daveBeamEfficiency.shape

p.plot(FeedCylinder.fAxis/1e6,(FeedCylinder.beamFeedAndDish.ellipticity),label='Cylinder Feed',color='k',ls='-',lw=2)
p.plot(FeedPickets.fAxis/1e6,(FeedPickets.beamFeedAndDish.ellipticity),label='Picket Feed',color='grey',ls='-',lw=2)
p.plot(FeedCylinderShortNoGap.fAxis/1e6,(FeedCylinderShortNoGap.beamFeedAndDish.ellipticity),label='Short Cylinder Feed\n No Gap',color='r',ls='-',lw=2)
p.plot(FeedHighBand.fAxis/1e6-50,(FeedHighBand.beamFeedAndDish.ellipticity),label='High Band (100-200MHz)',color='k',ls='--',lw=2)


#print beamCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#print beamNoCylinder.ellipticity[beamNoCylinder.fAxis==100e6]
#p.plot(beamNoCylinder.fAxis/1e6,beamNoCylinder.ellipticity,label='HERA Dish\nNo Cylinder',color='k',ls='--',lw=2)
#p.plot(beamBareDipole.fAxis/1e6,beamBareDipole.ellipticity,label='HERA Dish\nDipole Only',color='grey',lw=2,ls='-')
#p.plot(beamPaneledCylinder.fAxis/1e6,beamPaneledCylinder.ellipticity,label='HERA Dish\nPaneled Cylinder',color='grey',lw=2,ls='--')
#p.plot(beamCylinderV2.fAxis/1e6,beamCylinder.ellipticity,label='HERA Dish v2',color='k',ls=':',lw=2)



p.xlabel('f (MHz)',fontsize=20)
p.ylabel('Beam Ellipticity: $\\xi$',fontsize=20)
p.legend(loc='best',fontsize=10,ncol=2)
p.yscale('log')
#p.ylim(1e-4,1e0)

p.gca().tick_params('x',labelsize=16)
p.gca().tick_params('y',labelsize=16)
p.gcf().set_size_inches([8,6])
p.gca().yaxis.grid(which='minor')
p.grid()
p.savefig('../analysis/Polarization.pdf',bbox_inches='tight')


# In[ ]:



