#************************************************************
#library for reading cst and VNA data files
#************************************************************
import numpy as n
import healpy as hp
import numpy.fft as fft
import scipy.signal as signal

class MetaData():
    def __init__(self,device='',date='',dtype='',comment='',datarange=[]):
        self.device=device
        self.date=date
        self.dtype=dtype
        self.comment=comment
        self.datarange=datarange

#Read HP VNA data used in Greenbank measurements
def readVNAHP(fileName,comment=''):
    dataFile=open(fileName)
    dataLines=dataFile.readlines()
    FLOW=1e-9*float(dataLines[6].split()[1]);FHIGH=1e-9*float(dataLines[6].split()[2]);
    NDATA=int(dataLines[6].split()[3])
    data=n.loadtxt(fileName,skiprows=9,delimiter=',')
    device=dataLines[1][:-2]
    dtype=dataLines[4].split()[1]
    meta=MetaData(device=device,dtype=['FREQ',dtype],datarange=[FLOW,FHIGH,NDATA],comment=comment)
    fAxis=n.arange(NDATA)*(FHIGH-FLOW)/NDATA+FLOW
    print fAxis
    return fAxis,data[:,0]+1j*data[:,1],meta

#take ratio of fft of two inputs with padding
def fftRatio(convolved,kernel):
    nf=len(convolved)
    convolved_pad=n.pad(convolved,(nf/2,nf/2),mode='constant')
    kernel_pad=n.pad(kernel,(nf/2,nf/2),mode='constant')
    return fft.fftshift(fft.fft(convolved_pad)/fft.fft(kernel_pad))
    

#Read CST time trace file
def readCSTTimeTrace(fileName,comment=''):
    dataFile=open(fileName)
    dataLines=dataFile.readlines()
    header=dataLines[:2]
    if('ns' in header[0]):
        tFactor=1.
    if('ms' in header[0]):
        tFactor=1e6
    if('micro' in header[0]):
        tFactor=1e3
    if('sec' in header[0]):
        tFactor=1e9
        
    inputTrace=[]
    outputTrace1=[]
    outputTrace2=[]
    lNum=0
    while lNum <len(dataLines):
        if('o1' in dataLines[lNum]):
            thisTrace=outputTrace1
            lNum+=2
        elif('o2' in dataLines[lNum]):
            thisTrace=outputTrace2
            lNum+=2
        elif('i1' in dataLines[lNum]):
            thisTrace=inputTrace
            lNum+=2
            dtype='Terminal Excitation'
        elif('Plane wave' in dataLines[lNum]):
            thisTrace=inputTrace
            lNum+=2
            dtype='PlaneWave Excitation'
        else:
            entry=dataLines[lNum].split()
            if(len(entry)==2):
                thisTrace.append([float(entry[0]),float(entry[1])])
            lNum+=1
    inputTrace=n.array(inputTrace)
    outputTrace1=n.array(outputTrace1)
    outputTrace2=n.array(outputTrace2)
    inputTrace[:,0]*=tFactor
    outputTrace1[:,1]*=tFactor
    if(len(outputTrace2)>0):
        outputTrace2[:,1]*=tFactor
    meta=MetaData(device='CST',dtype=['TIME',dtype],datarange=[inputTrace[:,0].min(),inputTrace[:,0].max(),len(inputTrace[:,0])],comment=comment)
    return [inputTrace,outputTrace1,outputTrace2],meta
    
def readCSTS11(fileName,comment='',degrees=True):
    dB=False
    header=open(fileName+'_abs.txt').readlines()[:2]
    if('MHz' in header[0]):
        fFactor=1e-3
    elif('GHz' in header[0]):
        fFactor=1e0
    elif('kHz' in header[0]):
        fFactor=1e-6
    elif('Hz' in header[0]):
        fFactor=1e-9
    if('dB' in header[0]):
        dB=True
    amp=n.loadtxt(fileName+'_abs.txt',skiprows=2)
    fAxis=amp[:,0]
    amp=amp[:,1]
    if(dB):
        amp=10.**(amp/20.)

    pha=n.loadtxt(fileName+'_pha.txt',skiprows=2)[:,1]
    if(degrees):
        pha*=n.pi/180.
    data=amp*n.exp(1j*pha)
    meta=MetaData(device='CST',dtype=['FREQ','S11'],datarange=[fAxis.min(),fAxis.max(),len(fAxis)],comment=comment)
    return fFactor*fAxis,data,meta
    

FILETYPES=['CST_TimeTrace','CST_S11','VNAHP_S11']
class GainData():
    def __init__(self,fileName,fileType,fMin=None,fMax=None,windowFunction=None,comment='',filterNegative=False):
        assert fileType in FILETYPES
        if(windowFunction is None):
            windowFunction = 'blackman-harris'
        self.windowFunction=windowFunction
        if (fileType=='CST_TimeTrace'):
            [inputTrace,outputTrace,_],self.metaData=readCSTTimeTrace(fileName,comment=comment)
            self.fAxis=fft.fftshift(fft.fftfreq(len(inputTrace)*2,inputTrace[1,0]-inputTrace[0,0]))
            self.gainFrequency=fftRatio(outputTrace[:,1],inputTrace[:,1])
            
        elif(fileType=='CST_S11'):
            self.fAxis,self.gainFrequency,self.metaData=readCSTS11(fileName,comment=comment)
        elif(fileType=='VNAHP_S11'):
            self.fAxis,self.gainFrequency,self.metaData=readVNAHP(fileName,comment=comment)
        if(fMin is None):
            fMin=self.fAxis.min()            
        if(fMax is None):
            fMax=self.fAxis.max()
            
            
        selection=n.logical_and(self.fAxis>=fMin,self.fAxis<=fMax)
        self.fAxis=self.fAxis[selection]
        self.gainFrequency=self.gainFrequency[selection]
        if(windowFunction== 'blackman-harris'):
            wF=signal.blackmanharris(len(self.fAxis))
            wF/=n.sqrt(n.mean(wF**2.))
        self.tAxis=fft.fftshift(fft.fftfreq(len(self.fAxis),self.fAxis[1]-self.fAxis[0]))
        if(filterNegative):
            gainDelay=fft.fftshift(fft.ifft(fft.fftshift(self.gainFrequency)))
            gainDelay[self.tAxis<0.]=0.
            self.gainFrequency=fft.fftshift(fft.fft(fft.fftshift(gainDelay)))
        self.gainDelay=fft.fftshift(fft.ifft(fft.fftshift(self.gainFrequency*wF)))
    
