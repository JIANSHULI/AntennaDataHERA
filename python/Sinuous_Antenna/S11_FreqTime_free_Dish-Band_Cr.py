#!/usr/bin/python

import numpy as n
import matplotlib
from matplotlib import pylab as plt
import matplotlib.pyplot as p
import gainData as gainData

#fileNameTimeTraceCST='../cst/TallCylinderGapOverDish/TallCylinderGapOverDish_TerminalExcitation_timetrace.txt'
#fileNameS11VNA='../reflectometry/RichBradley_GreenBank/TallCylinderGapOverDish_S11_Greenbank_RichBradley.d1'
#fileNameS11CST='../cst/TallCylinderGapOverDish/TallCylinderGapOverDish_S11'

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



Growth_Rate=50
Outer_Diameter=275
Inner_Diameter=30
MP=8

PW = S11_Power = 1

Growth_Rate_List = [50]
Outer_Diameter_List = [275]
Inner_Diameter_List = []

for Growth_Rate in Growth_Rate_List:
	for Outer_Diameter in Outer_Diameter_List:

		#fileNameTimeTraceCST='../Results/Sinuous_Antenna/TimeDomain_0.65-30-325-MP9.txt'
		#fileNameS11CST='../Results/Sinuous_Antenna/S11_0.65-30-325-MP9'
		fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.%i-%i-%i_dish-band.txt' %(Growth_Rate, Inner_Diameter,Outer_Diameter)
		fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.%i-%i-%i_dish-band' %(Growth_Rate, Inner_Diameter,Outer_Diameter)
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
		p.title('S11_CST_Delay_0.%i-%i-%i_PW%i_dish-band' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power))
		#p.show()
		p.grid()
		#p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay.pdf',bbox_inches='tight')
		p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_0.%i-%i-%i_PW%i_Cr_dish-band.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power),bbox_inches='tight')
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
		p.title('S11_CST_Frequency_0.%i-%i-%i_PW%i_dish-band' %(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power)) 
		#p.show()
		p.grid()
		#p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency.pdf',bbox_inches='tight')
		p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_0.%i-%i-%i_PW%i_Cr_dish-band.pdf'%(Growth_Rate, Inner_Diameter,Outer_Diameter,S11_Power),bbox_inches='tight')
		p.close()