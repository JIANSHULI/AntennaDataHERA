import numpy as n
import matplotlib.pyplot as p
import gainData as gainData

#fileNameTimeTraceCST='../cst/TallCylinderGapOverDish/TallCylinderGapOverDish_TerminalExcitation_timetrace.txt'
#fileNameS11VNA='../reflectometry/RichBradley_GreenBank/TallCylinderGapOverDish_S11_Greenbank_RichBradley.d1'
#fileNameS11CST='../cst/TallCylinderGapOverDish/TallCylinderGapOverDish_S11'


#fileNameTimeTraceCST='../Results/Sinuous_Antenna/TimeDomain_0.65-30-325-MP9.txt'
#fileNameS11CST='../Results/Sinuous_Antenna/S11_0.65-30-325-MP9'
fileNameTimeTraceCST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/TimeDomain_0.6-30-325-MP8.txt'
fileNameS11CST='/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/S11_0.6-30-325-MP8'
#fileNameS11VNA='../reflectometry/RichBradley_GreenBank/TallCylinderGapOverDish_S11_Greenbank_RichBradley.d1'

FLOW=0.05
FHIGH=0.25


gainData_timeTrace=gainData.GainData(fileNameTimeTraceCST,
									 fileType='CST_TimeTrace',
									 fMin=FLOW,fMax=FHIGH,
									 comment='s11 derived from cst time domain data',filterNegative=True)
gainData_cst=gainData.GainData(fileNameS11CST,
							   fileType='CST_S11',
							   fMin=FLOW,fMax=FHIGH,
							   comment='s11 obtained directly from cst',filterNegative=True)
#gainData_vna=gainData.GainData(fileNameS11VNA,
#							   fileType='VNAHP_S11',
#							   fMin=FLOW,fMax=FHIGH,
#							   comment='s11 obtained from richs vna measurement',filterNegative=True)

print gainData_cst.gainFrequency.shape

#first make original plot comparing s11 of time trace and s11 of vna

#p.plot(gainData_vna.tAxis,10.*n.log10(n.abs(gainData_vna.gainDelay)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
p.plot(gainData_timeTrace.tAxis,10.*n.log10(n.abs(gainData_timeTrace.gainDelay)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
p.plot(gainData_cst.tAxis,10.*n.log10(n.abs(gainData_cst.gainDelay)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
p.xlim(-30,400)
p.ylim(-70,0)
p.ylabel('|$\widetilde{S}_{11}$|(dB)')
p.xlabel('delay (ns)')
p.legend(loc='best')
#p.show()
p.grid()
#p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Delay_NF.pdf',bbox_inches='tight')
p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Delay_NF.pdf',bbox_inches='tight')
p.close()

#p.plot(gainData_vna.fAxis,10.*n.log10(n.abs(gainData_vna.gainFrequency)),color='grey',ls='-',marker='o',label='VNA Measurement',markersize=4,markeredgecolor='none')
p.plot(gainData_timeTrace.fAxis,10.*n.log10(n.abs(gainData_timeTrace.gainFrequency)),color='k',ls='-',marker='o',label='CST timetrace',markersize=4,markeredgecolor='none')
p.plot(gainData_cst.fAxis,10.*n.log10(n.abs(gainData_cst.gainFrequency)),color='k',ls='--',marker='o',label='CST $S_{11}$',markersize=4,markeredgecolor='none')
p.xlim(.045,.255)
p.ylim(-25,0)
p.ylabel('|S$_{11}$|(dB)')
p.xlabel('f (GHz)')
p.legend(loc='best')
#p.show()
p.grid()
#p.savefig('../plots/s11_CST_vs_ReflectometryRich_TallCylinderGapFeedOnly_Frequency_NF.pdf',bbox_inches='tight')
p.savefig('/Users/JianshuLi/Documents/Miracle/Research/Cosmology/21cm Cosmology/Results/Sinuous_Antenna/Plots/S11_CST_Frequency_NF.pdf',bbox_inches='tight')
p.close()

