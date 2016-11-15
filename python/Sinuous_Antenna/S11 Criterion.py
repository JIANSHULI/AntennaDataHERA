import numpy as np 
#%matplotlib inline
import matplotlib
from matplotlib import pylab as plt


def calc_Tsky(nu, fq0=.180, T0=180., alpha=-2.5):
	'''Tsky = T0 * (nu/fq0)^alpha. nu, fq0 in GHz, T0 in K.''' 
	return T0 * (nu/fq0)**alpha

PW = 1

nu = np.linspace(.050,.250, 1024) 
Trx1 = 100 # K
Trx2 = 75 # K
Tsky = calc_Tsky(nu)
plt.semilogy(nu, Tsky, 'b', label='$T_{\\rm sky}$')
plt.semilogy(nu, Trx1 * np.ones_like(nu), 'r--', label='$T_{\\rm rx}=100$ K')
plt.semilogy(nu, Trx2 * np.ones_like(nu), 'k--', label='$T_{\\rm rx}=75$ K')
plt.xlim(.05,.25)
plt.ylim(50, 5000)
plt.xlabel('Frequency [GHz]')
plt.ylabel('Temperature [K]')
plt.legend(loc='best')
plt.show()



e_dish = 0.74
eT_rx = np.around(e_dish * (500 - calc_Tsky(.150)))

#print 'T_rx,eff =', eT_rx, 'K'

S11_150_1 = np.around(5*PW*np.log10(1 - Trx1/eT_rx))
S11_150_2 = np.around(5*PW*np.log10(1 - Trx2/eT_rx))
#print 'S_11 =', S11_150_1, 'dB and', S11_150_2,'dB, respectively.'

eT_rx = np.around(e_dish * (500 - Tsky))
plt.semilogy(nu, Tsky, 'b', label='$T_{\\rm sky}$')
plt.semilogy(nu, eT_rx, 'r', label='$T_{\\rm rx}$ target') 
plt.semilogy(nu, Trx1 * np.ones_like(nu), 'r--', label='$T_{\\rm rx}=100$ K') 
plt.semilogy(nu, Trx2 * np.ones_like(nu), 'k--', label='$T_{\\rm rx}=75$ K') 
plt.xlim(.05,.25); plt.ylim(50, 5000)
plt.xlabel('Frequency [GHz]')
plt.ylabel('Temperature [K]')
plt.legend(loc='best')
plt.show()

S11_1 = (1 - Trx1/eT_rx).clip(1e-10,1)
S11_2 = (1 - Trx2/eT_rx).clip(1e-10,1)
plt.plot(nu, 5*PW*np.log10(S11_1), 'r', label='$T_{\\rm rx}=100$ K') 
plt.plot(nu, 5*PW*np.log10(S11_2), 'b', label='$T_{\\rm rx}=75$ K') 
plt.xlim(.05,.25); plt.ylim(-10*PW, 0)
plt.xlabel('Frequency [GHz]')
plt.ylabel('$S_{11}$ [dB]')
plt.legend(loc='best')
plt.grid()
plt.show()

f_obs = 2
e_dish = 0.74
Tsys = np.sqrt(f_obs) * Tsky
Trx_target = e_dish * (Tsys - Tsky)
plt.semilogy(nu, Tsky, 'b', label='$T_{\\rm sky}$')
plt.semilogy(nu, Tsys, 'g', label='$T_{\\rm sys}$ target') 
plt.semilogy(nu, Trx_target, 'r', label='$T_{\\rm rx}$ target') 
plt.semilogy(nu, Trx1 * np.ones_like(nu), 'r--', label='$T_{\\rm rx}=100$ K') 
plt.semilogy(nu, Trx2 * np.ones_like(nu), 'k--', label='$T_{\\rm rx}=75$ K') 
plt.xlim(.05,.25); plt.ylim(5, 5000)
plt.xlabel('Frequency [GHz]')
plt.ylabel('Temperature [K]')
plt.legend(loc='best')
plt.grid()
plt.show()

S11_1 = (1 - Trx1 / Trx_target).clip(0,1)
S11_2 = (1 - Trx2 / Trx_target).clip(0,1)
plt.plot(nu, 5*PW*np.log10(S11_1), 'r', label='$T_{\\rm rx}=100$ K') 
plt.plot(nu, 5*PW*np.log10(S11_2), 'b', label='$T_{\\rm rx}=75$ K') 
plt.xlim(.05,.25)#; plt.ylim(50, 5000)
plt.xlabel('Frequency [GHz]')
plt.ylabel('$S_{11}$ [dB]')
plt.legend(loc='best')
plt.grid()
plt.show()

f_obs = 1.5
e_dish = 0.74
Trx_lotarg = np.around(e_dish * (500 - Tsky))
Tsys = np.sqrt(f_obs) * Tsky
Trx_hitarg = e_dish * (Tsys - Tsky)
Trx_targ = np.where(Trx_lotarg > Trx_hitarg, Trx_lotarg, Trx_hitarg).clip(0)
Trx_targ = np.where(Trx1 * np.ones_like(nu) > Trx_targ, Trx1 * np.ones_like(nu), Trx_targ).clip(0)
plt.semilogy(nu, Tsky, 'b', label='$T_{\\rm sky}$')
plt.semilogy(nu, Trx_lotarg, 'g:', label='$T_{\\rm rx}$ low-freq target') 
plt.semilogy(nu, Trx_hitarg, 'g--', label='$T_{\\rm rx}$ high-freq target') 
plt.semilogy(nu, Trx_targ, 'r', label='$T_{\\rm rx}$ target')
plt.semilogy(nu, Trx1 * np.ones_like(nu), 'r--', label='$T_{\\rm rx}=100$ K') 
plt.semilogy(nu, Trx2 * np.ones_like(nu), 'k--', label='$T_{\\rm rx}=75$ K') 
plt.xlim(.05,.25); plt.ylim(50, 5000)
plt.xlabel('Frequency [GHz]')
plt.ylabel('Temperature [K]')
plt.legend(loc='best')
plt.show()

S11_1 = (1 - Trx1 / Trx_targ).clip(0,1)
S11_2 = (1 - Trx2 / Trx_targ).clip(0,1)
S11_3 = (1 - 85 / Trx_targ).clip(0,1)
plt.plot(nu, 5*PW*np.log10(S11_1), 'r', label='$T_{\\rm rx}=100$ K') 
plt.plot(nu, 5*PW*np.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
plt.plot(nu, 5*PW*np.log10(S11_2), 'c', label='$T_{\\rm rx}=75$ K') 
plt.xlim(.05,.25); plt.ylim(-10*PW, 0)
plt.xlabel('Frequency [GHz]') 
plt.ylabel('$S_{11}$ [dB]; PW=%i'%PW)
plt.legend(loc='best')
plt.grid()
plt.show()


plt.plot(nu, 5*PW*np.log10(S11_3), 'b', label='$T_{\\rm rx}=85$ K') 
plt.xlim(.05,.25); plt.ylim(-5*PW, 0)
plt.xlabel('Frequency [GHz]') 
plt.ylabel('$S_{11}$ [dB]; PW=%i'%PW) 
plt.grid()
plt.show()