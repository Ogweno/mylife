from mudpy import view
from numpy import array,genfromtxt
from matplotlib import pyplot as plt

hf=genfromtxt('/Users/dmelgar/Coquimbo2015/BP/20s/HF.stf')
lf=genfromtxt('/Users/dmelgar/Coquimbo2015/BP/20s/lF.stf')
epicenter=array([-71.654,-31.570,29.8]) 
t,s=view.source_time_function(u'/Users/dmelgar/Slip_inv/Coquimbo_4s/output/inverse_models/models/gps_sm_tg_insar_9win_vel2.0_2insar_highWinsar.0009.inv',epicenter)


plt.figure()
plt.fill(t,s,'b',alpha=0.5)
plt.plot(hf[:,0],hf[:,1]*s.max(),'k')
plt.plot(lf[:,0],lf[:,1]*s.max(),'r')
plt.xlim([0,120])
plt.show()