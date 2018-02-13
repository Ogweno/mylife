from numpy import array,zeros
from matplotlib import pyplot as plt

num='0016'
path='/Users/dmelgar/Slip_inv/Amatrice_3Dfitsgeol_final1/output/inverse_models/models/_previous/'
root1='bigkahuna_vrtest3win_vr'
root2='.'+num+'.log'

vr=array([1.6,1.8,2.0,2.2,2.4,2.6])
vr_static=zeros(len(vr))
vr_insar=zeros(len(vr))
vr_velocity=zeros(len(vr))
for k in range(len(vr)):
    f=open(path+root1+str(vr[k])+root2,'r')
    while True:
        line=f.readline()
        if 'VR static' in line:
            vr_static[k]=float(line.split('=')[-1])
        elif 'VR velocity' in line:
            vr_velocity[k]=float(line.split('=')[-1])
        elif 'VR InSAR' in line:
            vr_insar[k]=float(line.split('=')[-1])
            break
    f.close()
    
plt.figure()
plt.plot(vr,vr_static+19)
plt.plot(vr,vr_velocity+44)
plt.plot(vr,vr_insar+14)
plt.legend(['GPS','SM','InSAR'],loc=3)
plt.xlabel('vr (km/s)')
plt.ylabel('VR (%)')
plt.show()