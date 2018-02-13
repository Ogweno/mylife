from matplotlib import pyplot as plt
from numpy import array

path=u'/Users/dmelgar/Slip_inv/Coquimbo_4s/output/inverse_models/models/'
models=['gps_sm_tg_insar_veltest_1.6.0000',
    'gps_sm_tg_insar_final1_veltest_1.8.0000',
    'gps_sm_tg_insar_final1_veltest_2.0.0000',
    'gps_sm_tg_insar_final1_veltest_2.2.0000',
    'gps_sm_tg_insar_final1_veltest_2.4.0000',
    'gps_sm_tg_insar_final1_veltest_2.6.0000',
    'gps_sm_tg_insar_final1_veltest_2.8.0000']
#models=['gps_sm_veltest_1.6.0000',
#    'gps_sm_veltest_1.8.0000',
#    'gps_sm_veltest_2.0.0000',
#    'gps_sm_veltest_2.2.0000',
#    'gps_sm_veltest_2.4.0000',
#    'gps_sm_veltest_2.6.0000',
#    'gps_sm_veltest_2.8.0000']
vr=[1.6,1.8,2.0,2.2,2.4,2.6,2.8]
    
vr_stat=[]
vr_disp=[]
vr_vel=[]
vr_insar=[]
vr_tsun=[]
for k in range(len(models)):
    f=open(path+models[k]+'.log')
    while True:
        line=f.readline()
        if 'static' in line:
            vr_stat.append(float(line.split()[-1]))
        if 'displacement' in line:
            vr_disp.append(float(line.split()[-1]))
        if 'VR velocity' in line:
            vr_vel.append(float(line.split()[-1]))
        if 'tsunami' in line:
            vr_tsun.append(float(line.split()[-1]))
        if 'InSAR' in line:
            vr_insar.append(float(line.split()[-1]))
        if line=='':
            break
vr_vel=array(vr_vel)
vr_tsun=array(vr_tsun)
vr_insar=array(vr_insar)
        
plt.figure(figsize=(4,8))
plt.plot(vr,vr_stat,lw=2,c='b')
plt.plot(vr,vr_disp,lw=2,c='#228B22')
plt.plot(vr,vr_vel+35,lw=2,c='r')
plt.plot(vr,vr_tsun+25,lw=2,c='#9932CC')
plt.plot(vr,vr_insar-10,lw=2,c='#808000')

plt.scatter(vr,vr_stat,marker='+',s=100,lw=1.5,c='k')
plt.scatter(vr,vr_disp,marker='+',s=100,lw=1.5,c='k')
plt.scatter(vr,vr_vel,marker='+',s=100,lw=1.5,c='k')
plt.scatter(vr,vr_tsun,marker='+',s=100,lw=1.5,c='k')#
plt.scatter(vr,vr_insar,marker='+',s=100,lw=1.5,c='k')

plt.annotate('InSAR',xy=(2.6,80),color='#808000')
plt.annotate('HR-GPS',xy=(1.9,75),color='#228B22')
plt.annotate('Strong motion',xy=(1.6,62),color='r')
plt.annotate('Tsunami',xy=(2.6,66.5),color='#9932CC')

#plt.legend(['Statics','Displacement','Velocity','Tsunami','InSAR'],loc=4,framealpha=0)
#plt.legend(['Displacement','Velocity','Tsunami'],loc=4)
plt.grid()
plt.xlabel('Maximum rupture speed (km/s)',fontsize=14)
plt.ylabel('Variance Reduction (%)',fontsize=14)

plt.show()