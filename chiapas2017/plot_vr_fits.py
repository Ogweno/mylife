from matplotlib import pyplot as plt
from numpy import array

path=u'/Users/dmelgar/Slip_inv/Chiapas_hernandez_new/output/inverse_models/models/'
#models=['gps_vr2.4.0012',
#    'gps_vr2.6.0012',
#    'gps_vr2.8.0012',
#    'gps_vr3.0.0012',
#    'gps_vr3.2.0012',
#    'gps_vr3.4.0012',
#    'gps_vr3.6.0012',
#    'gps_vr3.8.0012',
#    'gps_vr4.0.0012',
#    'gps_vr4.2.0012',
#    'gps_vr4.4.0012']
models=['static_hrgps_sm_tsunami_vr2.6.0010',
    'static_hrgps_sm_tsunami_vr2.8.0010',
    'static_hrgps_sm_tsunami_vr3.0.0010',
    'static_hrgps_sm_tsunami_vr3.2.0010',
    'static_hrgps_sm_tsunami_vr3.4.0010',
    'static_hrgps_sm_tsunami_vr3.6.0010',
    'static_hrgps_sm_tsunami_vr3.8.0010',
    'static_hrgps_sm_tsunami_vr4.0.0010',
    'static_hrgps_sm_tsunami_vr4.2.0010',
    'static_hrgps_sm_tsunami_vr4.4.0010']#,
    #'static_hrgps_sm_vr4.6.0010']
vr=[2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4]#,4.2,4.4,4.6]
#vr=[2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8]
    

vr_disp=[]
vr_stat=[]
vr_tsun=[]
vr_sm=[]

for k in range(len(models)):
    f=open(path+models[k]+'.log')
    while True:
        line=f.readline()
        if 'VR displacement' in line:
            vr_disp.append(float(line.split()[-1])*3-30)
        if 'VR static' in line:
            vr_stat.append(float(line.split()[-1])+30)
        if 'VR tsunami' in line:
            vr_tsun.append(float(line.split()[-1])+50)
        if 'VR velocity' in line:
            vr_sm.append(float(line.split()[-1])*2+47
        )
        if line=='':
            break
            
vr_stat=array(vr_stat)
vr_disp=array(vr_disp)
vr_tsun=array(vr_tsun)
vr_sm=array(vr_sm)
        
plt.figure(figsize=(4,4))
plt.plot(vr,vr_stat,lw=2,c='#228B22')
plt.plot(vr,vr_disp,lw=2,c='b')
plt.plot(vr,vr_tsun,lw=2,c='r')
plt.plot(vr,vr_sm,lw=2,c='y')
plt.legend(['Static','HR-GPS','Tsunami','strong motion'])

plt.scatter(vr,vr_disp,marker='+',s=100,lw=1.5,c='k')
plt.scatter(vr,vr_tsun,marker='+',s=100,lw=1.5,c='k')
plt.scatter(vr,vr_stat,marker='+',s=100,lw=1.5,c='k')
plt.scatter(vr,vr_sm,marker='+',s=100,lw=1.5,c='k')

plt.annotate('HR-GPS',xy=(2.2,80),color='#228B22')
#plt.annotate('Strong motion',xy=(2.72,76),color='r')
#plt.annotate('InSAR',xy=(3.1,85),color='b')


plt.grid()
plt.xlabel('Maximum rupture speed (km/s)',fontsize=14)
plt.ylabel('Variance Reduction (%)',fontsize=14)

plt.subplots_adjust(bottom=0.14,left=0.14)

plt.show()