import sys
sys.path.append('/Users/dmelgar/Slip_inv/run/')
#from nepal_ttests_fwd import batch
#from nepal_ttests_inv_batch import batch
from nepal_1s_batch_inv import batch
from numpy import arange

#project='Nepal_ttests_7'
#run=['7s_vr2.8','7s_vr3.0','7s_vr3.2','7s_vr3.4','7s_vr3.6']
#rupt=['ALOS_LRGPS_justalos.0021_vr2.8.rupt','ALOS_LRGPS_justalos.0021_vr3.0.rupt','ALOS_LRGPS_justalos.0021_vr3.2.rupt','ALOS_LRGPS_justalos.0021_vr3.4.rupt','ALOS_LRGPS_justalos.0021_vr3.6.rupt']
#vr=arange(2.8,3.8,0.2)
#for k in range(len(vr)):
#    print vr[k]
#    batch(project,run[k],rupt[k],vr[k])
   
   
#def move_files(trise):
#    from glob import glob
#    from obspy import read
#    files=glob(u'/Users/dmelgar/Slip_inv/Nepal_ttests_'+str(trise)+'/output/forward_models/*vr3.2*.vel*')
#    dirchop=u'/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/data/waveforms/'
#    dirout='/Users/dmelgar/Slip_inv/Nepal_ttests_inv/data/waveforms/'
#    for k in range(len(files)):
#        file_out=dirout+files[k].split('/')[-1].split('.')[2]+'.'+files[k].split('/')[-1].split('.')[3]+'.'+files[k].split('/')[-1].split('.')[4]
#        file_chop=dirchop+files[k].split('/')[-1].split('.')[2]+'.'+files[k].split('/')[-1].split('.')[3]+'.'+files[k].split('/')[-1].split('.')[4]
#        st=read(files[k])
#        stchop=read(file_chop)
#        st[0].trim(starttime=stchop[0].stats.starttime,endtime=stchop[0].stats.endtime)
#        st.write(file_out,format='SAC')
#     
#run_name=['6s_vr2.8','6s_vr3.4','6s_vr3.6','2s_vr3.2','4s_vr3.2','8s_vr3.2' ,'10s_vr3.2']
#trise=[6,6,6,2,4,8,10]
#rupt_speed=[2.8,3.4,3.6,3.2,3.2,3.2,3.2]
#for k in range(7):
#    move_files(trise[k])
#    batch(run_name[k],rupt_speed[k])

run_name=['review_3.1_low','review_3.3_low']
rupt_speed=[3.1,3.3]
for k in range(2):
    batch(run_name[k],rupt_speed[k])