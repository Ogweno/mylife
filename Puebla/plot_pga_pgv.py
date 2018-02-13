from numpy import genfromtxt,array,ones,zeros,r_,where,arange
from obspy import read
from pyproj import Geod
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz
from scipy.signal import butter,filtfilt
from mudpy.ruptfunctions import rotatedResponseSpectrum

#hypocenter
epi=array([-98.6878,18.3044,57.52])

fout='/Users/dmelgar/Puebla2017/strong_motion/synthesis/ground_motions_20s_highpass.txt'

#read station list
sta=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/really_final_stations.txt',usecols=1,dtype='S')
lonlat=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/really_final_stations.txt',usecols=[2,3])

#where's teh data?
path=u'/Users/dmelgar/Puebla2017/strong_motion/sac/'

#get event/station distances
P=Geod(ellps='WGS84')
az,baz,d=P.inv(lonlat[:,0],lonlat[:,1],ones(len(lonlat))*epi[0],ones(len(lonlat))*epi[1])
d=d/1000.

fcorner=1./20
order=2

oscFreqs=array([0.01,0.020,0.030,0.050,0.075,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.0])

#initalize
pga=zeros(len(sta))
pgv=zeros(len(sta))
SA=zeros((len(sta),len(oscFreqs)))



def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
        
    return data_filt



#get intenisyt parameters
for k in range(len(sta)):
    
    print k
    
    try:
        n=read(path+sta[k]+'.HLN.sac')
    except:
        n=read(path+sta[k]+'.HNN.sac')
    try:
        e=read(path+sta[k]+'.HLE.sac')
    except:
        e=read(path+sta[k]+'.HNE.sac')
    try:
        z=read(path+sta[k]+'.HLZ.sac')
    except:
        z=read(path+sta[k]+'.HNZ.sac')

    
    #get pga
    pga[k]=max(r_[n[0].data,e[0].data,z[0].data])
    
    #Integrate
    fsample=1./n[0].stats.delta
    vn=highpass(cumtrapz(n[0].data,n[0].times(),initial=0),fcorner,fsample,order)
    ve=highpass(cumtrapz(e[0].data,e[0].times(),initial=0),fcorner,fsample,order)
    vz=highpass(cumtrapz(z[0].data,z[0].times(),initial=0),fcorner,fsample,order)
    
    #get pgv
    pgv[k]=max(r_[vn,ve,vz])   
    
    #High pass filter for SA
    fsample=1./n[0].stats.delta
    n[0].data=highpass(n[0].data,fcorner,fsample,order)
    e[0].data=highpass(e[0].data,fcorner,fsample,order)
    z[0].data=highpass(z[0].data,fcorner,fsample,order)
    
    #Get SA
    dt=n[0].stats.delta
    #trim n and e to same starttimes
    t1=max([n[0].stats.starttime,e[0].stats.starttime])
    t2=min([n[0].stats.endtime,e[0].stats.endtime])
    n[0].trim(starttime=t1,endtime=t2)
    e[0].trim(starttime=t1,endtime=t2)
    sa=rotatedResponseSpectrum(dt, n[0].data, e[0].data, oscFreqs, oscDamping=0.05,
        percentiles=[50], angles=arange(0, 180, step=1)) 
    #unwrap
    for kf in range(len(oscFreqs)):
        SA[k,kf]=sa[0][kf][0]
    
    
  
#write to file
f=open(fout,'w')
f.write('# Ground motions for the M7.1 Puebla-Morelos EQ. rotd50 SAs are computed at 21 frequencies (in Hz) they are:\n# 0.010,0.020,0.030,0.050,0.075,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.0\n') 
f.write('#\n')
f.write('# station     lon             lat       depi(km)     pga(g)   pgv(cm/s)  SA in g at 21 frequencies\n')
for k in range(len(sta)):
    line='%s\t%12.4f\t%12.4f\t%7.1f\t%13.4e\t%.4f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n' % (sta[k],lonlat[k,0],lonlat[k,1],d[k],pga[k]/9.81,pgv[k]*100,SA[k,0]/9.81,SA[k,1]/9.81,SA[k,2]/9.81,SA[k,3]/9.81,SA[k,4]/9.81,SA[k,5]/9.81,SA[k,6]/9.81,SA[k,7]/9.81,SA[k,8]/9.81,SA[k,9]/9.81,SA[k,10]/9.81,SA[k,11]/9.81,SA[k,12]/9.81,SA[k,13]/9.81,SA[k,14]/9.81,SA[k,15]/9.81,SA[k,16]/9.81,SA[k,17]/9.81,SA[k,18]/9.81,SA[k,19]/9.81,SA[k,20]/9.81)  
    f.write(line)
f.close()
      
          
plt.figure(figsize=(12,6))
i_vm=[]
i_ig=[]
i_i=[]
for k in range(len(sta)):
    if 'VM' in sta[k]:
        i_vm.append(k)
    else:
        i_ig.append(k)


plt.subplot(231)
plt.loglog(d[i_ig],pga[i_ig]/9.81,marker='o',markersize=5,lw=0)
plt.loglog(d[i_vm],pga[i_vm]/9.81,marker='o',markersize=5,lw=0,c='#CC0066')
plt.legend(['IG or II','VM'])
plt.ylabel('PGA (g)')

plt.subplot(232)
plt.loglog(d[i_ig],pgv[i_ig]/100.,marker='o',markersize=5,lw=0)
plt.loglog(d[i_vm],pgv[i_vm]/100.,marker='o',markersize=5,lw=0,c='#CC0066')
plt.ylabel('PGV (cm/s)')

plt.subplot(233)
plt.loglog(d[i_ig],SA[i_ig,5]/9.81,marker='o',markersize=5,lw=0)
plt.loglog(d[i_vm],SA[i_vm,5]/9.81,marker='o',markersize=5,lw=0,c='#CC0066')
plt.ylabel('SA(T=10s) (g)')

plt.subplot(234)
plt.loglog(d[i_ig],SA[i_ig,11]/9.81,marker='o',markersize=5,lw=0)
plt.loglog(d[i_vm],SA[i_vm,11]/9.81,marker='o',markersize=5,lw=0,c='#CC0066')
plt.xlabel('Epicentral distance (km)')
plt.ylabel('SA(T=2.0s) (g)')

plt.subplot(235)
plt.loglog(d[i_ig],SA[i_ig,15]/9.81,marker='o',markersize=5,lw=0)
plt.loglog(d[i_vm],SA[i_vm,15]/9.81,marker='o',markersize=5,lw=0,c='#CC0066')
plt.xlabel('Epicentral distance (km)')
plt.ylabel('SA(T=0.5s) (g)')

plt.subplot(236)
plt.loglog(d[i_ig],SA[i_ig,18]/9.81,marker='o',markersize=5,lw=0)
plt.loglog(d[i_vm],SA[i_vm,18]/9.81,marker='o',markersize=5,lw=0,c='#CC0066')
plt.xlabel('Epicentral distance (km)')
plt.ylabel('SA(T=0.2s) (g)')

plt.subplots_adjust(left=0.06,right=0.97,top=0.95,bottom=0.1,wspace=0.25)




#plt.loglog(d,pga/9.81,marker='o',markersize=7,lw=0)
#plt.xlabel('Epicentral distance (km)')
#plt.ylabel('PGA (g)')
#for k in range(len(sta)):
#    plt.annotate(s=sta[k],xy=(d[k],pga[k]/9.81))
#
#
#
#plt.figure()
#plt.loglog(d,pgv*100,marker='o',markersize=7,lw=0)
#plt.xlabel('Epicentral distance (km)')
#plt.ylabel('PGV (cm/s)')
#for k in range(len(sta)):
#    plt.annotate(s=sta[k],xy=(d[k],pgv[k]*100))
#
#
#
#
##Comapre to II report
#staii=['PPIG','TLIG','MEIG','TPIG','CAIG','ZIIG','CMIG','MMIG','TUIG']
#pgaii=array([112.62,110.70,74.59,71.34,8.07,4.19,2.99,2.14,2.07])
#plt.figure()
#for k in range(len(staii)):
#    i=where(staii[k]==sta)[0]
#    plt.scatter(pga[i]*100,pgaii[k])
#    if k<4:
#        plt.annotate(s=staii[k],xy=(pga[i]*100+2,pgaii[k]+2))
#plt.xlabel('PGA formas de onda (cm/s/s)')
#plt.ylabel('PGA reporte II (cm/s/s)')
#plt.plot([0,120],[0,120])
#
#
#
##Comapre to II report
#pgvii=array([12.96,3.9,2.42,7.41,0.42,0.50,0.35,0.34,0.33])
#plt.figure()
#for k in range(len(staii)):
#    i=where(staii[k]==sta)[0]
#    plt.scatter(pgv[i]*100,pgvii[k])
#    if k<4:
#        plt.annotate(s=staii[k],xy=(pgv[i]*100+0.1,pgvii[k]+0.1))
#plt.xlabel('PGV formas de onda (cm/s)')
#plt.ylabel('PGV reporte II (cm/s)')
#plt.plot([0,14],[0,14])


plt.show()