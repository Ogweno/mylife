'''
Strong motion processing script for Napa event
First define functions that will be used repeatedly, batch run is further down in this file
'''

#Common modules
from numpy import genfromtxt,where,intersect1d,array,mean,float64
from obspy import read
from obspy.core import UTCDateTime
from glob import glob
from matplotlib import pyplot as plt
from scipy.signal import butter,filtfilt,cheby1
import re
from mudpy.forward import lowpass as lfilter
from scipy.integrate import cumtrapz
from datetime import timedelta


#Fucntion definitions


def apply_gain(acc_file,acc_out_path):
    '''
    Read acc file, apply gain and re-write.
    '''
    #Read data file
    st=read(acc_file)
    #Decide which gain file to use
    if st[0].stats.network=='BK':
        gain_file=gain_file1
    elif st[0].stats.network=='NC':
        gain_file=gain_file2
    elif st[0].stats.network=='NP':
        gain_file=gain_file3
    elif st[0].stats.network=='CE':
        gain_file=gain_file4
    #Parse gain_file for correct gain
    nets=genfromtxt(gain_file,usecols=0,dtype='S')
    stas=genfromtxt(gain_file,usecols=1,dtype='S')
    chans=genfromtxt(gain_file,usecols=2,dtype='S')
    gains=genfromtxt(gain_file,usecols=10)
    #Find correct row to read from
    inet=where(nets==st[0].stats.network)[0]
    ista=where(stas==st[0].stats.station)[0]
    ichan=where(chans==st[0].stats.channel)[0]
    i=intersect1d(inet,intersect1d(ista,ichan))
    if len(i)>1: #Several gains, keep latest
        i=i[-1]
        gain=gains[i]
    elif len(i)==1: #Only one gain, it's all good
        gain=gains[i]
    elif len(i)==0: #No simple gain, get it from sac pole-zero file, look at "sensitivity" field
        print '...no simple gain file, looking for SAC pole-zero file'
        sta=acc_files[0].split('/')[-1].split('.')[1]
        chan=acc_files[0].split('/')[-1].split('.')[3]
        try:
            pzfile=glob(sacpz_dir+'*'+sta+'*'+chan+'*')[0]
            #Read pzfile and find gain
            f=open(pzfile)
            for line in f:
                if re.match("(.*)SENSITIVITY(.*)",line):
                    gain=float(line.split(':')[1].split(' ')[1])
        except:
            print '... no gain information in this pole-zero file either, sorry.'
            return
    #Apply gain
    st[0].data=st[0].data/gain
    #Remove pre-event mean
    st[0].data=st[0].data-mean(st[0].data[0:500])
    st.write(acc_out_path+st[0].stats.network+'.'+st[0].stats.station+'.'+st[0].stats.channel+'.sac',format='SAC')

def stdecimate(st,factor,order):
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout
    
def bfilter(tr,fcorner,order):
    b, a = butter(order, fcorner,btype='bandpass')
    y = filtfilt(b, a, tr)
    return y
    
def lfilter(tr,fcorner,order):
    b, a = butter(order, fcorner,btype='lowpass')
    y = filtfilt(b, a, tr)
    return y
    
def hfilter(tr,fcorner,order):
    b, a = butter(order, fcorner,btype='highpass')
    y = filtfilt(b, a, tr)
    return y
    
def hcfilter(tr,fcorner,order,rp):
    b, a = cheby1(order,rp, fcorner,btype='highpass')
    y = filtfilt(b, a, tr)
    return y

################################################################################
#
#                         Batch processing section
#
################################################################################
    
#Places to find things
master_list='/Users/dmelgar/Napa2014/seis_only/stations.txt'
acc_dir='/Users/dmelgar/Napa2014/seis_only/raw/'
acc_out_path='/Users/dmelgar/Napa2014/seis_only/zeroth/'
filt_out_path='/Users/dmelgar/Napa2014/seis_only/filtered/'
plot_out_path='/Users/dmelgar/Napa2014/seis_only/plots/'
gain_file1='/Users/dmelgar/Napa2014/acc/_resp.txt'
gain_file2='/Users/dmelgar/Napa2014/acc/_resp.NC.txt'
gain_file3='/Users/dmelgar/Napa2014/acc/_resp.NP.txt'
gain_file4='/Users/dmelgar/Napa2014/acc/_resp.CE.txt'
    
#What do you want to do?
applygain_flag=False  #Comvert from counts to m/s^2
decimate=True
bandpass=False
trim=True
make_plot=True #Make and save plots to file

#Parameters
fcorner=array([1./50,0.5])
time_epi=UTCDateTime('2014-08-24T10:20:44')
t1=timedelta(seconds=10)
t2=timedelta(seconds=60)
    
#Read collocated stations list file
acc_list=genfromtxt(master_list,usecols=0,dtype='S')
trim_time=90

if applygain_flag==True:
    print 'Applying gain to accelerometer data...'
    for k in range(len(acc_list)):
        print '... working on '+acc_list[k]
        acc_files=glob(acc_dir+'*'+acc_list[k]+'*.HN*') #Should get 3 files here for NEU
        if len(acc_files)<3:
            print '... ERROR: no acceleroemter data found for station '+acc_list[k]
        else:
            apply_gain(acc_files[0],acc_out_path)
            apply_gain(acc_files[1],acc_out_path)
            apply_gain(acc_files[2],acc_out_path)
          
if bandpass==True:
        print 'Bandpass filtering ...'
        for k in range(len(acc_list)):
            print '...'+acc_list[k]
            #Read accelerometer data
            ae=read(glob(acc_out_path+acc_list[k]+'*.HNE*')[0])
            an=read(glob(acc_out_path+acc_list[k]+'*.HNN*')[0])
            az=read(glob(acc_out_path+acc_list[k]+'*.HNZ*')[0])
            #Copy to vel
            ve=ae.copy()
            vn=an.copy()
            vz=az.copy()
            #Integrate
            ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
            vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
            vz[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
            #Normalize corner
            fnyquist=1./(2*ve[0].stats.delta)
            Fc=fcorner/fnyquist
            #Filter
            ve[0].data=bfilter(ve[0].data,Fc,3)
            vn[0].data=bfilter(vn[0].data,Fc,3)
            vz[0].data=bfilter(vz[0].data,Fc,3)
            if trim==True:
                tstart=ve[0].stats.starttime
                tp=timedelta(seconds=float64(ve[0].stats.sac['t0'])) #pwave pick         
                ve[0].trim(starttime=tstart+tp-t1,endtime=tstart+tp+t2)
                vn[0].trim(starttime=tstart+tp-t1,endtime=tstart+tp+t2)
                vz[0].trim(starttime=tstart+tp-t1,endtime=tstart+tp+t2)
                
            #Write
            vn.write(filt_out_path+vn[0].stats.station+'.HYN.sac',format='SAC')
            ve.write(filt_out_path+ve[0].stats.station+'.HYE.sac',format='SAC')
            vz.write(filt_out_path+vz[0].stats.station+'.HYZ.sac',format='SAC')



if decimate==True:
        print 'Decimating and high-pass filtering ...'
        for k in range(len(acc_list)):
            print '...'+acc_list[k]
            #Read accelerometer data
            ae=read(glob(acc_out_path+acc_list[k]+'*.HNE*')[0])
            an=read(glob(acc_out_path+acc_list[k]+'*.HNN*')[0])
            az=read(glob(acc_out_path+acc_list[k]+'*.HNZ*')[0])
            #Copy to vel
            ve=ae.copy()
            vn=an.copy()
            vz=az.copy()
            #Integrate
            ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
            vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
            vz[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
            #Decimate
            if vn[0].stats.delta==0.005:
                vn=stdecimate(vn,5,13) ; vn=stdecimate(vn,5,13) ; vn=stdecimate(vn,2,13)
                ve=stdecimate(ve,5,13) ; ve=stdecimate(ve,5,13) ; ve=stdecimate(ve,2,13)
                vz=stdecimate(vz,5,13) ; vz=stdecimate(vz,5,13) ; vz=stdecimate(vz,2,13)
            if vn[0].stats.delta==0.01:
                vn=stdecimate(vn,5,13) ; vn=stdecimate(vn,5,13)
                ve=stdecimate(ve,5,13) ; ve=stdecimate(ve,5,13)
                vz=stdecimate(vz,5,13) ; vz=stdecimate(vz,5,13)
            #Normalize corner
            fnyquist=1./(2*ve[0].stats.delta)
            Fc=fcorner/fnyquist
            #Filter
            #ve[0].data=hfilter(ve[0].data,Fc,10)
            #vn[0].data=hfilter(vn[0].data,Fc,10)
            #vz[0].data=hfilter(vz[0].data,Fc,10)
            ve[0].data=bfilter(ve[0].data,Fc,4)
            vn[0].data=bfilter(vn[0].data,Fc,4)
            vz[0].data=bfilter(vz[0].data,Fc,4)
            if trim==True:
                tstart=ve[0].stats.starttime
                tp=timedelta(seconds=float64(ve[0].stats.sac['t0'])) #pwave pick         
                ve[0].trim(starttime=tstart+tp-t1,endtime=tstart+tp+t2,pad=True,fill_value=0)
                vn[0].trim(starttime=tstart+tp-t1,endtime=tstart+tp+t2,pad=True,fill_value=0)
                vz[0].trim(starttime=tstart+tp-t1,endtime=tstart+tp+t2,pad=True,fill_value=0)
            #Zero out things before p arrival
            #i=where(ve[0].times()<t1.seconds)[0]
            #ve[0].data[i]=0
            #i=where(vn[0].times()<t1.seconds)[0]
            #vn[0].data[i]=0
            #i=where(vz[0].times()<t1.seconds)[0]
            #vz[0].data[i]=0
            #Write
            vn.write(filt_out_path+vn[0].stats.station+'.HYN.sac',format='SAC')
            ve.write(filt_out_path+ve[0].stats.station+'.HYE.sac',format='SAC')
            vz.write(filt_out_path+vz[0].stats.station+'.HYZ.sac',format='SAC')
        



if make_plot == True:
        print 'Making plots ...'
        for k in range(len(acc_list)):
            print '...'+acc_list[k]
            #Read accelerometer data
            ve=read(glob(filt_out_path+acc_list[k].split('.')[1]+'*.HYE*')[0])
            vn=read(glob(filt_out_path+acc_list[k].split('.')[1]+'*.HYN*')[0])
            vz=read(glob(filt_out_path+acc_list[k].split('.')[1]+'*.HYZ*')[0])
            #Plot dispalcememnts
            plt.subplot(311)
            plt.plot(vn[0].times(),vn[0].data)
            plt.ylabel('North')
            plt.grid()
            plt.subplot(312)
            plt.plot(ve[0].times(),ve[0].data)
            plt.ylabel('East')
            plt.grid()
            plt.subplot(313)
            plt.plot(vz[0].times(),vz[0].data)
            plt.ylabel('Vertical')
            plt.grid()
            plt.suptitle(acc_list[k]+' velocity (m/s)')
            plt.savefig(plot_out_path+acc_list[k]+'.png',dpi=300)
            plt.close("all")
        
        

