'''
Diego Melgar, UC Berkeley, 01/2015
Pre-process Iquique strong motion data
'''

from obspy import read
from glob import glob
from numpy import float64,zeros,genfromtxt
from string import ljust
from matplotlib import pyplot as plt

applygain=False
getmeta=False
makeplots=True

path="/Users/dmelgar/Coquimbo2015/strong_motion/"
rawfiles=glob(path+'raw/*.SAC')
pathout=path+'proc/'

if applygain:
    #First gain correct
    for k in range(len(rawfiles)):
        #Read file get station name and channel
        st=read(rawfiles[k])
        sta=st[0].stats.station
        chan=st[0].stats.channel
        #Go find appropriate RESP file
        resp_file=glob(path+'resp/RESP.*'+sta+'..'+chan)[0]
        #Parse for channel sensitivity
        f=open(resp_file,'r')
        found_sensitivity=False
        stage_0=False
        infinite=0
        while True:
            line=f.readline()
            if 'Channel Sensitivity' in line:
                found_sensitivity=True
            if found_sensitivity:
                if 'Stage sequence number:                 0' in line:
                    stage_0=True
                if stage_0==True:
                    if 'Sensitivity' in line:
                        gain=float64(line.split()[-1])
                        print 'Gain is '+str(gain)+' for channel '+chan+' at station '+sta
                        break
            if infinite>5000:
                print 'Infinite loop'
                break
            infinite+=1
        f.close()
        #Apply gain
        st[0].data=st[0].data/gain
        #Correct crappy GFZ naming convention
        if chan=='HLE':
            chan='HNE'
        if chan=='HLN':
            chan='HNN'
        if chan=='HLZ':
            chan='HNZ'
        #Write to file
        st.write(pathout+sta+'.'+chan+'.sac',format='SAC')
    
#Now get metadata
if getmeta:
    rawfiles=glob(path+'proc/*HNE.sac')
    lat=zeros(len(rawfiles))
    lon=zeros(len(rawfiles))
    outfile=open(path+'strong_motion.sta','w')
    for k in range(len(rawfiles)):
        try:
            #Read file get station name and channel
            st=read(rawfiles[k])
            sta=st[0].stats.station
            #Go find appropriate PZ file
            pzfile=glob(path+'pz/*PZs*'+sta+'*')[0]
            f=open(pzfile,'r')
            while True:
                line=f.readline()
                if 'STATION' in line:
                    if k==0:
                        station=[line.split(':')[-1]]
                    else:
                        station.append(line.split(':')[-1])
                if 'DESCRIPTION' in line:
                    if k==0:
                        site=[line.split(':')[-1]]
                    else:
                        site.append(line.split(':')[-1])
                if 'LATITUDE' in line:
                    lat[k]=float64(line.split(':')[-1])
                if 'LONGITUDE' in line:
                    lon[k]=float64(line.split(':')[-1])
                    break
            f.close()
            #Write to file
            ll='%.6f\t%.6f'%(float64(lon[k]),float64(lat[k]))
            sta=ljust(station[k].strip(),5,' ')
            line=ll+'\t'+sta+'\t'+site[k].strip()+'\n'
            outfile.write(line)
        except:
            print 'Error ons tation '+sta
    outfile.close()

#Make plots of data       
if makeplots:
    sta=genfromtxt(u'/Users/dmelgar/Coquimbo2015/station_info/strong_motion.sta',usecols=2,dtype='S')
    for k in range(len(sta)):
        print sta[k]
        n=read('/Users/dmelgar/Coquimbo2015/strong_motion/proc/'+sta[k]+'.HNN.sac')
        e=read('/Users/dmelgar/Coquimbo2015/strong_motion/proc/'+sta[k]+'.HNE.sac')
        z=read('/Users/dmelgar/Coquimbo2015/strong_motion/proc/'+sta[k]+'.HNZ.sac')
        plt.figure()
        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'b')
        plt.grid()
        plt.ylabel(r'North $(m/s^2)$')
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'r')
        plt.grid()
        plt.ylabel(r'East $(m/s^2)$')
        plt.subplot(313)
        plt.plot(z[0].times(),z[0].data,'g')
        plt.grid()
        plt.ylabel(r'Up $(m/s^2)$')
        plt.xlabel(r'Time $(s)$')
        plt.suptitle('Station '+sta[k])
        plt.savefig('/Users/dmelgar/Coquimbo2015/strong_motion/plots/'+sta[k]+'.png')


 
    