from glob import glob
from numpy import zeros,array,r_,expand_dims
from obspy import Stream,Trace
from obspy import UTCDateTime
from matplotlib import pyplot as plt



files=glob(u'/Users/dmelgar/Chiapas2017/strong_motion/II/raw/*')
path=u'/Users/dmelgar/Chiapas2017/strong_motion/II/sac/'

for k in range(0,len(files)):
    with open(files[k]) as fp:
        for i, line in enumerate(fp):
            if i == 67:
                tstart=line.split()[-1]
                tstart=UTCDateTime('2017-09-08T'+tstart)
            elif i > 67:
                break
    f=open(files[k])
    lnum=0
    
    #What station am I on?
    sta=files[k].split('/')[-1][0:4]
    
    nsamples=1e6
    strout=''
    while True:
        line=f.readline()
        if line=='':
            break
        if lnum==16:
            strout+=line.split()[-1]
        if lnum==22:
            strout=strout+'\t'+line.split()[5]
        if lnum==23:
            strout=strout+'\t'+line.split()[1]
        if lnum==46:
            dt=float(line.split('/')[-1])
            #print dt
        if lnum==71:
            nsamples=int(line.split('/')[-1])
            #print nsamples
        elif lnum==107:
            components=line.split()
            print sta+' '+str(components)
        elif lnum==109:
            data=expand_dims(array(line.split()).astype('float'),0)
        elif lnum>109:
            data=r_[data,expand_dims(array(line.split()).astype('float'),0)]
        if lnum>=nsamples-1:
            break
        lnum+=1
        
    #print strout
    
    
    #put in streams
    n=Stream(Trace())
    e=Stream(Trace())
    z=Stream(Trace())
    
    for kcomp in range(3):
        if components[kcomp]=='ENE' or components[kcomp]=='N90E':
            e[0].data=data[:,kcomp]/100.
        elif components[kcomp]=='ENN' or components[kcomp]=='N00E':
            n[0].data=data[:,kcomp]/100.
        elif components[kcomp]=='ENZ' or components[kcomp]=='V':
            z[0].data=data[:,kcomp]/100.
        else:
            print ': Unknown component'
            
    n[0].stats.starttime=tstart ; n[0].stats.delta=dt
    e[0].stats.starttime=tstart ; e[0].stats.delta=dt
    z[0].stats.starttime=tstart ; z[0].stats.delta=dt
    
    n.write(path+sta+'.HNN.sac',format='SAC')
    e.write(path+sta+'.HNE.sac',format='SAC')
    z.write(path+sta+'.HNZ.sac',format='SAC')
    
    #plt.figure()
    #plt.subplot(311)
    #plt.plot(n[0].times(),n[0].data)
    #plt.ylabel('North (m/s/s)')
    #plt.subplot(312)
    #plt.plot(e[0].times(),e[0].data)
    #plt.ylabel('East (m/s/s)')
    #plt.subplot(313)
    #plt.plot(z[0].times(),z[0].data)
    #plt.ylabel('Up (m/s/s)')
    #plt.xlabel('Seconds')
    #plt.suptitle('Station '+sta)
    #plt.savefig('/Users/dmelgar/Chiapas2017/strong_motion/II/plots/'+sta+'.png')
    #plt.close()
    
    
    #Low pass filter and decimate
    
    
    
    
