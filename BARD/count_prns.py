from glob import glob
from numpy import arange,array,r_
import datetime
from matplotlib import pyplot as plt
from string import rjust

path= u'/Users/dmelgar/BARD/data/2015'
predirs=arange(244,331)

prns=array([])
t=array([])
for k in range(len(predirs)):
    print k
    f=glob(path+'.'+str(predirs[k])+'/brib*.15o')[0]
    brib=open(f,'r')
    #Get current date
    date=datetime.date(2015,1,1) + datetime.timedelta(predirs[k]-1)
    yr=date.year-2000
    mo=date.month
    da=date.day
    #
    search=str(yr)+' '+rjust(str(mo),2,' ')+' '+rjust(str(da),2,' ')
    #Go counting
    while True:
        line=brib.readline()
        if search in line:
            print line
            prns=r_[prns,int(line.split()[-1].split('G')[0])]
            epoch=line.split()
            t=r_[t,float(predirs[k])+float(epoch[3])/24+float(epoch[4])/1440+float(epoch[5])/86400]
        if line=='':
            break
    brib.close()
    
plt.figure()
plt.plot(t,prns,lw=0.5)
plt.xlabel('2015 day of year')
plt.ylabel('# SVs')
plt.xlim([244,332])
plt.plot([273,273],[0,15],'-k',lw=4)
plt.plot([328.5,328.5],[0,15],'-k',lw=4)
plt.title('BRIB, visible satellites')
plt.annotate('Installed Kestrels',xy=(274.5,12))
plt.annotate('New splitter',xy=(320,12))
plt.ylim([0,14])
plt.show()
        
            