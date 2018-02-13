from numpy import genfromtxt,zeros,savetxt
from matplotlib import pyplot as plt

sta=genfromtxt('/Users/dmelgar/Wenchuan2008/wenc.full.neu.sta',usecols=0,dtype='S')
lon_lat=genfromtxt('/Users/dmelgar/Wenchuan2008/wenc.full.neu.sta',usecols=[1,2])
neu_all=genfromtxt('/Users/dmelgar/Wenchuan2008/wenc.full.neu')
outpath='/Users/dmelgar/Slip_inv/wenc_2008/data/statics/'
#Intialize plotting variables
lonout=zeros(len(sta))
latout=zeros(len(sta))
eout=zeros(len(sta))
nout=zeros(len(sta))
#Make individual neu files
for k in range(len(sta)):
    neu=neu_all[3*k:3*k+3]
    lonout[k]=lon_lat[k,0]
    latout[k]=lon_lat[k,1]
    nout[k]=neu[0]
    eout[k]=neu[1]
    #Write individual .neu file
    savetxt(outpath+sta[k]+'.neu',neu)
#Plot to check everything is ok
fault=genfromtxt(u'/Users/dmelgar/Slip_inv/wenc_2008/data/model_info/wenc_unilateral.fault')
plt.figure()
plt.scatter(fault[:,1],fault[:,2])
plt.quiver(lonout,latout,eout,nout,2)
plt.show()
    
    
    