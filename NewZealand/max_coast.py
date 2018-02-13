from numpy import genfromtxt,zeros,where,c_,savetxt
from glob import glob
from matplotlib import pyplot as plt

path='/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/'
files=glob(path+'gauge01*')
fout='/Users/dmelgar/code/GMT/NewZealand/max_eta.txt'

# get coordiantes
lon=zeros(len(files))
lat=zeros(len(files))
max_eta=zeros(len(files))
for k in range(len(files)):
    
    if k%20==0:
        print k
    
    f=open(files[k])
    line=f.readline()
    lon[k]=float(line.split()[4])
    lat[k]=float(line.split()[5])
    f.close()
    
    #read full gauge
    g=genfromtxt(files[k])
    g0=g[0,5]
    t=g[:,1]
    eta=g[:,5]
    
    #Keep only stuff above level 3
    i=where(g[:,0]>=3)[0]
    t=t[i]
    eta=eta[i]
    aux=g[i,2]
    
    #Other filter
    i=where(aux>0)
    t=t[i]
    eta=eta[i]
    
    i=where(eta>0)[0]
    if len(i)>0:
        max_eta[k]=eta[i].max()
        


s=genfromtxt(u'/Users/dmelgar/NewZealand2016/tsunami/survey.txt')

plt.figure(figsize=(3,12))
plt.plot(max_eta,lat)
#plt.scatter(s[:,2],s[:,1])
#wellington
plt.scatter(0.4,-41.82,c='r')
plt.legend(['Model','Tide guges'],frameon=False,loc=4)
plt.annotate(s='Wellington',xy=(0.4,-41.82),xytext=(0.55,-41.86))
#Kaik
plt.scatter(2.3,-42.41,c='r')
plt.annotate(s='Kaikoura',xy=(2.3,-42.41),xytext=(1.8,-42.30))
#Kaik
plt.scatter(0.25,-40.69992,c='r')
plt.annotate(s='Castlepoint',xy=(0.15,-40.8992),xytext=(0.3,-40.65))

plt.ylabel('Latitude')
plt.xlabel('Max amplitude at coast (m)')
plt.subplots_adjust(left=0.2,bottom=0.06,top=0.98,right=0.98)
plt.show()

savetxt(fout,c_[max_eta,lat],fmt='%.3f')