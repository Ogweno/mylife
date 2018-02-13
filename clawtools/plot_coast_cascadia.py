# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


fgmax_file=u'/Users/dmelgar/Tsunamis/Cascadia/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/Cascadia/fort.FG1.aux1'

wet_tol=0.001
etaclip=50

#Get maximum amplitude first
lon=genfromtxt(fgmax_file,usecols=0)
lat=genfromtxt(fgmax_file,usecols=1)
amr=genfromtxt(fgmax_file,usecols=2)
H_in=genfromtxt(fgmax_file,usecols=3)
b_in=genfromtxt(aux_file)
unique_amr=unique(amr)
eta=zeros(len(H_in))
H=zeros(len(H_in))
b=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    if unique_amr[k]==0:
        eta[i]=0
    else:
        eta[i]=H_in[i]+b_in[i,unique_amr[k]+1]
        H[i]=H_in[i]
        b[i]=b_in[i,unique_amr[k]+1]
i=where(b>0)[0]
eta[i]=nan


dlat=0.25
#Make plot
plt.figure(figsize=(2, 6))
static=plt.plot(eta,lat,lw=0.5,c='k',zorder=3,label='Predicted')

xfill=array([0,0.3])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#32CD32',alpha=1.0)
xfill=array([0.3,1.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FFD700',alpha=1.0)
xfill=array([1.,3.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FF8C00',alpha=1.0)
xfill=array([3.,9.8])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#DC143C',alpha=1.0)

plt.ylim([39,51])
plt.xlim([0,7])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.show()

