# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,argsort,isnan,argmin
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})

fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_gps_tg_30/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/tohoku_gps_tg_30/_output/fort.FG1.aux1'
fgmax_file2=u'/Users/dmelgar/Tsunamis/tohoku_gps_tg_30_kinematic/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/tohoku_gps_tg_30_kinematic/_output/fort.FG1.aux1'

survey_file='/Users/dmelgar/Iquique2014/tsunami/iquique_survey.txt'

#Get maximum amplitude first
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
i=where(eta>etaclip)[0]
eta[i]=nan
eta_out=eta.copy()
lat_out=lat.copy()


#Get maximum amplitude first (again, second file)
lat2=genfromtxt(fgmax_file2,usecols=1)
amr=genfromtxt(fgmax_file2,usecols=2)
H_in=genfromtxt(fgmax_file2,usecols=3)
b_in=genfromtxt(aux_file2)
unique_amr=unique(amr)
eta2=zeros(len(H_in))
H=zeros(len(H_in))
b=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    if unique_amr[k]==0:
        eta2[i]=0
    else:
        eta2[i]=H_in[i]+b_in[i,unique_amr[k]+1]
        H[i]=H_in[i]
        b[i]=b_in[i,unique_amr[k]+1]
i=where(b>0)[0]
eta2[i]=nan
i=where(eta2>etaclip)[0]
eta2[i]=nan
eta_out2=eta2.copy()
lat_out2=lat2.copy()


#Sort things
i=argsort(lat_out)
lat_out=lat_out[i]
eta_out=eta_out[i]
i=argsort(lat_out2)
lat_out2=lat_out2[i]
eta_out2=eta_out2[i]



#Make plot
i=where(isnan(eta_out)==False)[0]
i2=where(isnan(eta_out2)==False)[0]
plt.figure(figsize=(6, 10))
static=plt.plot(eta_out,lat_out,lw=1.0,c='#6666FF',zorder=3,label='Instantaneous')
kinematic=plt.plot(eta_out2,lat_out2,lw=0.6,c='k',zorder=4,label='Kinematic')
kinematic=plt.plot(abs(eta_out2-eta_out),lat,lw=0.5,c='r',zorder=4,label='Difference')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),ncol=1, fancybox=True, shadow=True)
xfill=array([0,0.3])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#32CD32')
xfill=array([0.3,1.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FFD700')
xfill=array([1.,3.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FF8C00')
xfill=array([3.,20])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#DC143C')

#plt.ylim([lat.min(),lat.max()])
plt.xlim([0,20])
plt.ylim([34,43])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.show()

