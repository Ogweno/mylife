# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,r_
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


#fgmax_file=u'/Users/dmelgar/Tsunamis/maulong_gps/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maulong_gps/_output/fort.FG1.aux1'
fgmax_file3=u'/Users/dmelgar/Tsunamis/maulong_pgd/_output/fort.FG1.valuemax'
aux_file3='/Users/dmelgar/Tsunamis/maulong_pgd/_output/fort.FG1.aux1'
fgmax_file1=u'/Users/dmelgar/Tsunamis/maulong_gps_tg/_output/fort.FG1.valuemax'
aux_file1='/Users/dmelgar/Tsunamis/maulong_gps_tg/_output/fort.FG1.aux1'
fgmax_file4=u'/Users/dmelgar/Tsunamis/maulong_wphase/_output/fort.FG1.valuemax'
aux_file4='/Users/dmelgar/Tsunamis/maulong_wphase/_output/fort.FG1.aux1'
#fgmax_file2=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static/_output/fort.FG1.valuemax'
#aux_file2='/Users/dmelgar/Tsunamis/maule_gps_tg_static/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.aux1'


#fgmax_file=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maule_gps_tg_static/_output/fort.FG1.aux1'
#fgmax_file2=u'/Users/dmelgar/Tsunamis/maule_gps_tg_kinematic/_output/fort.FG1.valuemax'
#aux_file2='/Users/dmelgar/Tsunamis/maule_gps_tg_kinematic/_output/fort.FG1.aux1'
fgmax_file2=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.aux1'

survey_file='/Users/dmelgar/Maule2010/tsunami/survey.txt'
#gauges_file=
blend=False
wet_tol=0.001
etaclip=11

#Get maximum amplitude first
lon=genfromtxt(fgmax_file1,usecols=0)
lat1=genfromtxt(fgmax_file1,usecols=1)
amr=genfromtxt(fgmax_file1,usecols=2)
H_in=genfromtxt(fgmax_file1,usecols=3)
b_in=genfromtxt(aux_file1)
unique_amr=unique(amr)
eta1=zeros(len(H_in))
H=zeros(len(H_in))
b=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    if unique_amr[k]==0:
        eta1[i]=0
    else:
        eta1[i]=H_in[i]+b_in[i,unique_amr[k]+1]
        H[i]=H_in[i]
        b[i]=b_in[i,unique_amr[k]+1]
i=where(b>0)[0]
eta1[i]=nan

#Get maximum amplitude first
lon2=genfromtxt(fgmax_file2,usecols=0)-360
lat2=genfromtxt(fgmax_file2,usecols=1)
amr=genfromtxt(fgmax_file2,usecols=2)
H_in=genfromtxt(fgmax_file2,usecols=3)
b_in=genfromtxt(aux_file2)
unique_amr=unique(amr)
eta2=zeros(len(H_in))
b=zeros(len(H_in))
H=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    eta2[i]=H_in[i]+b_in[i,unique_amr[k]+1]
    H[i]=H_in[i]
    b[i]=b_in[i,unique_amr[k]+1]
i=where(H<wet_tol)[0]
eta2[i]=nan
#Now combine
latbounds=[-40,-32.5]
i=where(lat1<latbounds[0])[0] #THings south of the fine model
lat_out=lat1[i]
lon_out=lon[i]
eta_out=eta1[i]
#THings IN the fine model
lat_out=r_[lat_out,lat2]
lon_out=r_[lon_out,lon2]
eta_out=r_[eta_out,eta2]
#Things north of the fine model
i=where(lat1>latbounds[1])[0] #THings north of the fine model
lat_out=r_[lat_out,lat1[i]]
lon_out=r_[lon_out,lon[i]]
eta_out=r_[eta_out,eta1[i]]
#Ok done
lat1=lat_out
lon=lon_out
eta1=eta_out
i=where(eta1>etaclip)[0]
eta1[i]=nan

#Get maximum amplitude first
lon=genfromtxt(fgmax_file3,usecols=0)
lat3=genfromtxt(fgmax_file3,usecols=1)
amr=genfromtxt(fgmax_file3,usecols=2)
H_in=genfromtxt(fgmax_file3,usecols=3)
b_in=genfromtxt(aux_file3)
unique_amr=unique(amr)
eta3=zeros(len(H_in))
H=zeros(len(H_in))
b=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    if unique_amr[k]==0:
        eta3[i]=0
    else:
        eta3[i]=H_in[i]+b_in[i,unique_amr[k]+1]
        H[i]=H_in[i]
        b[i]=b_in[i,unique_amr[k]+1]
i=where(b>0)[0]
eta3[i]=nan

#Get maximum amplitude first
lon=genfromtxt(fgmax_file4,usecols=0)
lat4=genfromtxt(fgmax_file4,usecols=1)
amr=genfromtxt(fgmax_file4,usecols=2)
H_in=genfromtxt(fgmax_file4,usecols=3)
b_in=genfromtxt(aux_file4)
unique_amr=unique(amr)
eta4=zeros(len(H_in))
H=zeros(len(H_in))
b=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    if unique_amr[k]==0:
        eta4[i]=0
    else:
        eta4[i]=H_in[i]+b_in[i,unique_amr[k]+1]
        H[i]=H_in[i]
        b[i]=b_in[i,unique_amr[k]+1]
i=where(b>0)[0]
eta4[i]=nan

#Parse survey data
latsurvey=genfromtxt(survey_file,usecols=1)
etasurvey=genfromtxt(survey_file,usecols=3)
#    
#Get maximum gauge amplitudes
st=read(u'/Users/dmelgar/Maule2010/wave_gauges/ancu.tsun')
gaugemax=st[0].data.max()
st=read(u'/Users/dmelgar/Maule2010/wave_gauges/cald.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Maule2010/wave_gauges/coqu.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Maule2010/wave_gauges/corr.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Maule2010/wave_gauges/talc.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Maule2010/wave_gauges/valp.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
gaugelat=genfromtxt('/Users/dmelgar/Maule2010/wave_gauges/wave_gauge_list.sta',usecols=1)


#Make plot
plt.figure(figsize=(6.5, 10.5))
e1=plt.plot(eta1,lat1,lw=1.0,c='#000000',zorder=1,label='Land+ocean')
e3=plt.plot(eta3,lat3,lw=1.0,c='#6666FF',zorder=2,label='Magnitude')
e4=plt.plot(eta4,lat4,lw=1.0,c='#006400',zorder=3,label='Moment tensor')
tgauges=plt.scatter(gaugemax,gaugelat,marker='*',s=160,c='#1E90FF',zorder=4,label='Tide gauges')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2),
    ncol=1, fancybox=True, shadow=True)
xfill=array([0,0.3])
plt.fill_between(xfill,lat1.min(),lat1.max(),facecolor='#32CD32',alpha=0.5)
xfill=array([0.3,1.])
plt.fill_between(xfill,lat1.min(),lat1.max(),facecolor='#FFD700',alpha=0.5)
xfill=array([1.,3.])
plt.fill_between(xfill,lat1.min(),lat1.max(),facecolor='#FF8C00',alpha=0.7)
xfill=array([3.,99])
plt.fill_between(xfill,lat1.min(),lat1.max(),facecolor='#DC143C',alpha=0.8)

plt.ylim([lat1.min(),lat1.max()])
plt.xlim([0,10.5])
#plt.xlim([0,etasurvey.max()+0.1*etasurvey.max()])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.show()

