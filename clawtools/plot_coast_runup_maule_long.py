# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,r_,argmin
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


fgmax_file=u'/Users/dmelgar/Tsunamis/maulong_gps/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/maulong_gps/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/maulong_pgd/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maulong_pgd/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/maulong_gps_tg/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maulong_gps_tg/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/maulong_wphase/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maulong_wphase/_output/fort.FG1.aux1'
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
blend=True
wet_tol=0.001
etaclip=11

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

if blend==True:
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
    i=where(lat<latbounds[0])[0] #THings south of the fine model
    lat_out=lat[i]
    lon_out=lon[i]
    eta_out=eta[i]
    #THings IN the fine model
    lat_out=r_[lat_out,lat2]
    lon_out=r_[lon_out,lon2]
    eta_out=r_[eta_out,eta2]
    #Things north of the fine model
    i=where(lat>latbounds[1])[0] #THings north of the fine model
    lat_out=r_[lat_out,lat[i]]
    lon_out=r_[lon_out,lon[i]]
    eta_out=r_[eta_out,eta[i]]
    #Ok done
    lat=lat_out
    lon=lon_out
    eta=eta_out
else:
    eta_out=eta
i=where(eta_out>etaclip)[0]
eta_out[i]=nan

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

#Find maximum coastline comapred to gauges
gaugemax_pre=zeros(len(gaugemax[0]))
for k in range(len(gaugelat)):
    i=argmin(abs(gaugelat[k]-lat))
    gaugemax_pre[k]=eta[i]

dlat=0.25
#Make plot
plt.figure(figsize=(2, 6))
static=plt.plot(eta,lat,lw=0.5,c='k',zorder=3,label='Predicted')
tgauges=plt.scatter(gaugemax,gaugelat,marker='*',s=160,c='#1E90FF',zorder=4,label='Tide gauges')
#tgauges_obs=plt.barh(gaugelat+dlat,gaugemax[0],height=0.4,color='#1E90FF',zorder=2,label='TG obs.')
tgauges_pre=plt.barh(gaugelat-dlat,gaugemax_pre,height=0.6,color='#909090',zorder=2,label='TG pred.')
#kinematic=plt.plot(eta2,lat2,lw=1.0,c='#000000',zorder=4,label='Coarse (30s)')
#survey=plt.scatter(etasurvey,latsurvey,marker='o',s=40,c='#7FFF00',zorder=1,label='Survey')
#tgauges=plt.scatter(gaugemax,gaugelat,marker='*',s=160,c='#1E90FF',zorder=2,label='Tide gauges')
#plt.legend()
#plt.legend(['Predicted','Survey','Tide Gauges'],loc='upper center', bbox_to_anchor=(0.5, 1.27),
#    ncol=1, fancybox=True, shadow=True)
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),ncol=1, fancybox=True, shadow=True)
xfill=array([0,0.3])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#32CD32',alpha=1.0)
xfill=array([0.3,1.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FFD700',alpha=1.0)
xfill=array([1.,3.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FF8C00',alpha=1.0)
xfill=array([3.,10.5])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#DC143C',alpha=1.0)
#plt.scatter(etasurvey,latsurvey,marker='o',s=40,c='#7FFF00',zorder=5)
#plt.scatter(gaugemax,gaugelat,marker='*',s=160,c='#1E90FF',zorder=6)
plt.ylim([-45,-26.0125])
#plt.ylim([lat.min(),lat.max()])
plt.xlim([0,10.5])
#plt.xlim([0,etasurvey.max()+0.1*etasurvey.max()])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.show()

