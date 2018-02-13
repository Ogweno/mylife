# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


fgmax_file=u'/Users/dmelgar/Tsunamis/maule_gps_static/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/maule_gps_static/_output/fort.FG1.aux1'
fgmax_file2=u'/Users/dmelgar/Tsunamis/maule_gps_kinematic/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/maule_gps_kinematic/_output/fort.FG1.aux1'

#fgmax_file=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maule_gps_tg_static/_output/fort.FG1.aux1'
#fgmax_file2=u'/Users/dmelgar/Tsunamis/maule_gps_tg_kinematic/_output/fort.FG1.valuemax'
#aux_file2='/Users/dmelgar/Tsunamis/maule_gps_tg_kinematic/_output/fort.FG1.aux1'

#survey_file='/Users/dmelgar/Iquique2014/tsunami/iquique_survey.txt'
#gauges_file=

#Get maximum amplitude first
lat=genfromtxt(fgmax_file,usecols=1)
amr=genfromtxt(fgmax_file,usecols=2)
H=genfromtxt(fgmax_file,usecols=3)
b=genfromtxt(aux_file)
unique_amr=unique(amr)
eta=zeros(len(H))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    eta[i]=H[i]+b[i,unique_amr[k]+1]
i=where(eta>etaclip)[0]
eta[i]=nan

#Get maximum amplitude first (again, second file)
lat=genfromtxt(fgmax_file2,usecols=1)
amr=genfromtxt(fgmax_file2,usecols=2)
H=genfromtxt(fgmax_file2,usecols=3)
b=genfromtxt(aux_file2)
unique_amr=unique(amr)
eta2=zeros(len(H))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    eta2[i]=H[i]+b[i,unique_amr[k]+1]
i=where(eta2>etaclip)[0]
eta2[i]=nan

#Parse survey data
#latsurvey=genfromtxt(survey_file,usecols=5)
#etasurvey=genfromtxt(survey_file,usecols=7)
#    
##Get maximum gauge amplitudes
#st=read(u'/Users/dmelgar/Iquique2014/tsunami/gauge_data/aric.sac')
#gaugemax=st[0].data.max()
#st=read(u'/Users/dmelgar/Iquique2014/tsunami/gauge_data/iqui.sac')
#gaugemax=c_[gaugemax,st[0].data.max()]
#st=read(u'/Users/dmelgar/Iquique2014/tsunami/gauge_data/pata.sac')
#gaugemax=c_[gaugemax,st[0].data.max()]
#st=read(u'/Users/dmelgar/Iquique2014/tsunami/gauge_data/pisa.sac')
#gaugemax=c_[gaugemax,st[0].data.max()]
#gaugelat=genfromtxt('/Users/dmelgar/Iquique2014/station_info/tide_gauges.txt',usecols=2)


#Make plot
plt.figure(figsize=(6.5, 10.5))
static=plt.plot(eta,lat,lw=1.6,c='#6666FF',zorder=3,label='Static')
kinematic=plt.plot(eta2,lat,lw=1.0,c='#000000',zorder=4,label='Kinematic')
#survey=plt.scatter(etasurvey,latsurvey,marker='o',s=40,c='#7FFF00',zorder=1,label='Survey')
#tgauges=plt.scatter(gaugemax,gaugelat,marker='*',s=160,c='#1E90FF',zorder=2,label='Tide gauges')
plt.legend(['Instantaneous','Kinematic','Survey','Tide Gauges'])
#plt.legend(['Predicted','Survey','Tide Gauges'],loc='upper center', bbox_to_anchor=(0.5, 1.27),
#    ncol=1, fancybox=True, shadow=True)
xfill=array([0,0.3])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#32CD32',alpha=0.5)
xfill=array([0.3,1.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FFD700',alpha=0.5)
xfill=array([1.,3.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FF8C00',alpha=0.7)
xfill=array([3.,99])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#DC143C',alpha=0.8)
#plt.scatter(etasurvey,latsurvey,marker='o',s=40,c='#7FFF00',zorder=5)
#plt.scatter(gaugemax,gaugelat,marker='*',s=160,c='#1E90FF',zorder=6)
plt.ylim([lat.min(),lat.max()])
plt.xlim([0,8.5])
#plt.xlim([0,etasurvey.max()+0.1*etasurvey.max()])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.show()

