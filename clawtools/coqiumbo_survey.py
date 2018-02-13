# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,r_,argmin,squeeze,isnan
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_detailed/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/coquimbo_detailed/_output/fort.FG1.aux1'
gauge_data=u'/Users/dmelgar/Tsunamis/coquimbo_detailed/_output/fort.gauge'


blend=False
wet_tol=0.001
etaclip=4

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
eta_out=eta

#Wphase and slip inversion clipping
#iclip=where(lat>-25)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#i=where(isnan(eta)==False)[0]
#eta=eta[i]
#lat=lat[i]

#PGD clipping
etaclip=0.6
iclip=where(lat>-28.6)[0]
i=where(eta[iclip]>etaclip)[0]
eta[iclip[i]]=nan
i=where(isnan(eta)==False)[0]
eta=eta[i]
lat=lat[i]
iclip=where(lat<-34.7)[0]
i=where(eta[iclip]>etaclip)[0]
eta[iclip[i]]=nan
i=where(isnan(eta)==False)[0]
eta=eta[i]
lat=lat[i]


#Get tide gauge predictions
stations=genfromtxt('/Users/dmelgar/Coquimbo2015/station_info/tide_gauge_detailed.txt',usecols=3)
gaugelat_pred=genfromtxt('/Users/dmelgar/Coquimbo2015/station_info/tide_gauge_detailed.txt',usecols=1)
gaugedata=genfromtxt(gauge_data)
for kgauge in range(0,4):
    print stations[kgauge]
    i=where(gaugedata[:,0]==stations[kgauge])[0]
    g=gaugedata[i,6]
    if kgauge==0:
        gaugemax_pred=g.max()
    else:
        gaugemax_pred=c_[gaugemax_pred,g.max()]
gaugemax_pred=squeeze(gaugemax_pred)

#Get survey
survey_lat=genfromtxt('/Users/dmelgar/Coquimbo2015/tsunami/survey.txt',usecols=1)
survey_eta=genfromtxt('/Users/dmelgar/Coquimbo2015/tsunami/survey.txt',usecols=2) 
        
# Get tide gauge observations
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/valp.sac')
gaugemax=st[0].data.max()
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/quin.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/pich.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/coqu.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
gaugelat=array([-33.0273083,-32.7755556,-32.1356110,-29.9500660])


dlat=0.1
#Make plot
plt.figure(figsize=(2, 6.6))
i=where((lat>-33.5) & (lat<-29))[0]

#tgauges=plt.scatter(gaugemax,gaugelat+dlat,marker='*',s=160,c='#1E90FF',zorder=4,label='Tide gauges')
#tgauges_obs=plt.barh(gaugelat+dlat,gaugemax[0],height=0.4,color='#1E90FF',zorder=2,label='TG obs.')
survey=plt.barh(survey_lat,survey_eta,height=0.03,color='#FF6347',label='Survey',lw=0)
tgauges_pre=plt.barh(gaugelat_pred,gaugemax_pred,height=0.06,color='#008080',label='TG pred.',lw=0)
static=plt.plot(eta[i],lat[i],lw=1.0,c='#606060',label='Predicted')

#plt.legend()
#plt.legend(['Predicted','Survey','Tide Gauges'],loc='upper center', bbox_to_anchor=(0.5, 1.27),
#    ncol=1, fancybox=True, shadow=True)
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),ncol=1, fancybox=True, shadow=True)

plt.ylim([-33.5,-29])
#plt.ylim([lat.min(),lat.max()])
plt.xlim([0,11])
#plt.xlim([0,etasurvey.max()+0.1*etasurvey.max()])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off')        # ticks along the top edge are off
plt.show()

