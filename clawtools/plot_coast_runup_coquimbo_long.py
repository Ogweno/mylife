# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,r_,argmin,squeeze,isnan
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


#fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_3_gps_sm_tg/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/coquimbo_3_gps_sm_tg/_output/fort.FG1.aux1'

#fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_30_wphase/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/coquimbo_30_wphase/_output/fort.FG1.aux1'

fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_detailed/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/coquimbo_detailed/_output/fort.FG1.aux1'

#fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_30_pgd/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/coquimbo_30_pgd/_output/fort.FG1.aux1'

fgmax_file2=u'/Users/dmelgar/Tsunamis/coquimbo_30_gps_sm_tg/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/coquimbo_30_gps_sm_tg/_output/fort.FG1.aux1'

gauge_data_coarse=u'/Users/dmelgar/Tsunamis/coquimbo_30_wphase/_output/fort.gauge'
#gauge_data_coarse=u'/Users/dmelgar/Tsunamis/coquimbo_30_pgd/_output/fort.gauge'
#gauge_data_coarse=u'/Users/dmelgar/Tsunamis/coquimbo_30_gps_sm_tg/_output/fort.gauge'
gauge_data_fine=u'/Users/dmelgar/Tsunamis/coquimbo_3_gps_sm_tg/_output/fort.gauge'

blend=True
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

if blend==True:
    #Get maximum amplitude first
    lon2=genfromtxt(fgmax_file2,usecols=0)
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
    latbounds=[-35,-27]
    i=where(lat2<latbounds[0])[0] #THings south of the fine model
    lat_out=lat2[i]
    lon_out=lon2[i]
    eta_out=eta2[i]
    #THings IN the fine model
    lat_out=r_[lat_out,lat]
    lon_out=r_[lon_out,lon]
    eta_out=r_[eta_out,eta]
    #Things north of the fine model
    i=where(lat2>latbounds[1])[0] #THings north of the fine model
    lat_out=r_[lat_out,lat2[i]]
    lon_out=r_[lon_out,lon2[i]]
    eta_out=r_[eta_out,eta2[i]]
    #Ok done
    lat=lat_out
    lon=lon_out
    eta=eta_out
else:
    eta_out=eta

#Wphase and slip inversion clipping
#iclip=where(lat>-25)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#i=where(isnan(eta)==False)[0]
#eta=eta[i]
#lat=lat[i]

#PGD clipping
#etaclip=0.6
#iclip=where(lat>-28.6)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#i=where(isnan(eta)==False)[0]
#eta=eta[i]
#lat=lat[i]
#iclip=where(lat<-34.7)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#i=where(isnan(eta)==False)[0]
#eta=eta[i]
#lat=lat[i]


#Get tide gauge predictions
stations=genfromtxt('/Users/dmelgar/Coquimbo2015/station_info/tide_gauge_geoclaw_srtm3.txt',usecols=3)
gaugelat_pred=genfromtxt('/Users/dmelgar/Coquimbo2015/station_info/tide_gauge_geoclaw_srtm3.txt',usecols=1)
if blend==True: #Get 1000-1006 from fine file
    gaugedata=genfromtxt(gauge_data_fine)
    for kgauge in range(0,6):
        print stations[kgauge]
        i=where(gaugedata[:,0]==stations[kgauge])[0]
        g=gaugedata[i,6]
        if kgauge==0:
            gaugemax_pred=g.max()
        else:
            gaugemax_pred=c_[gaugemax_pred,g.max()]
    gaugedata=genfromtxt(gauge_data_coarse)
    for kgauge in range(6,len(stations)):
        print stations[kgauge]
        i=where(gaugedata[:,0]==stations[kgauge])[0]
        j=where(gaugedata[i,1]==3)[0]
        g=gaugedata[i[j],6]
        if len(g)==0:
            g=array([0,0])
        gaugemax_pred=c_[gaugemax_pred,g.max()]
else:
    gaugedata=genfromtxt(gauge_data_coarse)
    for kgauge in range(0,len(stations)):
        print stations[kgauge]
        i=where(gaugedata[:,0]==stations[kgauge])[0]
        j=where(gaugedata[i,1]==3)[0]
        g=gaugedata[i[j],6]
        if len(g)==0:
            g=array([0,0])
        if kgauge==0:
            gaugemax_pred=g.max()
        else:
            gaugemax_pred=c_[gaugemax_pred,g.max()]
gaugemax_pred=squeeze(gaugemax_pred)
    
        
        


    
# Get tide gauge observations
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/anto.sac')
gaugemax=st[0].data.max()
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/bahm.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/buca.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/cald.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/chan.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/cons.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/coqu.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/coro.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/corr.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/huas.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/lebu.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/meji.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/papo.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/pich.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/queu.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/quin.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/quir.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/talt.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/toco.sac')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read(u'/Users/dmelgar/Coquimbo2015/tsunami/sac/valp.sac')
gaugemax=c_[gaugemax,st[0].data.max()]

gaugelat=genfromtxt('/Users/dmelgar/Coquimbo2015/station_info/tide_gauge_noislands.txt',usecols=1)

#Find maximum coastline comapred to gauges
gaugemax_pre=zeros(len(gaugemax[0]))
for k in range(len(gaugelat)):
    i=argmin(abs(gaugelat[k]-lat))
    gaugemax_pre[k]=eta[i]

dlat=0.25
#Make plot
plt.figure(figsize=(2, 6))
static=plt.plot(eta,lat,lw=0.5,c='k',zorder=3,label='Predicted')
tgauges=plt.scatter(gaugemax+dlat,gaugelat,marker='*',s=160,c='#1E90FF',zorder=4,label='Tide gauges')
#tgauges_obs=plt.barh(gaugelat+dlat,gaugemax[0],height=0.4,color='#1E90FF',zorder=2,label='TG obs.')
tgauges_pre=plt.barh(gaugelat_pred,gaugemax_pred,height=0.3,color='#909090',zorder=2,label='TG pred.')


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
xfill=array([3.,9.8])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#DC143C',alpha=1.0)

plt.ylim([-40,-21])
#plt.ylim([lat.min(),lat.max()])
plt.xlim([0,9.8])
#plt.xlim([0,etasurvey.max()+0.1*etasurvey.max()])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.show()

