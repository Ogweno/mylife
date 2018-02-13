# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,argsort,isnan,argmin
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})

#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_gps_tg_30/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_gps_tg_30/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_gps_30/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_gps_30/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_pgd/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_pgd/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_wphase/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_wphase/_output/fort.FG1.aux1'
fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_glarms/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/tohoku_glarms/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_gps_tg_3/_output/fort.FG2.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_gps_tg_3/_output/fort.FG2.aux1'

fgmax_file2=u'/Users/dmelgar/Tsunamis/iquique_42_30_long/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/iquique_42_30_long/_output/fort.FG1.aux1'

survey_file='/Users/dmelgar/Iquique2014/tsunami/iquique_survey.txt'
blend=False

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

if blend==True:
    i=where(isnan(eta)==True)[0]
    #Now go one at a time and find closest one in other survey
    for k in range(len(i)):
        ifill=argmin(abs(lat[i[k]]-lat2))
        if eta2[ifill]!=nan:
            eta[i[k]]=eta2[ifill]
        else:
            eta[i[k]]=eta2[ifill+1]

        



#Parse survey data
latsurvey=genfromtxt(survey_file,usecols=5)
etasurvey=genfromtxt(survey_file,usecols=7)
    
#Get maximum gauge amplitudes
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/AKA.tsun')
gaugemax=st[0].data.max()
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/CHI.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/CHO.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/HAK.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/HAN.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/IRO.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/KUS.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/MAI.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/MER.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/MYE.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/OKA.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/ONA.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/OWA.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/SHM.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/SHR.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/TOB.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
st=read('/Users/dmelgar/Tohoku2011/Tsunami/gauges/URA.tsun')
gaugemax=c_[gaugemax,st[0].data.max()]
gaugelat=genfromtxt('/Users/dmelgar/Tohoku2011/Tsunami/gauges/gauges.sta',usecols=2)

gaugemax_pre=zeros(len(gaugemax[0]))
for k in range(len(gaugelat)):
    i=argmin(abs(gaugelat[k]-lat))
    print k
    print eta[i]
    if isnan(eta[i])==False:
        gaugemax_pre[k]=eta[i]
    else:
        gaugemax_pre[k]=eta[i+1]        




#Sort things
i=argsort(lat)
lat=lat[i]
eta=eta[i]

dlat=0.075
#Make plot
i=where(isnan(eta)==False)[0]
plt.figure(figsize=(2,4.5))
#static=plt.plot(eta2,lat,lw=1.6,c='#6666FF',zorder=3,label='Static')
kinematic=plt.plot(eta[i],lat[i],lw=0.5,c='k',zorder=4,label='Predicted')
#survey=plt.scatter(etasurvey,latsurvey,marker='o',s=40,c='#7FFF00',zorder=1,label='Survey')
#tgauges_obs=plt.barh(gaugelat+dlat,gaugemax[0],height=0.17,color='#1E90FF',zorder=2,label='TG obs.')
tgauges_pre=plt.barh(gaugelat-dlat,gaugemax_pre,height=0.2,color='#909090',zorder=2,label='TG pred.')
#plt.legend(['Instantaneous','Kinematic','Survey','Tide Gauges'])
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.23),
#    ncol=1, fancybox=True, shadow=True)
xfill=array([0,0.3])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#32CD32')
xfill=array([0.3,1.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FFD700')
xfill=array([1.,3.])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#FF8C00')
xfill=array([3.,99])
plt.fill_between(xfill,lat.min(),lat.max(),facecolor='#DC143C')
#plt.scatter(etasurvey,latsurvey,marker='o',s=40,c='#7FFF00',zorder=5)
plt.scatter(gaugemax,gaugelat,marker='*',s=160,c='#1E90FF',zorder=6)
plt.ylim([lat.min(),lat.max()])
plt.xlim([0,30])
#plt.xlim([0,etasurvey.max()+0.35*etasurvey.max()])
plt.xlabel('Maximum expected amplitude (m)')
plt.ylabel('Latitude ')
plt.subplots_adjust(left=0.22, right=0.95, bottom=0.1, top=0.8)
plt.show()

