from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


fgmax_file=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static_3_survey/_output/fort.FG2.valuemax'
aux_file='/Users/dmelgar/Tsunamis/maule_gps_tg_static_3_survey/_output/fort.FG2.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/maule_gps_static_3_survey/_output/fort.FG2.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/maule_gps_static_3_survey/_output/fort.FG2.aux1'


survey_file='/Users/dmelgar/Maule2010/tsunami/survey.txt'
latsurvey=genfromtxt(survey_file,usecols=1)
etasurvey=genfromtxt(survey_file,usecols=3)

#Get maximum amplitude first
lat=genfromtxt(fgmax_file,usecols=1)
amr=genfromtxt(fgmax_file,usecols=2)
Hi=genfromtxt(fgmax_file,usecols=3)
bi=genfromtxt(aux_file)
unique_amr=unique(amr)
eta=zeros(len(Hi))
H=zeros(len(Hi))
b=zeros(len(Hi))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    eta[i]=Hi[i]+bi[i,unique_amr[k]+1]
    H[i]=Hi[i]
    b[i]=bi[i,unique_amr[k]+1]
i=where(H==0)[0]
eta[i]=nan
H[i]=nan
b[i]=nan

plt.figure(figsize=(5.5,10))
plt.scatter(etasurvey,latsurvey,c='#66B2FF',label='Survey',s=40,lw=0)
plt.scatter(eta,lat,c='#FF4500',label='Modeled',s=40,lw=0)
plt.xlabel('Run-up (m)')
plt.ylim([-39.5,-32.5])
plt.xlim([0,35])
plt.legend()
plt.show()