from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,arange,argmin
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


fgmax_file=u'/Users/dmelgar/Tsunamis/iquique_42_3_long/_output/fort.FG3.valuemax'
aux_file='/Users/dmelgar/Tsunamis/iquique_42_3_long/_output/fort.FG3.aux1'



survey_file='/Users/dmelgar/Iquique2014/tsunami/iquique_survey.txt'
latsurvey=genfromtxt(survey_file,usecols=5)
etasurvey=genfromtxt(survey_file,usecols=7)

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
eta[i]=9999
H[i]=nan
b[i]=nan
eta_out=zeros(30)
lat_out=zeros(30)
for k in range(30):
    temp=[eta[k],eta[30+k],eta[60+k],eta[90+k]]
    lat_temp=[lat[k],lat[30+k],lat[60+k],lat[90+k]]
    itemp=argmin(temp)
    lat_out[k]=lat_temp[itemp]
    eta_out[k]=temp[itemp]
    

plt.figure(figsize=(5.5,10))
plt.scatter(etasurvey,latsurvey,c='#66B2FF',label='Survey',s=40,lw=0)
plt.scatter(eta_out,lat_out,c='#FF4500',label='Modeled',s=40,lw=0)
plt.xlabel('Run-up (m)')
#plt.ylim([-39.5,-32.5])
plt.xlim([0,5])
plt.legend()
plt.show()