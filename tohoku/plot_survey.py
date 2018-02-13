from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


survey=genfromtxt('/Users/dmelgar/Tohoku2011/Tsunami/inundation/survey_all.txt')
latsurvey=survey[:,0]
etasurvey=survey[:,1]

survey=genfromtxt('/Users/dmelgar/Tohoku2011/Tsunami/inundation/survey_yesinundated.txt')
lat=survey[:,0]
eta=survey[:,1]


plt.figure(figsize=(5.5,10))
plt.scatter(etasurvey,latsurvey,c='#66B2FF',label='Survey',s=20,lw=0)
plt.scatter(eta,lat,c='#FF4500',label='Modeled',s=20,lw=0)
plt.xlabel('Run-up (m)')
plt.ylim([36,41])
plt.xlim([0,40])
#plt.legend()
plt.show()