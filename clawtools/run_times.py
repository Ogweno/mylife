from matplotlib import pyplot as plt
from numpy import array

cores=array([1,2,3,4,5,6,9,12,15])
time_iq=array([145.9,73.9,51.7,39.6,32.9,28.5,23.1,22.8,22.5])/60
time_ma=array([483.2,248.2,169.6,129.7,106.5,89.8,68.3,64.1,60.1])/60
time_to=array([944.7,485.8,335.6,255.4,208.2,176.8,133.9,124.9,115.27])/60
time_co=array([399,205.4,141.5,109.5,89.3,75.92,59.5,55.7,53.3])/60

Kma=483.2/60
tma=Kma/cores

Kiq=145.9/60
tiq=Kiq/cores

Kto=944.67/60
tto=Kto/cores

Kco=399./60
tco=Kco/cores


plt.figure(figsize=(10,3))
plt.subplot(221)
plt.scatter(cores,time_ma,c='r',s=30)
plt.grid()
plt.ylabel('Avg. wall time (min)',fontsize=14)
plt.annotate('2010 Maule',xy=(7,7))
plt.plot(cores,tma,c='k')
plt.ylim([0,8.5])


plt.subplot(222)
plt.scatter(cores,time_to,c='r',s=30)
plt.grid()
plt.plot(cores,tto,c='k')
plt.ylim([0,16])
plt.annotate('2011 Tohoku-oki',xy=(3.8,13.5))


plt.subplot(223)
plt.scatter(cores,time_iq,c='r',s=30)
plt.grid()
plt.xlabel('No. CPUs',fontsize=14)
plt.plot(cores,tiq,c='k')
plt.ylim([0,2.5])
plt.annotate('2014 Iquique',xy=(6,2.1))
plt.ylabel('Avg. wall time (min)',fontsize=14)

plt.subplot(224)
plt.scatter(cores,time_co,c='r',s=30)
plt.grid()
plt.xlabel('No. CPUs',fontsize=14)
plt.plot(cores,tco,c='k')
plt.ylim([0,7])
plt.annotate('2015 Illapel',xy=(6,6))


plt.subplots_adjust(bottom=0.2)

plt.show()