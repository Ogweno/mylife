from matplotlib import pyplot as plt

vr=[2.6,2.8,3.0,3.1,3.2,3.3,3.4,3.6]
VR=[48.37,57.55,64.92,69.12,72.82,74.12,73.1,67.66]


plt.plot(vr,VR,'k',lw=1.5)
plt.scatter(vr,VR,marker='+',s=90,lw=1.5)
plt.grid()
plt.xlabel('Rupture speed (km/s)')
plt.ylabel('Variance reductio (%)')
plt.show()