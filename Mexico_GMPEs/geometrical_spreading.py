from matplotlib import pyplot as plt
from numpy import arange,where,zeros,sqrt,logspace


d=arange(1,400,1)

#coast
Gcoast=zeros(len(d))
Gcoast=1./d
i=where(d>50)[0]
Gcoast[i]=1./sqrt(50*d[i])


#land
Gland=zeros(len(d))
Gland=1./d
i=where((d>50) & (d<150))[0]
Gland[i]=1./50
i=where(d>=150)[0]
Gland[i]=sqrt(3)/sqrt(50*d[i])


plt.figure()
plt.loglog(d,Gcoast)
plt.loglog(d,Gland)
plt.xlabel('d (km)')
plt.ylabel('Geom. spreading')
plt.legend(['coast','land'])
plt.grid()
plt.show()


f=logspace(-2,2)
Qland=211*f**0.46
Qcoast=175*f**0.52

plt.figure()
plt.semilogx(f,Qcoast)
plt.semilogx(f,Qland)
plt.xlabel('frequency (Hz)')
plt.ylabel('Q')
plt.legend(['coast','land'])
plt.grid()
plt.show()