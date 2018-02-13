from numpy import savetxt,genfromtxt,arange,zeros,argmin,c_

g=genfromtxt('/Users/dmelgar/Cascadia_M9/glarms/glarms_P.txt')

t=arange(0,300.5,0.5)
mag=zeros(len(t))
for k in range(len(t)):
    i=argmin(abs(t[k]-g[:,0]))
    mag[k]=g[i,2]
    
savetxt(u'/Users/dmelgar/code/GMT/Cascadia_M9/glarms.mag',c_[t,mag],fmt='%.2f')