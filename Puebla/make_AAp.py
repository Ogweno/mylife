from numpy import genfromtxt,linspace,savetxt,c_


f=genfromtxt('/Users/dmelgar/Puebla2017/cross-section/AAp.txt')

x=linspace(f[0,0],f[1,0],100)
y=linspace(f[0,1],f[1,1],100)
out=c_[x,y]
savetxt('/Users/dmelgar/Puebla2017/cross-section/AAp.xy',out)