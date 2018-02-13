from numpy import genfromtxt,sqrt,zeros,argmin,c_,savetxt

rupt=genfromtxt('/Users/dmelgar/Tohoku2011/Minsons/rupt.txt')
f=genfromtxt('/Users/dmelgar/Tohoku2011/Minsons/minson_3km.fault')
fout='/Users/dmelgar/Tohoku2011/Minsons/minson_3km.rupt'
S=sqrt(rupt[:,8]**2+rupt[:,9]**2)

slip=zeros((len(f),1))
for k in range(len(f)):
    if k%100==0:
        print k
    dist=sqrt((f[k,1]-rupt[:,2])**2+(f[k,2]-rupt[:,1])**2)
    slip[k,0]=S[argmin(dist)]
    
# No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)
#1	-125.237101	 46.418952	 10.8050	 353.00	  10.00	 0.5	 0.00	 0.00	 0.00	  16926.90	  16926.90	 0.00	3.577392e+10
#f[:,4]=f[:,4]+180
out=c_[f[:,0:8],zeros((len(f),1)),slip,f[:,-2:],zeros((len(f),2))]
savetxt(fout,out,fmt='%10d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t')

#savetxt('/Users/dmelgar/Tohoku2011/Minsons/minson_3km.fault',f,fmt='%10d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f')