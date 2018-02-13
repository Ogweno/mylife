from numpy import genfromtxt,zeros

sta =genfromtxt('/Users/dmelgar/Slip_inv/Nepal_fwd/data/station_info/nepal.sta',usecols=0,dtype='S')
ll =genfromtxt('/Users/dmelgar/Slip_inv/Nepal_fwd/data/station_info/nepal.sta',usecols=[1,2])
f=open('/Users/dmelgar/Slip_inv/Nepal_fwd/output/forward_models/billham_offsets.txt','w')
f.write('lon,lat,N(m),E(m),U(m)\n')
for k in range(len(sta)):
    neu=genfromtxt(u'/Users/dmelgar/Slip_inv/Nepal_fwd/output/forward_models/'+sta[k]+'.static.neu')
    n=neu[0]
    e=neu[1]
    u=neu[2]
    out='%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n' %(sta[k],ll[k,0],ll[k,1],n,e,u)
    f.write(out)
f.close()
