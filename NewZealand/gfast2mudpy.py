from numpy import genfromtxt,savetxt,arange,c_,ones


#fin=genfromtxt('/Users/dmelgar/NewZealand2016/ff/gfast_kaikoura_slipmodel.txt')
#fout='/Users/dmelgar/NewZealand2016/ff/gfast__kaikoura.fault'
#rupt='/Users/dmelgar/NewZealand2016/ff/gfast_kaikoura.rupt'
#strike=212.93*ones(len(fin))
#dip=25.27*ones(len(fin))
#L=30.59e3*ones(len(fin))
#W=16.73e3*ones(len(fin))


fin=genfromtxt('/Users/dmelgar/NewZealand2016/ff/nz_hamling_40km_slipmodels.txt')
fout='/Users/dmelgar/NewZealand2016/ff/gfast_kaikoura_hamling.fault'
rupt='/Users/dmelgar/NewZealand2016/ff/gfast_kaikoura_hamling.rupt'
strike=236.77*ones(len(fin))
dip=69.56*ones(len(fin))
L=30.59e3*ones(len(fin))
W=16.73e3*ones(len(fin))



f=c_[arange(1,len(fin)+1),fin[:,1],fin[:,2],fin[:,3],strike,dip,0.5*ones(len(fin)),5.0*ones(len(fin)),L,W]
r=c_[arange(1,len(fin)+1),fin[:,1],fin[:,2],fin[:,3],strike,dip,0.5*ones(len(fin)),5.0*ones(len(fin)),fin[:,4],fin[:,5,],L,W,1.0*ones(len(fin)),30e9*ones(len(fin))]

savetxt(fout,f,fmt='%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f')
savetxt(rupt,r,fmt='%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2e')