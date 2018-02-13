from numpy import genfromtxt,arange,r_

nrows=52
ncols=420
strike=356.0
dip=14.0


fault=genfromtxt(u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.000006.rupt')
fout=u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.lowpass.000006.rupt'




#Find the root idnices for the averaging (the top left corner)
N=nrows*ncols/16
i=arange(0,420,4)
for k in range(nrows/4):
    if k==0:
        root=i.copy()
    else:
        root=r_[root,i+(k)*ncols*4]



#Select 
adjust=r_[arange(4),arange(4)+ncols,arange(4)+2*ncols,arange(4)+3*ncols]
f=open(fout,'w')
f.write('#\n')
for k in range(len(root)):
    j=adjust+root[k]
    #Get moment of fault patch
    M0=sum(fault[j,9]*fault[j,10]*fault[j,11]*fault[j,13])
    #Get area of current fault aptch
    area=sum(fault[j,10]*fault[j,11])
    #Get mean rigidity of fault patch
    mu=fault[j,13].mean()
    #get slip of new fault path
    slip_out=M0/(mu*area)
    #Get coordinates of new patch
    lon_out=fault[j,1].mean()
    lat_out=fault[j,2].mean()
    z_out=fault[j,3].mean()
    rise_time=fault[j,7].mean()
    onset_time=fault[j,12].mean()
    
    line='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e\n' %(k+1,lon_out,lat_out,z_out,strike,dip,0.5,rise_time,0.00,slip_out,area**0.5,area**0.5,onset_time,mu)
    f.write(line)
    
f.close()