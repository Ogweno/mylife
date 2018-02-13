from numpy import genfromtxt,arange,cos,sin,deg2rad,ones,zeros,c_,savetxt


f=genfromtxt('/Users/dmelgar/SStsunami/Wharton/Hill_et_al_2016_model6c.txt')

# No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)
no=arange(1,len(f)+1)
lon=f[:,0]
lat=f[:,1]
z=f[:,2]
strike=f[:,6]
dip=f[:,7]
rake=f[:,5]
slip=f[:,8]
length=f[:,3]*1000
width=f[:,4]*1000
ss_slip=slip*cos(deg2rad(rake))
ds_slip=slip*sin(deg2rad(rake))
dura=20*ones(len(f))
rise=0.5*ones(len(f))
rupt_time=zeros(len(f))
rigidity=60e9*ones(len(f))

out=c_[no,lon,lat,z,strike,dip,rise,dura,ss_slip,ds_slip,length,width,rupt_time,rigidity]

savetxt('/Users/dmelgar/SStsunami/Wharton/model6c.rupt',out,fmt='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e')