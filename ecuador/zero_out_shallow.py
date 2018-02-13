import matplotlib.path as path
from numpy import genfromtxt,where,savetxt

p=genfromtxt(u'/Users/dmelgar/Ecuador2016/tsunami/zero_outside.txt')
fault=genfromtxt('/Users/dmelgar/Slip_inv/Ecuador_insar/output/inverse_models/models/sm_insar_5win_vr3.4.0014.inv.total')
fout='/Users/dmelgar/Slip_inv/Ecuador_insar/output/inverse_models/models/sm_insar_5win_vr3.4.0014.noshallow.inv.total'


zero_path=path.Path(p)

i=where(zero_path.contains_points(fault[:,1:3])==False)[0]
fault[i,8]=0
fault[i,9]=0

#Write
savetxt(fout,fault,fmt='%d\t%10.4f\t%10.4f\t%8.4f\t%8.2f\t%6.2f\t%6.2f\t%6.2f\t%12.4e\t%12.4e\t%8.2f\t%8.2f\t%8.4e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
