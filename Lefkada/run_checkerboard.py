from numpy import genfromtxt,arange,zeros,expand_dims,load,eye,squeeze,savetxt
from mudpy.inverse import ds2rot,rot2ds
from scipy.optimize import nnls



home='/Users/dmelgar/Slip_inv/'
project_name='Lefkada65'
fault_name='lefkada65.fault' 
nfaults=(40,12)
num_windows=1



#Load model and roatet
#f=genfromtxt(u'/Users/dmelgar/Slip_inv/Lefkada_fwd/forward_models/checkerboard.rupt')
f=genfromtxt(u'/Users/dmelgar/Slip_inv/Lefkada_fwd/forward_models/checkerboard_1win.rupt')

iss=arange(0,2*len(f),2)
ids=arange(1,2*len(f),2)
m1=zeros((2*len(f),1))
m1[iss]=expand_dims(f[:,8],1)
m1[ids]=expand_dims(f[:,9],1)
m=ds2rot(m1,135)

#GFs
G=load("/Users/dmelgar/Slip_inv/Lefkada_fwd/GFs/matrices/gps_sm_insar_5win_vr2.6.npy")
G=G[:,0:960]
dsynth=G.dot(m)


#Smoothing
N=nfaults[0]*nfaults[1]*num_windows*2 #Get total no. of model parameters
Ls=eye(N) 
LsLs=Ls.transpose().dot(Ls)


#WHich data are we doing?

##GPS
#iselect=arange(0,459)
#weight=1.0
#lambda_spatial=0.0001
#fout=u'/Users/dmelgar/Slip_inv/Lefkada_fwd/forward_models/checkerboard_gps_ouput.rupt'

##SM
iselect=arange(459,2871)
weight=1.0
lambda_spatial=0.05
#lambda_spatial=0.0001
fout=u'/Users/dmelgar/Slip_inv/Lefkada_fwd/forward_models/checkerboard_sm_ouput.rupt'

##Insar
#iselect=arange(2871,len(dsynth))
#weight=1.0
#lambda_spatial=0.00001
#fout=u'/Users/dmelgar/Slip_inv/Lefkada_fwd/forward_models/checkerboard_insar_ouput.rupt'

##All
#iselect=arange(len(dsynth))
#weight=1.0
#lambda_spatial=0.03
#fout=u'/Users/dmelgar/Slip_inv/Lefkada_fwd/forward_models/checkerboard_all_ouput.rupt'

#select subset of data
dsynth=dsynth[iselect]
G=G[iselect,:]



#prepare matrices
W=eye(len(dsynth))/weight
WG=W.dot(G)
wd=W.dot(dsynth)
x=WG.transpose().dot(wd)
K=(WG.T).dot(WG)
Kinv=K+(lambda_spatial**2)*LsLs

#Invert
x=squeeze(x.T)
sol,res=nnls(Kinv,x)

#Output
mout=rot2ds(expand_dims(sol,1),135)
f[:,8]=mout[iss,0]
f[:,9]=mout[ids,0]
fmt='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
savetxt(fout,f,fmt)







