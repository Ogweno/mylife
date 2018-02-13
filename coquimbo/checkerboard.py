from numpy import load,empty,tile
from numpy import genfromtxt,arange,zeros,where,eye,logspace,squeeze,expand_dims,savetxt
from obspy import read
from scipy.optimize import nnls
from string import rjust
from mudpy import gmttools

run='joint'
G=load("/Users/dmelgar/Slip_inv/Coquimbo_4s/GFs/matrices/gps_sm_tg_insar_1win.npy")
f=genfromtxt(u'/Users/dmelgar/Slip_inv/Coquimbo_4s/output/forward_models/checker.rupt')
w=load("/Users/dmelgar/Slip_inv/Coquimbo_4s/GFs/matrices/data_weights.npy")
reg_spatial=logspace(-3,0,20)

nfaults=len(f)
iss=arange(0,nfaults*2,2)
ids=arange(1,nfaults*2,2)
m=zeros(nfaults*2)
m[iss]=f[:,8]
m[ids]=f[:,9]
#Make data vector
d=G.dot(m)

#Figure out what indices are which data
gf_file=u'/Users/dmelgar/Slip_inv/Coquimbo_4s/data/station_info/gps_sm_tg_insar.gflist'
GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
istatic=where(GF[:,0]==1)[0]
idisp=where(GF[:,1]==1)[0]
ivel=where(GF[:,2]==1)[0]
itsun=where(GF[:,3]==1)[0]
iinsar=where(GF[:,4]==1)[0]

GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],dtype='S') 
kstart=0
kend=0
#get indices
if len(istatic)>0:
    kend+=len(istatic)*3
    jstatic=arange(0,kend)
if len(idisp)>0: #read displacememnt waveforms
    kstart=kend
    for kfile in range(len(idisp)):
        n=read(GFfiles[idisp[kfile],1]+'.n')
        kend+=3*n[0].stats.npts
    jdisp=arange(kstart,kend)
if len(ivel)>0:
    kstart=kend
    for kfile in range(len(ivel)):
        n=read(GFfiles[ivel[kfile],2]+'.n')
        kend+=3*n[0].stats.npts
    jvel=arange(kstart,kend)
if len(itsun)>0:
    kstart=kend
    for kfile in range(len(itsun)):
        tsun=read(GFfiles[itsun[kfile],3])
        kend+=tsun[0].stats.npts
    jtsun=arange(kstart,kend)
if len(iinsar)>0:
    kstart=kend
    kend+=len(iinsar)
    jinsar=arange(kstart,kend)
    
#Decide which data
#i=jinsar
i=arange(0,kend)
d=d[i]
G=G[i,:]

#Data weights an dinversion if joint
W=empty(G.shape)
W=tile(w,(G.shape[1],1)).T
WG=empty(G.shape)
WG=W*G
wd=w*d.squeeze()
wd=expand_dims(wd,axis=1)
#Clear up extraneous variables
W=None
w=None
#Define inversion quantities
x=WG.transpose().dot(wd)
print 'Computing G\'G'
K=(WG.T).dot(WG)


##Inversion quantities if single data set
#K=(G.T).dot(G)
#x=G.transpose().dot(d)

#Get smoothing matrix
L=eye(len(m)) 
LL=L.transpose().dot(L)
for ks in range(len(reg_spatial)):
    print ks
    lambda_spatial=reg_spatial[ks]
    Kinv=K+(lambda_spatial**2)*LL
    x=squeeze(x.T)
    try:
        sol,res=nnls(Kinv,x)
    except:
        print '+++ WARNING: No solution found, writting zeros.'
        sol=zeros(G.shape[1])
        x=expand_dims(x,axis=1)
        sol=expand_dims(sol,axis=1)
    #Save
    number=rjust(str(ks),4,'0')
    fout='/Users/dmelgar/Slip_inv/Coquimbo_4s/output/forward_models/'+run+'.'+number+'.rupt'
    f[:,8]=sol[iss]
    f[:,9]=sol[ids]
    savetxt(fout,f,fmt='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e')
    #gmttools.make_total_model(fout,0)
