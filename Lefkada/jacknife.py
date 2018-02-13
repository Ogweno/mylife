import numpy as np
from obspy.core import UTCDateTime
from mudpy import inverse as inv
from numpy import zeros,dot,array,squeeze,expand_dims,empty,tile,eye,argsort,sort
from scipy.sparse import csr_matrix as sparse
from scipy.optimize import nnls
from numpy.random import rand

fraction_keep=0.8
iter=50
reg_spatial=0.17378008287493754
#reg_spatial=0.143844989
reg_temporal=None

home='/Users/dmelgar/Slip_inv/'
project_name='Lefkada65'
run_name='jacknife'

model_name='haslinger99.mod'   #Velocity model
fault_name='lefkada65.fault'    #Fault geometry
GF_list='gps_sm_insar.gflist'#What GFs are to be computed for each station
tgf_file=None#'seafloor.sta'
G_name='gps_sm_insar_5win_vr2.6' #run_name#'review_3.2_low' #Either name of GF matrix to load or name to save GF matrix with
# Displacement and velocity wcloseaveform parameters
NFFT=512; dt=0.2
#Tsunami deformation parameters
tsunNFFT=128 ; tsun_dt=2.0
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
custom_stf=None
################################################################################

#############               Inversion Parameters               ################# 
time_epi=UTCDateTime('2015-11-17T07:10:07')
epicenter=np.array([20.6002,38.6655,10.7]) #original depth is 25, 29 is to match subfault view
rupture_speed=2.6 #Fastest rupture allowed in km/s
num_windows=5
nstrike=40 ; ndip=12 ; nfaults=(nstrike,ndip) #set nstrike to total no. of faults and ndip to 1 if using Tikh
beta=135 #Rotational offset (in degrees) applied to rake (0 for normal)
Ltype=0 # 0 for Tikhonov and 2 for Laplacian
solver='nnls' # 'lstsq','nnls'
top='free' ; bottom='locked' ; left='locked' ; right='locked' #'locked' or 'free'
bounds=(top,bottom,left,right)
################################################################################e=
########      Run-time modifications to the time series             ############
weight=True
decimate=None  #Decimate by constant (=None for NO decimation)
bandpass=None #Corner frequencies in Hz =None if no filter is desired
################################################################################





dcomplete=inv.getdata(home,project_name,GF_list,decimate,bandpass=None)
#Get GFs
Gcomplete=inv.getG(home,project_name,fault_name,model_name,GF_list,True,G_name,epicenter,
            rupture_speed,num_windows,decimate,bandpass)

print 'Applying data weights'
w=inv.get_data_weights(home,project_name,GF_list,dcomplete,decimate)
W=empty(Gcomplete.shape)
W=tile(w,(Gcomplete.shape[1],1)).T
WG=empty(Gcomplete.shape)
WG=W*Gcomplete
wd=w*dcomplete.squeeze()
wd=expand_dims(wd,axis=1)
WGcomplete=WG.copy()
wdcomplete=wd.copy()
#Clear up extraneous variables
W=None
w=None       
            
for k in range(iter):
    
    kout=k
    print k
    #Throw the dice
    R=rand(len(dcomplete))
    keep=argsort(R)[0:int(len(dcomplete)*fraction_keep)]
    keep=sort(keep)
    d=dcomplete[keep,0]
    G=Gcomplete[keep,:]
    wd=wdcomplete[keep,0]
    WG=WGcomplete[keep,:]
    

    #Define inversion quantities
    x=WG.transpose().dot(wd)
    print 'Computd=ing G\'G'
    K=(WG.T).dot(WG)
    
    #Regularization
    static=False #Is it jsut a static inversion?
    if reg_spatial!=None:
        if Ltype==2: #Laplacian smoothing
            Ls=inv.getLs(home,project_name,fault_name,nfaults,num_windows,bounds)
        else: #Tikhonov smoothing
            N=nfaults[0]*nfaults[1]*num_windows*2 #Get total no. of model parameters
            Ls=eye(N) 
        Ninversion=1
    else:
        Ls=zeros(K.shape)
        reg_spatial=array([0.])
        Ninversion=1
    if reg_temporal!=None:
        Lt=inv.getLt(home,project_name,fault_name,num_windows)
        Ninversion=1
    else:
        Lt=zeros(K.shape)
        reg_temporal=0
        static=True
    #Make L's sparse
    Ls=sparse(Ls)
    Lt=sparse(Lt)
    #Get regularization tranposes for ABIC
    LsLs=Ls.transpose().dot(Ls)
    LtLt=Lt.transpose().dot(Lt)
    #inflate
    Ls=Ls.todense()
    Lt=Lt.todense()
    LsLs=LsLs.todense()
    LtLt=LtLt.todense()
    
    
    
    #Run inversion
    Kinv=K+(reg_spatial**2)*LsLs+(reg_temporal**2)*LtLt
    x=squeeze(x.T)
    sol,res=nnls(Kinv,x)
    x=expand_dims(x,axis=1)
    sol=expand_dims(sol,axis=1)

    #Get moment
    Mo,Mw=inv.get_moment(home,project_name,fault_name,model_name,sol)
    print Mw
    #If a rotational offset was applied then reverse it for output to file
    sol=inv.rot2ds(sol,beta)
    #Write log
    #Write output to file
    inv.write_model(home,project_name,'jacknife',fault_name,model_name,rupture_speed,num_windows,epicenter,sol,kout)
    
    
