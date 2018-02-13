from numpy import genfromtxt,ones,argmin,where,zeros


afters_thrust='/Volumes/Kanagawa/Slip_Inv/tohoku_10s/analysis/thrust_afters_select.txt'
afters_normal='/Volumes/Kanagawa/Slip_Inv/tohoku_10s/analysis/normal_afters_select.txt'
stress_file='/Users/dmelgar/code/GMT/tohoku/stress.xyz'

normal=genfromtxt(afters_normal)
thrust=genfromtxt(afters_thrust)
stress=genfromtxt(stress_file)

#Assign x variable (fault type)
x=ones(len(normal)+len(thrust))
x[0:len(normal)]=-1

#Now find if the stress drop is positive or negative at the lcoation of the aftershock
y=zeros(len(normal)+len(thrust))
for k in range(len(normal)):
    lon=normal[k,0]
    lat=normal[k,1]
    d=((stress[:,0]-lon)**2+(stress[:,1]-lat)**2)**0.5
    imin=argmin(d)
    if stress[imin,2]>0:
        y[k]=-1.
    else:
        y[k]=1.

i=0
for k in range(len(normal),len(normal)+len(thrust)):
    lon=thrust[i,0]
    lat=thrust[i,1]
    d=((stress[:,0]-lon)**2+(stress[:,1]-lat)**2)**0.5
    imin=argmin(d)
    if stress[imin,2]<0:
        y[k]=1.
    else:
        y[k]=-1.
    i+=1
    
#Compute phi
n11=len(where((x==1) & (y==1))[0])
n10=len(where((x==1) & (y==-1))[0])
n01=len(where((x==-1) & (y==1))[0])
n00=len(where((x==-1) & (y==-1))[0])

Nx1=n11+n01
Nx0=n10+n00
N1y=n11+n10
N0y=n01+n00

phi=(n11*n00-n10*n01)/((Nx1*Nx0*N1y*N0y)**0.5)
print phi