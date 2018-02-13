from numpy import genfromtxt,unique,zeros,where

f=genfromtxt(u'/Users/dmelgar/Slip_inv/Chiapas_hernandez_new/output/inverse_models/models/final_vr_3.4.0007.inv')

z=unique(f[:,3])

M0=zeros(len(z))
for k in range(len(z)):
    i=where(f[:,3]==z[k])[0]
    M0[k]=((f[i,8]**2+f[i,9]**2)**0.5*f[i,10]*f[i,11]*f[i,12]).sum()
    