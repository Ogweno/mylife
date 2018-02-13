from numpy import load,save,r_

fout='/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gsi_vr3.6.npy'

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gs_vr3.6.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/T62T83_4win.npy")

G=r_[G1,G2]

save(fout,G)