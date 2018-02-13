from numpy import load,save,r_



G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr2.2.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr2.2.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr2.2.npy",G)

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr2.4.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr2.4.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr2.4.npy",G)

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr2.6.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr2.6.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr2.6.npy",G)

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr2.8.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr2.8.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr2.8.npy",G)

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr3.0.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr3.0.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr3.0.npy",G)

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr3.2.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr3.2.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr3.2.npy",G)

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr3.4.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr3.4.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr3.4.npy",G)

G1=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_vr3.6.npy")
G2=load("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/sm_vr3.6.npy")
G=r_[G1,G2]
save("/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/matrices/gps_sm_vr3.6.npy",G)