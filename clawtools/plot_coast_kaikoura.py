# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_,r_,argmin,squeeze,isnan
from obspy import read
from matplotlib import rcParams

etaclip=50
rcParams.update({'font.size': 22})


fgmax_file=u'/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/fort.FG1.aux1'

wet_tol=0.001
etaclip=4
minamr=3

#Get maximum amplitude first
lon=genfromtxt(fgmax_file,usecols=0)
lat=genfromtxt(fgmax_file,usecols=1)
amr=genfromtxt(fgmax_file,usecols=2)
H_in=genfromtxt(fgmax_file,usecols=3)
b_in=genfromtxt(aux_file)
unique_amr=unique(amr)
#i=where(unique_amr>0)[0]
#unique_amr=unique_amr[i]
eta=zeros(len(H_in))
H=zeros(len(H_in))
b=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    if unique_amr[k]<minamr:
        eta[i]=0
    else:
        eta[i]=H_in[i]+b_in[i,int(unique_amr[k]+1)]
        H[i]=H_in[i]
        b[i]=b_in[i,int(unique_amr[k]+1)]

#i=where(b<0)[0]
#lon=lon[i]
#lat=lat[i]
#eta=eta[i]
#H=H[i]
#b=b[i]
#
#i=where(H<10)[0]
#lon=lon[i]
#lat=lat[i]
#eta=eta[i]
#H=H[i]
#b=b[i]
#remove onshore points



#Wphase and slip inversion clipping
#iclip=where(lat>-25)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#i=where(isnan(eta)==False)[0]
#eta=eta[i]
#lat=lat[i]

#PGD clipping
#etaclip=0.6
#iclip=where(lat>-28.6)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#i=where(isnan(eta)==False)[0]
#eta=eta[i]
#lat=lat[i]
#iclip=where(lat<-34.7)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#i=where(isnan(eta)==False)[0]
#eta=eta[i]
#lat=lat[i]

