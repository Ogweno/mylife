# -*- coding: utf-8 -*-
from numpy import genfromtxt,unique,zeros,where,nan,argmin,sqrt,array,r_,isnan
from glob import glob
    
#First read the simulation results
#

#fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_3_gps_sm_tg/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/coquimbo_3_gps_sm_tg/_output/fort.FG1.aux1'

#fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_30_wphase/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/coquimbo_30_wphase/_output/fort.FG1.aux1'

fgmax_file=u'/Users/dmelgar/Tsunamis/coquimbo_30_pgd/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/coquimbo_30_pgd/_output/fort.FG1.aux1'

fgmax_file2=u'/Users/dmelgar/Tsunamis/coquimbo_30_gps_sm_tg/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/coquimbo_30_gps_sm_tg/_output/fort.FG1.aux1'

fout='/Users/dmelgar/Coquimbo2015/tsunami/warnings/coquimbo_30_pgd.warn'

wet_tol=0.001 #Anything udner this is consdered dry
blend=False
etaclip=4


#Get maximum amplitude first
lon=genfromtxt(fgmax_file,usecols=0)
lat=genfromtxt(fgmax_file,usecols=1)
amr=genfromtxt(fgmax_file,usecols=2)
H_in=genfromtxt(fgmax_file,usecols=3)
b_in=genfromtxt(aux_file)
unique_amr=unique(amr)
eta=zeros(len(H_in))
b=zeros(len(H_in))
H=zeros(len(H_in))
for k in range(len(unique_amr)):
    i=where(amr==unique_amr[k])[0]
    eta[i]=H_in[i]+b_in[i,unique_amr[k]+1]
    H[i]=H_in[i]
    b[i]=b_in[i,unique_amr[k]+1]
i=where(H<wet_tol)[0]
eta[i]=nan

if blend==True:
    #Get maximum amplitude first
    lon2=genfromtxt(fgmax_file2,usecols=0)
    lat2=genfromtxt(fgmax_file2,usecols=1)
    amr=genfromtxt(fgmax_file2,usecols=2)
    H_in=genfromtxt(fgmax_file2,usecols=3)
    b_in=genfromtxt(aux_file2)
    unique_amr=unique(amr)
    eta2=zeros(len(H_in))
    b=zeros(len(H_in))
    H=zeros(len(H_in))
    for k in range(len(unique_amr)):
        i=where(amr==unique_amr[k])[0]
        eta2[i]=H_in[i]+b_in[i,unique_amr[k]+1]
        H[i]=H_in[i]
        b[i]=b_in[i,unique_amr[k]+1]
    i=where(H<wet_tol)[0]
    eta2[i]=nan
    #Now combine
    latbounds=[-35,-27]
    i=where(lat2<latbounds[0])[0] #THings south of the fine model
    lat_out=lat2[i]
    lon_out=lon2[i]
    eta_out=eta2[i]
    #THings IN the fine model
    lat_out=r_[lat_out,lat]
    lon_out=r_[lon_out,lon]
    eta_out=r_[eta_out,eta]
    #Things north of the fine model
    i=where(lat2>latbounds[1])[0] #THings north of the fine model
    lat_out=r_[lat_out,lat2[i]]
    lon_out=r_[lon_out,lon2[i]]
    eta_out=r_[eta_out,eta2[i]]
    #Ok done
    lat=lat_out
    lon=lon_out
    eta=eta_out
else:
    eta_out=eta
#Wphase clips
#iclip=where(lat>-25)[0]
#i=where(eta[iclip]>etaclip)[0]
#eta[iclip[i]]=nan
#PGD clipping
etaclip=0.4
iclip=where(lat>-28.6)[0]
i=where(eta[iclip]>etaclip)[0]
eta[iclip[i]]=nan
i=where(isnan(eta)==False)[0]
eta=eta[i]
lat=lat[i]
lon=lon[i]
iclip=where(lat<-34.7)[0]
i=where(eta[iclip]>etaclip)[0]
eta[iclip[i]]=nan
i=where(isnan(eta)==False)[0]
eta=eta[i]
lat=lat[i]
lon=lon[i]
iclip=where(lat>-23.1)[0]
eta[iclip]=0




#Define region class

class region:
  def __init__(self):
    self.lon = None
    self.lat = None
    self.eta = []
    self.color=None

def read_region(f):
    R=open(f,'r')
    lon=[]
    lat=[]
    for line in R:
        if '>' in line:
            pass
        else:
            lon.append(float(line.split()[0]))
            lat.append(float(line.split()[1]))
    lon=array(lon)
    lat=array(lat)
    return lon,lat

#Instantiate regions
anto=region()
arau=region()
atac=region()
ayse=region()
biob=region()
coqu=region()
losl=region()
maga=region()
maul=region()
ohig=region()
sant=region()
tara=region()
valp=region()

#Now read all regions

anto.lon,anto.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_antofagasta.txt')
arau.lon,arau.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_araucania.txt')
atac.lon,atac.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_atacama.txt')
ayse.lon,ayse.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_aysen.txt')
biob.lon,biob.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_biobio.txt')
coqu.lon,coqu.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_coquimbo.txt')
losl.lon,losl.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_loslagos.txt')
maga.lon,maga.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_magallanes.txt')
maul.lon,maul.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_maule.txt')
ohig.lon,ohig.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_ohiggins.txt')
sant.lon,sant.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_santiago.txt')
tara.lon,tara.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_tarapaca.txt')
valp.lon,valp.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_valparaiso.txt')

regions=[anto,arau,atac,ayse,biob,coqu,losl,maga,maul,ohig,sant,tara,valp]

#Loop over every coast point, see what region it belongs to
for k in range(len(lon)):
    d=[]
    if k==269:
        a=0
    print k
    for kregion in range(len(regions)):
        d.append(sqrt((lon[k]-regions[kregion].lon)**2+(lat[k]-regions[kregion].lat)**2).min())
    if argmin(d)==0:
        anto.eta.append(eta[k])
    if argmin(d)==1:
        arau.eta.append(eta[k])
    if argmin(d)==2:
        atac.eta.append(eta[k])
    if argmin(d)==3:
        ayse.eta.append(eta[k])
    if argmin(d)==4:
        biob.eta.append(eta[k])
    if argmin(d)==5:
        coqu.eta.append(eta[k])
    if argmin(d)==6:
        losl.eta.append(eta[k])
    if argmin(d)==7:
        maga.eta.append(eta[k])
    if argmin(d)==8:
        maul.eta.append(eta[k])
    if argmin(d)==9:
        ohig.eta.append(eta[k])
    if argmin(d)==10:
        sant.eta.append(eta[k])
    if argmin(d)==11:
        tara.eta.append(eta[k])
    if argmin(d)==12:
        valp.eta.append(eta[k])
        
#Assign wanring level
def warn_color(region):    
    eta=array(region.eta)
    i=where(isnan(eta)==True)[0]
    eta[i]=0
    red='220/20/60'
    orange='255/140/0'
    yellow='255/215/0'
    green='50/205/50'
    if len(eta) < 1:
        color=green
    elif max(eta) > 3:
        color=red
    elif max(eta) > 1:
        color=orange
    elif max(eta) > 0.3:
        color=yellow
    elif max(eta) <= 0.3:
        color=green
    return color

#Get warning levels
#[anto,arau,atac,ayse,biob,coqu,losl,maga,maul,ohig,sant,tara,valp]
anto.color=warn_color(anto)
arau.color=warn_color(arau)
atac.color=warn_color(atac)
ayse.color=warn_color(ayse)
biob.color=warn_color(biob)
coqu.color=warn_color(coqu)
losl.color=warn_color(losl)
maga.color=warn_color(maga)
maul.color=warn_color(maul)
ohig.color=warn_color(ohig)
sant.color=warn_color(sant)
tara.color=warn_color(tara)
valp.color=warn_color(valp)

#Write warning file
f=open(fout,'w')
f.write('Antofagasta='+anto.color+'\n')
f.write('Araucania='+arau.color+'\n')
f.write('Atacama='+atac.color+'\n')
f.write('Aysen='+ayse.color+'\n')
f.write('Bio Bio='+biob.color+'\n')
f.write('Coquimbo='+coqu.color+'\n')
f.write('Los Lagos='+losl.color+'\n')
f.write('Magallanes='+maga.color+'\n')
f.write('Maule='+maul.color+'\n')
f.write('O Higgins='+ohig.color+'\n')
f.write('Santiago=160/160/160\n')
f.write('Tarapaca='+tara.color+'\n')
f.write('Valparaiso='+valp.color)
f.close()