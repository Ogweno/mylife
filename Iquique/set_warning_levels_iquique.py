from numpy import genfromtxt,unique,zeros,where,nan,argmin,sqrt,array,r_,isnan
from glob import glob
    
#First read the simulation results

#fgmax_file=u'/Users/dmelgar/Tsunamis/iquique_42_3_long/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/iquique_42_3_long/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/iquique_pgd/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/iquique_pgd/_output/fort.FG1.aux1'
fgmax_file=u'/Users/dmelgar/Tsunamis/iquique_wphase/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/iquique_wphase/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/iquique_landonly_long/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/iquique_landonly_long/_output/fort.FG1.aux1'
fgmax_file2=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.aux1'
fout='/Users/dmelgar/Iquique2014/tsunami/warnings/wphase.warn'
wet_tol=0.001 #Anything udner this is consdered dry
blend=False
etaclip=6



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
i=where(eta>etaclip)[0]
eta[i]=nan

if blend==True:
    #Get maximum amplitude first
    lon2=genfromtxt(fgmax_file2,usecols=0)-360
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
    latbounds=[-40,-32.5]
    i=where(lat<latbounds[0])[0] #THings south of the fine model
    lat_out=lat[i]
    lon_out=lon[i]
    eta_out=eta[i]
    #THings IN the fine model
    lat_out=r_[lat_out,lat2]
    lon_out=r_[lon_out,lon2]
    eta_out=r_[eta_out,eta2]
    #Things north of the fine model
    i=where(lat>latbounds[1])[0] #THings south of the fine model
    lat_out=r_[lat_out,lat[i]]
    lon_out=r_[lon_out,lon[i]]
    eta_out=r_[eta_out,eta[i]]
    #Ok done
    lat=lat_out
    lon=lon_out
    eta=eta_out
else:
    lat_out=lat
    eta_out=eta
    lon_out=lon






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
arec=region()
atac=region()
icax=region()
moqu=region()
tacn=region()
tara=region()


#Now read all regions

anto.lon,anto.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_antofagasta.txt')
arec.lon,arec.lat=read_region('/Users/dmelgar/KMLs/peru/regions/PER_arequipa.txt')
atac.lon,atac.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_atacama.txt')
icax.lon,icax.lat=read_region('/Users/dmelgar/KMLs/peru/regions/PER_ica.txt')
moqu.lon,moqu.lat=read_region('/Users/dmelgar/KMLs/peru/regions/PER_moquegua.txt')
tacn.lon,tacn.lat=read_region('/Users/dmelgar/KMLs/peru/regions/PER_tacna.txt')
tara.lon,tara.lat=read_region('/Users/dmelgar/KMLs/chile/regions/CHL_tarapaca.txt')

regions=[anto,arec,atac,icax,moqu,tacn,tara]

#Loop over every coast point, see what region it belongs to
for k in range(len(lon_out)):
    d=[]
    print k
    for kregion in range(len(regions)):
        d.append(sqrt((lon_out[k]-regions[kregion].lon)**2+(lat_out[k]-regions[kregion].lat)**2).min())
    if argmin(d)==0:
        anto.eta.append(eta_out[k])
    if argmin(d)==1:
        arec.eta.append(eta_out[k])
    if argmin(d)==2:
        atac.eta.append(eta_out[k])
    if argmin(d)==3:
        icax.eta.append(eta_out[k])
    if argmin(d)==4:
        moqu.eta.append(eta_out[k])
    if argmin(d)==5:
        tacn.eta.append(eta_out[k])
    if argmin(d)==6:
        tara.eta.append(eta_out[k])
        
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
#[anto,arec,atac,icax,moqu,tacn,tara]
anto.color=warn_color(anto)
arec.color=warn_color(arec)
atac.color=warn_color(atac)
icax.color=warn_color(icax)
moqu.color=warn_color(moqu)
tacn.color=warn_color(tacn)
tara.color=warn_color(tara)

#Write warning file
f=open(fout,'w')
f.write('Antofagasta='+anto.color+'\n')
f.write('Arequipa='+arec.color+'\n')
f.write('Atacama='+atac.color+'\n')
f.write('Ica='+icax.color+'\n')
f.write('Moquegua='+moqu.color+'\n')
f.write('Tacna='+tacn.color+'\n')
f.write('Tarapaca='+tara.color+'\n')
f.close()