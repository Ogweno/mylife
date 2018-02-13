from numpy import genfromtxt,unique,zeros,where,nan,argmin,sqrt,array,r_,isnan
from glob import glob
    
#First read the simulation results

#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_gps_tg_30/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_gps_tg_30/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_gps_30/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_gps_30/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_pgd/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_pgd/_output/fort.FG1.aux1'
#fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_wphase/_output/fort.FG1.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/tohoku_wphase/_output/fort.FG1.aux1'
fgmax_file=u'/Users/dmelgar/Tsunamis/tohoku_glarms/_output/fort.FG1.valuemax'
aux_file='/Users/dmelgar/Tsunamis/tohoku_glarms/_output/fort.FG1.aux1'
fgmax_file2=u'/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.valuemax'
aux_file2='/Users/dmelgar/Tsunamis/maule_gps_tg_static_3/_output/fort.FG1.aux1'

fout='/Users/dmelgar/Tohoku2011/tsunami/warnings/glarms.warn'
wet_tol=0.001 #Anything udner this is consdered dry
blend=False
etaclip=40



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
aich=region()
akit=region()
aomo=region()
chib=region()
ehim=region()
fuki=region()
fuko=region()
fuks=region()
gifu=region()
gunm=region()
hiro=region()
hokk=region()
hyog=region()
ibar=region()
ishi=region()
iwat=region()
kaga=region()
kago=region()
kana=region()
koch=region()
kuma=region()
kyot=region()
miex=region()
miyg=region()
miyz=region()
naga=region()
naoa=region()
nara=region()
niig=region()
oita=region()
okay=region()
okin=region()
osak=region()
saga=region()
sait=region()
shig=region()
shim=region()
shiz=region()
toch=region()
toku=region()
toky=region()
tott=region()
toya=region()
waka=region()
yamg=region()
yamn=region()
yamu=region()



#Now read all regions

aich.lon,aich.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Aichi.txt')
akit.lon,akit.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Akita.txt')
aomo.lon,aomo.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Aomori.txt')
chib.lon,chib.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Chiba.txt')
ehim.lon,ehim.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Ehime.txt')
fuki.lon,fuki.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Fuki.txt')
fuko.lon,fuko.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Fukoka.txt')
fuks.lon,fuks.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Fukshima.txt')
gifu.lon,gifu.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Gifu.txt')
gunm.lon,gunm.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Gunma.txt')
hiro.lon,hiro.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Hiroshima.txt')
hokk.lon,hokk.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Hokkaido.txt')
hyog.lon,hyog.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Hyogo.txt')
ibar.lon,ibar.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Ibaraki.txt')
ishi.lon,ishi.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Ishikawa.txt')
iwat.lon,iwat.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Iwate.txt')
kaga.lon,kaga.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Kagawa.txt')
kago.lon,kago.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Kagoshima.txt')
kana.lon,kana.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Kanagawa.txt')
koch.lon,koch.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Kochi.txt')
kuma.lon,kuma.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Kumamoto.txt')
kyot.lon,kyot.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Kyoto.txt')
miex.lon,miex.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Miex.txt')
miyg.lon,miyg.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Miygi.txt')
miyz.lon,miyz.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Miyzaki.txt')
naga.lon,naga.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Nagano.txt')
naoa.lon,naoa.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Naoasaki.txt')
nara.lon,nara.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Nara.txt')
niig.lon,niig.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Niigata.txt')
oita.lon,oita.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Oita.txt')
okay.lon,okay.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Okayama.txt')
okin.lon,okin.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Okinawa.txt')
osak.lon,osak.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Osaka.txt')
saga.lon,saga.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Saga.txt')
sait.lon,sait.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Saitama.txt')
shig.lon,shig.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Shiga.txt')
shim.lon,shim.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Shimane.txt')
shiz.lon,shiz.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Shizuoka.txt')
toch.lon,toch.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Tochigi.txt')
toku.lon,toku.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Tokushima.txt')
toky.lon,toky.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Tokyo.txt')
tott.lon,tott.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Tottori.txt')
toya.lon,toya.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Toyama.txt')
waka.lon,waka.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Wakayama.txt')
yamg.lon,yamg.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Yamgata.txt')
yamn.lon,yamn.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Yamnashi.txt')
yamu.lon,yamu.lat=read_region('/Users/dmelgar/KMLs/japan/prefectures/JPN_Yamuchi.txt')


regions=[aich,akit,aomo,chib,ehim,fuki,fuko,fuks,gifu,gunm,hiro,hokk,hyog,ibar,ishi,iwat,kaga,kago,
    kana,koch,kuma,kyot,miex,miyg,miyz,naga,naoa,nara,niig,oita,okay,okin,osak,saga,sait,shig,shim,
    shiz,toch,toku,toky,tott,toya,waka,yamg,yamn,yamu]


#Loop over every coast point, see what region it belongs to
for k in range(len(lon_out)):
    d=[]
    print k
    for kregion in range(len(regions)):
        d.append(sqrt((lon_out[k]-regions[kregion].lon)**2+(lat_out[k]-regions[kregion].lat)**2).min())
    if argmin(d)==0:
   	aich.eta.append(eta_out[k])
    if argmin(d)==1:
   	akit.eta.append(eta_out[k])
    if argmin(d)==2:
   	aomo.eta.append(eta_out[k])
    if argmin(d)==3:
   	chib.eta.append(eta_out[k])
    if argmin(d)==4:
   	ehim.eta.append(eta_out[k])
    if argmin(d)==5:
   	fuki.eta.append(eta_out[k])
    if argmin(d)==6:
   	fuko.eta.append(eta_out[k])
    if argmin(d)==7:
   	fuks.eta.append(eta_out[k])
    if argmin(d)==8:
   	gifu.eta.append(eta_out[k])
    if argmin(d)==9:
   	gunm.eta.append(eta_out[k])
    if argmin(d)==10:
   	hiro.eta.append(eta_out[k])
    if argmin(d)==11:
   	hokk.eta.append(eta_out[k])
    if argmin(d)==12:
   	hyog.eta.append(eta_out[k])
    if argmin(d)==13:
   	ibar.eta.append(eta_out[k])
    if argmin(d)==14:
   	ishi.eta.append(eta_out[k])
    if argmin(d)==15:
   	iwat.eta.append(eta_out[k])
    if argmin(d)==16:
   	kaga.eta.append(eta_out[k])
    if argmin(d)==17:
   	kago.eta.append(eta_out[k])
    if argmin(d)==18:
   	kana.eta.append(eta_out[k])
    if argmin(d)==19:
   	koch.eta.append(eta_out[k])
    if argmin(d)==20:
   	kuma.eta.append(eta_out[k])
    if argmin(d)==21:
   	kyot.eta.append(eta_out[k])
    if argmin(d)==22:
   	miex.eta.append(eta_out[k])
    if argmin(d)==23:
   	miyg.eta.append(eta_out[k])
    if argmin(d)==24:
   	miyz.eta.append(eta_out[k])
    if argmin(d)==25:
   	naga.eta.append(eta_out[k])
    if argmin(d)==26:
   	naoa.eta.append(eta_out[k])
    if argmin(d)==27:
   	nara.eta.append(eta_out[k])
    if argmin(d)==28:
   	niig.eta.append(eta_out[k])
    if argmin(d)==29:
   	oita.eta.append(eta_out[k])
    if argmin(d)==30:
   	okay.eta.append(eta_out[k])
    if argmin(d)==31:
   	okin.eta.append(eta_out[k])
    if argmin(d)==32:
   	osak.eta.append(eta_out[k])
    if argmin(d)==33:
   	saga.eta.append(eta_out[k])
    if argmin(d)==34:
   	sait.eta.append(eta_out[k])
    if argmin(d)==35:
   	shig.eta.append(eta_out[k])
    if argmin(d)==36:
   	shim.eta.append(eta_out[k])
    if argmin(d)==37:
   	shiz.eta.append(eta_out[k])
    if argmin(d)==38:
   	toch.eta.append(eta_out[k])
    if argmin(d)==39:
   	toku.eta.append(eta_out[k])
    if argmin(d)==40:
   	toky.eta.append(eta_out[k])
    if argmin(d)==41:
   	tott.eta.append(eta_out[k])
    if argmin(d)==42:
   	toya.eta.append(eta_out[k])
    if argmin(d)==43:
   	waka.eta.append(eta_out[k])
    if argmin(d)==44:
   	yamg.eta.append(eta_out[k])
    if argmin(d)==45:
   	yamn.eta.append(eta_out[k])
    if argmin(d)==46:
   	yamu.eta.append(eta_out[k])
        
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
aich.color=warn_color(aich)
akit.color=warn_color(akit)
aomo.color=warn_color(aomo)
chib.color=warn_color(chib)
ehim.color=warn_color(ehim)
fuki.color=warn_color(fuki)
fuko.color=warn_color(fuko)
fuks.color=warn_color(fuks)
gifu.color=warn_color(gifu)
gunm.color=warn_color(gunm)
hiro.color=warn_color(hiro)
hokk.color=warn_color(hokk)
hyog.color=warn_color(hyog)
ibar.color=warn_color(ibar)
ishi.color=warn_color(ishi)
iwat.color=warn_color(iwat)
kaga.color=warn_color(kaga)
kago.color=warn_color(kago)
kana.color=warn_color(kana)
koch.color=warn_color(koch)
kuma.color=warn_color(kuma)
kyot.color=warn_color(kyot)
miex.color=warn_color(miex)
miyg.color=warn_color(miyg)
miyz.color=warn_color(miyz)
naga.color=warn_color(naga)
naoa.color=warn_color(naoa)
nara.color=warn_color(nara)
niig.color=warn_color(niig)
oita.color=warn_color(oita)
okay.color=warn_color(okay)
okin.color=warn_color(okin)
osak.color=warn_color(osak)
saga.color=warn_color(saga)
sait.color=warn_color(sait)
shig.color=warn_color(shig)
shim.color=warn_color(shim)
shiz.color=warn_color(shiz)
toch.color=warn_color(toch)
toku.color=warn_color(toku)
toky.color=warn_color(toky)
tott.color=warn_color(tott)
toya.color=warn_color(toya)
waka.color=warn_color(waka)
yamg.color=warn_color(yamg)
yamn.color=warn_color(yamn)
yamu.color=warn_color(yamu)


#Write warning file
f=open(fout,'w')
f.write('Aichi='+aich.color+'\n')
f.write('Akita='+akit.color+'\n')
f.write('Aomori='+aomo.color+'\n')
f.write('Chiba='+chib.color+'\n')
f.write('Ehime='+ehim.color+'\n')
f.write('Fuki='+fuki.color+'\n')
f.write('Fukoka='+fuko.color+'\n')
f.write('Fukshima='+fuks.color+'\n')
f.write('Gifu='+gifu.color+'\n')
f.write('Gunma='+gunm.color+'\n')
f.write('Hiroshima='+hiro.color+'\n')
f.write('Hokkaido='+hokk.color+'\n')
f.write('Hyogo='+hyog.color+'\n')
f.write('Ibaraki='+ibar.color+'\n')
f.write('Ishikawa='+ishi.color+'\n')
f.write('Iwate='+iwat.color+'\n')
f.write('Kagawa='+kaga.color+'\n')
f.write('Kagoshima='+kago.color+'\n')
f.write('Kanagawa='+kana.color+'\n')
f.write('Kochi='+koch.color+'\n')
f.write('Kumamoto='+kuma.color+'\n')
f.write('Kyoto='+kyot.color+'\n')
f.write('Miex='+miex.color+'\n')
f.write('Miygi='+miyg.color+'\n')
f.write('Miyzaki='+miyz.color+'\n')
f.write('Nagano='+naga.color+'\n')
f.write('Naoasaki='+naoa.color+'\n')
f.write('Nara='+nara.color+'\n')
f.write('Niigata='+niig.color+'\n')
f.write('Oita='+oita.color+'\n')
f.write('Okayama='+okay.color+'\n')
f.write('Okinawa='+okin.color+'\n')
f.write('Osaka='+osak.color+'\n')
f.write('Saga='+saga.color+'\n')
f.write('Saitama='+sait.color+'\n')
f.write('Shiga='+shig.color+'\n')
f.write('Shimane='+shim.color+'\n')
f.write('Shizuoka='+shiz.color+'\n')
f.write('Tochigi='+toch.color+'\n')
f.write('Tokushima='+toku.color+'\n')
f.write('Tokyo='+toky.color+'\n')
f.write('Tottori='+tott.color+'\n')
f.write('Toyama='+toya.color+'\n')
f.write('Wakayama='+waka.color+'\n')
f.write('Yamgata='+yamg.color+'\n')
f.write('Yamnashi='+yamn.color+'\n')
f.write('Yamuchi='+yamu.color+'\n')

f.close()