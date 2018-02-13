from numpy import cos,sin,deg2rad,r_,c_,savetxt,where,zeros,unique,array,cross,rad2deg,arctan,arcsin,arange
from numpy.linalg import norm
from gmsh_tools import llz2utm

f=open('/Users/dmelgar/Cascadia/schmalzle/gaus_flt_atr.gmt','r')
fout='/Users/dmelgar/Cascadia/schmalzle/gaus_fine.txt'
fout2='/Users/dmelgar/Cascadia/schmalzle/gaus_mudpy.txt'
fout3='/Users/dmelgar/Cascadia/schmalzle/gaus.rupt'

#f=open('/Users/dmelgar/Cascadia/schmalzle/wang_flt_atr.gmt','r')
#fout='/Users/dmelgar/Cascadia/schmalzle/wang_fine.txt'
#fout2='/Users/dmelgar/Cascadia/schmalzle/wang_mudpy.txt'
#fout3='/Users/dmelgar/Cascadia/schmalzle/wang.rupt'

filter=True
thresh=1.0 #in mm

def get_strike_and_dip(p,q,r):
    
    #convert points to utm
    x,y=llz2utm(p[0],p[1],projection_zone='10T')
    p[0]=x/1000
    p[1]=y/1000
    x,y=llz2utm(q[0],q[1],projection_zone='10T')
    q[0]=x/1000
    q[1]=y/1000
    x,y=llz2utm(r[0],r[1],projection_zone='10T')
    r[0]=x/1000
    r[1]=y/1000
    
    #form refernce vectors
    pq=q-p
    pr=r-p
    
    #get normal
    n=cross(pq,pr)
    
    
    #get strike and dip
    a=n[0]
    b=n[1]
    c=n[2]
    d=a*p[0]+b*p[1]+c*p[2]
    #Compute strike
    strike=rad2deg(arctan(-b/a))
    #Compute dip
    beta=deg2rad(strike+90)
    m=array([sin(beta),cos(beta),0]) #Points in dip direction
    dip=abs(rad2deg(arcsin(m.dot(n)/(norm(m)*norm(n)))))
    if strike<0:
        strike=360+strike
    
    return strike,dip
    





k=0
while True:
    line=f.readline()
    if line=='':
        break
    if '>' in line:
        lon=float(line.split()[-10])
        lat=float(line.split()[-9])
        depth=-float(line.split()[-8])
        slip=float(line.splif=t()[3])
        phi=float(line.split()[4])
        rake=float(line.split()[-11])
        #split into ds and ss
        ds=slip*sin(deg2rad(rake))
        ss=slip*cos(deg2rad(rake))
        nodex=int(line.split()[-7])
        nodez=int(line.split()[-6])
        segx=int(line.split()[-5])
        segz=int(line.split()[-4])
        area=float(line.split()[-3])
        if k==0:
            lon_out=lon
            lat_out=lat
            depth_out=depth
            phi_out=phi
            ss_out=ss
            ds_out=ds
            nodex_out=nodex
            nodez_out=nodez
            segx_out=segx
            segz_out=segz
            area_out=area
        else:
            lon_out=r_[lon_out,lon]
            lat_out=r_[lat_out,lat]
            depth_out=r_[depth_out,depth]
            phi_out=r_[phi_out,phi]
            ss_out=r_[ss_out,ss]
            ds_out=r_[ds_out,ds]
            nodex_out=r_[nodex_out,nodex]
            nodez_out=r_[nodez_out,nodez]
            segx_out=r_[segx_out,segx]
            segz_out=r_[segz_out,segz]
            area_out=r_[area_out,area]
            
    else: #read next 4 lines which are quadranlge points
        p=array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
        line=f.readline()
        q=array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
        line=f.readline()
        r=array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
        line=f.readline()
        #get strike and dip
        strike,dip=get_strike_and_dip(p,q,r)
        if k==0:
            strike_out=strike
            dip_out=dip
        else:
            strike_out=r_[strike_out,strike]
            dip_out=r_[dip_out,dip]
        k+=1
f.close()

#Save something I can use to make plots and maps
out=c_[lon_out,lat_out,depth_out,phi_out,ss_out,ds_out]
savetxt(fout,out,header='lon,lat,depth(km),phi,ss(mm/yr),ds(mm/yr)',fmt='%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f')





#Save something for mudpy (only main nodes)
Nx=len(unique(nodex_out))
Nz=len(unique(nodez_out))
N=Nx*Nz

lon_centroid=zeros(N)
lat_centroid=zeros(N)
depth_centroid=zeros(N)
ss_centroid=zeros(N)
ds_centroid=zeros(N)
strike_centroid=zeros(N)
dip_centroid=zeros(N)
total_area=zeros(N)

k=0
for ix in unique(nodex_out):
    for iz in unique(nodez_out):
        
        #Find points that make the node
        i=where((nodex_out==ix) & (nodez_out==iz))[0]
        N=len(i)
        
        #Get features of the node
        lon_centroid[k]=lon_out[i].sum()/N
        lat_centroid[k]=lat_out[i].sum()/N
        depth_centroid[k]=depth_out[i].sum()/N
        ss_centroid[k]=ss_out[i].sum()/N
        ds_centroid[k]=ds_out[i].sum()/N
        strike_centroid[k]=strike_out[i].sum()/N
        dip_centroid[k]=dip_out[i].sum()/N
        total_area[k]=area_out[i].sum()
        k+=1
       


#write the gorram file
if filter==True:
    slip=(ss_centroid**2+ds_centroid**2)**0.5
    i=where(slip>thresh)[0]
else: #Keep all
    i=arange(len(lon_centroid))

#Apply filter
lon_centroid=lon_centroid[i]
lat_centroid=lat_centroid[i]
depth_centroid=depth_centroid[i]
strike_centroid=strike_centroid[i]
dip_centroid=dip_centroid[i]
ss_centroid=ss_centroid[i]
ds_centroid=ds_centroid[i]
total_area=total_area[i]

f=open(fout2,'w')
f.write('# No,lon,lat,z[km],strike,dip,rtype,rise time,len(m),width(m)\n')
for k in range(len(i)):
    line='%d\t%.4f\t%12.4f\t%12.4f\t%8.2f\t%6.2f\t0.5\t1.0\t%.2f\t%.2f\n' % (k+1,lon_centroid[k],lat_centroid[k],depth_centroid[k],strike_centroid[k],dip_centroid[k],(total_area[k]**0.5)*1000,(total_area[k]**0.5)*1000)
    f.write(line)
f.close()  

  
#Now output in rupture format  

f=open(fout3,'w')
f.write('# No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)\n')
for k in range(len(lon_centroid)):
    line='%d\t%.4f\t%12.4f\t%12.4f\t%8.2f\t%6.2f\t0.5\t1.0\t%10.6f\t%10.6f%9.2f\t%9.2f\t0.0\t60e9\n' % (k+1,lon_centroid[k],lat_centroid[k],depth_centroid[k],strike_centroid[k],dip_centroid[k],-ss_centroid[k]/1000.,-ds_centroid[k]/1000.,(total_area[k]**0.5)*1000,(total_area[k]**0.5)*1000)
    f.write(line)
f.close() 