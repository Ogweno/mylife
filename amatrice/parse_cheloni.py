from numpy import array,c_,r_,savetxt,cos,sin,deg2rad,arctan,expand_dims,arange,ones,zeros
from pyproj import Proj
from mudpy import forward

##For GMT
#
#fault='/Users/dmelgar/Amatrice2016/2009events/fault.proj_F2T1_tt.gd'
#fout='/Users/dmelgar/Amatrice2016/2009events/aftershock_slip.txt'
#
#f=open(fault,'r')
#k=0
#while True:
#    
#    line=f.readline()
#    
#    if '>' in line:
#        
#        #Read vertices
#        line1=f.readline()
#        line2=f.readline()
#        line3=f.readline()
#        line4=f.readline()
#        vertex1=array([float(line1.split()[0]),float(line1.split()[1]),float(line1.split()[2])])
#        vertex2=array([float(line2.split()[0]),float(line2.split()[1]),float(line2.split()[2])])
#        vertex3=array([float(line3.split()[0]),float(line3.split()[1]),float(line3.split()[2])])
#        vertex4=array([float(line4.split()[0]),float(line4.split()[1]),float(line4.split()[2])])
#        centroid=(vertex1+vertex2+vertex3+vertex4)/4
#        
#        if k==0:
#            slip=array([float(line.split()[2])/1000])
#            lon=array([centroid[0]])
#            lat=array([centroid[1]])
#            k+=1
#        else:
#            slip=r_[slip,float(line.split()[2])/1000]
#            lon=r_[lon,centroid[0]]
#            lat=r_[lat,centroid[1]]
#
#    if line=='':
#        break
#        
#f.close()
#savetxt(fout,c_[lon,lat,slip],fmt='%.4f\t%.4f\t%.4f')
        
        
##For coulomb
#reference_point=[13.3,42.6]
#fault='/Users/dmelgar/Amatrice2016/2009events/fault.proj_F1T1_tt.gd'
#coulout='/Users/dmelgar/Amatrice2016/2009events/mainshock_coul.txt'
##coulout='/Users/dmelgar/Amatrice2016/2009events/aftershock_coul.txt'
#dip=47 #mainshock
##dip=60 #aftershock
#
#zone='33T'
#p = Proj(proj='utm',zone=zone,ellps='WGS84')
#reference_utm=p(reference_point[0],reference_point[1])
#reference_utm=array(reference_utm)/1000
#
#f=open(fault,'r')
#k=0
#while True:
#    
#    line=f.readline()
#    
#    if '>' in line:
#        
#        #Read vertices
#        line1=f.readline()
#        line2=f.readline()
#        line3=f.readline()
#        line4=f.readline()
#        vertex1=array([float(line1.split()[0]),float(line1.split()[1]),float(line1.split()[2])])
#        vertex2=array([float(line2.split()[0]),float(line2.split()[1]),float(line2.split()[2])])
#        vertex3=array([float(line3.split()[0]),float(line3.split()[1]),float(line3.split()[2])])
#        vertex4=array([float(line4.split()[0]),float(line4.split()[1]),float(line4.split()[2])])
#        #Convert to UTM
#        v1=p(vertex1[0],vertex1[1])
#        xstart=v1[0]/1000-reference_utm[0]
#        ystart=v1[1]/1000-reference_utm[1]
#        v2=p(vertex2[0],vertex2[1])
#        xfin=v2[0]/1000-reference_utm[0]
#        yfin=v2[1]/1000-reference_utm[1]
#        
#        #Get slip
#        rake=float(line.split()[-1])
#        slip=float(line.split()[2])/1000
#        ss=slip*cos(deg2rad(rake))
#        ds=slip*sin(deg2rad(rake))
#        
#        #Coords of top and bottom
#        top=vertex1[2]
#        bottom=vertex4[2]
#        
#        if k==0:
#            k+=1
#            coul=expand_dims(array([1,xstart,ystart,xfin,yfin,100,ss,ds,dip,top,bottom,k]),0)
#        else:
#            k+=1
#            coul=r_[coul,expand_dims(array([1,xstart,ystart,xfin,yfin,100,ss,ds,dip,top,bottom,k]),0)]
#
#    if line=='':
#        break
#        
#f.close()
#savetxt(coulout,coul,fmt='%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.1f\t%.4f\t%.4f\t%d')





#For Mudpy
#strike=140 #Main
#dip=47
strike=150 #After
dip=60
fault='/Users/dmelgar/Amatrice2016/2009events/fault.proj_F2T1_tt.gd'
fout='/Users/dmelgar/Amatrice2016/2009events/aftershock.rupt'

reference_point=[13.3,42.6]
zone='33T'
p = Proj(proj='utm',zone=zone,ellps='WGS84')
reference_utm=p(reference_point[0],reference_point[1])
reference_utm=array(reference_utm)/1000

f=open(fault,'r')
k=0
while True:
    
    line=f.readline()
    
    if '>' in line:
        
        #Read vertices
        line1=f.readline()
        line2=f.readline()
        line3=f.readline()
        line4=f.readline()
        vertex1=array([float(line1.split()[0]),float(line1.split()[1]),float(line1.split()[2])])
        vertex2=array([float(line2.split()[0]),float(line2.split()[1]),float(line2.split()[2])])
        vertex3=array([float(line3.split()[0]),float(line3.split()[1]),float(line3.split()[2])])
        vertex4=array([float(line4.split()[0]),float(line4.split()[1]),float(line4.split()[2])])
        centroid=(vertex1+vertex2+vertex3+vertex4)/4
        
        #Convert to UTM
        v1=p(vertex1[0],vertex1[1])
        xstart=v1[0]/1000-reference_utm[0]
        ystart=v1[1]/1000-reference_utm[1]
        v2=p(vertex2[0],vertex2[1])
        xfin=v2[0]/1000-reference_utm[0]
        yfin=v2[1]/1000-reference_utm[1]
        length=((xstart-xfin)**2+(ystart-yfin)**2)**0.5
        width=abs(vertex3[2]-vertex1[2])/sin(deg2rad(dip))
        
        if k==0:
            slip=array([float(line.split()[2])/1000])
            rake=array([float(line.split()[-1])])
            L=array([length])
            W=array([width])
            lon=array([centroid[0]])
            lat=array([centroid[1]])
            depth=array([centroid[2]])
            k+=1
        else:
            slip=r_[slip,float(line.split()[2])/1000]
            rake=r_[rake,float(line.split()[-1])]
            lon=r_[lon,centroid[0]]
            lat=r_[lat,centroid[1]]
            L=r_[L,length]
            W=r_[W,width]
            depth=r_[depth,array([centroid[2]])]
    if line=='':
        break
       
ss=slip*cos(deg2rad(rake))
ds=slip*sin(deg2rad(rake))
f.close()
out=c_[arange(1,len(lon)+1),lon,lat,depth,strike*ones(len(lon)),dip*ones(len(lon)),0.5*ones(len(lon)),ones(len(lon)),ss,ds,L*1000,W*1000,zeros(len(lon)),zeros(len(lon))]
fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e\t%10.1f\t%10.1f\t%8.4f\t%.4e'
savetxt(fout,out,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')





#Do slip inversion
#rupt='/Users/dmelgar/Slip_inv/Amatrice_3D_final1/output/inverse_models/models/big_kahuna_vrtest_epishift3win_vr2.4.0016.inv.total'
#slip_coul_out='/Users/dmelgar/Slip_inv/Amatrice_3D_final1/output/inverse_models/models/big_kahuna_vrtest_epishift3win_vr2.4.0016.inv.total.coul'
#forward.inv2coulomb(rupt,reference_point,slip_coul_out)