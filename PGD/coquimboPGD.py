from glob import glob
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.core import UTCDateTime
from obspy import read
from numpy import genfromtxt,r_,where,arange,zeros,savetxt,c_,ones,nan,unique,meshgrid
from matplotlib import pyplot as plt
from datetime import timedelta
from string import rjust
import scaling
from scipy.interpolate import griddata
from matplotlib import cm




#Things to do
get_pgd=True
compile_pgds=True
get_magnitude=True
plot=False


#Globals
norm=1
path='/Users/dmelgar/PGD/'
data_path='/Users/dmelgar/PGD/GPS/sac/Coquimbo2015'
earthquakes=glob(path+'/station_info/*')

suffix='all'
pgd_folders=['/Users/dmelgar/PGD/GPS/sac/Coquimbo2015',]


Sspeed=3.0 #S-wave mask
tmax=timedelta(seconds=5.*60)
minpgd=2 #in cm
maxdist=1000 #in km
minsta=2
Nboot=1000
hypo=[287.73,-31.22,18.00]
dl=0.25
L=1.5 #degrees
dz=3
Lz=12
#Make search grid
lon_grid=arange(hypo[0]-L,hypo[0]+L,dl)
lat_grid=arange(hypo[1]-L,hypo[1]+L,dl)
z_grid=array([18])#arange(hypo[2]-Lz,hypo[0]+Lz,dz)
lat_search=ones(len(lon_grid)**2*len(z_grid))
lon_search=ones(len(lon_grid)**2*len(z_grid))
z_search=ones(len(lon_grid)**2*len(z_grid))
k=0
for kz in range(len(z_grid)):
    for klon in range(len(lon_grid)):
        for klat in range(len(lat_grid)):
            lat_search[k]=lat_grid[klat]
            lon_search[k]=lon_grid[klon]
            z_search[k]=z_grid[kz]
            k+=1


#Get PGD as a function of time for all sites
if get_pgd:
    #initalize
    event='coquimbo'
    print event
    #Get hypocenter time
    hypo_time=UTCDateTime(genfromtxt(path+'event_info/'+event+'.hypo',dtype='S')[0])
    #Read station data
    stanames=genfromtxt(path+'station_info/coquimbo.sta',usecols=0,dtype='S')
    station_coords=genfromtxt(path+'station_info/coquimbo.sta',usecols=[1,2])
    for khypo in range(len(lon_search)): #Loops pver grid search nodes
        print 'Working on node '+str(khypo)+'/'+str(len(lon_search)) 
        iteration=rjust(str(khypo),4,'0')
        #Now loop over station data at each epoch
        for ksta in range(len(stanames)):
            #print '... fetching PGD for '+stanames[ksta]
            try:
                n=read(data_path+'/'+stanames[ksta]+'.LXN.sac')
                e=read(data_path+'/'+stanames[ksta]+'.LXE.sac')
                u=read(data_path+'/'+stanames[ksta]+'.LXZ.sac')
                if n[0].stats.npts>0:
                    #Trim to times of interest
                    n.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    e.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    u.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    #Remove mean of first ten seconds
                    #Get station-event distance
                    d,az,baz=gps2DistAzimuth(station_coords[ksta,1],station_coords[ksta,0],lat_search[khypo],lon_search[khypo])
                    hypo_dist=((d/1000)**2+hypo[2]**2)**0.5
                    #Now get PGD at that site
                    pgd_all_out=0
                    pgd_nov_out=0
                    time_to_all=0
                    time_to_nov=0
                    f=open(path+'PGD/Coquimbo2015/_gridsearch/'+stanames[ksta]+'.'+iteration+'.pgd','w')
                    f.write('# Seconds after OT , no horiz PGD [m] , all PGD [m] , seconds to horiz , seconds to all , hypo_dist [km]\n')
                    for kt in range(len(n[0].times())):
                        pgd_all=(n[0].data[kt]**2+e[0].data[kt]**2+u[0].data[kt]**2)**0.5
                        pgd_nov=(n[0].data[kt]**2+e[0].data[kt]**2)**0.5
                        if pgd_all>pgd_all_out:
                            pgd_all_out=pgd_all
                            time_to_all=kt
                        if pgd_nov>pgd_nov_out:
                            pgd_nov_out=pgd_nov
                            time_to_nov=kt
                        line='%5.1f\t%8.4f%8.4f\t%5.1f\t%5.1f\t%7.2f\n' %(kt,pgd_nov_out,pgd_all_out,time_to_nov,time_to_all,hypo_dist)
                        f.write(line)
                    f.close()
            except:
                print 'Station '+stanames[ksta]+' does not exist'

# Gather PGD's from all sites for each earthquake
if compile_pgds:
    event='Coquimbo2015'
    for khypo in range(len(lon_search)): #Loops pver grid search nodes
        print 'Working on node '+str(khypo)+'/'+str(len(lon_search)) 
        iteration=rjust(str(khypo),4,'0')
        #Get list of PGD files
        stations=glob(path+'PGD/'+event+'/_gridsearch/'+'*'+iteration+'.pgd')
        #Read PGD
        dall=[]
        dnov=[]
        R=[]
        #Open cleanup file
        cleansta=genfromtxt(path+'misc/'+event+'.clean',usecols=0,dtype='S')
        cleanvalue=genfromtxt(path+'misc/'+event+'.clean',usecols=1)
        for ksta in range(len(stations)):
            sta=stations[ksta].split('/')[-1].split('.')[0]
            iclean=where((sta.lower()==cleansta) | (sta.upper()==cleansta))[0]
            if cleanvalue[iclean]==1: #Station is to be used
                pgd=genfromtxt(stations[ksta])
                #Add to data vector
                dnov=r_[dnov,pgd[-1,1]]
                dall=r_[dall,pgd[-1,2]]
                R=r_[R,pgd[-1,5]]
            savetxt(path+'PGD/_gridsearch/'+event+'.'+iteration+'.obs.pgd',c_[R,dnov*100,dall*100],fmt='%10.3f\t%10.4f\t%10.4f',header='dist(km), pgd no vert(cm) , pgd all (cm)')




#Get magnitude
if get_magnitude:
    event='Coquimbo2015'
    print event
    for khypo in range(len(lon_search)): #Loops pver grid search nodes
        print 'Working on node '+str(khypo)+'/'+str(len(lon_search)) 
        iteration=rjust(str(khypo),4,'0')
        #Get list of PGD files
        stations=glob(path+'PGD/'+event+'/_gridsearch/'+'*'+iteration+'.pgd')
        #Now loop over times
        tpgd=arange(0,tmax.seconds,1)
        Mw_nov=zeros(len(tpgd))
        Mw_nov_minus=zeros(len(tpgd))
        Mw_nov_plus=zeros(len(tpgd))
        Mw_all=zeros(len(tpgd))
        Mw_all_minus=zeros(len(tpgd))
        Mw_all_plus=zeros(len(tpgd))
        res_nov=zeros(len(tpgd))
        res_all=zeros(len(tpgd))
        for kt in [290]:#(len(tpgd)):
            #read all sites
            dall=[]
            dnov=[]
            R=[]
            #print kt
            cleansta=genfromtxt(path+'misc/'+event+'.clean',usecols=0,dtype='S')
            cleanvalue=genfromtxt(path+'misc/'+event+'.clean',usecols=1)
            for ksta in range(len(stations)):
                smask=False
                sta=stations[ksta].split('/')[-1].split('.')[0]
                pgd=genfromtxt(stations[ksta])
                iclean=where((sta.lower()==cleansta) | (sta.upper()==cleansta))[0]
                #Check swave mask
                dist=pgd[-1,5]
                if dist>(Sspeed*tpgd[kt]): #Swave not there yet
                    smask=True
                if dist<maxdist:
                    dist_use=True
                if cleanvalue[iclean]==1 and smask==False and dist_use==True: #Station is to be used
                    i=where(pgd[:,0]==tpgd[kt])[0]
                    #Add to data vector
                    dnov=r_[dnov,pgd[i,1]]
                    dall=r_[dall,pgd[i,2]]
                    R=r_[R,pgd[i,5]]
            #Solve regression
            coeff_novert=genfromtxt(path+'coefficients/'+'novert.coeff')
            coeff_all=genfromtxt(path+'coefficients/'+'all.coeff')
            #Get coefficient  uncertainties
            sigma_novert=1.96*genfromtxt(path+'coefficients/novert.boot',usecols=1)
            sigma_novert[1]=-sigma_novert[1] ; sigma_novert[2]=-sigma_novert[2]
            sigma_all=1.96*genfromtxt(path+'coefficients/novert.boot',usecols=1)
            sigma_all[1]=-sigma_all[1] ; sigma_all[2]=-sigma_all[2]
            if len(dnov)>=minsta:
                #print "using "+str(len(dnov))+' stations'
                #Mw_nov[kt],res_nov[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert,weight=False,norm=2,residual=True)
                Mw_all,res_all=scaling.PGD(dall*100,R,coefficients=coeff_all,weight=False,norm=2,residual=True)
                #Now get upper and lower magnitude bound
                #Mw_nov_plus[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert+sigma_novert,weight=True,norm=2)
                #Mw_all_plus[kt]=scaling.PGD(dall*100,R,coefficients=coeff_all+sigma_all,weight=True,norm=2)
                #Mw_nov_minus[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert-sigma_novert,weight=True,norm=2)
                #Mw_all_minus[kt]=scaling.PGD(dall*100,R,coefficients=coeff_all-sigma_all,weight=True,norm=2)
                
            else:
                Mw_nov[kt]=nan
                Mw_all[kt]=nan
        Mw_nov=Mw_all
        Mw_nov_plus=Mw_all
        Mw_nov_minus=Mw_all
        Mw_all_plus=Mw_all
        Mw_all_minus=Mw_all
        res_nov=res_all
        #Write to file
        savetxt(path+'PGD/magnitudes/_gridsearch/'+event+'.'+iteration+'.mag',c_[tpgd[kt],Mw_nov,Mw_all,Mw_nov_minus,Mw_nov_plus,Mw_all_minus,Mw_all_plus,res_nov,res_all],fmt='%5.1f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f',header='t,Mw_nov,Mw_all,Mw_nov_minus,Mw_nov_plus,Mw_all_minus,Mw_all_plus,res_nov,res_all')


if plot:
    event='Coquimbo2015'
    N=143
    Mw_nov=zeros(N)
    Mw_all=zeros(N)
    res_nov=zeros(N)
    res_all=zeros(N)
    longrid=zeros(N)
    latgrid=zeros(N)
    for k in range(N):#len(lon_search)): #Loops pver grid search nodes
        iteration=rjust(str(k),4,'0')
        #Read magnitude file
        #,c_[tpgd,Mw_nov,Mw_all,Mw_nov_minus,Mw_nov_plus,Mw_all_minus,Mw_all_plus,res_nov,res_all]
        mag=genfromtxt(path+'PGD/magnitudes/_gridsearch/'+event+'.'+iteration+'.mag')
        Mw_nov[k]=mag[1]
        Mw_all[k]=mag[2]
        res_nov[k]=mag[7]
        res_all[k]=mag[7]
        longrid[k]=lon_search[k]
        latgrid[k]=lat_search[k]
    #Convert to grid
    X,Y=meshgrid(arange(longrid.min(),longrid.max(),0.05),arange(latgrid.min(),latgrid.max(),0.05))
    Mw_nov_grid=griddata((longrid,latgrid),Mw_nov,(X,Y))
    Mw_all_grid=griddata((longrid,latgrid),Mw_all,(X,Y))
    res_nov_grid=griddata((longrid,latgrid),res_nov,(X,Y))
    res_all_grid=griddata((longrid,latgrid),res_all,(X,Y))
    #plot
    plt.figure()
    plt.contourf(X-360,Y,Mw_all_grid,80,cmap=cm.pink_r)
    cb=plt.colorbar()
    cb.ax.set_xlabel('Mw')
    plt.contour(X-360,Y,res_all_grid,50,cmap=cm.rainbow)
    cb2=plt.colorbar()
    cb2.ax.set_xlabel('Residual')
    plt.show()
    
    
    
