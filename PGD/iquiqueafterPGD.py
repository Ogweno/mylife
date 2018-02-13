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
get_pgd=False
compile_pgds=False
get_magnitude=True
plot=True


#Globals
norm=1
path='/Users/dmelgar/PGD/'
data_path='/Users/dmelgar/PGD/GPS/sac/IquiqueAfter2014'
earthquakes=glob(path+'/station_info/*')

suffix='all'
pgd_folders=['/Users/dmelgar/PGD/GPS/sac/IquiqueAfter2014',]


Sspeed=6000 #S-wave mask
tmax=timedelta(seconds=100)
minpgd=0.05 #in m
maxdist=500 #in km
minsta=2
Nboot=1000
hypo=[-70.493,-20.571,22.4]
delta=0.05

#Get PGD as a function of time for all sites
if get_pgd:
    #initalize
    event='iquiqueafter'
    print event
    #Get hypocenter time
    hypo_time=UTCDateTime(genfromtxt(path+'event_info/'+event+'.hypo',dtype='S')[0])
    #Read station data
    stanames=genfromtxt(path+'station_info/iquique_after.sta',usecols=0,dtype='S')
    station_coords=genfromtxt(path+'station_info/iquique_after.sta',usecols=[1,2])
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
                d,az,baz=gps2DistAzimuth(station_coords[ksta,1],station_coords[ksta,0],hypo[1],hypo[0])
                hypo_dist=((d/1000)**2+hypo[2]**2)**0.5
                #Now get PGD at that site
                pgd_all_out=0
                pgd_nov_out=0
                time_to_all=0
                time_to_nov=0
                f=open(path+'PGD/IquiqueAfter2014/'+stanames[ksta]+'.pgd','w')
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
    event='IquiqueAfter2014'
    #Get list of PGD files
    stations=glob(path+'PGD/'+event+'/*.pgd')
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
        savetxt(path+'PGD/'+event+'.obs.pgd',c_[R,dnov*100,dall*100],fmt='%10.3f\t%10.4f\t%10.4f',header='dist(km), pgd no vert(cm) , pgd all (cm)')




#Get magnitude
if get_magnitude:
    event='IquiqueAfter2014'
    print event
    #Get list of PGD files
    stations=glob(path+'PGD/'+event+'/*.pgd')
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
    for kt in range(len(tpgd)):
        print kt
        #read all sites
        dall=[]
        dnov=[]
        R=[]
        #print kt
        cleansta=genfromtxt(path+'misc/'+event+'.clean',usecols=0,dtype='S')
        cleanvalue=genfromtxt(path+'misc/'+event+'.clean',usecols=1)
        for ksta in range(len(stations)):
            print stations[ksta]
            smask=False
            min_pgd_ok=False
            sta=stations[ksta].split('/')[-1].split('.')[0]
            pgd=genfromtxt(stations[ksta])
            iclean=where((sta.lower()==cleansta) | (sta.upper()==cleansta))[0]
            #Check swave mask
            dist=pgd[-1,5]
            if dist>(Sspeed*tpgd[kt]): #Swave not there yet
                smask=True
            if dist<maxdist:
                dist_use=True
            if pgd[kt,2]>minpgd:
                min_pgd_ok=True
            if cleanvalue[iclean]==1 and smask==False and dist_use==True and min_pgd_ok==True: #Station is to be used
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
            Mw_nov[kt]=scaling.PGD((dnov+delta)*100,R,coefficients=coeff_novert,weight=False,norm=norm)
            Mw_all[kt]=scaling.PGD((dall+delta)*100,R,coefficients=coeff_all,weight=False,norm=norm)
            #Now get upper and lower magnitude bound
            Mw_nov_plus[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert+sigma_novert,weight=True,norm=norm)
            Mw_all_plus[kt]=scaling.PGD(dall*100,R,coefficients=coeff_all+sigma_all,weight=True,norm=norm)
            Mw_nov_minus[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert-sigma_novert,weight=True,norm=norm)
            Mw_all_minus[kt]=scaling.PGD(dall*100,R,coefficients=coeff_all-sigma_all,weight=True,norm=norm)
            
        else:
            Mw_nov[kt]=nan
            Mw_all[kt]=nan
    #Write to file
    savetxt(path+'PGD/magnitudes/'+event+'.mag',c_[tpgd,Mw_nov,Mw_all,Mw_nov_minus,Mw_nov_plus,Mw_all_minus,Mw_all_plus],fmt='%5.1f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f',header='t,Mw_nov,Mw_all,Mw_nov_minus,Mw_nov_plus,Mw_all_minus,Mw_all_plus')


if plot:
    event='IquiqueAfter2014'
    M=genfromtxt('/Users/dmelgar/PGD/PGD/magnitudes/IquiqueAfter2014.mag')
    plt.figure()
    merr=0.25*ones(len(M))
    #plt.errorbar(M[:,0],M[:,1],yerr=[M[:,1]-M[:,3],M[:,4]-M[:,1]],fmt='o',color='#FF6347',ecolor='#FF6347',markersize=3,label=r'No vertical')
    plt.errorbar(M[:,0],M[:,1],yerr=[merr,merr],fmt='o',color='#FF6347',ecolor='#FF6347',markersize=3,label=r'No vertical')
    plt.errorbar(M[:,0],M[:,2],yerr=[merr,merr],fmt='o',color='#6495ED',ecolor='#6495ED',markersize=3,label=r'With vertical')
    plt.legend(loc=5)
    plt.ylim([6,8.0])
    plt.grid()
    plt.title('PGD magnitude, 2014 Mw7.7 Iquique Aftershock',fontsize=16)
    plt.ylabel(r'$M_w$',fontsize=16)
    plt.xlabel(r'Seconds after origin time',fontsize=16)
    #plt.plot([5,5],[0,10],'k--',lw=2)
    #plt.annotate('5s after OT',xy=(6,5.3))
    plt.show()
    
    
    
