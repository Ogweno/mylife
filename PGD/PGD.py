'''
Diego Melgar 02/2015
Perform PGD regressions and magnitude computations
'''

from glob import glob
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.core import UTCDateTime
from obspy import read
from numpy import genfromtxt,r_,where,arange,zeros,savetxt,c_,ones,log10,expand_dims,array,logspace,nan,diag,setxor1d
from matplotlib import pyplot as plt
from datetime import timedelta
from string import rjust
import scaling
from numpy.linalg import lstsq
from scipy.linalg import norm as vecnorm
from l1 import l1
from cvxopt import matrix
from numpy.random import randint



#Things to do
get_pgd=False
compile_pgds=False
run_regression=False
get_magnitude=False
bootstrap=False
magnitude_error=False
#Things to plot
plot_sites_and_pgd=False
plot_all=False
plot_magnitudes=False
plot_pgd_dist=False
plot_pgd_dist_all=False
#Extras
clean=False

#Globals
norm=2
path='/Users/dmelgar/PGD/'
data_path='/Users/dmelgar/PGD/GPS/sac/'
earthquakes=glob(path+'/station_info/*')
# Hack for event ommitance testing
#pgd_folders=glob(path+'GPS/sac/*')
#suffix='Ignore_Tokachi'
#pgd_folders=['/Users/dmelgar/PGD/GPS/sac/Tohoku2011',
# '/Users/dmelgar/PGD/GPS/sac/Aegean2014',
# '/Users/dmelgar/PGD/GPS/sac/ElMayor2010',
# '/Users/dmelgar/PGD/GPS/sac/Iquique2014',
# '/Users/dmelgar/PGD/GPS/sac/Napa2014',
# '/Users/dmelgar/PGD/GPS/sac/Maule2010',
# '/Users/dmelgar/PGD/GPS/sac/Mentawai2010',
# '/Users/dmelgar/PGD/GPS/sac/Nicoya2012',
# '/Users/dmelgar/PGD/GPS/sac/Parkfield2004']
#pgd_folders_mag=['/Users/dmelgar/PGD/GPS/sac/Tokachi2003']
# End hack

suffix='all'
pgd_folders=['/Users/dmelgar/PGD/GPS/sac/Aegean2014',
 '/Users/dmelgar/PGD/GPS/sac/Coquimbo2015',
 '/Users/dmelgar/PGD/GPS/sac/ElMayor2010',
 '/Users/dmelgar/PGD/GPS/sac/Iquique2014',
 '/Users/dmelgar/PGD/GPS/sac/Maule2010',
 '/Users/dmelgar/PGD/GPS/sac/Mentawai2010',
 '/Users/dmelgar/PGD/GPS/sac/Napa2014',
 '/Users/dmelgar/PGD/GPS/sac/Nicoya2012',
 '/Users/dmelgar/PGD/GPS/sac/Parkfield2004',
 '/Users/dmelgar/PGD/GPS/sac/Tohoku2011',
 '/Users/dmelgar/PGD/GPS/sac/Tokachi2003']


Sspeed=3.0 #S-wave mask
tmax=timedelta(seconds=5.*60)
minpgd=1 #in cm
maxdist=1000 #in km
minsta=2
Nboot=1000




#Get PGD as a function of time for all sites
if get_pgd:
    #initalize
    for k in [1]:#range(len(earthquakes)):
        event=earthquakes[k].split('/')[-1].split('.')[0]
        print event
        #Get hypocenter
        hypo_time=UTCDateTime(genfromtxt(path+'event_info/'+event+'.hypo',dtype='S')[0])
        hypo=genfromtxt(path+'event_info/'+event+'.hypo')
        #Read station data
        stanames=genfromtxt(earthquakes[k],usecols=0,dtype='S')
        station_coords=genfromtxt(earthquakes[k],usecols=[1,2])
        #Now loop over station data at each epoch
        event=pgd_folders[k].split('/')[-1]
        for ksta in range(len(stanames)):
            print 'Fetching PGD for '+stanames[ksta]
            try:
                if event=='Tohoku2011' or event=='Tokachi2003':
                    n=read(pgd_folders[k]+'/'+rjust(stanames[ksta],4,'0')+'.LXN.sac')
                    e=read(pgd_folders[k]+'/'+rjust(stanames[ksta],4,'0')+'.LXE.sac')
                    u=read(pgd_folders[k]+'/'+rjust(stanames[ksta],4,'0')+'.LXZ.sac')
                else: 
                    n=read(pgd_folders[k]+'/'+stanames[ksta]+'.LXN.sac')
                    e=read(pgd_folders[k]+'/'+stanames[ksta]+'.LXE.sac')
                    u=read(pgd_folders[k]+'/'+stanames[ksta]+'.LXZ.sac')
                if n[0].stats.npts>0:
                    #Trim to times of interest
                    n.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    e.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    u.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    #decimate if necessary (5Hz data)
                    if n[0].stats.delta==0.2:
                        n.decimate(5, no_filter=True)
                        e.decimate(5, no_filter=True)
                        u.decimate(5, no_filter=True)
                    #Remove mean of first ten seconds
                    #Get station-event distance
                    d,az,baz=gps2DistAzimuth(station_coords[ksta,1],station_coords[ksta,0],hypo[2],hypo[1])
                    hypo_dist=((d/1000)**2+hypo[3]**2)**0.5
                    #Now get PGD at that site
                    pgd_all_out=0
                    pgd_nov_out=0
                    time_to_all=0
                    time_to_nov=0
                    f=open(path+'PGD/'+event+'/'+stanames[ksta]+'.pgd','w')
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
    #initalize
    for k in [1]:#range(0,len(pgd_folders)):
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        print event
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
            if event=='Tohoku2011' or event=='Tokachi2003':
                sta=rjust(sta,4,'0')
            iclean=where((sta.lower()==cleansta) | (sta.upper()==cleansta))[0]
            if cleanvalue[iclean]==1: #Station is to be used
                pgd=genfromtxt(stations[ksta])
                #Add to data vector
                dnov=r_[dnov,pgd[-1,1]]
                dall=r_[dall,pgd[-1,2]]
                R=r_[R,pgd[-1,5]]
            savetxt(path+'PGD/'+event+'.obs.pgd',c_[R,dnov*100,dall*100],fmt='%10.3f\t%10.4f\t%10.4f',header='dist(km), pgd no vert(cm) , pgd all (cm)')

if run_regression:
    #initalize
    for k in range(0,len(pgd_folders)):
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        #Assign magnitudes
        if event=='ElMayor2010':
            Mcurrent=7.18
        if event=='Iquique2014':
            Mcurrent=8.19
        if event=='Maule2010':
            Mcurrent=8.85
        if event=='Mentawai2010':
            Mcurrent=7.68
        if event=='Napa2014':
            Mcurrent=6.11
        if event=='Parkfield2004':
            Mcurrent=5.92
        if event=='Tohoku2011':
            Mcurrent=9.09
        if event=='Tokachi2003':
            Mcurrent=8.25
        if event=='Nicoya2012':
            Mcurrent=7.57
        if event=='Aegean2014':
            Mcurrent=6.87
        pgd_temp=genfromtxt(path+'PGD/'+event+'.obs.pgd')
        ifilt_novert=where((pgd_temp[:,0]<maxdist) & (pgd_temp[:,1]>minpgd))[0]
        ifilt_all=where((pgd_temp[:,0]<maxdist) & (pgd_temp[:,2]>minpgd))[0]
        if k==0:
            pgd_novert=pgd_temp[ifilt_novert,:].copy()
            Mw_novert=ones(len(pgd_novert))*Mcurrent
            W_novert=ones(len(pgd_novert))/vecnorm(log10(pgd_novert[:,1]))
            pgd_all=pgd_temp[ifilt_all,:].copy()
            Mw_all=ones(len(pgd_all))*Mcurrent
            W_all=ones(len(pgd_all))/vecnorm(log10(pgd_all[:,2]))
        else:
            pgd_novert=r_[pgd_novert,pgd_temp[ifilt_novert]]
            Mw_novert=r_[Mw_novert,ones(len(pgd_temp[ifilt_novert,:]))*Mcurrent]
            W_novert=r_[W_novert,ones(len(pgd_novert[ifilt_novert,:]))/vecnorm(log10(pgd_novert[ifilt_novert,1]))]
            pgd_all=r_[pgd_all,pgd_temp[ifilt_all]]
            Mw_all=r_[Mw_all,ones(len(pgd_temp[ifilt_all,:]))*Mcurrent]
            W_all=r_[W_all,ones(len(pgd_all[ifilt_all,:]))/vecnorm(log10(pgd_all[ifilt_all,2]))]
    #Define regression quantities
    dnovert=log10(pgd_novert[:,1])
    dall=log10(pgd_all[:,2])
    #make matrix of event weights
    W_all=diag(W_all)
    W_novert=diag(W_novert)
    #Make matrix of data weights
    iall=ones((len(dall),1))
    inovert=ones((len(dnovert),1))
    G_all=c_[iall,expand_dims(Mw_all,1)*iall,expand_dims(Mw_all*log10(pgd_all[:,0]),1)]
    G_novert=c_[inovert,expand_dims(Mw_novert,1)*inovert,expand_dims(Mw_novert*log10(pgd_novert[:,0]),1)]
    #Run regression
    # log(PGD)=A+B*Mw+C*Mw*log(R)
    if norm==2:
        coefficients_novert=lstsq(W_novert.dot(G_novert),W_novert.dot(dnovert))[0]
        coefficients_all=lstsq(W_all.dot(G_all),W_all.dot(dall))[0]
    elif norm==1:
        P=matrix(W_novert.dot(G_novert))
        q=matrix(W_novert.dot(dnovert))
        coefficients_novert=array(l1(P,q))
        P=matrix(W_all.dot(G_all))
        q=matrix(W_all.dot(dall))
        coefficients_all=array(l1(P,q))
    f=open(path+'coefficients/'+suffix+'.novert.coeff','w')
    f.write('#A,B,C\n%.6f\n%.6f\n%.6f' %(coefficients_novert[0],coefficients_novert[1],coefficients_novert[2]))
    f.close()
    f=open(path+'coefficients/'+suffix+'.all.coeff','w')
    f.write('#A,B,C\n%.6f\n%.6f\n%.6f' %(coefficients_all[0],coefficients_all[1],coefficients_all[2]))
    f.close()
    print 'All coefficients '+str(coefficients_all)
    print 'Novert coefficients '+str(coefficients_novert)
        

#Determine parameter uncertainties
if bootstrap:
    Aall=zeros(Nboot)
    Ball=zeros(Nboot)
    Call=zeros(Nboot)
    Anovert=zeros(Nboot)
    Bnovert=zeros(Nboot)
    Cnovert=zeros(Nboot)
    for kboot in range(0,Nboot):
        print 'kboot = '+str(kboot)
        #initalize
        for k in range(0,len(pgd_folders)):
            event=pgd_folders[k].split('/')[-1].split('.')[0]
            #Assign magnitudes
            if event=='ElMayor2010':
                Mcurrent=7.18
            if event=='Iquique2014':
                Mcurrent=8.19
            if event=='Maule2010':
                Mcurrent=8.85
            if event=='Mentawai2010':
                Mcurrent=7.68
            if event=='Napa2014':
                Mcurrent=6.11
            if event=='Parkfield2004':
                Mcurrent=5.92
            if event=='Tohoku2011':
                Mcurrent=9.09
            if event=='Tokachi2003':
                Mcurrent=8.25
            if event=='Nicoya2012':
                Mcurrent=7.57
            if event=='Aegean2014':
                Mcurrent=6.87
            pgd_temp=genfromtxt(path+'PGD/'+event+'.obs.pgd')
            ifilt_novert=where((pgd_temp[:,0]<maxdist) & (pgd_temp[:,1]>minpgd))[0]
            ifilt_all=where((pgd_temp[:,0]<maxdist) & (pgd_temp[:,2]>minpgd))[0]
            if k==0:
                pgd_novert=pgd_temp[ifilt_novert,:].copy()
                Mw_novert=ones(len(pgd_novert))*Mcurrent
                W_novert=ones(len(pgd_novert))/vecnorm(log10(pgd_novert[:,1]))
                pgd_all=pgd_temp[ifilt_all,:].copy()
                Mw_all=ones(len(pgd_all))*Mcurrent
                W_all=ones(len(pgd_all))/vecnorm(log10(pgd_all[:,2]))
            else:
                pgd_novert=r_[pgd_novert,pgd_temp[ifilt_novert]]
                Mw_novert=r_[Mw_novert,ones(len(pgd_temp[ifilt_novert,:]))*Mcurrent]
                W_novert=r_[W_novert,ones(len(pgd_novert[ifilt_novert,:]))/vecnorm(log10(pgd_novert[ifilt_novert,1]))]
                pgd_all=r_[pgd_all,pgd_temp[ifilt_all]]
                Mw_all=r_[Mw_all,ones(len(pgd_temp[ifilt_all,:]))*Mcurrent]
                W_all=r_[W_all,ones(len(pgd_all[ifilt_all,:]))/vecnorm(log10(pgd_all[ifilt_all,2]))]
        #Define regression quantities
        dnovert=log10(pgd_novert[:,1])
        dall=log10(pgd_all[:,2])
        iall=ones((len(dall),1))
        inovert=ones((len(dnovert),1))
        G_all=c_[iall,expand_dims(Mw_all,1)*iall,expand_dims(Mw_all*log10(pgd_all[:,0]),1)]
        G_novert=c_[inovert,expand_dims(Mw_novert,1)*inovert,expand_dims(Mw_novert*log10(pgd_novert[:,0]),1)]
        #Remove 10% of data
        remove=randint(0,len(dnovert),len(dnovert)/10)
        alli=arange(len(dnovert))
        keep=setxor1d(alli,remove)
        dnovert=dnovert[keep]
        dall=dall[keep]
        G_novert=G_novert[keep,:]
        G_all=G_all[keep,:]
        #make matrix of weights
        W_all=diag(W_all[keep])
        W_novert=diag(W_novert[keep])
        #Run regression
        P=matrix(W_novert.dot(G_novert))
        q=matrix(W_novert.dot(dnovert))
        coefficients_novert=array(l1(P,q))
        P=matrix(W_all.dot(G_all))
        q=matrix(W_all.dot(dall))
        coefficients_all=array(l1(P,q))
        Anovert[kboot]=coefficients_novert[0]
        Bnovert[kboot]=coefficients_novert[1]
        Cnovert[kboot]=coefficients_novert[2]
        Aall[kboot]=coefficients_all[0]
        Ball[kboot]=coefficients_all[1]
        Call[kboot]=coefficients_all[2]
        f=open(path+'coefficients/novert.boot','w')
        f.write('# mean , std dev (A,B,C)\n%.6f\t%.6f\n%.6f\t%.6f\n%.6f\t%.6f' %(Anovert.mean(),Anovert.std(),Bnovert.mean(),Bnovert.std(),Cnovert.mean(),Cnovert.std()))
        f.close()
        f=open(path+'coefficients/all.boot','w')
        f.write('# mean , std dev (A,B,C)\n%.6f\t%.6f\n%.6f\t%.6f\n%.6f\t%.6f' %(Aall.mean(),Aall.std(),Ball.mean(),Ball.std(),Call.mean(),Call.std()))
        f.close()
        print 'All coefficients '+str(coefficients_all)
        print 'Novert coefficients '+str(coefficients_novert)


#Get magnitude
if get_magnitude:
    #initalize
    for k in [1]:#range(0,len(pgd_folders)): #### pgd_folders_mag is for ommitance testing!!!!!
        event=pgd_folders[k].split('/')[-1].split('.')[0]
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
        for kt in range(len(tpgd)):
            #read all sites
            dall=[]
            dnov=[]
            R=[]
            print kt
            cleansta=genfromtxt(path+'misc/'+event+'.clean',usecols=0,dtype='S')
            cleanvalue=genfromtxt(path+'misc/'+event+'.clean',usecols=1)
            for ksta in range(len(stations)):
                smask=False
                sta=stations[ksta].split('/')[-1].split('.')[0]
                if event=='Tohoku2011' or event=='Tokachi2003':
                    sta=rjust(sta,4,'0')
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
            #coeff_novert=genfromtxt(path+'coefficients/'+suffix+'.novert.coeff')
            #coeff_all=genfromtxt(path+'coefficients/'+suffix+'.all.coeff')
            coeff_novert=genfromtxt(path+'coefficients/'+'novert.coeff')
            coeff_all=genfromtxt(path+'coefficients/'+'all.coeff')
            #Get coefficient  uncertainties
            sigma_novert=1.96*genfromtxt(path+'coefficients/novert.boot',usecols=1)
            sigma_novert[1]=-sigma_novert[1] ; sigma_novert[2]=-sigma_novert[2]
            sigma_all=1.96*genfromtxt(path+'coefficients/novert.boot',usecols=1)
            sigma_all[1]=-sigma_all[1] ; sigma_all[2]=-sigma_all[2]
            if len(dnov)>=minsta:
                print "using "+str(len(dnov))+' stations'
                Mw_nov[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert,weight=True,norm=2)
                Mw_all[kt]=scaling.PGD(dall*100,R,coefficients=coeff_all,weight=True,norm=2)
                #Now get upper and lower magnitude bound
                Mw_nov_plus[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert+sigma_novert,weight=True,norm=2)
                Mw_all_plus[kt]=scaling.PGD(dall*100,R,coefficients=coeff_all+sigma_all,weight=True,norm=2)
                Mw_nov_minus[kt]=scaling.PGD(dnov*100,R,coefficients=coeff_novert-sigma_novert,weight=True,norm=2)
                Mw_all_minus[kt]=scaling.PGD(dall*100,R,coefficients=coeff_all-sigma_all,weight=True,norm=2)
            else:
                Mw_nov[kt]=nan
                Mw_all[kt]=nan
        #Write to file
        savetxt(path+'PGD/magnitudes/'+event+'.'+suffix+'.mag',c_[tpgd,Mw_nov,Mw_all,Mw_nov_minus,Mw_nov_plus,Mw_all_minus,Mw_all_plus],fmt='%5.1f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f')


#Get magnitude error from the residuals
if magnitude_error==True:
    #Get coefficeints
    coeff_novert=genfromtxt(path+'coefficients/novert.coeff')
    coeff_all=genfromtxt(path+'coefficients/all.coeff')
    #Assemble all PGD's
    for k in range(0,len(pgd_folders)):
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        #Assign magnitudes
        if event=='ElMayor2010':
                Mcurrent=7.18
        if event=='Iquique2014':
            Mcurrent=8.19
        if event=='Maule2010':
            Mcurrent=8.85
        if event=='Mentawai2010':
                Mcurrent=7.68
        if event=='Napa2014':
            Mcurrent=6.11
        if event=='Parkfield2004':
            Mcurrent=5.92
        if event=='Tohoku2011':
            Mcurrent=9.09
        if event=='Tokachi2003':
            Mcurrent=8.25
        if event=='Nicoya2012':
            Mcurrent=7.57
        if event=='Aegean2014':
            Mcurrent=6.87
        pgd_temp=genfromtxt(path+'PGD/'+event+'.obs.pgd')
        ifilt_novert=where((pgd_temp[:,0]<maxdist) & (pgd_temp[:,1]>minpgd))[0]
        ifilt_all=where((pgd_temp[:,0]<maxdist) & (pgd_temp[:,2]>minpgd))[0]
        if k==0:
            pgd_novert=pgd_temp[ifilt_novert,:].copy()
            Mw_novert=ones(len(pgd_novert))*Mcurrent
            W_novert=ones(len(pgd_novert))/vecnorm(log10(pgd_novert[:,1]))
            pgd_all=pgd_temp[ifilt_all,:].copy()
            Mw_all=ones(len(pgd_all))*Mcurrent
            W_all=ones(len(pgd_all))/vecnorm(log10(pgd_all[:,2]))
        else:
            pgd_novert=r_[pgd_novert,pgd_temp[ifilt_novert]]
            Mw_novert=r_[Mw_novert,ones(len(pgd_temp[ifilt_novert,:]))*Mcurrent]
            W_novert=r_[W_novert,ones(len(pgd_novert[ifilt_novert,:]))/vecnorm(log10(pgd_novert[ifilt_novert,1]))]
            pgd_all=r_[pgd_all,pgd_temp[ifilt_all]]
            Mw_all=r_[Mw_all,ones(len(pgd_temp[ifilt_all,:]))*Mcurrent]
            W_all=r_[W_all,ones(len(pgd_all[ifilt_all,:]))/vecnorm(log10(pgd_all[ifilt_all,2]))]
    #Get pgd and distance
    dnovert=pgd_novert[:,1]
    dall=pgd_all[:,2]
    rnovert=pgd_novert[:,0]
    rall=pgd_all[:,0]
    #Forward calcualte magnitudes
    A=coeff_all[0]
    B=coeff_all[1]
    C=coeff_all[2]
    Mf_all=(log10(dall)-A)/(B+C*log10(rall))
    s=Mf_all-Mw_all
    print 'Standard deviation using all components is '+str(s.std())
    A=coeff_novert[0]
    B=coeff_novert[1]
    C=coeff_novert[2]
    Mf_novert=(log10(dnovert)-A)/(B+C*log10(rnovert))
    s=Mf_novert-Mw_novert
    print 'Standard deviation with no vertical is '+str(s.std())

#This plots waveforms and their PGD for QCing
if plot_sites_and_pgd:
    #initalize
    for k in range(len(earthquakes)):
        event=event=pgd_folders[k].split('/')[-1].split('.')[0]
        event_hypo=earthquakes[k].split('/')[-1].split('.')[0]
        hypo_time=UTCDateTime(genfromtxt(path+'event_info/'+event_hypo+'.hypo',dtype='S')[0])
        print event
        #Read station data
        stanames=genfromtxt(earthquakes[k],usecols=0,dtype='S')
        #Now loop over station data at each epoch
        event=pgd_folders[k].split('/')[-1]
        for ksta in range(len(stanames)):
            #Read waveforms
            print stanames[ksta]
            try:
                if event=='Tohoku2011':
                    n=read(pgd_folders[k]+'/'+rjust(stanames[ksta],4,'0')+'.LXN.sac')
                    e=read(pgd_folders[k]+'/'+rjust(stanames[ksta],4,'0')+'.LXE.sac')
                    u=read(pgd_folders[k]+'/'+rjust(stanames[ksta],4,'0')+'.LXZ.sac')
                else: 
                    n=read(pgd_folders[k]+'/'+stanames[ksta]+'.LXN.sac')
                    e=read(pgd_folders[k]+'/'+stanames[ksta]+'.LXE.sac')
                    u=read(pgd_folders[k]+'/'+stanames[ksta]+'.LXZ.sac')
                if n[0].stats.npts>0:
                    #Trim to times of interest
                    n.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    e.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    u.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                    #decimate if necessary (5Hz data)
                    if n[0].stats.delta==0.2:
                        n.decimate(5, no_filter=True)
                        e.decimate(5, no_filter=True)
                        u.decimate(5, no_filter=True)
                #Read PGDs
                pgd_values=genfromtxt(path+'PGD/'+event+'/'+stanames[ksta]+'.pgd')
                tnovert=pgd_values[-1,3]
                dnovert=pgd_values[-1,1]
                tall=pgd_values[-1,4]
                dall=pgd_values[-1,2]
                plt.close("all")
                plt.figure()
                plt.subplot(311)
                plt.title('Station '+stanames[ksta])
                plt.plot(n[0].times(),n[0].data)
                plt.plot(r_[tnovert,tnovert],r_[-dnovert,dnovert],c='r')
                plt.plot(r_[tall,tall],r_[-dall,dall],c='g')
                plt.legend(['','no vert','all'])
                plt.ylabel('North (m)')
                plt.subplot(312)
                plt.plot(e[0].times(),e[0].data)
                plt.plot(r_[tnovert,tnovert],r_[-dnovert,dnovert],c='r')
                plt.plot(r_[tall,tall],r_[-dall,dall],c='g')
                plt.ylabel('East (m)')
                plt.subplot(313)
                plt.plot(u[0].times(),u[0].data)
                plt.plot(r_[tnovert,tnovert],r_[-dnovert,dnovert],c='r')
                plt.plot(r_[tall,tall],r_[-dall,dall],c='g')
                plt.ylabel('Up (m)')
                plt.xlabel('Seconds after origin time')
                plt.savefig(path+'plots/waveforms/'+event+'/'+stanames[ksta]+'.png')
                
            except:
                print stanames[ksta]+' not found'
                        
# This makes a loglog scatter plot of all pgd's
if plot_all:
    #initalize
    #No verticals plot
    plt.figure()
    #THingy for legend
    #for k in [7,2,8,1,3,5,0,4,6]:
    for k in [9,4,1,10,3,5,7,2,0,6,8]:
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        print event
        pgd_novert=genfromtxt(path+'PGD/'+event+'.obs.pgd')
        i=where((pgd_novert[:,0]<maxdist) & (pgd_novert[:,1]>minpgd))[0]
        if event=='ElMayor2010':
            plt.loglog(1e-9,1e-9,'s',lw=0,c='#FFD700',markersize=4.5)
        if event=='Coquimbo2015':
            plt.loglog(1e-9,1e-9,'p',lw=0,c='#008000',markersize=6)
        if event=='Iquique2014':
            plt.loglog(1e-9,1e-9,'v',lw=0,c='#FF4500',markersize=6)
        if event=='Maule2010':
            plt.loglog(1e-9,1e-9,'o',lw=0,c='#1E90FF',markersize=6)
        if event=='Mentawai2010':
            plt.loglog(1e-9,1e-9,'D',lw=0,c='#32CD32',markersize=4.5)
        if event=='Napa2014':
            plt.loglog(1e-9,1e-9,'^',lw=0,c='#0000CD',markersize=6)
        if event=='Parkfield2004':
            plt.loglog(1e-9,1e-9,'p',lw=0,c='#DAA520',markersize=6)
        if event=='Tohoku2011':
            plt.loglog(1e-9,1e-9,'s',lw=0,c='#9ACD32',markersize=4.5)
        if event=='Tokachi2003':
            plt.loglog(1e-9,1e-9,'o',lw=0,c='#7B68EE',markersize=4.5)
        if event=='Nicoya2012':
            plt.loglog(1e-9,1e-9,'o',lw=0,c='#8B008B',markersize=6)
        if event=='Aegean2014':
            plt.loglog(1e-9,1e-9,'p',lw=0,c='#FF69B4',markersize=7)
    #end legend
    #get coefficients
    coeff=genfromtxt(path+'coefficients/novert.coeff')
    A=coeff[0]
    B=coeff[1]
    C=coeff[2]
    plt.suptitle('PGD scaling, no verticals')
    plt.title(r'log(PGD) = %.3f + %.3f$M_w$ %.3f$M_w$log(R)' %(A,B,C),fontsize=12)
    plt.legend([r'Tohoku-oki ($M_w9.09$)',r'Maule ($M_w8.85$)',r'Coquimbo ($M_w8.41$)',r'Tokachi-oki ($M_w8.25$)',r'Iquique ($M_w8.19$)',r'Mentawai ($M_w7.68$)',r'Nicoya ($M_w7.57$)',r'El Mayor ($M_w7.18$)',r'Aegean ($M_w6.87$)',r'Napa ($M_w6.11$)',r'Parkfield ($M_w5.92$)'],numpoints=1,loc=2,fontsize=10)
    plt.xlabel('Hypocentral Distance (km)')
    plt.ylabel('PGD (cm)')
    Mminor=arange(1,10,0.1)
    dist=10*logspace(-5,5)
    for kplot in range(len(Mminor)):
        plt.loglog(dist,10**(A+B*Mminor[kplot]+C*Mminor[kplot]*log10(dist)),color='#808080',lw=0.5)
    Mmajor=arange(1,11,1)
    dist=10*logspace(-5,5)
    for kplot in range(len(Mmajor)):
        plt.loglog(dist,10**(A+B*Mmajor[kplot]+C*Mmajor[kplot]*log10(dist)),color='k',lw=1)
    plt.xlim([8,1000])
    plt.ylim([1,1000])
    #Make magnitude annotations
    xy=(19.02, 5.42)
    plt.annotate(r"$M_w6$", xy=xy,xytext=xy,rotation=-23)
    xy=(60, 14.33)
    plt.annotate(r"$M_w7$", xy=xy,xytext=xy,rotation=-27)    
    xy=(47.04, 126.35)
    plt.annotate(r"$M_w8$", xy=xy,xytext=xy,rotation=-30)
    xy=(50.36, 777)
    plt.annotate(r"$M_w9$", xy=xy,xytext=xy,rotation=-32)
    #for k in [7,2,8,1,3,5,0,4,6]:
    for k in [9,4,1,10,3,5,7,2,0,6,8]:
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        print event
        pgd_novert=genfromtxt(path+'PGD/'+event+'.obs.pgd')
        i=where((pgd_novert[:,0]<maxdist) & (pgd_novert[:,1]>minpgd))[0]
        if event=='ElMayor2010':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'s',lw=0,c='#FFD700',markersize=4.5)
        if event=='Coquimbo2015':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'p',lw=0,c='#008000',markersize=6)
        if event=='Iquique2014':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'v',lw=0,c='#FF4500',markersize=6)
        if event=='Maule2010':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'o',lw=0,c='#1E90FF',markersize=6)
        if event=='Mentawai2010':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'D',lw=0,c='#32CD32',markersize=4.5)
        if event=='Napa2014':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'^',lw=0,c='#0000CD',markersize=6)
        if event=='Parkfield2004':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'p',lw=0,c='#DAA520',markersize=6)
        if event=='Tohoku2011':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'s',lw=0,c='#9ACD32',markersize=4.5)
        if event=='Tokachi2003':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'o',lw=0,c='#7B68EE',markersize=4.5)
        if event=='Nicoya2012':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'o',lw=0,c='#8B008B',markersize=6)
        if event=='Aegean2014':
            plt.loglog(pgd_novert[i,0],pgd_novert[i,1],'p',lw=0,c='#FF69B4',markersize=7)

    
    plt.figure()
    #THingy for legend
    #for k in [7,2,8,1,3,5,0,4,6]:
    for k in [8,3,9,2,4,6,1,0,5,7]:
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        print event
        pgd_novert=genfromtxt(path+'PGD/'+event+'.obs.pgd')
        i=where((pgd_novert[:,0]<maxdist) & (pgd_novert[:,1]>minpgd))[0]
        if event=='ElMayor2010':
            plt.loglog(1e-9,1e-9,'s',lw=0,c='#FFD700',markersize=4.5)
        if event=='Iquique2014':
            plt.loglog(1e-9,1e-9,'v',lw=0,c='#FF4500',markersize=6)
        if event=='Maule2010':
            plt.loglog(1e-9,1e-9,'o',lw=0,c='#1E90FF',markersize=6)
        if event=='Mentawai2010':
            plt.loglog(1e-9,1e-9,'D',lw=0,c='#32CD32',markersize=4.5)
        if event=='Napa2014':
            plt.loglog(1e-9,1e-9,'^',lw=0,c='#0000CD',markersize=6)
        if event=='Parkfield2004':
            plt.loglog(1e-9,1e-9,'p',lw=0,c='#DAA520',markersize=6)
        if event=='Tohoku2011':
            plt.loglog(1e-9,1e-9,'s',lw=0,c='#9ACD32',markersize=4.5)
        if event=='Tokachi2003':
            plt.loglog(1e-9,1e-9,'o',lw=0,c='#7B68EE',markersize=4.5)
        if event=='Nicoya2012':
            plt.loglog(1e-9,1e-9,'o',lw=0,c='#8B008B',markersize=6)
        if event=='Aegean2014':
            plt.loglog(1e-9,1e-9,'p',lw=0,c='#FF69B4',markersize=7)
    #end legend
    #get coefficients
    coeff=genfromtxt(path+'coefficients/all.coeff')
    A=coeff[0]
    B=coeff[1]
    C=coeff[2]
    plt.suptitle('PGD scaling, all components')
    plt.title(r'log(PGD) = %.3f + %.3f$M_w$ %.3f$M_w$log(R)' %(A,B,C),fontsize=12)
    plt.legend([r'Tohoku-oki ($M_w9.09$)',r'Maule ($M_w8.85$)',r'Tokachi-oki ($M_w8.25$)',r'Iquique ($M_w8.19$)',r'Mentawai ($M_w7.68$)',r'Nicoya ($M_w7.57$)',r'El Mayor ($M_w7.18$)',r'Aegean ($M_w6.87$)',r'Napa ($M_w6.11$)',r'Parkfield ($M_w5.92$)'],numpoints=1,loc=2,fontsize=10)
    plt.xlabel('Hypocentral Distance (km)')
    plt.ylabel('PGD (cm)')
    Mminor=arange(1,10,0.1)
    dist=10*logspace(-5,5)
    for kplot in range(len(Mminor)):
        plt.loglog(dist,10**(A+B*Mminor[kplot]+C*Mminor[kplot]*log10(dist)),color='#808080',lw=0.5)
    Mmajor=arange(1,11,1)
    dist=10*logspace(-5,5)
    for kplot in range(len(Mmajor)):
        plt.loglog(dist,10**(A+B*Mmajor[kplot]+C*Mmajor[kplot]*log10(dist)),color='k',lw=1)
    plt.xlim([8,1000])
    plt.ylim([1,1000])
    xy=(20.17, 6.4)
    plt.annotate(r"$M_w6$", xy=xy,xytext=xy,rotation=-23)
    xy=(46.87, 19.45)
    plt.annotate(r"$M_w7$", xy=xy,xytext=xy,rotation=-27)    
    xy=(46.70, 126.34)
    plt.annotate(r"$M_w8$", xy=xy,xytext=xy,rotation=-30)
    xy=(57.7, 673)
    plt.annotate(r"$M_w9$", xy=xy,xytext=xy,rotation=-32)
    #for k in [7,2,8,1,3,5,0,4,6]:
    for k in [8,3,9,2,4,6,1,0,5,7]:
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        pgd_all=genfromtxt(path+'PGD/'+event+'.obs.pgd')
        i=where((pgd_all[:,0]<maxdist) & (pgd_all[:,2]>minpgd))[0]
        if event=='ElMayor2010':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'s',lw=0,c='#FFD700',markersize=4.5)
        if event=='Iquique2014':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'v',lw=0,c='#FF4500',markersize=6)
        if event=='Maule2010':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'o',lw=0,c='#1E90FF',markersize=6)
        if event=='Mentawai2010':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'D',lw=0,c='#32CD32',markersize=4.5)
        if event=='Napa2014':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'^',lw=0,c='#0000CD',markersize=6)
        if event=='Parkfield2004':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'p',lw=0,c='#DAA520',markersize=6)
        if event=='Tohoku2011':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'s',lw=0,c='#9ACD32',markersize=4.5)
        if event=='Tokachi2003':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'o',lw=0,c='#7B68EE',markersize=4.5)
        if event=='Nicoya2012':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'o',lw=0,c='#8B008B',markersize=6)
        if event=='Aegean2014':
            plt.loglog(pgd_all[i,0],pgd_all[i,1],'p',lw=0,c='#FF69B4',markersize=7)
    
    plt.show()
                               
        
if plot_magnitudes:
    fig, axarr = plt.subplots(5,2) 
    ycount=0
    xcount=0
    #for k in [7,2,8,1,3,5,0,4,6]:
    for k in [8,3,9,2,4,6,1,0,5,7]:
        #Define which axes is beig used
        ax=axarr[ycount,xcount]
        ax2=ax.twinx()
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        #load magnitude data
        #M=genfromtxt(path+'PGD/magnitudes/swave3/'+event+'.mag')
        M2=genfromtxt(path+'PGD/magnitudes/2sta_baduncert/swave2/'+event+'.mag')
        M3=genfromtxt(path+'PGD/magnitudes/2sta_baduncert/swave3/'+event+'.mag')
        M4=genfromtxt(path+'PGD/magnitudes/2sta_baduncert/swave4/'+event+'.mag')
        #Make plot
        t=arange(0,1000)
        y=ones(len(t))
        if event=='ElMayor2010':
            #tscale
            tlims=[0,120]
            #do moment rate
            mr=genfromtxt(path+'STFs/elmayor.stf')
            scale=1e18
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*7.18,'--',lw=1.5,c='r')
            xy=(75,5.7)
            ax2.annotate(r"El Mayor $M_w7.18$", xy=xy,xytext=xy)
            ax2.set_ylim((5.5,7.5))
            ax2.set_ylabel(r'($\times 10^{18}$Nm/s)',labelpad=25)
        if event=='Iquique2014':
            #tscale
            tlims=[0,150]
            #do moment rate
            mr=genfromtxt(path+'STFs/iquique.stf')
            scale=1e26
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*8.19,'--',lw=1.5,c='r')
            xy=(96,6.4)
            ax2.annotate(r"Iquique $M_w8.19$", xy=xy,xytext=xy)
            ax2.set_ylim((6.0,8.5))
            ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=20)
        if event=='Maule2010':
            #tscale
            tlims=[0,250]
            #do moment rate
            mr=genfromtxt(path+'STFs/maule.stf')
            scale=1e27
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            #magnitude
            ax2.plot(t,y*8.85,'--',lw=1.5,c='r')
            xy=(160,6.5)
            ax2.annotate(r"Maule $M_w8.85$", xy=xy,xytext=xy)
            ax2.set_ylim((6,9.5))
            ax2.set_xlim(tlims)
            ax2.set_ylabel(r'($\times 10^{20}$Nm/s)',labelpad=25)
        if event=='Mentawai2010':
            #tscale
            tlims=[0,150]
            #do moment rate
            mr=genfromtxt(path+'STFs/mentawai.stf')
            scale=1e26
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*7.68,'--',lw=1.5,c='r')
            xy=(90,5.7)
            ax2.annotate(r"Mentawai $M_w7.68$", xy=xy,xytext=xy)
            ax2.set_ylim((5.5,8.0))
            ax.set_xlabel('Seconds after OT')
            ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=25)
        if event=='Napa2014':
            #tscale
            tlims=[0,40]
            #do moment rate
            mr=genfromtxt(path+'STFs/napa.stf')
            scale=1
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*6.11,'--',lw=1.5,c='r')
            xy=(25,5.25)
            ax2.annotate(r"Napa $M_w6.11$", xy=xy,xytext=xy)
            ax2.set_ylim((5.0,6.5))
            ax2.set_ylabel(r'($\times 10^{17}$Nm/s)',labelpad=25)
        if event=='Parkfield2004':
            #tscale
            tlims=[0,40]
            #do moment rate
            mr=genfromtxt(path+'STFs/parkfield.stf')
            scale=1
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*5.92,'--',lw=1.5,c='r')
            xy=(25,5.15)
            ax2.annotate(r"Parkfield $M_w5.92$", xy=xy,xytext=xy)
            ax2.set_ylim((5.0,6.0))
            ax2.set_ylabel(r'($\times 10^{17}$Nm/s)',labelpad=25)
            ax.set_xlabel('Seconds after OT')
        if event=='Tohoku2011':
            #Stuff for legend
            ax2.plot(-10,-10,lw=0,marker='s',color='#6495ED',markersize=3,label=r'$\beta=4km/s$')
            ax2.plot(-10,-10,lw=0,marker='s',color='#FF6347',markersize=3,label=r'$\beta=3km/s$')
            ax2.plot(-10,-10,lw=0,marker='s',color='#32CD32',markersize=3,label=r'$\beta=2km/s$')
            #tscale
            tlims=[0,250]
            #do moment rate
            mr=genfromtxt(path+'STFs/tohoku.stf')
            scale=1
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*9.09,'--',lw=1.5,c='r')
            xy=(150,5.9)
            ax2.annotate(r"Tohoku-oki $M_w9.09$", xy=xy,xytext=xy)
            ax2.set_ylim((5.5,9.8))
            ax2.set_ylabel(r'($\times 10^{20}$Nm/s)',labelpad=25)
            ax2.legend(bbox_to_anchor=(0., 1.10, 1., .102), loc=3,
                ncol=3, mode="expand", borderaxespad=0.,fontsize=12,numpoints=1)
        if event=='Tokachi2003':
            #tscale
            tlims=[0,150]
            #do moment rate
            mr=genfromtxt(path+'STFs/tokachi.stf')
            scale=1e26
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*8.25,'--',lw=1.5,c='r')
            xy=(90,5.8)
            ax2.annotate(r"Tokachi-oki $M_w8.25$", xy=xy,xytext=xy)
            ax2.set_ylim((5.5,8.5))
            ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=25)
        if event=='Nicoya2012':
            #tscale
            tlims=[0,80]
            #do moment rate
            mr=genfromtxt(path+'STFs/nicoya.stf')
            scale=1e26
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*7.57,'--',lw=1.5,c='r')
            xy=(50,6.3)
            ax2.annotate(r"Nicoya $M_w7.57$", xy=xy,xytext=xy)
            ax2.set_ylim((6.,8.0))
            ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=25)
        if event=='Aegean2014':
            #tscale
            tlims=[0,80]
            #do moment rate
            mr=genfromtxt(path+'STFs/aegean.stf')
            scale=1e18
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*6.87,'--',lw=1.5,c='r')
            xy=(50,5.3)
            ax2.annotate(r"Aegean $M_w6.87$", xy=xy,xytext=xy)
            ax2.set_ylim((5,7.5))
            ax2.set_ylabel(r'($\times 10^{18}$Nm/s)',labelpad=25)
        #Plot magnitudes
        #ax.plot(M[:,0],M[:,1],lw=0,marker='s',color='#6495ED',markersize=4)
        ax2.errorbar(M4[:,0],M4[:,2],yerr=[M4[:,2]-M4[:,5],M4[:,6]-M4[:,2]],fmt='o',color='#6495ED',ecolor='#6495ED',markersize=3,label=r'$\beta=4km/s$')
        ax2.errorbar(M3[:,0],M3[:,2],yerr=[M3[:,2]-M3[:,5],M3[:,6]-M3[:,2]],fmt='o',color='#FF6347',ecolor='#FF6347',markersize=3,label=r'$\beta=3km/s$')
        ax2.errorbar(M2[:,0],M2[:,2],yerr=[M2[:,2]-M2[:,5],M2[:,6]-M2[:,2]],fmt='o',color='#32CD32',ecolor='#32CD32',markersize=3,label=r'$\beta=2km/s$')
        #ax2.errorbar(M4[:,0],M4[:,1],yerr=[M4[:,1]-M4[:,3],M4[:,4]-M4[:,1]],fmt='o',color='#6495ED',ecolor='#6495ED',markersize=3,label=r'$\beta=4km/s$')
        #ax2.errorbar(M3[:,0],M3[:,1],yerr=[M3[:,1]-M3[:,3],M3[:,4]-M3[:,1]],fmt='o',color='#FF6347',ecolor='#FF6347',markersize=3,label=r'$\beta=3km/s$')
        #ax2.errorbar(M2[:,0],M2[:,1],yerr=[M2[:,1]-M2[:,3],M2[:,4]-M2[:,1]],fmt='o',color='#32CD32',ecolor='#32CD32',markersize=3,label=r'$\beta=2km/s$')
        ax2.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.set_ylabel(r'$M_w$',labelpad=26)
        #Adjust subplot counters
        ycount+=1
        if ycount==5:
            ycount=0
            xcount+=1
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.9, wspace=0.3, hspace=0.3)
    #Clear extraneous axes
    #fig.delaxes(axarr[4,1])
    plt.show()

#PGD move out plt for one event at a time
if plot_pgd_dist:
    for k in range(len(pgd_folders)):
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        print event
        #Get PGD times
        pgd_files=glob(path+'PGD/'+event+'/*.pgd')
        distance=zeros(len(pgd_files))
        time2pgd=zeros(len(pgd_files))
        cleansta=genfromtxt(path+'misc/'+event+'.clean',usecols=0,dtype='S')
        cleanvalue=genfromtxt(path+'misc/'+event+'.clean',usecols=1)
        for ksta in range(len(pgd_files)):
            sta=pgd_files[ksta].split('/')[-1].split('.')[0]
            if event=='Tohoku2011' or event=='Tokachi2003':
                sta=rjust(sta,4,'0')
            iclean=where((sta.lower()==cleansta) | (sta.upper()==cleansta))[0]
            if cleanvalue[iclean]==1: #Station is to be used
                pgd=genfromtxt(pgd_files[ksta])
                distance[ksta]=pgd[-1,5]
                time2pgd[ksta]=pgd[-1,3] #Using no vert for now
        #Make plot
        plt.figure()
        plt.plot(distance,time2pgd,lw=0,marker='o')
        plt.grid()
        plt.xlabel('Hypocentral distance (km)')
        plt.ylabel('Time to PGD (s)')
        
        #Plot reference lines
        for j in range(5):
            d=arange(0,maxdist+10,1)
            t=d/(1*j)
            plt.plot(d,t,'--',c='k')
        plt.xlim([0,maxdist])
        plt.ylim([0,time2pgd.max()])
        plt.title(event)
        plt.savefig(path+'plots/time2pgd/'+event+'.time2pgd.horiz.png')
    plt.close("all")
  
#Plot  PGD moveout for all events at the same time      
if plot_pgd_dist_all:
    plt.figure()
    plt.xlabel('Hypocentral distance (km)')
    plt.ylabel('Time to PGD (s)')
    #Legend stuff
    for k in [6,2,7,1,3,0,4,5]:
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        if event=='ElMayor2010':
            plt.plot(-100,-100,'s',lw=0,c='#FFD700',markersize=4.5)
        if event=='Coquimbo2015':
            plt.plot(-100,-100,'v',lw=0,c='#9ACD32',markersize=4.5)
        if event=='Iquique2014':
            plt.plot(-100,-1009,'v',lw=0,c='#FF4500',markersize=6)
        if event=='Maule2010':
            plt.plot(-100,-100,'o',lw=0,c='#1E90FF',markersize=6)
        if event=='Mentawai2010':
            plt.plot(-100,-100,'D',lw=0,c='#32CD32',markersize=4.5)
        if event=='Napa2014':
            plt.plot(1-100,-100,'^',lw=0,c='#0000CD',markersize=6)
        if event=='Parkfield2004':
            plt.plot(-100,-100,'p',lw=0,c='#DAA520',markersize=6)
        if event=='Tohoku2011':
            plt.plot(-100,-100,'s',lw=0,c='#9ACD32',markersize=4.5)
        if event=='Tokachi2003':
            plt.plot(-100,-100,'o',lw=0,c='#7B68EE',markersize=4.5)
    plt.legend(['Tohoku','Maule','Tokachi','Iquique','Mentawai','El Mayor','Napa','Parkfield'],loc=4)
    #Plot reference lines
    for j in range(5):
        d=arange(0,maxdist+10,1)
        t=d/(1*j)
        plt.plot(d,t,'--',c='k')
    plt.xlim([0,maxdist])
    plt.ylim([0,300])
    plt.title('All events')
    plt.grid()
    
    for k in [6,2,7,1,3,0,4,5]:
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        print event
        #Get PGD times
        pgd_files=glob(path+'PGD/'+event+'/*.pgd')
        distance=zeros(len(pgd_files))
        time2pgd=zeros(len(pgd_files))
        cleansta=genfromtxt(path+'misc/'+event+'.clean',usecols=0,dtype='S')
        cleanvalue=genfromtxt(path+'misc/'+event+'.clean',usecols=1)
        for ksta in range(len(pgd_files)):
            sta=pgd_files[ksta].split('/')[-1].split('.')[0]
            if event=='Tohoku2011' or event=='Tokachi2003':
                sta=rjust(sta,4,'0')
            iclean=where((sta.lower()==cleansta) | (sta.upper()==cleansta))[0]
            if cleanvalue[iclean]==1: #Station is to be used
                pgd=genfromtxt(pgd_files[ksta])
                distance[ksta]=pgd[-1,5]
                time2pgd[ksta]=pgd[-1,3] #Using no vert for now
        #Make plot
        if event=='ElMayor2010':
            plt.plot(distance,time2pgd,'s',lw=0,c='#FFD700',markersize=4.5)
        if event=='Iquique2014':
            plt.plot(distance,time2pgd,'v',lw=0,c='#FF4500',markersize=6)
        if event=='Maule2010':
            plt.plot(distance,time2pgd,'o',lw=0,c='#1E90FF',markersize=6)
        if event=='Mentawai2010':
            plt.plot(distance,time2pgd,'D',lw=0,c='#32CD32',markersize=4.5)
        if event=='Napa2014':
            plt.plot(distance,time2pgd,'^',lw=0,c='#0000CD',markersize=6)
        if event=='Parkfield2004':
            plt.plot(distance,time2pgd,'p',lw=0,c='#DAA520',markersize=6)
        if event=='Tohoku2011':
            plt.plot(distance,time2pgd,'s',lw=0,c='#9ACD32',markersize=4.5)
        if event=='Tokachi2003':
            plt.plot(distance,time2pgd,'o',lw=0,c='#7B68EE',markersize=4.5)
        #Make mmoveout annotations
        xy=(490, 143)
        plt.annotate(r"$4km/s$", xy=xy,xytext=xy,rotation=33)
        xy=(440, 173)
        plt.annotate(r"$3km/s$", xy=xy,xytext=xy,rotation=41)
        xy=(395, 230)
        plt.annotate(r"$2km/s$", xy=xy,xytext=xy,rotation=51)
        xy=(175, 225)
        plt.annotate(r"$1km/s$", xy=xy,xytext=xy,rotation=66)
        plt.savefig(path+'plots/time2pgd/all_time2pgd.horiz.png')
    plt.close("all")
        
        
        
#Inspect all waveforms and cleanup
if clean:
    #initalize
    for k in [1]:#range(len(earthquakes)):
        #Now loop over station data at each epoch
        event=pgd_folders[k].split('/')[-1]
        #Read station data
        stanames=glob(path+'GPS/sac/'+event+'/*LXN.sac')
        #Open cleanup tracking file
        f=open(path+'misc/'+event+'.clean','w')
        f.write('#Station name, 1 is keep 0 is ignore\n')
        print event
        event_hypo=earthquakes[k].split('/')[-1].split('.')[0]
        hypo_time=UTCDateTime(genfromtxt(path+'event_info/'+event_hypo+'.hypo',dtype='S')[0])
        for ksta in range(len(stanames)):
            #Read waveforms
            print 'Station '+str(ksta)+' of '+str(len(stanames))
            station=stanames[ksta].split('/')[-1].split('.')[0]
            if event=='Tohoku2011':
                n=read(pgd_folders[k]+'/'+rjust(station,4,'0')+'.LXN.sac')
                e=read(pgd_folders[k]+'/'+rjust(station,4,'0')+'.LXE.sac')
                u=read(pgd_folders[k]+'/'+rjust(station,4,'0')+'.LXZ.sac')
            else: 
                n=read(pgd_folders[k]+'/'+station+'.LXN.sac')
                e=read(pgd_folders[k]+'/'+station+'.LXE.sac')
                u=read(pgd_folders[k]+'/'+station+'.LXZ.sac')
            if n[0].stats.npts>0:
                #Trim to times of interest
                n.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                e.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                u.trim(starttime=hypo_time,endtime=hypo_time+tmax)
                #decimate if necessary (5Hz data)
                if n[0].stats.delta==0.2:
                    n.decimate(5, no_filter=True)
                    e.decimate(5, no_filter=True)
                    u.decimate(5, no_filter=True)
                plt.close("all")
                plt.figure(1)
                plt.subplot(311)
                plt.title('Station '+stanames[ksta])
                plt.plot(n[0].times(),n[0].data)
                plt.ylabel('North (m)')
                plt.subplot(312)
                plt.plot(e[0].times(),e[0].data)
                plt.ylabel('East (m)')
                plt.subplot(313)
                plt.plot(u[0].times(),u[0].data)
                plt.ylabel('Up (m)')
                plt.xlabel('Seconds after origin time')
                plt.show()
                print 'Keep station (1/0)?'
                keep = raw_input()
                f.write(station+'\t'+keep+'\n')
        f.close()        
        
        
        
            