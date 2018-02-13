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
plot_magnitudes=True


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


if plot_magnitudes:
    fig, axarr = plt.subplots(5,2) 
    ycount=0
    xcount=0
    #for k in [7,2,8,1,3,5,0,4,6]:
    for k in [8,3,9,2,4,6,1,0,5,7]:
        #Define which axes is beig used
        event=pgd_folders[k].split('/')[-1].split('.')[0]
        #load magnitude data
        #M=genfromtxt(path+'PGD/magnitudes/swave3/'+event+'.mag')
        #M2=genfromtxt(path+'PGD/magnitudes/2sta_baduncert/swave2/'+event+'.mag')
        #M3=genfromtxt(path+'PGD/magnitudes/2sta_baduncert/swave3/'+event+'.mag')
        #M4=genfromtxt(path+'#PGD/magnitudes/2sta_baduncert/swave4/'+event+'.mag')
        #Make plot
        t=arange(0,1000)
        y=ones(len(t))
        print event
        if event=='Iquique2014':
            ax=axarr[1,0]
            ax2=ax.twinx()
            #tscale
            tlims=[0,300]
            #do moment rate
            mr=genfromtxt(path+'STFs/iquique.stf')
            scale=1e26
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*8.19,'--',lw=1.5,c='r')
            xy=(150,6.4)
            ax2.annotate(r"Iquique $M_w8.19$", xy=xy,xytext=xy)
            ax2.set_ylim((6.0,9.0))
            ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=20)
            ax2.yaxis.tick_left()
            ax.yaxis.tick_right()
            ax.set_ylabel(r'$M_w$',labelpad=26)
        if event=='Coquimbo2015':
            ax=axarr[0,0]
            ax2=ax.twinx()
            #tscale
            tlims=[0,300]
            #do moment rate
            mr=genfromtxt(path+'STFs/illapel.stf')
            scale=1e26
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*8.19,'--',lw=1.5,c='r')
            xy=(150,6.4)
            ax2.annotate(r"Illapel $M_w8.19$", xy=xy,xytext=xy)
            ax2.set_ylim((6.0,9.0))
            ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=20)
            ax2.yaxis.tick_left()
            ax.yaxis.tick_right()
            ax.set_ylabel(r'$M_w$',labelpad=26)
        if event=='Nicoya2012':
            ax=axarr[2,0]
            ax2=ax.twinx()
            #tscale
            tlims=[0,300]
            #do moment rate
            mr=genfromtxt(path+'STFs/iquique_after.stf')
            scale=1e26
            ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
            ax.set_xlim(tlims)
            ax2.plot(t,y*7.54,'--',lw=1.5,c='r')
            xy=(130,6.4)
            ax2.annotate(r"Iquique after $M_w7.54$", xy=xy,xytext=xy)
            ax2.set_ylim((6,9.0))
            ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=20)
            ax2.yaxis.tick_left()
            ax.yaxis.tick_right()
            ax.set_ylabel(r'$M_w$',labelpad=26)
        #if event=='Maule2010':
        #    #tscale
        #    tlims=[0,250]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/maule.stf')
        #    scale=1e27
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    #magnitude
        #    ax2.plot(t,y*8.85,'--',lw=1.5,c='r')
        #    xy=(160,6.5)
        #    ax2.annotate(r"Maule $M_w8.85$", xy=xy,xytext=xy)
        #    ax2.set_ylim((6,9.5))
        #    ax2.set_xlim(tlims)
        #    ax2.set_ylabel(r'($\times 10^{20}$Nm/s)',labelpad=25)
        #if event=='Mentawai2010':
        #    #tscale
        #    tlims=[0,150]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/mentawai.stf')
        #    scale=1e26
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    ax2.plot(t,y*7.68,'--',lw=1.5,c='r')
        #    xy=(90,5.7)
        #    ax2.annotate(r"Mentawai $M_w7.68$", xy=xy,xytext=xy)
        #    ax2.set_ylim((5.5,8.0))
        #    ax.set_xlabel('Seconds after OT')
        #    ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=25)
        #if event=='Napa2014':
        #    #tscale
        #    tlims=[0,40]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/napa.stf')
        #    scale=1
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    ax2.plot(t,y*6.11,'--',lw=1.5,c='r')
        #    xy=(25,5.25)
        #    ax2.annotate(r"Napa $M_w6.11$", xy=xy,xytext=xy)
        #    ax2.set_ylim((5.0,6.5))
        #    ax2.set_ylabel(r'($\times 10^{17}$Nm/s)',labelpad=25)
        #if event=='Parkfield2004':
        #    #tscale
        #    tlims=[0,40]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/parkfield.stf')
        #    scale=1
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    ax2.plot(t,y*5.92,'--',lw=1.5,c='r')
        #    xy=(25,5.15)
        #    ax2.annotate(r"Parkfield $M_w5.92$", xy=xy,xytext=xy)
        #    ax2.set_ylim((5.0,6.0))
        #    ax2.set_ylabel(r'($\times 10^{17}$Nm/s)',labelpad=25)
        #    ax.set_xlabel('Seconds after OT')
        #if event=='Tohoku2011':
        #    #Stuff for legend
        #    ax2.plot(-10,-10,lw=0,marker='s',color='#6495ED',markersize=3,label=r'$\beta=4km/s$')
        #    ax2.plot(-10,-10,lw=0,marker='s',color='#FF6347',markersize=3,label=r'$\beta=3km/s$')
        #    ax2.plot(-10,-10,lw=0,marker='s',color='#32CD32',markersize=3,label=r'$\beta=2km/s$')
        #    #tscale
        #    tlims=[0,250]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/tohoku.stf')
        #    scale=1
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    ax2.plot(t,y*9.09,'--',lw=1.5,c='r')
        #    xy=(150,5.9)
        #    ax2.annotate(r"Tohoku-oki $M_w9.09$", xy=xy,xytext=xy)
        #    ax2.set_ylim((5.5,9.8))
        #    ax2.set_ylabel(r'($\times 10^{20}$Nm/s)',labelpad=25)
        #    ax2.legend(bbox_to_anchor=(0., 1.10, 1., .102), loc=3,
        #        ncol=3, mode="expand", borderaxespad=0.,fontsize=12,numpoints=1)
        #if event=='Tokachi2003':
        #    #tscale
        #    tlims=[0,150]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/tokachi.stf')
        #    scale=1e26
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    ax2.plot(t,y*8.25,'--',lw=1.5,c='r')
        #    xy=(90,5.8)
        #    ax2.annotate(r"Tokachi-oki $M_w8.25$", xy=xy,xytext=xy)
        #    ax2.set_ylim((5.5,8.5))
        #    ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=25)
        #if event=='Nicoya2012':
        #    #tscale
        #    tlims=[0,80]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/nicoya.stf')
        #    scale=1e26
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    ax2.plot(t,y*7.57,'--',lw=1.5,c='r')
        #    xy=(50,6.3)
        #    ax2.annotate(r"Nicoya $M_w7.57$", xy=xy,xytext=xy)
        #    ax2.set_ylim((6.,8.0))
        #    ax2.set_ylabel(r'($\times 10^{19}$Nm/s)',labelpad=25)
        #if event=='Aegean2014':
        #    #tscale
        #    tlims=[0,80]
        #    #do moment rate
        #    mr=genfromtxt(path+'STFs/aegean.stf')
        #    scale=1e18
        #    ax.fill_between(mr[:,0],0,mr[:,1]/scale,facecolor='#FFE4E1')
        #    ax.set_xlim(tlims)
        #    ax2.plot(t,y*6.87,'--',lw=1.5,c='r')
        #    xy=(50,5.3)
        #    ax2.annotate(r"Aegean $M_w6.87$", xy=xy,xytext=xy)
        #    ax2.set_ylim((5,7.5))
        #    ax2.set_ylabel(r'($\times 10^{18}$Nm/s)',labelpad=25)
        #Plot magnitudes
        #ax.plot(M[:,0],M[:,1],lw=0,marker='s',color='#6495ED',markersize=4)
        #ax2.errorbar(M4[:,0],M4[:,2],yerr=[M4[:,2]-M4[:,5],M4[:,6]-M4[:,2]],fmt='o',color='#6495ED',ecolor='#6495ED',markersize=3,label=r'$\beta=4km/s$')
        #ax2.errorbar(M3[:,0],M3[:,2],yerr=[M3[:,2]-M3[:,5],M3[:,6]-M3[:,2]],fmt='o',color='#FF6347',ecolor='#FF6347',markersize=3,label=r'$\beta=3km/s$')
        #ax2.errorbar(M2[:,0],M2[:,2],yerr=[M2[:,2]-M2[:,5],M2[:,6]-M2[:,2]],fmt='o',color='#32CD32',ecolor='#32CD32',markersize=3,label=r'$\beta=2km/s$')
        #ax2.errorbar(M4[:,0],M4[:,1],yerr=[M4[:,1]-M4[:,3],M4[:,4]-M4[:,1]],fmt='o',color='#6495ED',ecolor='#6495ED',markersize=3,label=r'$\beta=4km/s$')
        #ax2.errorbar(M3[:,0],M3[:,1],yerr=[M3[:,1]-M3[:,3],M3[:,4]-M3[:,1]],fmt='o',color='#FF6347',ecolor='#FF6347',markersize=3,label=r'$\beta=3km/s$')
        #ax2.errorbar(M2[:,0],M2[:,1],yerr=[M2[:,1]-M2[:,3],M2[:,4]-M2[:,1]],fmt='o',color='#32CD32',ecolor='#32CD32',markersize=3,label=r'$\beta=2km/s$')
        #Adjust subplot counters
        ycount+=1
        if ycount==5:
            ycount=0
            xcount+=1
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.9, wspace=0.3, hspace=0.3)
    #Clear extraneous axes
    #fig.delaxes(axarr[4,1])
    plt.show()