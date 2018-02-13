from scipy.interpolate import LSQBivariateSpline
from numpy import linspace,genfromtxt,where,arange,array,zeros,ones
import matplotlib.path as mplPath
from scipy.optimize import minimize

decimate=500
g=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Oct26-30_combined.txt')
sta=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Oct26-30_combined.txt',usecols=0,dtype='S')
asc=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M5.8_M6.6/T117_Italy/T117_Italy.lltnde')
desc=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M5.8_M6.6/T22_Italy/T22_Italy.lltnde')
overlap_poly=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/T22_T117_overlap.txt')

#Path for overalp are of all interferos
bbPath = mplPath.Path(overlap_poly)

#Decimate descending
i=where(bbPath.contains_points(desc[:,0:2])==True)[0]
j=arange(0,len(i),decimate)
i=i[j]
Ndesc=len(i)
desc=desc[i,:]

x=desc[:,0]
y=desc[:,1]
f=desc[:,6]
look=desc[:,3:6].ravel()
u_test=ones(look.T.shape)

G=zeros((Ndesc,3*Ndesc))
for k in range(Ndesc):
    if k==0:
        print 'Working on pixel No:'
    if k%100==0:
        print '..'+str(k)
    
    lookE_desc=desc[k,3]
    lookN_desc=desc[k,4]
    lookU_desc=desc[k,5]
    
    #Form inversion quantities
    Gtemp=array([[lookE_desc,lookN_desc,lookU_desc]])
    G[k,3*k:3*k+3]=Gtemp

xknots = linspace(13, 13.5, 10)
yknots = linspace(42.6, 43.1, 10)

#spline=LSQBivariateSpline(x, y, f, xknots, yknots,w=None)
#X,Y=meshgrid(linspace(13,13.5),linspace(42.6,43.1))
#Z=spline.ev(X,Y)
#plt.figure()
#pcolormesh(X,Y,Z,vmin=-300,vmax=300)
#colorbar()



#Minimization test

def func(u):
    from numpy import meshgrid,linspace
    from scipy.linalg import norm
    
    coords=meshgrid(linspace(13,13.5),linspace(42.6,43.1))
    xknots = linspace(13, 13.5, 10)
    yknots = linspace(42.6, 43.1, 10)
    spline_fit=LSQBivariateSpline(x, y, G.dot(u), xknots, yknots,w=None).ev(coords[0],coords[1])
    spline_fit=spline_fit.ravel()
    Z=norm(spline_fit-f)
    return Z




#deg=2
#vander=polynomial.polyvander2d(x,y,[deg,deg])
#c = np.linalg.lstsq(vander, f)[0]
#X,Y=meshgrid(linspace(13,13.5),linspace(42.6,43.1))
#Z=polynomial.polyval2d(X, Y, c)
#pcolormesh(X,Y,Z,vmin=-300,vmax=300)
#scatter(x,y,f,lw=0,vmin=-300,vmax=300,lw=0)