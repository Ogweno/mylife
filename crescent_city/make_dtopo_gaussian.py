from numpy import arange,exp,meshgrid,c_,zeros,ones,r_,savetxt

fout='/Users/dmelgar/tidegauge_noise/data/cres/_dtopo/gauss_0.2deg_south.dtopo'

width=0.2 #in degrees
sigma=width/3. #99.7% under the curve
max_amplitude=1.0

lon_center=-126.0
lat_center=39.5

# one degree max
lon=arange(-1,1,0.005)
lat=arange(-1,1,0.005)
lat=lat[::-1] #for proper output format


# make actual grid
lon=lon+lon_center
lat=lat+lat_center
X,Y=meshgrid(lon,lat)
lon=X.ravel()
lat=Y.ravel()

#Make gaussian
amplitude=max_amplitude*exp(-((lon-lon_center)**2)/(2*sigma**2)-((lat-lat_center)**2)/(2*sigma**2))

out1=c_[zeros(len(lon)),lon,lat,zeros(len(lon))]
out2=c_[ones(len(lon)),lon,lat,amplitude]
out=r_[out1,out2]

savetxt(fout,out,fmt='%.1f\t%.3f\t%.3f\t%.3f')