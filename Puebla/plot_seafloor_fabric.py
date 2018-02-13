from numpy import where,isnan,histogram,genfromtxt,deg2rad,pi
import matplotlib.cm as cm
from matplotlib.pyplot import figure, show, rc


f=genfromtxt('/Users/dmelgar/code/GMT/Puebla/strike.xyz')

# force square figure and square axes looks better for polar, IMO
fig = figure(figsize=(7,7))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

N = 20
i=where(isnan(f[:,2])==False)[0]
radii,theta= histogram(f[i,2],40)
theta=deg2rad(theta)
theta=theta[0:-1]
bars = ax.bar(theta, radii, bottom=0.0,width=0.15,color='#DAA520')

ax.bar(deg2rad(299),radii.max()*1.2,width=0.05,color='#00008B',alpha=0.5)
ax.bar(deg2rad(299-180),radii.max()*1.05,width=0.05,color='#00008B',alpha=0.5)
ax.set_ylim([0,16000])
ax.legend(['Seafloor lineaments','Fault strike'],bbox_to_anchor=(0.47, 1.15))
ax.set_theta_offset(0.5*pi)
ax.set_theta_direction('clockwise')




#for r,bar in zip(radii, bars):
#    bar.set_facecolor( cm.jet(r/10.))
#    bar.set_alpha(0.5)

show()