from mudpy import view
from matplotlib import pyplot as plt

view.plot_grd(u'/Users/dmelgar/Slab_Models/sam_slab1.0_strclip.grd',[-20,20],plt.cm.jet)
#CSN hypocenter
plt.scatter(360-74.391,-43.517,c='w',s=50)
plt.ylim([-46,-37])
plt.xlim([280,295])
plt.title('Strike')

view.plot_grd(u'/Users/dmelgar/Slab_Models/sam_slab1.0_dipclip.grd',[-40,0],plt.cm.jet)
#CSN hypocenter
plt.scatter(360-74.391,-43.517,c='w',s=50)
plt.ylim([-46,-37])
plt.xlim([280,295])
plt.title('Dip')

plt.show()