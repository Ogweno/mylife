import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import where,genfromtxt,squeeze

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

afin=genfromtxt('/Users/dmelgar/Chiapas2017/afters/afters_3days_clean.txt')
af=afin.squeeze()
i=where(af[:,4]>-97)[0]
af=af[i,:]
ax.scatter(af[:,4],af[:,3],-af[:,5])

plt.show()