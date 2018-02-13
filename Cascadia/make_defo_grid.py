from mudpy import forward
from numpy import genfromtxt,where
import matplotlib.path as mplPath
from string import rjust

path=genfromtxt('/Users/dmelgar/Cascadia/stations/deformation/deformation_area.txt')
forward.make_grid(-130,-120,40,52,0.1,0.1,u'/Users/dmelgar/Cascadia/stations/deformation/deformation_grid.txt')
grid=genfromtxt(u'/Users/dmelgar/Cascadia/stations/deformation/deformation_grid.txt')

#Make path
path = mplPath.Path(path)

#filter
in_path=path.contains_points(grid[:,1:])
i=where(in_path==True)[0]

grid=grid[i,1:]

#Save
f=open(u'/Users/dmelgar/Cascadia/stations/deformation/deformation.sta','w')
for k in range(len(i)):
    line='CA%s\t%.4f\t%.4f\n' %(rjust(str(k+1),4,'0'),grid[k,0],grid[k,1])
    f.write(line)
f.close()