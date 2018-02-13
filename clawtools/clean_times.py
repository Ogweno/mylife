import sys
from numpy import genfromtxt,where,savetxt,nan

fgmax_file=sys.argv[1]
fg=genfromtxt(fgmax_file)
i=where(fg[:,-1]<0)[0]
fg[i,-1]=nan
savetxt(fgmax_file,fg,fmt='%.12e\t%.12e\t%d\t%.12e\t%.12e\t%.12e')