from glob import glob
from os import remove


f=glob('/Volumes/Illapel/FQ/illapel/GFs/dynamic/*sub*')

for k in range(len(f)):
    print f[k]
    f2=glob(f[k]+'/*.[e,n,z]')
    for j in range(len(f2)):
        remove(f2[j])