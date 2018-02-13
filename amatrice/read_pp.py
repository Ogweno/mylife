from numpy import r_,array

gps_file='/Users/dmelgar/Downloads/ppp_5/pos_2016303.txt'
sta='LNSS'

f=open(gps_file,'r')
kpos=0
while True:
    line=f.readline()
    if sta in line:
        line=f.readline()
        if kpos==0:
            x=array([line.split()[1]])
            y=array([line.split()[2]])
            z=array([line.split()[3]])
            kpos+=1
        else:
            x=r_[x,array([line.split()[1]])]
            y=r_[x,array([line.split()[2]])]
            z=r_[x,array([line.split()[3]])]
    if line == '':
        break
        