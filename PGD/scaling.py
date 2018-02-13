'''
Scaling regression tools
Original by B. Crowell U.Washington
Modified by D. Melgar U.C. Berkeley
'''

from numpy import exp,log10,power,expand_dims,ones,where,r_,array
from numpy.linalg import lstsq
from scipy.linalg import norm as vecnorm
from l1 import l1
from cvxopt import matrix

def PGD(pgd, r,coefficients,weight=True,norm=2,residual=False):
    '''
    Run PGD regression
    
    IN:
        pgd is pgd data
        r are station/event distances
        coefficients is a 3 vector with A,B and C
        weight is boolean to apply distance weighting
        norm switches between L1 and L2 norm minimizing solver
    '''
    mintol=1e-5 #This value is used isntad of zero
    #Define regression coefficients
    A=coefficients[0]
    B=coefficients[1]
    C=coefficients[2] 
    #Decide whether to apply distance weights
    if weight:
        #W = expand_dims(exp(-power(1./log10(pgd),2)/2/power(2*min(1./log10(pgd)),2)),1)
        W = expand_dims(exp(-power(r,2)/2/power(2*min(r),2)),1)
        #W = expand_dims(exp(-power(r,3)/2/power(2*min(r),3)),1)
    else:
        W=expand_dims(ones(len(r)),1)
    # Green's functions
    G = expand_dims(B+C*(log10(r)),1)
    # Clean data for small values
    i=where(pgd<mintol)[0]
    pgd[i]=mintol
    #Define data vector
    b = expand_dims(log10(pgd)-A,1)
    # L1 or L2 norm minimizing solver
    if norm==2:
        P=W*G
        q=W*b
        M = lstsq(P,q)[0]
        R=P.dot(M)-q
        res = vecnorm(R)
    elif norm==1:
        P=matrix(W*G)
        q=matrix(W*b)
        print P
        print q
        M=l1(P,q)[0]
        #get residuals
        res=sum(abs(P*matrix(M)-q))
    if residual:
        return M,res
    else:
        return M


def PD(d, r):
	A = -0.811
	B = 0.546
	C = -1.708
	M = (math.log10(d)-C*(math.log10(r))-A)/B
	return(M)


