from mudpy import gmttools
from string import rjust

path='/Users/dmelgar/Slip_inv/Amatrice_3Dfitsgeol_final1/output/inverse_models/models/'
iter=200
for k in range(iter):
    mod=rjust(str(k),4,'0')
    gmttools.make_total_model(u'/Users/dmelgar/Slip_inv/Amatrice_3Dfitsgeol_final1/output/inverse_models/models/jacknife.'+mod+'.inv',0)