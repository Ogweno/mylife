from mudpy import forward

rupt=u'/Users/dmelgar/FakeQuakes/bbp_test/output/ruptures/bbp.000003.rupt'
log_file=u'/Users/dmelgar/FakeQuakes/bbp_test/output/ruptures/bbp.000003.log'
outfile=u'/Users/dmelgar/FakeQuakes/bbp_test/output/ruptures/bbp.000003.srf'

forward.mudpy2srf(rupt,log_file)