from pyrocko import cake
from numpy import arange,array

km = 1000.

# Load builtin 'prem-no-ocean' model ('.m': medium resolution variant)
#model = cake.load_model('prem-no-ocean.m')
model = cake.load_model('/Users/dmelgar/code/pyrocko/build/lib.macosx-10.6-x86_64-2.7/pyrocko/data/earthmodels/Nocal.nd')

# Source depth [m].
#source_depth = 300. * km
source_depth= 5.8

# Distances as a numpy array [deg].
distances = array([5,10,20,50,100])*km * cake.m2d

# Define the phase to use.
Phase = cake.PhaseDef('P')

# calculate distances and arrivals and print them:
print 'distance [km]      time [s]'
for arrival in model.arrivals(distances, phases=Phase, zstart=source_depth):
    print '%13g %13g' % (arrival.x*cake.d2m/km, arrival.t)