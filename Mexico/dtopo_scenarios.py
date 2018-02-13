from mudpy import forward
from numpy import arange


x=arange(-106,-95,0.05)+360
y=arange(14,20,0.05)

#forward.move_seafloor_okada(u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000046.rupt',u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000046.dtopo',x,y)
#forward.move_seafloor_okada(u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000057.rupt',u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000057.dtopo',x,y)
#forward.move_seafloor_okada(u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000062.rupt',u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000062.dtopo',x,y)
forward.move_seafloor_okada(u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000042.rupt',u'/Users/dmelgar/FakeQuakes/Mexico_scenarios/output/ruptures/all_CA.000042.dtopo',x,y)
