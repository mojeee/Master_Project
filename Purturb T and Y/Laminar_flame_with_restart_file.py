#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='run'),
    InitialCondition(restartFile='prof000009.h5'),
    StrainParameters(initial=0,final=0),
    General(fixedLeftLocation=True,continuityBC='fixedLeft',flameGeometry='planar',twinFlame=True,nThreads=3),
    Chemistry(transportModel='Approx', kineticsModel='standard',threshold=1e-05),    
    Grid(vtol=0.1, dvtol=0.15, gridMin=5e-6, gridMax=0.001),
    PositionControl(proportionalGain=2000, xInitial=0.005, xFinal=0.005),
    TerminationCondition(tEnd = 0.0005),
    Times(profileStepInterval=50,globalTimestep = 2e-09,currentStateStepInterval = 1),
    Debug(adaptation = False,regridding = False))

if __name__ == '__main__':
   
    conf.run()

