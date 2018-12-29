#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='run'),
    InitialCondition(flameType = 'premixed',pressure = 101325,fuel='CH4:1.0',oxidizer='O2:1, N2:3.76',
    nPoints = 2000,Tu = 300,equivalenceRatio=1.0,xLeft=0.0,xRight=0.01),
    StrainParameters(initial=0,final=0),
    General(fixedLeftLocation = True,continuityBC='fixedLeft',flameGeometry='planar',twinFlame=True,nThreads=5),
    Chemistry(transportModel='Approx', kineticsModel='standard',threshold=1e-05),    
    Grid(gridMax=5e-06,gridMin=5e-06),
    PositionControl(proportionalGain=2000, xInitial=0.007, xFinal=0.007),
    TerminationCondition(tEnd = 0.008),
    Times(profileStepInterval=50,globalTimestep = 2e-09,currentStateStepInterval = 1),
    Debug(adaptation = False,regridding = False))

if __name__ == '__main__':
   
    conf.run()

