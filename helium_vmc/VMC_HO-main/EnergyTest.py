import numpy as np
from VMCHelp import *
from VMC_HO import *
 
import stats as CalcStatistics
 
HO=HOClass()
alpha=[0.33]
HO.SetParams(alpha)
R=np.zeros( (1,1) )
R[0,0]=0.5
HO.SetCoords(R)
local_energy=alpha[0] + ( 0.5 - 2.*alpha[0]**2 )*R[0,0]**2
print ("The Local Energy should be",local_energy,HO.LocalEnergy(R))
np.random.seed()
energyList=VMC(HO,1000)
print ("The correct answer is approximately 0.543")
print ("A short run gives ",CalcStatistics.Stats(np.array(energyList)))
energyList=VMC(HO,100000)
print ("A longer run gives ",CalcStatistics.Stats(np.array(energyList)))
plt.plot(energyList)
plt.savefig("EnergyGraph.pdf")
plt.show()
