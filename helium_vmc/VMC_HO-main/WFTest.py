from VMCHelp import *
import numpy
import random
from VMC_HO import *
 
HO=HOClass()
HO.SetParams([0.33])
coords=numpy.zeros((1,1),float)
coords[0,0]=0.8
energyVal=HO.WaveFunction(coords)
print ("Your wavefunction is working if this is ",numpy.exp(-0.33*0.8**2),energyVal)
