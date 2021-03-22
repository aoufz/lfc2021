import numpy as np
from VMCHelp import *
from VMC_HO import *
import stats as CalcStatistics
 
#set up wave function as always
 
HO=HOClass()
HO.SetParams([0.8])
R=np.zeros((1,1))
R[0,0]=np.random.uniform(-1.5,1.5,1) 
print (R[0,0])
HO.SetCoords(R)
(optimizeList,bestAlpha,bestEnergy)=HO.Optimize()
print ("The optimum alpha is ",bestAlpha)
fig=plt.figure(2, figsize=(9, 4), dpi=100)
s1 = fig.add_subplot(1, 2, 1)
s2 = fig.add_subplot(1, 2, 2)
s1.errorbar(optimizeList[:,0],optimizeList[:,1],optimizeList[:,2], marker='.', mec='red', mfc='red', elinewidth=2, ecolor='red',  mew=2)
s2.plot(optimizeList[:,0],optimizeList[:,2], marker='.', mec='red', mfc='red', mew=2)
#plt.axis([0., 1., bestEnergy-0.1, 1.3*bestEnergy])
plt.savefig("optimize.pdf")
plt.show()
