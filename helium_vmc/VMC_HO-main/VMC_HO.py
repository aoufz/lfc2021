import numpy as np
from VMCHelp import *
import stats
 
 
class HOClass:
  def __init__(self):
    self.coord = np.zeros(1)
    self.SetParams([0.32])
    self.name = 'HO_WF'

  def SetParams(self,params):
     self.params=np.copy(params)
    
  def SetCoords(self,coords):
     self.coords=coords.copy()
    
  def WaveFunction(self,R):
     alpha=self.params[0]
     return   np.exp(-alpha*R**2)
    
  def LocalEnergy(self,R):
     KE=-0.5*LaplacianPsiOverPsi(R,self.WaveFunction)
     #change this to actually calculate V
     V=  0.5*R**2
     return V+KE
 
  def VMC(self,numSteps):
     EnergyList=np.zeros(numSteps)
     myHist=Histogram(-5,5,500)
     #pseudoCode:
     delta=1.618034*2.
     #Looping over steps:
     movesAttempted=0.0
     movesAccepted=0.0
     R=np.zeros((1,1) )
     R=self.coords.copy()
     nR=np.zeros((1,1))
     Psi=self.WaveFunction(R)
     energy=self.LocalEnergy(R)
     for step in range(numSteps):
       #choose new coordinates
       dR = np.random.uniform(-delta,delta,1)
       nR = R + dR
       #evaluate the wave function on these new cooridinates
       newPsi=self.WaveFunction(nR)
       if ( newPsi**2/Psi**2 > np.random.uniform() ):
         #accept the move making sure coords and oldWaveFunction has the up to date information
          R=nR.copy()
          Psi=newPsi
          movesAccepted+=1.
          energy=self.LocalEnergy(R)
          EnergyList[step]=energy
       else:
        #reject the move making sure coords and oldWaveFunction has the up to date information
          EnergyList[step]=energy
       movesAttempted+=1.
       myHist.add(R[0,0])
     print ("Acceptance ratio: ", movesAccepted/movesAttempted   )
     myHist.plotMe("density-"+str(self.params)+".pdf")
     return EnergyList
 
  def Optimize(self,npoints=15):
     import stats as CalcStatistics
     #add things here to loop over alpha and do optimization
     bestAlpha=0.
     bestEnergy=0.
     oldene=100.
     alpha = np.linspace(0.15,0.85, num=npoints, endpoint=True,)
     optimizeList=np.zeros( (npoints,3) )
     for i in range(npoints):
        self.SetParams([alpha[i]])
        energyList=self.VMC(100000)
        energy, variance, error, corrtime =CalcStatistics.Stats(np.array(energyList))
        optimizeList[i,0] = alpha[i]
        optimizeList[i,1] = energy
        optimizeList[i,2] = error
        print (" alpha= %g , energy= %g , var = %g , err = %g , ac-time= %g" % (alpha[i],energy,variance,error,corrtime))
        if energy<=oldene :
           bestAlpha=alpha[i]
           bestEnergy=energy
           oldene=energy
     return (optimizeList,bestAlpha,bestEnergy)
