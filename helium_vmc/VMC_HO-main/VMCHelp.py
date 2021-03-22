import numpy as np
import matplotlib.pyplot as plt

def dist(x,y):
    return np.sqrt(np.dot(x-y,x-y))

#calcola numericamente il contributo cinetico all'energia locale, approssimando la derivata seconda con approssimazione centrale
def LaplacianPsiOverPsi(coords,WaveFunction):
    total=0.0
    delta=0.0001
    tempVal3=WaveFunction(coords)
    for i in range(0,len(coords)):
        for j in range(0,len(coords[0])):
            coords[i,j]=coords[i,j]+delta
            tempVal=WaveFunction(coords)
            coords[i,j]=coords[i,j]-2*delta
            tempVal2=WaveFunction(coords)
            coords[i,j]=coords[i,j]+delta
            total +=(tempVal+tempVal2)-2.0*tempVal3
    return total/(delta*delta*tempVal3)

class Histogram:
        min=0.0
        max=5.0
        delta=0.05
        numSamples=0
        bins=np.array([0])
        def add(self,val):
                if self.min<val and val<self.max:
                        index=int((val-self.min)/self.delta)
                        self.bins[index]=self.bins[index]+1
                self.numSamples=self.numSamples+1
        def printMe(self):
                for i in range(0,len(self.bins)):
                        print (self.min+self.delta*i,self.bins[i]/self.numSamples)

        def plotMe(self,fileName=""):
                print ("plotting")
                plt.clf()
                self.bins=self.bins/self.numSamples
                xCoord=[self.min+self.delta*i for i in range(0,len(self.bins))]
                plt.plot(xCoord,self.bins)
                plt.gca().xaxis.major.formatter.set_scientific(False)
                if not(fileName==""):
                   plt.savefig(fileName)
                else:
                   plt.show()
                
        def plotMeNorm(self,fileName):
                print ("plotting")
                plt.clf()
                self.bins=self.bins/self.numSamples
                xCoord=np.array([self.min+self.delta*i for i in range(0,len(self.bins))])
                plt.plot(xCoord,self.bins/xCoord**2)
                plt.gca().xaxis.major.formatter.set_scientific(False)
                plt.savefig(fileName)
                plt.show()

        def __init__(self,min,max,numBins):
                self.min=min
                self.max=max
                self.delta=(max-min)/numBins
                self.bins=np.zeros(numBins)
                numSamples=0
