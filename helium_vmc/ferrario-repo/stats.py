import numpy as np

def UnweightedAvg(meanList,errorList):
    mean=np.sum(meanList)/len(meanList)
    error=0.0;
    for e in errorList:
        error=error+e*e
    error=sqrt(error)/len(errorList)
    return (mean,error)

def WeightedAvg (means, errors):
    zeroErrors = False
    for i in errors:
        if i == 0.0:
            zeroErrors = True
    if (not zeroErrors):
        weights = map (lambda x: 1.0/(x*x), errors)
        norm = 1.0/np.sum(weights)
        weights = map(lambda x: x*norm, weights)
        avg = 0.0
        error2 = 0.0
        for i in range (0,len(means)):
            avg = avg + means[i]*weights[i]
            error2 = error2 + weights[i]**2*errors[i]*errors[i]
        return (avg, np.sqrt(error2))
    else:
        return (np.sum(means)/len(means), 0.0)

def MeanErrorString (mean, error):
     if (mean!=0.0):
          meanDigits = math.floor(np.log(np.abs(mean))/np.log(10))
     else:
          meanDigits=2
     if (error!=0.0):
          rightDigits = -np.floor(np.log(error)/np.log(10))+1
     else:
          rightDigits=2
     if (rightDigits < 0):
          rightDigits = 0
     formatstr = '%1.' + '%d' % rightDigits + 'f'
     meanstr  = formatstr % mean
     errorstr = formatstr % error
     return (meanstr, errorstr)

def c(i,x,mean,var):
    N=len(x)
    if var==0:#if the variance is 0 return an effectively infinity corr
        return 1e100
#    print (len(x([0:N-1])),len(x([i:N])))
    corr=1.0/var*1.0/(N-i)*np.sum((x[0:N-i]-mean)*(x[i:N]-mean))
    return corr
                         
def Stats(x):
    N=len(x)
    mean=np.sum(x)/N
    xSquared=x*x
    var=np.sum(xSquared)/N-mean*mean
    i=0          
    tempC=0.5
    kappa=0.0
    while (tempC>0 and i<(N-1)):
        kappa=kappa+2.0*tempC
        i=i+1
        tempC=c(i,x,mean,var) 
    if kappa == 0.0:
        kappa = 1.0
    Neff=(N+0.0)/(kappa+0.0)
    error=np.sqrt(var/Neff)
    return (mean,var,error,kappa)
