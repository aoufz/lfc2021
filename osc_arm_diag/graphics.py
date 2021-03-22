import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import re

#nomi file e parametri conosciuti
aval_txt = 'eval.txt'
avec_txt = 'evec.txt'
M = 10

#inizializzazione
n = np.array(dtype = float,object=[])
l = np.array(dtype = float,object=[])
x = np.array(dtype = float,object=[])
evalues_err = np.array(dtype = float,object=[])
evec_data = []
evec_nl_index = []
evec_target = []

#leggiamo i file riorganizzando i dati opportunamente-----------
with open(aval_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    if index % 2 == 0: 
        line_e = re.findall(r'\S+',line)
        n = np.append(n,float(line_e[0]))
        l = np.append(l,float(line_e[1]))
        evalues_err = np.append(evalues_err,float(line_e[2]))
       
    #scegliamo di ignorare la riga contentente gli errori sui singoli autovalori, per ora non utile
    
    
f.close()


with open(avec_txt) as f:
    con = f.readlines()

for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    for i,item in enumerate(line_e):
        line_e[i] = float(item)

    if(index % (M+2) == 0):
        evec_nl_index.append([line_e[0],line_e[1]])
    else:
        evec_data.append(line_e)
    
f.close()
#------------------------


#consideriamo i valori unici di n,l 
nu = np.unique(n)
lu = np.unique(l)

#costruiamo la matrice z, funzione di nu e lu
z = np.reshape(evalues_err, (len(lu),len(nu)))

#costruiamo le matrici x e y 
x = np.matrix(np.ones_like(lu)).T * np.matrix(nu)
y = np.matrix(np.ones_like(nu)).T * np.matrix(lu)

#grafichiamo i risultati in scala logaritmica in z
fig = plt.figure(figsize = (8,6))
ax = plt.axes(projection='3d',autoscale_on = True, xlabel = 'N', ylabel = 'L', zlabel = 'log_10(error)')
ax.set_title("Dipendenza dell'errore da N ed L")
ax.plot_surface(x,y.T,np.log10(z), cmap = plt.cm.coolwarm)

plt.savefig('error_graph.jpg')

#troviamo il minimo dell'errore: ricerca di Nstar e Lstar automatica -------
l_ind,n_ind = np.where(z == np.amin(z))

Nstar = nu[n_ind[0]]
Lstar = lu[l_ind[0]]
ch = False #variabile per capire se è stato settato un Nstar, Lstar manuale
#---------------------z

#Se si vuole settare un N,L manuale, lasciare la sezione non commentata ------
Nstar = 2901
Lstar = 10
n_a = np.logical_not(np.all(x-Nstar,axis=0))
l_a = np.logical_not(np.all(y-Lstar,axis=0))
n_ind = np.where(n_a[0,:])
l_ind = np.where(l_a[0,:])
ch = True
#---------------


#troviamo la posizione di Nstar e Lstar nell'array degli indici

for i in range(0,int(len(evec_nl_index))):
    if (evec_nl_index[i][0] == Nstar and evec_nl_index[i][1] == Lstar):
        position = i

#estraiamo dai dati raw la suddivisione dell'intervallo desiderata
x_evec_target = evec_data[(M+1)*position]

#estraiamo gli autovalori 
for i in range((M+1)*position+1, (M+1)*position + M+1):
    evec_target.append(evec_data[i])

if ch:
    inc = 1
#plottiamo le M autofunzioni in funzione di x
fig1, ax1 = plt.subplots(2,int(M/2),sharex= True, sharey = True, figsize = (8,6))
fig1.suptitle("Prime M="+str(M)+" autofunzioni per N="+str(Nstar)+" ed L="+str(Lstar)+" (Error = "+str(z[l_ind[inc],n_ind[inc]])+")")

for axi in ax1.flat:
    axi.label_outer()

for i in range(0,M):
    if i < M/2:
        ax1[0,i].plot(x_evec_target,evec_target[i], color = 'red')
    else:
        ax1[1,i-int(M/2)].plot(x_evec_target,evec_target[i], color = 'red')


plt.savefig('evec_graph.jpg')
