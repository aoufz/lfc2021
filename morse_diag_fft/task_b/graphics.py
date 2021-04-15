import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import re

input_txt = 'input.txt'
aval_txt = 'eval.txt'
avec_txt = 'evec.txt'
L = 40
N = 2**9
a=1.5
x0 = 5
M = 6
astart = 1.5
astop = 1.5
astep = 1

'''
#input file read
with open(input_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    if index == 0:
        L = float(line_e[0])
        tol = float(line_e[1])
    elif index == 1:
        Nstart = int(line_e[0])
        Nstop = int(line_e[1])
        Nstep = int(line_e[2])
        M = int(line_e[3])
    elif index == 2:
        astart = float(line_e[0])
        astop = float(line_e[1])
        astep = float(line_e[2])
        x0 = float(line_e[3])
    elif index == 3:
        aval_txt = line_e[0][1:-1]
        avec_txt = line_e[1][1:-1]
        err_txt = line_e[2][1:-1]
        
f.close()

#creation of intervals
n = np.arange(Nstart,Nstop,Nstep, dtype=int)
a = np.arange(astart,astop,astep, dtype=float)
'''
evalues = []
arr_info = []
x = []
evec = []

#data read
with open(aval_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    for i,item in enumerate(line_e):
        line_e[i] = float(item)

    if index % 2 == 0: 
        evalues.append(line_e)
    else:
        arr_info.append(line_e)
        


f.close()


with open(avec_txt) as f:
    con = f.readlines()

for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    for i,item in enumerate(line_e):
        line_e[i] = float(item)

    if(index % (M+1) == 0):
        x.append(line_e)
    else:
        evec.append(line_e)
    
f.close()
#------------------------



a = list(np.arange(astart,astop,astep))
if (astart-astop) % astep == 0: #fix per aggiungere l'estremo
    a.append(astop)

#per ogni valore di a faccio un grafico nella cartella figures
for index, aeval in enumerate(evalues):
    aevec = evec[index*M:(index+1)*M]
    fig1, ax1 = plt.subplots(2,int(M/2),sharex= True, sharey = True, figsize = (8,6))
    fig1.suptitle("Prime M="+str(M)+" autofunzioni per N="+str(arr_info[index][1])+", L="+str(L)+" e a="+str(a[index]))
    
    for axi in ax1.flat:
        axi.label_outer()

    for i in range(0,M):
        if i < M/2:
            ax1[0,i].plot(x[index],[0 for i in range(int(arr_info[index][1]))],linestyle='--')
            ax1[0,i].plot(x[index],aevec[i], color = 'red')
            ax1[0,i].title.set_text("aval: "+str(aeval[i]))
        else:
            ax1[1,i-int(M/2)].plot(x[index],[0 for i in range(int(arr_info[index][1]))],'--')
            ax1[1,i-int(M/2)].plot(x[index],aevec[i], color = 'red')
            ax1[1,i-int(M/2)].title.set_text("aval: "+str(aeval[i]))


    plt.savefig('./figures/evec_graph_'+str(a[index])+'.jpg')
    


"""

#costruiamo la matrice z, funzione di nu e lu


z = np.reshape(evalues, (len(l),len(n)))

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

"""