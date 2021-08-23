import matplotlib.pyplot as plt
import numpy as np
import re
import os
from os import listdir

#pulizia folder "figures"
path = "./figures/"
for file_name in listdir(path):
    if file_name.endswith('.jpg'):
        os.remove(path+file_name)


#input file read
input_txt = "input.txt"
with open(input_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    if index == 0:
        Nstart = int(line_e[0])
        Nstop = int(line_e[1])
        Nstep = int(line_e[2])
        M = int(line_e[3])
    elif index == 1:
        L = float(line_e[0])
        astart = float(line_e[1])
        astop = float(line_e[2])
        astep = float(line_e[3])
        x0 = float(line_e[4])
    elif index == 2:
        tol = float(line_e[0])
        aval_txt = line_e[1][1:-1]
        avec_txt = line_e[2][1:-1]
        err_txt = line_e[3][1:-1]

        
f.close()

#creation of intervals
n = np.arange(Nstart,Nstop,Nstep, dtype=int)
a = np.arange(astart,astop,astep, dtype=float)

evalues = []
arr_info = []
x = []
evec = []
err = []

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

with open(err_txt) as f:
    con = f.readlines()

for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    for i,item in enumerate(line_e):
        line_e[i] = float(item)
    err.append(line_e[0])

f.close()
#------------------------



a = list(np.arange(astart,astop+astep,astep))
if (astart-astop) % astep == 0: #fix per aggiungere l'estremo
    a.append(astop+astep)

#per ogni valore di a faccio un grafico nella cartella figures
for index, aeval in enumerate(evalues):
    aevec = evec[index*M:(index+1)*M]
    fig1, ax1 = plt.subplots(2,int(M/2),sharex= True, sharey = True, figsize = (10,8))
    fig1.suptitle(t = "Prime M="+str(M)+" autofunzioni per N="+str(int(arr_info[index][1]))+", L="+str(L)+" e a="+str('%.2f' % a[index]), fontfamily='Serif',fontsize='x-large',fontweight='bold')
    
    for axi in ax1.flat:
        axi.label_outer()

    for i in range(0,M):
        if i < M/2:
            ax1[0,i].plot(x[index],[0 for i in range(int(arr_info[index][1]))],linestyle='--')
            ax1[0,i].plot(x[index],aevec[i], color = 'red')
            ax1[0,i].title.set_text("Autovalore: "+str(aeval[i]))
        else:
            ax1[1,i-int(M/2)].plot(x[index],[0 for i in range(int(arr_info[index][1]))],'--')
            ax1[1,i-int(M/2)].plot(x[index],aevec[i], color = 'red')
            ax1[1,i-int(M/2)].title.set_text("Autovalore: "+str(aeval[i]))


    plt.savefig('./figures/evec_graph_'+str('%.2f' % a[index])+'.jpg')

    cycles = int(arr_info[index][0])
   
    err_index = 0
    for i in range(index+1):
        err_index = err_index + int(arr_info[i][0]) + 1 
    

    fig2,ax2 = plt.subplots(figsize=(10,8))
    fig2.suptitle(t="Errore in funzione di N per a = "+str('%.2f' % a[index]),fontfamily='Serif',fontsize='x-large',fontweight='bold')
    ax2.set_xlabel("N")
    ax2.set_ylabel("Error")
    ax2.set_yscale("log")

    ax2.plot(n[1:cycles+1],err[err_index-cycles:err_index],color='red')

    plt.savefig('./figures/err_graph_'+str('%.2f' % a[index])+'.jpg')
    plt.close('all')


