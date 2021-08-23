import matplotlib.pyplot as plt
import numpy as np
import re
import os 
from os import listdir
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm



#LETTURA INPUT
input_txt = 'input.txt'
with open(input_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)

    if index == 0: 
        output_txt = line_e[1][1:-1]


#PULIZIA FOLDER "imgs"
path = "./imgs/"
for file_name in listdir(path):
    if file_name.endswith('err.JPG'):
        os.remove(path+file_name)

free = []
N = []
L = []
M = []
LT = []
devE = []
veldiff = []
sigmadiff = []
normdiff = []

#LETTURA OUTPUT
with open(output_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    if index % 2 != 0:
        free.append(int(line_e[0]))
        N.append(int(line_e[1]))
        L.append(float(line_e[2]))
        M.append(int(line_e[3]))
        LT.append(float(line_e[4]))
        devE.append(float(line_e[5]))
        veldiff.append(float(line_e[6]))
        sigmadiff.append(float(line_e[7]))
        normdiff.append(float(line_e[8]))

f.close()

#------------------------

dxdt = [L[i]/N[i] * M[i]/LT[i] for i in range(0,len(N))] 

#GRAFICO L'ERRORE IN FUNZIONE DI DELTAX/DELTAT  
fig1, ax1 = plt.subplots(figsize = (8,6))
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('deltax/deltat')
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))
eplot, = ax1.plot(dxdt,devE, color = 'red')
eplot.set_label('Errore su energia')
ax1.title.set_text("Errore su integrali del moto in funzione di dx/dt")
ax1.legend()
plt.savefig('./imgs/eerr.JPG')




#SE LIBERA, GRAFICO ANCHE XDIFF E SIGMADIFF
if free[0] == 1:
    fig2, ax2 = plt.subplots(figsize = (8,6))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('deltax/deltat')
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    xplot, = ax2.plot(dxdt,veldiff,color = 'blue')
    xplot.set_label('Errore sulla velocità')
    ax2.title.set_text("Errore su integrali del moto in funzione di dx/dt")
    ax2.legend()
    plt.savefig('./imgs/velerr.JPG')

    fig3, ax3 = plt.subplots(figsize = (8,6))
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('deltax/deltat')
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    sigmaplot, = ax3.plot(dxdt,sigmadiff,color = 'green')
    sigmaplot.set_label('Errore su sigma^2')
    ax3.title.set_text("Errore su integrali del moto in funzione di dx/dt")
    ax3.legend()
    plt.savefig('./imgs/sigmaerr.JPG')

fig4, ax4 = plt.subplots(figsize = (8,6))
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlabel('deltax/deltat')
ax4.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))
normplot, = ax4.plot(dxdt,normdiff, color = 'purple')
normplot.set_label('Errore sulla norma della funzione')
ax4.title.set_text("Errore su integrali del moto in funzione di dx/dt")
ax4.legend()
plt.savefig('./imgs/normerr.JPG')

