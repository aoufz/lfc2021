import matplotlib.pyplot as plt
import numpy as np
import math as m
import cmath as cm
import re

input_txt='fort_input.txt'
params_txt='genparams.txt'

#LETTURA INPUT
with open(input_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    if index == 0:
        N = int(line_e[0])
    elif index == 1:
        L = float(line_e[0])
    elif index == 2:
        filenamepacket = line_e[2][1:-1]

f.close()

with open(params_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    if index == 0:
        k0 = float(line_e[0])
        sigma = float(line_e[1])
    
f.close()


#ULTERIORI PARAMETRI FUNZIONE
x0 = 2*L/3


#COSTRUZIONE VETTORI
x = np.array([i*L/N for i in range(0,N)])
packet = [0 for i in range(0,N)]
impacket = [0 for i in range(0,N)]
###

#COSTRUZIONE FUNZIONE
for i in range(0,N):
    packet[i] = (1/m.pow(2*m.pi*sigma**2,1/4))*m.exp(-((x[i]-x0)**2)/(4*sigma**2))*m.cos(k0*x[i])
    impacket[i] =(1/m.pow(2*m.pi*sigma**2,1/4))*m.exp(-(x[i]-x0)**2/(4*sigma**2))*m.sin(k0*x[i])

f = open(filenamepacket,"w")
for i in range(0,N):
    f.write("(%.5f,%.5f) " % (packet[i], impacket[i]))
f.close()

#PLOT FUNZIONE PER DEBUG E/O CONFRONTO
fig,ax = plt.subplots()
ax.plot(x,packet)
fig.suptitle("Starting packet")
plt.savefig("debugpacket.jpg")