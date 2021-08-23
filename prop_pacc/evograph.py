import matplotlib.pyplot as plt
import numpy as np
import re
import os 
from os import listdir
import imageio 


input_txt = 'fort_input.txt'

#PULIZIA FOLDER "imgs"
path = "./imgs/"
for file_name in listdir(path):
    if file_name.endswith('.jpg'):
        os.remove(path+file_name)

#LETTURA INPUT
with open(input_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    if index == 0:
        N = int(line_e[0])
        free = int(line_e[2])
    elif index == 1:
        L = float(line_e[0])
    elif index == 2:
        filenamemonitor = line_e[0][1:-1]
        filenameevo = line_e[1][1:-1]

f.close()

E = []
sigma = []
xbar = []
vel = []

#LETTURA ENERGIE
with open(filenamemonitor) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    if index != 0:
        E.append(float(line_e[0]))
        if free == 1:
            sigma.append(float(line_e[1]))
            xbar.append(float(line_e[2]))
            vel.append(float(line_e[4]))

   
f.close()


x = []
packs = []
pot = []

#LETTURA PACCHETTI PER OGNI FRAME
with open(filenameevo) as f:
    con = f.readlines()

for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)
    for i,item in enumerate(line_e):
        line_e[i] = float(item)

    if(index == 0):
        x.append(line_e)
    elif(index==1):
        pot.append(line_e)
    else:
        packs.append(line_e)
    
f.close()
#------------------------


imgsfile = []
fig1, ax1 = plt.subplots(figsize = (8,6))
#PER OGNI FRAME GRAFICO L'AUTOFUNZIONE
for index in range(0,len(packs)):
    pack = packs[index]
   
    packplot, = ax1.plot(x[0],pack, color = 'red')
    packplot.set_label('Packet Plot')
    potplot, = ax1.plot(x[0],pot[0],color = 'blue')
    potplot.set_label('Potenziale')

    ax1.title.set_text(str(index))
    ax1.set_ylim([-0.05,max(packs[0])*3])
    ax1.legend()
    ax1.set_xlabel('Posizione')
    props = dict(boxstyle='round',facecolor='white',alpha=0.5)
    ax1.text(0.05,0.95,"Energy: %.4f" % E[index],transform=ax1.transAxes,fontsize = 11,verticalalignment='top',horizontalalignment='left',bbox = props)
    if free == 1:
        ax1.text(0.05,0.85,r"Sigma^2: %.2f, " % sigma[index] + r"xbar: %.2f" % xbar[index],transform=ax1.transAxes,fontsize = 11,verticalalignment='top',horizontalalignment='left', bbox=props)
        ax1.text(0.05,0.75,"vel: %.4f " % vel[index] ,transform=ax1.transAxes,fontsize = 11,verticalalignment='top',horizontalalignment='left', bbox=props)

    imgsfile.append('./imgs/evo_'+str(index)+'.jpg')
    plt.savefig('./imgs/evo_'+str(index)+'.jpg')
    plt.cla()

#CREO GIF
with imageio.get_writer('evolution_ani.gif', mode='I') as writer:
    for filename in imgsfile:
        image = imageio.imread(filename)
        writer.append_data(image)

