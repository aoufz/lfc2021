import subprocess
import os 
import re

input_txt = 'input.txt'
input_fort_txt = 'fort_input.txt'

N = []
M = []
free = []

#LETTURA INPUT
with open(input_txt) as f:
    con = f.readlines()
    
for index,line in enumerate(con):
    line_e = re.findall(r'\S+',line)

    if index == 0: 
        mode = int(line_e[0])
        output_txt = line_e[1][1:-1]
    elif index == 1: 
        L = float(line_e[0])
        LT = float(line_e[1])
        phi = int(line_e[2])
        E1 = float(line_e[3])
        E2 = float(line_e[4])
    elif index == 2:
        filenamemonitor = line_e[0][1:-1]
        filenameevo = line_e[1][1:-1]
        filenamepacket = line_e[2][1:-1]
    else:
        if mode == 0:
            N.append(line_e[0])
            M.append(line_e[1])
            free.append(line_e[2])
        elif mode == 1:
            Nstart = int(line_e[0])
            Nstop = int(line_e[1])
            Nstep = int(line_e[2])
            Mstart = int(line_e[3])
            Mstop = int(line_e[4])
            Mstep = int(line_e[5])
            freeauto = int(line_e[6])
            
f.close()

#ELIMINIAMO IL FILE DI OUTPUT SE ESISTE DA UNA PRECEDENTE RUN
if (os.path.exists(output_txt)):
    os.remove(output_txt)


if mode == 0: #MODALITA' MANUALE

    for i in range (0,len(N)):

        #SCRIVIAMO SU INPUT_FORT
        f = open(input_fort_txt,"w")
    
        f.write(str(N[i])+" "+str(M[i])+" "+str(free[i])+"\n")
        f.write(str(L)+" "+str(LT)+" "+str(phi)+" "+ str(E1)+" "+str(E2)+"\n")
        f.write("'"+filenamemonitor+"' '"+filenameevo+"' '"+filenamepacket+"' '"+output_txt+"' \n")

        f.close()

        #ESEGUIAMO IL FORTRAN PER OGNI ELEMENTO
        bashcomm1 = "python3 gen.py"
        bashcomm2 = "make"
        bashcomm3 = "./main.exe"
        process = subprocess.run(bashcomm1.split(), stdout=subprocess.PIPE)
        print(process)
        process = subprocess.run(bashcomm2.split(), stdout=subprocess.PIPE)
        print(process)
        process = subprocess.run(bashcomm3.split(), stdout=subprocess.PIPE)
        print(process)

elif mode == 1: #MODALITA' AUTOMATICA

    N = [i for i in range(Nstart,Nstop+1,Nstep)]
    M = [i for i in range(Mstart,Mstop+1,Mstep)]
    for i in range (0,len(N)):
        for j in range(0,len(M)):

            #SCRIVIAMO SU INPUT_FORT
            f = open(input_fort_txt,"w")

            f.write(str(N[i])+" "+str(M[j])+" "+str(freeauto)+"\n")
            f.write(str(L)+" "+str(LT)+" "+str(phi)+" "+ str(E1)+" "+str(E2)+"\n")
            f.write("'"+filenamemonitor+"' '"+filenameevo+"' '"+filenamepacket+"' '"+output_txt+"' \n")

            f.close()

            #ESEGUIAMO IL FORTRAN PER OGNI ELEMENTO
            bashcomm1 = "python3 gen.py"
            bashcomm2 = "make"
            bashcomm3 = "./main.exe"
            process = subprocess.run(bashcomm1.split(), stdout=subprocess.PIPE)
            print(process)
            process = subprocess.run(bashcomm2.split(), stdout=subprocess.PIPE)
            print(process)
            process = subprocess.run(bashcomm3.split(), stdout=subprocess.PIPE)
            print(process)

else: 

    print("Scegli una modalità di esecuzione compatibile tra 0 (manuale) e 1 (automatica)")
    exit()