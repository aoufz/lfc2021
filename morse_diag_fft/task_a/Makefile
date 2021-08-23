MF= Makefile

FC= gfortran

FFLAGS= -g 
#FFLAGS= -g -I/usr/include
LFLAGS= $(FFLAGS) -L/usr/local/lib -llapack
LFLAGS= $(FFLAGS) -L/usr/lib -llapack

EXE= taska.exe

SRC= \
	main.f90 

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ= $(SRC:.f90=.o)
MOD= $(SRC:.f90=.mod)

.f90.o:
	$(FC) $(FFLAGS) -c mod_costr.f90 mod_resolve.f90 $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) mod_costr.o mod_resolve.o $(LFLAGS)

$(OBJ): $(MF)

mainFFTW.o: mainFFTW.f90 

clean:
	rm -f $(OBJ) $(MOD) $(EXE) core
