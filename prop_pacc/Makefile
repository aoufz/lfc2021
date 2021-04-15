MF= Makefile

FC= gfortran

FFLAGS= -g -I/usr/local/include 
#FFLAGS= -g -I/usr/include
LFLAGS= $(FFLAGS) -L/usr/local/lib -lfftw3 -llapack
LFLAGS= $(FFLAGS) -L/usr/lib -lfftw3 -llapack

EXE = main.exe
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
	$(FC) $(FFLAGS) -c mod1.f90 $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC)  $(FFLAGS) -o $(EXE) mod1.o $(OBJ) $(LFLAGS)

$(OBJ): $(MF)


clean:
	rm -f $(OBJ) $(MOD) $(EXE) core
