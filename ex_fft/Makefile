MF= Makefile

FC= gfortran

FFLAGS= -g -I/usr/local/include
#FFLAGS= -g -I/usr/include
LFLAGS= $(FFLAGS) -L/usr/local/lib -lfftw3
LFLAGS= $(FFLAGS) -L/usr/lib -lfftw3

EXE= mainFFTW.exe

SRC= \
	mainFFTW.f90 

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ= $(SRC:.f90=.o)
MOD= $(SRC:.f90=.mod)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ): $(MF)

mainFFTW.o: mainFFTW.f90 

clean:
	rm -f $(OBJ) $(MOD) $(EXE) core
