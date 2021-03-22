PROGRAM main
	USE mod_result
	IMPLICIT NONE

	!PARAMETRI CHE VERRANNO LETTI DA INPUT
	REAL(KIND = KIND(0.0d0)) :: Lstart, Lstop, Lstep , tol, L
	INTEGER  ::  Nstart, Nstop, Nstep, M, ueval = 10, uevec = 11, unitinput = 12, N
	CHARACTER(LEN=10) :: filenameeval, filenameevec, filenameinput = 'input.txt'

	OPEN(UNIT = unitinput, FILE = filenameinput, ACTION = 'READ')

	READ(unitinput, *) Lstart, Lstop, Lstep, tol
	READ(unitinput, *) Nstart, Nstop, Nstep, M 
	READ(unitinput, *) filenameeval, filenameevec 

	CLOSE(UNIT = unitinput)

	!CHECK SUI VALORI SCELTI
	IF (MOD(Nstart,2) == 0 .or. MOD(Nstep,2) == 1) THEN
		WRITE (*,*) "Nstart is even or Nstep is odd: every iteration must have an odd N to exploit the hamiltonian's parity."
		STOP 
	END IF

	!ITERATION
	OPEN(UNIT=ueval, FILE = filenameeval)
	OPEN(UNIT=uevec, FILE = filenameevec)
	
	L = Lstart
   	DO WHILE (L >= Lstart .and. L<=Lstop)
		DO N = Nstart,Nstop,Nstep
    		CALL resolve(N,L,M,ueval,uevec,tol) 
		END DO
		L = L + Lstep
    END DO
         
    CLOSE(ueval)
    CLOSE(uevec)


END PROGRAM
