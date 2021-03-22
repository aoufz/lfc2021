PROGRAM main
	USE mod_result
	IMPLICIT NONE

	!PARAMETRI CHE VERRANNO LETTI DA INPUT
	REAL(KIND = KIND(0.0d0)) :: Lstart, Lstop, Lstep , tol, L, astart, astop, astep, a, x0
	INTEGER  ::  Nstart, Nstop, Nstep, M, ueval = 10, uevec = 11, unitinput = 12, N
	CHARACTER(LEN=10) :: filenameeval, filenameevec, filenameinput = 'input.txt'

	OPEN(UNIT = unitinput, FILE = filenameinput, ACTION = 'READ')

	READ(unitinput, *) Lstart, Lstop, Lstep, tol
	READ(unitinput, *) Nstart, Nstop, Nstep, M 
	READ(unitinput, *) astart, astop, astep, x0
	READ(unitinput, *) filenameeval, filenameevec 

	CLOSE(UNIT = unitinput)

	!ITERATION
	OPEN(UNIT=ueval, FILE = filenameeval)
	OPEN(UNIT=uevec, FILE = filenameevec)
	
	a = astart
	DO WHILE (a >= astart .and. a<=astop)
		L = Lstart
   		DO WHILE (L >= Lstart .and. L<=Lstop)
			DO N = Nstart,Nstop,Nstep
    			CALL resolve(N,L,M,ueval,uevec,tol,a,x0) 
			END DO
			L = L + Lstep
    	END DO
		a = a + astep
	END DO
         
    CLOSE(ueval)
    CLOSE(uevec)


END PROGRAM
