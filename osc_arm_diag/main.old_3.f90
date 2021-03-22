PROGRAM main
    USE mod_result
	IMPLICIT NONE	

	!PARAMETRI PER LA SCELTA DELLA GRIGLIA E DEI FILE DI INPUT/OUTPUT
	REAL(KIND = KIND(0.0d0)), PARAMETER :: Lstart = 10, Lstop = 30, Lstep = 2, tol = 1.0d-5
	INTEGER, PARAMETER ::  Nstart = 51, Nstop = 3001, Nstep = 150, M = 10, ueval = 11, uevec = 12
	CHARACTER(LEN=*),PARAMETER :: filename1 = 'eval.txt', filename2 = 'evec.txt'
	!-------
	REAL(KIND = KIND(0.0d0)) :: L 
	INTEGER :: N

	IF (MOD(Nstart,2) == 0 .or. MOD(Nstep,2) == 1) THEN
		WRITE (*,*) "Nstart is even or Nstep is odd: every iteration must have an odd N to exploit the hamiltonian's parity."
		STOP 
	END IF

    OPEN(UNIT=ueval, FILE = filename1)
	OPEN(UNIT=uevec, FILE = filename2)
	!Impostiamo le iterazioni da Nstart a Nstop aumentando di Nstep, uguale per L
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
