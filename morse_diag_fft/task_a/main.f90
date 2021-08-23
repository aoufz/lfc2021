PROGRAM main
	USE mod_resolve
	IMPLICIT NONE

	!INIT
	REAL(KIND = KIND(0.0d0)) ::  tol, L, astart, astop, astep, a, x0
	INTEGER  ::  Nstart, Nstop, Nstep, M, ueval = 10, uevec = 11, unitinput = 12, uerr = 15, N
	CHARACTER(LEN=10) :: filenameeval, filenameevec, filenameerr, filenameinput = 'input.txt'

	!LAST EVALS, ISDONE E CONTEGGIO STEP
	REAL(KIND=KIND(0.d0)), ALLOCATABLE, DIMENSION(:) :: eval_last
	LOGICAL :: isdone = .false.
	INTEGER :: step_counter

	!LETTURA DA FILE
	OPEN(UNIT = unitinput, FILE = filenameinput, ACTION = 'READ')
	
	READ(unitinput, *) L, tol
	READ(unitinput, *) Nstart, Nstop, Nstep, M 
	READ(unitinput, *) astart, astop, astep, x0
	READ(unitinput, *) filenameeval, filenameevec, filenameerr

	CLOSE(UNIT = unitinput)

	ALLOCATE(eval_last(M))

	!ITERAZIONE
	OPEN(UNIT=ueval, FILE = filenameeval)
	OPEN(UNIT=uevec, FILE = filenameevec)
	OPEN(UNIT=uerr, FILE = filenameerr)

	a = astart
	DO WHILE (a >= astart .and. a<=astop)
		eval_last(:) = 0
		step_counter = 0
		DO N = Nstart,Nstop,Nstep
    		CALL resolve(N,L,M,ueval,uevec,uerr,tol,a,x0,eval_last,isdone) 
			IF (isdone) THEN
				WRITE(ueval,*)  step_counter, N, a
				EXIT
			END IF 
			step_counter = step_counter + 1				
		END DO

		!CONTROLLO CONVERGENZA 
		IF (isdone) THEN
			isdone = .false.
		ELSE			
			WRITE(*,*) "Nessuna convergenza nel range di N fornito per a = ",a
		ENDIF

		a = a + astep

	END DO
    
    CLOSE(ueval)
    CLOSE(uevec)
	CLOSE(uerr)
	DEALLOCATE(eval_last)
		
END PROGRAM
