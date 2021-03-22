PROGRAM main
    USE mod_costr
	IMPLICIT NONE	

	!PARAMETRI PER LA SCELTA DELLA GRIGLIA E DEI FILE DI INPUT/OUTPUT
	REAL(KIND = KIND(0.0d0)), PARAMETER :: Lstart = 10, Lstop = 30, Lstep = 2
	INTEGER, PARAMETER ::  Nstart = 51, Nstop = 3001, Nstep = 150, M = 10 
	CHARACTER(LEN=*),PARAMETER :: filename1 = 'eval.txt', filename2 = 'evec.txt'
	!-------


    INTEGER :: N, F, ns, info, i, j, k, c, p, o, ueval=11, uevec=12
	CHARACTER(LEN=15) :: fm_eval, fm_evec
	REAL(KIND = KIND(0.0d0)) :: L, h, norm, tol = 1.0d-5
   	REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: x, ham_i, ham_i_1, work, work1, eval_true, eval, error
   	REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:,:) :: evec
   	INTEGER, ALLOCATABLE, DIMENSION(:) :: ib, is, iwork, ifail


	IF (MOD(Nstart,2) == 0 .or. MOD(Nstep,2) == 1) THEN
		WRITE (*,*) "Nstart is even or Nstep is odd: every iteration must have an odd N to exploit the hamiltonian's parity."
		STOP 
	END IF

    OPEN(UNIT=ueval, FILE = filename1)
	OPEN(UNIT=uevec, FILE = filename2)
	L = Lstart
   	DO WHILE (L >= Lstart .and. L<=Lstop)
		DO N = Nstart,Nstop,Nstep
    		ALLOCATE(x(N))
		    ALLOCATE(ham_i(N))
    		ALLOCATE(ham_i_1(N-1))
    		ALLOCATE(eval(N))
    		ALLOCATE(eval_true(M))
    	   	ALLOCATE(error(M))
    		ALLOCATE(ib(N))
    		ALLOCATE(is(N))
    		ALLOCATE(work(4*N))
    		ALLOCATE(work1(5*N))
    		ALLOCATE(iwork(3*N))
    		ALLOCATE(ifail(M))
    		ALLOCATE(evec(N,M))
    		 
    		! costruzione griglia ed hamiltoniana
    	 	CALL grid(L,N,h,x)
    		CALL ham_cr(h,x,ham_i,ham_i_1)

		    !trova autovalori
    	 	!CALL DSTEBZ(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,NSPLIT,W,IBLOCK,ISPLIT,WORK,IWORK,INFO), vedi documentazione lapack per info su  	parametri in ingresso e uscita
   		    CALL DSTEBZ('I', 'E', N, 1, M, 1, M, tol, ham_i, ham_i_1, F, ns, eval, ib,is,work, iwork, info)
   	
   	 	    !teniamo conto del fatto che gli eval calcolati sono 2*E
        	eval = eval/2
        
        	!calcoliamo il vettore che contiene i veri autovalori 
        	!calcoliamo il vettore che contiene gli errori sugli autovalori calcolati tramite diagonalizzazione di H
        	!calcoliamo la norma 1 del vettore errore (tenendo conto che le sue componenti sono già in v. ass.)
            DO k=1,M            	 
            	eval_true(k) = (k - 1.0) + 1/2.0
                error(k) = abs(eval_true(k) - eval(k))
			END DO 
			norm = SUM(error)
        	
        	!scriviamo su file i risultati ottenuti per gli autovalori
			WRITE(fm_eval,'(a,i2,a)') "(",M,"f10.5)"  !impostiamo il formato usando M come variabile di ripetizione
			WRITE(fm_evec,'(a,i4,a)') "(",N,"f10.5)"  !impostiamo il formato usando N come variabile di ripetizione
			WRITE(ueval,'(i5, 2f10.5)') N, L, norm
        	!WRITE(ueval, fm_eval) (eval(c), c=1,M) 
        	 
        	!trova autovettori
            !CALL DSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAIL,INFO), anche qui vedi documentazione lapack per info sui parametri
            CALL DSTEIN(N,ham_i,ham_i_1,M,eval,ib,is,evec,N,work1,iwork,ifail,info)  
                 
            !scriviamo su file i risultati ottenuti per gli autovettori
            WRITE(uevec, '(i5, f10.5)') N,L
			WRITE(uevec, fm_evec) (x(o), o = 1,N)
            DO p = 1,M  !indice di riga (scrive un autovettore per riga, )
				WRITE(uevec, fm_evec) (evec(o,p), o = 1,N)
            END DO
        	 
           	DEALLOCATE(x)
		    DEALLOCATE(ham_i)
    		DEALLOCATE(ham_i_1)
    		DEALLOCATE(eval)
    		DEALLOCATE(eval_true)
    	   	DEALLOCATE(error)
    		DEALLOCATE(ib)
    		DEALLOCATE(is)
    		DEALLOCATE(work)
    		DEALLOCATE(work1)
    		DEALLOCATE(iwork)
    		DEALLOCATE(ifail)
    		DEALLOCATE(evec)       
    		 
		END DO
		L = L + Lstep
    END DO
         
    CLOSE(ueval)
    CLOSE(uevec)
END PROGRAM
