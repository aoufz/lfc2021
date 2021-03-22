MODULE mod_result
    USE mod_costr
	IMPLICIT NONE	

	CONTAINS
	!Risolve il problema agli autovettori e autovalori dati i 3 parametri N,L,M, le unità di scrittura dei risultati e la tolleranza dell'algoritmo utilizzato
	SUBROUTINE resolve(N,L,M,ueval,uevec,tol)
		!INPUT 
		INTEGER, INTENT(IN) ::  M, N, ueval, uevec
		REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L, tol
        !------

    	INTEGER :: F, ns, info, i, j, k, c, p, o
		CHARACTER(LEN=15) :: fm_eval, fm_evec
		REAL(KIND = KIND(0.0d0)) :: h, norm
   		REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: x, ham_i, ham_i_1, work, work1, eval_true, eval, error
   		REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:,:) :: evec
   		INTEGER, ALLOCATABLE, DIMENSION(:) :: ib, is, iwork, ifail

   		
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
    		 
    	! costruzione griglia ed hamiltoniana----------------------
    	CALL grid(L,N,h,x)
		CALL ham_cr(h,x,ham_i,ham_i_1)
		


		!autovalori---------------------------
    	!CALL DSTEBZ(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,NSPLIT,W,IBLOCK,ISPLIT,WORK,IWORK,INFO), vedi documentazione lapack per info su  	parametri in ingresso e uscita
   		CALL DSTEBZ('I', 'E', N, 1, M, 1, M, tol, ham_i, ham_i_1, F, ns, eval, ib,is,work, iwork, info)
   	
   	 	!teniamo conto del fatto che gli eval calcolati sono 2*E
        eval = eval/2
        
        !calcoliamo il vettore che contiene gli autovalori analitici
        !calcoliamo il vettore che contiene gli errori sugli autovalori calcolati tramite diagonalizzazione di H
        !calcoliamo l'errore da associare alla computazione degli autovalori
                	 
        eval_true = (/ ( (k - 1.0) + 1/2.0, k = 1,M) /)
		error(:) = abs(eval_true(:) - eval(:))
		norm = max(error)
        	
        !scriviamo su file i risultati ottenuti per gli autovalori
		WRITE(fm_eval,'(a,i2,a)') "(",M,"f10.5)"  !impostiamo il formato usando M come variabile di ripetizione
		WRITE(fm_evec,'(a,i4,a)') "(",N,"f10.5)"  !impostiamo il formato usando N come variabile di ripetizione
		WRITE(ueval,'(i5, 2f10.5)') N, L, norm
        WRITE(ueval, fm_eval) (error(c), c=1,M)  !opzionale: stampa anche gli autovalori
			 
		

        !autovettori---------------------------------
        !CALL DSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAIL,INFO), anche qui vedi documentazione lapack per info sui parametri
        CALL DSTEIN(N,ham_i,ham_i_1,M,2*eval,ib,is,evec,N,work1,iwork,ifail,info)  
                 
        !scriviamo su file i risultati ottenuti per gli autovettori
        WRITE(uevec, '(i5, f10.5)') N,L
		WRITE(uevec, fm_evec) (x(o), o = 1,N)
        DO p = 1,M  !indice di riga (scrive un autovettore per riga)
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
    
	END SUBROUTINE 

END MODULE
