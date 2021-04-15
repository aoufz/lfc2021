MODULE mod_result
    USE mod_costr_morse
	IMPLICIT NONE	


	CONTAINS
	SUBROUTINE resolve(N,L,M,ueval,uevec,uerr,tol,a,x0,eval_last,isdone)

		!INPUT 
		INTEGER, INTENT(IN) ::  M, N, ueval, uevec, uerr
		REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L, tol, a, x0
		LOGICAL :: isdone
   		REAL(KIND = KIND(0.0d0)), DIMENSION(:) ::  eval_last
        !------

    	INTEGER :: F, ns, info, i, j, k, c, p, o
		CHARACTER(LEN=15) :: fm_eval, fm_evec
		REAL(KIND = KIND(0.0d0)) :: h, err
   		REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: va, x, ham_i, ham_i_1, work, work1 , eval
   		REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:,:) :: evec
   		INTEGER, ALLOCATABLE, DIMENSION(:) :: ib, is, iwork, ifail

   		!setup array
    	ALLOCATE(x(N))
		ALLOCATE(va(N))
		ALLOCATE(ham_i(N))
    	ALLOCATE(ham_i_1(N-1))
    	ALLOCATE(eval(N))
    	ALLOCATE(ib(N))
    	ALLOCATE(is(N))
    	ALLOCATE(work(4*N))
    	ALLOCATE(work1(5*N))
    	ALLOCATE(iwork(3*N))
    	ALLOCATE(ifail(M))
    	ALLOCATE(evec(N,M))
    		 
    	! costruzione griglia ed hamiltoniana
    	CALL grid(L,N,h,x)
		va = v(a,x,x0)
		CALL ham_cr(h,va,ham_i,ham_i_1)
		

		!autovalori----
    	!CALL DSTEBZ(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,NSPLIT,W,IBLOCK,ISPLIT,WORK,IWORK,INFO), vedi documentazione lapack per info su  	parametri in ingresso e uscita
   		CALL DSTEBZ('I', 'E', N, 1, M, 1, M, tol, ham_i, ham_i_1, F, ns, eval, ib,is,work, iwork, info)
   	
   	 	!teniamo conto del fatto che gli eval calcolati sono 2*E
        eval = eval/2
		err = maxval(abs(eval-eval_last))
		WRITE(uerr,*) err

		!Se abbiamo raggiunto la convergenza, stampiamo su file
		IF (err <= tol) THEN
			WRITE(uerr,*) "done for a=",a 
			isdone = .true.
        	!scriviamo su file i risultati ottenuti per gli autovalori
			WRITE(fm_eval,'(a,i2,a)') "(",M,"f10.5)"  !impostiamo il formato usando M come variabile di ripetizione
			WRITE(fm_evec,'(a,i4,a)') "(",N,"f10.5)"  !impostiamo il formato usando N come variabile di ripetizione
  	    	WRITE(ueval, fm_eval) (eval(c), c=1,M)  ! inserisce autovalori nel file
			 
		
        	!CALL DSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAIL,INFO), anche qui vedi documentazione lapack per info sui parametri
        	CALL DSTEIN(N,ham_i,ham_i_1,M,2*eval,ib,is,evec,N,work1,iwork,ifail,info)  
                 
        	!scriviamo su file i risultati ottenuti per gli autovettori
			WRITE(uevec, fm_evec) (x(o), o = 1,N)
        	DO p = 1,M  !indice di riga (scrive un autovettore per riga)
				WRITE(uevec, fm_evec) (evec(o,p), o = 1,N)
        	END DO

			RETURN
		END IF

		!setup per il prossimo ciclo di esecuzione
		eval_last = eval

        DEALLOCATE(x)
		DEALLOCATE(va)
		DEALLOCATE(ham_i)
    	DEALLOCATE(ham_i_1)
    	DEALLOCATE(eval)
    	DEALLOCATE(ib)
    	DEALLOCATE(is)
    	DEALLOCATE(work)
		DEALLOCATE(work1)
   		DEALLOCATE(iwork)
   		DEALLOCATE(ifail)
   		DEALLOCATE(evec)       
    
	END SUBROUTINE 

END MODULE
