MODULE mod_resolve
    USE mod_costr
	IMPLICIT NONE	


	CONTAINS
	SUBROUTINE resolve(N,L,M,ueval,uevec,uerr,tol,a,x0,eval_last,isdone)

		
		INTEGER, INTENT(IN) ::  M, N, ueval, uevec, uerr
		REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L, tol, a, x0
		LOGICAL :: isdone
   		REAL(KIND = KIND(0.0d0)), DIMENSION(:) ::  eval_last
        

    	INTEGER :: F, ns, info, i, j, k, c, p, o
		CHARACTER(LEN=15) :: fm_eval, fm_evec
		REAL(KIND = KIND(0.0d0)) :: h, err, norm
   		REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: va, x, ham_i, ham_i_1, work, work1 , eval
   		REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:,:) :: evec
   		INTEGER, ALLOCATABLE, DIMENSION(:) :: ib, is, iwork, ifail

    	ALLOCATE(x(N))
		ALLOCATE(va(N))
		ALLOCATE(ham_i(N))
    	ALLOCATE(ham_i_1(N-1))
    	ALLOCATE(eval(M))
    	ALLOCATE(ib(N))
    	ALLOCATE(is(N))
    	ALLOCATE(work(4*N))
    	ALLOCATE(work1(5*N))
    	ALLOCATE(iwork(3*N))
    	ALLOCATE(ifail(M))
    	ALLOCATE(evec(N,M))
    		 
    	!COSTRUZIONE GRIGLIA ED HAMILTONIANA
    	CALL grid(L,N,h,x)
		va = v(a,x,x0)
		CALL ham_cr(h,va,ham_i,ham_i_1)
		

    	!DSTEBZ(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,NSPLIT,W,IBLOCK,ISPLIT,WORK,IWORK,INFO)
   		CALL dstebz('I', 'E', N, 1, M, 1, M, tol, ham_i, ham_i_1, F, ns, eval, ib,is,work, iwork, info)
		
		IF (info .ne. 0) THEN
			WRITE(*,*) "Error in LAPACK call"
			STOP
		ENDIF

   	 	!EVAL CALCOLATI = 2*EVAL
        eval = eval/2
		err = maxval(abs(eval-eval_last))
		WRITE(uerr,*) err

		!CONVERGENZA
		IF (err <= tol) THEN

			isdone = .true.
        	
			WRITE(fm_eval,'(a,i2,a)') "(",M,"f10.5)"  !formati
			WRITE(fm_evec,'(a,i4,a)') "(",N,"f10.5)" 
			
			
  	    	WRITE(ueval, fm_eval) eval
			 
        	!DSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAIL,INFO)
        	CALL dstein(N,ham_i,ham_i_1,M,2*eval,ib,is,evec,N,work1,iwork,ifail,info)  
                 
			WRITE(uevec, fm_evec) (x(o), o = 1,N)

        	DO p = 1,M !CALCOLO NORMA
				norm = 0  
				DO i = 2,N 
					norm =  norm + L/(2*dble(N))*(evec(i-1,p)**2+evec(i,p)**2)
				END DO
				evec(:,p) = evec(:,p) / sqrt(norm) !NORMALIZZAZIONE AUTOVETTORI
				
				WRITE(uevec, fm_evec) evec(:,p)
        	END DO

			RETURN

		END IF

		!SETUP PROSSIMO CICLO
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
