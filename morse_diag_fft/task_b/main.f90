PROGRAM main
    USE fourier_mod
    USE eigen_mod
    IMPLICIT NONE

    !INIT
    INTEGER  :: Nstart, Nstop, Nstep, N, M, ueval=11, uevec = 12, unitinput = 13, uerr = 14
    REAL(KIND=KIND(0.0d0)) :: L, a, astart, astop, astep, x0, tol
    CHARACTER(LEN=10) :: filenameeval, filenameevec, filenameerr,filenameinput = 'input.txt'
    !---

    !CONTROLLO
    LOGICAL :: is_done
    REAL(KIND=KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: eval_last
    !---

    INTEGER :: i,step_counter
    REAL(KIND=KIND(0.d0)) :: err
    CHARACTER(LEN=15):: fm_eval, fm_evec
    REAL(KIND=KIND(0.0d0)), PARAMETER :: e = exp(1.d0), pi = 4*atan(1.d0)
    COMPLEX(KIND=KIND(0.0d0)), PARAMETER :: im = cmplx(0,1)
    COMPLEX(KIND=KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: tf 
    REAL(KIND=KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: x, k, f
    REAL(KIND=KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: eval
    COMPLEX(KIND=KIND(0.0d0)), ALLOCATABLE, DIMENSION(:,:) :: H, evec, coeff 


    !LETTURA DA FILE
    OPEN(UNIT = unitinput, FILE = filenameinput, ACTION = 'READ')

    READ(unitinput, *) Nstart, Nstop, Nstep, M
    READ(unitinput, *) L, astart, astop, astep, x0
    READ(unitinput, *) tol, filenameeval, filenameevec,filenameerr

    CLOSE(unitinput)

    ALLOCATE(eval_last(M))

    OPEN(ueval, file=filenameeval, action='WRITE')
    OPEN(uevec, file=filenameevec, action='WRITE')
    OPEN(uerr, file=filenameerr, action='WRITE')
 
    !ITERAZIONE
    a = astart
    DO WHILE (a>=astart .and. a<=astop)
        eval_last(:) = 100
        step_counter = 0
        is_done = .false.
    

        DO N=Nstart,Nstop,Nstep
            ALLOCATE(tf(0:N-1))
            ALLOCATE(x(0:N-1))
            ALLOCATE(k(0:N-1))
            ALLOCATE(f(0:N-1))
            ALLOCATE(H(0:N-1,0:N-1))
            ALLOCATE(eval(M))
            ALLOCATE(evec(M,0:N-1))
            ALLOCATE(coeff(0:N-1,M))
        
            
            !CREA GRIGLIE
            CALL grids(L,N,x,k)

            f = v(a,x,x0)

            !CALCOLO TRASFORMATA
            CALL fourier(f,N,tf)

            !COSTRUZIONE HAMILTONIANA
            CALL ham(H,k,tf,N,L)

            !RISOLVE CON ZHEEVR
            CALL resolve(N,H,coeff,eval,M,tol)
 
            !CALCOLO AUTOVETTORI
            CALL eigenvec(N,L,k,x,coeff,evec,M)

            !CALCOLO ERRORI
            err = maxval(abs(eval_last - eval))
            WRITE(uerr,*) err 
            
            !CONVERGENZA
            IF (err <= tol) THEN
                WRITE(fm_eval,'(a,i2,a)') "(",M,"f10.5)"
                WRITE(ueval,fm_eval) eval/2
                WRITE(ueval,*) step_counter,N,a

                WRITE(fm_evec,'(a,i4,a)') "(",N,"f10.5)"
                WRITE(uevec,fm_evec) x
                DO i = 1,M
                    WRITE(uevec,fm_evec) real(evec(i,:))
                END DO

                is_done = .true.

                DEALLOCATE(tf)
                DEALLOCATE(x)
                DEALLOCATE(k)
                DEALLOCATE(f)
                DEALLOCATE(H)
                DEALLOCATE(eval)
                DEALLOCATE(evec)
                DEALLOCATE(coeff)

                EXIT

            END IF

            !SETUP PROSSIMO CICLO
            step_counter = step_counter + 1
            eval_last = eval

            DEALLOCATE(tf)
            DEALLOCATE(x)
            DEALLOCATE(k)
            DEALLOCATE(f)
            DEALLOCATE(H)
            DEALLOCATE(eval)
            DEALLOCATE(evec)
            DEALLOCATE(coeff)
        END DO

        IF (.not. is_done) THEN
            WRITE(*,*) "Nessuna convergenza nel range di N fornito per a =",a 
        END IF
        a = a + astep
        
    END DO

    CLOSE(ueval)
    CLOSE(uevec)
    CLOSE(uerr)
    DEALLOCATE(eval_last)

END PROGRAM
