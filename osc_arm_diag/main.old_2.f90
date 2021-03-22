program main
    USE mod_costr
    IMPLICIT NONE 

    ! i parametri del problema sono N,M e L
    INTEGER  :: N=1001,M=20,F,ns,info
    REAL(KIND = KIND(0.0d0)) :: L=30, tol = 0.00001
    ! fine dei parametri

    
    REAL(KIND = KIND(0.0d0)) :: h 
    REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:) :: x, ham_i, ham_i_1, eval, work
    REAL(KIND = KIND(0.0d0)), ALLOCATABLE, DIMENSION(:,:) :: evec
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ib, is, iwork, ifail

    
    ALLOCATE(x(N))
    ALLOCATE(ham_i(N))
    ALLOCATE(ham_i_1(N-1))
    ALLOCATE(eval(N))
    ALLOCATE(ib(N))
    ALLOCATE(is(N))
    ALLOCATE(work(4*N))
    ALLOCATE(iwork(3*N))
    ALLOCATE(evec(N,M))
    ALLOCATE(ifail(M))

    ! costruzione griglia ed hamiltoniana
    CALL grid(L,N,h,x)
    CALL ham_cr(h,x,ham_i,ham_i_1)

    !WRITE(*,*) ham_i
    !WRITE(*,*) ham_i_1

    ! ---------------------
    !trova autovalori
    !CALL DSTEBZ(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,NSPLIT,W,IBLOCK,ISPLIT,WORK,IWORK,INFO), vedi documentazione per info sui parametri in ingresso e uscita
    CALL DSTEBZ('I', 'E', N, 1, M, 1, M, tol, ham_i, ham_i_1, F, ns, eval, ib,is, work, iwork, info)

    !teniamo conto del fatto che eval_true = eval / 2
    eval = eval/2
    WRITE(*,*) eval(1:M), info 
    ! ---------------------



    ! --------------------
    !trova autovettori
    !CALL DSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAIL,INFO), anche qui vedi documentazione per info sui parametri
    CALL DSTEIN(N,ham_i,ham_i_1,M,eval,ib,is,evec,N,work,iwork,info)

    ! WRITE(*,*) evec, info
    ! -------------------

    

end program 