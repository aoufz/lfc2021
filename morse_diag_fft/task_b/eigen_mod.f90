MODULE eigen_mod
IMPLICIT NONE

CONTAINS

    !calcola autovalori e autovettori di una matrice H hermitiana di dimensione N
    SUBROUTINE resolve(N,H,evec,eval,M)
        INTEGER, INTENT(IN) :: N,M
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(:,:), INTENT(IN) :: H
        REAL(KIND=KIND(0.d0)), DIMENSION(:), INTENT(OUT) :: eval
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(:,:), INTENT(OUT) :: evec

        CHARACTER :: jobz = 'V', uplo = 'l', range = 'i'
        INTEGER :: lda, lwork, info, mout, lrwork, liwork
        INTEGER, DIMENSION(:),ALLOCATABLE :: isuppz,iwork
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(:), ALLOCATABLE :: work
        REAL(KIND=KIND(0.d0)), DIMENSION(:), ALLOCATABLE :: rwork
        
        lrwork=max(1,24*N)
        lda = max(1,N)
        lwork = max(1,2*N)
        liwork = max(1,10*N)
        
        ALLOCATE(iwork(max(1,liwork)))
        ALLOCATE(isuppz(2*max(1,M)))
        ALLOCATE(work(max(1,lwork)))
        ALLOCATE(rwork(max(1,lrwork)))
        !evec = H
        CALL zheevr(jobz,'i',uplo,N,H,lda,0.d0,100.d0,1,M,0.00000001,&
        & mout,eval,evec,N,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
        !CALL zheev(jobz,uplo,N,evec,lda,eval,work,lwork,rwork,info)
        print *, mout
        IF (info .ne. 0) THEN
            PRINT *, "ERROR IN LAPACK CALL"
            STOP
        END IF
        
    END SUBROUTINE

    SUBROUTINE eigenvec(N,L,k,x,coeff,evec)
        
        INTEGER :: i, j
        COMPLEX(KIND=KIND(0.d0)) :: img = dcmplx(0,1)

        INTEGER, INTENT(IN) :: N 
        REAL(KIND=KIND(0.0d0)), INTENT(IN) :: L
        REAL(KIND=KIND(0.0d0)), DIMENSION(0:N-1), INTENT(IN) :: x,k 
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(0:N-1,0:N-1), INTENT(IN) :: coeff
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(0:N-1,0:N-1), INTENT(OUT) :: evec
        
        DO i = 0,N-1
            DO j = 0,N-1
                evec(i,j) = 1/sqrt(L)*sum(coeff(:,i)*exp(img*k(:)*x(j)))
            END DO
        END DO

    

    END SUBROUTINE




END MODULE