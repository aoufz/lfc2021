MODULE eigen_mod
IMPLICIT NONE

CONTAINS

    !CALCOLO AUTOVALORI E AUTOVETTORI DI UNA MATRICE HERMITIANA DI DIMENSIONE N
    SUBROUTINE resolve(N,H,evec,eval,M,tol)
        INTEGER, INTENT(IN) :: N,M
        REAL(KIND=KIND(0.d0)), INTENT(IN) :: tol
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(:,:), INTENT(IN) :: H
        REAL(KIND=KIND(0.d0)), DIMENSION(:), INTENT(OUT) :: eval
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(:,:), INTENT(OUT) :: evec

        CHARACTER :: jobz = 'V', uplo = 'l', range = 'i'
        INTEGER :: lda, lwork, info, mout, lrwork, liwork, i
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

        !ZHEEVR(JOBZ,RANGE,UPLO,N,H,LDA,VL,VU,IL,IU,ABSTOL,
        !    M,W,Z,LDZ,ISUPPZ,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
        CALL zheevr(jobz,range,uplo,N,H,lda,0.d0,100.d0,1,M,tol,&
        & mout,eval,evec,N,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

        IF (info .ne. 0) THEN
            WRITE(*,*) "ERROR IN LAPACK CALL"
            STOP
        END IF
        
    END SUBROUTINE

    !COMBINAZIONE LINEARE CON ONDE PIANE
    SUBROUTINE eigenvec(N,L,k,x,coeff,evec,M)
        
        INTEGER :: i, j
        REAL(KIND=KIND(0.d0)) :: norm
        COMPLEX(KIND=KIND(0.d0)) :: img = dcmplx(0,1)

        INTEGER, INTENT(IN) :: N, M 
        REAL(KIND=KIND(0.0d0)), INTENT(IN) :: L
        REAL(KIND=KIND(0.0d0)), DIMENSION(0:N-1), INTENT(IN) :: x,k 
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(0:N-1,M), INTENT(IN) :: coeff
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(M,0:N-1), INTENT(OUT) :: evec
        
        
        DO i = 1,M 
            DO j = 0,N-1 !CL DI ONDE PIANE
                evec(i,j) = 1/sqrt(L)*sum(coeff(:,i)*exp(img*k(:)*x(j)))
            END DO

            norm = 0  !NORMALIZZAZIONE AUTOVETTORI
            DO j = 1, N-1 
                norm = norm + L/(2*dble(N)) * (conjg(evec(i,j-1))*evec(i,j-1) + conjg(evec(i,j))*evec(i,j))
            END DO
            evec(i,:) = evec(i,:) / sqrt(norm)
        END DO
    

    END SUBROUTINE




END MODULE