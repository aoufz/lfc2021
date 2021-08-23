MODULE fourier_mod
IMPLICIT NONE 
INCLUDE 'fftw3.f'


CONTAINS 

!CREAZIONE GRIGLIE
SUBROUTINE grids(L,N,x,k)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER, INTENT(IN) :: N
    REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L
    REAL(KIND = KIND(0.0d0)) :: dx, dk
    REAL(KIND = KIND(0.0d0)), DIMENSION(0:N-1), INTENT(OUT) :: x,k 
    
    dx = L/(dble(N))
    dk = 2*4*atan(1.0d0)/L

    x(:) = (/ (dble(i)*dx, i = 0,N-1) /)
    
    DO i = 0, N/2-1
        k(i) = dble(i)*dk
    END DO
    DO i = -N/2, -1
        k(i+N) = dble(i)*dk
    END DO
    

END SUBROUTINE

!FUNZIONE POTENZIALE MORSE
FUNCTION v(a,x,x0)
	IMPLICIT NONE
	REAL(KIND = KIND(0.0d0)) :: a, x0
	REAL(KIND = KIND(0.0d0)), DIMENSION(:):: x
	REAL(KIND = KIND(0.0d0)), DIMENSION(SIZE(x)) :: w
	REAL(KIND = KIND(0.0d0)), DIMENSION(SIZE(x)):: v

    
	w(:) = -a*(x(:)-x0)
	v(:) = exp(2*w(:)) - 2*exp(w(:))

    !Potenziale armonico di prova
    !v(:) = (x-20)**2

END FUNCTION
 
!CALCOLO TRASFORMATA DI FOURIER
SUBROUTINE fourier(f,N,trf)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(KIND=KIND(0.0d0)), DIMENSION(0:N-1), INTENT(IN) :: f
    COMPLEX(KIND=KIND(0.0d0)), DIMENSION(0:N-1), INTENT(OUT) :: trf
    INTEGER :: plan, i

    CALL dfftw_plan_dft_1d(plan,N,dcmplx(f,0),trf,FFTW_FORWARD,FFTW_ESTIMATE)
    CALL dfftw_execute_dft(plan,dcmplx(f,0),trf)
    CALL dfftw_destroy_plan(plan)

    trf(0:N-1) = (1/dble(N))*trf(0:N-1) !NORMALIZZAZIONE TRASFORMATA
    
END SUBROUTINE

!CREAZIONE HAMILTONIANA
SUBROUTINE ham(H,k,trf,N,L)

    REAL(KIND=KIND(0.d0)), DIMENSION(0:N-1), INTENT(IN) :: k
    INTEGER :: i, j
    REAL(KIND=KIND(0.d0)), INTENT(IN) :: L
    INTEGER, INTENT(IN) :: N
    COMPLEX(KIND=KIND(0.d0)), DIMENSION(0:N-1), INTENT(IN) :: trf
    COMPLEX(KIND=KIND(0.d0)), DIMENSION(0:N-1,0:N-1), INTENT(OUT):: H

    
    DO i = 0,N-1
        DO j = 0,N-1
            IF (i .ne. j) THEN
                IF (i-j > 0) THEN
                    H(i,j) = trf(i-j)
                ELSE
                    H(i,j) = trf(N-(j-i)) 
                ENDIF
            ELSE   
                H(i,j) = k(i)**2 + trf(0) 
            END IF
        END DO
    END DO
    
END SUBROUTINE

END MODULE