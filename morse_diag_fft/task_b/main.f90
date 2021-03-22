PROGRAM main
IMPLICIT NONE
!INCLUDE 'fftw3.f'

INTEGER, PARAMETER :: N = 2**10
REAL(KIND=KIND(0.0d0)), PARAMETER:: L = 20, pi=4*ATAN(1.d0)
REAL(KIND=KIND(0.0d0)) :: deltax, deltak
REAL(KIND=KIND(0.0d0)), DIMENSION(0:N-1,0:N-1) ::  H, kt
COMPLEX(KIND=KIND(0.0d0)), DIMENSION(0:N-1) :: f, T
REAL(KIND=KIND(0.0d0)), DIMENSION(0:N-1) :: x, k
INTEGER:: i, j, plan, o, p, q

deltax=L/dble(N)
deltak= pi/L

!costruzione funzione f(x) = x**2, griglia e vettori d'onda
DO i = 0, N-1
    x(i) = -L/2 + i*deltax
    k(i) = i*deltak
END DO

DO o = 0:N-1
    DO p = 0:N-1
        kt(o,p) = k(o) - k(p)
    END DO
END DO

f(:) = x(:)**2

DO q = 0:N-1
CALL dfftw_plan_dft_1d(plan,N,f,T,FFTW_forward,fftw_estimate)
CALL dfftw_execute(plan,f,T)
CALL dfftw_destroy_plan(plan)


END PROGRAM