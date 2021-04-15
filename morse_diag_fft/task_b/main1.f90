PROGRAM main1
USE fourier_mod
USE eigen_mod
IMPLICIT NONE
    INTEGER :: N=2, M=2
    COMPLEX(KIND=KIND(0.0d0)), DIMENSION(2,2):: H, coeff
    REAL(KIND=KIND(0.0d0)), DIMENSION(2):: eval

    H(1,1) = 1
    H(1,2) = dcmplx(0,1)
    H(2,1) = -dcmplx(0,1)
    H(2,2) = 2

    CALL resolve(N,H,coeff,eval,M)

    PRINT *, eval


END PROGRAM