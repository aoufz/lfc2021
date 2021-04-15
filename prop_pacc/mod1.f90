MODULE mod1
IMPLICIT NONE

CONTAINS 
    
    SUBROUTINE grids(x,t,L,T,N,M)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: N, M
        REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L, T
        REAL(KIND = KIND(0.0d0)), DIMENSION(0:N-1), INTENT(OUT) :: x,t
        REAL(KIND = KIND(0.0d0)) :: dx, dt


        dx = L/dble(N)
        dt = T/dble(M)
        DO i = 0,N-1
            x(i) = dble(i)*dx
            t(i) = dble(i)*dt
        END DO
    

    END SUBROUTINE





END MODULE
