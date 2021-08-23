MODULE mod_costr
    IMPLICIT NONE   
    
    CONTAINS

    !CREAZIONE HAMILTONIANA 
    SUBROUTINE ham_cr(h,v,ham_i,ham_i_1)
        IMPLICIT NONE 
        INTEGER :: i
        REAL(KIND = KIND(0.0d0)), INTENT(IN) :: h
        REAL(KIND = KIND(0.0d0)), INTENT(IN), DIMENSION(:) :: v
        REAL(KIND = KIND(0.0d0)), INTENT(OUT), DIMENSION(SIZE(v)) :: ham_i
        REAL(KIND = KIND(0.0d0)), INTENT(OUT), DIMENSION(SIZE(v)-1) :: ham_i_1


        ham_i(:) = 2/(h**2) + v(:)
        ham_i_1(:) = -1/(h**2) 


    END SUBROUTINE

    !CREAZIONE GRIGLIA
    SUBROUTINE grid(L,N,h,x)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: N
        REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L
        REAL(KIND = KIND(0.0d0)), INTENT(OUT) :: h
        REAL(KIND = KIND(0.0d0)), DIMENSION(1:N), INTENT(OUT) :: x

        h = L/REAL(N-1)
	    x =(/ ( real(i-1)*h, i = 1,N ) /)

               
   END SUBROUTINE

   !FUNZIONE POTENZIALE DI MORSE
    FUNCTION v(a,x,x0)
        IMPLICIT NONE
	    REAL(KIND = KIND(0.0d0)) :: a, x0
	    REAL(KIND = KIND(0.0d0)), DIMENSION(:):: x
	    REAL(KIND = KIND(0.0d0)), DIMENSION(SIZE(x)) :: w
	    REAL(KIND = KIND(0.0d0)), DIMENSION(SIZE(x)):: v

	    w(:) = -a*(x(:)-x0)
	    v(:) = exp(2*w(:)) - 2*exp(w(:))

    END FUNCTION
 

END MODULE
