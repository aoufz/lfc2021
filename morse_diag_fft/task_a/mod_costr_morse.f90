MODULE mod_costr_morse
    IMPLICIT NONE   
    
    CONTAINS
    
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

    SUBROUTINE grid(L,N,h,x)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: N
        REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L
        REAL(KIND = KIND(0.0d0)), INTENT(OUT) :: h
        REAL(KIND = KIND(0.0d0)), DIMENSION(:), INTENT(OUT) :: x

        h = L/REAL(N-1)
        x(1) = 0
	    x(2:) =(/ (x(1) + real(i-1)*h, i = 2,N) /)

               
   END SUBROUTINE

   FUNCTION v(a,x,x0)
	IMPLICIT NONE
	REAL(KIND = KIND(0.0d0)) :: a, x0
	REAL(KIND = KIND(0.0d0)), DIMENSION(:):: x
	REAL(KIND = KIND(0.0d0)), DIMENSION(SIZE(x)) :: w
	REAL(KIND = KIND(0.0d0)), DIMENSION(SIZE(x)):: v

	w(:) = -a*(x(:)-x0)
	v(:) = (1 - exp(w(:)))**2

   END FUNCTION
 

END MODULE
