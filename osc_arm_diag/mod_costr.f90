MODULE mod_costr
    IMPLICIT NONE   
    
    CONTAINS
    

    !COSTRUZIONE HAMILTONIANA
    SUBROUTINE ham_cr(h,x,ham_i,ham_i_1)
        IMPLICIT NONE 
        INTEGER :: i
        REAL(KIND = KIND(0.0d0)), INTENT(IN) :: h
        REAL(KIND = KIND(0.0d0)), INTENT(IN), DIMENSION(:) :: x 
        REAL(KIND = KIND(0.0d0)), INTENT(OUT), DIMENSION(SIZE(x)) :: ham_i
        REAL(KIND = KIND(0.0d0)), INTENT(OUT), DIMENSION(SIZE(x)-1) :: ham_i_1
 
        ham_i(:) = 2/(h**2) + x(:)**2
        ham_i_1(:) = -1/(h**2) 

    END SUBROUTINE 


    !COSTRUZIONE GRIGLIA
    SUBROUTINE grid(L,N,h,x)
        IMPLICIT NONE
        INTEGER :: i

        INTEGER, INTENT(IN) :: N
        REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L
        REAL(KIND = KIND(0.0d0)), INTENT(OUT) :: h
        REAL(KIND = KIND(0.0d0)), DIMENSION(:), INTENT(OUT) :: x

        
        h = L/REAL(N-1)
        x(1) = -L/2
        x(2:) = (/ (x(1) + real(i-1)*h , i = 2,N) /)
        
          
    END SUBROUTINE

END MODULE
