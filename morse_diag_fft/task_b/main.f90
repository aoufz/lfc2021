PROGRAM main
USE fourier_mod
USE eigen_mod
IMPLICIT NONE

INTEGER, PARAMETER :: N = 2**9, M = 6
REAL(KIND=KIND(0.0d0)), PARAMETER:: L = 40, pi=4*ATAN(1.d0), a = 1.5, x0 = 5, e = exp(1.d0)
COMPLEX(KIND=KIND(0.0d0)), PARAMETER :: im = cmplx(0,1)
COMPLEX(KIND=KIND(0.0d0)), DIMENSION(0:N-1) :: tf , tftemp
REAL(KIND=KIND(0.0d0)), DIMENSION(0:N-1) :: x, k, f
INTEGER:: i, j, plan, o, p
COMPLEX(KIND=KIND(0.0d0)), DIMENSION(0:N-1,0:N-1) :: H
REAL(KIND=KIND(0.0d0)), DIMENSION(0:N-1) :: eval
COMPLEX(KIND=KIND(0.d0)),DIMENSION(0:N-1,0:N-1) :: evec
COMPLEX(KIND=KIND(0.0d0)), DIMENSION(0:N-1,0:N-1):: coeff, VT

!crea le griglie
CALL grids(L,N,x,k)

f = v(a,x,x0)

!calcola la trasformata
CALL fourier(f,N,tf)

!CALCOLO VALORE TEORICO
DO i = 0,N-1
    DO j = 0,N-1
        IF (i .ne. j) THEN
            VT(i,j) = (4*exp(-20*im*(k(i)-k(j)))*((200*(k(i)-k(j))**2-1)*sin(20*(k(i)-k(j)))+20*(k(i)-k(j))&
            &*cos(20*(k(i)-k(j)))))/((k(i)-k(j))**3*40)
        else
            VT(i,i) = L**2/12 + k(i)**2
        ENDIF
    END DO
END DO

!costruisce hamiltoniana
CALL ham(H,k,tf,N,L)

!risolve con zheevr
CALL resolve(N,H,coeff,eval,M)

!calcola autovettori
CALL eigenvec(N,L,k,x,coeff,evec)



!SCRITTURA SU FILE DEI DATI
OPEN(16, file="output.txt", action='write')
WRITE(16,*)"            x            k             f            RF               IF           "
DO i= 0,N-1
    WRITE(16,*) x(i), k(i), f(i), tf(i), VT(i,:)
END DO
CLOSE(16)

OPEN(17, file="hamiltonian.txt",action='write')
DO i = 0,N-1
    WRITE(17,*) H(i,:)
END DO
CLOSE(17)

OPEN(18, file="eval.txt", action="write")
WRITE(18,"(6f10.5)") eval(0:M-1)/2
WRITE(18,*) 0, N
CLOSE(18)

OPEN(19, file="evec.txt", action='write')
WRITE(19,"(8192f10.5)") x
DO i = 0,M-1
    WRITE(19,"(8192f10.5)") real(evec(i,:))
END DO
CLOSE(19)

END PROGRAM
