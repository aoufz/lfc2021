 program main
 !use subr
 implicit none
 include 'fftw3.f'
!
 real (kind=8), parameter :: a=1.d0
 real (kind=8), parameter :: pi=3.141592653589793d0
 integer, parameter:: N=2**10
 real(kind=8):: L, deltax,deltaf,fc,mu1
 integer::i,k,m
 complex(kind=8), dimension(0:N-1)::fd1, fd2, fd3
 real(kind=8), dimension(0:N-1):: x
 real(kind=8), dimension(0:N-1)::th, f
 integer(kind=8) :: plan1,plan2
!
 L=20.d0/a
 deltax=2*L/dble(N)
 deltaf=pi/L
 fc=N*deltaf
 mu1=0.d0
! costruzione della funzione f = e^(-|x|)
! do i= 0,N-1
!     x(i)=-L/2+i*deltax
!     if (x(i)<0.d0) then
!         fd1(i)=exp(a*x(i))
!     else
!         fd1(i)=exp(-a*x(i))
!    endif
! enddo
!
do i = 0,N-1
    x(i) = -L + i*deltax
end do

fd1 = x**2

 do m=0,N-1
     f(m)=m*deltaf 
     th(m)= (((f(m)**2)*(L**2)-2)*sin(f(m)*L)+2*f(m)*L*cos(f(m)*L))/(L*(f(m)**3))
 enddo
 th(0) = L**2/3
!
! call cdft(N,0,fd1,fd2)
! call cdft(N,1,fd2,fd3)
! 

! trasformata
call dfftw_plan_dft_1d(plan1,N,fd1,fd2,FFTW_forward,fftw_estimate)
call dfftw_execute(plan1,fd1,fd2)
call dfftw_destroy_plan(plan1)
fd2 = fd2/dble(N)
!antitrasfromata
call dfftw_plan_dft_1d(plan2,N,fd2,fd3,FFTW_backward,fftw_estimate)
call dfftw_execute(plan2,fd2,fd3)
call dfftw_destroy_plan(plan2)
!
open(unit=2,file='ft.dat')
!
write(2,*) '#    i     x              ref(x)   imf(x)    k   th(k)    reF(k)    imF(k)    reF-1(x)   imF-1(x)' 
 do k=0,N-1
    write(2,fmt='(i7,g13.5,2g15.7,g13.4,5g15.7)') k,x(k),fd1(k),f(k),th(k),fd2(k),fd3(k) !nota: L*fd2 è la trasformata
 enddo

open(unit=1,file='ft.txt')
do k=0,N-1
    write(1,*) fd2(k), fd3(k)
end do

end program
