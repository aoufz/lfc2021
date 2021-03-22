 program main
 use subr 
 implicit none
!
 real (kind=8), parameter :: a=1.d0
 real (kind=8), parameter :: pi=3.141592653589793d0
 integer, parameter:: N=8192
 real(kind=8):: L, deltax,deltaf,fc,mu1
 integer::i,k,m
 complex(kind=8), dimension(0:N-1)::fd1, fd2, fd3
 real(kind=8), dimension(0:N-1):: x
 real(kind=8), dimension(0:N-1)::th, f
!
 L=20.d0/a
 deltax=L/dble(N)
 deltaf=2.d0*pi/L
 fc=N*deltaf
 mu1=0.d0
!
 do i= 0,N-1
     x(i)=-L/2+i*deltax
     if (x(i)<0.d0) then
         fd1(i)=exp(a*x(i))
     else
         fd1(i)=exp(-a*x(i))
     endif
 enddo
!
 do m=0,N-1
    f(m)=m*deltaf -fc/2.d0
    th(m)=a*2.d0/(a**2 + f(m)**2)
 enddo
!
 call cdft(N,0,fd1,fd2)
 call cdft(N,1,fd2,fd3)
!
 open(unit=1,file='ft.dat')
!
write(1,*) '#    i     x    ref(x)   imf(x)    k   th(k)    reF(k)    imF(k)    reF-1(x)   imF-1(x)' 
 do k=0,N-1
    write(1,fmt='(i7,g13.5,2g15.7,g13.4,5g15.7)') k,x(k),fd1(k),f(k),th(k),L*fd2(k),fd3(k)
 enddo
!
 end program
