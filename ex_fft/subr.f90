module subr

contains

 subroutine cdft(N,switch, fd1, fd2)
 implicit none
 real(kind=8), parameter :: pi=3.141592653589793d0
 integer, intent (in)::N, switch
 integer :: j,k
 complex (kind=8), dimension (-N/2:N/2-1), intent(in) :: fd1
 complex (kind=8), dimension (-N/2:N/2), intent(out):: fd2
 complex (kind=8) :: w,wk
!
 if (switch==0) then
    w = dcmplx( cos(2.d0*pi/dble(N)) , sin(-2.d0*pi/dble(N)) )
 elseif (switch==1) then
    w = dcmplx( cos(2.d0*pi/dble(N)) , sin( 2.d0*pi/dble(N)) )
 endif
!
 do k=-N/2,N/2-1
    fd2(k)=(0.,0.)
    wk= w**k
    do j=-N/2,N/2
        fd2(k)=fd2(k)+fd1(j)*wk**j
    end do
    if (switch==0) then
        fd2(k)=fd2(k)/dble(N)
    endif
 end do
!
 end subroutine cdft

end module subr 
