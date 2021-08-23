PROGRAM main
USE mod1
IMPLICIT NONE

!INIT
INTEGER :: unitinput = 13, uevo = 14, umonitor = 15, upacket = 16, uint = 17
INTEGER :: N, M, phi, free
REAL(KIND=KIND(0.d0)) :: L , LT, E1, E2
CHARACTER(LEN=15) :: filenameinput = 'fort_input.txt', filenamepacket, filenameevo, filenamemonitor, filenameint 
!---

INTEGER :: i, counter, counterlim
REAL(KIND=KIND(0.d0)) :: dt, dx, dk, stddevE, veldiff, sigmadiff, normdiff
REAL(KIND=KIND(0.d0)), ALLOCATABLE, DIMENSION(:) :: x,k,V,t
REAL(KIND=KIND(0.d0)), ALLOCATABLE, DIMENSION(:) :: E,sigma,xbar,tsample,norm,vel
COMPLEX(KIND=KIND(0.d0)), ALLOCATABLE, DIMENSION(:) :: packet, trpacket
COMPLEX(KIND=KIND(0.d0)), ALLOCATABLE, DIMENSION(:) :: propV, propT
COMPLEX(KIND=KIND(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: evolution

!LETTURA INPUT
OPEN (unit = unitinput, file = filenameinput, action = 'READ')

READ(unitinput,*) N, M, free
READ(unitinput,*) L, LT, phi, E1, E2
READ(unitinput,*) filenamemonitor, filenameevo, filenamepacket, filenameint

CLOSE(unitinput)

!CHECK SE FREE HA SENSO
IF ((free .ne. 0) .and. (free .ne. 1)) THEN
    WRITE(*,*) "Scegli se eseguire per pacchetto d'onda libero (1) oppure in presenza di un potenziale (0)"
    STOP 
ENDIF

counterlim = floor(LT/phi) 
ALLOCATE(x(0:N-1))
ALLOCATE(k(0:N-1))
ALLOCATE(V(0:N-1))
ALLOCATE(propT(0:N-1))
ALLOCATE(propV(0:N-1))
ALLOCATE(t(0:M))
ALLOCATE(tsample(0:counterlim))
ALLOCATE(norm(0:counterlim))
ALLOCATE(vel(0:counterlim))
ALLOCATE(E(0:counterlim))
ALLOCATE(sigma(0:counterlim))
ALLOCATE(xbar(0:counterlim))
ALLOCATE(evolution(0:counterlim,0:N-1))

!COSTRUZIONE GRIGLIA E PROPAGATORI PRECALCOLATI
CALL grids(x,t,k,L,LT,N,M,dt,dx,dk)
CALL makeV(x,dx,dt,L,N,E1,E2,free,propV,V)
propT = makeT(k,dt,N)

ALLOCATE(packet(0:N-1))
ALLOCATE(trpacket(0:N-1))

!LETTURA PACCHETTO
OPEN(upacket,file=filenamepacket,action = 'READ')
READ(upacket,*) packet
CLOSE(upacket)

!PASSO 0
evolution(0,:) = packet(:)
CALL monitor(N,L,k,x,packet,V,E(0),sigma(0),xbar(0),vel(0),norm(0))
tsample(0) = 0
counter = 1

!ITERAZIONE
DO i = 1,M
    !AZIONE 1 IN SPAZIO REALE PROPV
    packet(:) = propV(:)*packet(:)  
    !AZIONE 2 IN SPAZIO K PROPT
    CALL fourier(packet,N,trpacket,0)
    trpacket(:) = propT(:)*trpacket(:)
    CALL fourier(trpacket,N,packet,1)
    !AZIONE 3 IN SPAZIO REALE PROPV
    packet(:) = propV(:)*packet(:)
    !CONTROLLO SE E' TEMPO DI MONITORARE
    IF ( t(i) >= phi*counter ) THEN
        tsample(counter) = t(i)
        CALL monitor(N,L,k,x,packet,V,E(counter),sigma(counter),xbar(counter),vel(counter),norm(counter))
        evolution(counter,:) = packet(:)
        counter = counter + 1
    ENDIF

END DO

!CONTROLLO INTEGRALI DEL MOTO 
CALL mintegrals(counterlim,L,free,E,vel,sigma,norm,tsample,counter,phi,stddevE,veldiff,sigmadiff,normdiff)

!SCRITTURA EVOLUZIONE
OPEN(uevo,file=filenameevo)
WRITE(uevo,*) x
WRITE(uevo,*) V
DO i = 0,counter-1
    WRITE(uevo,*) real(evolution(i,:)*conjg(evolution(i,:)))
END DO 
CLOSE(uevo)

!SCRITTURA MONITOR
OPEN(umonitor,file=filenamemonitor)
WRITE(umonitor,*) "         Energy(t)              Sigma(t)                xbar(t)                  norm(t)            vel(t)"
DO i = 0,counter-1
    WRITE(umonitor,*) E(i), sigma(i), xbar(i), norm(i), vel(i)
END DO
CLOSE(umonitor)

!SCRITTURA OUTPUT
OPEN(uint,file=filenameint,position='append')
WRITE(uint,*) "         free        N            L             M               LT                         stddevE             &
& veldiff               sigmadiff                  normdiff"
WRITE(uint,*) free, N,L,M,LT,stddevE, veldiff, sigmadiff,normdiff
CLOSE(uint)

DEALLOCATE(x)
DEALLOCATE(k)
DEALLOCATE(V)
DEALLOCATE(propT)
DEALLOCATE(propV)
DEALLOCATE(t)
DEALLOCATE(tsample)
DEALLOCATE(E)
DEALLOCATE(sigma)
DEALLOCATE(xbar)
DEALLOCATE(vel)
DEALLOCATE(evolution)
DEALLOCATE(packet)
DEALLOCATE(trpacket)
DEALLOCATE(norm)

END PROGRAM