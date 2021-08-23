MODULE mod1
IMPLICIT NONE
INCLUDE 'fftw3.f'

CONTAINS 
    
    !COSTRUZIONE GRIGLE 
    SUBROUTINE grids(x,t,k,L,LT,N,M,dt,dx,dk)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: N, M
        REAL(KIND = KIND(0.0d0)), INTENT(IN) :: L, LT
        REAL(KIND = KIND(0.0d0)), DIMENSION(0:N-1), INTENT(OUT) :: x,k
        REAL(KIND = KIND(0.0d0)), DIMENSION(0:M), INTENT(OUT) :: t
        REAL(KIND = KIND(0.0d0)), INTENT(OUT) :: dt, dx, dk

        dx = L/dble(N)
        dt = LT/dble(M)
        dk = 2*4*atan(1.0d0)/L

        x(:) = (/ (dble(i)*dx , i = 0,N-1) /)
        t(:) = (/ (dble(i)*dt , i = 0,M) /)
    
        DO i = 0, N/2-1
            k(i) = dble(i)*dk
        END DO

        DO i = -N/2, -1
            k(i+N) = dble(i)*dk
        END DO
    
    END SUBROUTINE

    !COSTRUZIONE PROPAGATORE V
    SUBROUTINE makeV(x,dx,dt,L,N,E1,E2,free,propV,V)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: N, free
        COMPLEX, PARAMETER :: im = dcmplx(0,1)
        REAL(KIND=KIND(0.0d0)),DIMENSION(0:N-1), INTENT(IN):: x
        REAL(KIND=KIND(0.0d0)),DIMENSION(0:N-1), INTENT(OUT):: V
        COMPLEX(KIND=KIND(0.0d0)),DIMENSION(0:N-1),INTENT(OUT):: propV
        REAL(KIND=KIND(0.0d0)), INTENT(IN) :: dx, dt, L, E1, E2
        REAL(KIND=KIND(0.0d0)) :: l1, l2

        !SE LIBERO, CALCOLA E RITORNA
        IF (free == 1) THEN
            V(:) = 0
            propV(:) = 1
            RETURN
        ENDIF  

        !SPESSORE BARRIERA
        l1 = L/137
        l2 = L/137

        !CREAZIONE POTENZIALE
        DO i = 0,N-1
            IF (L/3 - l1/2 < x(i) .and. x(i) < L/3 + l1/2) THEN
                V(i) = E1
            ELSE IF (L/2 - l2/2 < x(i) .and. x(i) < L/2 + l2/2) THEN 
                V(i) = E2
            ELSE 
                V(i) = 0
            END IF
        END DO
        
        propV(:) = exp(-im*(dt/2)*V(:))
        
    END SUBROUTINE

    !CREAZIONE PROPAGATORE T
    FUNCTION makeT(k,dt,N)
        IMPLICIT NONE 
        INTEGER :: N
        COMPLEX(KIND=KIND(0.d0)) :: im = dcmplx(0,1)
        REAL(KIND=KIND(0.d0)), DIMENSION(0:N-1) :: k
        REAL(KIND=KIND(0.d0)) :: dt 
        COMPLEX(KIND=KIND(0.d0)), DIMENSION(0:N-1) :: makeT

        makeT(:) = exp(-im * (dt/2) * k(:)**2)

    END FUNCTION

    !TRASFORMATA FOURIER, INVERSA OPPURE NO 
    SUBROUTINE fourier(f,N,trf,info)
        IMPLICIT NONE
        INTEGER, INTENT(IN) ::  N, info
        COMPLEX(KIND=KIND(0.0d0)), DIMENSION(0:N-1), INTENT(IN) :: f
        COMPLEX(KIND=KIND(0.0d0)), DIMENSION(0:N-1), INTENT(OUT) :: trf
        INTEGER :: plan, i

        
        IF (info == 0) THEN
            CALL dfftw_plan_dft_1d(plan,N,f,trf,FFTW_FORWARD,FFTW_ESTIMATE)
            CALL dfftw_execute_dft(plan,f,trf)
            CALL dfftw_destroy_plan(plan)
        ELSE IF (info == 1) THEN
            CALL dfftw_plan_dft_1d(plan,N,f,trf,FFTW_BACKWARD,FFTW_ESTIMATE)
            CALL dfftw_execute_dft(plan,f,trf)
            CALL dfftw_destroy_plan(plan)
        ENDIF
        
        trf(0:N-1) = trf(0:N-1)/sqrt((dble(N))) !NORMALIZZAZIONE
        
    END SUBROUTINE

    !MONITORA QUANTITA' RILEVANTI
    SUBROUTINE monitor(N,L,k,x,packet,V,E,sigma,xbar,vel,norm)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N 
        REAL(KIND=KIND(0.d0)), INTENT(IN) :: L
        REAL(KIND = KIND(0.0d0)), DIMENSION(0:N-1), INTENT(IN) :: k, x, V
        COMPLEX(KIND = KIND(0.0d0)), DIMENSION(0:N-1), INTENT(IN) :: packet
        COMPLEX(KIND = KIND(0.0d0)), DIMENSION(0:N-1) :: trpacket
        REAL(KIND=KIND(0.d0)) :: T,U, normtr
        REAL(KIND=KIND(0.d0)), INTENT(OUT) :: E, sigma, xbar, norm, vel
        INTEGER :: j

        
        norm = 0  !NORMALIZZAZIONE AUTOVETTORI
        DO j = 1, N-1 
            norm = norm + L/(2*dble(N)) * (conjg(packet(j-1))*packet(j-1) + conjg(packet(j))*packet(j))
        END DO
        

        U = 0 !VALORE ASPETTAZIONE ENERGIA POTENZIALE
        DO j = 1, N-1
            U = U + L/(2*dble(N)) * (V(j-1)*conjg(packet(j-1))*packet(j-1) + V(j)*conjg(packet(j))*packet(j))
        END DO
        U=U/norm
        
        CALL fourier(packet,N,trpacket,0)

        normtr = 0  !NORMALIZZAZIONE TRASFORMATA
        DO j = 1, N-1 
            normtr = normtr + (4*atan(1.d0)/L) * (conjg(trpacket(j-1))*trpacket(j-1) + conjg(trpacket(j))*trpacket(j))
        END DO
        
        T = 0 !VALORE ASPETTAZIONE ENERGIA CINETICA
        DO j = 1, N-1
            T = T +  (4*atan(1.d0)/L) * (k(j-1)**2/2*conjg(trpacket(j-1))*trpacket(j-1) + k(j)**2/2*conjg(trpacket(j))*trpacket(j))
        END DO
        T = T/normtr
        
        vel = 0 !CALCOLO VELOCITA' 
        DO j = 1,N-1
            vel = vel + (4*atan(1.d0)/L) * (k(j-1)*conjg(trpacket(j-1))*trpacket(j-1) + k(j)*conjg(trpacket(j))*trpacket(j))
        END DO
        vel = vel/normtr 


        xbar = 0 !VALORE MEDIO X
        DO j = 1, N-1
            xbar = xbar + L/(2*dble(N)) * (x(j-1)*conjg(packet(j-1))*packet(j-1) + x(j)*conjg(packet(j))*packet(j))
        END DO
        xbar = xbar/norm

        sigma = 0 !SIGMA^2
        DO j = 1, N-1
            sigma = sigma + L/(2*dble(N)) * ((x(j-1)-xbar)**2*conjg(packet(j-1))*packet(j-1) + (x(j)-xbar)**2*&
            & conjg(packet(j))*packet(j))
        END DO
        sigma = sigma/norm

        E = T + U

    END SUBROUTINE

    !CONTROLLO INTEGRALI DEL MOTO
    SUBROUTINE mintegrals(counterlim,L,free,E,vel,sigma,norm,tsample,counter,phi,stddevE,veldiff,sigmadiff,normdiff)
        INTEGER, INTENT(IN) :: free, counter, phi, counterlim
        REAL(KIND=KIND(0.d0)), DIMENSION(0:counterlim), INTENT(IN) :: E,vel,sigma,tsample,norm
        REAL(KIND=KIND(0.d0)), INTENT(IN) :: L
        REAL(KIND=KIND(0.d0)), INTENT(OUT) :: stddevE, veldiff, sigmadiff,normdiff

        INTEGER :: i, genunit = 20
        CHARACTER (LEN=15) :: params_txt = 'genparams.txt'
        REAL(KIND=KIND(0.d0)) :: Ei,vi,sigmai

        !LETTURA PARAMETRI INIZIALI PACCHETTO
        OPEN (unit = genunit, file = params_txt, action = 'READ')
        READ(genunit,*) vi, sigmai
        CLOSE(genunit)

        Ei = vi**2 / 2.d0

        !CALCOLO ERRORE SU ENERGIA
        stddevE = maxval(abs(Ei-E(1:counter-1)))

        !CALCOLO ERRORE SU NORMA
        normdiff = maxval(abs(1-norm(1:counter-1)))
        
        IF (free == 1) THEN
            
            !CALCOLO ERRORE SU VELOCITA'
            veldiff = maxval( abs(vi-vel(1:counter-1) ))

            !CALCOLO ERRORE SU SIGMA^2
            sigmadiff = maxval( abs( sigma(1:counter-1) - sigmai**2 - tsample(1:counter-1)**2/ (4*sigmai**2 ) ))

        ENDIF 

    END SUBROUTINE 

END MODULE
