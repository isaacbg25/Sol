PROGRAM solterra
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    real(kind=DP), parameter :: m_t=5.972168e24, m_s=1.9885e30, g=6.674302e-11, r_0=152098232e3, v_0=29.29e3, t_n = 365.25*3600*24
    integer, parameter :: n=300, d=n/10
    real(kind=DP), parameter :: L=r_0*v_0*m_t
    real(kind=DP) :: k_1, c_1
    real(kind=DP), parameter :: betha=(L/m_t)**2, k=g*m_s, alpha=betha/k, v_=k/(betha)**(1/2), dt = t_n/n!(t_n*v_)/(n*alpha)
    integer :: i,t 
    real(kind=DP) :: v(0:n),r(0:n),theta(0:n)
    v(0)=0
    r(0)=r_0
    do i = 0,n-1 !perquè ens dona en el temps i+1
        k_1 = -(1/r(i)**2) + ((1/r(i)**3))
        v(i+1)=v(i)+K_1*dt
        c_1 = v(i)
        r(i+1)= r(i)+dt*c_1 
    end do
    
    theta(0)=0
    do i = 0,n-1
        theta(i+1)=theta(i)+dt*L/(m_t*r(i)**2)
    end do
    
    open(unit=11, file="radi.txt", status="replace")
    do i= 0, d
        write(11, '(F30.15, 1X, F30.15)') dt*i*10, r(i*10)
    end do
    close(11)
    
    open(unit=12, file="orbita.txt", status="replace")
    do i= 0, d
        write(12, '(F30.15, 1X, F30.15)') r(i*10)*COS(theta(i*10)), r(i*10)*SIN(theta(i*10))
    end do
    close(12)
    
    
    
    END PROGRAM solterra