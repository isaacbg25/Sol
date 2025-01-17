PROGRAM solterra
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    real(kind=DP), parameter :: m_t=5.972168e24, m_s=1.9885e30, g=6.674302e-11, r_0=147098070e3, v_0=30.29e3, t_n = 365*3600*24
    integer, parameter :: n=365*24*60 !discretització temporal
    real(kind=DP), parameter :: L=r_0*v_0*m_t
    real(kind=DP) :: k_1, c_1 !variables que utilitzarem a Euler
    real(kind=DP), parameter :: betha=(L/m_t)**2, k=g*m_s, alpha=betha/k, v_=k/(betha)**(0.5) !constants de normalització
    integer :: i
    real(kind=DP) :: v(0:n),r(0:n),theta(0:n) !velocitat radial, posició i angle
    real(kind=dp), parameter :: t_n_norm = t_n*v_/alpha, dt=t_n_norm/n
    
    ! ci de les edos en r
    v(0)=0
    r(0)=r_0/alpha
    
    !solució de les equacions per Euler
    do i = 0,n-1 !perquè ens dona en el temps i+1
        k_1 = -(1/r(i)**2) + ((1/r(i)**3))
        v(i+1)=v(i)+K_1*dt
        c_1 = v(i)
        r(i+1)= r(i)+dt*c_1 
    end do
    
    theta(0)=0
    do i = 0,n-1
        theta(i+1)=theta(i)+dt*(r(i)**2)
    end do
    
    !escrivim els resultats a un .txt desnormalitzant i en cartesianes
    open(unit=12, file="orbita.txt", status="replace")
        do i= 0, n
            write(12, '(F30.15, 1X, F30.15)') r(i)*alpha*COS(theta(i)), r(i)*alpha*SIN(theta(i))
        end do
    close(12)
    
    END PROGRAM solterra