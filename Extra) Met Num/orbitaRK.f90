PROGRAM orbitark

    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    real(kind=DP), parameter :: m_t=5.972168e24, m_s=1.9885e30, g=6.674302e-11, r_0=147098070e3, v_0=30.29e3, t_n = 365*3600*24
    integer, parameter :: n=3000 !discretització temporal
    real(kind=DP), parameter :: L=r_0*v_0*m_t
    real(kind=DP) :: k1_v, k2_v, k1_r, k2_r, c1_b, c2_b, c3_b , c4_b, c1_f, c2_f, c3_f , c4_f !variables que utilitzarem a RK
    real(kind=DP), parameter :: betha=(L/m_t)**2, k=g*m_s, alpha=betha/k, v_=k/(betha)**(0.5) !constants de normalització
    integer :: i
    real(kind=DP) :: v(0:n),r(0:n),theta_RK2(0:n) !velocitat radial, posició i angle RK2
    real(kind=DP) :: b(0:n),f(0:n),theta_RK4(0:n) !velocitat radial, posició i angle RK4
    real(kind=dp), parameter :: t_n_norm = t_n*v_/alpha, dt=t_n_norm/n

    ! ci de les edos en r RK2
    v(0)=0
    r(0)=r_0/alpha

    ! ci de les edos en r RK4
    b(0)=0
    f(0)=r_0/alpha

    !solució de les equacions per cada MN
    do i = 0,n-1 !perquè ens dona en el temps i+1
        !RK-2
        !primer pas
        K1_v = dt*(-(1/r(i)**2) + ((1/r(i))**3))
        K1_r = dt*v(i)

        !segon pas
        K2_v = dt*(-(1/(r(i) + K1_v*0.5)**2) + ((1/(r(i) + K1_v*0.5)**3)))
        K2_r = dt*(v(i)+ K1_r*0.5)

        !actualitzacio
        v(i+1)= v(i)+K2_v
        r(i+1)= r(i)+K2_r

        !RK-4
        C1_b = (-(1/f(i)**2) + ((1/f(i)**3)))
        C1_f = b(i)
        
        C2_b = (-(1/(f(i)+C1_b*dt*0.5)**2) + ((1/(f(i)+C1_b*dt*0.5)**3)))
        C2_f = b(i)+C1_f*dt*0.5
        
        C3_b = (-(1/(f(i)+C2_b*dt*0.5)**2) + ((1/(f(i)+C2_b*dt*0.5)**3)))
        C3_f = b(i)+C2_f*dt*0.5
        
        C4_b = (-(1/(f(i)+C3_b*dt)**2) + ((1/(f(i)+C3_b*dt)**3)))
        C4_f = b(i)+C3_f*dt

        !actualitzacio
        b(i+1)= b(i)+(dt/6)*(C1_b + 2*C2_b + 2*C3_b + C4_b)
        f(i+1)= f(i)+(dt/6)*(C1_f + 2*C2_f + 2*C3_f + C4_f)

    end do

    theta_RK2(0)=0
    theta_RK4(0)=0
    do i = 0,n-1
        theta_RK2(i+1)=theta_RK2(i)+dt*(r(i)**2)
        theta_RK4(i+1)=theta_RK4(i)+dt*(f(i)**2)
    end do

    !escrivim els resultats a un .txt desnormalitzant i en cartesianes
    open(unit=12, file="orbitaRK2.txt", status="replace")
        do i= 0, n
            write(12, '(F30.15, 1X, F30.15)') r(i)*alpha*COS(theta_RK2(i)), r(i)*alpha*SIN(theta_RK2(i))
        end do
    close(12)

        open(unit=12, file="orbitaRK4.txt", status="replace")
        do i= 0, n
            write(12, '(F30.15, 1X, F30.15)') f(i)*alpha*COS(theta_RK4(i)), f(i)*alpha*SIN(theta_RK4(i))
        end do
    close(12)

    END PROGRAM orbitark
