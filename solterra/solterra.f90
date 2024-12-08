PROGRAM solterra
INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
real(kind=DP), parameter :: m_t=5.972168e24, m_s=1.9885e30, g=6.674302e-11, r_0=152098232e3, v_0=29.29e3, t_n = 365.25
integer, parameter :: n=300
real(kind=DP), parameter :: L=r_0*v_0*m_t
real(kind=DP) :: k_1, k_2, k_3, k_4, c_1, c_2, c_3, c_4
real(kind=DP), parameter :: betha=(L/m_t)**2, k=g*m_s, alpha=betha/k, v_=k/(betha)**(1/2), dt = (t_n*v_)/(n*alpha)
integer :: i,t 
real(kind=DP) :: v(0:n),r(0:n)
v(0)=v_0/v_
r(0)=r_0/alpha
do i = 0,n-1 !perquè ens dona en el temps i+1
    t = i*dt !crec que no cal perquè no apareix a les equs
    k_1 = (1/r(i)**2) - (1/r(i)**3)
    K_2 = (1/(r(i)+K_1*dt*0.5)**2) - (1/(r(i)+K_1*dt*0.5)**3)
    k_3 = (1/(r(i)+K_2*dt*0.5)**2) - (1/(r(i)+K_2*dt*0.5)**3)
    k_4 = (1/(r(i)+K_3*dt)**2) - (1/(r(i)+K_3*dt)**3)
    v(i+1)=v(i)+(dt/6)*(k_1+2*K_2+2*k_3+k_4)
    c_1 = v(i)
    c_2 =v(i)+0.5*dt*c_1
    c_3 = v(i)+0.5*dt*c_2
    c_4 = v(i)+dt*c_3
    r(i+1)= r(i)+ (dt/6)*(c_1+2*c_2+2*c_3+c_4)
end do
print*, r(0:20), v(0:20), r(100)
    

    











END PROGRAM solterra