program angles
implicit none
real, parameter :: pi = 3.1416
real, parameter :: r_t = 6.378e6,r_s_mod = 1.5e11, betha=(23.44*2*pi)/360, w_rot=(2*pi)/96
real, parameter ::alpha=(41.468*2*pi)/360, w_trans=(2*pi)/35040
integer,parameter :: hours=96, n=3, days= 365
integer :: i,j,k, valors(4)
character(len=100):: nom_fitxer
real :: r_0(n), phi, rho(n), r_s(n), r_phi(n,n), r_pla(n), r_betha(n,n)
real ::  a_v(0:96,0:365), a_h(0:96,0:365), rho_0(n), r_betha_inv(n,n), r_gamma(n,n)

!trobar el radi  i l'angle de l'orbita per cada temps (apartat 1)
INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
real(kind=DP), parameter :: m_t=5.972168e24, m_s=1.9885e30, g=6.674302e-11, r_0_t=147098290e3 , v_0=29.29e3, t_n = 365*3600*24
integer, parameter :: n_temps=365*96 !calculem cada 15min durant tot l'any
real(kind=DP), parameter :: L=r_0_t*v_0*m_t
real(kind=DP) :: k_1, k_2, k_3, k_4, c_1, c_2, c_3, c_4
real(kind=DP), parameter :: betha_par=(L/m_t)**2, k_norm=g*m_s, alpha_norm=betha_par/k_norm
integer :: t 
real(kind=DP), parameter :: v_=k_norm/(betha_par)**(1/2), dt = t_n/n_temps!(t_n*v_)/(n_temps*alpha_norm)
real(kind=DP) :: gamma,v(0:96,0:365),r(0:96,0:365),theta(0:96,0:365), r_xy(0:96,0:365,3) !x=quardores y=dia z=coordenada
v(0,0)=0
r(0,0)=r_0_t
do j= 0,365
    do i = 0,95 !perquè ens dona en el temps i+1
        t = i*dt !crec que no cal perquè no apareix a les equs
        k_1 = -(g*m_s/r(i,j)**2) + ((l**2)/(m_t**(2)*r(i,j)**3))
        !K_2 = (1/(r(i,j)+K_1*dt*0.5)**2) + (1/(r(i,j)+K_1*dt*0.5)**3)
        !k_3 = (1/(r(i,j)+K_2*dt*0.5)**2) + (1/(r(i,j)+K_2*dt*0.5)**3)
        !k_4 = (1/(r(i,j)+K_3*dt)**2) + (1/(r(i,j)+K_3*dt)**3)
        v(i+1,j)=v(i,j)+K_1*dt!(dt/6)*(k_1+2*K_2+2*k_3+k_4)
        c_1 = v(i,j)
        !c_2 =v(i)+0.5*dt*c_1
        !c_3 = v(i)+0.5*dt*c_2
        !c_4 = v(i)+dt*c_3
        r(i+1,j)= r(i,j)+dt*c_1 !(dt/6)*(c_1+2*c_2+2*c_3+c_4)
    end do

    theta(0,0)=0
    do i = 0,95
        theta(i+1,j)=theta(i,j)+dt*L/(m_t*r(i,j)**2)
    end do
end do

do j = 0, 365
    do i = 0,96 
        r_xy(i,j,1)=r(i,j)*COS(theta(i,j))
        r_xy(i,j,2)=r(i,j)*SIN(theta(i,j))
        r_xy(i,j,3)=0.0
    end do
end do

!calcular l'angle gamma entre el solstici d'hivern i el 3 de gener(periheli)
gamma = (r_xy(0,0,1)*r_xy(0,352,1)+r_xy(0,0,2)*r_xy(0,352,2)+r_xy(0,0,3)*r_xy(0,352,3))
gamma = acos(gamma/(sqrt(sum(r_xy(0,0,:)**2))*sqrt(sum(r_xy(0,352,:)**2))))

r_betha = reshape([1.0,0.0,0.0,0.0,cos(betha),sin(betha), 0.0, -sin(betha), cos(betha)], shape(r_betha))
r_betha_inv = reshape([1.0,0.0,0.0,0.0,cos(-1*betha),sin(-1*betha), 0.0, -sin(-1*betha), cos(-1*betha)], shape(r_betha))
valors = (/0,91,183,274/)
r_gamma = reshape([1.0_8,0.0_8,0.0_8,0.0_8,cos(gamma),sin(gamma), 0.0_8, -sin(gamma), cos(gamma)], shape(r_gamma))

do j=0,365
    
    do i= 0,96
        !r_s=r_s_mod*(/cos(-w_trans*(96*j+i)),sin(-w_trans*(j*96+i)),0.0/)
        r_s=(/r(i,j)*COS(theta(i,j)), r(i,j)*SIN(theta(i,j)),0.0_8/)
        r_0=r_t*(/ cos(alpha)*cos(w_rot*i), cos(alpha)*sin(w_rot*i), sin(alpha) /)
        r_pla = matmul(r_0, r_betha)
        rho =r_s + r_pla
        rho_0= matmul(rho, r_betha_inv)
        a_h(i,j) =pi-acos((rho_0(1)*r_0(1)+rho_0(2)*r_0(2))/(sqrt(rho_0(1)**2+rho_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
        a_v(i,j) = -pi/2 + acos((r_pla(1)*rho(1)+r_pla(2)*rho(2)+r_pla(3)*rho(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho**2))))
    end do


    open(unit=11, file = 'angles.txt', status='old', position='append')
        do k = 0,96
            if (a_v(k,j)>0) then
            write(11,*) a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
            end if
        end do
        do k = 0,96
            if (a_v(k,j)>0) then
            write(11,*) -1*a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
            end if
        end do
        write(11,*)""
    close(11)
end do


do i= 1, size(valors)
    write(nom_fitxer, '("equinoci",I0,".txt")') i
    open(unit=i, file =trim(nom_fitxer), status='replace')
        j = valors(i)
        do k = 0,96
            if (a_v(k,j)>0) then
                write(i,*) a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
            end if
        end do
        do k = 0,96
            if (a_v(k,j)>0) then
                write(i,*) -1*a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
            end if
        end do
        write(i,*)""
    close(unit=i)
end do
print*, gamma*180/pi, r(2,1), dt
end program angles

