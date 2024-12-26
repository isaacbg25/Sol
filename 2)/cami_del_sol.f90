program cami_del_sol
implicit none
real, parameter :: pi = 3.1416
integer,parameter ::   quarts=96, n=3, dies= 365
real, parameter :: r_t = 6.378e6, betha=(23.44*2*pi)/360, w_rot=(2*pi)/quarts
real, parameter ::alpha=(41.468*2*pi)/360, w_trans=(2*pi)/35040
integer :: i,j,k, valors(4)
character(len=100):: nom_fitxer
real(8) :: r_0(n), rho_0(n), rho_pla(n), r_s_pla(n), r_pla(n), rho_terra(n), r_terra(n)
real(8) ::  a_v(0:quarts,0:dies), a_h(0:quarts,0:dies),r_betha(n,n), r_betha_inv(n,n), r_gamma(n,n), r_gamma_inv(n,n)


!---------------------trobar el radi  i l'angle de l'orbita per cada temps (apartat 1)---------------------------------------------------

INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
real(kind=DP), parameter :: m_t=5.972168e24, m_s=1.9885e30, g=6.674302e-11, r_0_t=147098290e3 , v_0=30.29e3, t_n = dies*3600*24 !v_0 velocitat al periheli
integer, parameter :: n_temps=dies*quarts !calculem cada 15min durant tot l'any
real(kind=DP), parameter :: L=r_0_t*v_0*m_t
real(kind=DP) :: k_1,c_1
real(kind=DP), parameter :: betha_par=(L/m_t)**2, k_norm=g*m_s, alpha_norm=betha_par/k_norm
integer :: t 
real(kind=DP), parameter :: v_=k_norm/(betha_par)**(1/2), dt = t_n/n_temps!(t_n*v_)/(n_temps*alpha_norm)
real(kind=DP) :: gamma,v(0:quarts,0:dies),r(0:quarts,0:dies),theta(0:quarts,0:dies), r_xy(0:quarts,0:dies,3) !x=quarts d'hora y=dia z=coordenada
v(0,0)=0!ci velocitat en r(periheli)
r(0,0)=r_0_t!ci radi inicial (periheli)

do j= 0,dies !dies
    do i = 0,quarts-1 !quarts d'hora en cada dia
        k_1 = -(g*m_s/r(i,j)**2) + ((l**2)/(m_t**(2)*r(i,j)**3))
        v(i+1,j)=v(i,j)+K_1*dt
        c_1 = v(i,j)
        r(i+1,j)= r(i,j)+dt*c_1 
    end do
    
    theta(0,0)=0 !ci per l'angle
    do i = 0,quarts-1
        theta(i+1,j)=theta(i,j)+dt*L/(m_t*r(i,j)**2)
    end do
    if (j==dies) then
    exit
    end if
    theta(0,j+1)=theta(quarts,j)+dt*L/(m_t*r(quarts,j)**2)
    k_1 = -(g*m_s/r(quarts,j)**2) + ((l**2)/(m_t**(2)*r(quarts,j)**3))
    v(0,j+1)=v(quarts,j)+K_1*dt!(dt/6)*(k_1+2*K_2+2*k_3+k_4)
    c_1 = v(quarts,j)
    r(0,j+1)= r(quarts,j)+dt*c_1 !(dt/6)*(c_1+2*c_2+2*c_3+c_4)
end do

do j = 0, dies
    do i = 0,quarts 
        r_xy(i,j,1)=r(i,j)*COS(theta(i,j))
        r_xy(i,j,2)=r(i,j)*SIN(theta(i,j))
        r_xy(i,j,3)=0.0
    end do
end do

!calcular l'angle gamma entre el vector posició al solstici d'hivern i al periheli
gamma = (r_xy(0,0,1)*r_xy(0,352,1)+r_xy(0,0,2)*r_xy(0,352,2)+r_xy(0,0,3)*r_xy(0,352,3))
gamma = acos(gamma/(sqrt(sum(r_xy(0,0,:)**2))*sqrt(sum(r_xy(0,352,:)**2))))


!----------------------Calcular els angles del sol amb la nostra casa--------------------------------------------------------

!matrius per fer tots els canvis de coordenades
r_betha = reshape([1.0,0.0,0.0,0.0,cos(betha),sin(betha), 0.0, -sin(betha), cos(betha)], shape(r_betha))!sistema amb z orientat amb l'eix de la terra (sist ref 0) a sistema amb z perpendicular al pla d'orbita (sist ref terra)
r_betha_inv = reshape([1.0,0.0,0.0,0.0,cos(-1*betha),sin(-1*betha), 0.0, -sin(-1*betha), cos(-1*betha)], shape(r_betha))
valors = (/352,80,172,266/)
r_gamma = reshape([cos(gamma),sin(gamma), 0.0_8, -sin(gamma), cos(gamma),0.0_8,0.0_8,0.0_8,1.0_8], shape(r_gamma))!sistema terra amb l'angle de l'eix de la terra al pla yz(sist ref terra) a sistema amb eixos x i y orientats amb l'lel·lipse de l'orbita (sist ref pla)
r_gamma_inv = reshape([cos(-1*gamma),sin(-1*gamma), 0.0_8, -sin(-1*gamma), cos(-1*gamma),0.0_8,0.0_8,0.0_8,1.0_8], shape(r_gamma))

do j=0,dies
    
    do i= 0,quarts
        r_s_pla=(/r(i,j)*COS(theta(i,j)), r(i,j)*SIN(theta(i,j)),0.0_8/)!r_s_pla és el vector entre el centre de la terra i el del sol
        r_0=r_t*(/ cos(alpha)*cos(w_rot*i), cos(alpha)*sin(w_rot*i), sin(alpha) /)! r_0 és el vector entre el centre de la terra i la casa
        r_terra = matmul(r_0, r_betha) !r_terra és el vector de r_0 al sistema ref terra
        r_pla = matmul(r_terra,r_gamma)
        rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
        rho_terra = matmul(rho_pla,r_gamma_inv)
        rho_0= matmul(rho_terra, r_betha_inv)
        a_h(i,j) =pi-acos((rho_0(1)*r_0(1)+rho_0(2)*r_0(2))/(sqrt(rho_0(1)**2+rho_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
        a_v(i,j) = acos((r_pla(1)*rho_pla(1)+r_pla(2)*rho_pla(2)+r_pla(3)*rho_pla(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho_pla**2))))
        a_v(i,j)=a_v(i,j) - pi/2
    end do

!escrivim els angles trobats a un txt
    open(unit=11, file = 'angles.txt', status='old', position='append')
        do k = 0,quarts
            if (a_v(k,j)>0) then
            write(11,*) a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
            end if
        end do
        do k = 0,quarts
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
        do k = 0,quarts
            if (a_v(k,j)>0) then
                write(i,'(F10.2,1x,F10.2)') a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
            end if
        end do
        do k = 0,quarts
            if (a_v(k,j)>0) then
                write(i,'(F10.2,1x,F10.2)') -1*a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
            end if
        end do
        write(i,*)""
    close(unit=i)
end do

end program cami_del_sol