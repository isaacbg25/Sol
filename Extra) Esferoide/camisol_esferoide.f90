program cami_sol_esferoide
    implicit none
    real, parameter :: pi = 3.1416
    integer,parameter ::   minuts=1439, n=3, dies= 364 !minuts al dia i dies a l'any començant pel zero
    real, parameter :: r_tes = 6.378e6, betha=(23.44*2*pi)/360, w_rot=(2*pi)/minuts
    real, parameter ::alpha_es=(41.468*2*pi)/360, w_trans=(2*pi)/35040
    !introduim parametres esferoide
    real, parameter :: a = 6378136.6, b = 6356751.9 !semieix major (radi equatorial), semieix menor (radi polar)
    real, parameter :: alpha_Tel = atan((b/a)*tan(alpha_es))
    real, parameter :: num_frac = ((((a**2)*cos(alpha_es))**2)+(((b**2)*sin(alpha_es))**2))
    real, parameter :: denom_frac = (((a*cos(alpha_es))**2)+((b*sin(alpha_es))**2))
    real, parameter :: r_Tel = sqrt(num_frac/denom_frac)
    integer :: i,j,k, valors(4)
    character(len=100):: nom_fitxer
    real :: phi
    real(8) :: r_0(n), rho_0(n), rho_pla(n), r_s_pla(n), r_pla(n), rho_terra(n), r_terra(n), r_phi(n,n)
    real(8) ::  a_v(0:minuts,0:dies), a_h(0:minuts,0:dies),r_betha(n,n), r_betha_inv(n,n), r_gamma(n,n), r_gamma_inv(n,n)
    real(8) :: r_s_terra(n), r_s_0(n), mod_rho(0:minuts,0:dies)
    
    
    !---------------------trobar el radi  i l'angle de l'orbita per cada temps (apartat 1)---------------------------------------------------
    
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    real(kind=DP), parameter :: m_t=5.972168e24, m_s=1.9885e30, g=6.67430e-11, r_0_t=147095000e3 , v_0=30.29e3, t_n = 365*3600*24 !v_0 velocitat al periheli
    integer, parameter :: n_temps=dies*minuts !calculem cada 15min durant tot l'any
    real(kind=DP), parameter :: L=r_0_t*v_0*m_t
    real(kind=DP) :: k_1,c_1
    real(kind=DP), parameter :: betha_par=(L/m_t)**2, k_norm=g*m_s, alpha_norm=betha_par/k_norm
    integer :: t 
    real(kind=DP), parameter :: v_=k_norm/(betha_par)**(0.5), t_n_norm = t_n*v_/alpha_norm, dt=t_n_norm/n_temps
    real(kind=DP) :: gamma,v(0:minuts,0:dies),r(0:minuts,0:dies),theta(0:minuts,0:dies), r_xy(0:minuts,0:dies,3) !x=minuts d'hora y=dia z=coordenada
    v(0,0)=0!ci velocitat en r(periheli)
    r(0,0)=r_0_t/alpha_norm!ci radi inicial (periheli)
    
    do j= 0,dies !dies
        do i = 0,minuts-1 !minuts d'hora en cada dia
            k_1 = -(1/r(i,j)**2) + (1/(r(i,j)**3))
            v(i+1,j)=v(i,j)+k_1*dt
            c_1 = v(i,j)
            r(i+1,j)= r(i,j)+dt*c_1 
        end do
        
        theta(0,0)=0 !ci per l'angle
        do i = 0,minuts-1
            theta(i+1,j)=theta(i,j)+dt/(r(i,j)**2)
        end do
        if (j==dies) then
        exit
        end if
        theta(0,j+1)=theta(minuts,j)+dt/(r(minuts,j)**2)
        k_1 = -(1/r(minuts,j)**2) + ((1)/(r(minuts,j)**3))
        v(0,j+1)=v(minuts,j)+K_1*dt!(dt/6)*(k_1+2*K_2+2*k_3+k_4)
        c_1 = v(minuts,j)
        r(0,j+1)= r(minuts,j)+dt*c_1 !(dt/6)*(c_1+2*c_2+2*c_3+c_4)
    end do
    !vectors sol-centre terra al llarg de l'any
    do j = 0, dies
        do i = 0,minuts 
            r_xy(i,j,1)=r(i,j)*alpha_norm*SIN(theta(i,j))
            r_xy(i,j,2)=r(i,j)*alpha_norm*COS(theta(i,j))
            r_xy(i,j,3)=0.0
        end do
    end do
    
    !calcular l'angle gamma entre el vector posició al solstici d'hivern i al periheli
    gamma = (r_xy(0,0,1)*r_xy(0,350,1)+r_xy(0,0,2)*r_xy(0,350,2)+r_xy(0,0,3)*r_xy(0,350,3))
    gamma = acos(gamma/(sqrt(sum(r_xy(0,0,:)**2))*sqrt(sum(r_xy(0,351,:)**2))))
    
    
    !----------------------Calcular els angles del sol amb la nostra casa--------------------------------------------------------
    
!matrius per fer tots els canvis de coordenades
r_betha = reshape([1.0,0.0,0.0,0.0,cos(betha),-sin(betha), 0.0, sin(betha), cos(betha)], shape(r_betha))!sistema amb z orientat amb l'eix de la terra (sist ref 0) a sistema amb z perpendicular al pla d'orbita (sist ref terra)
r_betha_inv = reshape([1.0,0.0,0.0,0.0,cos(-1*betha),-sin(-1*betha), 0.0, sin(-1*betha), cos(-1*betha)], shape(r_betha))
valors = (/350,75,168,262/)
r_gamma= reshape([cos(gamma),sin(gamma), 0.0_8, -sin(gamma), cos(gamma),0.0_8,0.0_8,0.0_8,1.0_8], shape(r_gamma))!sistema terra amb l'angle de l'eix de la terra al pla yz(sist ref terra) a sistema amb eixos x i y orientats amb l'lel·lipse de l'orbita (sist ref pla)
r_gamma_inv=reshape([cos(-1*gamma),sin(-1*gamma), 0.0_8, -sin(-1*gamma), cos(-1*gamma),0.0_8,0.0_8,0.0_8,1.0_8], shape(r_gamma))

!fem el calcul per tres epoques de l'any començant a j=0 que coincideix amb el 3 de gener (terra al periheli)
!hem de separar el calcul per sumar els angles correctament (separem l'any entre solsticis(t=168 i t=250) però comencem al periheli, per tant 3 divisions)
do j=0,167

    do i= 0,minuts
        if (i==0) then
            r_s_pla=(/r(i,j)*alpha_norm*SIN(theta(i,j)),r(i,j)*alpha_norm*COS(theta(i,j)),0.0_8/)!r_s_pla és el vector entre el centre de la terra i el del sol
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2), cos(alpha_Tel)*sin(w_rot*i+pi/2), sin(alpha_Tel) /)! r_0 és el vector entre el centre de la terra i la casa
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            rho_terra = matmul(r_gamma_inv,rho_pla)
            rho_0= matmul( r_betha_inv,rho_terra)
            phi=-acos((rho_0(1)*r_0(1)+rho_0(2)*r_0(2))/(sqrt(rho_0(1)**2+rho_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            !r_phi =reshape([cos(phi),-sin(phi), 0.0_8, sin(phi), cos(phi),0.0_8,0.0_8,0.0_8,1.0_8], shape(r_gamma))
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2+phi), cos(alpha_Tel)*sin(w_rot*i+pi/2+phi), sin(alpha_Tel) /)
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            r_s_terra = matmul(r_gamma_inv,rho_pla)
            r_s_0= matmul( r_betha_inv,r_s_terra)
            a_h(i,j) =pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            a_v(i,j)=((r_pla(1)*rho_pla(1)+r_pla(2)*rho_pla(2)+r_pla(3)*rho_pla(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho_pla**2))))
            a_v(i,j)=acos(a_v(i,j)) - pi/2
        else 
            r_s_pla=(/r(i,j)*alpha_norm*SIN(theta(i,j)),r(i,j)*alpha_norm*COS(theta(i,j)),0.0_8/)!r_s_pla és el vector entre el centre de la terra i el del sol
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2+phi), cos(alpha_Tel)*sin(w_rot*i+pi/2+phi), sin(alpha_Tel) /)
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            r_s_terra = matmul(r_gamma_inv,rho_pla)
            r_s_0= matmul( r_betha_inv,r_s_terra)
            if (i<720) then
                a_h(i,j)=pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            else
                a_h(i,j)=-(pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2))))
            end if
            a_v(i,j)=((r_pla(1)*rho_pla(1)+r_pla(2)*rho_pla(2)+r_pla(3)*rho_pla(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho_pla**2))))
            a_v(i,j)=acos(a_v(i,j)) - pi/2    
        end if
        mod_rho(i,j)=sqrt(sum(rho_pla**2))

    end do

end do

do j=168,349

    do i= 0,minuts
        if (i==0) then
            r_s_pla=(/r(i,j)*alpha_norm*SIN(theta(i,j)),r(i,j)*alpha_norm*COS(theta(i,j)),0.0_8/)!r_s_pla és el vector entre el centre de la terra i el del sol
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2), cos(alpha_Tel)*sin(w_rot*i+pi/2), sin(alpha_Tel) /)! r_0 és el vector entre el centre de la terra i la casa
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            rho_terra = matmul(r_gamma_inv,rho_pla)
            rho_0= matmul( r_betha_inv,rho_terra)
            phi=acos((rho_0(1)*r_0(1)+rho_0(2)*r_0(2))/(sqrt(rho_0(1)**2+rho_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            !r_phi =reshape([cos(phi),-sin(phi), 0.0_8, sin(phi), cos(phi),0.0_8,0.0_8,0.0_8,1.0_8], shape(r_gamma))
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2+phi), cos(alpha_Tel)*sin(w_rot*i+pi/2+phi), sin(alpha_Tel) /)
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            r_s_terra = matmul(r_gamma_inv,rho_pla)
            r_s_0= matmul( r_betha_inv,r_s_terra)
            if (i<720) then
                a_h(i,j)=pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            else
                a_h(i,j)=-(pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2))))
            end if
            a_v(i,j)=((r_pla(1)*rho_pla(1)+r_pla(2)*rho_pla(2)+r_pla(3)*rho_pla(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho_pla**2))))
            a_v(i,j)=acos(a_v(i,j)) - pi/2    
        else 
            r_s_pla=(/r(i,j)*alpha_norm*SIN(theta(i,j)),r(i,j)*alpha_norm*COS(theta(i,j)),0.0_8/)!r_s_pla és el vector entre el centre de la terra i el del sol
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2+phi), cos(alpha_Tel)*sin(w_rot*i+pi/2+phi), sin(alpha_Tel) /)
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            r_s_terra = matmul(r_gamma_inv,rho_pla)
            r_s_0= matmul( r_betha_inv,r_s_terra)
            if (i<720) then
                a_h(i,j)=pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            else
                a_h(i,j)=-(pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2))))
            end if
            a_v(i,j)=((r_pla(1)*rho_pla(1)+r_pla(2)*rho_pla(2)+r_pla(3)*rho_pla(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho_pla**2))))
            a_v(i,j)=acos(a_v(i,j)) - pi/2       
        end if
        mod_rho(i,j)=sqrt(sum(rho_pla**2))

    end do

end do

do j=350,364

    do i= 0,minuts
        if (i==0) then
            r_s_pla=(/r(i,j)*alpha_norm*SIN(theta(i,j)),r(i,j)*alpha_norm*COS(theta(i,j)),0.0_8/)!r_s_pla és el vector entre el centre de la terra i el del sol
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2), cos(alpha_Tel)*sin(w_rot*i+pi/2), sin(alpha_Tel) /)! r_0 és el vector entre el centre de la terra i la casa
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            rho_terra = matmul(r_gamma_inv,rho_pla)
            rho_0= matmul( r_betha_inv,rho_terra)
            phi=-acos((rho_0(1)*r_0(1)+rho_0(2)*r_0(2))/(sqrt(rho_0(1)**2+rho_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            !r_phi =reshape([cos(phi),-sin(phi), 0.0_8, sin(phi), cos(phi),0.0_8,0.0_8,0.0_8,1.0_8], shape(r_gamma))
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2+phi), cos(alpha_Tel)*sin(w_rot*i+pi/2+phi), sin(alpha_Tel) /)
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            r_s_terra = matmul(r_gamma_inv,rho_pla)
            r_s_0= matmul( r_betha_inv,r_s_terra)
            if (i<720) then
                a_h(i,j)=pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            else
                a_h(i,j)=-(pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2))))
            end if
            a_v(i,j)=((r_pla(1)*rho_pla(1)+r_pla(2)*rho_pla(2)+r_pla(3)*rho_pla(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho_pla**2))))
            a_v(i,j)=acos(a_v(i,j)) - pi/2    
        else 
            r_s_pla=(/r(i,j)*alpha_norm*SIN(theta(i,j)),r(i,j)*alpha_norm*COS(theta(i,j)),0.0_8/)!r_s_pla és el vector entre el centre de la terra i el del sol
            r_0=r_Tel*(/ cos(alpha_Tel)*cos(w_rot*i+pi/2+phi), cos(alpha_Tel)*sin(w_rot*i+pi/2+phi), sin(alpha_Tel) /)
            r_terra = matmul(r_betha,r_0) !r_terra és el vector de r_0 al sistema ref terra
            r_pla = matmul(r_gamma,r_terra)
            rho_pla =r_s_pla + r_pla !rho_pla és el vector entre el sol i la casa al sistema de ref orbita
            r_s_terra = matmul(r_gamma_inv,rho_pla)
            r_s_0= matmul( r_betha_inv,r_s_terra)
            if (i<720) then
                a_h(i,j)=pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2)))
            else
                a_h(i,j)=-(pi-acos((r_s_0(1)*r_0(1)+r_s_0(2)*r_0(2))/(sqrt(r_s_0(1)**2+r_s_0(2)**2)*sqrt(r_0(1)**2+r_0(2)**2))))
            end if
            a_v(i,j)=((r_pla(1)*rho_pla(1)+r_pla(2)*rho_pla(2)+r_pla(3)*rho_pla(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho_pla**2))))
            a_v(i,j)=acos(a_v(i,j)) - pi/2      
        end if
        mod_rho(i,j)=sqrt(sum(rho_pla**2))

        

    end do



   
end do
    
    !escrivim els angles trobats a un txt
open(unit=11, file = 'angles_esferoide.txt', status='replace')
    do j=0,dies
        do k = 0,minuts 
            write(11,'(F10.4,1x,F10.4)') a_h(k,j)*(180/pi), a_v(k,j)*(180/pi) 
        end do
    end do
close(11)
    
    do i= 1, size(valors)
        write(nom_fitxer, '("eqesferoide",I0,".txt")') i
        open(unit=i, file =trim(nom_fitxer), status='replace')
            j = valors(i)
            do k = 0,minuts
                if (a_v(k,j)>0) then
                    write(i,'(F10.2,1x,F10.2)') a_h(k,j)*(180/pi), a_v(k,j)*(180/pi)
                end if
            end do
            write(i,*)""
        close(unit=i)
    end do

open(unit=12, file = 'dist_sol_esferoide.txt', status='replace')
    do j=0,dies
        do k = 0,minuts 
            write(12,*) mod_rho(k,j)
        end do
    end do
close(12)
end program cami_sol_esferoide
