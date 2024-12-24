program angles
implicit none
real, parameter :: pi = 3.1416
real, parameter :: r_t = 6.378e6,r_s_mod = 1.5e11, betha=(23.44*2*pi)/360, w_rot=(2*pi)/96
real, parameter ::alpha=(41.468*2*pi)/360, w_trans=(2*pi)/35040
integer,parameter :: hours=96, n=3, days= 365
integer :: i,j,k, valors(4)
character(len=100):: nom_fitxer
real :: r_0(n), phi, r(n), rho(n), r_s(n), r_phi(n,n), r_pla(n), r_betha(n,n)
real ::  a_v(0:96,0:365), a_h(0:96,0:365), rho_0(n), r_betha_inv(n,n)

r_betha = reshape([1.0,0.0,0.0,0.0,cos(betha),sin(betha), 0.0, -sin(betha), cos(betha)], shape(r_betha))
r_betha_inv = reshape([1.0,0.0,0.0,0.0,cos(-1*betha),sin(-1*betha), 0.0, -sin(-1*betha), cos(-1*betha)], shape(r_betha))
valors = (/0,91,183,274/)

do j=0,365
    
    do i= 0,96
        r_s=r_s_mod*(/cos(-w_trans*(96*j+i)),sin(-w_trans*(j*96+i)),0.0/)
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

end program angles

