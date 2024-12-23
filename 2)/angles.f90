program angles
implicit none
real, parameter :: pi = 3.1416
real, parameter :: r_t = 1,r_s_mod = 5, betha=(23.44*2*pi)/360, w=(2*pi)/24, alpha=(41.468*2*pi)/360, w_dies=(2*pi)/365
integer,parameter :: hours=24, n=3, days= 365
integer :: i,j
real :: r_0(n), phi, r(n), rho(n), r_s(n), r_phi(n,n), r_pla(n), r_betha(n,n), a_3(0:24,0:365), a_5(0:24,0:365)
real :: a_4(0:24,0:365), a_7(0:24,0:365)
r_betha = reshape([1.0,0.0,0.0,0.0,cos(betha),sin(betha), 0.0, -sin(betha), cos(betha)], shape(r_betha))


do j=0,365
    r_s=r_s_mod*(/cos(w_dies*j),sin(w_dies*j),0.0/)
    do i= 0,24 ! days
        r_0=r_t*(/ cos(alpha)*cos(w*i), cos(alpha)*sin(w*i), sin(alpha) /)
        r_pla = matmul(r_0, r_betha)
        !r_s= r_s_mod*(/0.0,-1.0,0.0/)
        rho =r_s + r_pla
        a_3(i,j)=acos((rho(1)*r_s(1)+rho(2)*r_s(2))/(sqrt(rho(1)**2+rho(2)**2)*sqrt(r_s(1)**2+r_s(2)**2)))
        if (r_pla(1)>0) then
            if (i==1) then
            a_5(i,j) = pi-acos((r_s(1)*r_pla(1)+r_s(2)*r_pla(2))/(sqrt(r_s(1)**2+r_s(2)**2)*sqrt(r_pla(1)**2+r_pla(2)**2)))
            else
            a_5(i,j) = (pi-acos((r_s(1)*r_pla(1)+r_s(2)*r_pla(2))/(sqrt(r_s(1)**2+r_s(2)**2)*sqrt(r_pla(1)**2+r_pla(2)**2))))
            end if
        else 
            if (i==1) then
            a_5(i,j) = -(pi-acos((r_s(1)*r_pla(1)+r_s(2)*r_pla(2))/(sqrt(r_s(1)**2+r_s(2)**2)*sqrt(r_pla(1)**2+r_pla(2)**2))))
            else
            a_5(i,j) = -(pi-acos((r_s(1)*r_pla(1)+r_s(2)*r_pla(2))/(sqrt(r_s(1)**2+r_s(2)**2)*sqrt(r_pla(1)**2+r_pla(2)**2))))
            end if
        end if
        a_4(i,j)=a_5(i,j)+a_3(i,j)
        a_7(i,j) = -pi/2 + acos((r_pla(1)*rho(1)+r_pla(2)*rho(2)+r_pla(3)*rho(3))/(sqrt(sum(r_pla**2))*sqrt(sum(rho**2))))
    end do


    open(unit=11, file = 'angles.txt', status='old', position='append')
        do i = 0,24
            write(11,*) a_4(i,j)*(180/pi), a_7(i,j)*(180/pi)
        end do
        write(11,*)""
    close(11)
end do
!contains    
    !function ang(x,y) result(angle)
    !implicit ::

    !if catx>caty
end program angles

