program cel
implicit none
real, parameter :: pi = 3.1416
real, parameter :: r_t = 1,r_s_mod = 5, betha=(23.44*2*pi)/360, w=(2*pi)/24, alpha=(41.468*2*pi)/360, w_dies=(2*pi)/365
integer,parameter :: hours=24, n=3, days= 365
integer :: i,j
real :: r_0(n), phi, r(n), rho(n), r_s(n), r_phi(n,n), r_pla(n), r_betha(n,n)

r_0=r_t*(/ 0.0, cos(alpha), sin(alpha) /)
phi = betha !acos(r_0(3)/r_t)

r_betha = reshape([1.0,0.0,0.0,0.0,cos(-1*betha),sin(-1*betha), 0.0, -sin(-1*betha), cos(-1*betha)], shape(r_betha))
r=[0.0,0.0,-r_t]

do i= 1,365, 90! days
    r_0=r_t*(/ cos(alpha)*cos(w_dies*i), cos(alpha)*sin(w_dies*i), sin(alpha) /)
    r_pla = matmul(r_betha,r_0)
    phi=acos(r_pla(3)/r_t)
    r_phi = reshape([1.0,0.0,0.0,0.0,cos(phi),sin(phi), 0.0, -sin(phi), cos(phi)], shape(r_phi))

    do j = 0,hours
        r_s=matmul([r_s_mod*cos(w*j), r_s_mod*sin(w*j),0.0],r_phi)
        rho = r + r_s
        open(unit=11, file = 'data.txt', status='old', position='append')
            write(11,*) rho
        close(11)

    end do
end do

end program cel
