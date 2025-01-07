program opt

  implicit none(type, external)
  integer, parameter :: dp = selected_real_kind(15, 300), punts_per_dia = 1440, n = 365 * punts_per_dia
  real(kind=dp), dimension(n) :: integrand_a, integrand_b, integrand_c
  real(kind=dp) :: a, b, c, h, gamma, beta, pi = 2.0*acos(0.0)
  real(kind=dp), parameter :: t_0 = 86400.0_dp, P_0 = 400.0_dp
  integer :: i

  ! Obrim l'arxiu que conté els 3 integrands a cada minut de l'any
  open(unit=10, file='../3)/integrands.txt', status='old', action='read')
  do i = 1, n
    read(10, *) integrand_a(i), integrand_b(i), integrand_c(i)
  end do
  close(10)

  ! Normalitzem els integrands
  integrand_a = integrand_a / (P_0*t_0)
  integrand_b = integrand_b / (P_0*t_0)
  integrand_c = integrand_c / (P_0*t_0)

  ! Definim l'interval temporal normalitzat
  h = 1.0_dp / real(punts_per_dia-1, kind=dp)

  ! Claculem numèricament les integrals a, b i c utilitzant Simpson 1/3 
  a = integrand_a(1) + integrand_a(n)
  b = integrand_b(1) + integrand_b(n)
  c = integrand_c(1) + integrand_c(n)

  do i = 2, n - 1, 2
    a = a + 4.0_dp * integrand_a(i)
    b = b + 4.0_dp * integrand_b(i)
    c = c + 4.0_dp * integrand_c(i)
  end do

  do i = 3, n - 2, 2
    a = a + 2.0_dp * integrand_a(i)
    b = b + 2.0_dp * integrand_b(i)
    c = c + 2.0_dp * integrand_c(i)
  end do

  a = (a * h / 3.0_dp)
  b = (b * h / 3.0_dp)
  c = (c * h / 3.0_dp)

  ! Calculem l'angle gamma
  gamma = atan(c/b)

  ! Calculem l'angle beta
  beta = atan((b*cos(gamma)+c*sin(gamma))/a)

  ! Imprimim a la terminal els resultats amb un decimal i en graus
  write(*,'(A,F6.1,A)') "L'angle gamma òptim és:", gamma*180.0_dp/pi, "°"
  write(*,'(A,F6.1,A)') "L'angle beta òptim és:", beta*180.0_dp/pi, "°"

end program opt

  
