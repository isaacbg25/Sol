program energia

  implicit none(type,external)
  integer, parameter :: dp = selected_real_kind(15, 300), punts_per_dia = 1440, n = 365 * punts_per_dia, num_cols = 11
  integer :: i, dia, j
  real(kind=dp) :: potencies(n, num_cols), energies(365, num_cols), h, integral, dummy
  !Definim el paràmetre alpha resultant de la normalització
  real(kind=dp), parameter :: P_0 = 400.0_dp, t_0 = 86400.0_dp

  !Llegim les dades del fitxer "any_pot_dades.txt" ignorant la primera columna i les afegim a la matriu potencies(,)
  open(unit=10, file='any_pot_dades.txt', status='old', action='read')
  do i = 1, n
      read(10, *) dummy, (potencies(i, j), j = 1, num_cols)
  end do
  close(10)

  !Normalitzem la potència
  potencies = potencies / P_0

  !Definim l'interval temporal normalitzat
  h = 1.0_dp / real(punts_per_dia - 1, kind=dp)

  !Claculem numèricament l'energia utilitzant Simpson 1/3 per diferents combinacions d'angles beta i gamma
  !i guardem els resultats a la matriu energies(,)
  do j = 1, num_cols
      do dia = 1, 365
          integral = 0.0_dp

          integral = potencies((dia - 1) * punts_per_dia + 1, j) + potencies(dia * punts_per_dia, j)

          do i = (dia - 1) * punts_per_dia + 2, dia * punts_per_dia - 1, 2
              integral = integral + 4.0_dp * potencies(i, j)
          end do

          do i = (dia - 1) * punts_per_dia + 3, dia * punts_per_dia - 2, 2
              integral = integral + 2.0_dp * potencies(i, j)
          end do

          integral = integral * (h / 3.0_dp)

          energies(dia, j) = integral
      end do
  end do

  !Guardem els resultats desnormalitzats al fitxer "energies.txt"
  open(unit=20, file='energies.txt', status='unknown', action='write')
  do dia = 1, 365
      write(20, *) dia, ((energies(dia, j) * P_0*t_0) / (3.6 * (10**6)), j = 1, num_cols)
  end do
  close(20)

end program energia
