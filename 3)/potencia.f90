program potenc
    implicit none
 
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,300)
 
    ! Declaració de variables
    real(KIND=DP) :: thetaz, gamma_s, costheta, d, potencia, diacoma, hora, pot_caracteristica, r, gamma, beta, d_caracteristica
    real(KIND=DP) :: intensitat, inta, intb, intc
    integer :: i, j, dia, total_intervals, g, b
    real(KIND=DP), dimension(:,:), allocatable :: angles_horitz
    real(KIND=DP), dimension(:,:), allocatable :: angles_inc
    real(KIND=DP), dimension(:), allocatable :: distancies
    real(KIND=DP), dimension(3) :: betas
    real(KIND=DP), dimension(5) :: gammas
 
    ! Declaració de paràmetres
    ! Constants numèriques genèriques
    real(KIND=DP), parameter :: pi = 2.0 * ACOS(0.0)
    real(KIND=DP), parameter :: sigmaSB = 5.6704e-8 ! Constant de Stefan-Boltzmann (W/m^2 K^4)
    real(KIND=DP), parameter :: alpha = 1.0 - 0.294 ! Percentatge de llum solar (%) aborbida per la Terra
    ! 0.294 és el "Bond albedo" de "Earth Fact Sheet"
    real(KIND=DP), parameter :: RS = 695700.0e3     ! Radi del Sol (m) ("Volumetric mean radius" de "Sun Fact Sheet")
    real(KIND=DP), parameter :: TS = 5772.0         ! Temperatura superfície del Sol (K) ("Effective temperature" de "Sun Fact Sheet")
 
    ! Discretització temporal
    real(KIND=DP), parameter :: num_intervals_per_hora = 60    ! Predeterminat segons el nombre d'angles per dia
    integer, parameter :: total_dies = 365                     ! Nombre de dies en un any
 
    ! Constants numèriques associades a la placa fotovoltàica
    real(KIND=DP), parameter :: A_placa = 2.0          ! Àrea de la placa (m^2)
    real(KIND=DP), parameter :: pot_el_max = 400.0  ! Màxim d'electricitat que genera (W) quan
    real(KIND=DP), parameter :: irr_max = 1.0e3     ! irradiació solar (W) incideix perpendicularment sobre la placa
    r = pot_el_max / (irr_max * A_placa)   ! Rendiment (%)
    ! Angles de la inclinació i orientació de la placa que volem provar; introduir manualment
    betas = [0.0, 35.0, 44.5]                            ! Inclinacions (°)
    gammas = [-60.0, -15.0, 0.0, 15.0, 60.0]          ! Angles azimutals (°) 
 
    ! Constants caracterítiques del sistema que utilitzarem per a normalitzar les variables
    pot_caracteristica = pot_el_max
    d_caracteristica = sqrt((r * alpha * (RS**2) * sigmaSB * (TS**4) * A_placa) / (pot_caracteristica))
 
    total_intervals = num_intervals_per_hora * 24   ! Nombre d'intervals per dia
 
    ! Allocate les matrius dels angles i el vector de distàncies
    allocate(angles_horitz(total_dies, total_intervals))
    allocate(angles_inc(total_dies, total_intervals))
    allocate(distancies(total_dies * total_intervals))
 
    ! A continuació hi afegim les dades corresponents dels fitxers .txt obtinguts prèviament
 
    ! Obrim el fitxer de distàncies i guardem en "distancies" les dades llegides
    open(unit=11, file='../2)/dist_sol.txt', status='unknown')
    i = 1
    do while (.true.)
        read(11, *, iostat=j) distancies(i)
        if (j /= 0) exit
        if (distancies(i) <= 0.0) then
            print *, "Distància invàlida detectada en l'índex ", i
            stop
        end if
        i = i + 1
    end do
    close(11)
 
    ! Obrim el fitxer dels angles i assignem els valors a les matrius
    open(unit=10, file='../2)/angles.txt', status='unknown')
    dia = 1
    j = 1
    do while (.true.)
       read(10, *, iostat=i) angles_horitz(dia, j), angles_inc(dia, j)
       if (i /= 0) exit
       j = j + 1
       if (j > total_intervals) then ! Dia ha acabat, passem al següent i tornem a començar
          j = 1
          dia = dia + 1
       end if
    end do
    close(10)
 
    ! Obrim fitxers .txt de sortida on escriurem els resultats de la potència
    open(unit=12, file='any_pot_dades.txt', status='unknown')
    open(unit=13, file='dia1_pot_dades.txt', status='unknown')
    open(unit=14, file='dia170_pot_dades.txt', status='unknown')
    ! Calcularem també els integrands de l'integral de l'energia, que seran necessaris en el programa per a l'optimització dels angles
    open(unit=15, file='integrands.txt', status='unknown')
 
    ! Calcularem la potència total_intervals vegades cada dia, durant 365 dies
    do dia = 1, total_dies
       do j = 1, total_intervals
          ! En primer lloc, en la primera columna de cada fitxer escrivim els temps corresponents a cada potència que calcularem
          diacoma = REAL(dia) - 1.0 + REAL(j) / (num_intervals_per_hora * 24.0)  ! Temps en dies
          write(12, '(ES12.5,1X)', advance='no') diacoma
          if (dia == 1) then ! En aquest fitxer guardarem les dades d'un únic dia (el 1r)
             hora = REAL(j) / num_intervals_per_hora  ! Temps en hores
             write(13, '(ES12.5,1X)', advance='no') hora
          else if (dia == 170) then ! En aquest fitxer guardarem les dades d'un únic dia (el 170è)
            hora = REAL(j) / num_intervals_per_hora  ! Temps en hores
            write(14, '(ES12.5,1X)', advance='no') hora
          end if
          
          
         ! Calculem els integrands
         ! Com la potència és nul·la si no es compleix que thetaz >= -pi/2 .and. thetaz <= pi/2,
         ! imposem la mateixa condició sobre els integrands
          if (thetaz >= -pi/2 .and. thetaz <= pi/2) then
            intensitat = (alpha * (RS**2) * sigmaSB * (TS**4)) / & 
            ((distancies((dia - 1) * total_intervals + j)**2) * d_caracteristica**2)
            inta = r * intensitat * A_placa * COS(thetaz)
            intb = r * intensitat * A_placa * SIN(thetaz) * COS(gamma_s)
            intc = r * intensitat * A_placa * SIN(thetaz) * SIN(gamma_s)
         else
            inta = 0.0
            intb = 0.0
            intc = 0.0
         end if
         write(15, '(ES12.5, 1X, ES12.5, 1X, ES12.5)') inta, intb, intc


         ! Calcularem la potència per cada beta i gamma
         ! També registrarem el temps per a poder fer-ne un gràfic en funció del temps
          DO b = 1, SIZE(betas)
            beta = betas(b) * pi / 180.0
            DO g = 1, SIZE(gammas)
               gamma = gammas(g) * pi / 180.0
               
               ! Si beta és 0 el resultat és independent de gamma així que ho fem només per gamma=0
               IF(beta==0.0) THEN
                  gamma = 0.0
               END IF
 
                ! Seleccionem els valors dels angles del Sol en el moment en què estem, i els convertim a radians
                thetaz = ( 90.0 - angles_inc(dia, j) ) * pi / 180.0
                gamma_s = ( angles_horitz(dia, j) ) * pi / 180.0
 
                ! Calculem cos(theta)
                IF (thetaz >= -pi/2 .AND. thetaz <= pi/2) THEN
                   costheta = MAX(COS(thetaz) * COS(beta) + SIN(thetaz) * SIN(beta) * COS(gamma_s - gamma), 0.0)
                    ! Agafem només valors positius del cosinus, és a dir, imposem que -pi/2 <= theta <= pi/2 perquè en cas contrari el Sol
                    ! es trobaria darrere la placa, és a dir, que la mateixa placa bloquejaria la llum solar
                ELSE
                   costheta = 0.0 ! Perquè el Sol està per sota l'horitzó i per tant no arriba llum a la placa
                END IF
 
                ! Seleccionem la distància en el temps corresponent i la normalitzem
                d = distancies((dia - 1) * total_intervals + j) / d_caracteristica
 
                ! Calculem la potència
                potencia = (costheta / (d**2)) * pot_caracteristica
 
                ! Guardem el resultat de la potència
                write(12, '(1X, ES12.5)', advance='no') potencia
                
                if (b==1 .or. b==2) then
                    if (dia == 1) then
                    write(13, '(1X, ES12.5)', advance='no') potencia
                    else if (dia == 170) then
                        write(14, '(1X, ES12.5)', advance='no') potencia
                    end if
                end if
                
                ! Si beta és zero, aturem la iteració (només ho fem una vegada)
                IF(beta==0) THEN
                  exit
                END IF

             END DO
          END DO
          write(12, *)
          if (dia == 1) then
            write(13, *)
          else if (dia==170) then
            write(14, *)
          end if
       end do
    end do
 
    close(12)
    close(13)
    close(14)
    close(15)
 
 end program potenc
 
