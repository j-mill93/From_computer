! A function which calculates the temperature field for a cylindrical aerogel
! crucible furnace with two boron-nitride heating elements (one at the top
! of the rod and one at the bottom). The bottom of the sample is coupled to a
! molybdenum rod which acts as a heat sink.
!
! As the furnace is radially adiabatic through the centre (360 degree radial
! symmetry) a 2-dimensional calculation in cylindrical co-ordinates is used to
! obtain the thermal field through the furnace.
!
! ASCII representation of the furnace:
!
!		     ______________________
!				[______________________] Top of crucible insulation (T_inf = 20C
!  [ o ]|										   |  								= 293K)
!  [ o ]|			 						   	 | Top heater ([ o ]) (T_heater1 = )
!  [ o ]|						 			  	 |
!  [ o ]|						   				 |
!  [~~~]|						   				 |
!  [~~~]|						 				   |
!  [~~~]|						  				 | Al-Cu sample
!  [~~~]|						 				   | Thermal insulation ([~~~]) (T_termal = 20C
!  [~~~]|						  				 | 									= 293K)
!  [~~~]|						   				 |Centre of the sample
!  [ o ]|						  				 |
!  [ o ]|						 				   | Bottom heater ([ o ]) (T_heater2 = )
!  [ o ]|						 				   |
!  [ o ][______________________] Molybdenum-Al-Cu contact point
!  [~~~][						   				 ]
!  [~~~]|						   				 |
!  [~~~]|						   				 |
!  [~~~]|						   				 |
!  [~~~]|						   				 | Molybdenum rod
!  [~~~]|						   				 |
!  [~~~]|						   				 |
!  [~~~]|						  				 |
!  [~~~]|						  				 |
!  [~~~][______________________] (0,0,0) origin (radius, length, angle)
!		Bottom of the molybdenum rod (T_bottom_rod = T_inf)
!
!
! Equations used to obtain the finite difference heat transfer through the cell:
!
!

! begin program
program Aerogel_Furnace_Temp_Sim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialise parameters, 8-bit real numbers and 8-bit integers                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Al-Cu parameters
real(kind = 8) :: radius_rod! Radius of the Al-Cu rod [m]
real(kind = 8) :: length_rod! Length of the Al-Cu rod [m]
real(kind = 8) :: C_AlCu! Heat capacity of the Al-Cu rod [J/kgK]
real(kind = 8) :: dens_AlCu! Density of Al-Cu [kg/m^3]
real(kind = 8) :: k_Al_Cu! Thermal conductivity of Al-Cu [W/Km]

! Furnace parameters
real(kind = 8) :: T_heater1! Temperature of the top heater [K]
real(kind = 8) :: T_heater2! Temperature of the bottom heater [K]
real(kind = 8) :: T_inf! Temperature at infinity, room temperature [K]
real(kind = 8) :: T_0! Temperature of the Al-Cu rod and molybdenum heat
! sink at the time t = 0s [K]
real(kind = 8) :: L_heater1! Length of top heater [m]
real(kind = 8) :: L_heater2! Length of the bottom heater [m]
real(kind = 8) :: L_aerogel! Length of the aerogel crucible (L_AlCu - L_heater1 -
! L_heater2 = L_aerogel) [m]
real(kind = 8) :: R_AlCu_Heater! Thermal resistivity between the boron-nitride
! heater and the Al-Cu rod [K/W]
real(kind = 8) :: T_new! New temperature at time t + Delta t

! Molybdenum parameters
real(kind = 8) :: length_moly! length of the molybdenum heat sink. The radius of
! the rod is equal to the radius of the Al-Cu rod
real(kind = 8) :: C_moly! Heat capacity of the molybdenum rod [J/kgK]
real(kind = 8) :: dens_moly! Density of molybdenum [kg/m^3]
real(kind = 8) :: k_moly! Thermal conductivity of molybdenum [W/Km]
real(kind = 8) :: R_AlCu_Moly! Thermal resistivity between the molybdenum rod and
! the Al-Cu rod [K/W]

! Finite difference parameters
real(kind = 8) :: t_total! Total time taken for the heating/cooling
real(kind = 8) :: Delta_t! The time step [s]
integer(kind = 8) :: no_time_steps! Total number of time steps
real(kind = 8) :: Delta_r! The radius step, used to calculate heat across
! the small volumes during the finite difference
! calculations
real(kind = 8) :: Delta_L! The length step, used to calculate heat across small
! volumes during the finite difference calculations
real(kind = 8) :: K_conduct_N! Thermal conductance, North [W/K]
real(kind = 8) :: K_conduct_E! Thermal conductance, East [W/K]
real(kind = 8) :: K_conduct_S! Thermal conductance, South [W/K]
real(kind = 8) :: K_conduct_W! Thermal conductance, West [W/K]
real(kind = 8) :: Q_N! Thermal flux North
real(kind = 8) :: Q_E! Thermal flux East
real(kind = 8) :: Q_S! Thermal flux South
real(kind = 8) :: Q_W! Thermal flux West

! Real numbers for the changing temperature of the heaters 1 and 2
real(kind = 8) :: T_centre_start ! when the average temperature of the central
                                 ! points of the rod are greater than this value
                                 ! then the heaters will begin to change temperature
real(kind = 8) :: T_centre ! Average temperature in the centre of the rod
real(kind = 8) :: T_heater1_new  ! The new heater1 temperature
real(kind = 8) :: T_heater2_new  ! The new heater2 temperature
real(kind = 8) :: Time_for_T_change ! Time for the change from the initial to final temperatures
                                    ! in the heaters.
real(kind = 8) :: Thermal_gradient ! thermal gradient between the two heaters
                                   ! T_heater1/T_heater2
real(kind = 8) :: T_heater1_Final ! The final temperature of heater 1.
real(kind = 8) :: dTdt ! change in temperature each second
! 1D arrays of radius and length values used in the calculation of the temperature
real(kind = 8), dimension(:), allocatable :: radius_array
real(kind = 8), dimension(:), allocatable :: length_array
real(kind = 8), dimension(:), allocatable :: centre_rod

! Array the data will be occuring over

real(kind = 8), dimension(:, :), allocatable :: Temp_array
real(kind = 8), dimension(:, :), allocatable :: Temp_array_new
real(kind = 8), dimension(:, :), allocatable :: Resistance_array

! Intergers for counting in loops
integer(kind = 8) :: i ! for looping along the rows of the array (length of furnace)
integer(kind = 8) :: j ! for looping along the columns of the array (radius of furnace)
integer(kind = 8) :: t_step ! for looping along the time step of the function
integer(kind = 8) :: n_col_array ! total number of columns in the array
integer(kind = 8) :: n_row_array ! total number of rows in the array
integer(kind = 8) :: heater1_array_length ! length of heater1 in terms of array integers length
integer(kind = 8) :: heater2_array_length ! length of heater2 in terms of array integer length
integer(kind = 8) :: aerogel_array_length ! length of the aerogel crucible in terms of array integer length
integer(kind = 8) :: moly_array_length ! length of the molybdenum heat sink in terms of array integer length
integer(kind = 8) :: rod_array_length ! length of the Cu 10wt% Al rod in terms of array integer length
integer(kind = 8) :: rod_array_radius ! radius of the Cu 10wt% Al rod and molybdenum rod in terms of
                            ! array integer length
integer(kind = 8) :: row ! for looping through variables to print out to a text file
integer(kind = 8) :: col ! for looping through variables to print out to a text file
integer(kind = 8) :: T_to_write
character*32:: filename

! real numbers used to calculate thermal conductivity parameters
real(kind = 8) :: a_poly ! coeffcient a of the polynomial equation
real(kind = 8) :: b_poly ! coeffcient b of the polynomial equation
real(kind = 8) :: c_poly ! coeffcient c of the polynomial equation
real(kind = 8) :: k_ij   ! thermal conductivity at ij
real(kind = 8) :: k_N    ! thermal conductivity at North of ij
real(kind = 8) :: k_E    ! thermal conductivity at East of ij
real(kind = 8) :: k_S    ! thermal conductivity at South of ij
real(kind = 8) :: k_W    ! thermal conductivity at West of ij

!IMPLICIT NONE
INTEGER :: counter, NR
! Assign variables to the parameters !

! Declare constant Pi
real(kind = 8), parameter :: pi= 3.1415927D0

! Non-advancing status counter:
!Shows how to place a non-advancing status counter...




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign the parameters                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NR = 100 ! For the counter

radius_rod = 5D-3 ! Radius of the rod [m]
length_rod = 0.2D0 ! length of the rod in [m]
C_AlCu = 900D0 ! heat capacity of Al-10wt%Cu, this needs to change to be dependent upon temperature.
dens_AlCu = 2550D0! Density of pure Al at T_melt [kg/m^3]
k_Al_Cu = 213D0 ! thermal diffusivity

T_heater1 = 1500D0 ! Temperature of the top heater (heater 1)
T_heater2 = 1500D0 ! Temperature of the bottom heater (heater 2)
T_inf = 0.0D0 ! Room temperature i.e. 20C
T_0 = 273D0 ! Initial temperature of the system (Al-Cu and molybdenum)

T_centre_start = 900D0 ! Temperature for when the cooling curve controlling the
                       ! heaters starts
Time_for_T_change = 1000D0! Time for the change from the initial to final temperatures
                                                           ! in the heaters. [s]
Thermal_gradient = 7D0! thermal gradient between the two heaters
                                                          ! T_heater1/T_heater2
T_heater1_Final = 300D0 ! final heater temperature

dTdt = (T_heater1_Final-T_heater1)/Time_for_T_change

L_heater1 = 0.05D0 ! length of the first heater
L_heater2 = 0.05D0 ! length of the second heater
L_aerogel = 0.05D0 ! length of the aerogel between the two heaters
R_AlCu_Heater = 0D0 ! Resistance between the heaters and the Al-Cu

length_moly = 0.01D0 ! length of the molybdenum heat sink
C_moly = 280D0 ! heat capacity of the molybdenum heat sink
dens_moly = 10220D0 ! density of the molybdenum
k_moly = 138D0 ! thermal conductivity of the molybdenum
R_AlCu_Moly = 0D0 ! Resistance between the molybdenum and the Al-Cu

t_total = 1D3 ! total time the analysis is run over [s]
Delta_t = 1D-1 ! time step
Delta_r = 1D-4 ! The radius step
Delta_L = 1D-3 ! The length step

a_poly = -1.4D-7 ! coefficient a of the polynomial equation
b_poly = 2.12D-3 ! coefficient b of the polynomial equation
c_poly = 1.154   ! coefficient c of the polynomial equation

! Calculate the lengths of the furnace in terms of array points
no_time_steps = (t_total/Delta_t) ! total number of time steps
print* , "no time steps", no_time_steps
n_col_array = (radius_rod/Delta_r) ! Total number of columns in Al-Cu rod
n_row_array = ((length_rod + length_moly)/Delta_L) ! Total number of rows in
! Al-Cu + moly rods

! Calculate the length of the parts of the heater in terms of integer array length
heater1_array_length = (L_heater1/Delta_L)
heater2_array_length = (L_heater2/Delta_L)
aerogel_array_length = (L_aerogel/Delta_L)
rod_array_length = (length_rod/Delta_L)
moly_array_length = (length_moly/Delta_L)
rod_array_radius = (radius_rod/Delta_r)

! Print the values to the console
print*, "number of rows", n_row_array
print*, 'number of columns', n_col_array
print*, 'length of heater1', heater1_array_length
print*, 'length of heater2', heater2_array_length
print*, 'length of aerogel crucible', aerogel_array_length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign the areas of the Temperature array with variables                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(Temp_array(1:(n_row_array + 1), 1:(n_col_array + 1))) ! size of array

  Temp_array = T_inf ! initially set all values to T_inf (Temperature at infinity).
                     ! This is the temperature in the aerogel and top and bottom
                     ! left boundaries
  Temp_array(2:heater1_array_length, 1) = T_heater1 ! Set heater 1 temperature
  Temp_array((heater1_array_length + aerogel_array_length + 1) &
  :(heater1_array_length + aerogel_array_length + 1 + heater2_array_length) &
  , 1) = T_heater2 ! Set heater 2 temperature
  Temp_array(2:n_row_array, 2:n_col_array) = T_0 ! Set initial temperature in
                   ! molybdenum and Al_Cu

  ! Save te first Temp_array as a .csv file
  open(unit=8, file='temp/Temp_array0.csv', ACTION="write", STATUS="replace")
    write(8, *) "X",",", "Y",",","Z",",", "T"
    do row = 1, n_row_array
      do col = 1, n_col_array ! loop along the columns
        T_to_write = Temp_array(row,col)
        write(8, *) real(row*Delta_L),",", real(col*Delta_r), &
        ",",real(0),"," ,real(T_to_write) !, 0, T_to_write) ! Write each row to file
      enddo
  enddo ! end loop along the rows


! Assign area of Temp_array_new used in the loop
allocate(Temp_array_new(1:(n_row_array + 1), 1:(n_col_array + 1))) ! size of array

  Temp_array_new = Temp_array

! Assign resistance to areas which are coupled together
allocate(Resistance_array(1:(n_row_array + 1), 1:(n_col_array + 1))) ! size of array
  Resistance_array(:,:) = 0.0D0

  Resistance_array(2:heater1_array_length, 1) = 0.0D0!R_AlCu_Heater ! Resistance between

  Resistance_array(2:rod_array_length, 1) = 0.0D0!1D100 ! Resistance between
    ! the heater 1 and the Aluminium-Copper
  Resistance_array((heater1_array_length + aerogel_array_length + 1) &
    :heater2_array_length, 1) = 0.0D0 !R_AlCu_Heater ! Resistance between the heater 2
    ! and the aluminium copper

  Resistance_array((rod_array_length + 1), 2:rod_array_radius) = R_AlCu_Moly
    ! Resistance between the Al-Cu rod and the molybdenum heat sink

! Allocate the values of length and radius for the finite difference calculuations
! for every Delta_L and Delta_r at each Delta_t*t_step
allocate(length_array(1:(moly_array_length + rod_array_length)))
  do i = 1,(moly_array_length + rod_array_length)
    length_array(i) = ((length_moly+length_rod)) - i*Delta_L
  enddo

allocate(radius_array(1:rod_array_radius))
  do j = 1, rod_array_radius
    radius_array(j) = (radius_rod) - j*Delta_r
  enddo

print*, "radius_array", size(radius_array)
!print*, "length_array", size(length_rod)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Looping to calculate the new temperature values at each point i,j in the     !
! array using finite difference method.                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do t_step = 1, no_time_steps   ! start time loop

    ! Step counter. Outputs percentage to completion to the console
    write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
    & " Percent Complete: ", (real(t_step)/real(no_time_steps))*100.0, "%"
    ! initialise a new temperature array the same as the previous loops array
  Temp_array_new = Temp_array

  centre_rod = Temp_array_new(2:rod_array_length,size(radius_array))

  T_centre = sum(centre_rod)/(max(1,size(centre_rod)))

  !print*,T_centre,T_heater1,Temp_array_new(2,2),Temp_array_new(2,size(radius_array)-2),Temp_array_new(2,size(radius_array))
  ! If the temperature in the centre of the rod is greater than the value of T_centre_start
  ! then a new temperature for the heaters is to be calculated dependent upon the timestep
  ! and initial and final temperatures of the heaters.

    if (T_centre .gt. T_centre_start &
        .and. T_heater1 .ge. T_heater1_Final) then ! start changing temperatures of the heaters

      T_heater1_new = T_heater1 + dTdt
      T_heater2_new = T_heater1_new - Thermal_gradient

      T_heater1 = T_heater1_new
      T_heater2 = T_heater2_new

      Temp_array_new(2:heater1_array_length, 1) = T_heater1 ! Set heater 1 temperature
      Temp_array_new((heater1_array_length + aerogel_array_length + 1) &
      :(heater1_array_length + aerogel_array_length + 1 + heater2_array_length)&
      , 1) = T_heater2 ! Set heater 2 temperature
                       ! molybdenum and Al_Cu
    endif ! stop the changing temperatures of the heaters

      do i = 2, (size(length_array) - 1)  ! loop over length of the array
        do j = 2, (size(radius_array) - 1) ! loop over radius (width) of the array


    !Assign variables needed for the loop to calculate the new Temperature, T_new
        radius_j = radius_array(j) ! the radius at the point i,j
        radius_East = radius_array(j + 1) ! radius to point to east of i, j
        radius_West = radius_array(j - 1) ! radius to point west of i, j
        Temp_ij = Temp_array(i,j) ! Temperature of point i, j
        Temp_N = Temp_array(i + 1, j) ! Temperature of point to north of i, j
        Temp_E = Temp_array(i,j + 1) ! Temperature of point to East of i, j
        Temp_S = Temp_array(i - 1, j) ! Temperature of point to South of i, j
        Temp_W = Temp_array(i, j - 1) ! Temperature of point to West of i, j


      ! Calculate thermal conductivities

      ! If the the point lies in the Aluminium-Copper region its thermal conductivity
      ! is given by:
      ! Calculate the thermal conductivity of all of the relevent points, i.e.
      ! k_ij, k_N, k_E, k_S,k_W. The data is calculated from a polynomial fitted
      ! curve of thermal conductivity data for Al-10wt%Cu from 200K to 800K
      ! data from NIST:
      ! http://www.nist.gov/data/PDFfiles/jpcrd123.pdf
      if (i .gt. 1 .and. i .le.rod_array_length) then
        k_ij = a_poly*Temp_ij + b_poly*Temp_ij + c_poly
        k_N = a_poly*Temp_N + b_poly*Temp_N + c_poly
        k_E = a_poly*Temp_E + b_poly*Temp_E + c_poly
        k_S = a_poly*Temp_S + b_poly*Temp_S + c_poly
        k_W = a_poly*Temp_W + b_poly*Temp_W + c_poly
      endif

      ! If the the point lies in the molybdenum region its thermal conductivity
      ! is given by:
      if (i .gt. rod_array_length .and. i .le.rod_array_length &
            + moly_array_length) then
        k_ij = a_poly*Temp_ij + b_poly*Temp_ij + c_poly
        k_N = a_poly*Temp_N + b_poly*Temp_N + c_poly
        k_E = a_poly*Temp_E + b_poly*Temp_E + c_poly
        k_S = a_poly*Temp_S + b_poly*Temp_S + c_poly
        k_W = a_poly*Temp_W + b_poly*Temp_W + c_poly
      endif



      ! Density calculuations for Aluminium-Copper and the molybdenum

      ! Calculate the density of the Aluminium-Copper. The data for the calculuations
      ! was taken from the Thermocalc database TCFe7. Three curves have been fitted
      ! to the data which correspond to the temperature ranges 150K - 715K, 715K
      ! - 840K, and 840K+. The graphs for the corresponding regions are included
      ! in the file which contains this program.

      ! Calculate the density of the first temperature region: 150K - 715K
      if (Temp_ij .lt. 715) then
          dens_AlCu = (-0.0000726)*(Temp_ij**2) + (-0.149)*Temp_ij + 2.93d3
        ! Calculate the density of the second temperature region 715K - 840K
      elseif (Temp_ij .ge. 715 .and. Temp_ij .lt. 840) then
          dens_AlCu = (7.77d6) + (-5.08d4)*Temp_ij + 133*(Temp_ij**2) + (-0.174) &
            *(Temp_ij**3) + 0.000114*(Temp_ij**4) + (-2.98d-8)*(Temp_ij**5)
      elseif (Temp_ij .ge. 840) then
          dens_AlCu = (-0.230)*Temp_ij + 2.76d3
      endif

      ! Calculate the heat capacity of the Aluminium-Copper. The data for the
      ! calculations was taken from the thermocalc database TCFe7. Three curves
      ! have been fitted to the data which correspond to the temperature ranges
      ! 80K to 820K, 820K to 905K, and 950K+.

      ! Calculate the heat capacity of the first region: 80K to 820K
      if (Temp_ij .lt. 820) then
        C_AlCu = (-56.3) + (1.2)*Temp_ij + (-0.00730)*(Temp_ij**2) &
          + (0.0000228)*(Temp_ij**3) + (-3.78d-8)*(Temp_ij**4) &
          + (3.15d-11)*(Temp_ij**5) + (-1.04d-14)*(Temp_ij**6)
      ! Calculate the heat capacity of the second region: 820K to 905K
      elseif (Temp_ij .ge. 820 .and. Temp_ij .lt. 905) then
      !  C_AlCu = (1.81d7) + ((-8.52d4)*Temp_ij) + ((150)*(Temp_ij**2)) &
      !    + ((-0.118)*(Temp_ij**3)) + ((0.0000345)*(Temp_ij**4))
C_AlCu = 300
          !print*(C_AlCu)
      ! Calculte the heat capacity of the third region: 905K+
      elseif (Temp_ij .ge. 905) then
        C_AlCu = 31.2
      endif
       ! Calculate the conductivity across the edges of the finite volumes

       ! North and South vectors
       K_conduct_N = (2*pi*radius_array(j)*Delta_r) &
       /(((0.5*Delta_L)/k_N) + ((0.5*Delta_L)/(k_ij)) &
       + (Resistance_array(i, j)))

       if (Temp_N .lt. T_inf + 1) then
         K_conduct_N = 0D0
       endif

       K_conduct_S = (2*pi*radius_array(j)*Delta_r)/((((0.5*Delta_L) &
         /k_S) + ((0.5*Delta_L)/k_ij)) + (Resistance_array(i, j)))

       if (Temp_S .lt. T_inf + 1) then
           K_conduct_S = 0D0
       endif
       ! East and West vectors

       K_conduct_E = (Delta_L) &
       /(((1/(2*pi*k_E))*log((radius_array(j)+(Delta_r/2)&
       /radius_array(j + 1)))) + ((1/(2*pi*k_ij)) &
       *log((radius_array(j))/(radius_array(j) + (Delta_r/2)))) &
      + (Resistance_array(i,j)/(2*pi*(radius_array(j) + (Delta_r/2)))))

      if (Temp_E .lt. T_inf + 1) then
        K_conduct_E = 0D0
      endif

       K_conduct_W = (Delta_L) &
       /(((1/(2*pi*k_W))*log((radius_array(j)-(Delta_r/2)) &
       /radius_array(j - 1))) + ((1/(2*pi*k_ij))*log((radius_array(j)) &
      /radius_array(j) - (Delta_r/2))) + ((Resistance_array(i,j) &
       /(2*pi*(radius_array(j) - (Delta_r/2))))))

       if (Temp_W .lt. T_inf + 1) then
         K_conduct_W = 0D0
       endif
      ! Calculate the fluxes
       Q_N = K_conduct_N*(Temp_N - Temp_ij)
       Q_E = K_conduct_E*(Temp_E - Temp_ij)
       Q_S = K_conduct_S*(Temp_S - Temp_ij)
       Q_W = K_conduct_W*(Temp_W - Temp_ij)
      ! Calculate the new Temperature at the point i,j and time t_step

       if (Q_S .gt. Q_N .and. Q_E .gt. Q_W) then
          T_new = (((-Q_W + Q_E + Q_S - Q_N)*Delta_t) &
             / (2*pi*(C_AlCu*dens_AlCu)*radius_array(j)*Delta_r)) + Temp_ij

          Temp_array_new(i,j) = T_new
        elseif (Q_S .gt. Q_N .and. Q_E .lt. Q_W) then
             T_new = (((Q_W - Q_E + Q_S - Q_N)*Delta_t) &
                / (2*pi*(C_AlCu*dens_AlCu)*radius_array(j)*Delta_r)) + Temp_ij
        elseif (Q_S .gt. Q_N .and. Q_E .eq. Q_W) then
             T_new = (((Q_S - Q_N)*Delta_t) &
                / (2*pi*(C_AlCu*dens_AlCu)*radius_array(j)*Delta_r)) + Temp_ij
             Temp_array_new(i,j) = T_new
       elseif (Q_N .lt. Q_S .and. Q_E .gt. Q_W) then
          T_new = (((- Q_W + Q_E - Q_S + Q_N)*Delta_t) &
             / (2*pi*(C_AlCu*dens_AlCu)*radius_array(j)*Delta_r)) + Temp_ij
          Temp_array_new(i,j) = T_new
        elseif (Q_N .lt. Q_S .and. Q_E .eq. Q_W) then
           T_new = (((- Q_S + Q_N)*Delta_t) &
              / (2*pi*(C_AlCu*dens_AlCu)*radius_array(j)*Delta_r)) + Temp_ij
           Temp_array_new(i,j) = T_new
        elseif (Q_N .lt. Q_S .and. Q_E .lt. Q_W) then
           T_new = (((- Q_W + Q_E - Q_S + Q_N)*Delta_t) &
              / (2*pi*(C_AlCu*dens_AlCu)*radius_array(j)*Delta_r)) + Temp_ij
           Temp_array_new(i,j) = T_new
       else


         T_new = (((-Q_W + Q_E + Q_S - Q_N)*Delta_t) &
            / (2*pi*(C_AlCu*dens_AlCu)*radius_array(j)*Delta_r)) + Temp_ij

         Temp_array_new(i,j) = T_new

        endif

    end do   ! end loop over radius

    Temp_array_new(i,size(radius_array)) = Temp_array(i,size(radius_array)-1)-10

  end do   ! end loop over length

  ! Assign new Temperature array to the old one.
  Temp_array = Temp_array_new

  ! Save the new Temperature array as a .txt file from 2 -> the last one.
  ! The first file was saved above earlier.
  if (MOD(t_step,10000) .eq. 0) then
  write (filename, '( "temp/Temp_array", i0,".csv" )' ) t_step
  open(unit=t_step,file=filename) ! open file
      write(t_step, *) "X",",", "Y",",","Z",",", "T"
      do row = 1, n_row_array
      do col = 1, n_col_array ! loop along the columns
        T_to_write = Temp_array(row,col)
      write(t_step, *) real(row*Delta_L),",", real(col*Delta_r), &
      ",",real(0),"," ,real(T_to_write)!, 0, T_to_write) ! Write each row to file
      enddo
      enddo ! end loop along the rows
  close (t_step) ! close the file Temp_array"t_step".txt
  endif
!enddo ! end time loop

enddo ! End the status counter
end program Aerogel_Furnace_Temp_Sim
