!-----------------------------------------------------------------------
!Module: physics
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!!
!! Description: contains calculations of pressure, energy density, etc.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_p_e_Z
!! populate_arrays
!! TOV
!! 
!-----------------------------------------------------------------------
module physics
use types
use math, only : linear_interpolation

implicit none

private
public :: populate_arrays, calculate_p_e_Z, TOV, read_p_e

contains

!-----------------------------------------------------------------------
!! Subroutine: populate_arrays
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!!
!! Description: This subroutine populates parameter arrays.
!!----------------------------------------------------------------------
!
! Input:
!
!   flavor_max      - Number of flavors/leptons to consider
!!----------------------------------------------------------------------
! Output:
!
!   m_f             - Array containing quark masses for each flavor
!   k_f             - Array containing Fermi momentum per quark flavor
!   g               - 1D integer array degeneracy factor
!   q               - charge of each species
!
!-----------------------------------------------------------------------
subroutine populate_arrays(flavor_max, m_f, k_f, g, q)
    implicit none
    integer, intent(in) :: flavor_max
    real(dp), intent(out), allocatable :: m_f(:), k_f(:), q(:)
    integer, intent(out), allocatable :: g(:)
    integer :: AllocateStatus

    ! allocate space
    allocate(m_f(flavor_max), stat=AllocateStatus)
    if(AllocateStatus /= 0) stop "*** Not enough memory ***"
    allocate(k_f(flavor_max), stat=AllocateStatus)
    if(AllocateStatus /= 0) stop "*** Not enough memory ***"
    allocate(g(flavor_max), stat=AllocateStatus)
    if(AllocateStatus /= 0) stop "*** Not enough memory ***"
    allocate(q(flavor_max), stat=AllocateStatus)
    if(AllocateStatus /= 0) stop "*** Not enough memory ***"

    ! masses
    m_f(1) = 2.55_dp                            ! mass of up quark in MeV              
    m_f(2) = 4.865_dp                           ! mass of down quark in MeV            
    m_f(3) = 101.0_dp                           ! mass of strange quark in MeV
    m_f(4) = 0.510998_dp                        ! mass of electron in MeV
    if (flavor_max > 4) m_f(5) = 105.6583755_dp ! mass of muon in MeV
    if (flavor_max > 5) m_f(6) = 1270.0_dp      ! mass of CHARM QUARK

    ! Fermi momenta
    k_f(1) =  1.0_dp * conversion_factor        ! k_up in MeV
    k_f(4) = k_f(1) * ( m_f(4) / m_f(1) )       ! k_electron as a fraction of k_up --> k_up * ratio of masses
    ! other Fermi momenta not defined

    ! degeneracy factors
    g(1) = 6 ! 6 for quarks
    g(2) = 6
    g(3) = 6
    g(4) = 2 ! 2 for leptons
    if (flavor_max > 4) g(5) = 2 ! muon
    if (flavor_max > 5) g(6) = 6 ! charm quark

    ! charges
    q(1) = 2.0_dp/3.0_dp    ! +2/3 for up quarks
    q(2) = -1.0_dp/3.0_dp   ! -1/3 for down quarks
    q(3) = -1.0_dp/3.0_dp   ! -1/3 for strange quarks
    q(4) = -1.0_dp  ! -1 for electrons
    if (flavor_max > 4) q(5) = -1.0_dp  ! -1 for muons
    if (flavor_max > 5) q(6) = 2.0_dp/3.0_dp ! +2/3 for charm quarks

end subroutine populate_arrays

!-----------------------------------------------------------------------
!! Subroutine: calculate_p_e_Z
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!
!  This subroutine takes 6 values:
!   
!   k                   - Constant multiplying k_f to update its value for each iteration
!   k_f                 - Array containing Fermi momentum per quark flavor
!   m_f                 - Array containing quark masses for each flavor
!   g                   - 1D integer array degeneracy factor
!   q                   - charge of each species
!   flavors             - The number of quark flavors/types of particles to be considered 
!   mu_electron         - Electron chemical potential value
!               
!
! And calculates pressure and energy density using these values:
!
!   E                   - Energy density
!   P                   - Pressure
!
! This subroutine also calculates the number densities and chemical potentials of
! each species, as well as the total number density.
!
!   mu_f                - Array containing individual chemical potentials
!   n                   - Array containing individual number densities
!   n_tot               - Total number density
!   chem_EQ             - Array containing chemical equilibrium conditions
!   charge_neutrality   - Value of total charge
!
!
!-----------------------------------------------------------------------
subroutine calculate_p_e_z(k, k_f, m_f, g, q, flavors, mu_electron, epsilon, P, mu_f, n, n_tot, chem_eq, charge_neutrality)
    implicit none
    real(dp), intent(in) :: k, m_f(:), k_f(:), mu_electron, q(:)
    integer, intent(in) :: flavors, g(:)
    real(dp), intent(out) :: n_tot, chem_eq(1:2), charge_neutrality, epsilon, P
    real(dp), allocatable :: z(:), mu_f(:), n(:)
    real(dp) :: sum_p, sum_e, p_constant, p_term_1, p_term_2, e_constant, e_term_1, e_term_2, b, epsilon_charm
    real(dp) :: n_e_diff, n_mu_diff
    integer :: i, AllocateStatus, flavors_used

    flavors_used = flavors

    ! allocate space
    allocate(z(flavors_used), stat=AllocateStatus)
    if(AllocateStatus /= 0) stop "*** Not enough memory ***"


    ! populate chemical potential array mu_f

    mu_f(1) = sqrt(m_f(1)**2.0_dp + (k*k_f(1))**2.0_dp)       ! Set mu_up based off of k from iteration
    !mu_f(4) = mu_electron                                    ! mu_e (OLD)

    mu_f(4) = sqrt(m_f(4)**2.0_dp + (k*k_f(4))**2.0_dp)       ! NEW

    mu_f(2) = mu_f(1) + mu_f(4)                               ! mu_down = mu_up + mu_e
    mu_f(3) = mu_f(2)                                         ! mu_strange = mu_down
    if (flavors_used > 4) mu_f(5) = mu_f(4)                   ! mu_mu = mu_e
    if (flavors_used > 5) mu_f(6) = mu_f(1)                   ! mu_charm = mu_up


    ! Initializations to zero
    sum_p = 0.0 ; sum_e = 0.0 ; n_tot = 0.0 ; n = 0.0

    ! sum flavor contributions for pressure, energy density, number density, and chemical potential
    do i = 1, flavors_used

        z(i) = m_f(i) / mu_f(i)

        if (z(i) <= 1.0_dp) then    ! only contribute to the sums if the chem pot is beyond a certain threshold
                                    ! this allows charm quarks, etc. to only come into play beyond a certain energy density
            ! pressure sum
            p_constant = (g(i) * mu_f(i)**4.0_dp) / (24.0_dp * (pi**2.0_dp))
            p_term_1 = (sqrt(abs(1.0_dp - (z(i)**2.0_dp)))) * (1.0_dp - (2.5_dp * (z(i)**2.0_dp)))
            p_term_2 = 1.5_dp * (z(i)**4.0_dp) * log((1.0_dp + sqrt(abs(1.0_dp - (z(i)**2.0_dp)))) / z(i))

            sum_p = sum_p + p_constant * (p_term_1 + p_term_2)

            ! energy density sum
            e_constant = (g(i) * (mu_f(i)**4.0_dp)) / (8.0_dp * (pi**2.0_dp))
            e_term_1 = (sqrt(abs(1.0_dp - (z(i)**2.0_dp)))) * (1.0_dp - (0.5_dp * (z(i)**2.0_dp)))
            e_term_2 = 0.5_dp * (z(i)**4.0_dp) * log((1.0_dp + sqrt(abs(1.0_dp - (z(i)**2.0_dp)))) / z(i))

            sum_e = sum_e + e_constant * (e_term_1 - e_term_2)

            ! number density
            n(i) = ((g(i) * mu_f(i))**3.0_dp) * (abs(1.0_dp - (z(i)**2.0_dp)))**(3.0_dp/2.0_dp) / (6.0_dp * (pi**2.0_dp))

            n_tot = n_tot + n(i)
        
        endif

    end do

    P = -bag * conversion_factor**3.0_dp + sum_p

    ! solve for energy density
    epsilon = bag * conversion_factor**3.0_dp + sum_e

    ! charge neutrality (should be zero)
    charge_neutrality = 0.0_dp
    do i = 1, flavors_used
        charge_neutrality = charge_neutrality + q(i) * ( n(i) / conversion_factor**3.0_dp )
    end do


!!________________________________________________________________________________________
    !! NOT YET FINISHED
    ! ! neutralize charge by adding electrons:

    if (z(5) >= 1.0_dp) then        ! if muons are not yet present




        n_e_diff = ( charge_neutrality / abs( q(4) ) ) * conversion_factor**3.0_dp
        ! Define the number of electrons to relate to how off the charge is


        

        if (charge_neutrality > 0) n(4) = n(4) + n_e_diff ! if charge neutrality is greater than zero
        if (charge_neutrality < 0) n(4) = n(4) - n_e_diff ! if charge neutrality is greater than zero



        ! set chemical potential of electrons after this charge neutralization

        ! n(i) = ((g(i) * mu_f(i))**3.0_dp) * (abs(1.0_dp - (z(i)**2.0_dp)))**(3.0_dp/2.0_dp) / (6.0_dp * (pi**2.0_dp))

        mu_f(4) = (  n(4) * (6.0_dp * (pi**2.0_dp)) / ( g(4) *  (abs(1.0_dp - (z(4)**2.0_dp)))**(3.0_dp/2.0_dp) )  )&
         ** (1.0_dp / 3.0_dp)


    end if

                                    


    ! if muons are present:


    if (z(5) < 1.0_dp) then    



        n_mu_diff = ( charge_neutrality / abs( q(5) ) ) * conversion_factor**3.0_dp 
        ! Define the number of electrons to relate to how off the charge is

        

        if (charge_neutrality > 0) n(5) = n(5) + n_mu_diff ! if charge neutrality is greater than zero
        if (charge_neutrality < 0) n(5) = n(5) - n_mu_diff ! if charge neutrality is greater than zero


        ! set chemical potential of electrons after this charge neutralization

        ! n(i) = ((g(i) * mu_f(i))**3.0_dp) * (abs(1.0_dp - (z(i)**2.0_dp)))**(3.0_dp/2.0_dp) / (6.0_dp * (pi**2.0_dp))

        mu_f(5) = (  n(5) * (6.0_dp * (pi**2.0_dp)) / ( g(5) *  (abs(1.0_dp - (z(5)**2.0_dp)))**(3.0_dp/2.0_dp) )  )&
         ** (1.0_dp / 3.0_dp)


    end if

!!________________________________________________________________________________________


! Set charge neutrality again

    ! charge neutrality (should be zero)

    charge_neutrality = 0.0_dp
    do i = 1, flavors_used
        charge_neutrality = charge_neutrality + q(i) * ( n(i) / conversion_factor**3.0_dp )
    end do






    ! chemical equilibrium (both should be zero)
    chem_eq(1) = mu_f(2) - mu_f(1) - mu_f(4)    ! mu_down = mu_up + mu_electron
    chem_eq(2) = mu_f(2) - mu_f(3)              ! mu_down = mu_strange






   

end subroutine calculate_p_e_z

!-----------------------------------------------------------------------
!! Subroutine: TOV
!-----------------------------------------------------------------------
!! By: David Wilkins, Luke Glass, and Fridolin Weber
!!
!! Description: 
!! This subroutine solves the equation of hydrostatic
!! equilibrium, outputting mass and radius values into output files.
!! This subroutine takes two arrays corresponding to y and x values from
!! given file, referenced in MAIN.f90
!!----------------------------------------------------------------------
!
! Input:
!
! p_work            1D array containing pressure values from file (MeV/fm^3)
!
! e_work            1D array containing energy density values from file (MeV/fm^3)
!
!-----------------------------------------------------------------------
subroutine TOV(p_work, e_work)
    implicit none
    real(dp) :: e_c, p_c, e, p, f, pdr, deltaec, RNS_km, MNS_Msun, z
    real(dp) :: cf, mass, r, dr
    integer :: choice
    real(dp), intent(in) :: p_work(:), e_work(:)

    open (unit = 11, file ='MvsR_new.dat', status ='unknown')
    open (unit = 22, file ='zvsM_new.dat', status ='unknown')

    ! write file header

    write(11,*) 'MvsR data'
    write(11,'(3a15)') 'Radius (km)', 'Mass (M_sun)', 'Central Density (MeV/fm^3)'

    write(*,*) 'Newtonian (1) or General Relativistic (2) calculation?'
    read(*,*) choice

    if (choice == 1) then
        write(*,*) 'Newtonian stellar model'
        else if (choice == 2) then
        write(*,*) 'General Relativistic stellar model'
        else
        write(*,*) 'Input error -> computation terminated!'
        stop
    end if

    cf = msun_km * c18 / msun_mev               ! Compute conversion factor
    e_c = 4.2_dp * bag                          ! Initialize central density of 1st stellar model
                                                ! e_c in MeV/fm^3
    ! deltaec = bag / 50.0_dp                     ! Initialize increase in central density of next
                                                ! stellar model 




    do while (e_c <= (250000.0_dp * bag))           ! Compute entire sequence of stellar models 
        




        ! * Original p_c
        ! ** REMOVE THIS FOR NEW MODEL THAT INCLUDES INTERPOLTION
        ! p_c = (e_c - 4.0_dp * bag) / 3.0_dp     ! Initialize central pressure


        ! **** call interpolation subroutine to find pressure value for density
        call linear_interpolation(e_work, p_work, e_c, p_c)




        mass = 0.0_dp                           ! Initialize star's mass
        r = 1.0e16                              ! Initialize radial distance in fermi 
        dr = 1.0e16                             ! Initialize step size in fermi
        e = e_c                                 ! Initialize e
        p = p_c                                 ! Initialize p
                                    
        pdr = 1.0_dp

        do while (pdr > 0.0_dp)                 ! Integrate TOV equation
            if (choice == 1) then               ! Choice 1 = Newtonian
                f = e * mass * cf / (r * r)
                else if (choice == 2) then      ! Choice 2 = Einsteinian
                f = (e + p) * (4.0_dp * pi * r**3.0_dp * p + mass) * cf / (r * r * (1.0_dp - 2.0_dp * mass * cf / r))
            end if
            mass = mass + 4.0_dp * pi * r * r * e * dr    
                                                ! mass comes from integrating energy density over the surface 
                                                ! (A=4pi*r^2) with respect to r. Energy density times volume gives 
                                                ! total energy (mass)

            pdr  = p - f * dr               
            p = pdr                               ! update p for next iteration
        ! e = 3.0_dp * pdr + 4.0_dp * bag         ! update e for next iteration 
                                                ! ** REMOVE THIS FOR NEW MODEL
                                                ! THAT INCLUDES INTERPOLTION


        ! **** call interpolation subroutine to find density value for given pressure
        call linear_interpolation(p_work, e_work, p, e)

            r = r + dr                      
                                                ! pdr in units of F/m
                                                ! Energy: kg*m^2/s^2, force: kg*m/s^2, force*distance
        end do                                  ! Energy density in F/m = kg/s^2, E/m^2? Not per unit volume?


        ! First, the force of gravity is computed (depending on choice 1 or 2, Newtonian or Einsteinian 
        ! gravity is considered).
        ! Second, the mass is computed from the energy density and the volume inside the shell

!         write(*,99) r / c18, e_c, mass / msun_mev
! 99      format(' R=',F7.2,' km',4x,'e_c=',F7.2,' MeV/fm^3', 4x,'M/M_sol=',F7.4) 
        
        RNS_km   = r / c18
        MNS_Msun = mass/msun_mev
        write(11,20) RNS_km, MNS_Msun, e_c           ! Output r in km, mass of NS in M_sun
20      format(2x, f8.4, 4x, f9.4, 4x, f12.2)
        z = 1.0_dp / sqrt(1.0_dp - 2.0_dp * MNS_Msun * msun_km/ RNS_km) - 1.0_dp
        write(22,20) mass / msun_mev, z         ! Output mass in M_sun, grav. redshift

        deltaec = 0.001_dp * e_c
        e_c = e_c + deltaec                     ! Central density of next stellar model

    end do                                      ! Compute next stellar model

    close(unit = 22)
    close(unit = 11)

end subroutine TOV
!!---------------------------------------------------------------------------------------


!-----------------------------------------------------------------------
! Subroutine: read_p_e
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!!
! Description: 
!
! This subroutine reads pressure from the first column and energy
! density from the second column of a data file into respective arrays.
!
!!----------------------------------------------------------------------
!
! Input:
!
!   filename       character filename corresponding to the input file
!
!!----------------------------------------------------------------------
! Output:
!
!   p              1D array containing value of pressure (y value)
!
!   e              1D array containing value of energy density (x value)
!
!
!-----------------------------------------------------------------------
subroutine read_p_e(p, e, filename)
    implicit none
    character(len=50), intent(in) :: filename
    real(dp), allocatable, intent(out) :: p(:), e(:)
    integer :: file_size, i, stat
    real(dp) :: dummy_real
    character(len=50) :: dummy_character

    ! open file to read p and e from
    open (unit=100, file=filename, status='old')

    ! skip first two header lines
    read(100,*) dummy_character
    read(100,*) dummy_character

    ! Determine size of file
    file_size = 0   ! Initialize file_size index

40  read(100,*, end=25) dummy_real ! read dummy variable on each line
    file_size = file_size + 1 ! update file size

        goto 40 ! Reiterate

25  continue ! exit implicit loop

    ! rewind file
    rewind(100)

    ! skip first two header lines again
    read(100,*) dummy_character
    read(100,*) dummy_character

    ! allocate arrays
    allocate(p(1:file_size))
    allocate(e(1:file_size))

    ! read p and e from file
    do i = 1, file_size

        read(100,*) e(i), p(i)

    end do

    print*, "file length = ", file_size

    ! close file
    close(100)

end subroutine read_p_e

!--------------------------------------------------------------------------------




end module physics