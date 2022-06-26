!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!!
!! Description: This module writes the necessary values to output files.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! write_bag_model
!! write_new_eq_of_state
!! 
!-----------------------------------------------------------------------
module read_write
use types
use physics, only : calculate_p_e_Z, read_p_e

implicit none

private
public :: write_bag_model, write_new_eq_of_state, write_new_eq_of_state_phase_change, write_phase_change

contains

!-----------------------------------------------------------------------
!Subroutine: write_bag_model
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!!
!! Description: This subroutine writes the bag model not considering
!! quarks with mass.
!! 
!-----------------------------------------------------------------------
subroutine write_bag_model(filename)
    implicit none
    real(dp) :: epsilon_c, P_c
    character(len=50), intent(out) :: filename
    real(dp), allocatable :: p(:), e(:)
    
    filename = 'bag_model.dat'

    ! open output file
    open (unit=20, file=filename, status='unknown')

    ! write file header
    write(20,*) 'Bag Model Equation of State'
    write(20,'(2a36)') 'central energy density (MeV/fm^3)', 'central pressure (MeV/fm^3)'

    epsilon_c = 4.0_dp * bag           ! Initialize central energy density (in MeV/fm^3)

    do while (epsilon_c <= (25.0_dp*bag))           
        P_c = (epsilon_c - 4.0_dp*bag) / 3.0_dp          ! determine central pressure
        write(20,'(2g36.16)') epsilon_c, P_c
        epsilon_c = epsilon_c + 0.05_dp * bag            ! update central energy density
    end do

    close(20)

end subroutine write_bag_model

!-----------------------------------------------------------------------
!Subroutine: write_new_eq_of_state
!-----------------------------------------------------------------------
!By: David Wilkins and Luke Glass
!
!  This subroutine takes 6 values:
!
!   flavors             - The number of quark flavors/types of particles to be considered 
!   k_f                 - Array containing Fermi momentum per quark flavor
!   q                   - charge of each species
!   m_f                 - Array containing quark masses for each flavor
!   g                   - 1D integer array degeneracy factor
!
!              
! And writes pressure and energy density using the above values in the 
! caluclate_p_e_Z subroutine :
!
!   E                   - Energy density
!   P                   - Pressure
!
! This subroutine also writes the number densities and chemical potentials of
! each species, as well as the total number density.
!
!   mu_f                - Array containing individual chemical potentials
!   n                   - Array containing individual number densities
!   n_tot               - Total number density
!   chem_EQ             - Array containing chemical equilibrium conditions
!   charge_neutrality   - Value of total charge
!
!-----------------------------------------------------------------------
subroutine write_new_eq_of_state(flavor_max, k_f, q, m_f, g)
    implicit none
    real(dp), intent(in) :: m_f(:), k_f(:), q(:)
    integer, intent(in) :: flavor_max, g(:)
    real(dp) :: epsilon_c, P_c, n_tot, chem_eq(1:2), charge_neutrality, k, k_step, k_max
    real(dp) :: mu_electron!, mu_electron_step, mu_electron_max, mu_up, mu_up_step, mu_up_max
    real(dp), allocatable :: mu_f(:), n(:)
    integer :: i, AllocateStatus

    ! allocate space
    allocate(mu_f(flavor_max), stat=AllocateStatus)
    if(AllocateStatus /= 0) stop "*** Not enough memory ***"
    allocate(n(flavor_max), stat=AllocateStatus)
    if(AllocateStatus /= 0) stop "*** Not enough memory ***"
    
    ! open output files

    open (unit=30, file='new_eq_of_state.dat', status='unknown')
    open (unit=40, file='chem_pot.dat', status='unknown')
    open (unit=41, file='num_dens.dat', status='unknown')
    open (unit=50, file='charge_neutrality.dat', status='unknown')
    open (unit=60, file='chem_eq.dat',status='unknown')
    open (unit=110, file='muon_chem_pot_vs_num_dens.dat',status='unknown')
    open (unit=120, file='charm_chem_pot_vs_num_dens.dat',status='unknown')

    ! write file headers

    write(30,*) 'New Equation of State Model'
    write(30,'(2a36)') 'central energy density (MeV/fm^3)', 'central pressure (MeV/fm^3)'

    write(40,*) 'Chemical Potentials (MeV) vs Number Densities (1/fm^3)'
    write(40, '(6a28)') 'Mu_Up', 'Mu_Down', 'Mu_Strange', 'Mu_Electron', 'Mu_Muon', 'Mu_Charm'

    write(41, *) 'Number Densities (1/fm^3)'
    write(41, '(7a28)') 'N_Up', 'N_Down', 'N_Strange', 'N_Electron', 'N_Muon', 'N_Charm', 'N_Tot'

    write(50,*) "Charge Neutrality"
    write(50,'(2a28)') 'iteration', 'charge'

    write(60,*) "Chemical equilibrium"
    write(60,'(3a28)') 'iteration', 'mu_d - mu_u + mu_e = 0', 'mu_d - mu_s = 0'

    write(110,*) 'Muon Chemical Potential (MeV) vs Number Densities (1/fm^3)'
    write(110, '(4a28)') 'Mu_Muon', 'muon num dens', 'total num dens', 'energy density'

    write(120,*) 'Charm Quark Chemical Potential (MeV) vs Number Densities (1/fm^3)'
    write(120, '(4a28)') 'Mu_Charm', 'charm num dens', 'total num dens', 'energy density'



    k = 0.41_dp                                  ! loop iterates over k
    k_step = 0.0001_dp
    k_max = 17.5_dp                            

    i = 1

    do while (k <= k_max)

        ! solve for central pressure and energy density
        call calculate_p_e_Z(k, k_f, m_f, g, q, flavor_max, mu_electron, epsilon_c, P_c,  &
                            & mu_f, n, n_tot, chem_EQ, charge_neutrality)

        ! write values to output files

        ! eq of state ; epsilon_c and P_c in MeV/fm^3
        write(30,'(2g28.16)') epsilon_c/conversion_factor**3.0_dp, P_c/conversion_factor**3.0_dp

        ! individual chem pots (mu in MeV)
        write(40,'(6g28.16)') mu_f(1), mu_f(2), mu_f(3), mu_f(4), mu_f(5), mu_f(6)

        ! individual and total num densities (n in 1/fm^3)
        write(41,'(7g28.16)' ) n(1)/conversion_factor**3.0_dp, &
                            & n(2)/conversion_factor**3.0_dp, n(3)/conversion_factor**3.0_dp, &
                            & n(4)/conversion_factor**3.0_dp, n(5)/conversion_factor**3.0_dp, &
                              n(6)/conversion_factor**3.0_dp, n_tot/conversion_factor**3.0_dp
    
        ! charge neutrality 
        write(50, '(2g28.16)') real(i), charge_neutrality

        ! chem eq
        write(60, '(3g28.16)') real(i), chem_eq(1), chem_eq(2)

        ! muon chemical potential vs number densities
        if (flavor_max > 4) write(110, '(4g28.16)') mu_f(5), n(5)/conversion_factor**3.0_dp, &
                        & n_tot/conversion_factor**3.0_dp, epsilon_c/conversion_factor**3.0_dp

        ! charm chemical potential vs number densities
        if (flavor_max > 5) write(120, '(4g28.16)') mu_f(6), n(6)/conversion_factor**3.0_dp, &
                        & n_tot/conversion_factor**3.0_dp, epsilon_c/conversion_factor**3.0_dp

        ! update k
        k = k + k_step

        ! update iteration count
        i = i + 1

    end do


!   Close output files
    close(30)
    close(40)
    close(41)
    close(50)
    close(60)

end subroutine write_new_eq_of_state


!-----------------------------------------------------------------------
!Subroutine: write_phase_change
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!
! This subroutine will write the new equation of state including phase changes
! where pressure is constant for a range of energy density
!
subroutine write_phase_change(output_file)
    implicit none
    real(dp), allocatable :: radius(:), mass(:), e_c(:), p_c(:), energy_density_MvsR(:)
    real(dp) :: mass_1, mass_2, epsilon_0, epsilon_step, epsilon_1, delta_epsilon, pressure, dummy_real
    integer :: i, j, k, l, n, index_MvsR, index_eq_of_state, AllocateStatus
    character(len=50), intent(out) :: output_file
    character(len=50) :: mass_radius_file, eq_of_state_input_file, dummy_character

    ! open input equation of state file
    open (unit=88, file='MvsR_new.dat', status='unknown')
    open (unit=89, file='new_eq_of_state.dat', status='unknown')

    ! open output file
    open (unit=90, file='eq_of_state_phase_change.dat', status='unknown')

    ! write file header
    write(90,*) 'New Equation of State Model Including Phase Change'
    write(90,'(2a36)') 'central energy density (MeV/fm^3)', 'central pressure (MeV/fm^3)'

    output_file = 'eq_of_state_phase_change.dat' ! unit=90

    mass_radius_file = 'MvsR_new.dat' ! unit = 88

    eq_of_state_input_file = 'new_eq_of_state.dat' ! unit = 89

    call read_p_e(mass, radius, mass_radius_file) ! put radius and mass into arrays

    ! allocate energy_density_MvsR array
    allocate(energy_density_MvsR(size(mass)), stat=AllocateStatus)
        if(AllocateStatus /= 0) stop "*** Not enough memory ***"

    ! Set iteration number i

    i = 1

    ! initialize mass comparisons
    mass_1 = mass(i)
    mass_2 = mass(i+1)

    do while (mass_2 > mass_1 )

        mass_1 = mass(i)
        mass_2 = mass(i+1)

        ! update iteration number
        i = i+1

    end do

    index_MvsR = i ! set the index at which mass begins to decrease in the MvsR plot

    print*, 'index_MvsR = ', index_MvsR

    ! write the third column (energy density) from MvsR into an array:

    ! skip first two header lines
    read(88,*) dummy_character
    read(88,*) dummy_character

    do k = 1, size(mass)

        read(88,*) dummy_real, dummy_real, energy_density_MvsR(k) ! energy_density_MvsR in MeV^4

    end do 

    print*, "energy density at index_MvsR is ", energy_density_MvsR(index_MvsR), " MeV/fm^3"

    call read_p_e(p_c, e_c, eq_of_state_input_file) ! put p and e into arrays
                        ! p and e in MeV/fm^3

    ! find index at which the value of energy density lies in this file

    l = 1

    do while (e_c(l) < energy_density_MvsR(index_MvsR))

        ! energy_density_MvsR(index_MvsR) is the energy density at which mass begins to decrease in the MvsR plot
        ! energy density in MeV/fm^3

        l = l + 1 ! update index

    end do 

    index_eq_of_state = l

    print*, 'index_eq_of_state = ', index_eq_of_state

    do j = 1, index_eq_of_state 

        ! eq of state ; e_c and p_c in MeV/fm^3
        write(90,'(2g28.16)') e_c(j), p_c(j)

    end do

    ! Now write constant pressure section
    ! First, set energy density iteration parameters

    epsilon_0 = e_c(index_eq_of_state)
    epsilon_step = epsilon_0 / 100.0_dp
    epsilon_1 = 10000.0_dp * epsilon_0 

    ! width of this: 
    delta_epsilon = epsilon_1 - epsilon_0

    do while (epsilon_0 <= epsilon_1)

        pressure = p_c(index_eq_of_state) ! set pressure (constant value)

        write(90,'(2g28.16)') epsilon_0, pressure

        ! update epsilon

        epsilon_0 = epsilon_0 + epsilon_step

    end do ! constant pressure segment complete

    ! Finish new eq of state file

    do n = index_eq_of_state, size(p_c) ! run from the index to the end of the original file for the
                             ! number of points to add to the file

        write(90,'(2g28.16)') ( e_c(n) + delta_epsilon ), p_c(n)

    end do

    ! close files
    close(88)
    close(89)
    close(90)

end subroutine write_phase_change




end module read_write