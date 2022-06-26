! Program: program
! By: David Wilkins and Luke Glass
!-----------------------------------------------------------------------------
! Description:
!
! This program solves for various quantities related to the equation of state 
! of compact matter involving both quarks and leptons. The program writes the
! bag model equation of state not considering quarks with mass, as well as an 
! adjusted model now considering the mass of quarks.
!
! The equation of state is written to an output file, and 
! is iterated over different values of electron chemical potential.
! As well as equation of state values (pressure and central energy density),
! the program calculates chemical potentials of each species, as well as 
! individual number densities and total number density (per fm^3).
! 
!-----------------------------------------------------------------------------
program program
use types
use physics, only : calculate_p_e_Z, populate_arrays, TOV, read_p_e
use read_write, only : write_bag_model, write_new_eq_of_state, write_new_eq_of_state_phase_change, write_phase_change

implicit none

real(dp), allocatable :: m_f(:), k_f(:), q(:), p_work(:), e_work(:)
integer :: flavor_max
integer, allocatable :: g(:)
character(len=50) :: filename, bag_filename, output_file

! set flavor_max
flavor_max = 6 ! up, down, strange, electron, (muon) --> after a certain energy density, consider flavor_max = 6 for charm quarks

! equation of state, no phase change
call populate_arrays(flavor_max, m_f, k_f, g, q)
call write_new_eq_of_state(flavor_max, k_f, q, m_f, g)

! original bag model:
call write_bag_model(filename)

! TOV part...

filename = 'new_eq_of_state.dat'

call read_p_e(p_work, e_work, filename)
call TOV(p_work, e_work)
call write_phase_change(output_file) ! uses new_eq_of_state.dat and MvsR_new.dat to write the output file
call TOV(p_work, e_work)


end program program