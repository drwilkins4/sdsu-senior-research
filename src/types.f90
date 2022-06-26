!-----------------------------------------------------------------------
!Module: types
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! Description: Data types and global parameters are defined in this
!! module. More may be added as needed.
!-----------------------------------------------------------------------
module types

    use iso_fortran_env
    
    implicit none
    
    integer, parameter :: sp = REAL32 ! single precision
    integer, parameter :: dp = REAL64 ! double precision
    integer, parameter :: qp = REAL128! quadruple precision
    real(dp), parameter :: pi = acos(-1.0_dp)! pi
    real(dp), parameter :: conversion_factor = 197.329_dp ! MeV*fm
    real(dp), parameter :: bag = 20.0_dp ! Bag constant in MeV/fm^3
                                          ! Previously was 57 MeV/fm^3
    real(dp), parameter :: msun_mev = 1.115829d60  ! Mass of Sun in MeV
    real(dp), parameter :: msun_km = 1.475  ! Mass of Sun in km
    real(dp), parameter :: c18 =1.0d18         ! 1.0 x 10^18
    
end module types