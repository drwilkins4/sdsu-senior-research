!-----------------------------------------------------------------------
!Module: math
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!!
!! Description: contains subprograms that carry out mathematical methods
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_p_e_Z
!! populate_arrays
!! TOV
!! 
!-----------------------------------------------------------------------
module math
use types

implicit none

private
public :: linear_interpolation

contains

!-----------------------------------------------------------------------
!Subroutine: linear_interpolation
!-----------------------------------------------------------------------
!! By: David Wilkins and Luke Glass
!!
!! Description:
!!
!! This subroutine takes a given x value, checks where the value falls inside
!! a given file, constructs the slope between the closest points to the given
!! x value in the file, and creates a linear interpolation between those points.
!! The given x value is then fed into that linear equation, and outputs a y value
!! to be returned corresponding to the given x value
!! 
!!----------------------------------------------------------------------
!
! Input:
!
!   x_array         1D array of x values from file
!
!   y_array         1D array of y values from file
!
!   x_0             inputted x value
!
!!----------------------------------------------------------------------
! Output:
!
! 	y_0				output value
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine linear_interpolation(x_array, y_array, x_0, y_0)
    implicit none
    real(dp), intent(in) :: x_array(:), y_array(:), x_0
    real(dp), intent(out) :: y_0
    real(dp) :: x_lower, x_upper, slope
    integer :: i, array_size, exit_stat, lower_index, upper_index
 
    array_size = size(x_array) ! size(x_array) should equal size(y_array)

    ! find values at which x_0 is between

    ! lower bound

    i = 1 ! initialize index number

    do while (x_0 <= x_lower) ! when the x_0 input becomes greater than the lower bound,
                              ! exit the loop and mark this index

        x_lower = x_array(i) ! set lower value equal to the array value of i

        i = i + 1 ! update i 

    end do 

    lower_index = i ! set lower index

    ! upper bound

    i = 1 ! reinitialize index number

    do while (x_0 >= x_upper) ! when the x_0 input becomes lower than the lower bound, 
                              ! exit the loop and mark this index

        x_upper = x_array(array_size - i) ! set upper value equal to the array value of 
                                          ! array_size - i

        i = i + 1 ! update i 

    end do 

    upper_index = array_size - i ! set lower index

    ! check if x_0 is within the bounds
    if (x_0 < x_lower) then

        print*, "x_0 is lower than the low bound. Stop"
        call exit

    end if

    if (x_0 > x_upper) then

        print*, "x_0 is greater than the high bound. Stop"
        call exit

    end if

    ! Construct slope
    ! slope = (y2-y1)/(x2-x1)
    slope = ( y_array(upper_index) - y_array(lower_index) ) /&
            ( x_array(upper_index) - x_array(lower_index) )

    ! construct line equation in point slope form to solve for y_0 using
    ! one of the indexed points

    y_0 = slope * (x_0 - x_array(lower_index) ) + y_array(lower_index)    

end subroutine linear_interpolation

!--------------------------------------------------------------------------------



end module math