! =====================================================================================
! ASCIDIAN
!
! Copyright (c) Tristan Salles (The University of Sydney)
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3.0 of the License,or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful,but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not,write to the Free Software Foundation,Inc.,59 Temple
! Place,Suite 330,Boston,MA 02111-1307 USA
! =====================================================================================
! =====================================================================================
!
!       Filename:  Interpolation.f90
!
!    Description:  Interpolation routine.
!
!        Version:  1.0
!        Created:  19/09/15 15:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module interpolation

  implicit none

  public

contains

  ! ============================================================================
  function binarysearch(length,array,value,delta)

    integer,intent(in)::length
    real,dimension(length),intent(in)::array
    real,intent(in)::value
    real,intent(in),optional::delta
    integer::binarysearch
    integer::left,middle,right
    real::d

    d=1e-9
    if(present(delta)) d=delta

    left=1
    right=length
    do
      if(left>right) exit
      middle=nint((left+right)/2.0)
      if(abs(array(middle)-value)<=d)then
        binarysearch=middle
        return
      elseif(array(middle)>value)then
        right=middle-1
      else
        left=middle+1
      endif
    enddo
    binarysearch=right

  end function binarysearch
  ! ============================================================================
  real function interpolate(x_len,x_array,y_len,y_array,f,x,y)

    integer,intent(in)::x_len,y_len
    real,dimension(x_len),intent(in)::x_array
    real,dimension(y_len),intent(in)::y_array
    real,dimension(x_len,y_len),intent(in)::f
    real,intent(in)::x,y
    real::denom,x1,x2,y1,y2
    integer::i,j

    i=binarysearch(x_len,x_array,x)
    j=binarysearch(y_len,y_array,y)
    x1=x_array(i)
    x2=x_array(i+1)
    y1=y_array(j)
    y2=y_array(j+1)
    denom=(x2-x1)*(y2-y1)
    interpolate=(f(i,j)*(x2-x)*(y2-y)+f(i+1,j)*(x-x1)*(y2-y)+ &
      f(i,j+1)*(x2-x)*(y-y1)+f(i+1,j+1)*(x-x1)*(y-y1))/denom

  end function interpolate
  ! ============================================================================

end module interpolation
