!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_constants.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glad_constants

  use glimmer_global, only: dp
  use glimmer_physcon, only: pi
  
  implicit none
  
  ! ------------------------------------------------------------
  ! global parameters
  ! ------------------------------------------------------------

  real(dp), parameter :: lapse = 0.0065_dp   ! atm lapse rate, deg/m
  real(dp),parameter :: days2hours = 24.d0
  real(dp),parameter :: hours2seconds = 3600.d0   !> Hours to seconds conversion factor

  real(dp), parameter :: default_diy = 360.d0                    !> Default number of days in year
  real(dp), parameter :: default_y2h = days2hours*default_diy     !> Default years to hours conversion

  ! Constants set at run-time

  integer  :: days_in_year = default_diy        !> The number of days in a year  
  real(dp) :: years2hours = default_y2h        !> Years to hours conversion factor
  real(dp) :: hours2years = 1.d0/default_y2h !> Hours to years conversion factor

  private :: default_diy, default_y2h
  
  ! Minimum thickness of ice, at or below which a point is considered bare land for upscaling/
  ! downscaling purposes. Values other than 0 can result in odd behavior - e.g., a value
  ! greater than 0 means that CLM would consider a point to have become icesheet, and so
  ! would send positive SMB, but if this SMB didn't reach the min_thck threshold, then
  ! CISM would effectively tell CLM, "no, it's not actually icesheet yet - it's still
  ! bare ground".
  real(dp), parameter :: min_thck = 0.d0

contains

  subroutine glad_set_year_length(daysinyear)

    integer, intent(in) :: daysinyear

    days_in_year = daysinyear
    years2hours = days2hours*days_in_year 
    hours2years = 1.d0/years2hours      

  end subroutine glad_set_year_length

end module glad_constants
