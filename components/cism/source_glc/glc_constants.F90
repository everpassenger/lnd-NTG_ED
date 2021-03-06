!=======================================================================
!BOP
!
! !MODULE: glc_constants - constants used by glc modules
!
  module glc_constants

! !DESCRIPTION:
!
! This module contains constants used by glc modules.
!
! Note that many of the required parameters are contained
! in glimmer_physcon and glimmer_params.
! Many of the parameters defined here are standard constants in POP.
!
! !REVISION HISTORY:
!  Author: William Lipscomb, LANL

! !USES:

  use glc_kinds_mod
  use shr_const_mod, only: radius=> SHR_CONST_REARTH,&
        	 	   tkfrz=>  SHR_CONST_TKFRZ

!lipscomb - Previously, stdout was defined in glc_constants.
!           Moved to glimmer_paramets so that it can be accessed from
!            glimmer source code as well as glc source code.
!           Glimmer does most of its standard output by calling the
!            write_log subroutine, which has a private output index
!            called glimmer_unit, but it is convenient sometimes to
!            write diagnostics directly to stdout.  
!           In CESM runs, glimmer_unit is set to stdout at initialization. 

  use glimmer_paramets, only: stdout
!EOP
!=======================================================================

  implicit none
  public

  include 'netcdf.inc'

   !-----------------------------------------------------------------
   ! elevation class info
   !-----------------------------------------------------------------

  logical, parameter :: verbose = .false.

  logical ::   &
     glc_smb              ! if true, get surface mass balance from CLM via coupler
                          ! (in multiple elevation classes)
                          ! if false, use PDD scheme in GLIMMER
                          ! set in glc_cpl_indices_set

   !-----------------------------------------------------------------
   !  common formats for formatted output
   !-----------------------------------------------------------------

   integer (i4), public :: &
      nml_in,            &! reserved unit for namelist input
!!      stdout,            &! reserved unit for standard output
                            ! see note above
      stderr              ! reserved unit for standard error

   character (1), parameter, public :: &
      char_delim = ','
 
   character (9), parameter, public :: &
      delim_fmt  = "(72('-'))",         &
      ndelim_fmt = "(72('='))"

   character (5), parameter, public :: &
      blank_fmt = "(' ')"

   character (char_len), public ::  &
      char_blank          ! empty character string

   !-----------------------------------------------------------------
   ! numbers
   !-----------------------------------------------------------------
 
   real (r8), parameter, public :: &
      c0     =    0.0_r8   ,&
      c1     =    1.0_r8

!EOP
!

!------------------------------------------------------------------------

  end module glc_constants

!------------------------------------------------------------------------
