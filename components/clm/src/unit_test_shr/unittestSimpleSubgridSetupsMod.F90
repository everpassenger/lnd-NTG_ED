module unittestSimpleSubgridSetupsMod

  ! This module provides wrappers to unittestSubgridMod, which give you a variety of
  ! simple subgrid setups.
  !
  ! Note that these routines do everything needed with the subgrid setup. So once you
  ! call these routines, you cannot add any more gridcells, landunits, etc.

  use unittestSubgridMod
  use shr_kind_mod , only : r8 => shr_kind_r8
  use landunit_varcon, only : istsoil

  implicit none
  private
  save

  ! ------------------------------------------------------------------------
  ! Routines that do everything needed with the subgrid setup, including the begin & end
  ! call. Once you call these routines, you cannot add any more gridcells, landunits, etc.
  ! ------------------------------------------------------------------------

  ! Create a grid that has a single gridcell with a single vegetated patch
  public :: setup_single_veg_patch

  ! Create a grid that has N grid cells, each with a single vegetated patch
  public :: setup_ncells_single_veg_patch

  ! ------------------------------------------------------------------------
  ! Routines that create a single grid cell with certain properties. You can do other
  ! subgrid setup (creating other grid cells) before and after this.
  ! ------------------------------------------------------------------------

  ! Create a grid cell that is 100% natural veg, with a single patch
  public :: create_gridcell_single_veg_patch

contains

  ! ========================================================================
  ! Routines that do everything needed with the subgrid setup, including the begin & end
  ! call. Once you call these routines, you cannot add any more gridcells, landunits, etc.
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine setup_single_veg_patch(pft_type)
    !
    ! !DESCRIPTION:
    ! Create a grid that has a single gridcell with a single vegetated patch, with veg
    ! type given by the pft_type argument
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: pft_type  ! the type of the single vegetated patch
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'setup_single_veg_patch'
    !-----------------------------------------------------------------------
    
    call setup_ncells_single_veg_patch(ncells=1, pft_type=pft_type)

  end subroutine setup_single_veg_patch

  !-----------------------------------------------------------------------
  subroutine setup_ncells_single_veg_patch(ncells, pft_type)
    !
    ! !DESCRIPTION:
    ! Create a grid that has ncells grid cells, each with a single vegetated patch. All
    ! vegetated patches have the same type, given by pft_type.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: ncells    ! number of grid cells
    integer, intent(in) :: pft_type  ! pft type
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'setup_ncells_single_veg_patch'
    !-----------------------------------------------------------------------

    call unittest_subgrid_setup_start()
    do i = 1, ncells
       call create_gridcell_single_veg_patch(pft_type = pft_type)
    end do
    call unittest_subgrid_setup_end()

  end subroutine setup_ncells_single_veg_patch


  ! ========================================================================
  ! Routines that create a single grid cell with certain properties. You can do other
  ! subgrid setup (creating other grid cells) before and after this.
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine create_gridcell_single_veg_patch(pft_type)
    !
    ! !DESCRIPTION:
    ! Create a grid cell that is 100% natural veg, with a single patch
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: pft_type  ! the type of the single vegetated patch
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'create_gridcell_single_veg_patch'
    !-----------------------------------------------------------------------

    call unittest_add_gridcell()
    call unittest_add_landunit(my_gi=gi, ltype=istsoil, wtgcell=1.0_r8)
    call unittest_add_column(my_li=li, ctype=1, wtlunit=1.0_r8)
    call unittest_add_patch(my_ci=ci, ptype=pft_type, wtcol=1.0_r8)

  end subroutine create_gridcell_single_veg_patch


end module unittestSimpleSubgridSetupsMod
