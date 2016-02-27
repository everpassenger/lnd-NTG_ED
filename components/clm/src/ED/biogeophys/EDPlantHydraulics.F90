module EDPlantHydraulicsMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! 
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8;
  use clm_varctl       , only : use_ed_planthydraulics
  use pftconMod        , only : pftcon
  use EDTypesMod       , only : ed_patch_type, ed_cohort_type, numpft_ed
  use EDEcophysContype , only : EDecophyscon
  use GridcellType     , only : grc
  use ColumnType       , only : col
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Hydraulics_Drive
  ! added by Brad Christoffersen Jan 2016 for use in ED hydraulics
  !    van Genuchten (1980)-specific functions for the swc (soil water characteristic)
  !    and for the kunsat (unsaturated hydraulic conductivity) curves
  !
  ! BOC...for quick implementation avoided JT's abstract interface,
  !    but these should be converted to interfaces in the future
  public :: swc_psi_VG          !returns water potential, psi
  public :: swc_satfrac_VG      !returns saturation fraction, s
  public :: swc_dpsidth_VG      !derivative of psi wrt theta
  public :: unsatk_flc_VG       !returns k/ksat ('fractional loss of conductivity' flc)
  public :: unsatk_dflcdpsi_VG  !derivative of k/ksat (flc) wrt psi
  !
  !
  !type(ed_cohort_type), pointer   :: currentCohort   ! current cohort   !! are these needed here ???
  !type(ed_patch_type) , pointer   :: currentPatch    ! current patch    !! are these needed here ???
  integer             , parameter :: npool_leaf  = 1                      ! 
  integer             , parameter :: npool_stem  = 1                      ! 
  integer             , parameter :: npool_troot = 1                      ! 
  integer             , parameter :: npool_aroot = 1                      ! 
  integer             , parameter :: npool_ag    = npool_leaf+npool_stem  ! number of aboveground plant water storage nodes
  integer             , parameter :: npool_bg    = npool_troot            ! number of belowground plant water storage nodes (except absorbing roots)
  integer             , parameter :: nshell      = 11                     ! number of concentric soil cylinders surrounding absorbing root
  integer             , parameter :: npool_tot   = npool_ag + 2 + nshell
  integer             , parameter :: n_porous_media = 5
  integer                         :: porous_media(npool_tot)              ! 1=leaf, 2=stem, 3=troot, 4=aroot, 5=soil
  
  ! 01/18/16: Created by Brad Christoffersen
  !------------------------------------------------------------------------------
   
contains 
  !------------------------------------------------------------------------------

  subroutine hydraulics_initialize()
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    !
    integer :: k !index
    
    do k = 1,npool_tot
       if(k <= npool_leaf) then
          porous_media(k) = 1
       else if(k <= (npool_leaf+npool_stem)) then
          porous_media(k) = 2
       else if(k <= (npool_leaf+npool_stem+npool_troot)) then
          porous_media(k) = 3
       else if(k <= (npool_leaf+npool_stem+npool_troot+npool_aroot)) then
          porous_media(k) = 4
       else
          porous_media(k) = 5
       end if
    enddo
    
  end subroutine hydraulics_initialize

  !------------------------------------------------------------------------------
  subroutine hydraulics_drive( bounds, ed_allsites_inst, &
             soilstate_inst, waterflux_inst, waterstate_inst, &
	     temperature_inst, energyflux_inst)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use shr_const_mod      , only : shr_const_pi  
    use decompMod          , only : bounds_type
    use clm_time_manager   , only : get_step_size
    use clm_varpar         , only : nlevgrnd
    use clm_varctl         , only : iulog
    use clm_varcon         , only : tfrz, denice, denh2o
    use SoilStateType      , only : soilstate_type
    use WaterFluxType      , only : waterflux_type
    use WaterStateType     , only : waterstate_type
    use TemperatureType    , only : temperature_type
    use EnergyFluxType     , only : energyflux_type
    use PatchType          , only : patch
    use EDTypesMod         , only : ed_site_type, map_clmpatch_to_edpatch 
    !
    ! !ARGUMENTS    
    type(bounds_type)      , intent(in)            :: bounds  ! clump bounds
    type(ed_site_type)     , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    type(waterflux_type)   , intent(inout)         :: waterflux_inst
    type(waterstate_type)  , intent(in)            :: waterstate_inst
    type(temperature_type) , intent(in)            :: temperature_inst
    type(energyflux_type)  , intent(inout)         :: energyflux_inst
    !

    SELECT CASE (use_ed_planthydraulics)

       CASE (1)
	       
	  call hydraulics_BC(bounds, ed_allsites_inst(bounds%begg:bounds%endg), &
               soilstate_inst, waterflux_inst, waterstate_inst, temperature_inst, energyflux_inst)

       CASE (2)
	       
          !call Hydraulics_CX()
		  
       CASE DEFAULT
	     
    end SELECT
	     
  end subroutine Hydraulics_Drive
  
  !-------------------------------------------------------------------------------!
  subroutine hydraulics_bc ( bounds, ed_allsites_inst, &
             soilstate_inst, waterflux_inst, waterstate_inst, &
	     temperature_inst, energyflux_inst)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use shr_const_mod      , only : shr_const_pi  
    use decompMod          , only : bounds_type
    use clm_time_manager   , only : get_step_size
    use clm_varpar         , only : nlevgrnd, nlevsoi
    use clm_varctl         , only : iulog
    use clm_varcon         , only : tfrz, denice, denh2o, grav
    use SoilStateType      , only : soilstate_type
    use WaterFluxType      , only : waterflux_type
    use WaterStateType     , only : waterstate_type
    use TemperatureType    , only : temperature_type
    use EnergyFluxType     , only : energyflux_type
    use PatchType          , only : patch
    use EDTypesMod         , only : ed_site_type, ed_patch_type, ed_cohort_type, map_clmpatch_to_edpatch 
    !
    ! !ARGUMENTS    
    type(bounds_type)      , intent(in)            :: bounds  ! clump bounds
    type(ed_site_type)     , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    type(waterflux_type)   , intent(inout)         :: waterflux_inst
    type(waterstate_type)  , intent(in)            :: waterstate_inst
    type(temperature_type) , intent(in)            :: temperature_inst
    type(energyflux_type)  , intent(inout)         :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: iv !leaf layer
    integer :: g  !gridcell
    integer :: c  !column
    integer :: p  !patch
    integer :: j  !soil layer
    integer :: k  !1D plant-soil continuum array
    integer :: ft ! plant functional type index
    !----------------------------------------------------------------------

    type (ed_patch_type),  pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    integer, parameter    :: nv = 5           ! Number of canopy layers

    ! hydraulics global constants
    real(r8), parameter :: thresh         = 1.e-5_r8     ! threshold for water balance error                    [mm h2o]
    real(r8), parameter :: rs1            = 0.001_r8     ! mean absorbing fine root radius                      [m]
    real(r8), parameter :: rs2            = 0.001_r8     ! mean fine root radius for biomass estimates          [m]
    real(r8), parameter :: aral           = 1._r8        ! root:leaf area ratio                                 [m2 m-2]
    real(r8), parameter :: resisfrac_stem = 0.625_r8     ! fraction of total resistance in stem                 [-]
    real(r8), parameter :: rhoc_stem      = 1.54_r8      ! dry matter (or cell wall) density                    [g cm-3]
    
    ! hydraulics timestep adjustments for acceptable water balance error
    integer  :: maxiter        = 10           ! maximum iterations for timestep reduction                       [-]
    integer  :: imult          = 2            ! iteration index multiplier                                      [-]
    real(r8) :: wh_tot_err                    ! 1D plant-soil continuum water error                             [kgh2o]
    real(r8) :: wh_tot_err_1l                 ! 1D plant-soil continuum water error (single-layer soln)         [kgh2o]
    
    ! hydraulics outputs from 1D solver
    real(r8) :: dwh_leaf                      ! del(total water storage) in leaves                              [kgh2o/indiv]
    real(r8) :: dwh_stem                      ! del(total water storage) in stem                                [kgh2o/indiv]
    real(r8) :: dwh_troot                     ! del(total water storage) in transporting root                   [kgh2o/indiv]
    real(r8) :: dwh_aroot_tot                 ! del(total water storage) in absorbing roots summed across layers[kgh2o/indiv]
    real(r8) :: dwh_rhiz_tot(nshell)          ! del(total water storage) in 'representative' rhizosphere shell  [kgh2o/indiv/shell]
    real(r8) :: dw_soil_tot_coh               ! del(total water storage) across all rhiz shells (by cohort)     [kgh2o/indiv/timestep]

    ! cohort-specific arrays to hold 1D hydraulics geometric & state variables for entire continuum (leaf,stem,root,soil)
    integer  :: tisstype(     npool_tot)      ! tissue type (1=leaf, 2=stem, 3=troot, 4=aroot)                  [-]
    real(r8) :: z_node(       npool_tot)      ! nodal height of water storage compartments                      [m]
    real(r8) :: z_node_1l(    npool_tot)      ! nodal height of water storage compartments (single-layer soln)  [m]
    real(r8) :: z_upper(      npool_tot)      ! upper boundary height (more proximal to atmosphere)             [m]
    real(r8) :: z_lower(      npool_tot)      ! lower boundary height (more distal to atmosphere)               [m]
    real(r8) :: v_node(       npool_tot)      ! volume of water storage compartments                            [m3]
    real(r8) :: v_node_1l(    npool_tot)      ! volume of water storage compartments (single-layer soln)        [m3]
    real(r8) :: w_node(       npool_tot)      ! mass of water in water storage compartments                     [kgh2o/indiv]
    real(r8) :: psi_node(     npool_tot)      ! water potential in water storage compartments                   [MPa]
    real(r8) :: ths_node(     npool_tot)      ! saturated volumetric water in water storage compartments        [m3 m-3]
    real(r8) :: ths_node_1l(  npool_tot)      ! saturated volumetric water in water storage compartments (single-layer soln) [m3 m-3]
    real(r8) :: thr_node(     npool_tot)      ! residual volumetric water in water storage compartments         [m3 m-3]
    real(r8) :: thr_node_1l(  npool_tot)      ! residual volumetric water in water storage compartments (single-layer soln) [m3 m-3]
    real(r8) :: the_node(     npool_tot)      ! error resulting from supersaturation or below-residual th_node  [m3 m-3]
    real(r8) :: the_node_1l(  npool_tot)      ! like the_node(:) but for specific single soil layer             [m3 m-3]
    real(r8) :: th_node(      npool_tot)      ! volumetric water in water storage compartments                  [m3 m-3]
    real(r8) :: th_node_1l(   npool_tot)      ! volumetric water in water storage compartments (single-layer soln) [m3 m-3]
    real(r8) :: dth_node(     npool_tot)      ! change in volumetric water in water storage compartments        [m3 m-3]
    real(r8) :: dth_node_1l(  npool_tot)      ! like dth_node_1l(:) but for specific single soil layer          [m3 m-3]
    real(r8) :: kmax_bound(   npool_tot)      ! lower boundary maximum hydraulic conductance of compartments    [kg s-1 MPa-1]
    real(r8) :: kmax_bound_1l(npool_tot)      ! lower boundary maximum hydraulic conductance of compartments (single-layer soln) [kg s-1 MPa-1]
    real(r8) :: l_aroot_tot_coh               ! total length of absorbing roots across all soil layers (cohort) [m]
    
    ! column-specific arrays to hold rhizosphere geometric & state variables
    real(r8) :: dz_tot                        ! total soil depth (to bottom of bottom layer)                    [m]
    real(r8) :: l_aroot_tot_col               ! total length of absorbing roots across all soil layers          [m]
    real(r8) :: r_out_shell_1D(nshell)        ! outer radius of rhizosphere compartment                         [m]
    real(r8) :: r_node_shell_1D(nshell)       ! nodal radius of rhizosphere compartment                         [m]
    real(r8) :: dth_layershell_col(nlevsoi,nshell) ! accumulated water content change over all cohorts in a column   [m3 m-3]
    real(r8) :: ths_shell_1D(nshell)          ! saturated water content of rhizosphere compartment              [m3 m-3]
    real(r8) :: thr_shell_1D(nshell)          ! residual water content of rhizosphere compartment               [m3 m-3]
    real(r8) :: v_shell_1D(nshell)            ! shell volume of rhizosphere compartment                         [m3]
    real(r8) :: kmax_bound_shell_1D(nshell)   ! maximum conductance at upper (closer to atm) shell boundaries   [kg s-1 MPa-1]
    real(r8) :: kmax_bound_shell_1l(nshell)   ! like kmax_bound_shell_1D(:) but for specific single soil layer  [kg s-1 MPa-1]
    real(r8) :: psi_node_shell_1D(nshell)     ! soil matric potential of rhizosphere compartment                [MPa]
    real(r8) :: ths_aroot_1D                  ! saturated water content of 1D representation of fine roots      [m3 m-3]
    real(r8) :: thr_aroot_1D                  ! residual water content of 1D representation of fine roots       [m3 m-3]
    real(r8) :: vtot_aroot_1D                 ! sum of fine root volume across soil layers                      [m3]
    real(r8) :: psi_node_aroot_1D             ! water potential of absorbing root                               [MPa]

    ! hydraulics architecture & allometry variables
    real(r8) :: leafmass_tot                  ! cohort total leaf dry mass                                      [kg]
    real(r8) :: leafarea_tot                  ! cohort total leaf area                                          [m2]
    real(r8) :: sapwoodarea                   ! sapwood area                                                    [m2]
    real(r8) :: rootarea_tot                  ! cohort total fine root surface absorbing area                   [m2]
    real(r8) :: rootlength_tot                ! cohort total fine root length                                   [kg]
    real(r8) :: rootlength(nlevgrnd)          ! fine root length in each layer                                  [m]
    real(r8) :: sf_stem_to_leaf               ! conductivity to conductance scale factor for stem-canopy pathway[m]
    real(r8) :: sf_troot_to_stem              ! conductivity to conductance scale factor for troot-stem pathway [m]
    
    ! hydraulics conductances
    real(r8) :: kmax_bound_bylayershell(nlevsoi,nshell) ! maximum conductance at shell boundaries in each layer [kg s-1 MPa-1]
    real(r8) :: ksoil_bylayer(nlevsoi)        ! total rhizosphere conductance (over all shells) by soil layer   [MPa]
    real(r8) :: ksoil_tot                     ! total rhizosphere conductance (over all shells and soil layers  [MPa]
    real(r8) :: kmax_aroot_radial_in          ! maximum aroot (radial,incoming) conductance                     [kg s-1 MPa-1]
    real(r8) :: kmax_aroot_radial_out         ! maximum aroot (radial,outgoing) conductance                     [kg s-1 MPa-1]
    real(r8) :: kmax_aroot_radial             ! maximum aroot (radial,bidirectional) conductance                [kg s-1 MPa-1]
    real(r8) :: kmax_stem                     ! maximum whole-stem (above troot to leaf) conductance            [kg s-1 MPa-1]
    real(r8) :: kmax_aroot_leaf               ! maximum whole-tree (surface of aroot to leaf) conductance       [kg s-1 MPa-1]
    
    ! hydraulics other
    real(r8) :: qflx_tran_veg_indiv           ! individiual transpiration rate                                  [kgh2o indiv-1 s-1]
    real(r8) :: gpp_clm_patch                 ! sum of gpp_clm across all cohorts within a patch                [kgC patch-1 timestep-1]
    real(r8) :: dtime                         ! timestep size                                                   [s]
    integer  :: cprev                         ! previous column index
    integer  :: ncoh_col                      ! number of cohorts across all non-veg patches within a column

    !------------------------------------------------------------------------------

    associate(& 
         dz                  => col%dz                              , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         z                   => col%z                               , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         v_shell             => col%v_shell                         , & ! Input:  [real(r8) (:,:,:)] volume of rhizosphere compartment (m)
         r_node_shell        => col%r_node_shell                    , & ! Input:  [real(r8) (:,:,:)] nodal radius of rhizosphere compartment (m)

         smpso               => pftcon%smpso                        , & ! Input:  soil water potential at full stomatal opening (mm)
         smpsc               => pftcon%smpsc                        , & ! Input:  soil water potential at full stomatal closure (mm)

         sucsat              => soilstate_inst%sucsat_col           , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm) 
         watsat              => soilstate_inst%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         watdry              => soilstate_inst%watdry_col           , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=0
         watopt              => soilstate_inst%watopt_col           , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran = 1
         bsw                 => soilstate_inst%bsw_col              , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b" 
         soilbeta            => soilstate_inst%soilbeta_col         , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity 
         sand                => soilstate_inst%sandfrac_patch       , & ! Input:  [real(r8) (:)   ]  % sand of soil 
         rootr               => soilstate_inst%rootr_patch          , & ! Output: [real(r8) (:,:) ]  Fraction of water uptake in each layer

         h2osoi_ice          => waterstate_inst%h2osoi_ice_col      , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_vol          => waterstate_inst%h2osoi_vol_col      , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] 
         h2osoi_liq          => waterstate_inst%h2osoi_liq_col      , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) 
         h2osoi_vol_shell    => waterstate_inst%h2osoi_vol_shell    , & ! Input:  [real(r8) (:,:,:)] volumetric water in rhizosphere compartment (m3/m3)

         t_soisno            => temperature_inst%t_soisno_col       , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)

         btran               => energyflux_inst%btran_patch         , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)
         btran2              => energyflux_inst%btran2_patch        , & ! Output: [real(r8) (:)   ] 
         rresis              => energyflux_inst%rresis_patch        , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd) 

         qflx_tran_veg_col   => waterflux_inst%qflx_tran_veg_col    , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)  
         qflx_tran_veg_patch => waterflux_inst%qflx_tran_veg_patch    & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         )

       !set timestep & initialization
       dtime           = get_step_size()
       wh_tot_err      = 0.0_r8
       cprev           = -9999

       ! do f = 1, fn
       !    p = filterp(f)
       !    if (patch%is_veg(p)) then
       !       c = patch%column(p)
       !       if((c \= cprev)) then
       !       end if
       
       do c = bounds%begc, bounds%endc
          
   	  ! AVERAGE ROOT WATER UPTAKE (BY RHIZOSPHERE SHELL) ACROSS ALL COHORTS WITHIN A COLUMN
	  dth_layershell_col(:,:)   = 0._r8
	  ncoh_col                  = 0
	  
	  do p = col%patchi(c), col%patchf(c)

             if (patch%is_veg(p)) then

                g            =  patch%gridcell(p)
                currentPatch => map_clmpatch_to_edpatch(ed_allsites_inst(g), p) 
	 
    	        gpp_clm_patch = 0.0_r8
	        currentCohort=>currentPatch%tallest
                do while(associated(currentCohort))
	           gpp_clm_patch =  gpp_clm_patch + currentCohort%gpp_clm * currentCohort%n    ! kgC/indiv/timestep * n indiv = kgC/cohort/timestep
                   currentCohort => currentCohort%shorter
                enddo !cohort
	     
	        currentCohort=>currentPatch%tallest
    	        do while(associated(currentCohort))
	           ft       = currentCohort%pft
	           ncoh_col = ncoh_col + 1
	    
	           qflx_tran_veg_indiv = qflx_tran_veg_patch(p) * currentPatch%area * currentCohort%gpp_clm/gpp_clm_patch
		  
		   ! BUCKET APPROXIMATION OF THE SOIL-ROOT HYDRAULIC GRADIENT (weighted average across layers)
                   call map2d_to_1d_shells(soilstate_inst, waterstate_inst, g, c, rs1, dz_tot, currentCohort%l_aroot_layer, r_out_shell_1D, r_node_shell_1D, ths_shell_1D, &
		                           thr_shell_1D, v_shell_1D, psi_node_shell_1D, &
					   ksoil_bylayer, ksoil_tot, kmax_bound_bylayershell)
                   call boundary_kmax_rhizosphere(soilstate_inst, c, 1, 1, nshell, currentCohort%l_aroot_tot, &
		                                  rs1, r_node_shell_1D, kmax_bound_shell_1D)
	           psi_node(  (npool_tot-nshell+1):npool_tot) = psi_node_shell_1D(:)
         
		   ! REPRESENTATIVE SINGLE FINE ROOT POOL (weighted average across layers)
		   call map2d_to_1d_aroot(ft, currentCohort, kmax_bound_bylayershell, ths_aroot_1D, thr_aroot_1D, vtot_aroot_1D, psi_node_aroot_1D)
		   psi_node(npool_ag+npool_bg+1)            = psi_node_aroot_1D

		   ! SET NODE HEIGHTS AND VOLUMES
		   z_node(                   1 : npool_ag)           = currentCohort%z_node_ag(:)
		   z_node(         (npool_ag+1):(npool_ag+npool_bg)) = currentCohort%z_node_bg(:)
		   z_node((npool_ag+npool_bg+1): npool_tot)          = z_node(npool_ag+npool_bg)
		   v_node(                   1 : npool_ag)           = currentCohort%v_ag(:)
		   v_node(         (npool_ag+1):(npool_ag+npool_bg)) = currentCohort%v_bg(:)
		   v_node((npool_ag+npool_bg+1)                    ) = vtot_aroot_1D
	           v_node( (npool_tot-nshell+1): npool_tot)          = v_shell_1D(:)/dz_tot*currentCohort%l_aroot_tot

		   ! SET SATURATED & RESIDUAL WATER CONTENTS (leaf, stem, troot)
	           ths_node(  (npool_tot-nshell+1):npool_tot) = ths_shell_1D(:)
	           thr_node(  (npool_tot-nshell+1):npool_tot) = thr_shell_1D(:)
		   ths_node(  (npool_tot-nshell)            ) = ths_aroot_1D
		   thr_node(  (npool_tot-nshell)            ) = thr_aroot_1D
		   do k=1,npool_ag+npool_bg
		      ths_node(k) = pftcon%thetas_node(ft,porous_media(k))
		      thr_node(k) = pftcon%thetas_node(ft,porous_media(k)) * pftcon%resid_node(ft,porous_media(k))
		   enddo

		   ! SET BOUNDARY MAX CONDUCTANCES
		   !! assign cohort-level conductances to the 1D array
		   kmax_bound(          1 : npool_ag)         = currentCohort%kmax_bound
		   !! estimate troot-aroot and aroot-radial components as a residual:
		   !! 25% each of total (surface of aroots to leaves) resistance
		   kmax_bound((npool_ag+1):(npool_ag+2))      = 2._r8 * currentCohort%kmax_treebg_tot
		   kmax_bound(              npool_ag+2)       = 1._r8/(1._r8/kmax_bound(npool_ag+2) + 1._r8/kmax_bound_shell_1D(1))
		   kmax_bound((npool_tot-nshell+2):npool_tot) = kmax_bound_shell_1D(2:nshell)
		     
		   ! CONVERT WATER POTENTIALS TO WATER CONTENTS FOR THE NEW 'BUCKET' RHIZOSPHERE (fine roots and rhizosphere shells)
		   do k = (npool_tot - nshell), npool_tot
		      call th_from_psi(soilstate_inst, ft, porous_media(k), psi_node(k), th_node(k), c)
		   enddo !aroot thru outer rhiz shell

		   ! MAP REMAINING WATER CONTENTS (leaf, stem, troot) TO THE 1D ARRAY
		   th_node(          1 : npool_ag          ) = currentCohort%th_ag
		   th_node((npool_ag+1):(npool_ag+npool_bg)) = currentCohort%th_bg
		
                   ! AGGREGATE 1-D THETA-BASED SOLUTION TO RICHARDS' EQUATION
	   	   call Hydraulics_1DSolve(soilstate_inst, ft, c, z_node, v_node, ths_node, thr_node, kmax_bound, &
		                           th_node, qflx_tran_veg_indiv, &
                                           thresh, maxiter, imult, dtime, &
		                           dth_node, the_node, wh_tot_err)

		   ! UPDATE WATER CONTENT & POTENTIAL IN LEAVES, STEM, AND TROOT (COHORT-LEVEL)
		   do k=1,npool_ag
		      currentCohort%th_ag(k)          = th_node(k)
		      call psi_from_th(soilstate_inst, ft, porous_media(k), currentCohort%th_ag(k), currentCohort%psi_ag(k), c)
		   enddo
		   do k=(npool_ag+1),(npool_ag+npool_bg)
		      currentCohort%th_bg(k-npool_ag) = th_node(k)
		      call psi_from_th(soilstate_inst, ft, porous_media(k), currentCohort%th_bg(k), currentCohort%psi_bg(k), c)
		   enddo
		
		   ! UPDATE BTRAN BASED ON NEW LEAF WATER POTENTIAL (COHORT-LEVEL)
		   currentCohort%btran(:) = (1._r8 + (currentCohort%psi_ag(1)/pftcon%p50_gs(ft))**pftcon%avuln_gs(ft))**(-1._r8)
		  
		   ! LAYER-SPECIFIC 1-D THETA-BASED SOLUTION TO RICHARDS' EQUATION
		   !! partition dw_soil_layer_coh by shell in proportion to each shell's contribution to dW_richards_1layer,
		   !! where dW_richards_1layer is the solution to the 1D richards equation for 1 layer
		   do j=1,nlevsoi
		      z_node_1l(                       1 :(npool_ag+npool_bg)) = z_node(1:(npool_ag+npool_bg))
		      z_node_1l(    (npool_ag+npool_bg+1): npool_tot         ) = z(c,j)
		   
		      v_node_1l(                       1 :(npool_ag+npool_bg)) = v_node(1:(npool_ag+npool_bg))
		      v_node_1l(    (npool_ag+npool_bg+1)                    ) = currentCohort%v_aroot_layer(j)
		      v_node_1l(     (npool_tot-nshell+1): npool_tot         ) = v_shell(c,j,:)/dz(c,j)*currentCohort%l_aroot_layer(j)
		   
		      !! FIX (BOC): make these layer-specific for the soil nodes
		      ths_node_1l(                     1 : npool_tot         ) = ths_node(1:npool_tot)
		      thr_node_1l(                     1 : npool_tot         ) = thr_node(1:npool_tot)
		   
                      call boundary_kmax_rhizosphere(soilstate_inst, c, j, 1, nshell, currentCohort%l_aroot_layer(j), &
		                                     rs1, r_node_shell(c,j,:), kmax_bound_shell_1l)
		      kmax_bound_1l(                   1 : npool_ag          ) = kmax_bound(1:npool_ag)
		      kmax_bound_1l((npool_ag+1)         :(npool_ag+2)       ) = 2._r8 * currentCohort%kmax_treebg_layer(j)
		      kmax_bound_1l((npool_ag+2)                             ) = 1._r8/(1._r8/kmax_bound_1l(npool_ag+2) + 1._r8/kmax_bound_shell_1l(1))
		      kmax_bound_1l((npool_tot-nshell+2) : npool_tot         ) = kmax_bound_shell_1l(2:nshell)
		   
		      th_node_1l(                      1 :(npool_ag+npool_bg)) = th_node(1:(npool_ag+npool_bg))
		      th_node_1l(   (npool_ag+npool_bg+1)                    ) = currentCohort%th_aroot(j)
		      th_node_1l(   (npool_ag+npool_bg+2): npool_tot         ) = h2osoi_vol_shell(c,j,:)
		   
 	  	      call Hydraulics_1DSolve(soilstate_inst, ft, c, z_node_1l, v_node_1l, ths_node_1l, thr_node_1l, kmax_bound_1l, &
		                              th_node_1l, qflx_tran_veg_indiv*ksoil_bylayer(j)/ksoil_tot, &
                                              thresh, maxiter, imult, dtime, &
	  	                              dth_node_1l, the_node_1l, wh_tot_err_1l)
                   
		      ! UPDATE WATER CONTENT & POTENTIAL IN AROOT (COHORT-LEVEL)
		      currentCohort%th_aroot(j) = th_node_1l(npool_ag+npool_bg+1)
		      call psi_from_th(soilstate_inst, ft, porous_media(npool_ag+npool_bg+1), currentCohort%th_aroot(j), currentCohort%psi_aroot(j), c)

		      ! ACCUMULATE WATER CONTENT CHANGE IN SOIL (COLUMN-LEVEL)
		      dth_layershell_col(j,:)   = dth_layershell_col(j,:) + dth_node_1l((npool_tot-nshell+1):npool_tot)
		   enddo

                   currentCohort => currentCohort%shorter

                enddo !cohort 

             end if !vegetated patch

          enddo !patch
	  
          ! AVERAGE THE ACCUMULATED WATER CONTENT CHANGE OVER ALL COHORTS IN A COLUMN
          dth_layershell_col(:,:) = dth_layershell_col(:,:) / ncoh_col
	  h2osoi_vol_shell(c,:,:) = h2osoi_vol_shell(c,:,:) + dth_layershell_col(:,:)

       enddo !column
      
    end associate

  end subroutine Hydraulics_BC
  
  !-------------------------------------------------------------------------------!
  subroutine Hydraulics_1DSolve(soilstate_inst, ft, c, z_node, v_node, ths_node, thr_node, kmax_bound, &
                                th_node, qtop, &
                                thresh, maxiter, imult, dtime, &
				dth_node, the_node, wh_tot_err)
    !
    ! !DESCRIPTION: 
    !
    ! !USES:
    use clm_varcon         , only : denh2o
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    integer     , intent(in)      :: ft                   ! PFT index
    integer     , intent(in)      :: c                    ! column index
    real(r8)    , intent(in)      :: z_node(:)            ! nodal height of water storage compartments                      [m]
    real(r8)    , intent(in)      :: v_node(:)            ! volume of water storage compartments                            [m3]
    real(r8)    , intent(in)      :: ths_node(:)          ! saturated volumetric water in water storage compartments        [m3 m-3]
    real(r8)    , intent(in)      :: thr_node(:)          ! residual volumetric water in water storage compartments         [m3 m-3]
    real(r8)    , intent(in)      :: kmax_bound(:)        ! lower boundary maximum hydraulic conductance of compartments    [kg s-1 MPa-1]
    real(r8)    , intent(inout)   :: th_node(:)           ! volumetric water in water storage compartments                  [m3 m-3]
    real(r8)    , intent(in)      :: qtop                 ! evaporative flux from canopy                                    [kgh2o indiv-1 s-1]
    integer     , intent(in)      :: maxiter              ! maximum iterations for timestep reduction                       [-]
    integer     , intent(in)      :: imult                ! iteration index multiplier                                      [-]
    real(r8)    , intent(in)      :: thresh               ! threshold for water balance error                               [mm h2o]
    real(r8)    , intent(in)      :: dtime                ! timestep size                                                   [s]
    real(r8)    , intent(out)     :: dth_node(:)          ! change in volumetric water in water storage compartments        [m3 m-3]
    real(r8)    , intent(out)     :: the_node(:)          ! error resulting from supersaturation or below-residual th_node  [m3 m-3]
    real(r8)    , intent(out)     :: wh_tot_err           ! 1D plant-soil continuum water error                             [kgh2o]
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: small_num      = 1.e-13_r8    ! keeps th_node within [thr+small_num, ths-small_num]  [m3 m-3]
    integer  :: k                             ! 1D plant-soil continuum array
    integer  :: iterh1, iterh2                ! iteration indices                                               [-]
    real(r8) :: wh_tot_beg, wh_tot_end        ! 1D plant-soil continuum total water storage                     [kgh2o]
    integer  :: dt_fac                        ! timestep divisor                                                [-]
    integer  :: dt_new                        ! new timestep                                                    [s]
    real(r8) :: th_node_init( npool_tot)      ! initial volumetric water in water storage compartments          [m3 m-3]
    real(r8) :: psi_node(     npool_tot)      ! water potential in water storage compartments                   [MPa]
    real(r8) :: dpsidth_node( npool_tot)      ! derivative of water potential wrt to theta                      [MPa]
    real(r8) :: flc_node(     npool_tot)      ! fractional loss of conductivity at water storage nodes          [-]
    real(r8) :: dflcdpsi_node(npool_tot)      ! derivative of fractional loss of conductivity wrt psi           [MPa-1]
    real(r8) :: k_bound(      npool_tot)      ! lower boundary hydraulic conductance of compartments            [kg s-1 MPa-1]
    real(r8) :: q_bound(      npool_tot)      ! lower boundary flux rate                                        [kg s-1]
    real(r8) :: hdiff_bound(  npool_tot)      ! total water potential difference across lower boundary          [MPa-1]
    real(r8) :: dkbounddpsi0( npool_tot)      ! derivative of lower boundary conductance wrt psi above          [kg s-1 MPa-2]
    real(r8) :: dkbounddpsi1( npool_tot)      ! derivative of lower boundary conductance wrt psi below          [kg s-1 MPa-2]
    real(r8) :: dqbounddpsi0( npool_tot)      ! derivative of lower boundary flux rate wrt psi above            [kg s-1 MPa-1]
    real(r8) :: dqbounddpsi1( npool_tot)      ! derivative of lower boundary flux rate wrt psi below            [kg s-1 MPa-1]
    real(r8) :: dqbounddth0(  npool_tot)      ! derivative of lower boundary flux rate wrt theta above          [kg s-1 m3 m-3]
    real(r8) :: dqbounddth1(  npool_tot)      ! derivative of lower boundary flux rate wrt theta below          [kg s-1 m3 m-3]
    real(r8) :: amx(npool_tot)                ! "a" left off diagonal of tridiagonal matrix                     [kg s-1]
    real(r8) :: bmx(npool_tot)                ! "b" diagonal of tridiagonal matrix                              [kg s-1]
    real(r8) :: cmx(npool_tot)                ! "c" right off diagonal of tridiagonal matrix                    [kg s-1]
    real(r8) :: rmx(npool_tot)                ! "r" forcing term of tridiagonal matrix                          [kg s-1]
    real(r8) :: th_prev                       ! temporary                                                       [m3 m-3]
    real(r8) :: dth_prev                      ! temporary                                                       [m3 m-3]
    !----------------------------------------------------------------------

    ! STORE INITIAL STATES
    !! in case timestep needs to be chopped in half to balance water
    th_node_init(:) = th_node(:)
		
    ! OUTER DO-WHILE LOOP
    !! cuts timestep in half until all sub-timesteps (inner do-while loop) balance the water budget
    iterh1 = 0
    do while(  ((iterh1 == 0) .or. (abs(wh_tot_err) > thresh)) .and. iterh1 <= maxiter )
       dt_fac = max(imult*iterh1,1)
       dt_new = dtime/dt_fac

       !! restore initial states for a fresh attempt using new sub-timesteps
       if(iterh1 .gt. 0) then
          th_node(:) = th_node_init(:)
       end if
		     
       ! INNER DO-WHILE LOOP
       !! repeats, for dt_frac times, the removal of 1/dt_fac * transpiration (the top boundary flux condition)
       !! stops and returns to outer loop if at any point the water budget doesn't balance, so the timestep can be chopped in half again
       iterh2 = 0
       do while(  (iterh2 .lt. dt_fac) .and. (abs(wh_tot_err) .lt. thresh)  )
          iterh2 = iterh2 + 1

          ! SET DERIVED STATE VARIABLES OVER ALL NODES
          do k = 1, npool_tot
             call psi_from_th(soilstate_inst, ft, porous_media(k), th_node(k), psi_node(k), c)
             call dpsidth_from_th(soilstate_inst, ft, porous_media(k), th_node(k), dpsidth_node(k), c)
             call flc_from_psi(soilstate_inst, ft, porous_media(k), psi_node(k), flc_node(k), c)
             call dflcdpsi_from_psi(soilstate_inst, ft, porous_media(k), psi_node(k), dflcdpsi_node(k), c)
          enddo
		  
          ! SET BOUNDARY PRESSURE DIFFERENCES & CONDUCTANCES
          !! compute water potential differences + conductances and their derivatives wrt water potential
          call boundary_hdiff_and_k(z_node, psi_node, flc_node, dflcdpsi_node, kmax_bound, &
                                    hdiff_bound, k_bound, dkbounddpsi0, dkbounddpsi1)
		  
          ! SET BOUNDARY FLUX TERMS
          !! compute flux terms and their derivatives wrt water content
          q_bound(      1:(npool_tot-1)) =   -1._r8 * k_bound(1:(npool_tot-1)) * hdiff_bound(1:(npool_tot-1))
          dqbounddpsi0( 1:(npool_tot-1)) = ( -1._r8 * k_bound(1:(npool_tot-1)) - dkbounddpsi0(1:(npool_tot-1)) ) * dpsidth_node(1:(npool_tot-1))
          dqbounddpsi1( 1:(npool_tot-1)) = ( -1._r8 * k_bound(1:(npool_tot-1)) - dkbounddpsi1(1:(npool_tot-1)) ) * dpsidth_node(2:npool_tot)
          dqbounddth0(  1:(npool_tot-1)) = ( -1._r8 * k_bound(1:(npool_tot-1)) - dkbounddpsi0(1:(npool_tot-1)) * hdiff_bound(1:(npool_tot-1)) ) &
                                           * dpsidth_node(1:(npool_tot-1))
          dqbounddth1(  1:(npool_tot-1)) = (          k_bound(1:(npool_tot-1)) - dkbounddpsi1(1:(npool_tot-1)) * hdiff_bound(1:(npool_tot-1)) ) &
                                           * dpsidth_node(2: npool_tot   )
          !! zero-flux outer soil shell boundary condition
          q_bound(         npool_tot)    =    0._r8
          dqbounddpsi0(    npool_tot)    =    0._r8
          dqbounddpsi1(    npool_tot)    =    0._r8
          dqbounddth0(     npool_tot)    =    0._r8
          dqbounddth1(     npool_tot)    =    0._r8

          ! STORE BEGINNING WATER BALANCE
          wh_tot_beg = sum( th_node(:)*v_node(:)*denh2o )

          ! SET UP TRIDIAGONAL MATRIX
          !! upper (leaf) layer
          k = 1
          rmx(k)    =  qtop - q_bound(k)
          amx(k)    =  0._r8
          bmx(k)    =  dqbounddth0(k) - 0._r8 - v_node(k)*denh2o/dt_new
          cmx(k)    =  dqbounddth1(k)
          !! intermediate nodes (plant and soil)
          do k=2,(npool_tot-1)
             rmx(k) =  q_bound(k-1) - q_bound(k)
	     amx(k) = -1._r8 * dqbounddth0(k-1)
	     bmx(k) =  dqbounddth0(k) - dqbounddth1(k-1) - v_node(k)*denh2o/dt_new
	     cmx(k) =  dqbounddth1(k)
          enddo
          !! outermost rhizosphere shell
          k = npool_tot
          rmx(k)    =  q_bound(k-1)
          amx(k)    = -1._r8 * dqbounddth0(k-1)
          bmx(k)    = -dqbounddth1(k-1) - v_node(k)*denh2o/dt_new
          cmx(k)    =  0._r8
		      
          ! SOLVE TRIDIAGONAL MATRIX
          call Hydraulics_Tridiagonal(amx, bmx, cmx, rmx, dth_node)
	      
          ! UPDATE WATER BUDGETS
          do k=1,npool_tot
             th_prev     = th_node(k)
	     dth_prev    = dth_node(k)
	     th_node(k)  = max(min(th_node(k)+dth_node(k),ths_node(k)-small_num),thr_node(k)+small_num)
	     dth_node(k) = th_node(k) - th_prev
	     the_node(k) = dth_node(k) - dth_prev
          enddo
			
          ! UPDATE ERROR TERM
          wh_tot_end = sum( th_node(:)*v_node(:)*denh2o )
          wh_tot_err = wh_tot_end - wh_tot_beg + qtop*dt_new
			
       end do ! loop over sub-timesteps

       iterh1 = iterh1 + 1
		     
    end do ! loop to get a timestep divisor that balances water
    
    dth_node(:) = th_node(:) - th_node_init(:)
    
  end subroutine Hydraulics_1DSolve

  !-------------------------------------------------------------------------------!
  subroutine Hydraulics_Tridiagonal(a, b, c, r, u)
    !
    ! !DESCRIPTION: An abbreviated version of biogeophys/TridiagonalMod.F90
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8), intent(in)    :: a(:)           ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: b(:)           ! "b" diagonal column of tridiagonal matrix
    real(r8), intent(in)    :: c(:)           ! "c" right off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: r(:)           ! "r" forcing term of tridiagonal matrix
    real(r8), intent(out)   :: u(:)           ! solution
    !
    ! !LOCAL VARIABLES:
    real(r8) :: bet                           ! temporary
    real(r8) :: gam(npool_tot)                ! temporary
    integer  :: k                             ! index
    !----------------------------------------------------------------------

    bet = b(1)
    do k=1,npool_tot
       if(k == 0) then
          u(k)   = r(k) / bet
       else
          gam(k) = c(k-1) / bet
	  bet    = b(k) - a(k) * gam(k)
	  u(k)   = (r(k) - a(k)*u(k-1)) / bet
       end if
    enddo
  
    do k=npool_tot-1,1,-1
          u(k)   = u(k) - gam(k+1) * u(k+1)
    enddo
    
  end subroutine Hydraulics_Tridiagonal

  !-------------------------------------------------------------------------------!
  subroutine boundary_kmax_rhizosphere(soilstate_inst, c, j, begk, endk, l_aroot, rs1, r_node_shell, kmax_bound)
    !
    ! !DESCRIPTION: 
    !
    ! !USES:
    use clm_varcon           , only : grav, rpi
    use clm_varpar          , only: nlevsoi
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    integer  , intent(in)  :: c               ! current column index
    integer  , intent(in)  :: j               ! current soil layer
    integer  , intent(in)  :: begk, endk      ! shell start and end indices
    real(r8) , intent(in)  :: l_aroot         ! absorbing root length                                           [m]
    real(r8) , intent(in)  :: rs1             ! fine root absorbing radius                                      [m]
    real(r8) , intent(in)  :: r_node_shell(:) ! nodal radius of rhizosphere shells                              [m]
    real(r8) , intent(out) :: kmax_bound(:)   ! maximum conductance at shell boundaries (NOTE:                  [kg s-1 MPa-1] 
                                              ! these are the boundaries closer to atmosphere, not boundaries
					      ! further from atmosphere, as is the case in Hydraulics_BC)
    !
    ! !LOCAL VARIABLES:
    integer    :: k                           ! shell index
    real(r8)   :: hksat_s(nlevsoi)            ! hksat converted to units of 10^6sec which is equiv to           [kg s-1 MPa-1]
    !----------------------------------------------------------------------

    associate(& 
         hksat             =>    soilstate_inst%hksat_col             & ! Input:  [real(r8) (:,:)  ] hydraulic conductivity at saturation (mm H2O /s)
         )

    hksat_s(:)  = hksat(c,:) * 1.e-3_r8 * 1/grav * 1.e6_r8

    do k=begk,endk
       if(k == 1) then
	  kmax_bound(k) = 2*rpi*l_aroot*log(r_node_shell(k)/rs1)*hksat_s(j)
       else
	  kmax_bound(k) = 2*rpi*l_aroot*log(r_node_shell(k)/r_node_shell(k-1))*hksat_s(j)
       end if
    enddo
  
    end associate 

  end subroutine boundary_kmax_rhizosphere
  
  !-------------------------------------------------------------------------------!
  subroutine boundary_hdiff_and_k(z_node, psi_node, flc_node, dflcdpsi_node, kmax_bound, hdiff_bound, k_bound, dkbounddpsi0, dkbounddpsi1)
    !
    ! !DESCRIPTION: 
    !
    ! !USES:
    use clm_varcon         , only : denh2o, grav
    !
    ! !ARGUMENTS
    real(r8)    , intent(in)   :: z_node(:)           ! height of node                                                  [m]
    real(r8)    , intent(in)   :: psi_node(:)         ! water potential in water storage compartments                   [MPa]
    real(r8)    , intent(in)   :: flc_node(:)         ! fractional loss of conductivity at water storage nodes          [-]
    real(r8)    , intent(in)   :: dflcdpsi_node(:)    ! derivative of fractional loss of conductivity wrt psi           [MPa-1]
    real(r8)    , intent(in)   :: kmax_bound(:)       ! lower boundary maximum hydraulic conductance of compartments    [kg s-1 MPa-1]
    real(r8)    , intent(out)  :: hdiff_bound(:)      ! total water potential difference across lower boundary          [MPa-1]
    real(r8)    , intent(out)  :: k_bound(:)          ! lower boundary hydraulic conductance of compartments            [kg s-1 MPa-1]
    real(r8)    , intent(out)  :: dkbounddpsi0(:)     ! derivative of lower boundary conductance wrt psi above          [kg s-1 MPa-2]
    real(r8)    , intent(out)  :: dkbounddpsi1(:)     ! derivative of lower boundary conductance wrt psi below          [kg s-1 MPa-2]
    !
    ! !LOCAL VARIABLES:
    integer  :: k                                     ! shell index
    !----------------------------------------------------------------------

    do k = 1, (npool_tot-1)
       hdiff_bound(k) = 1.e-6_r8*denh2o*grav*(z_node(k) - z_node(k+1)) + (psi_node(k) - psi_node(k+1))
       ! examine direction of water flow; use the upstream node's k for the boundary k. (as suggested by Ethan Coon, LANL)
       if(hdiff_bound(k) < 0._r8) then
	  k_bound(k)   = kmax_bound(k) * flc_node(k+1)  ! water moving towards atmosphere
	  dkbounddpsi0 = 0._r8
          dkbounddpsi1 = k_bound(k) * dflcdpsi_node(k+1)
       else                                           
	  k_bound(k)   = kmax_bound(k) * flc_node(k)    ! water moving towards soil
	  dkbounddpsi0 = k_bound(k) * dflcdpsi_node(k)
          dkbounddpsi1 = 0._r8
       end if
    enddo
  
  end subroutine boundary_hdiff_and_k
  
  !-------------------------------------------------------------------------------!
  subroutine map2d_to_1d_shells(soilstate_inst, waterstate_inst, g, c, rs1, dz_tot, l_aroot_layer, r_out_shell_1D, r_node_shell_1D, ths_shell_1D, &
                                thr_shell_1D, v_shell_1D, psi_node_shell_1D, &
				ksoil_bylayer, ksoil_tot, kmax_bound_bylayershell)
    !
    ! !DESCRIPTION: Converts the two-dimensional system (soil layer x rhizosphere shell) into
    ! an approximate 1-dimensional system (single soil layer with 1 set of rhizosphere shells)
    !
    ! !USES:
    !use EDMainMod            , only : shellGeom
    use clm_varcon           , only : grav
    use clm_varpar           , only : nlevsoi
    use shr_const_mod        , only : SHR_CONST_PI
    use SoilStateType        , only : soilstate_type
    use WaterStateType       , only : waterstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)  :: soilstate_inst
    type(waterstate_type)  , intent(in)     :: waterstate_inst
    integer                , intent(in)     :: g                            ! gridcell index
    integer                , intent(in)     :: c                            ! column index
    real(r8)               , intent(in)     :: rs1                          ! fine root absorbing radius (m)
    real(r8)               , intent(out)    :: dz_tot                       ! total soil thickness (m)
    real(r8)               , intent(in)     :: l_aroot_layer(:)             ! fine root length by soil layer (m)
    real(r8)               , intent(out)    :: r_out_shell_1D(:)            ! outer radii of rhizosphere compartments (m)
    real(r8)               , intent(out)    :: r_node_shell_1D(:)           ! nodal radii of rhizosphere compartments (m)
    real(r8)               , intent(out)    :: ths_shell_1D(:)              ! saturated water content of rhizosphere compartments (m3 m-3)
    real(r8)               , intent(out)    :: thr_shell_1D(:)              ! residual water content of rhizosphere compartments (m3 m-3)
    real(r8)               , intent(out)    :: v_shell_1D(:)                ! shell volume of rhizosphere compartments (m3)
    real(r8)               , intent(out)    :: psi_node_shell_1D(:)         ! soil matric potential of rhizosphere compartments (MPa)
    real(r8)               , intent(out)    :: ksoil_bylayer(:)             ! total rhizosphere conductance (over all shells) in each soil layer (MPa)
    real(r8)               , intent(out)    :: ksoil_tot                    ! total rhizosphere conductance (over all shells and soil layers) (MPa)
    real(r8)               , intent(out)    :: kmax_bound_bylayershell(:,:) ! maximum conductance at shell boundaries (NOTE: at
                                                                            ! boundaries closer to atmosphere, boundaries further
						                            ! from atmosphere, as is the case in Hydraulics_BC)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: area_col                                 ! column area [m2]
    real(r8) :: l_aroot_tot                              ! total fine root length across all soil layers    [m]
    integer  :: j,k,begk,endk                            ! soil layer and shell indices
    ! UNAGGREGATED VARIABLES (2D - by layer and shell)
    real(r8) :: psi_node_bylayershell(nlevsoi,nshell)    !
    real(r8) :: flc_node_bylayershell(nlevsoi,nshell)    !
    real(r8) :: hdiff_bylayershell(nlevsoi,nshell)       ! note location of boundaries as indicated by NOTE above for kmax_bound_bylayershell
    real(r8) :: k_bound_bylayershell(nlevsoi,nshell)     !
    ! AGGREGATED VARIABLES (1D - by shell only)
    real(r8) :: r_out_shell(nshell)                      !
    real(r8) :: r_node_shell(nshell)                     ! 

    real(r8) :: z_node(nshell)                           ! dummy; not needed for this routine but is an output of boundary_hdiff_and_k
    real(r8) :: dflcdpsi_node(nshell)                    ! dummy; not needed for this routine but is an output of boundary_hdiff_and_k
    real(r8) :: dkbounddpsi0(nshell)                     ! dummy; not needed for this routine but is an output of boundary_hdiff_and_k
    real(r8) :: dkbounddpsi1(nshell)                     ! dummy; not needed for this routine but is an output of boundary_hdiff_and_k
    real(r8) :: kmax_bound(   npool_tot)                 ! lower boundary maximum hydraulic conductance of compartments    [kg s-1 MPa-1]
    real(r8) :: k_bound(      npool_tot)                 ! lower boundary hydraulic conductance of compartments            [kg s-1 MPa-1]
    real(r8) :: hdiff_bound(  npool_tot)                 ! total water potential difference across lower boundary          [MPa-1]
    
    real(r8) :: soilresis_bylayer(nlevsoi)               !
    !----------------------------------------------------------------------

    associate(& 
         area_gcell        =>    grc%area                           , & ! Input:  [real(r8) (:)    ] gridcell area (km2)
         wtgcell           =>    col%wtgcell                        , & ! Input:  [real(r8) (:)    ] weight (relative to gridcell)
         dz                =>    col%dz                             , & ! Input:  [real(r8) (:,:)  ] layer thickness (m)                             
         r_out_shell       =>    col%r_out_shell                    , & ! Input:  [real(r8) (:,:,:)] outer radius of rhizosphere compartment (m)
         r_node_shell      =>    col%r_node_shell                   , & ! Input:  [real(r8) (:,:,:)] nodal radius of rhizosphere compartment (m)
         v_shell           =>    col%v_shell                        , & ! Input:  [real(r8) (:,:,:)] volume of rhizosphere compartment (m)
         h2osoi_vol_shell  =>    waterstate_inst%h2osoi_vol_shell   , & ! Input:  [real(r8) (:,:,:)] mass of water in rhizosphere compartment (kg)         h2osoi_vol_shell  =>    waterstate_inst%h2osoi_vol_shell   , & ! Input:  [real(r8) (:,:,:)] volumetric water in rhizosphere compartment (m3/m3)
         hksat             =>    soilstate_inst%hksat_col           , & ! Input:  [real(r8) (:,:)  ] hydraulic conductivity at saturation (mm H2O /s)
         watsat            =>    soilstate_inst%watsat_col          , & ! Input:  [real(r8) (:,:)  ] col volumetric soil water at saturation (porosity)
         watres            =>    soilstate_inst%watres_col          , & ! Input:  [real(r8) (:,:)  ] col volumetric residual soil water
         alpha_VG          =>    soilstate_inst%alpha_VG_col        , & ! Input:  [real(r8) (:,:)  ] col inverse of air-entry pressure (MPa-1)
         n_VG              =>    soilstate_inst%n_VG_col            , & ! Input:  [real(r8) (:,:)  ] col pore-size distribution index
         m_VG              =>    soilstate_inst%m_VG_col            , & ! Input:  [real(r8) (:,:)  ] = 1 - 1/n_VG
         l_VG              =>    soilstate_inst%l_VG_col              & ! Input:  [real(r8) (:,:)  ] col pore tortuosity parameter
         )

    area_col           = area_gcell(g) * wtgcell(c) * 1.e6_r8
    z_node(:)          = 0._r8   ! trivial for the calculations in this subroutine
    dflcdpsi_node(:)   = 0._r8   ! trivial for the calculations in this subroutine

    ! ------------------------------------------------------------------------------------------
    ! MULTIPLE LAYERS X SHELLS: CALCULATIONS TO GET LAYER CONTRIBUTION TO BULK SOIL RHIZOSPHERE PSI PROFILE
    ! ------------------------------------------------------------------------------------------

    ! get the maximum conductance at the boundaries between rhizosphere shells for each soil layer
    ! independent of soil moisture but dependent on fine root length
    do j=1,nlevsoi
       begk = 1
       endk = nshell
       call boundary_kmax_rhizosphere(soilstate_inst, c, j, begk, endk, l_aroot_layer(j), rs1, r_node_shell(c,j,:), kmax_bound_bylayershell(j,:))
    enddo

    ! get the water potentials at the nodes and fractional loss of conductivity at the boundaries of rhizosphere shells for each soil layer
    do j=1, nlevsoi
       do k=1, nshell
          call psi_from_th(soilstate_inst, -9999, 5, h2osoi_vol_shell(c,j,k), psi_node_bylayershell(j,k), c)   ! first argument is a dummy for ft
          call flc_from_psi(soilstate_inst, -9999, 5, psi_node_bylayershell(j,k), flc_node_bylayershell(j,k), c)  ! first argument is a dummy for ft
       enddo
    enddo
    
    ! get the water potential differences and conductances at the boundaries of rhizosphere shells for each layer
    do j=1, nlevsoi
       call boundary_hdiff_and_k(z_node, psi_node_bylayershell(j,:), flc_node_bylayershell(j,:), &
                                 dflcdpsi_node, kmax_bound_bylayershell(j,:), &
				 hdiff_bound, k_bound_bylayershell(j,:), dkbounddpsi0, dkbounddpsi1)
    enddo
    
    ! get each layer's soil conductance (sum across shells) and the total soil conductance (sum across layers)
    soilresis_bylayer(:) = 0._r8
    ksoil_tot            = 0._r8
    do j=1, nlevsoi
       do k=1, nshell
          soilresis_bylayer(j) = soilresis_bylayer(j) + 1._r8/k_bound_bylayershell(j,k)
       enddo
       ksoil_bylayer(j) = 1._r8/soilresis_bylayer(j)
       ksoil_tot        = ksoil_tot + ksoil_bylayer(j)
    enddo
    
    ! get single layer rhizosphere shell-distribution of psi by weighting each layer's shell by that layer's soil conductance
    !      optional: replace total soil conductance with maximum total soil conductance.
    psi_node_shell_1D(:) = 0._r8
    do k=1, nshell
       do j=1, nlevsoi
          psi_node_shell_1D(k) = psi_node_bylayershell(j,k) * ksoil_bylayer(j) / ksoil_tot
       enddo
    enddo
       
    ! ------------------------------------------------------------------------------------------
    ! SINGLE LAYER X SHELLS: GET THE GEOMETRY (SHELL RADII AND VOLUMES)
    ! ------------------------------------------------------------------------------------------

    ! get the single soil layer (bucket) rhizosphere geometry by assuming roots are distributed evenly distributed over the entire profile
    dz_tot = 0._r8
    l_aroot_tot    = 0._r8
    do j=1,nlevsoi
       l_aroot_tot = l_aroot_tot + l_aroot_layer(j)
       dz_tot      = dz_tot + dz(c,j)
    enddo
    call shellGeom(l_aroot_tot, rs1, area_col, dz_tot, r_out_shell_1D, r_node_shell_1D, v_shell_1D)
    
    ! ------------------------------------------------------------------------------------------
    ! SINGLE LAYER X SHELLS: SET THE SATURATED & RESIDUAL WATER CONTENTS
    ! NOTE, FIX: This assumes soil hydraulic properties are constant over layers.  A suitable averaging scheme
    ! should be developed, or this can be avoided altogether when a pressure-based solution gets implemented.
    ! ------------------------------------------------------------------------------------------

    ths_shell_1D(:) = watsat(c,1)
    thr_shell_1D(:) = watres(c,1)
    
    end associate 
  
  end subroutine map2d_to_1d_shells
  
  !-------------------------------------------------------------------------------!
  subroutine map2d_to_1d_aroot(ft, cc_p, kmax_bound_bylayershell, ths_aroot, thr_aroot, vtot_aroot, psi_node_aroot)
    !
    ! !DESCRIPTION: Converts the forked system (troot to multiple aroot by soil layer) into
    ! an approximate 1-dimensional system (troot connecting to single aroot)
    !
    ! !USES:
    use clm_varpar           , only : nlevsoi
    !
    ! !ARGUMENTS
    integer              , intent(in)             :: ft                            ! PFT index
    type(ed_cohort_type) , intent(inout), target  :: cc_p                          ! current cohort pointer
    real(r8)             , intent(in)             :: kmax_bound_bylayershell(:,:)  ! total volume of fine roots (cohort-level)              [m3]
    real(r8)             , intent(out)            :: ths_aroot                     ! saturated water content of fine roots (cohort-level)   [m3]
    real(r8)             , intent(out)            :: thr_aroot                     ! residual water content of fine roots (cohort-level)    [m3]
    real(r8)             , intent(out)            :: vtot_aroot                    ! total volume of fine roots (cohort-level)              [m3]
    real(r8)             , intent(out)            :: psi_node_aroot                ! aggregate water potential of fine roots (cohort-level) [MPa]
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type) , pointer                :: currentCohort                 !
    real(r8)                                      :: kmaxsoil_bylayer(nlevsoi)     !
    real(r8)                                      :: kmaxsoil_tot                  !
    real(r8)                                      :: soilminresis_bylayer(nlevsoi) !
    integer                                       :: j, k                          ! indices
    !----------------------------------------------------------------------

    associate(& 
         thetas   => pftcon%thetas_node        , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content
         resid    => pftcon%resid_node           & ! Input: [real(r8) (:,:) ] P-V curve: residual water fraction
         )

    currentCohort => cc_p

    ! get each layer's maximum soil conductance (sum across shells) and the total maximum soil conductance across all layers
    soilminresis_bylayer(:) = 0._r8
    kmaxsoil_tot            = 0._r8
    do j=1, nlevsoi
       do k=1, nshell
          soilminresis_bylayer(j) = soilminresis_bylayer(j) + 1._r8/kmax_bound_bylayershell(j,k)
       enddo
       kmaxsoil_bylayer(j)  = 1._r8/soilminresis_bylayer(j)
       kmaxsoil_tot         = kmaxsoil_tot + kmaxsoil_bylayer(j)
    enddo
    
    psi_node_aroot = 0._r8
    vtot_aroot     = 0._r8
    do j=1,nlevsoi
       psi_node_aroot = psi_node_aroot + currentCohort%psi_aroot(j) * kmaxsoil_bylayer(j) / kmaxsoil_tot
       vtot_aroot     = vtot_aroot + currentCohort%v_aroot_layer(j)
    enddo
    
    ths_aroot      = pftcon%thetas_node(ft,4)
    thr_aroot      = pftcon%thetas_node(ft,4) * pftcon%resid_node(ft,4)

    end associate

  end subroutine map2d_to_1d_aroot
  
  !-------------------------------------------------------------------------------!
  subroutine flc_from_psi(soilstate_inst, ft, pm, psi_node, flc_node, c)
    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to a fractional loss of conductivity
    !
    ! !USES:
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: psi_node    ! water potential                  [MPa]
    real(r8)         , intent(out)    :: flc_node    ! fractional loss of conductivity  [-]
    integer, optional, intent(in)     :: c           ! column index
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         avuln    => pftcon%avuln_node            , & ! Input: [real(r8) (:,:) ] PLC curve: vulnerability curve shape parameter          [-]
         p50      => pftcon%p50_node              , & ! Input: [real(r8) (:,:) ] PLC curve: water potential at 50% loss of conductivity  [Pa]
         watsat   => soilstate_inst%watsat_col    , & ! Input: [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  [m3 m-3]
         watres   => soilstate_inst%watres_col    , & ! Input: [real(r8) (:,:) ] volumetric residual soil water                  [m3 m-3]
         alpha    => soilstate_inst%alpha_VG_col  , & ! Input: [real(r8) (:,:) ] inverse of air-entry pressure                   [MPa-1]
         n        => soilstate_inst%n_VG_col      , & ! Input: [real(r8) (:,:) ] pore-size distribution index                    [-]
         m        => soilstate_inst%m_VG_col      , & ! Input: [real(r8) (:,:) ] = 1 - 1/n_VG                                    [-]
         l        => soilstate_inst%l_VG_col        & ! Input: [real(r8) (:,:) ] col pore tortuosity parameter                   [-]
         )
    
    if(pm <= 4) then
       flc_node = 1._r8/(1._r8 + (psi_node/p50(ft,pm))**avuln(ft,pm))
    else
       call unsatk_flc_VG(psi_node, watsat(c,1), watres(c,1), alpha(c,1), n(c,1), m(c,1), l(c,1), flc_node)
    end if
	     
    end associate

  end subroutine flc_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine dflcdpsi_from_psi(soilstate_inst, ft, pm, psi_node, dflcdpsi_node, c)
    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to a fractional loss of conductivity
    !
    ! !USES:
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    integer          , intent(in)     :: ft             ! PFT index
    integer          , intent(in)     :: pm             ! porous media index
    real(r8)         , intent(in)     :: psi_node       ! water potential                  [MPa]
    real(r8)         , intent(out)    :: dflcdpsi_node  ! fractional loss of conductivity  [-]
    integer, optional, intent(in)     :: c              ! column index
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         avuln    => pftcon%avuln_node            , & ! Input: [real(r8) (:,:) ] PLC curve: vulnerability curve shape parameter          [-]
         p50      => pftcon%p50_node              , & ! Input: [real(r8) (:,:) ] PLC curve: water potential at 50% loss of conductivity  [Pa]
         watsat   => soilstate_inst%watsat_col    , & ! Input: [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  [m3 m-3]
         watres   => soilstate_inst%watres_col    , & ! Input: [real(r8) (:,:) ] volumetric residual soil water                  [m3 m-3]
         alpha    => soilstate_inst%alpha_VG_col  , & ! Input: [real(r8) (:,:) ] inverse of air-entry pressure                   [MPa-1]
         n        => soilstate_inst%n_VG_col      , & ! Input: [real(r8) (:,:) ] pore-size distribution index                    [-]
         m        => soilstate_inst%m_VG_col      , & ! Input: [real(r8) (:,:) ] = 1 - 1/n_VG                                    [-]
         l        => soilstate_inst%l_VG_col        & ! Input: [real(r8) (:,:) ] col pore tortuosity parameter                   [-]
         )
    
    if(pm <= 4) then
       dflcdpsi_node = -1._r8 * (1._r8 + (psi_node/p50(ft,pm))**avuln(ft,pm))**(-2._r8) * &
                                avuln(ft,pm)/p50(ft,pm)*(psi_node/p50(ft,pm))**(avuln(ft,pm)-1._r8)
    else
       call unsatk_dflcdpsi_VG(psi_node, watsat(c,1), watres(c,1), alpha(c,1), n(c,1), m(c,1), l(c,1), dflcdpsi_node)
    end if
	     
    end associate

  end subroutine dflcdpsi_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine th_from_psi(soilstate_inst, ft, pm, psi_node, th_node, c)
    ! 
    ! !DESCRIPTION: calls necessary routines (plant vs. soil) for converting
    ! plant tissue or soil water potentials to volumetric water contents
    !
    ! !USES:
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: psi_node    ! water potential   [MPa]
    real(r8)         , intent(out)    :: th_node     ! water content     [m3 m-3]
    integer, optional, intent(in)     :: c           ! column index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lower                ! lower bound of initial estimate         [m3 m-3]
    real(r8) :: upper                ! upper bound of initial estimate         [m3 m-3]
    real(r8) :: xtol                 ! error tolerance for x-variable          [m3 m-3]
    real(r8) :: ytol                 ! error tolerance for y-variable          [MPa]
    real(r8) :: satfrac              ! soil saturation fraction                [0-1]
    !----------------------------------------------------------------------
  
    associate(& 
         thetas   => pftcon%thetas_node           , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content
         resid    => pftcon%resid_node            , & ! Input: [real(r8) (:,:) ] P-V curve: residual water fraction
         corrint  => pftcon%corrint_node          , & ! Input: [real(r8) (:,:) ] P-V curve: correction for nonzero psi0
         watsat   => soilstate_inst%watsat_col    , & ! Input: [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  [m3 m-3]
         watres   => soilstate_inst%watres_col    , & ! Input: [real(r8) (:,:) ] volumetric residual soil water                  [m3 m-3]
         alpha    => soilstate_inst%alpha_VG_col  , & ! Input: [real(r8) (:,:) ] inverse of air-entry pressure                   [MPa-1]
         n        => soilstate_inst%n_VG_col      , & ! Input: [real(r8) (:,:) ] pore-size distribution index                    [-]
         m        => soilstate_inst%m_VG_col      , & ! Input: [real(r8) (:,:) ] = 1 - 1/n_VG                                    [-]
         l        => soilstate_inst%l_VG_col        & ! Input: [real(r8) (:,:) ] col pore tortuosity parameter                   [-]
         )
   
    if(pm <= 4) then
       lower  = thetas(ft,pm)*(resid(ft,pm) + 0.0001_r8)/corrint(ft,pm)
       upper  = thetas(ft,pm)
       xtol   = 1.e-5_r8
       ytol   = 1.e-6_r8
       call bisect_pv(ft, pm, lower, upper, xtol, ytol, psi_node, th_node)
    else
       call swc_satfrac_VG(psi_node, alpha(c,1), n(c,1), m(c,1), l(c,1), satfrac)
       th_node = watsat(c,1) - satfrac*(watsat(c,1) - watres(c,1))
    end if
	     
    end associate

  end subroutine th_from_psi
  
  !-------------------------------------------------------------------------------!
  subroutine bisect_pv(ft, pm, lower, upper, xtol, ytol, psi_node, th_node)
    ! 
    ! !DESCRIPTION: Bisection routine for getting the inverse of the plant PV curve.
    !  An analytical solution is not possible because quadratic smoothing functions
    !  are used to remove discontinuities in the PV curve.
    !
    ! !USES:
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(inout)  :: lower       ! lower bound of estimate           [m3 m-3]
    real(r8)      , intent(inout)  :: upper       ! upper bound of estimate           [m3 m-3]
    real(r8)      , intent(in)     :: xtol        ! error tolerance for x-variable    [m3 m-3]
    real(r8)      , intent(in)     :: ytol        ! error tolerance for y-variable    [MPa]
    real(r8)      , intent(in)     :: psi_node    ! water potential                   [MPa]
    real(r8)      , intent(out)    :: th_node     ! water content                     [m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x_new                  ! new estimate for x in bisection routine
    real(r8) :: y_lo                   ! corresponding y value at lower
    real(r8) :: f_lo                   ! y difference between lower bound guess and target y
    real(r8) :: y_hi                   ! corresponding y value at upper
    real(r8) :: f_hi                   ! y difference between upper bound guess and target y
    real(r8) :: y_new                  ! corresponding y value at x.new
    real(r8) :: f_new                  ! y difference between new y guess at x.new and target y
    real(r8) :: chg                    ! difference between x upper and lower bounds (approach 0 in bisection)
    type(soilstate_type) :: soilstate_inst ! dummy var
    !----------------------------------------------------------------------
  
    call psi_from_th(soilstate_inst, ft, pm, lower, y_lo)
    call psi_from_th(soilstate_inst, ft, pm, lower, y_hi)
    f_lo  = y_lo - psi_node
    f_hi  = y_hi - psi_node
    chg   = upper - lower
    do while(abs(chg) .gt. xtol)
       x_new = 0.5_r8*(lower + upper)
       call psi_from_th(soilstate_inst, ft, pm, x_new, y_new)
       f_new = y_new - psi_node
       if(abs(f_new) .le. ytol) then
          EXIT
       end if
       if((f_lo * f_new) .lt. 0._r8) upper = x_new
       if((f_hi * f_new) .lt. 0._r8) lower = x_new
       chg = upper - lower
    end do
    
    th_node = x_new
	     
  end subroutine bisect_pv
  
  !-------------------------------------------------------------------------------!
  subroutine psi_from_th(soilstate_inst, ft, pm, th_node, psi_node, c)
    ! 
    ! !DESCRIPTION: evaluates the plant PV curve (returns water potential, psi)
    ! at a given water content (th)
    !
    ! !USES:
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: th_node     ! water content     [m3 m-3]
    real(r8)         , intent(out)    :: psi_node    ! water potential   [MPa]
    integer, optional, intent(in)     :: c           ! column index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: satfrac                  ! saturation fraction [0-1]
    !----------------------------------------------------------------------
  
    associate(& 
         corrint  => pftcon%corrint_node          , & ! Input: [real(r8) (:,:) ] P-V curve: correction for nonzero psi0 (pft x porous_media)
         watsat   => soilstate_inst%watsat_col    , & ! Input: [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  [m3 m-3]
         watres   => soilstate_inst%watres_col    , & ! Input: [real(r8) (:,:) ] volumetric residual soil water                  [m3 m-3]
         alpha    => soilstate_inst%alpha_VG_col  , & ! Input: [real(r8) (:,:) ] inverse of air-entry pressure                   [MPa-1]
         n        => soilstate_inst%n_VG_col      , & ! Input: [real(r8) (:,:) ] pore-size distribution index                    [-]
         m        => soilstate_inst%m_VG_col      , & ! Input: [real(r8) (:,:) ] = 1 - 1/n_VG                                    [-]
         l        => soilstate_inst%l_VG_col        & ! Input: [real(r8) (:,:) ] col pore tortuosity parameter                   [-]
         )
   
    if(pm <= 4) then       ! plant

       call tq2(ft, pm, th_node*corrint(ft,pm), psi_node)

    else if(pm == 5) then  ! soil

!! NOTE. FIX: The below sidesteps the problem of averaging potentially variable soil hydraulic properties with depth
!!        and simply assigns the bulk soil (bucket) approximation of hydraulic properties as equal to the top soil layer.
       satfrac = (watsat(c,1) - th_node)/(watsat(c,1) - watres(c,1))
       call swc_psi_VG(satfrac, alpha(c,1), n(c,1), m(c,1), l(c,1), psi_node)

    end if
	     
    end associate

  end subroutine psi_from_th
  
  !-------------------------------------------------------------------------------!
  subroutine dpsidth_from_th(soilstate_inst, ft, pm, th_node, y, c)
    ! 
    ! !DESCRIPTION: evaluates the plant PV curve (returns water potential, psi)
    ! at a given water content (th)
    !
    ! !USES:
    use SoilStateType      , only : soilstate_type
    !
    ! !ARGUMENTS
    type(soilstate_type)   , intent(inout)         :: soilstate_inst
    integer          , intent(in)     :: ft          ! PFT index
    integer          , intent(in)     :: pm          ! porous media index
    real(r8)         , intent(in)     :: th_node     ! water content                            [m3 m-3]
    real(r8)         , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    integer, optional, intent(in)     :: c           ! column index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: satfrac                  ! saturation fraction [0-1]
    !----------------------------------------------------------------------
  
    associate(& 
         corrint  => pftcon%corrint_node          , & ! Input: [real(r8) (:,:) ] P-V curve: correction for nonzero psi0
         watsat   => soilstate_inst%watsat_col    , & ! Input: [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  [m3 m-3]
         watres   => soilstate_inst%watres_col    , & ! Input: [real(r8) (:,:) ] volumetric residual soil water                  [m3 m-3]
         alpha    => soilstate_inst%alpha_VG_col  , & ! Input: [real(r8) (:,:) ] inverse of air-entry pressure                   [MPa-1]
         n        => soilstate_inst%n_VG_col      , & ! Input: [real(r8) (:,:) ] pore-size distribution index                    [-]
         m        => soilstate_inst%m_VG_col      , & ! Input: [real(r8) (:,:) ] = 1 - 1/n_VG                                    [-]
         l        => soilstate_inst%l_VG_col        & ! Input: [real(r8) (:,:) ] col pore tortuosity parameter                   [-]
         )
   
    if(pm <= 4) then       ! plant
       call dtq2dth(ft, pm, th_node*corrint(ft,pm), y)
    else if(pm == 5) then  ! soil
       satfrac = (watsat(c,1) - th_node)/(watsat(c,1) - watres(c,1))
       call swc_dpsidth_VG(satfrac, watsat(c,1), watres(c,1), alpha(c,1), n(c,1), m(c,1), l(c,1), y)
    end if
	     
    end associate

  end subroutine dpsidth_from_th
  
  !-------------------------------------------------------------------------------!
  subroutine tq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: smoothing function for elastic-to-cavitation region of the
    !  plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq2                  ! returned y (psi) value from bq2()
    real(r8) :: y_cq2                  ! returned y (psi) value from cq2()
    real(r8) :: beta2=0.99_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    call bq2(ft, pm, x, y_bq2)
    call cq2(ft, pm, x, y_cq2)
    y = (-y_bq2 + sqrt(y_bq2*y_bq2 - 4._r8*beta2*y_cq2))/(2*beta2)

  end subroutine tq2
  
  !-------------------------------------------------------------------------------!
  subroutine dtq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: smoothing function for elastic-to-cavitation region of the
    !  plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq2                  ! returned y (psi) value from bq2()
    real(r8) :: y_cq2                  ! returned y (psi) value from cq2()
    real(r8) :: dydth_bq2              ! returned derivative from dbq2dth()
    real(r8) :: dydth_cq2              ! returned derivative from dcq2dth()
    real(r8) :: beta2=0.99_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    call bq2(ft, pm, x, y_bq2)
    call cq2(ft, pm, x, y_cq2)
    call dbq2dth(ft, pm, x, dydth_bq2)
    call dcq2dth(ft, pm, x, dydth_cq2)
    y = 1._r8/(2._r8*beta2)*(-dydth_bq2 + 0.5_r8*((y_bq2*y_bq2 - 4._r8*beta2*y_cq2)**(-0.5_r8)) * &
                                                 (2._r8*y_bq2*dydth_bq2 - 4._r8*beta2*dydth_cq2))

  end subroutine dtq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine bq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    y = -1._r8*(y_tq1 + y_cavitation)
	     
  end subroutine bq2
  
  !-------------------------------------------------------------------------------!
  subroutine dbq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dydth_tq1              ! returned derivative from dtq1dth()
    real(r8) :: dcavdth                ! returned derivative from dcavitationdth()
    !----------------------------------------------------------------------
  
    call dtq1dth(ft, pm, x, dydth_tq1)
    call dcavitationPVdth(ft, pm, x, dcavdth)
    y = -1._r8*(dydth_tq1 + dcavdth)

  end subroutine dbq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine cq2(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for elastic-to-cavitation region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    y = y_tq1*y_cavitation
	     
  end subroutine cq2
  
  !-------------------------------------------------------------------------------!
  subroutine dcq2dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cq2() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_tq1                  ! returned y (psi) value from tq1()
    real(r8) :: y_cavitation           ! returned y (psi) value from cavitationPV()
    real(r8) :: dydth_tq1              ! returned derivative from dtq1dth()
    real(r8) :: dcavdth                ! returned derivative from dcavitationdth()
    !----------------------------------------------------------------------
  
    call tq1(ft, pm, x, y_tq1)
    call cavitationPV(ft, pm, x, y_cavitation)
    call dtq1dth(ft, pm, x, dydth_tq1)
    call dcavitationPVdth(ft, pm, x, dcavdth)
    y = y_tq1*dcavdth + dydth_tq1*y_cavitation
	     
  end subroutine dcq2dth
  
  !-------------------------------------------------------------------------------!
  subroutine tq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: either calls the elastic region of the PV curve (leaves) or
    ! does a smoothing function for capillary-to-elastic region of the plant PV
    ! curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq1                  ! returned y (psi) value from bq1()
    real(r8) :: y_cq1                  ! returned y (psi) value from cq1()
    real(r8) :: y_elastic              ! returned y (psi) value from elasticPV()
    real(r8) :: beta1=0.80_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    if(pm == 1) then            ! leaves have no capillary region in their PV curves
       call elasticPV(ft, pm, x, y_elastic)
       y = y_elastic
    else if(pm <= 4) then       ! sapwood has a capillary region
       call bq1(ft, pm, x, y_bq1)
       call cq1(ft, pm, x, y_cq1)
       y = (-y_bq1 - sqrt(y_bq1*y_bq1 - 4._r8*beta1*y_cq1))/(2*beta1)
    end if !porous media

  end subroutine tq1
  
  !-------------------------------------------------------------------------------!
  subroutine dtq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of tq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_bq1                  ! returned y (psi) value from bq1()
    real(r8) :: y_cq1                  ! returned y (psi) value from cq1()
    real(r8) :: dydth_bq1              ! returned derivative from dbq1dth()
    real(r8) :: dydth_cq1              ! returned derivative from dcq1dth()
    real(r8) :: delasticdth            ! returned derivative from delasticPVdth()
    real(r8) :: beta1=0.80_r8          ! smoothing factor
    !----------------------------------------------------------------------
  
    if(pm == 1) then            ! leaves have no capillary region in their PV curves
       call delasticPVdth(ft, pm, x, delasticdth)
       y = delasticdth
    else if(pm <= 4) then       ! sapwood has a capillary region
       call bq1(ft, pm, x, y_bq1)
       call cq1(ft, pm, x, y_cq1)
       call dbq1dth(ft, pm, x, dydth_bq1)
       call dcq1dth(ft, pm, x, dydth_cq1)
       y = 1._r8/(2._r8*beta1)*(-dydth_bq1 - 0.5_r8*((y_bq1*y_bq1 - 4._r8*beta1*y_cq1)**(-0.5_r8)) * &
                                                    (2._r8*y_bq1*dydth_bq1 - 4._r8*beta1*dydth_cq1))
    end if

  end subroutine dtq1dth
  
  !-------------------------------------------------------------------------------!
  subroutine bq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for capillary-to-elastic region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    y = -1._r8*(y_capillary + y_elastic)
	     
  end subroutine bq1
  
  !-------------------------------------------------------------------------------!
  subroutine dbq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of bq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dcapdth               ! returned derivative from dcapillaryPVdth()
    real(r8) :: delasticdth           ! returned derivative from delasticPVdth()
    !----------------------------------------------------------------------
  
    call dcapillaryPVdth(ft, pm, x, dcapdth)
    call delasticPVdth(ft, pm, x, delasticdth)
    y = -1._r8*(delasticdth + dcapdth)
	     
  end subroutine dbq1dth
  
  !-------------------------------------------------------------------------------!
  subroutine cq1(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: component smoothing function for capillary-to-elastic region
    ! of the plant PV curve where a discontinuity exists
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    y = y_capillary*y_elastic
	     
  end subroutine cq1
  
  !-------------------------------------------------------------------------------!
  subroutine dcq1dth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cq1() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_capillary           ! returned y (psi) value from capillaryPV()
    real(r8) :: y_elastic             ! returned y (psi) value from elasticPV()
    real(r8) :: dcapdth               ! returned derivative from dcapillaryPVdth()
    real(r8) :: delasticdth           ! returned derivative from delasticPVdth()
    !----------------------------------------------------------------------
  
    call capillaryPV(ft, pm, x, y_capillary)
    call elasticPV(ft, pm, x, y_elastic)
    call dcapillaryPVdth(ft, pm, x, dcapdth)
    call delasticPVdth(ft, pm, x, delasticdth)
    y = y_elastic*dcapdth + delasticdth*y_capillary
	     
  end subroutine dcq1dth
  
  !-------------------------------------------------------------------------------!
  subroutine cavitationPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the elastic region of the plant PV
    ! curve as the sum of both solute and elastic components.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_solute           ! returned y (psi) value from solutepsi()
    !----------------------------------------------------------------------
  
    call solutepsi(ft, pm, x, y_solute)
    y = y_solute
	     
  end subroutine cavitationPV
  
  !-------------------------------------------------------------------------------!
  subroutine dcavitationPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of cavitationPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dsoldth           ! returned derivative from dsolutepsidth()
    !----------------------------------------------------------------------
  
    call dsolutepsidth(ft, pm, x, dsoldth)
    y = dsoldth
	     
  end subroutine dcavitationPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine elasticPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the elastic region of the plant PV
    ! curve as the sum of both solute and elastic components.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: y_solute           ! returned y (psi) value from solutepsi()
    real(r8) :: y_pressure         ! returned y (psi) value from pressurepsi()
    !----------------------------------------------------------------------
  
    call solutepsi(ft, pm, x, y_solute)
    call pressurepsi(ft, pm, x, y_pressure)
    y = y_solute + y_pressure
	     
  end subroutine elasticPV
  
  !-------------------------------------------------------------------------------!
  subroutine delasticPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of elasticPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dsoldth           ! returned derivative from dsolutepsidth()
    real(r8) :: dpressdth         ! returned derivative from dpressurepsidth()
    !----------------------------------------------------------------------
  
    call dsolutepsidth(ft, pm, x, dsoldth)
    call dpressurepsidth(ft, pm, x, dpressdth)
    y = dsoldth + dpressdth
	     
  end subroutine delasticPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine solutepsi(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes solute water potential (negative) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         pinot   => pftcon%pinot_node     , & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => pftcon%thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 rwcft   => pftcon%rwcft_node     , & ! Input: [real(r8) (:,:) ] P-V curve: total RWC @ which elastic drainage begins     [-]
	 resid   => pftcon%resid_node       & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         )
    
    y = pinot(ft,pm)*thetas(ft,pm)*(rwcft(ft,pm) - resid(ft,pm)) / &
        (x - thetas(ft,pm)*resid(ft,pm))
	     
    end associate

  end subroutine solutepsi
  
  !-------------------------------------------------------------------------------!
  subroutine dsolutepsidth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of solutepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         pinot   => pftcon%pinot_node     , & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => pftcon%thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 rwcft   => pftcon%rwcft_node     , & ! Input: [real(r8) (:,:) ] P-V curve: total RWC @ which elastic drainage begins     [-]
	 resid   => pftcon%resid_node       & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         )
    
    y = -1._r8*thetas(ft,pm)*pinot(ft,pm)*(rwcft(ft,pm) - resid(ft,pm)) / &
        ((x - thetas(ft,pm)*resid(ft,pm))**2._r8)
	     
    end associate

  end subroutine dsolutepsidth
  
  !-------------------------------------------------------------------------------!
  subroutine pressurepsi(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes pressure water potential (positive) as a function of
    !  water content for the plant PV curve.
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         pinot   => pftcon%pinot_node     , & ! Input: [real(r8) (:,:) ] P-V curve: osmotic potential at full turgor              [MPa]
         thetas  => pftcon%thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 rwcft   => pftcon%rwcft_node     , & ! Input: [real(r8) (:,:) ] P-V curve: total RWC @ which elastic drainage begins     [-]
	 resid   => pftcon%resid_node     , & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         epsil   => pftcon%epsil_node       & ! Input: [real(r8) (:,:) ] P-V curve: bulk elastic modulus                          [MPa]
         )
   
    y = epsil(ft,pm) * (x - thetas(ft,pm)*rwcft(ft,pm)) / &
                      (thetas(ft,pm)*(rwcft(ft,pm)-resid(ft,pm))) - pinot(ft,pm)
	     
    end associate

  end subroutine pressurepsi
  
  !-------------------------------------------------------------------------------!
  subroutine dpressurepsidth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of pressurepsi() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
  
    associate(& 
         thetas  => pftcon%thetas_node    , & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
	 rwcft   => pftcon%rwcft_node     , & ! Input: [real(r8) (:,:) ] P-V curve: total RWC @ which elastic drainage begins     [-]
	 resid   => pftcon%resid_node     , & ! Input: [real(r8) (:,:) ] P-V curve: residual fraction                             [-]
         epsil   => pftcon%epsil_node       & ! Input: [real(r8) (:,:) ] P-V curve: bulk elastic modulus                          [MPa]
         )
   
    y = epsil(ft,pm)/(thetas(ft,pm)*(rwcft(ft,pm) - resid(ft,pm)))
	     
    end associate

  end subroutine dpressurepsidth
  
  !-------------------------------------------------------------------------------!
  subroutine capillaryPV(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: computes water potential in the capillary region of the plant
    !  PV curve (sapwood only)
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content     [m3 m-3]
    real(r8)      , intent(out)    :: y           ! water potential   [MPa]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    associate(& 
         intercept => pftcon%intercept_node     , & ! Input: [real(r8) (:,:) ] P-V curve: intercept of capillary region of curve (sapwood only)
         slp       => pftcon%slp_node             & ! Input: [real(r8) (:,:) ] P-V curve: slope of capillary region of curve (sapwood only)
         )
   
    y = intercept(ft,pm) + slp(ft,pm)*x
	     
    end associate

  end subroutine capillaryPV
  
  !-------------------------------------------------------------------------------!
  subroutine dcapillaryPVdth(ft, pm, x, y)
    ! 
    ! !DESCRIPTION: returns derivative of capillaryPV() wrt theta
    !
    ! !USES:
    !
    ! !ARGUMENTS
    integer       , intent(in)     :: ft          ! PFT index
    integer       , intent(in)     :: pm          ! porous media index
    real(r8)      , intent(in)     :: x           ! water content                            [m3 m-3]
    real(r8)      , intent(out)    :: y           ! derivative of water potential wrt theta  [MPa m3 m-3]
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    associate(& 
         slp       => pftcon%slp_node           , & ! Input: [real(r8) (:,:) ] P-V curve: slope of capillary region of curve (sapwood only)
         thetas    => pftcon%thetas_node          & ! Input: [real(r8) (:,:) ] P-V curve: saturated volumetric water content for node   [m3 m-3]
         )
   
    y = slp(ft,pm)/thetas(ft,pm)
	     
    end associate

  end subroutine dcapillaryPVdth
  
  !-------------------------------------------------------------------------------!
  subroutine swc_psi_VG(satfrac, alpha, n, m, l, psi)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns water potential given saturation fraction and shape parameters
  !
  !USES
  use shr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: satfrac  !saturation fraction            [0-1]
  real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
  real(r8), intent(in)            :: n        !pore-size distribution index   [-]
  real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
  real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
  real(r8), intent(out)           :: psi      !soil matric potential          [MPa]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  psi = -1._r8/alpha*(satfrac**(-1._r8/m)-1._r8)**(1._r8/n)

  end subroutine swc_psi_VG

  !-----------------------------------------------------------------------
  subroutine swc_satfrac_VG(psi, alpha, n, m, l, satfrac)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns saturation fraction given water potential and shape parameters
  !
  !USES
  use shr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! !ARGUMENTS:
  real(r8), intent(in)            :: psi      !soil matric potential          [MPa]
  real(r8), intent(in)            :: alpha    !inverse of air-entry pressure  [MPa-1]
  real(r8), intent(in)            :: n        !pore-size distribution index   [-]
  real(r8), intent(in)            :: m        != 1 - 1/n_VG                   [-]
  real(r8), intent(in)            :: l        !pore tortuosity parameter      [-]
  real(r8), intent(out)           :: satfrac  !soil saturation fraction       [0-1]
  !
  ! !LOCAL VARIABLES:
  !------------------------------------------------------------------------------

  satfrac = (1._r8/(1._r8 + (alpha*abs(psi))**n))**m

  end subroutine swc_satfrac_VG

  !-----------------------------------------------------------------------
  subroutine swc_dpsidth_VG(satfrac, watsat, watres, alpha, n, m, l, dpsidth)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given saturation fraction and shape parameters
  !
  !USES
  use shr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: satfrac  !saturation fraction                            [0-1]
  real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)  :: watres   !volumetric residual soil water                 [m3 m-3]
  real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
  real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
  real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
  real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
  real(r8), intent(out) :: dpsidth  !derivative of psi wrt theta                    [MPa/m3m-3]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: temp0    !temporary
  real(r8)              :: temp1    !temporary
  real(r8)              :: temp2    !temporary
  real(r8)              :: temp3    !temporary
  !------------------------------------------------------------------------------

  temp0   = 1._r8/(m*n*alpha*(watsat-watres))
  temp1   = satfrac**(-1._r8/m) - 1._r8
  temp2   = temp1**(1._r8/n - 1._r8)
  temp3   = satfrac**(-1._r8/m - 1._r8)
  dpsidth = temp0*temp2*temp3

  end subroutine swc_dpsidth_VG

  !-----------------------------------------------------------------------
  subroutine unsatk_flc_VG(psi, watsat, watres, alpha, n, m, l, flc)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given saturation fraction and shape parameters
  !
  !USES
  use shr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
  real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)  :: watres   !volumetric residual soil water                 [m3 m-3]
  real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
  real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
  real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
  real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
  real(r8), intent(out) :: flc      !k/ksat ('fractional loss of conductivity')     [-]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: temp          !temporary
  real(r8)              :: fac1a         !temporary
  real(r8)              :: fac1b         !temporary
  real(r8)              :: fac1          !temporary
  real(r8)              :: fac2          !temporary
  !------------------------------------------------------------------------------

  temp       = ( alpha*abs(psi)      ) ** (n)
  fac1a      = ( alpha*abs(psi)      ) ** (n-1._r8)
  fac1b      = ( 1._r8 + temp        ) ** (-1._r8*m)
  fac1       = ( 1._r8 - fac1a*fac1b ) ** (2._r8)
  fac2       = ( 1._r8 + temp        ) ** (-0.5_r8*m)
  
  flc        =   fac1 * fac2

  end subroutine unsatk_flc_VG

  !-----------------------------------------------------------------------
  subroutine unsatk_dflcdpsi_VG(psi, watsat, watres, alpha, n, m, l, dflcdpsi)
  !
  ! DESCRIPTION
  ! van Genuchten (1980) soil water characteristic (retention) curve
  ! returns derivative of water water potential with respect to water content
  ! given saturation fraction and shape parameters
  !
  !USES
  use shr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! !ARGUMENTS:
  real(r8), intent(in)  :: psi      !soil matric potential                          [MPa]
  real(r8), intent(in)  :: watsat   !volumetric soil water at saturation (porosity) [m3 m-3]
  real(r8), intent(in)  :: watres   !volumetric residual soil water                 [m3 m-3]
  real(r8), intent(in)  :: alpha    !inverse of air-entry pressure                  [MPa-1]
  real(r8), intent(in)  :: n        !pore-size distribution index                   [-]
  real(r8), intent(in)  :: m        != 1 - 1/n_VG                                   [-]
  real(r8), intent(in)  :: l        !pore tortuosity parameter                      [-]
  real(r8), intent(out) :: dflcdpsi !derivative of k/ksat (flc) wrt psi             [MPa-1]
  !
  ! !LOCAL VARIABLES:
  real(r8)              :: temp          !temporary
  real(r8)              :: fac1a         !temporary
  real(r8)              :: fac1b         !temporary
  real(r8)              :: fac1          !temporary
  real(r8)              :: fac2          !temporary
  real(r8)              :: dtemp         !temporary
  real(r8)              :: dfac1adpsi    !temporary
  real(r8)              :: dfac1bdpsi    !temporary
  real(r8)              :: dfac1dpsi     !temporary
  real(r8)              :: dfac2dpsi     !temporary
  !------------------------------------------------------------------------------

  temp       = ( alpha*abs(psi)      ) ** (n)
  fac1a      = ( alpha*abs(psi)      ) ** (n-1._r8)
  fac1b      = ( 1._r8 + temp        ) ** (-1._r8*m)
  fac1       = ( 1._r8 - fac1a*fac1b ) ** (2._r8)
  fac2       = ( 1._r8 + temp        ) ** (-0.5_r8*m)
  
  dtemp      =   n * alpha * ( alpha*abs(psi) ) ** (n-1._r8)
  dfac1adpsi = ( n-1._r8 ) * alpha * ( alpha*abs(psi) ) ** (n-2._r8)
  dfac1bdpsi = ( -1._r8 ) * m * dtemp * ( 1._r8 + temp ) ** (-1._r8*m - 1._r8)
  dfac1dpsi  = ( 2._r8 ) * ( 1._r8 - fac1a*fac1b ) * ( -1._r8*dfac1bdpsi*fac1a - dfac1adpsi*fac1b )
  dfac2dpsi  = ( -0.5_r8 ) * m * dtemp * (1._r8 + temp)**(-0.5_r8*m-1._r8)
  
  dflcdpsi   = ( -1._r8 ) * ( dfac2dpsi*fac1 + dfac1dpsi*fac2 )    ! BOC... mult by -1 because unsatk eqn is based on abs(psi)

  end subroutine unsatk_dflcdpsi_VG

  !-------------------------------------------------------------------------------!
  subroutine shellGeom(l_aroot, rs1, area, dz, r_out_shell, r_node_shell, v_shell)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil surrounding fine roots remains
    ! the same.  
    !
    ! !USES:
    !use EDPlantHydraulicsMod , only : nshell                             !BOC...# of rhizosphere 'shells'
    use shr_const_mod        , only : SHR_CONST_PI
    !
    ! !ARGUMENTS:
    real(r8)     , intent(in)             :: l_aroot
    real(r8)     , intent(in)             :: rs1
    real(r8)     , intent(in)             :: area
    real(r8)     , intent(in)             :: dz
    real(r8)     , intent(out)            :: r_out_shell(:)
    real(r8)     , intent(out)            :: r_node_shell(:)
    real(r8)     , intent(out)            :: v_shell(:)                 ! volume of a single rhizosphere shell
    !
    ! !LOCAL VARIABLES:
    real(r8)                       :: pi                                !
    integer                        :: k                                 ! rhizosphere shell indicies
    !-----------------------------------------------------------------------
    
    pi                      = SHR_CONST_PI

    ! update outer radii of column-level rhizosphere shells (same across patches and cohorts)
    r_out_shell(nshell) = (pi*l_aroot/(area*dz))**(-0.5_r8)                  ! eqn(8) S98
    do k = 1,nshell-1
       r_out_shell(k)   = rs1*(r_out_shell(nshell)/rs1)**(k/nshell)          ! eqn(7) S98
    enddo

    ! set nodal (midpoint) radii of these shells
    do k = 1,nshell
       if(k == 1) then
          ! BOC...not doing this as it requires PFT-specific fine root thickness, but this is at column level
          ! r_node_shell(k) = 0.5_r8*(rs1 + r_out_shell(k))
          r_node_shell(k) = 0.5_r8*(r_out_shell(k))
       else
          r_node_shell(k) = 0.5_r8*(r_out_shell(k-1) + r_out_shell(k))
       end if
    enddo

    ! update volumes
    do k = 1,nshell
       if(k == 1) then
	  ! BOC...not doing this as it requires PFT-specific fine root thickness but this is at column level
          ! v_shell(k)   = pi*(r_out_shell(k)**2._r8 - rs1**2._r8)*dz
          v_shell(k)     = pi*(r_out_shell(k)**2._r8)*dz
       else
          v_shell(k)     = pi*(r_out_shell(k)**2._r8 - r_out_shell(k-1)**2._r8)*dz
       end if
    enddo

  end subroutine shellGeom
  
  !-----------------------------------------------------------------------

end module EDPlantHydraulicsMod
