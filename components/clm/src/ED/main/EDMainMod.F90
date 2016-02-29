module EDMainMod

  ! ===========================================================================
  ! Main ED module.    
  ! ============================================================================

  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc
  use decompMod            , only : bounds_type
  use clm_varctl           , only : iulog
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type  
  use TemperatureType      , only : temperature_type
  use WaterStateType       , only : waterstate_type
  use EDCohortDynamicsMod  , only : allocate_live_biomass, terminate_cohorts, fuse_cohorts, sort_cohorts, count_cohorts
  use EDPatchDynamicsMod   , only : disturbance_rates, fuse_patches, spawn_patches, terminate_patches
  use EDPhysiologyMod      , only : canopy_derivs, non_canopy_derivs, phenology, recruitment, trim_canopy
  use SFMainMod            , only : fire_model
  use EDtypesMod           , only : ncwd, n_sub, numpft_ed, udata
  use EDtypesMod           , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDPhenologyType      , only : ed_phenology_type
  use EDCLMLinkMod         , only : ed_clm_type

  implicit none
  private

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: ed_driver
  public  :: ed_update_site
  integer             , parameter :: nshell      = 11                     ! number of concentric soil cylinders surrounding absorbing root
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ed_ecosystem_dynamics
  private :: ed_integrate_state_variables
  private :: ed_total_balance_check
  private :: hydraulics_TreeHydGeom               ! BC added for hydraulics
  private :: hydraulics_RhizHydGeom               ! BC added for hydraulics
  
  logical :: DEBUG  = .false.
  !
  ! 10/30/09: Created by Rosie Fisher
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ed_driver( bounds, ed_allsites_inst, ed_clm_inst, ed_phenology_inst, &
       atm2lnd_inst, soilstate_inst, temperature_inst, waterstate_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Main ed model routine containing gridcell loop   
    !
    ! !USES:
    use clm_time_manager     , only : get_days_per_year, get_curr_date
    use clm_time_manager     , only : get_ref_date, timemgr_datediff 
    use CanopySTateType      , only : canopystate_type
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)            :: bounds           
    type(ed_site_type)      , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    type(ed_clm_type)       , intent(inout)         :: ed_clm_inst
    type(ed_phenology_type) , intent(inout)         :: ed_phenology_inst
    type(atm2lnd_type)      , intent(in)            :: atm2lnd_inst
    type(soilstate_type)    , intent(in)            :: soilstate_inst
    type(temperature_type)  , intent(in)            :: temperature_inst
    type(waterstate_type)   , intent(inout)         :: waterstate_inst
    type(canopystate_type)  , intent(inout)         :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_site_type), pointer :: currentSite
    real(r8) :: dayDiff                  ! day of run
    integer  :: dayDiffInt               ! integer of day of run
    integer  :: g                        ! gridcell  
    integer  :: yr                       ! year (0, ...)
    integer  :: mon                      ! month (1, ..., 12)
    integer  :: day                      ! day of month (1, ..., 31)
    integer  :: sec                      ! seconds of the day
    integer  :: ncdate                   ! current date
    integer  :: nbdate                   ! base date (reference date)
    !-----------------------------------------------------------------------

    call ed_clm_inst%SetValues( bounds, 0._r8 )

    ! timing statements. 
    n_sub = get_days_per_year()
    udata%deltat = 1.0_r8/n_sub !for working out age of patches in years        
    if(udata%time_period == 0)then             
       udata%time_period = n_sub
    endif
    
    call get_curr_date(yr, mon, day, sec)
    ncdate = yr*10000 + mon*100 + day
    call get_ref_date(yr, mon, day, sec)
    nbdate = yr*10000 + mon*100 + day

    call timemgr_datediff(nbdate, 0, ncdate, sec, dayDiff)

    dayDiffInt = floor(dayDiff)
    udata%time_period = mod( dayDiffInt , n_sub )

    ! where most things happen
    do g = bounds%begg,bounds%endg
       if (ed_allsites_inst(g)%istheresoil) then
          currentSite => ed_allsites_inst(g)
          call ed_ecosystem_dynamics(currentSite, &
               ed_clm_inst, ed_phenology_inst, atm2lnd_inst, &
               soilstate_inst, temperature_inst, waterstate_inst)

          call ed_update_site( ed_allsites_inst(g))
       endif
    enddo

    ! link to CLM structures
    call ed_clm_inst%ed_clm_link( bounds, ed_allsites_inst(bounds%begg:bounds%endg),  &
         ed_phenology_inst, waterstate_inst, canopystate_inst)

    if (masterproc) then
      write(iulog, *) 'clm: leaving ED model', bounds%begg, bounds%endg, dayDiffInt
    end if

  end subroutine ed_driver

  !-------------------------------------------------------------------------------!
  subroutine ed_ecosystem_dynamics(currentSite, &
       ed_clm_inst, ed_phenology_inst, atm2lnd_inst, &
       soilstate_inst, temperature_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    !  Core of ed model, calling all subsequent vegetation dynamics routines         
    !
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), pointer :: currentSite
    type(ed_phenology_type) , intent(in)             :: ed_phenology_inst
    type(ed_clm_type)       , intent(in)             :: ed_clm_inst
    type(atm2lnd_type)      , intent(in)             :: atm2lnd_inst
    type(soilstate_type)    , intent(in)             :: soilstate_inst
    type(temperature_type)  , intent(in)             :: temperature_inst
    type(waterstate_type)   , intent(inout)          :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: currentPatch
    !-----------------------------------------------------------------------

    !**************************************************************************
    ! Fire, growth, biogeochemistry. 
    !**************************************************************************

    !FIX(SPM,032414) take this out.  On startup these values are all zero and on restart it
    !zeros out values read in the restart file
   
    call ed_total_balance_check(currentSite, 0)
    
    call phenology(currentSite, ed_phenology_inst, temperature_inst, waterstate_inst)

    call fire_model(currentSite, atm2lnd_inst, temperature_inst)

    ! Calculate disturbance and mortality based on previous timestep vegetation.
    call disturbance_rates(currentSite)

    ! Integrate state variables from annual rates to daily timestep
    call ed_integrate_state_variables(currentSite, soilstate_inst, temperature_inst, waterstate_inst) 

    !******************************************************************************
    ! Reproduction, Recruitment and Cohort Dynamics : controls cohort organisation 
    !******************************************************************************

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))                 

       ! adds small cohort of each PFT
       call recruitment(0,currentPatch)                

       currentPatch => currentPatch%younger
    enddo
       
    call ed_total_balance_check(currentSite,1)

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))

       ! kills cohorts that are too small
       call terminate_cohorts(currentPatch)       

       ! puts cohorts in right order
       call sort_cohorts(currentPatch)            

       ! fuses similar cohorts
       call fuse_cohorts(currentPatch)            

       currentPatch => currentPatch%younger
    enddo
   
    call ed_total_balance_check(currentSite,2)

    !*********************************************************************************
    ! Patch dynamics sub-routines: fusion, new patch creation (spwaning), termination.
    !*********************************************************************************

    ! make new patches from disturbed land
    call spawn_patches(currentSite)       
   
    call ed_total_balance_check(currentSite,3)

    ! fuse on the spawned patches.
    call fuse_patches(currentSite)        
   
    call ed_total_balance_check(currentSite,4)

    ! kill patches that are too small
    call terminate_patches(currentSite)   
   
    call ed_total_balance_check(currentSite,5)

  end subroutine ed_ecosystem_dynamics

  !-------------------------------------------------------------------------------!
  subroutine ed_integrate_state_variables(currentSite, soilstate_inst, temperature_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    ! FIX(SPM,032414) refactor so everything goes through interface
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(in)    :: currentSite
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(waterstate_type)  , intent(inout) :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort

    integer  :: c                     ! Counter for litter size class 
    integer  :: p                     ! Counter for PFT
    real(r8) :: small_no              ! to circumvent numerical errors that cause negative values of things that can't be negative
    real(r8) :: cohort_biomass_store  ! remembers the biomass in the cohort for balance checking
    !-----------------------------------------------------------------------

    small_no = 0.0000000000_r8  ! Obviously, this is arbitrary.  RF - changed to zero

    currentPatch => currentSite%youngest_patch

    do while(associated(currentPatch))

       currentPatch%age = currentPatch%age + udata%deltat
       ! FIX(SPM,032414) valgrind 'Conditional jump or move depends on uninitialised value'
       if( currentPatch%age  <  0._r8 )then
          write(iulog,*) 'negative patch age?',currentSite%clmgcell, currentPatch%age, &
               currentPatch%patchno,currentPatch%area
       endif

       ! Find the derivatives of the growth and litter processes. 
       call canopy_derivs(currentPatch)     ! BC cohort hite is updated here (so TreeHydGeom must be below this line)
       
       ! Update Canopy Biomass Pools
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 

          cohort_biomass_store  = (currentCohort%balive+currentCohort%bdead+currentCohort%bstore)
          currentCohort%dbh    = max(small_no,currentCohort%dbh    + currentCohort%ddbhdt    * udata%deltat )
          currentCohort%balive = currentCohort%balive + currentCohort%dbalivedt * udata%deltat 
          currentCohort%bdead  = max(small_no,currentCohort%bdead  + currentCohort%dbdeaddt  * udata%deltat )
          if ( DEBUG ) then
             write(iulog,*) 'EDMainMod dbstoredt I ',currentCohort%bstore, &
                  currentCohort%dbstoredt,udata%deltat
          end if
          currentCohort%bstore = currentCohort%bstore + currentCohort%dbstoredt * udata%deltat 
          if ( DEBUG ) then
             write(iulog,*) 'EDMainMod dbstoredt II ',currentCohort%bstore, &
                  currentCohort%dbstoredt,udata%deltat
          end if

          if( (currentCohort%balive+currentCohort%bdead+currentCohort%bstore)*currentCohort%n<0._r8)then
            write(iulog,*) 'biomass is negative', currentCohort%n,currentCohort%balive, &
                 currentCohort%bdead,currentCohort%bstore
          endif

          if(abs((currentCohort%balive+currentCohort%bdead+currentCohort%bstore+udata%deltat*(currentCohort%md+ &
               currentCohort%seed_prod)-cohort_biomass_store)-currentCohort%npp_acc) > 1e-8_r8)then
             write(iulog,*) 'issue with c balance in integration', abs(currentCohort%balive+currentCohort%bdead+ &
                  currentCohort%bstore+udata%deltat* &
                 (currentCohort%md+currentCohort%seed_prod)-cohort_biomass_store-currentCohort%npp_acc)
          endif  
          !do we need these any more?        
          currentCohort%npp_acc  = 0.0_r8
          currentCohort%gpp_acc  = 0.0_r8
          currentCohort%resp_acc = 0.0_r8
          
          call allocate_live_biomass(currentCohort,1)

          ! BOC...update tree 'hydraulic geometry' (size --> heights of elements --> hydraulic path lengths --> maximum node-to-node conductances)
          call hydraulics_TreeHydGeom(currentCohort)

          currentCohort => currentCohort%taller

       enddo
      
       if ( DEBUG ) then
          write(6,*)'DEBUG18: calling non_canopy_derivs with pno= ',currentPatch%clm_pno
       endif

       call non_canopy_derivs( currentPatch, temperature_inst, soilstate_inst, waterstate_inst )

       !update state variables simultaneously according to derivatives for this time period. 
       do p = 1,numpft_ed
          currentPatch%seed_bank(p) = currentPatch%seed_bank(p) + currentPatch%dseed_dt(p)*udata%deltat
       enddo

       do c = 1,ncwd
          currentPatch%cwd_ag(c) =  currentPatch%cwd_ag(c) + currentPatch%dcwd_ag_dt(c)* udata%deltat
          currentPatch%cwd_bg(c) =  currentPatch%cwd_bg(c) + currentPatch%dcwd_bg_dt(c)* udata%deltat
       enddo

       do p = 1,numpft_ed
          currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + currentPatch%dleaf_litter_dt(p)* udata%deltat
          currentPatch%root_litter(p) = currentPatch%root_litter(p) + currentPatch%droot_litter_dt(p)* udata%deltat
       enddo

       ! Check for negative values. Write out warning to show carbon balance. 
       do p = 1,numpft_ed
          if(currentPatch%seed_bank(p)<small_no)then
            write(iulog,*) 'negative seedbank', currentPatch%seed_bank(p)
            currentPatch%seed_bank(p) = small_no
          endif
       enddo

       do c = 1,ncwd
          if(currentPatch%cwd_ag(c)<small_no)then
            write(iulog,*) 'negative CWD_AG', currentPatch%cwd_ag(c),CurrentSite%lat,currentSite%lon
            currentPatch%cwd_ag(c) = small_no
          endif
          if(currentPatch%cwd_bg(c)<small_no)then
            write(iulog,*) 'negative CWD_BG', currentPatch%cwd_bg(c),CurrentSite%lat,CurrentSite%lon
            currentPatch%cwd_bg(c) = small_no
          endif
       enddo

       do p = 1,numpft_ed
          if(currentPatch%leaf_litter(p)<small_no)then
            write(iulog,*) 'negative leaf litter numerical error', currentPatch%leaf_litter(p),CurrentSite%lat,CurrentSite%lon,&
            currentPatch%dleaf_litter_dt(p),currentPatch%leaf_litter_in(p),currentPatch%leaf_litter_out(p),currentpatch%age
            currentPatch%leaf_litter(p) = small_no
          endif
          if(currentPatch%root_litter(p)<small_no)then
               write(iulog,*) 'negative root litter numerical error', currentPatch%root_litter(p), &
               currentPatch%droot_litter_dt(p)* udata%deltat, &
               CurrentSite%lat,CurrentSite%lon
            currentPatch%root_litter(p) = small_no
          endif
       enddo

     
       ! update cohort number. This needs to happen after the CWD_input and seed_input calculations as they 
       ! assume the pre-mortality currentCohort%n. 
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 
         currentCohort%n = max(small_no,currentCohort%n + currentCohort%dndt * udata%deltat )  
         currentCohort => currentCohort%taller
       enddo

       currentPatch => currentPatch%older

    enddo

    ! BOC...update 'rhizosphere geometry' (column-level root biomass + rootfr --> root length density --> node radii and volumes)
    call hydraulics_RhizHydGeom(currentSite, soilstate_inst, waterstate_inst)

  end subroutine ed_integrate_state_variables

  !-------------------------------------------------------------------------------!
  subroutine ed_update_site( currentSite )
    !
    ! !DESCRIPTION:
    ! Calls routines to consolidate the ED growth process.
    ! Canopy Structure to assign canopy layers to cohorts
    ! Canopy Spread to figure out the size of tree crowns
    ! Trim_canopy to figure out the target leaf biomass. 
    ! Extra recruitment to fill empty patches.  
    !
    ! !USES:
    use EDCanopyStructureMod , only : canopy_spread, canopy_structure
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer :: currentPatch   
    integer :: cohort_number ! To print out the number of cohorts.  
    integer :: g             ! Counter for sites
    !-----------------------------------------------------------------------

    call canopy_spread(currentSite)

    call ed_total_balance_check(currentSite,6)

    call canopy_structure(currentSite)

    call ed_total_balance_check(currentSite,7)

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))

       call terminate_cohorts(currentPatch) 

       ! FIX(SPM,040314) why is this needed for BFB restarts? Look into this at some point
       cohort_number = count_cohorts(currentPatch)  
       if ( DEBUG ) then
          write(iulog,*) 'tempCount ',cohort_number
       endif

       ! Note (RF)
       ! This breaks the balance check, but if we leave it out, then 
       ! the first new patch that isn't fused has no cohorts at the end of the spawn process
       ! and so there are radiation errors instead. 
       ! Fixing this would likely require a re-work of how seed germination works which would be tricky. 
       if(currentPatch%countcohorts < 1)then
          !write(iulog,*) 'ED: calling recruitment for no cohorts',currentPatch%siteptr%clmgcell,currentPatch%patchno
          !call recruitment(1,currentPatch)
          ! write(iulog,*) 'patch empty',currentPatch%area,currentPatch%age
       endif

       currentPatch => currentPatch%younger    

    enddo

    ! FIX(RF,032414). This needs to be monthly, not annual
    if((udata%time_period == N_SUB-1))then 
       write(iulog,*) 'calling trim canopy' 
       call trim_canopy(currentSite)  
    endif

  end subroutine ed_update_site

  !-------------------------------------------------------------------------------!
  subroutine ed_total_balance_check (currentSite, call_index )
    !
    ! !DESCRIPTION:
    ! This routine looks at the carbon in and out of the ED model and compares it to 
    ! the change in total carbon stocks. 
    ! Fluxes in are NPP. Fluxes out are decay of CWD and litter into SOM pools.  
    ! ed_allsites_inst%flux_out and ed_allsites_inst%flux_in are set where they occur 
    ! in the code. 
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout) :: currentSite
    integer            , intent(in)    :: call_index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: biomass_stock   ! total biomass   in KgC/site
    real(r8) :: litter_stock    ! total litter    in KgC/site
    real(r8) :: seed_stock      ! total seed mass in KgC/site
    real(r8) :: total_stock     ! total ED carbon in KgC/site
    real(r8) :: change_in_stock ! Change since last time we set ed_allsites_inst%old_stock in this routine.  KgC/site
    real(r8) :: error           ! How much carbon did we gain or lose (should be zero!) 
    real(r8) :: net_flux        ! Difference between recorded fluxes in and out. KgC/site

    ! nb. There is no time associated with these variables 
    ! because this routine can be called between any two 
    ! arbitrary points in code, even if no time has passed. 
    ! Also, the carbon pools are per site/gridcell, so that 
    ! we can account for the changing areas of patches. 

    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort
    !-----------------------------------------------------------------------

    change_in_stock = 0.0_r8
    biomass_stock   = 0.0_r8
    litter_stock    = 0.0_r8
    seed_stock      = 0.0_r8

    if (currentSite%istheresoil) then
       currentPatch => currentSite%oldest_patch 
       do while(associated(currentPatch))

          litter_stock = litter_stock + currentPatch%area * (sum(currentPatch%cwd_ag)+ &
               sum(currentPatch%cwd_bg)+sum(currentPatch%leaf_litter)+sum(currentPatch%root_litter))
          seed_stock   = seed_stock   + currentPatch%area * sum(currentPatch%seed_bank)
          currentCohort => currentPatch%tallest;

          do while(associated(currentCohort))

             biomass_stock =  biomass_stock + (currentCohort%bdead + currentCohort%balive + &
                  currentCohort%bstore) * currentCohort%n
             currentCohort => currentCohort%shorter;

          enddo !end cohort loop 

          currentPatch => currentPatch%younger

       enddo !end patch loop

    endif

    total_stock     = biomass_stock + seed_stock +litter_stock
    change_in_stock = total_stock - currentSite%old_stock  
    net_flux        = currentSite%flux_in - currentSite%flux_out
    error           = abs(net_flux - change_in_stock)   

    if ( abs(error) > 10e-6 ) then
       write(iulog,*) 'total error:in,out,net,dstock,error',call_index, currentSite%flux_in, & 
            currentSite%flux_out,net_flux,change_in_stock,error
       write(iulog,*) 'biomass,litter,seeds', biomass_stock,litter_stock,seed_stock
       write(iulog,*) 'lat lon',currentSite%lat,currentSite%lon
    endif

    currentSite%flux_in   = 0.0_r8
    currentSite%flux_out  = 0.0_r8  
    currentSite%old_stock = total_stock

  end subroutine ed_total_balance_check

  !-------------------------------------------------------------------------------!
  subroutine hydraulics_TreeHydGeom(cc_p)
    !
    ! !DESCRIPTION: Updates absorbing root length (total and its vertical distribution)
    !   as well as the consequential change in the size of the 'representative' rhizosphere
    !   shell radii, volumes
    !
    ! !USES:
    use pftconMod          , only : pftcon
    use shr_const_mod      , only : SHR_CONST_PI
    use EDEcophysConType   , only : EDecophyscon
    use clm_varpar         , only : nlevsoi
    use clm_varcon         , only : zisoi
    use EDTypesMod         , only : npool_leaf, npool_stem, npool_troot, npool_ag, npool_bg, nlevsoi_hyd
    !
    ! !ARGUMENTS:
    type(ed_cohort_type)   , intent(inout), target  :: cc_p ! current cohort pointer
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: currentCohort
    type(ed_patch_type),  pointer :: currentPatch
    integer  :: i,j,k,FT                     ! indices
    real(r8) :: pi                           ! pi
    real(r8) :: b_tot_carb                   ! total individual biomass in carbon units                              [kgC/indiv]
    real(r8) :: b_bg_carb                    ! belowground biomass (coarse + fine roots) in carbon units             [kgC/indiv]
    ! TRANSPORTING ROOT QUANTITIES
    real(r8) :: dcumul_rf                    ! cumulative root distribution discretization                           [-]
    real(r8) :: cumul_rf                     ! cumulative root distribution where depth is determined                [-]
    real(r8) :: z_cumul_rf                   ! depth at which cumul_rf occurs                                        [m]
    real(r8) :: b_troot_carb                 ! transporting root biomass in carbon units                             [kgC/indiv]
    real(r8) :: b_troot_biom                 ! transporting root biomass in dry wt units                             [kg/indiv]
    real(r8) :: v_troot                      ! transporting root volume                                              [m3/indiv]
    ! CANOPY or LEAF QUANTITIES
    real(r8) :: sla                          ! specific leaf area                                                    [cm2/g]
    real(r8) :: depth_canopy                 ! crown (canopy) depth                                                  [m]
    real(r8) :: dz_canopy                    ! vertical canopy discretization                                        [m]
    real(r8) :: a_leaf_tot                   ! total leaf area                                                       [m2/indiv]
    real(r8) :: b_canopy_carb                ! total leaf (canopy) biomass in carbon units                           [kgC/indiv]
    real(r8) :: b_canopy_biom                ! total leaf (canopy) biomass in dry wt units                           [kg/indiv]
    real(r8) :: v_canopy                     ! total leaf (canopy) volume                                            [m3/indiv]
    real(r8) :: denleaf                      ! leaf dry mass per unit fresh leaf volume                              [kg/m3]
    ! STEM OR SAPWOOD QUANTITIES
    real(r8) :: a_sapwood                    ! sapwood area                                                          [m2]
    real(r8) :: v_sapwood                    ! sapwood volume                                                        [m3]
    real(r8) :: z_stem                       ! tree height, minus any crown depth                                    [m]
    real(r8) :: dz_stem                      ! vertical stem discretization                                          [m]
    real(r8) :: b_stem_carb                  ! total aboveground stem biomass in carbon units                        [kgC/indiv]
    real(r8) :: b_stem_biom                  ! total agoveground stem biomass in dry wt units                        [kg/indiv]
    real(r8) :: v_stem                       ! total aboveground stem volume                                         [m3/indiv]
    ! HYDRAULIC MAXIMUM CONDUCTANCES and assoc vars
    real(r8) :: p=1._r8/3._r8                ! Savage et al. (2010) xylem taper exponent                             [-]
    real(r8) :: rint_jansenchoat=22._r8      ! conduit radius at branch location where kmax measured, tropical mean  [um]
    real(r8) :: rint_petiole=10._r8          ! petiole conduit radius (assumed invariant, sensu Savage et al. 2010)  [um]
    real(r8) :: kmax_node_petiole            ! maximum hydraulic conductivity at petiole                             [kg m-1 s-1 MPa-1]
    real(r8) :: kmax_cum_notaper(npool_ag-1) ! cumulative kmax, petiole to node k, conduit taper effects excluded    [kg s-1 MPa-1]
    real(r8) :: chi_cum_tapnotap(npool_ag-1) ! ratio of cumulative kmax with taper effects included to that without  [-]
    real(r8) :: kmax_stem                    ! total stem (petiole to transporting root node) hydraulic conductance  [kg s-1 MPa-1]
    real(r8) :: kmax_tot                     ! total tree (leaf to root tip) hydraulic conductance                   [kg s-1 MPa-1]
    real(r8) :: dz                           ! distance between any given two nodes                                  [m]
    ! INITIAL VOLUMES
    real(r8) :: v_ag_init(npool_ag)          ! initial aboveground compartment volumes                               [m3]
    real(r8) :: v_bg_init(npool_bg)          ! initial belowground compartment volumes                               [m3]
    real(r8) :: v_aroot_layer_init(nlevsoi_hyd) ! initial absorbing root compartment volumes, by soil layer          [m3]
    !-----------------------------------------------------------------------

    pi                               =  SHR_CONST_PI
    currentCohort                    => cc_p
    currentPatch                     => currentCohort%patchptr
    FT                               =  currentCohort%pft
    v_ag_init(:)                     =  currentCohort%v_ag(:)
    v_bg_init(:)                     =  currentCohort%v_bg(:)
    v_aroot_layer_init(:)            =  currentCohort%v_aroot_layer(:)
    
    ! CANOPY HEIGHTS
    !in special case where npool_leaf = 1, the node height of the canopy water pool is
    !1/2 the distance from the bottom of the canopy to the top of the tree
    !depth_canopy                     = exp(-1.169_r8)*currentCohort%hite**1.098_r8    !! crown depth from Poorter, Bongers & Bongers
    depth_canopy                     = 0._r8
    dz_canopy                        = depth_canopy / npool_leaf
    do k=1,npool_leaf
       currentCohort%z_lower_ag(k)   = currentCohort%hite - dz_canopy*k
       currentCohort%z_node_ag(k)    = currentCohort%z_lower_ag(k) + 0.5_r8*dz_canopy
       currentCohort%z_upper_ag(k)   = currentCohort%z_lower_ag(k) + dz_canopy
    enddo

    ! TOTAL LEAF VOLUME
    b_canopy_carb                    = currentCohort%bl
    b_canopy_biom                    = b_canopy_carb / pftcon%ccontent(FT)
    sla                              = 1.e4_r8 / pftcon%lma(FT)
    denleaf                          = -2.3231_r8*sla + 781.899_r8    !! empirical regression data from leaves at Caxiuana (~ 8 spp)
    v_canopy                         = b_canopy_biom / denleaf
    
    ! STEM HEIGHTS
    !in special case where npool_stem = 1, the node height of the stem water pool is
    !1/2 the height from the ground to the bottom of the canopy
    z_stem                           = currentCohort%hite - depth_canopy
    dz_stem                          = z_stem / npool_stem
    do k=npool_leaf+1,npool_ag
       currentCohort%z_upper_ag(k)   = (npool_ag - k + 1)*dz_stem
       currentCohort%z_node_ag(k)    = currentCohort%z_upper_ag(k) - 0.5_r8*dz_stem
       currentCohort%z_lower_ag(k)   = currentCohort%z_upper_ag(k) - dz_stem
    enddo

    ! STEM & SAPWOOD VOLUME
    b_stem_carb                      = b_tot_carb - b_bg_carb - b_canopy_carb
    b_stem_biom                      = b_stem_carb / pftcon%ccontent(FT)
    v_stem                           = b_stem_biom / (EDecophyscon%wood_density(FT)*1.e3_r8)
    a_leaf_tot                       = b_canopy_biom * sla * 1.e3_r8 / 1.e4_r8
    a_sapwood                        = a_leaf_tot / pftcon%latosa(FT) / 1.e4_r8
    v_sapwood                        = a_sapwood * z_stem
    
    ! TRANSPORTING ROOT DEPTHS
    !in special case where npool_troot = npool_bg = 1, the node depth of the single troot pool
    !is the depth at which 50% total root distribution is attained
    dcumul_rf                        = 1._r8/npool_troot
    do k=1,npool_troot
       cumul_rf                      = dcumul_rf*k
       call bisect_rootfr(pftcon%roota_par(FT), pftcon%roota_par(FT), 0._r8, abs(zisoi(nlevsoi)), &
                          0.001_r8, 0.001_r8, cumul_rf, z_cumul_rf)
       currentCohort%z_lower_bg(k)   = -z_cumul_rf
       call bisect_rootfr(pftcon%roota_par(FT), pftcon%roota_par(FT), 0._r8, abs(zisoi(nlevsoi)), &
                          0.001_r8, 0.001_r8, cumul_rf-0.5_r8*dcumul_rf, z_cumul_rf)
       currentCohort%z_node_bg(k)    = -z_cumul_rf
       currentCohort%z_upper_bg(k)   = 2._r8*currentCohort%z_node_bg(k) - currentCohort%z_lower_bg(k)
    enddo
    
    ! TRANSPORTING ROOT VOLUME
    !Determine total belowground biomass as a function of total (live and dead) aboveground biomass
    !then subtract out the fine root biomass to get coarse (transporting) root biomass
    b_tot_carb                       = currentCohort%bsw + currentCohort%bdead + currentCohort%bl + currentCohort%br
    b_bg_carb                        = pftcon%rootshoot(FT)/(1._r8 + pftcon%rootshoot(FT)) * b_tot_carb
    b_troot_carb                     = b_bg_carb - currentCohort%br
    b_troot_biom                     = b_troot_carb / pftcon%ccontent(FT)
    v_troot                          = b_troot_biom / (EDecophyscon%wood_density(FT)*1.e3_r8)
    
    ! COHORT-LEVEL HYDRAULIC CONTINUUM VOLUMES
    currentCohort%v_ag(1:npool_leaf)            = v_canopy / npool_leaf
    currentCohort%v_ag((npool_leaf+1):npool_ag) = v_sapwood / npool_stem
    currentCohort%v_bg(:)                       = v_troot / npool_troot    !! BOC not sure if/how we should multiply this by the sapwood fraction
    currentCohort%v_aroot_tot                   = currentCohort%br/pftcon%ccontent(FT)/pftcon%rootdens(FT)
    currentCohort%l_aroot_tot                   = currentCohort%v_aroot_tot/(pi*pftcon%rs1(FT)**2)
    currentCohort%l_aroot_layer(:)              = currentPatch%rootfr_ft(FT,:)*currentCohort%l_aroot_tot
    currentCohort%v_aroot_layer(:)              = currentPatch%rootfr_ft(FT,:)*currentCohort%v_aroot_tot
    
    ! MAXIMUM (SIZE-DEPENDENT) HYDRAULIC CONDUCTANCES
    ! first estimate cumulative (petiole to node k) conductances without taper as well as the chi taper function
    kmax_node_petiole                     = pftcon%kmax_node(FT,2) * (rint_petiole/rint_jansenchoat) ** 4._r8
    do k=npool_leaf,npool_ag
       if(k < npool_ag) then
          dz                              = currentCohort%z_node_ag(1) - currentCohort%z_node_ag(k+1)
       else
          dz                              = currentCohort%z_node_ag(1) - currentCohort%z_node_bg(1)
       end if
       kmax_cum_notaper(k)                = kmax_node_petiole * a_sapwood / dz
       chi_cum_tapnotap(k)                = xylemtaper(p, dz)
    enddo
    ! then calculate the conductances at node boundaries as the difference of cumulative conductances
    do k=npool_leaf,npool_ag
       if(k == npool_leaf) then
          currentCohort%kmax_bound(k)     = kmax_cum_notaper(k) * chi_cum_tapnotap(k)
       else
          currentCohort%kmax_bound(k)     = ( 1._r8/(kmax_cum_notaper(k)  *chi_cum_tapnotap(k)  ) - &
	                                      1._r8/(kmax_cum_notaper(k-1)*chi_cum_tapnotap(k-1))     ) ** (-1._r8)
       end if
    enddo
    ! finally, estimate the remaining tree conductance belowground as a residual
    kmax_stem                             = sum(currentCohort%kmax_bound(npool_leaf:npool_ag))
    kmax_tot                              = pftcon%rfrac_stem(FT) * kmax_stem
    currentCohort%kmax_treebg_tot         = ( 1._r8/kmax_tot - 1._r8/kmax_stem ) ** (-1._r8)
    do j=1,nlevsoi
       currentCohort%kmax_treebg_layer(j) = currentCohort%kmax_treebg_tot * currentPatch%rootfr_ft(FT,j)
    enddo

    ! UPDATE WATER CONTENTS (assume water for growth comes from within tissue itself -- apply water mass conservation)
    do k=1,npool_ag
       currentCohort%th_ag(k)             = currentCohort%th_ag(k)    * v_ag_init(k)          / currentCohort%v_ag(k)
    enddo
    do k=1,npool_bg
       currentCohort%th_bg(k)             = currentCohort%th_bg(k)    * v_bg_init(k)          / currentCohort%v_bg(k)
    enddo
    do j=1,nlevsoi_hyd
       currentCohort%th_aroot(j)          = currentCohort%th_aroot(j) * v_aroot_layer_init(j) / currentCohort%v_aroot_layer(j)
    enddo
    
    ! UPDATES OF WATER POTENTIALS ARE DONE PRIOR TO RICHARDS' SOLUTION WITHIN EDPLANTHYDRAULICSMOD.F90

  end subroutine hydraulics_TreeHydGeom

  !-------------------------------------------------------------------------------!
  subroutine bisect_rootfr(a, b, lower_init, upper_init, xtol, ytol, crootfr, x_new)
    ! 
    ! !DESCRIPTION: Bisection routine for getting the inverse of the cumulative root
    !  distribution. No analytical soln bc crootfr ~ exp(ax) + exp(bx).
    !
    ! !USES:
    !
    ! !ARGUMENTS
    real(r8)      , intent(in)     :: a, b        ! pft root distribution constants
    real(r8)      , intent(in)     :: lower_init  ! lower bound of initial x estimate [m]
    real(r8)      , intent(in)     :: upper_init  ! upper bound of initial x estimate [m]
    real(r8)      , intent(in)     :: xtol        ! error tolerance for x_new         [m]
    real(r8)      , intent(in)     :: ytol        ! error tolerance for crootfr       [-]
    real(r8)      , intent(in)     :: crootfr     ! cumulative root fraction at x_new [-]
    real(r8)      , intent(out)    :: x_new       ! soil depth                        [m]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lower                  ! lower bound x estimate [m]
    real(r8) :: upper                  ! upper bound x estimate [m]
    real(r8) :: y_lo                   ! corresponding y value at lower
    real(r8) :: f_lo                   ! y difference between lower bound guess and target y
    real(r8) :: y_hi                   ! corresponding y value at upper
    real(r8) :: f_hi                   ! y difference between upper bound guess and target y
    real(r8) :: y_new                  ! corresponding y value at x.new
    real(r8) :: f_new                  ! y difference between new y guess at x.new and target y
    real(r8) :: chg                    ! difference between x upper and lower bounds (approach 0 in bisection)
    !----------------------------------------------------------------------
    
    lower = lower_init
    upper = upper_init
    f_lo  = zeng2001_crootfr(a, b, lower) - crootfr
    f_hi  = zeng2001_crootfr(a, b, upper) - crootfr
    chg   = upper - lower
    do while(abs(chg) .gt. xtol)
       x_new = 0.5_r8*(lower + upper)
       f_new = zeng2001_crootfr(a, b, x_new) - crootfr
       if(abs(f_new) .le. ytol) then
          EXIT
       end if
       if((f_lo * f_new) .lt. 0._r8) upper = x_new
       if((f_hi * f_new) .lt. 0._r8) lower = x_new
       chg = upper - lower
    end do
  end subroutine bisect_rootfr
  
  !-------------------------------------------------------------------------------!
  subroutine hydraulics_RhizHydGeom(currentSite, soilstate_inst, waterstate_inst)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil tapped by fine roots remains
    ! the same.  
    !
    ! !USES:
    use pftconMod            , only : pftcon
    use clm_varpar           , only : nlevsoi
    use clm_varcon           , only : grav, denh2o
    use ColumnType           , only : col
    use GridcellType         , only : grc
    use shr_const_mod        , only : SHR_CONST_PI
    !use EDPlantHydraulics    , only : nshell                             !BOC...# of rhizosphere 'shells'
    !use EDPlantHydraulics    , only : swc_psi_VG, swc_satfrac_VG
    !
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(in)             :: currentSite
    type(soilstate_type)   , intent(in)             :: soilstate_inst
    type(waterstate_type)  , intent(inout)          :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort
    real(r8)                       :: l_aroot_col_init(nlevsoi)         ! initial individual total fine root length                   [m/indiv]
    real(r8)                       :: v_shell_init(nlevsoi,nshell)      ! initial volume of rhizosphere compartment                   [m3]
    real(r8)                       :: v_rhiz(nlevsoi)                   ! updated volume of all rhizosphere compartments              [m3]
    real(r8)                       :: r_delta                           ! change in radius of innermost rhizosphere compartment       [m]
    real(r8)                       :: dpsidr                            ! water potential gradient near root surface                  [MPa/m]
    real(r8)                       :: r_node_shell_init(nlevsoi,nshell) ! initial nodal radius of each rhizosphere compartment        [m]
    real(r8)                       :: w_shell_new(nlevsoi,nshell)       ! updated water mass in rhizosphere compartment               [kg]
    real(r8)                       :: w_layer_init(nlevsoi)             ! initial water mass by layer                                 [kg]
    real(r8)                       :: w_layer_interp(nlevsoi)           ! water mass after interpolating to new rhizosphere           [kg]
    real(r8)                       :: w_layer_new(nlevsoi)              ! water mass by layer after interpolation and fudging         [kg]
    real(r8)                       :: h2osoi_liq_col_new(nlevsoi)       ! water mass per area after interpolating to new rhizosphere  [kg/m2]
    real(r8)                       :: s_shell_init(nlevsoi,nshell)      ! initial saturation fraction in rhizosphere compartment      [0-1]
    real(r8)                       :: s_shell_interp(nlevsoi,nshell)    ! interpolated saturation fraction in rhizosphere compartment [0-1]
    real(r8)                       :: psi_shell_init(nlevsoi,nshell)    ! initial water potential in rhizosphere compartment          [MPa]
    real(r8)                       :: psi_shell_interp(nlevsoi,nshell)  ! interpolated psi_shell to new r_node_shell                  [MPa]
    real(r8)                       :: rh_out_rhiz(nlevsoi,nshell)       ! outer radius of rhizosphere shells                          [m]
    real(r8)                       :: delta_s(nlevsoi)                  ! change in saturation fraction needed to ensure water bal    [0-1]
    real(r8)                       :: errh2o(nlevsoi)                   ! water budget error after updating                           [kg/m2]
    real(r8)                       :: pi                                !
    real(r8)                       :: area_col                          ! column area                                                 [m2]
    real(r8)                       :: hksatMPa                          ! hksat converted to units for plant hydraulics               [kg m-1 s-1 MPa-1]
    integer                        :: g,c,j,k                           ! gridcell, column, soil layer, rhizosphere shell indicies
    integer                        :: indexc,indexj                     ! column and layer indices where there is a water balance error
    logical                        :: found                             ! flag in search loop
    !-----------------------------------------------------------------------
    
    associate(& 
         area_gcell        =>    grc%area                           , & ! Input:  [real(r8) (:)    ] gridcell area (km2)
         wtgcell           =>    col%wtgcell                        , & ! Input:  [real(r8) (:)    ] weight (relative to gridcell)
         dz                =>    col%dz                             , & ! Input:  [real(r8) (:,:)  ] layer thickness (m)                             
         r_out_shell       =>    col%r_out_shell                    , & ! Input:  [real(r8) (:,:,:)] outer radius of rhizosphere compartment (m)
         r_node_shell      =>    col%r_node_shell                   , & ! Input:  [real(r8) (:,:,:)] nodal radius of rhizosphere compartment (m)
         v_shell           =>    col%v_shell                        , & ! Input:  [real(r8) (:,:,:)] volume of rhizosphere compartment (m)
         l_aroot_col       =>    col%l_aroot_col                    , & ! Input:  [real(r8) (:,:)  ] total (across cohorts) absorbing root length by layer (m)
         h2osoi_vol_shell  =>    waterstate_inst%h2osoi_vol_shell   , & ! Input:  [real(r8) (:,:,:)] volumetric water in rhizosphere compartment (m3/m3)
         h2osoi_liq_col    =>    waterstate_inst%h2osoi_liq_col     , & ! Input:  [real(r8) (:,:)  ] col liquid water (kg/m2)
         hksat             =>    soilstate_inst%hksat_col           , & ! Input:  [real(r8) (:,:)  ] hydraulic conductivity at saturation (mm H2O /s)
         watsat            =>    soilstate_inst%watsat_col          , & ! Input:  [real(r8) (:,:)  ] col volumetric soil water at saturation (porosity)
         watres            =>    soilstate_inst%watres_col          , & ! Input:  [real(r8) (:,:)  ] col volumetric residual soil water
         alpha_VG          =>    soilstate_inst%alpha_VG_col        , & ! Input:  [real(r8) (:,:)  ] col inverse of air-entry pressure (MPa-1)
         n_VG              =>    soilstate_inst%n_VG_col            , & ! Input:  [real(r8) (:,:)  ] col pore-size distribution index
         m_VG              =>    soilstate_inst%m_VG_col            , & ! Input:  [real(r8) (:,:)  ] = 1 - 1/n_VG
         l_VG              =>    soilstate_inst%l_VG_col              & ! Input:  [real(r8) (:,:)  ] col pore tortuosity parameter
         )

    pi                      = SHR_CONST_PI
    g                       = currentSite%clmgcell
    c                       = currentSite%clmcolumn
    area_col                = area_gcell(g) * wtgcell(c) * 1.e6_r8
    l_aroot_col_init(:)     = l_aroot_col(c,:)
    r_node_shell_init(:,:)  = r_node_shell(c,:,:)
    v_shell_init(:,:)       = v_shell(c,:,:)
    s_shell_init(:,:)       = 0._r8
    psi_shell_init(:,:)     = 0._r8
    psi_shell_interp(:,:)   = 0._r8
    s_shell_interp(:,:)     = 0._r8

    ! update cohort-level root length density and accumulate it across cohorts and patches to the column level
    l_aroot_col(c,:) = 0._r8
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))
          l_aroot_col(c,:)               = l_aroot_col(c,:) + currentCohort%l_aroot_tot
          currentCohort => currentCohort%shorter
       enddo !cohort
       currentPatch => currentPatch%older
    enddo !patch
    
    ! update outer radii of column-level rhizosphere shells (same across patches and cohorts)
    do j = 1,nlevsoi
       ! proceed only if l_aroot_coh has changed
       if( l_aroot_col(c,j) /= l_aroot_col_init(j) ) then
          ! BOC...not doing this as it requires PFT-specific fine root thickness but this is at column level
          !call shellGeom(l_aroot_col(c,j), pftcon%rs1(FT), area_col, dz(c,j), r_out_shell(c,j,:), r_node_shell(c,j,:), v_shell(c,j,:))
          call shellGeom(l_aroot_col(c,j), 0.001_r8, area_col, dz(c,j), r_out_shell(c,j,:), r_node_shell(c,j,:), v_shell(c,j,:))
       end if !has l_aroot_coh changed?
    enddo
       
    ! calculate initial s, psi by layer and shell
    do j = 1,nlevsoi
       ! proceed only if l_aroot_coh has changed
       if( l_aroot_col(c,j) /= l_aroot_col_init(j) ) then
          do k = 1,nshell
             s_shell_init(j,k)     = (h2osoi_vol_shell(c,j,k) - watres(c,j)) / (watsat(c,j) - watres(c,j))
	     call swc_psi_VG(s_shell_init(j,k),alpha_VG(c,j),n_VG(c,j),m_VG(c,j),l_VG(c,j),psi_shell_init(j,k))
          enddo
       end if !has l_aroot_coh changed?
    enddo
    
    ! interpolate initial psi values by layer and shell
    ! BOC...To-Do: need to constrain psi to be within realistic limits (i.e., < 0)
    do j = 1,nlevsoi
       ! proceed only if l_aroot_coh has changed
       if( l_aroot_col(c,j) /= l_aroot_col_init(j) ) then
          if(r_node_shell(c,j,nshell) < r_node_shell_init(j,nshell)) then   ! fine root length increased, thus shrinking the rhizosphere size
             r_delta               = r_node_shell(c,j,1) - r_node_shell_init(j,1)
	     dpsidr                = (psi_shell_init(j,2) - psi_shell_init(j,1)) / (r_node_shell_init(j,2) - r_node_shell_init(j,1))
             psi_shell_interp(j,1) = dpsidr * r_delta
             do k = 2,nshell
  	        r_delta               = r_node_shell(c,j,k) - r_node_shell_init(j,k)
                dpsidr                = (psi_shell_init(j,k) - psi_shell_init(j,k-1)) / (r_node_shell_init(j,k) - r_node_shell_init(j,k-1))
                psi_shell_interp(j,k) = dpsidr * r_delta
             enddo
          else                                                              ! fine root length decreased, thus increasing the rhizosphere size
             do k = 1,(nshell-1)
	        r_delta               = r_node_shell(c,j,k) - r_node_shell_init(j,k)
                dpsidr                = (psi_shell_init(j,k+1) - psi_shell_init(j,k)) / (r_node_shell_init(j,k+1) - r_node_shell_init(j,k))
                psi_shell_interp(j,k) = dpsidr * r_delta
             enddo
             r_delta               = r_node_shell(c,j,nshell) - r_node_shell_init(j,nshell)
	     dpsidr                = (psi_shell_init(j,nshell) - psi_shell_init(j,nshell-1)) / (r_node_shell_init(j,nshell) - r_node_shell_init(j,nshell-1))
             psi_shell_interp(j,k) = dpsidr * r_delta
          end if
       end if !has l_aroot_coh changed?
    enddo
    
    ! 1st guess at new s based on interpolated psi
    do j = 1,nlevsoi
       ! proceed only if l_aroot_coh has changed
       if( l_aroot_col(c,j) /= l_aroot_col_init(j) ) then
          do k = 1,nshell
             call swc_satfrac_VG(psi_shell_interp(j,k),alpha_VG(c,j),n_VG(c,j),m_VG(c,j),l_VG(c,j),s_shell_interp(j,k))
          enddo
       end if !has l_aroot_coh changed?
    enddo
    
    ! accumlate water across shells for each layer (initial and interpolated)
    do j = 1,nlevsoi
       ! proceed only if l_aroot_coh has changed
       if( l_aroot_col(c,j) /= l_aroot_col_init(j) ) then
          w_layer_init(j)      = 0._r8
          w_layer_interp(j)    = 0._r8
          v_rhiz(j)            = 0._r8
          do k = 1,nshell
             w_layer_init(j)   = w_layer_init(j) + denh2o/dz(c,j)*( l_aroot_col_init(j)*v_shell_init(j,k)*h2osoi_vol_shell(c,j,k) )
             w_layer_interp(j) = w_layer_init(j) + denh2o/dz(c,j)*( l_aroot_col(c,j)*v_shell(c,j,k)*(s_shell_interp(j,k)*(watsat(c,j)-watres(c,j))-watres(c,j)) )
             v_rhiz(j)         = v_rhiz(j) + v_shell(c,j,k)
          enddo
       end if !has l_aroot_coh changed?
    enddo
    
    ! estimate delta_s across all shells needed to ensure total water in each layer doesn't change
    ! BOC...FIX: need to handle special cases where delta_s causes s_shell to go above or below 1 or 0, respectively.
    do j = 1,nlevsoi
       ! proceed only if l_aroot_coh has changed
       if( l_aroot_col(c,j) /= l_aroot_col_init(j) ) then
          delta_s(j) = ( w_layer_init(j) - w_layer_interp(j) )/( v_rhiz(j) * denh2o*l_aroot_col(c,j)/dz(c,j) )
       end if !has l_aroot_coh changed?
    enddo
    
    ! update h2osoi_vol_shell and h2osoi_liq_shell
    do j = 1,nlevsoi
       ! proceed only if l_aroot_coh has changed
       if( l_aroot_col(c,j) /= l_aroot_col_init(j) ) then
          w_layer_new(j)             = 0._r8
          do k = 1,nshell
             s_shell_interp(j,k)     = s_shell_interp(j,k) + delta_s(j)
   	     h2osoi_vol_shell(c,j,k) = s_shell_interp(j,k) * ( watsat(c,j)-watres(c,j) ) + watres(c,j)
	     w_shell_new(j,k)        = h2osoi_vol_shell(c,j,k) * v_shell(c,j,k) * denh2o
             w_layer_new(j)          = ( w_layer_new(j) + w_shell_new(j,k) )/( v_rhiz(j)/dz(c,j) )
          enddo
          h2osoi_liq_col_new(j)      = w_layer_new(j)/( v_rhiz(j)/dz(c,j) )
       end if !has l_aroot_coh changed?
    enddo
    
    ! balance check
    do j = 1,nlevsoi
       errh2o(j) = h2osoi_liq_col_new(j) - h2osoi_liq_col(c,j)
       if (abs(errh2o(j)) > 1.e-9_r8) then
          found = .true.
	  indexc = c
          indexj = j
          if( found ) then
                write(iulog,*)'WARNING:  water balance error ',&
                     ' local indexc= ',indexc,&
                     ' local indexj= ',indexj,&
                     ' errh2o= ',errh2o(indexj)
          end if
       end if
    enddo
    
    end associate 

  end subroutine hydraulics_RhizHydGeom
  
  !-------------------------------------------------------------------------------!
  subroutine shellGeom(l_aroot, rs1, area, dz, r_out_shell, r_node_shell, v_shell)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil surrounding fine roots remains
    ! the same.  
    !
    ! !USES:
    !use EDPlantHydraulics    , only : nshell                             !BOC...# of rhizosphere 'shells'
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
  
  !-------------------------------------------------------------------------------!
  function zeng2001_crootfr(a, b, z) result(crootfr)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: a,b    ! pft parameters
    real(r8) , intent(in) :: z      ! soil depth (m)
    !
    ! !RESULT
    real(r8) :: crootfr                         ! cumulative root fraction
    !
    !------------------------------------------------------------------------
    crootfr      = .5_r8*(exp(-a*z) + exp(-b*z))

    return

  end function zeng2001_crootfr

  !-------------------------------------------------------------------------------!
  function xylemtaper(p, dz) result(chi_tapnotap)

    ! !ARGUMENTS:
    real(r8) , intent(in) :: p      ! Savage et al. (2010) taper exponent                                                                [-]
    real(r8) , intent(in) :: dz     ! hydraulic distance from petiole to node of interest                                                [m]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: atap,btap           ! scaling exponents for total conductance ~ tree size (ratio of stem radius to terminal twig radius)
    real(r8) :: anotap,bnotap       ! same as atap, btap, but not acounting for xylem taper (Savage et al. (2010) p = 0)
                                    ! NOTE: these scaling exponents were digitized from Fig 2a of Savage et al. (2010)
				    ! Savage VM, Bentley LP, Enquist BJ, Sperry JS, Smith DD, Reich PB, von Allmen EI. 2010.
				    !    Hydraulic trade-offs and space filling enable better predictions of vascular structure
				    !    and function in plants. Proceedings of the National Academy of Sciences 107(52): 22722-22727.
    real(r8) :: lN=0.04_r8          ! petiole length                                                                                     [m]
    real(r8) :: little_n=2._r8      ! number of daughter branches per parent branch, assumed constant throughout tree (self-similarity)  [-]
    real(r8) :: big_n               ! number of branching levels (allowed here to take on non-integer values): increases with tree size  [-]
    real(r8) :: ktap                ! hydraulic conductance along the pathway, accounting for xylem taper                                [kg s-1 MPa-1]
    real(r8) :: knotap              ! hydraulic conductance along the pathway, not accounting for xylem taper                            [kg s-1 MPa-1]
    real(r8) :: num                 ! temporary
    real(r8) :: den                 ! temporary
    !
    ! !RESULT
    real(r8) :: chi_tapnotap        ! ratio of total tree conductance accounting for xylem taper to that without, over interval dz
    !
    !------------------------------------------------------------------------
    
    anotap  = 7.19903e-13_r8
    bnotap  = 1.326105578_r8
    if (p >= 1.0_r8) then
       btap  = 2.00586217_r8
       atap  = 1.82513E-12_r8
    else if (p >= (1._r8/3._r8) .AND. p < 1._r8) then
       btap  = 1.854812819_r8
       atap  = 6.66908E-13_r8
    else if (p >= (1._r8/6._r8) .AND. p < (1._r8/3._r8)) then
       btap  = 1.628179741_r8
       atap  = 6.58345E-13_r8
    else
       btap  = bnotap
       atap  = anotap
    end if

    num          = 3._r8*log(1._r8 - dz/lN * (1._r8-little_n**(1._r8/3._r8)))
    den          = log(little_n)
    big_n        = num/den - 1._r8
    ktap         = atap   * (little_n**(big_N*  btap/2._r8))
    knotap       = anotap * (little_n**(big_N*bnotap/2._r8))
    chi_tapnotap = ktap / knotap

    return

  end function xylemtaper

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

end module EDMainMod
