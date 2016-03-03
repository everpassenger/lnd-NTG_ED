module FuncPedotransferMod
!
!DESCRIPTIONS:
!module contains different pedotransfer functions to
!compute the mineral soil hydraulic properties.
!currenty, only the Clapp-Hornberg formulation is used.
!HISTORY:
!created by Jinyun Tang, Mar.1st, 2014
!pedotransf_tomasella_hodnett() added by Brad Christoffersen Jan 12th, 2016
!get_textural_class() added by Brad Christoffersen Jan 12th, 2016
implicit none
  private
  public :: pedotransf
  public :: get_ipedof
  public :: get_textural_class    ! BOC...added 01/12/16
  public :: init_pedof

  integer, parameter :: cosby_1984_table5 = 0       !by default uses this form
  integer, parameter :: cosby_1984_table4 = 1
  integer, parameter :: noilhan_lacarrere_1995 = 2
  integer, parameter :: hodnett_tomasella = 3       !BOC...Amazon-specific
  integer :: ipedof0
contains
   
!------------------------------------------------------------------------------------------
   subroutine init_pedof()
   !
   !DESCRIPTIONS
   !initialize the default pedotransfer function
   use clm_varctl    , only : use_ed_planthydraulics
   implicit none
   
   
   ipedof0 = cosby_1984_table5                     !the default pedotransfer function
   if (use_ed_planthydraulics == 1) then
      ipedof0 = hodnett_tomasella                  !BOC...overwrite pedotransfer function if Brad's hydraulics turned on
   end if
   end subroutine init_pedof
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf(ipedof, sand, clay, watsat, bsw, sucsat, xksat, watres, alpha_VG, n_VG, m_VG, l_VG)
   !pedotransfer function to compute hydraulic properties of mineral soil
   !based on input soil texture
   
   use shr_kind_mod         , only : r8 => shr_kind_r8
   use abortutils    , only : endrun         
   implicit none
   integer,  intent(in) :: ipedof   !type of pedotransfer function, use the default pedotransfer function  
   real(r8), intent(in) :: sand     !% sand
   real(r8), intent(in) :: clay     !% clay
   real(r8), intent(out):: watsat   !v/v saturate moisture         [m3 m-3] (for C&H and van Genuchten SWC)
   real(r8), intent(out):: bsw      !b shape parameter
   real(r8), intent(out):: sucsat   !mm, soil matric potential
   real(r8), intent(out):: xksat    !mm/s, saturated hydraulic conductivity
   real(r8), optional, intent(out):: watres   !v/v residual moisture         [m3 m-3] (for van Genuchten SWC only)
   real(r8), optional, intent(out):: alpha_VG !inverse of air-entry pressure [MPa-1]  (for van Genuchten SWC only)
   real(r8), optional, intent(out):: n_VG     !pore-size distribution index  [-]      (for van Genuchten SWC only)
   real(r8), optional, intent(out):: m_VG     != 1 - 1/n_VG                  [-]      (for van Genuchten SWC only)
   real(r8), optional, intent(out):: l_VG     !pore tortuosity parameter     [-]      (for van Genuchten SWC only)
   
   character(len=32) :: subname = 'pedotransf'  ! subroutine name 
   select case (ipedof)
   case (cosby_1984_table4)
      call pedotransf_cosby1984_table4(sand, clay, watsat, bsw, sucsat, xksat)
   case (hodnett_tomasella)
      call pedotransf_hodnett_tomasella(sand, clay, watsat, watres, alpha_VG, n_VG, m_VG, l_VG, bsw, sucsat, xksat)
   case (noilhan_lacarrere_1995)
      call pedotransf_noilhan_lacarrere1995(sand, clay, watsat, bsw, sucsat, xksat) 
   case (cosby_1984_table5)  
      call pedotransf_cosby1984_table5(sand, clay, watsat, bsw, sucsat, xksat)
   case default
      call endrun(subname // ':: a pedotransfer function must be specified!')  
   end select

   end subroutine pedotransf
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf_cosby1984_table4(sand, clay, watsat, bsw, sucsat, xksat)
   !
   !DESCRIPTIONS
   !compute hydraulic properties based on functions derived from Table 4 in cosby et al, 1984
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   real(r8), intent(in) :: sand   !% sand
   real(r8), intent(in) :: clay   !% clay
   real(r8), intent(out):: watsat !v/v saturate moisture
   real(r8), intent(out):: bsw    !b shape parameter
   real(r8), intent(out):: sucsat !mm, soil matric potential
   real(r8), intent(out):: xksat  !mm/s, saturated hydraulic conductivity
   
   !Cosby et al. Table 4
   watsat = 0.505_r8-0.00142_r8*sand-0.00037*clay
   bsw = 3.10+0.157*clay-0.003*sand            
   sucsat  = 10._r8 * ( 10._r8**(1.54_r8-0.0095_r8*sand+0.0063*(100._r8-sand-clay)))            
   xksat         = 0.0070556 *(10.**(-0.60+0.0126*sand-0.0064*clay) )     !mm/s now use table 4.
      
   end subroutine pedotransf_cosby1984_table4
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf_hodnett_tomasella(sand, clay, watsat, watres, alpha_VG, n_VG, m_VG, l_VG, bsw, sucsat, xksat)
   !
   !created by Brad Christoffersen
   !
   !DESCRIPTIONS
   !Compute hydraulic properties for van Genuchten soil water characteristic
   !based on values given in Table 6 of Hodnett & Tomasella 2002 Geoderma.
   !
   !As such, this is not a true pedotransfer function; the PTF advocated in H&T requires input
   !data on soil organic carbon and pH currently not assimilated as surface data in CLM.
   !
   !For tropical South American applications, an alternative is to read in maps of VG hydraulic
   !parameters directly using the dataset of Marthews et al. (2014) GMD (15 arcsec resolution)
   !
   ! !USES:
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   !
   ! !ARGUMENTS:
   real(r8), intent(in) :: sand        !% sand
   real(r8), intent(in) :: clay        !% clay
   real(r8), intent(out):: watsat      !v/v saturate moisture
   real(r8), intent(out):: watres      !v/v residual moisture             [m3 m-3] (for van Genuchten SWC only)
   real(r8), intent(out):: alpha_VG    !inverse of air-entry pressure     [MPa-1]  (for van Genuchten SWC only)
   real(r8), intent(out):: n_VG        !pore-size distribution index      [-]      (for van Genuchten SWC only)
   real(r8), intent(out):: m_VG        != 1 - 1/n_VG                      [-]      (for van Genuchten SWC only)
   real(r8), intent(out):: l_VG        !pore tortuosity parameter         [-]      (for van Genuchten SWC only)
   real(r8), intent(out):: bsw         !b shape parameter
   real(r8), intent(out):: sucsat      !mm, soil matric potential
   real(r8), intent(out):: xksat       !saturated hydraulic conductivity  [mm/s]

   ! !LOCAL VARIABLES:
   integer    :: texture_class
   
   !first determine soil textural class.
   texture_class = get_textural_class(sand, clay)

   !Table 6 Hodnett & Tomasella 2002 
   SELECT CASE (texture_class)
      CASE (1)
         alpha_VG = 0.380_r8
	 n_VG     = 2.474_r8
	 watsat   = 0.410_r8
	 watres   = 0.037_r8
      CASE (2)
         alpha_VG = 0.837_r8
	 n_VG     = 1.672_r8
	 watsat   = 0.438_r8
	 watres   = 0.062_r8
      CASE (3)
         alpha_VG = 0.396_r8
	 n_VG     = 1.553_r8
	 watsat   = 0.461_r8
	 watres   = 0.111_r8
      CASE (4)
         alpha_VG = 0.246_r8
	 n_VG     = 1.461_r8
	 watsat   = 0.521_r8
	 watres   = 0.155_r8
      CASE (5)
         alpha_VG = 0.191_r8
	 n_VG     = 1.644_r8
	 watsat   = 0.601_r8
	 watres   = 0.223_r8
      CASE (6)                  !BOC...Silt NA in Table 6. tc_silt_loam is the adjacent class on the triangle.
         alpha_VG = 0.191_r8
	 n_VG     = 1.644_r8
	 watsat   = 0.601_r8
	 watres   = 0.223_r8
      CASE (7)
         alpha_VG = 0.644_r8
	 n_VG     = 1.535_r8
	 watsat   = 0.413_r8
	 watres   = 0.149_r8
      CASE (8)
         alpha_VG = 0.392_r8
	 n_VG     = 1.437_r8
	 watsat   = 0.519_r8
	 watres   = 0.226_r8
      CASE (9)
         alpha_VG = 0.298_r8
	 n_VG     = 1.513_r8
	 watsat   = 0.586_r8
	 watres   = 0.267_r8
      CASE (10)
         alpha_VG = 0.509_r8
	 n_VG     = 1.396_r8
	 watsat   = 0.460_r8
	 watres   = 0.199_r8
      CASE (11)
         alpha_VG = 0.258_r8
	 n_VG     = 1.466_r8
	 watsat   = 0.570_r8
	 watres   = 0.278_r8
      CASE (12)
         alpha_VG = 0.463_r8
	 n_VG     = 1.514_r8
	 watsat   = 0.546_r8
	 watres   = 0.267_r8
      CASE DEFAULT
   end SELECT
   
   m_VG   = 1._r8 - 1._r8/n_VG
   l_VG   = 0.5_r8
      
   !BOC...NOTE: Current plant hydraulics routine solves for root water uptake by updating soil water contents, but uses a different pedotransfer function.
   !   As of yet, I have not implemented Jinyun's interface for using the same pedotransfer function for the vertical fluxes in SoilWater()
   !   Thus it is necessary to make sure the watsat values are consistent between the two (inconsistent!) pedotransfer functions used in plant+soil hydraulics and the CLM soil hydrology
   !   The below does 3 things:
   !      overwrites the above H&T watsat value with Cosby value
   !      ensures that bsw and sucsat parameters necessary for doing vertical water fluxes get assigned
   !      assigns xksat, for which no tropical pedotransfer function exists
   !BOC... This is not an ideal implementation and should/will be replaced with consistent pedotransfer functions for root water uptake and vertical water fluxes
   call pedotransf_cosby1984_table5(sand, clay, watsat, bsw, sucsat, xksat)

   end subroutine pedotransf_hodnett_tomasella
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf_cosby1984_table5(sand, clay, watsat, bsw, sucsat, xksat)
   !
   !DESCRIPTIONS
   !compute hydraulic properties based on functions derived from Table 5 in cosby et al, 1984
   
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   real(r8), intent(in) :: sand   !% sand
   real(r8), intent(in) :: clay   !% clay
   real(r8), intent(out):: watsat !v/v saturate moisture
   real(r8), intent(out):: bsw    !b shape parameter
   real(r8), intent(out):: sucsat !mm, soil matric potential
   real(r8), intent(out):: xksat  !mm/s, saturated hydraulic conductivity
   
   !Cosby et al. Table 5     
   watsat = 0.489_r8 - 0.00126_r8*sand
   bsw    = 2.91 + 0.159*clay
   sucsat = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )            
   xksat         = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s, from table 5 
      
   end subroutine pedotransf_cosby1984_table5
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf_noilhan_lacarrere1995(sand, clay, watsat, bsw, sucsat, xksat)
   !
   !DESCRIPTIONS
   !compute hydraulic properties based on functions derived from Noilhan and Lacarrere, 1995
   
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   real(r8), intent(in) :: sand   !% sand
   real(r8), intent(in) :: clay   !% clay
   real(r8), intent(out):: watsat !v/v saturate moisture
   real(r8), intent(out):: bsw    !b shape parameter
   real(r8), intent(out):: sucsat !mm, soil matric potential
   real(r8), intent(out):: xksat  !mm/s, saturated hydraulic conductivity
   
   !Noilhan and Lacarrere, 1995
   watsat = -0.00108*sand+0.494305
   bsw = 0.137*clay + 3.501
   sucsat = 10._r8**(-0.0088*sand+2.85)
   xksat = 10._r8**(-0.0582*clay-0.00091*sand+0.000529*clay**2._r8-0.0001203*sand**2._r8-1.38)
   end subroutine pedotransf_noilhan_lacarrere1995
!------------------------------------------------------------------------------------------
   function get_ipedof(soil_order)result(ipedof)
   !
   ! DESCRIPTION
   ! select the pedotransfer function to be used
   implicit none
   integer, intent(in) :: soil_order
   
   integer :: ipedof
   
   if(soil_order==0)then
      ipedof=ipedof0
   endif
   
   end function get_ipedof   
!------------------------------------------------------------------------------------------
   function get_textural_class(sand, clay) result(texture_class)
   !
   !created by Brad Christoffersen
   ! DESCRIPTION
   !identify where within the USDA texture triangle a particular combo of sand and clay resides
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   real(r8), intent(in) :: sand
   real(r8), intent(in) :: clay
   !
   ! !LOCAL VARIABLES:
   real(r8)           ::   silt
   integer, parameter ::   tc_sand            = 1
   integer, parameter ::   tc_loamy_sand      = 2
   integer, parameter ::   tc_sandy_loam      = 3
   integer, parameter ::   tc_loam            = 4
   integer, parameter ::   tc_silt_loam       = 5
   integer, parameter ::   tc_silt            = 6
   integer, parameter ::   tc_sandy_clay_loam = 7
   integer, parameter ::   tc_clay_loam       = 8
   integer, parameter ::   tc_silty_clay_loam = 9
   integer, parameter ::   tc_sandy_clay      = 10
   integer, parameter ::   tc_silty_clay      = 11
   integer, parameter ::   tc_clay            = 12
   integer            ::   texture_class
   
   silt   = 100._r8 - sand - silt                !% silt

   !using calculations from http://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/?cid=nrcs142p2_054167
   if ((silt + 1.5_r8*clay) < 15._r8) then
      texture_class = tc_sand
   else if (((silt + 1.5_r8*clay) >= 15._r8) .and. ((silt + 2._r8*clay) < 30._r8)) then
      texture_class = tc_loamy_sand
   else if ((clay >= 7._r8 .and. clay < 20._r8) .and. (sand > 52._r8) .and. ((silt + 2._r8*clay) >= 30._r8) &
       .or. (clay < 7._r8 .and. silt < 50._r8 .and. (silt+2._r8*clay) >= 30._r8)) then
      texture_class = tc_sandy_loam
   else if ((clay >= 7._r8 .and. clay < 27._r8) .and. (silt >= 28._r8 .and. silt < 50._r8) .and. (sand <= 52._r8)) then
      texture_class = tc_loam
   else if ((silt >= 50._r8 .and. (clay >= 12._r8 .and. clay < 27._r8)) &
       .or. ((silt >= 50._r8 .and. silt < 80._r8) .and. clay < 12._r8)) then
      texture_class = tc_silt_loam
   else if (silt >= 80._r8 .and. clay < 12._r8) then
      texture_class = tc_silt
   else if ((clay >= 20._r8 .and. clay < 35._r8) .and. (silt < 28._r8) .and. (sand > 45._r8)) then
      texture_class = tc_sandy_clay_loam
   else if ((clay >= 27._r8 .and. clay < 40._r8) .and. (sand > 20._r8 .and. sand <= 45._r8)) then
      texture_class = tc_clay_loam
   else if ((clay >= 27._r8 .and. clay < 40._r8) .and. (sand  <= 20._r8)) then
      texture_class = tc_silty_clay_loam
   else if (clay >= 35._r8 .and. sand > 45._r8) then
      texture_class = tc_sandy_clay
   else if (clay >= 40._r8 .and. silt >= 40._r8) then
      texture_class = tc_silty_clay
   else if (clay >= 40._r8 .and. sand <= 45._r8 .and. silt < 40._r8) then
      texture_class = tc_clay
   end if

   end function get_textural_class   
end module FuncpedotransferMod
