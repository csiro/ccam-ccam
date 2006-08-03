! carbon.f90
!
! Carbon store routines source file for CABLE, CSIRO land surface model
!
! The flow of carbon between the vegetation compartments and soil is described
! by a simple carbon pool model of Dickinson et al 1998, J. Climate, 11, 2823-2836.
! Model implementation by Eva Kowalczyk, CSIRO Marine and Atmospheric Research.
! 
! Fortran-95 coding by Harvey Davies and Gab Abramowitz,
! bugs to gabsun@gmail.com.

MODULE carbon_module
  USE define_types
  IMPLICIT NONE
  private
  ! Do these really need to be module variables
  ! clitt is used in eva_output, rest are just local in carbon_pl
  REAL(r_1), DIMENSION(:), allocatable, public, save :: clitt
  REAL(r_1), DIMENSION(:), allocatable, public, save :: coef_cd ! total stress coeff. for vegetation (eq. 6)
  REAL(r_1), DIMENSION(:), allocatable, public, save :: coef_cold  ! coeff. for the cold stress (eq. 7)
  REAL(r_1), DIMENSION(:), allocatable, public, save :: coef_drght ! coeff. for the drought stress (eq. 8)
  PUBLIC carbon_pl, soilcarb
CONTAINS

  SUBROUTINE carbon_pl(dt, soil, ssoil, veg, canopy, bgc)
!  SUBROUTINE carbon_pl(dt, ktau, soil, ssoil, veg, canopy, bgc)
    TYPE (soil_parameter_type), INTENT(IN)		:: soil
    TYPE (soil_snow_type), INTENT(IN)	                :: ssoil
    TYPE (veg_parameter_type), INTENT(IN)		:: veg
    TYPE (canopy_type), INTENT(IN)	         	:: canopy
    TYPE (bgc_pool_type), INTENT(INOUT)	                :: bgc
    REAL(r_1), INTENT(IN)			        :: dt ! integration time step (s)
!    INTEGER(i_d), INTENT(IN)                           :: ktau  ! integration step number
    REAL(r_1), PARAMETER        :: beta = 0.9
    REAL(r_1), DIMENSION(mp)	:: cfsf     ! fast soil carbon turnover
    REAL(r_1), DIMENSION(mp)	:: cfrts    ! roots turnover
    REAL(r_1), DIMENSION(mp)	:: cfwd     ! wood turnover 
    REAL(r_1), DIMENSION(mp)	:: fcl ! fraction of assimilated carbon that goes to the
    !					construction of leaves  (eq. 5)
    REAL(r_1), DIMENSION(mp)	:: fr
    REAL(r_1), PARAMETER, DIMENSION(13) :: rw = (/16., 8.7, 12.5, 16., 18., 7.5, 6.1, .84, &
                 10.4, 15.1, 9., 5.8, 0.001 /) ! approximate ratio of wood to nonwood carbon
    !	         				 inferred from observations 
    REAL(r_1), PARAMETER, DIMENSION(13) :: tfcl = (/0.248, 0.345, 0.31, 0.42, 0.38, 0.35, &
         0.997,	0.95, 2.4, 0.73, 2.4, 0.55, 0.9500/)         ! leaf allocation factor
    REAL(r_1), PARAMETER, DIMENSION(13) :: trnl = 3.17e-8    ! leaf turnover rate 1 year
    REAL(r_1), PARAMETER, DIMENSION(13) :: trnr = 4.53e-9    ! root turnover rate 7 years
    REAL(r_1), PARAMETER, DIMENSION(13) :: trnsf = 1.057e-10 ! soil transfer rate coef. 30 years
    REAL(r_1), PARAMETER, DIMENSION(13) :: trnw = 6.342e-10  ! wood turnover 50  years
    REAL(r_1), PARAMETER, DIMENSION(13) :: tvclst = (/ 283., 278., 278., 235., 268., &
                                           278.0, 278.0, 278.0, 278.0, 235., 278., &
                                           278., 268. /) ! cold stress temp. below which 
!                                                         leaf loss is rapid
    REAL(r_1), DIMENSION(mp)	:: wbav ! water stress index

    if ( .not. allocated(clitt) ) then
       allocate ( clitt(mp), coef_cd(mp), coef_cold(mp), coef_drght(mp) )
    end if

    !
    ! coef_cold = EXP(-(canopy%tv - tvclst(veg%iveg)))   ! cold stress
    ! Limit size of exponent to avoif overflow when tv is very cold
    coef_cold = EXP(min(50., -(canopy%tv - tvclst(veg%iveg))))   ! cold stress
    wbav = SUM(soil%froot * ssoil%wb, 2)
    coef_drght = exp(10.*( min(1., max(1.,wbav**(2-soil%ibp2)-1.) / & ! drought stress
         (soil%swilt**(2-soil%ibp2) - 1.)) - 1.))
    coef_cd = ( coef_cold + coef_drght ) * 2.0e-7
    !
    ! CARBON POOLS
    !
    fcl = exp(-tfcl(veg%iveg) * veg%vlai)  ! fraction of assimilated carbon that goes
    !                                         to the construction of leaves  (eq. 5)

    !							 LEAF
    ! resp_lfrate is omitted below as fpn represents photosythesis - leaf transpiration
    ! calculated by the CBM 
    !
    clitt = (coef_cd + trnl(veg%iveg)) * bgc%cplant(:,1)
    bgc%cplant(:,1) = bgc%cplant(:,1) - dt * (canopy%fpn * fcl + clitt)
    !
    !							 WOOD
    !	                           fraction of photosynthate going to roots, (1-fr) to wood, eq. 9
    fr = min(1., exp(- rw(veg%iveg) * beta * bgc%cplant(:,3) / max(bgc%cplant(:,2), 0.01)) / beta)
    !
    !                                            
    cfwd = trnw(veg%iveg) * bgc%cplant(:,2)
    bgc%cplant(:,2) = bgc%cplant(:,2) - dt * (canopy%fpn * (1.-fcl) * (1.-fr) + canopy%frpw + cfwd )

    !							 ROOTS
    !				
    cfrts = trnr(veg%iveg) * bgc%cplant(:,3)
    bgc%cplant(:,3) = bgc%cplant(:,3) - dt * (canopy%fpn * (1. - fcl) * fr + cfrts + canopy%frpr )
    !
    !							 SOIL
    !			                                	fast carbon 
    cfsf = trnsf(veg%iveg) * bgc%csoil(:,1)
    bgc%csoil(:,1) = bgc%csoil(:,1) + dt * (0.98 * clitt + 0.9 * cfrts + cfwd  - cfsf &
         - 0.98 * canopy%frs)
    !			                                	slow carbon 
    bgc%csoil(:,2) = bgc%csoil(:,2) + dt * (0.02 * clitt  + 0.1 * cfrts + cfsf &
         - 0.02 * canopy%frs)

    bgc%cplant(:,1)  = max(0.001, bgc%cplant(:,1))
    bgc%cplant(:,2)  = max(0.001, bgc%cplant(:,2))
    bgc%cplant(:,3) = max(0.001, bgc%cplant(:,3))
    bgc%csoil(:,1) = max(0.001, bgc%csoil(:,1))
    bgc%csoil(:,2) = max(0.001, bgc%csoil(:,2))
  END SUBROUTINE carbon_pl

  SUBROUTINE soilcarb(soil, ssoil, veg, bgc, met, canopy)
    TYPE (soil_parameter_type), INTENT(IN) :: soil
    TYPE (soil_snow_type), INTENT(IN)	:: ssoil
    TYPE (veg_parameter_type), INTENT(IN)  :: veg
    TYPE (bgc_pool_type), INTENT(IN)	:: bgc
    TYPE (met_type), INTENT(IN)		:: met	
    TYPE (canopy_type), INTENT(OUT)	:: canopy
    REAL(r_1), DIMENSION(mp)		:: den ! sib3
    INTEGER(i_d)			:: k
    REAL(r_1), DIMENSION(mp)		:: rswc
    REAL(r_1), DIMENSION(mp)		:: sss
    REAL(r_1), DIMENSION(mp)		:: e0rswc
    REAL(r_1), DIMENSION(mp)		:: ftsoil
    REAL(r_1), DIMENSION(mp)		:: ftsrs
    REAL(r_1), PARAMETER, DIMENSION(13)	:: rswch = 0.16
    REAL(r_1), PARAMETER, DIMENSION(13)	:: soilcf = 1.0
    REAL(r_1), PARAMETER		:: t0 = -46.0
    REAL(r_1), DIMENSION(mp)		:: tref
    REAL(r_1), DIMENSION(mp)		:: tsoil
    REAL(r_1), PARAMETER, DIMENSION(13)	:: vegcf = &
         (/ 1.95, 1.5, 1.55, 0.91, 0.73, 2.8, 2.75, 0.0, 2.05, 0.6, 0.4, 2.8, 0.0 /)

    den = soil%sfc - soil%swilt     
    rswc =  max(0.0001, soil%froot(:,1)*(ssoil%wb(:,2) - soil%swilt)) / den
!    rswc = soil%froot(:, 1) * max(0.0001, ssoil%wb(:,2) - soil%swilt) / den
    tsoil = soil%froot(:,1) * ssoil%tgg(:,2) - 273.15
    tref = max(t0 + 1.,ssoil%tgg(:,ms) - 273.1)

    do k = 2,ms  ! start from 2nd index for less dependency on the surface condition
       rswc = rswc +  &
            max(0.0001, soil%froot(:,k) * (ssoil%wb(:,k) - soil%swilt)) / den
       tsoil = tsoil + soil%froot(:,k) * ssoil%tgg(:,k)
    enddo
    rswc = min(1.,rswc)
    tsoil = max(t0 + 2., tsoil)
    e0rswc = 52.4 + 285. * rswc
    ftsoil=1./(tref - t0) - 1./(tsoil - t0)
    sss = min(15.,e0rswc * ftsoil)
    ftsrs=exp(sss)
    !    ftsrs=exp(e0rswc * ftsoil)
    !        soiltref=soilcf(soil%isoilm)*min(1.,1.4*max(.3,.0278*tsoil+.5))
    !        rpsoil=vegcf(veg%iveg)*soiltref* ftsrs * frswc
    !     &              * 12./44.*12./1.e6 * 
!    print *,'rswc',rswc
!    print *,'sss',sss
!    print *,'ftsrs',ftsrs
    canopy%frs = vegcf(veg%iveg) * (144.0 / 44.0e6)  &
         * soilcf(soil%isoilm) * min(1.,1.4 * max(.3,.0278 * tsoil + .5)) &
         * ftsrs &
         * rswc / (rswch(soil%isoilm) + rswc)
!         * exp((52.4 + 285. * rswc) * (1. / (tref - t0) - 1. / (tsoil - t0))) &
 
  END SUBROUTINE soilcarb

END MODULE carbon_module
