
MODULE carbon_module
  USE define_types
  USE define_dimensions, ONLY: r_1,i_d,mp,ms,mvtype,mstype
  IMPLICIT NONE
  PRIVATE
  PUBLIC carbon_pl, soilcarb, plantcarb
CONTAINS

  SUBROUTINE carbon_pl(dels, soil, ssoil, veg, canopy, bgc)
      use cable_common_module, only : cable_runtime, cable_user
      use cable_diag_module, only : cable_stat
   implicit none
    REAL(r_1), INTENT(IN)                 :: dels     ! integration time step (s)
    TYPE (soil_parameter_type), INTENT(IN):: soil   ! soil parameters
    TYPE (soil_snow_type), INTENT(IN)     :: ssoil  ! soil/snow variables
    TYPE (veg_parameter_type), INTENT(IN) :: veg    ! vegetation parameters
    TYPE (canopy_type), INTENT(IN)        :: canopy ! canopy/veg variables
    TYPE (bgc_pool_type), INTENT(INOUT)   :: bgc    ! biogeochemistry variables
!    INTEGER(i_d), INTENT(IN)              :: mvtype  ! Number of vegetation types
    REAL(r_1), PARAMETER     :: beta = 0.9
    REAL(r_1), DIMENSION(mp) :: cfsf     ! fast soil carbon turnover
    REAL(r_1), DIMENSION(mp) :: cfrts    ! roots turnover
    REAL(r_1), DIMENSION(mp) :: cfwd     ! wood turnover
    REAL(r_1), DIMENSION(mp) :: fcl ! fraction of assimilated carbon that
                                ! goes to the construction of leaves  (eq. 5)
    REAL(r_1), DIMENSION(mp) :: fr
    REAL(r_1), DIMENSION(mp) :: clitt
    REAL(r_1), DIMENSION(mp) :: coef_cd ! total stress coeff. for veg (eq. 6)
    REAL(r_1), DIMENSION(mp) :: coef_cold  ! coeff. for cold stress (eq. 7)
    REAL(r_1), DIMENSION(mp) :: coef_drght ! coeff. for drought stress (eq. 8)
    REAL(r_1), DIMENSION(:), ALLOCATABLE :: rw, tfcl, tvclst
!jhan:see earliervn for comment
    REAL(r_1), DIMENSION(:), ALLOCATABLE :: trnl, trnr, trnsf, trnw

    REAL(r_1), DIMENSION(mp) :: wbav ! water stress index 

      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('carbon_pl')

    ALLOCATE( rw(mvtype), tfcl(mvtype), tvclst(mvtype) )

    ALLOCATE( trnl(mvtype), trnr(mvtype), trnsf(mvtype), trnw(mvtype) )
    trnl = 3.17e-8
    trnr = 4.53e-9
    trnsf = 1.057e-10
    trnw = 6.342e-10

    SELECT CASE (mvtype)
      CASE (13)     ! CASA vegetation types
        rw   = (/ 16., 8.7, 12.5, 16., 18., 7.5, &
                & 6.1, .84, 10.4, 15.1, 9., 5.8, 0.001 /)
        tfcl = (/ 0.248, 0.345, 0.31, 0.42, 0.38, 0.35, &
                & 0.997, 0.95, 2.4, 0.73, 2.4, 0.55, 0.9500 /)
        tvclst = (/ 283., 278., 278., 235., 268., 278.0, &
                  & 278.0, 278.0, 278.0, 235., 278., 278., 268. /)
      CASE (15)     ! CSIRO types for UM
        rw   = (/ 16., 16., 18., 8.7, 10.4, 6.1, 6.1, 6.1, &
                  5.8, 5.8, 0.001, 9.0, 0.001, 0.001, 0.001 /)
        tfcl = (/ 0.42, 0.248, 0.38, 0.345, 2.4, 0.997, 0.997, 0.997, &
                  0.55, 0.55, 0.9500, 2.4, 0.9500, 0.9500, 0.9500 /)
        tvclst = (/ 235., 283., 268., 278., 278.0, 278.0, 278.0, 278.0, &
                    278., 278., 278.0, 278., 278., 278., 268. /)
      CASE (16)     ! IGBP vegetation types without water bodies
        rw   = (/ 16., 16., 18., 8.7, 12.5, 15.1, 10.4, 7.5, &
                & 6.1, 6.1, 0.001, 5.8, 0.001, 5.8, 0.001, 9.0 /)
        tfcl = (/ 0.42, 0.248, 0.38, 0.345, 0.31, 0.73, 2.4, 0.35, &
                & 0.997, 0.997, 0.9500, 0.55, 0.9500, 0.55, 0.9500, 2.4 /)
        tvclst = (/ 235., 283., 268., 278., 278., 235., 278.0, 278.0, &
                  & 278.0, 278.0, 278.0, 278., 278., 278., 268., 278. /)
      CASE (17)     ! IGBP vegetation types with water bodies
!! rml: may not be the best values for our current 17 types, but will be superceeded by 
!! CASA-CNP anyway
        rw   = (/ 16., 16., 18., 8.7, 12.5, 15.1, 10.4, 7.5, &
                & 6.1, 6.1, 0.001, 5.8, 0.001, 5.8, 0.001, 9.0, 0.001 /)
        tfcl = (/ 0.42, 0.248, 0.38, 0.345, 0.31, 0.73, 2.4, 0.35, 0.997, &
                & 0.997, 0.9500, 0.55, 0.9500, 0.55, 0.9500, 2.4, 0.9500 /)
        tvclst = (/ 235., 283., 268., 278., 278., 235., 278.0, 278.0, &
                  & 278.0, 278.0, 278.0, 278., 278., 278., 268., 278., 278. /)
      CASE DEFAULT
        PRINT *, 'Error! Dimension not compatible with CASA or CSIRO or IGBP types!'
        PRINT *, 'Dimension =', mvtype
        PRINT *, 'At the rw section.'
        STOP
    END SELECT

      ! Limit size of exponent to avoif overflow when tv is very cold
      coef_cold = EXP(MIN(1., -(canopy%tv - tvclst(veg%iveg)))) 
      wbav = REAL(SUM(veg%froot * ssoil%wb, 2),r_1)
      wbav = max(0.01,wbav)  ! EAK Jan2011
      ! drought stress
      coef_drght = EXP(5.*( MIN(1., MAX(1.,wbav**(2-soil%ibp2)-1.) / & 
            (soil%swilt**(2-soil%ibp2) - 1.)) - 1.))
      coef_cd = ( coef_cold + coef_drght ) * 2.0e-7
      !
      ! CARBON POOLS
      !
      fcl = EXP(-tfcl(veg%iveg) * veg%vlai)  ! fraction of assimilated carbon
                             ! that goes to the construction of leaves (eq. 5)
      !	LEAF
      ! resp_lfrate is omitted below as fpn represents photosythesis - leaf
      ! transpiration calculated by the CBM 
      !
      clitt = (coef_cd + trnl(veg%iveg)) * bgc%cplant(:,1)
      bgc%cplant(:,1) = bgc%cplant(:,1) - dels * (canopy%fpn * fcl + clitt)
      !
      !	WOOD
      ! fraction of photosynthate going to roots, (1-fr) to wood, eq. 9
      fr = MIN(1., EXP(- rw(veg%iveg) * beta * bgc%cplant(:,3) &
            & / MAX(bgc%cplant(:,2), 0.01)) / beta)
      !
      !                                            
      cfwd = trnw(veg%iveg) * bgc%cplant(:,2)
      bgc%cplant(:,2) = bgc%cplant(:,2) - dels * (canopy%fpn * (1.-fcl) &
                    & * (1.-fr) + canopy%frpw + cfwd )

      ! ROOTS
      !				
      cfrts = trnr(veg%iveg) * bgc%cplant(:,3)
      bgc%cplant(:,3) = bgc%cplant(:,3) - dels * (canopy%fpn * (1. - fcl) * fr &
                    & + cfrts + canopy%frpr )
      !
      ! SOIL
      !     fast carbon 
      cfsf = trnsf(veg%iveg) * bgc%csoil(:,1)
      bgc%csoil(:,1) = bgc%csoil(:,1) + dels * (0.98 * clitt + 0.9 * cfrts &
                   & + cfwd - cfsf - 0.98 * canopy%frs)
      !     slow carbon 
      bgc%csoil(:,2) = bgc%csoil(:,2) + dels * (0.02 * clitt  + 0.1 * cfrts &
                   & + cfsf - 0.02 * canopy%frs)
       
! rml 17/1/11 change minimum pool size from 0.001 to 0. (since want 0. for vegtype=ice)
      bgc%cplant(:,1)  = MAX(0.00, bgc%cplant(:,1))
      bgc%cplant(:,2)  = MAX(0.00, bgc%cplant(:,2))
      bgc%cplant(:,3) = MAX(0.00, bgc%cplant(:,3))
      bgc%csoil(:,1) = MAX(0.00, bgc%csoil(:,1))
      bgc%csoil(:,2) = MAX(0.00, bgc%csoil(:,2))
  
    DEALLOCATE( rw, tfcl, tvclst )
    DEALLOCATE( trnl, trnr, trnsf, trnw )

  END SUBROUTINE carbon_pl

  SUBROUTINE soilcarb( soil, ssoil, veg, bgc, met, canopy)

    use physical_constants
    use cable_common_module   

    TYPE (soil_parameter_type), INTENT(IN)   :: soil
    TYPE (soil_snow_type), INTENT(IN)        :: ssoil
    TYPE (veg_parameter_type), INTENT(IN)    :: veg
    TYPE (bgc_pool_type), INTENT(IN)         :: bgc
    TYPE (met_type), INTENT(IN)              :: met 
    TYPE (canopy_type), INTENT(OUT)          :: canopy

    REAL(r_1), DIMENSION(mp)  :: den ! sib3  
    INTEGER(i_d)                             :: k
    REAL(r_1), DIMENSION(mp)  :: rswc   
    REAL(r_1), DIMENSION(mp)  :: sss    
    REAL(r_1), DIMENSION(mp)  :: e0rswc 
    REAL(r_1), DIMENSION(mp)  :: ftsoil 
    REAL(r_1), DIMENSION(mp)  :: ftsrs  
    REAL(r_1), PARAMETER                     :: t0 = -46.0
    REAL(r_1), DIMENSION(mp)  :: tref   
    REAL(r_1), DIMENSION(mp)  :: tsoil 
    REAL(r_1), DIMENSION(mstype)             :: rswch
    REAL(r_1), DIMENSION(mstype)             :: soilcf
    REAL(r_1), DIMENSION(mp)            :: avgtrs !root weighted mean soil temperature
    REAL(r_1), DIMENSION(mp)            :: avgwrs !root weighted mean soil moisture

    if (cable_user%DIAG_SOIL_RESP == 'off' .OR. &
         cable_user%DIAG_SOIL_RESP == 'OFF'    )  then
      avgwrs = sum(veg%froot * ssoil%wb,2)
      avgtrs = max(0.0,sum(veg%froot * ssoil%tgg,2)-tfrz)
      canopy%frs = veg%rs20 * min(1.0, max(0.0, min(&
           -0.0178+0.2883*avgwrs+5.0176*avgwrs*avgwrs-4.5128*avgwrs*avgwrs*avgwrs, &
           0.3320+22.6726*exp(-5.8184*avgwrs)))) &
           * min(1.0, max(0.0, min( 0.0104*(avgtrs**1.3053), 5.5956-0.1189*avgtrs)))
      canopy%frs = canopy%frs &
           * sum(spread(bgc%ratecs,1,mp) * bgc%csoil,2) &
           /(365.0*24.0*3600.0)   !convert 1/year to 1/second
      WHERE (ssoil%snowd > 1.)
        canopy%frs = canopy%frs / max(0.001,min(100.,ssoil%snowd))
      END WHERE

    else

      rswch = 0.16
      soilcf = 1.0

      den = max(0.07,soil%sfc - soil%swilt)
      rswc = MAX(0.0001, veg%froot(:,1)*(REAL(ssoil%wb(:,2),r_1) - soil%swilt))&
         & / den
      tsoil = veg%froot(:,1) * ssoil%tgg(:,2) - tfrz
    
      tref = MAX(0.,ssoil%tgg(:,ms) - (tfrz-.05) )
    
      DO k = 2,ms 
         rswc = rswc + MAX(0.0001, veg%froot(:,k) &
            & * (REAL(ssoil%wb(:,k),r_1) - soil%swilt)) / den
         tsoil = tsoil + veg%froot(:,k) * ssoil%tgg(:,k)
      ENDDO
      rswc = MIN(1.,rswc)
      tsoil = MAX(t0 + 2., tsoil)
      e0rswc = 52.4 + 285. * rswc
      ftsoil=min(0.0015,1./(tref - t0) - 1./(tsoil - t0))
      sss = MAX(-15.,MIN(1.,e0rswc * ftsoil))
      ftsrs=EXP(sss)
      canopy%frs = veg%vegcf * (144.0 / 44.0e6)  &
         & * soilcf(soil%isoilm) * MIN(1.,1.4 * MAX(.3,.0278 * tsoil + .5)) &
         & * ftsrs * rswc / (rswch(soil%isoilm) + rswc)
    endif
  END SUBROUTINE soilcarb

! rml 17/1/11 made plant respiration calculation into a subroutine
  SUBROUTINE plantcarb(veg, bgc, met, canopy)
   use physical_constants
    TYPE (veg_parameter_type), INTENT(IN)    :: veg
    TYPE (bgc_pool_type), INTENT(IN)         :: bgc
    TYPE (met_type), INTENT(IN)              :: met
    TYPE (canopy_type), INTENT(OUT)          :: canopy

    REAL(r_1), DIMENSION(mp) :: poolcoef1 ! non-leaf carbon turnover rate * non-leaf pool size
    REAL(r_1), DIMENSION(mp) :: poolcoef1w ! wood carbon turnover rate * wood pool size
    REAL(r_1), DIMENSION(mp) :: poolcoef1r ! root carbon turnover rate * root pool size
    REAL(r_1), DIMENSION(mp) :: tmp1,tmp2,tmp3


    poolcoef1=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) - &
         bgc%ratecp(1)*bgc%cplant(:,1))
    poolcoef1w=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -  &
         bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(3)*bgc%cplant(:,3))
    poolcoef1r=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -  &
         bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(2)*bgc%cplant(:,2))

    tmp1(:) = 3.22 - 0.046 * (met%tk(:)-tfrz)
    tmp2(:) = 0.1 * (met%tk(:)-tfrz-20.0)
    tmp3(:) = tmp1(:) ** tmp2(:)
    canopy%frp  = veg%rp20 * tmp3 * poolcoef1  / (365.0*24.0*3600.0)
    canopy%frpw = veg%rp20 * tmp3 * poolcoef1w / (365.0*24.0*3600.0)
    canopy%frpr = veg%rp20 * tmp3 * poolcoef1r / (365.0*24.0*3600.0)

    return
  END SUBROUTINE plantcarb



END MODULE carbon_module
