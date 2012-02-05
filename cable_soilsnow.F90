!
!soil_snow.f90
!
! Soil and snow routines source file for CABLE, CSIRO land surface model
!
! Science development by Eva Kowalczyk, CSIRO Marine and Atmospheric Research
! 
! Fortran-95 coding by Harvey Davies, Gab Abramowitz, Martin Dix and Bernard Pak
! Please report bugs to bernard.pak@csiro.au
!
! Major rewrite to modular format in Nov 2007 by Eva Kowalczyk
!
! This file contains the soil_snow_module only
! which has the following subroutines:
!   trimb,
!   smoisturev,
!   snowdensity,
!   snow_melting,
!   snow_accum,
!   surfbv,
!   snow_albedo,
!   stempv,
!   stempvsn,
!   snowcheck,
!   snowl_adjust,
!   soilfreeze,
!   remove_trans, and
!   soil_snow
!
MODULE soil_snow_module
  USE physical_constants ! from cable_variables.f90
  USE define_types       ! from cable_variables.f90
  IMPLICIT NONE
  PRIVATE
  REAL(r_1), PARAMETER :: cgsnow = 2090.0 ! specific heat capacity for snow
  REAL(r_1), PARAMETER :: csice = 2.100e3 ! specific heat capacity for ice
  REAL(r_1), PARAMETER :: cswat = 4.218e3 ! specific heat capacity for water
  REAL(r_1), PARAMETER :: hl = 2.5104e6   ! latent heat of evaporation
  REAL(r_1), PARAMETER :: hlf = 0.335e6   ! latent heat of fusion
  REAL(r_1), PARAMETER :: cp = 1004.64    ! specific heat capacity for air
  REAL(r_1), PARAMETER :: rhowat = 1000.0 ! density of water
!  REAL(r_1), PARAMETER :: snmin = 0.11    ! for 3-layer;
!  REAL(r_1), PARAMETER :: snmin =  100000.0 ! for 1-layer;
!  REAL(r_1), PARAMETER :: snmin =  20.0 ! for 3-layer in land-ice areas;
  REAL(r_1), PARAMETER :: snmin =  1.0 ! for 3-layer in land-ice areas;

  ! This module contains the following subroutines:
  PUBLIC soil_snow ! must be available outside this module
  PRIVATE trimb, smoisturev, snow_accum, stempv, stempvsn 
  PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust 
!  PRIVATE soilfreeze, remove_trans, snow_albedo
   PRIVATE soilfreeze, remove_trans
!  INTEGER(i_d), PARAMETER, PRIVATE  :: idjd = 9985
!  INTEGER(i_d), PARAMETER, PRIVATE  :: idjd = 6488
!  INTEGER(i_d), PARAMETER, PRIVATE  :: idjd = 10280 
  INTEGER(i_d), PARAMETER, PRIVATE  :: idjd = 2654 

CONTAINS

  !----------------------------------------------------------------------
  ! SUBROUTINE trimb
  !
  !      this routine solves the system
  !	   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kmax-1
  !	   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)	       for k=1
  !	   and	 a(k)*u(k-1)+b(k)*u(k)=rhs(k)	       for k=kmax
  !
  !	 the Thomas algorithm is used for solving sets of linear equation
  !	 rhs initially contains rhs; leaves with answer (jlm)
  !	 n.b. this one does not assume b = 1-a-c
  !
  SUBROUTINE trimb (a, b, c, rhs, kmax)
    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: a ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: b ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: c ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(INOUT)  :: rhs ! right hand side of eq
    INTEGER(i_d), INTENT(IN)                  :: kmax ! no. of discrete layers
    INTEGER(i_d)                              :: k ! do lloop counter
    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: e 
    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: g 
    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: temp 

    e(:,1) = c(:,1) / b(:,1)
    DO k = 2, kmax - 1
      temp(:,k) = 1. / (b(:,k) - a(:,k) * e(:,k-1) )
      e(:,k) = c(:,k) * temp(:,k)
    END DO
    g(:,1) = rhs(:,1) / b(:,1)
    DO k = 2, kmax - 1
      g(:,k) = (rhs(:,k) - a(:,k) * g(:,k-1) ) * temp(:,k)
    END DO

    ! do back substitution to give answer now
    rhs(:,kmax) = (rhs(:,kmax) - &
       & a(:,kmax) * g(:,kmax-1)) / (b(:,kmax) - a(:,kmax) * e(:,kmax-1))
    DO k = kmax - 1, 1, - 1
      rhs(:,k) = g(:,k) - e(:,k) * rhs(:,k + 1)
    END DO
  END SUBROUTINE trimb

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoisturev (fwtop,dels,ktau,ssoil,soil)
  !      Solves implicit soil moisture equation
  !      Science development by Eva Kowalczyk and John McGregor, CMAR
  !
  SUBROUTINE smoisturev (dels,ktau,ssoil,soil)
    REAL(r_1), INTENT(IN)                     :: dels    ! time step size (s)
    INTEGER(i_d), INTENT(IN)                  :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    INTEGER(i_d), PARAMETER                   :: ntest = 0 ! 2 for funny pre-set
                                                           ! for idjd
    ! nmeth selects the solution method
    INTEGER(i_d), PARAMETER                   :: nmeth = - 1 ! preferred method
    !                                  Values as follows:
    !                                   -1 for simple implicit D
    !                                    1 for fully implicit solution
    !                                    2 for simpler implicit
    !                                    3 for simple implicit D, explicit K 
    !                                    4 for simple implicit D, implicit K
    !                                    0 for simple implicit D, new jlm TVD K
    REAL(r_2), DIMENSION(mp,3*ms) :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,3*ms) :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,3*ms) :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp)      :: fact
    REAL(r_2), DIMENSION(mp)      :: fact2
    REAL(r_2), DIMENSION(mp)      :: fluxhi
    REAL(r_2), DIMENSION(mp)      :: fluxlo
    REAL(r_2), DIMENSION(mp)      :: hydss  ! hydraulic conductivity
                                            ! adjusted for ice
    INTEGER(i_d)                  :: k
    REAL(r_2), DIMENSION(mp)      :: phi
    REAL(r_2), DIMENSION(mp)      :: pwb
!    REAL(r_2), DIMENSION(:), ALLOCATABLE, SAVE :: pwb_min
    REAL(r_2), DIMENSION(mp)      :: rat
    REAL(r_2), DIMENSION(mp)      :: speed_k
    REAL(r_2), DIMENSION(mp)      :: ssatcurr_k
    REAL(r_1), DIMENSION(mp)      :: wblfmn
    REAL(r_1), DIMENSION(mp)      :: wblfmx
    REAL(r_2), DIMENSION(mp,ms+1) :: wbh
    REAL(r_2), DIMENSION(mp,ms+1) :: z1
    REAL(r_2), DIMENSION(mp,ms+1) :: z2
    REAL(r_2), DIMENSION(mp,ms+1) :: z3
    REAL(r_1), DIMENSION(mp,ms+1) :: z1mult
    REAL(r_2), DIMENSION(mp,0:ms) :: fluxh
    REAL(r_2), DIMENSION(mp,0:ms) :: delt
    REAL(r_2), DIMENSION(mp,0:ms) :: dtt
    REAL(r_2), DIMENSION(mp)      :: pwb_wbh
    REAL(r_2), DIMENSION(mp,ms)   :: ssatcurr
    REAL(r_1), DIMENSION(mp)      :: totwba ! diagnostic
    REAL(r_1), DIMENSION(mp)      :: totwbb
    REAL(r_1), DIMENSION(mp)      :: totwbc
    REAL(r_1), DIMENSION(mp)      :: totwblb
    REAL(r_1), DIMENSION(mp)      :: totwblc
    REAL(r_1), DIMENSION(mp)      :: wbficemx
    REAL(r_2), DIMENSION(mp)      :: wbh_k
    REAL(r_2), DIMENSION(mp)      :: wbl_k
    REAL(r_2), DIMENSION(mp)      :: wbl_kp
    REAL(r_2), DIMENSION(mp)      :: wh
    REAL(r_2), DIMENSION(mp)      :: z3_k
    LOGICAL :: is_open     ! Is file open?
    INTEGER :: u           ! I/O unit
    !
!    IF (ktau == 1) THEN
!      ! For some sorts of nested models, this might have been allocated on
!      ! another nest
!      IF ( ALLOCATED(pwb_min) ) DEALLOCATE(pwb_min)
!      ALLOCATE(pwb_min(mp))
!      ! Set working variable:
!      pwb_min = (soil%swilt / soil%ssat ) **soil%ibp2
!    END IF
    ! Block below for testing purposes only: - - - - - - - - - - - -
!        print *, 'wb   smoist1',(ssoil%wb(idjd,k) , k = 1, ms)
!        print *, 'wblf smoist1',(ssoil%wblf(idjd,k),k = 1, ms)
    IF (ntest > 0) THEN
      PRINT * , 'entering smoisturev fwtop1,i2bp3,swilt,sfc,ssat: ', &
        & ssoil%fwtop1(idjd), soil%i2bp3(idjd), soil%swilt(idjd), &
        & soil%sfc(idjd), soil%ssat(idjd)
      u = 97
      inquire (u, opened=is_open)
      IF (.NOT. is_open) THEN
        open (u, file='f97.txt', status='replace')
        write (u, *) 'ktau', ' fwtop', ' i2bp3', ' swilt', ' sfc', ' ssat'
      END IF
      write (u, *) ktau,ssoil%fwtop1,soil%i2bp3,soil%swilt,soil%sfc,soil%ssat
!      IF (ntest == 2) THEN ! just to test conservation
!        IF (ktau == 1) ssoil%wb(:,ms) = soil%swilt
!        ssoil%fwtop1 = 0.0
!        ssoil%fwtop2 = 0.0
!        ssoil%fwtop3 = 0.0
!      END IF
      WRITE (6, " ('wb   ', 6f8.3) ")    (ssoil%wb(idjd,k) , k = 1, ms)
      WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k) , k = 1, ms)
      totwba = 0.0
      DO k = 1, ms
        totwba = totwba + soil%zse(k) * REAL(ssoil%wb(:,k),r_1)
      END DO
    END IF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
    at = 0.0
    bt = 1.0
    ct = 0.0
    z1mult(:,1) = 0.0       ! corresponds to 2b+3
    z1mult(:,ms+1) = 0.0    ! corresponds to 2b+3
    z1(:,1) = 0.0           ! i.e. K(.5),    value at surface
    z1(:,ms+1) = 0.0        ! i.e. K(ms+.5), value at bottom
    ! nmeth: equation solution technique
    IF (nmeth <= 0) THEN
      ! jlm split TVD version
      ! all land points
      delt(:,0) = 0.0
      fluxh(:,0) = 0.0
      fluxh(:,ms) = 0.0
      DO k = 1, ms-1
        ! Calculate amount of liquid soil water:
        wbl_k = MAX(  0.01_r_2, ssoil%wb(:,k)   - ssoil%wbice(:,k)   )
        wbl_kp = MAX( 0.01_r_2, ssoil%wb(:,k+1) - ssoil%wbice(:,k+1) )
        ! Calculate difference in liq soil water b/w consecutive layers:
        delt(:,k) = wbl_kp - wbl_k
        ! especially to allow for isolated frozen layers, use min speed
        wh = MIN(wbl_k, wbl_kp)
        where(ssoil%wbice(:,k).gt.0.05.or.ssoil%wbice(:,k+1).gt.0.01) &
           wh = 0.9*wbl_k + 0.1*wbl_kp
        ! with 50% wbice, reduce hyds by 1.e-5
        ! Calculate hyd conductivity adjusted for ice:
        hydss = soil%hyds 
!        hydss = soil%hyds * (1.0 - MIN(2.0_r_2 * ssoil%wbice(:,k) &
!              & / MAX(0.01_r_2,  ssoil%wb(:,k) ), 0.99999_r_2) )
        speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
        ! update wb by TVD method
        rat = delt(:,k - 1) / (delt(:,k)+SIGN(REAL(1.0e-20,r_2), delt(:,k)))
        phi = MAX(0.0, MIN(1.0, 2.0 * rat), MIN(2., rat) ) ! 0 for -ve rat
!        phi = MAX(0._r_2, MIN(1._r_2, 2._r_2 * rat), &
!            & MIN(2._r_2, rat) ) ! 0 for -ve rat
        fluxhi = wh
        fluxlo = wbl_k
        ! Block below for testing purposes only:
        IF (ntest > 0) THEN
          PRINT * , 'in TVD for k= ', k
          PRINT * , 'wbl,wh,hydss ', wbl_k(idjd), wh(idjd), hydss(idjd)
          PRINT * , 'speeda,speedb,fluxhi,fluxlo,delt,rat,phi ', &
             & speed_k(idjd), 0.5*soil%zse(k)/dels, fluxhi(idjd), &
             & fluxlo(idjd), delt(idjd,k), rat(idjd), phi(idjd)
        END IF
        ! scale speed to grid lengths per dt & limit speed for stability
        ! 1. OK too for stability
        speed_k = MIN(speed_k, 0.5 * soil%zse(k) / dels)
        fluxh(:,k) = speed_k * (fluxlo + phi * (fluxhi - fluxlo) )
      END DO
      ! calculate drainage (this code replaces the code in the surfb)
      k = ms 
!      WHERE( ssoil%wb(:,ms) > soil%sfc(:))
!      WHERE( soil%albsoil(:) .lt. 0.25 )
!      WHERE( soil%albsoil(:) .lt. 0.30 )
!      WHERE( ssoil%wb(:,ms) > soil%sfc(:))
!      WHERE( ssoil%wb(:,ms) > max(soil%sfc(:),0.14) )
      WHERE( ssoil%wb(:,ms) > soil%sfc(:) )
!      WHERE( ssoil%wb(:,ms) > 0.5*(soil%sfc(:)+soil%ssat(:)) )
        wbl_k = MAX(0.001_r_2, ssoil%wb(:,ms) - ssoil%wbice(:,ms) )
        wbl_kp = MAX(0.001_r_2, soil%ssat(:) - ssoil%wbice(:,ms) )
        !wh = 0.9*wbl_k + 0.1*wbl_kp
        wh = MIN(wbl_k, wbl_kp)
        !where(ssoil%wbice(:,ms).gt. 0.05) wh = 0.9*wbl_k + 0.1*wbl_kp
        where(ssoil%wbice(:,ms).gt. 0.05) wh = 0.8*wbl_k + 0.2*wbl_kp
        ! Calculate hyd conductivity adjusted for ice:
        hydss = soil%hyds 
        !hydss = soil%hyds * ( 1. - MIN( 2.0_r_2 * ssoil%wbice(:,ms) &
        !      & / MAX(0.01_r_2, ssoil%wb(:,ms)), 0.99999_r_2 ) )
        speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
        !speed_k =  0.5*speed_k / (1. - min(0.5,10.*ssoil%wbice(:,ms)))
        fluxlo = wbl_k
        ! scale speed to grid lengths per dt & limit speed for stability
        speed_k = MIN(speed_k, 0.5 * soil%zse(ms) / dels)
        fluxh(:,ms) = MAX(0.0,speed_k * fluxlo )
       END WHERE
!      ELSEWHERE
!       WHERE( ssoil%wb(:,ms) > soil%sfc(:))
!        wbl_k = MAX(0.01_r_2, ssoil%wb(:,ms) - ssoil%wbice(:,ms) )
!        wbl_kp = MAX(0.01_r_2, soil%ssat(:) - ssoil%wbice(:,ms) )
!        wh = 0.9*wbl_k + 0.1*wbl_kp
!        ! Calculate hyd conductivity adjusted for ice:
!        hydss = soil%hyds * ( 1. - MIN( 2.0_r_2 * ssoil%wbice(:,ms) &
!              & / MAX(0.01_r_2, ssoil%wb(:,ms)), 0.99999_r_2 ) )
!        speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
!        fluxlo = wbl_k
!        ! scale speed to grid lengths per dt & limit speed for stability
!!        speed_k = MIN(speed_k, 0.5 * soil%zse(ms) / dels)
!        speed_k = MIN(0.5*speed_k, 0.5 * soil%zse(ms) / dels)
!        fluxh(:,ms) = MAX(0.0,speed_k * fluxlo )
!       END WHERE
!      END WHERE

      ! update wb by TVD method
      DO k = ms, 1, -1
        IF (nmeth == -1) THEN ! each new wb constrained by ssat
          fluxh(:,k-1) = MIN(fluxh (:,k-1), (soil%ssat &
                 - ssoil%wb(:,k) ) * soil%zse(k) / dels + fluxh(:,k) )
        END IF
        ! fluxh (:,ms) is drainage
        ssoil%wb(:,k) = ssoil%wb(:,k) + dels * (fluxh(:,k-1) - fluxh(:,k)) &
                                         & / soil%zse(k)
        ! re-calculate wblf
        ssatcurr_k = soil%ssat - ssoil%wbice(:,k)
        dtt(:,k) = dels / (soil%zse(k) * ssatcurr_k)
        ! this defn of wblf has different meaning from previous one in surfbv
        ! N.B. are imposing wbice<wb, so wblf <1
        ssoil%wblf(:,k) = (ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr_k
      END DO
      ssoil%rnof2 = dels * REAL(fluxh(:,ms),r_1) * 1000.0

      ! wbh_k represents wblf(k-.5)
      DO k = 2, ms
        ssatcurr_k = REAL(soil%ssat,r_2) - ssoil%wbice(:,k)
        wbh_k = ( soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1) &
              & * ssoil%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
        ! i.e. wbh**(bch+1)

        fact = wbh_k** (soil%ibp2 - 1)
        ! with 50% wbice, reduce hbsh by 1.e-5
!        pwb_wbh = (soil%hsbh * (1. - MIN(2. * MIN(0.99_r_2, MAX( &
!        pwb_wbh = (soil%hsbh * (1. - MIN(2. * MIN(0.2_r_2, MAX( &
!        pwb_wbh = (soil%hsbh * (1. - MIN(2. * MIN(0.1_r_2, MAX( &
        pwb_wbh = (soil%hsbh * (1. - MIN(2. * MIN(0.1_r_2, MAX( &
                & ssoil%wbice(:,k-1) / MAX(0.01_r_2, ssoil%wb(:,k-1)), &
                & ssoil%wbice(:,k)   / MAX(0.01_r_2, ssoil%wb(:,k)  ) )) &
                & , 0.1_r_2) )) &
                & * MAX(soil%pwb_min, wbh_k * fact)
!                & , 0.2_r_2) )) &
!                & , 0.99999_r_2) )) &
        ! moisture diffusivity (D) is  wbh*pwb; hsbh includes b
        ! i.e. D(k-.5)/soil%zshh(k)
        z3_k = pwb_wbh / soil%zshh (k)
        ! PRINT * , 'z3_k ', z3_k
        ! where dtt=dels/(soil%zse(k)*ssatcurr_k)
        at (:,k) = - dtt(:,k) * z3_k
        ct (:,k-1) = - dtt(:,k-1) * z3_k
      END DO
      bt = 1. - at - ct
      ! Block below for testing purposes only:
!        print *, 'wb   smoist2',(ssoil%wb(idjd,k) , k = 1, ms)
!        print *, 'wblf smoist2',(ssoil%wblf(idjd,k) , k = 1, ms)
      IF (ntest > 0) THEN
!        PRINT * , 'midway through nmeth<=0'
!        PRINT * , 'fluxh ', (fluxh(idjd,k) , k = 1, ms)
        WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k) , k = 1, ms)
        WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k) , k = 1, ms)
        totwbb = 0.
        totwblb = 0.
        DO k = 1, ms
          totwbb = totwbb + soil%zse(k) * REAL(ssoil%wb(:,k),r_1) ! diagnostic
          totwblb=totwblb + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1) ! diagnostic
        END DO
!        PRINT * , 'nmeth, b+2, 2b+3: ',nmeth, soil%ibp2(idjd), soil%i2bp3(idjd)
        WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k) , k = 1, ms)
        WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k) , k = 1, ms)
        WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k) , k = 1, ms)
!        PRINT * , 'zse ', soil%zse
!        PRINT * , 'zshh ', soil%zshh
!        PRINT * , 'dtt ', (dtt(idjd,k) , k = 1, ms)
!        PRINT * , 'at ', (at(idjd,k) , k = 1, ms)
!        PRINT * , 'bt ', (bt(idjd,k) , k = 1, ms)
!        PRINT * , 'ct ', (ct(idjd,k) , k = 1, ms)
      END IF
      ssoil%wblf(:,1) = ssoil%wblf(:,1) + dtt(:,1) * ssoil%fwtop1 / rhowat
      ssoil%wblf(:,2) = ssoil%wblf(:,2) + dtt(:,2) * ssoil%fwtop2 / rhowat
      ssoil%wblf(:,3) = ssoil%wblf(:,3) + dtt(:,3) * ssoil%fwtop3 / rhowat
    END IF

    IF (nmeth > 0) THEN
      wbficemx = 0.0
      DO k = 1, ms
        ssatcurr(:,k) = REAL(soil%ssat,r_2) - ssoil%wbice(:,k)
        ! this defn of wblf has different meaning from previous one in surfbv
        ! N.B. are imposing wbice<wb, so wblf <1
        ssoil%wblf(:,k) = (ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr(:,k)
        ssoil%wbfice(:,k) = REAL(ssoil%wbice(:,k),r_1) / soil%ssat
        wbficemx = MAX(wbficemx, ssoil%wbfice(:,k) )
        dtt(:,k) = dels / (soil%zse(k) * ssatcurr(:,k) )
      END DO
      IF (nmeth == 1) THEN ! full implicit method
        DO k = 2, ms
          ! wbh(k)=MIN(1.,ww(k)*wblf(:,k-1)+(1.-ww(k))*wblf(:,k))
          ! jlm: this is same as:
          wbh(:,k) = (soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1) &
                   & * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k-1) )
          fact = wbh(:,k) ** (soil%ibp2 - 1) ! i.e. wbh**(bch+1)
          fact2 = fact * fact
          pwb = soil%hsbh * fact
          ! moisture diffusivity (D) is  wbh*pwb
          ! other term (K) is wbh*soil%hyds*fact2
          z1(:,k) = wbh(:,k) * ( (soil%i2bp3 - 1) * soil%hyds * fact2 - &
                  & soil%ibp2 * pwb * (ssoil%wblf(:,k) - ssoil%wblf(:,k-1) ) &
                  & / soil%zshh (k) )
          z2(:,k) = - soil%i2bp3 * soil%hyds * fact2 + soil%ibp2 * pwb &
                  & * (ssoil%wblf(:,k) - ssoil%wblf(:,k-1) ) / soil%zshh (k)
          z3(:,k) = pwb * wbh(:,k) / soil%zshh (k)
          at(:,k) = dtt(:,k) * (z2(:,k) * 0.5 * soil%zse(k) / soil%zshh (k) &
                  & - z3(:,k) )
        END DO
        DO k = 1, ms - 1
          ct(:,k) = dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k) &
                  & / soil%zshh (k+1) - z3(:,k+1) )
          bt(:,k) = 1.0 + dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k+1) &
                  & / soil%zshh (k+1) + z2(:,k) * 0.5 * soil%zse(MAX(k-1,1)) &
                  & / soil%zshh (k) + z3(:,k+1) + z3(:,k) )
        END DO
        bt(:,ms) = 1.0 + dtt(:,ms) * (z2(:,ms) * 0.5 * soil%zse(ms) &
                  & / soil%zshh (ms) + z3(:,ms) )
        DO k = 1, ms
          ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) * (z1(:,k+1) - z1(:,k) )
        END DO
      END IF
      IF (nmeth >= 2) THEN ! part implicit method
        DO k = 2, ms
          z1mult(:,k) = soil%i2bp3 ! corresponds to 2b+3
        END DO
        DO k = 2, ms ! wbh(k) represents wblf(k-.5)
          wbh(:,k) = (soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1) &
                   & * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k-1) )
          fact = wbh(:,k) ** (soil%ibp2 - 1) ! i.e. wbh**(bch+1)
          IF (nmeth == 2) pwb_wbh = soil%hsbh * wbh(:,k) * fact
          IF (nmeth >= 3) pwb_wbh = soil%hsbh * MAX(soil%pwb_min,wbh(:,k)*fact)
          fact2 = fact * fact
          ! moisture diffusivity (D) is  wbh*pwb
          ! other term (K) is wbh*soil%hyds*fact2
          z1(:,k) = soil%hyds * fact2 !  i.e. K(k-.5)/wbh(:,k)
          z3(:,k) = pwb_wbh / soil%zshh(k) !  i.e. D(k-.5)/soil%zshh(k)
          at(:,k) = - dtt(:,k) * z3(:,k)
          ct(:,k-1) = - dtt(:,k-1) * z3(:,k)
        END DO
        bt = 1. - at - ct
        IF (nmeth == 4) THEN ! for simple implicit D, implicit K
          bt(:,1) = bt(:,1) + dtt(:,1) * z1mult(:,1+1) * z1(:,1+1) &
                  & * soil%zse(1+1) / (soil%zse(1) + soil%zse(1+1) )
          DO k = 2, ms
            at(:,k)   = at(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k) &
                      & * soil%zse(k) / (soil%zse(k) + soil%zse(k-1) )
            ct(:,k-1) = ct(:,k-1) + dtt(:,k-1) * z1mult(:,k) * z1(:,k) &
                      & * soil%zse(k-1) / (soil%zse(k) + soil%zse(k-1) )
            bt(:,k) = bt(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k) &
                      & * soil%zse(k-1) / (soil%zse(k) + soil%zse(k-1) ) &
                      & + dtt(:,k) * z1mult(:,k+1) * z1(:,k+1) &
                      & * soil%zse(k+1) / (soil%zse(k) + soil%zse(k+1) )
          END DO
        ! (nmeth == 4)
        END IF
        DO k = 2, ms
          ! i.e. now K(k-.5)
          z1(:,k) = wbh(:,k) * z1(:,k)
        END DO

        ! the following top & bottom b.c.'s will preserve a uniform column
        !     z1(1) =z1(2)   ! simple dk/dz=0
        !     z1(ms+1)=z1(ms) ! simple dk/dz=0
        ! N.B. z1 are here +ve
        z1(:,1) = MIN(z1(:,2), z1(:,ms) )
        z1(:,ms + 1) = z1(:,1)
        ! no gravit. term if too much ice 11/12/00
        DO k = 1, ms
          IF (nmeth == 4) THEN
            WHERE (wbficemx < 0.75)
              ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) &
                              & * ( (z1mult(:,k+1) - 1.0) * z1(:,k+1) &
                              & - (z1mult(:,k) - 1.0) * z1(:,k) )
            END WHERE
          ELSE
            WHERE (wbficemx < 0.75)
              ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) &
                              & * (z1(:,k) - z1(:,k+1) )
            END WHERE
          END IF
        END DO
      END IF
      IF (ntest > 0) THEN
        wblfmx = MAXVAL(REAL(ssoil%wblf,r_1), 2)
        wblfmn = MINVAL(REAL(ssoil%wblf,r_1), 2)
      END IF
      ! Block below for testing purposes only:
      IF (ntest > 0) THEN
        totwbb = 0.0
        totwblb = 0.0
        DO k = 1, ms
          ! diagnostic
          totwbb = totwbb + soil%zse(k) * REAL(ssoil%wb(:,k),r_1)
          ! diagnostic
          totwblb = totwblb + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1)
        END DO
!        PRINT * , 'nmeth, b+2, 2b+3: ',nmeth, soil%ibp2(idjd), soil%i2bp3(idjd)
        WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k), k = 1, ms)
        WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k), k = 1, ms)
        WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k), k = 1, ms)
        WRITE (6, " ('wbh  ', 7f8.3) ") wbh(idjd,:)
        WRITE (6, " ('ssatcurr', 6f8.3) ") ssatcurr(idjd,:)
!        PRINT * , 'pwb_wbh,soil%pwb_min* for ms ', pwb_wbh(idjd), &
!                                 & soil%hsbh(idjd) * soil%pwb_min(idjd)
!        PRINT * , 'wblfmx,wblfmn', wblfmx(idjd), wblfmn(idjd)
!        PRINT * , 'zse ', soil%zse
!        PRINT * , 'zshh ', soil%zshh
!        PRINT * , 'at ', (at(idjd,k), k = 1, ms)
!        PRINT * , 'bt ', (bt(idjd,k), k = 1, ms)
!        PRINT * , 'ct ', (ct(idjd,k), k = 1, ms)
      END IF
      IF (nmeth == 3) THEN
        ! artificial fix applied here for safety (explicit nmeth only)
        DO k = 1, ms
          ssoil%wblf(:,k) = MAX(0.0_r_2, MIN(ssoil%wblf(:,k), 1.0_r_2) )
        END DO
        ! (nmeth == 3)
      END IF
      ssoil%wblf(:,1) = ssoil%wblf(:,1) + dtt(:,1) * ssoil%fwtop1 / rhowat
      ssoil%wblf(:,2) = ssoil%wblf(:,2) + dtt(:,2) * ssoil%fwtop2 / rhowat
      ssoil%wblf(:,3) = ssoil%wblf(:,3) + dtt(:,3) * ssoil%fwtop3 / rhowat
    END IF  ! IF (nmeth > 0)

    CALL trimb(at, bt, ct, ssoil%wblf, ms)
    DO k = 1, ms
      ssatcurr(:,k) = soil%ssat - ssoil%wbice(:,k)
      ssoil%wb(:,k) = ssoil%wblf(:,k) * ssatcurr(:,k) + ssoil%wbice(:,k)
      ssoil%wbice(:,k) = MIN(ssoil%wbice(:,k), 0.85 * ssoil%wb(:,k) )
    END DO
    ! Block below for testing purposes only:
!      print *, 'wb   3', (ssoil%wb(idjd,k), k = 1, ms)
!      print *, 'wbice3', (ssoil%wbice(idjd,k), k = 1, ms)
!      print *, 'wblf 3', (ssoil%wblf(idjd,k), k = 1, ms)
    IF (ntest > 0) THEN
!      PRINT * , 'at end of smoisturev,fwtop ', ssoil%fwtop(idjd)
!      PRINT * , 'tgg ', (ssoil%tgg(idjd,k), k = 1, ms)
      WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k), k = 1, ms)
      WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k), k = 1, ms)
      WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k), k = 1, ms)
      totwbc = 0.
      totwblc = 0.
      DO k = 1, ms
        totwbc = totwbc + soil%zse(k) * REAL(ssoil%wb(:,k),r_1)
        totwblc = totwblc + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1)
      END DO
!      PRINT * , 'totwba,totwbb,totwbc ', totwba(idjd), totwbb(idjd),totwbc(idjd)
!      PRINT * , 'totwblb,totwblc ', totwblb(idjd), totwblc(idjd)
    END IF
  END SUBROUTINE smoisturev


  !-------------------------------------------------------------------------
  SUBROUTINE snowdensity (dels, ssoil, soil)
    REAL(r_1), INTENT(IN)   :: dels   ! integration time step (s)
    TYPE(soil_snow_type),      INTENT(INOUT) :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    WHERE (ssoil%snowd > 0.1 .and. ssoil%isflag == 0)
      ssoil%ssdn(:,1) = MIN(750.0,MAX(120.0, ssoil%ssdn(:,1) + dels &
          & * ssoil%ssdn(:,1) * 3.1e-6 * EXP( -0.03 &
          & * (tfrz - MIN(tfrz, ssoil%tgg(:,1) )) &
          & - merge(0.046, 0.0, ssoil%ssdn(:,1) >= 150.0) &
          & * (ssoil%ssdn(:,1) - 150.0) ) ))
      ssoil%ssdn(:,1) = MIN(750.0,ssoil%ssdn(:,1) + dels * 9.806 &
          & * ssoil%ssdn(:,1) * 0.75 * ssoil%snowd &
          & / (3.0e7 * EXP(0.021 * ssoil%ssdn(:,1) + 0.081 &
          & * (tfrz - MIN(tfrz, ssoil%tgg(:,1))))))
      WHERE (soil%isoilm /= 9) ssoil%ssdn(:,1) = MIN(450.0,ssoil%ssdn(:,1))
      ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
                                                        & + 0.074, 2.51) )
!                                                        & + 0.074, 1.0) )
      ssoil%sconds(:,2) = ssoil%sconds(:,1) 
      ssoil%sconds(:,3) = ssoil%sconds(:,1) 
      ssoil%ssdnn = ssoil%ssdn(:,1)
      ssoil%ssdn(:,2) = ssoil%ssdn(:,1)
      ssoil%ssdn(:,3) = ssoil%ssdn(:,1)
    END WHERE
    WHERE (ssoil%isflag == 1)
      ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * ssoil%ssdn(:,1) * 3.1e-6 &
          & * EXP( -0.03 * (tfrz - MIN(tfrz, ssoil%tggsn(:,1))) &
          & - merge(0.046, 0.0, ssoil%ssdn(:,1) >= 150.0) &
          & * (ssoil%ssdn(:,1) - 150.0) )
      ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * ssoil%ssdn(:,2) * 3.1e-6 &
          & * EXP( -0.03 * (tfrz - MIN(tfrz, ssoil%tggsn(:,2))) &
          & - merge(0.046, 0.0, ssoil%ssdn(:,2) >= 150.0) &
          & * (ssoil%ssdn(:,2) - 150.0) )
      ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * ssoil%ssdn(:,3) * 3.1e-6 &
          & * EXP( -0.03 * (tfrz - MIN(tfrz, ssoil%tggsn(:,3))) &
          & - merge(0.046, 0.0, ssoil%ssdn(:,3) >= 150.0) &
          & * (ssoil%ssdn(:,3) - 150.0) )
      ssoil%ssdn(:,1) = min( ssoil%ssdn(:,1) + dels * 9.806 * ssoil%ssdn(:,1) &
          & * ssoil%t_snwlr(:)*ssoil%ssdn(:,1) &
          & / (3.0e7 * EXP(.021 * ssoil%ssdn(:,1) + 0.081 &
          & * (tfrz - MIN(tfrz, ssoil%tggsn(:,1))))) , 750.)
      ssoil%ssdn(:,2) = min( ssoil%ssdn(:,2) + dels * 9.806 * ssoil%ssdn(:,2) &
          & * (ssoil%t_snwlr * ssoil%ssdn(:,1) + 0.5 * ssoil%smass(:,2) ) &
          & / (3.0e7 * EXP(.021 * ssoil%ssdn(:,2) + 0.081 &
          & * (tfrz - MIN(tfrz, ssoil%tggsn(:,2))))) , 750.)
      ssoil%ssdn(:,3) = min( ssoil%ssdn(:,3) + dels * 9.806 * ssoil%ssdn(:,3) &
          & * (ssoil%t_snwlr*ssoil%ssdn(:,1) + ssoil%smass(:,2) + 0.5*ssoil%smass(:,3)) &
          & / (3.0e7 * EXP(.021 * ssoil%ssdn(:,3) + 0.081 &
          & * (tfrz - MIN(tfrz, ssoil%tggsn(:,3))))) , 750.)
      WHERE (soil%isoilm /= 9) ssoil%ssdn(:,1) = MIN(450.0,ssoil%ssdn(:,1))
      WHERE (soil%isoilm /= 9) ssoil%ssdn(:,2) = MIN(450.0,ssoil%ssdn(:,2))
      WHERE (soil%isoilm /= 9) ssoil%ssdn(:,3) = MIN(450.0,ssoil%ssdn(:,3))
      ssoil%sdepth(:,1) =  ssoil%smass(:,1) / ssoil%ssdn(:,1) 
      ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
      ssoil%sdepth(:,3) =  ssoil%smass(:,3) / ssoil%ssdn(:,3) 
      ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%smass(:,1) + ssoil%ssdn(:,2) &
             & * ssoil%smass(:,2) + ssoil%ssdn(:,3) * ssoil%smass(:,3) ) &
             & / ssoil%snowd
      ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
                                                        & + 0.074, 2.51) )
      ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2) ** 2 &
                                                        & + 0.074, 2.51) )
      ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3) ** 2 &
                                                        & + 0.074, 2.51) )
    END WHERE
!    print *,'ssdn & scnd',ssoil%ssdn
!    print *,'ssdn & scnda',ssoil%ssdnn
!    print *,'ssdn & scnd1',ssoil%snowd
!    print *,'ssdn & scnd2',ssoil%sconds
!    print *,'ssdn & scnd3',ssoil%ssdn,ssoil%ssdnn,ssoil%snowd,ssoil%sconds
!    print *,'ssdn & scnd',mp
!    print *,'snowdensity',ssoil%ssdnn
  END SUBROUTINE snowdensity

  !-------------------------------------------------------------------------
  SUBROUTINE snow_melting (dels, snowmlt, ktau, ssoil, soil )
    REAL(r_1), INTENT(IN)                 :: dels   ! integration time step (s)
    REAL(r_1), DIMENSION(mp), INTENT(OUT) :: snowmlt ! snow melt   
    INTEGER(i_d), INTENT(IN)              :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT)   :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
    INTEGER(i_d)                 :: k
    REAL(r_1), DIMENSION(mp)     :: osm
    REAL(r_1), DIMENSION(mp)     :: sgamm
    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
    REAL(r_1), DIMENSION(mp)     :: snowflx
   
    snowmlt= 0.0
    smelt1 = 0.0

!    ssoil%dtmlt = 0.0
    WHERE (ssoil%snowd > 0.0 .and. ssoil%isflag == 0 &
                       & .and. ssoil%tgg(:,1) >= tfrz )
      ! snow covered land
      ! following done in sflux  via  ga= ... +cls*egg + ...
      ! ** land,snow,melting
      snowflx = (ssoil%tgg(:,1) - tfrz) * ssoil%gammzz(:,1)
!      ssoil%dtmlt(:,1) =ssoil%dtmlt(:,1) + max(0.,ssoil%tgg(:,1) - tfrz)
      ! prevent snow depth going negative
      snowmlt = MIN(snowflx / hlf, ssoil%snowd )
      ssoil%dtmlt(:,1) = ssoil%dtmlt(:,1) + snowmlt * hlf / ssoil%gammzz(:,1)
!      ssoil%qfsrf = ssoil%qfsrf + snowmlt*hlf/dels
      ssoil%snowd = ssoil%snowd - snowmlt
      ssoil%tgg(:,1) = ssoil%tgg(:,1) - snowmlt * hlf / ssoil%gammzz(:,1)
    END WHERE

    smelt1(:,0) = 0.0
    DO k = 1, 3
      WHERE (ssoil%snowd > 0.0 .and. ssoil%isflag > 0)
        sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
        ! snow melt refreezing
!        snowflx = smelt1(:,k-1) * hlf / dels
!!        ssoil%tggsn(:,k) = ssoil%tggsn(:,k) + snowflx * dels / sgamm
!        ssoil%tggsn(:,k) = ssoil%tggsn(:,k) + ( snowflx * dels +  &
!                         & smelt1(:,k-1)*cswat*(tfrz-ssoil%tggsn(:,k)) ) / &
!                         & ( sgamm + cswat*smelt1(:,k-1) )
!        ! increase density due to snowmelt
!        osm = ssoil%smass(:,k)
!        ssoil%smass(:,k) = ssoil%smass(:,k) + smelt1(:,k-1)
!        ssoil%ssdn(:,k) = MAX(120.0,MIN(ssoil%ssdn(:,k) * osm/ssoil%smass(:,k) &
!                        & + rhowat*(1.0-osm/ssoil%smass(:,k)), 750.0))
!        WHERE (soil%isoilm /= 9) ssoil%ssdn(:,k) = MIN(450.0,ssoil%ssdn(:,k))
!        ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
!        sgamm = ssoil%smass(:,k) * cgsnow
!        smelt1(:,k-1) = 0.0
        smelt1(:,k) = 0.0
!        ! snow melting
!        ssoil%dtmlt(:,k) = 0.0
        WHERE (ssoil%tggsn(:,k) > tfrz)
          snowflx = (ssoil%tggsn(:,k) - tfrz) * sgamm
!          ssoil%dtmlt(:,k) = ssoil%dtmlt(:,k)+ max(0.,ssoil%tggsn(:,k) - tfrz)
          smelt1(:,k) = MIN(snowflx / hlf, 0.6 * ssoil%smass(:,k) )
          ssoil%dtmlt(:,k) = ssoil%dtmlt(:,k) + smelt1(:,k) * hlf / sgamm
          osm = ssoil%smass(:,k)
          ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
          ssoil%tggsn(:,k) = ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm
!          ssoil%tggsn(:,k) = MIN(ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm, &
!                               & tfrz)
          ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
        END WHERE
      END WHERE
    END DO
    WHERE (ssoil%snowd > 0.0 .and. ssoil%isflag > 0)
      snowmlt = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
!      ssoil%qfsrf = ssoil%qfsrf + snowmlt*hlf/dels
      ssoil%snowd = ssoil%snowd - snowmlt
    END WHERE
!    print *,'snowmelt routine',ssoil%snowd,snowmlt,smelt1(1,:)
!    print *,'snowmelt routine',ssoil%tggsn(1,:)
  END SUBROUTINE snow_melting


  !-------------------------------------------------------------------------
  SUBROUTINE snow_accum (dels,  ktau, canopy, met, ssoil, soil)
    REAL(r_1), INTENT(IN)                    :: dels   ! integration time step (s)
    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
    TYPE(canopy_type), INTENT(INOUT)         :: canopy ! vegetation variables
    TYPE(met_type), INTENT(INOUT)            :: met   ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    INTEGER(i_d), PARAMETER  :: ntest = 0 ! for snow diag prints
!    REAL(r_1), DIMENSION(mp) :: evapsn   ! now put in soil_snow_type EAK aug08
    INTEGER(i_d)             :: k
    REAL(r_1), DIMENSION(mp) :: osm
    REAL(r_1), DIMENSION(mp) :: sgamm
!    REAL(r_1), DIMENSION(mp) :: snowmlt
    REAL(r_1), DIMENSION(mp) :: xxx
    !
    IF (ntest > 0) THEN
      PRINT * , 'surfb1'
      PRINT * , 'entering snow_accum', canopy%precis(idjd)
      PRINT * , 'osnowd,snowd,isflag', ssoil%osnowd(idjd), &
                  & ssoil%snowd(idjd), ssoil%isflag(idjd)
      PRINT * , 'tggsn ', (ssoil%tggsn(idjd,k), k = 1, 3)
      PRINT * , 'tgg ', (ssoil%tgg(idjd,k), k = 1, ms)
      PRINT * , 'wb ', (ssoil%wb(idjd,k), k = 1, ms)
      PRINT * , 'wbice ', (ssoil%wbice(idjd,k), k = 1, ms)
      PRINT * , 'gammzz ', (ssoil%gammzz(idjd,k), k = 1, ms)
      PRINT * , 'albnew ', (ssoil%albsoilsn(idjd,k), k = 1, 2)
    END IF

!    snowmlt =0.0

!    WHERE (canopy%precis > 0.)  ! precis is both liquid and snow
    WHERE (canopy%precis > 0.0 .and. ssoil%isflag == 0)
      ssoil%snowd = MAX(ssoil%snowd + met%precip_sn, 0.0) ! accumulate solid part
      canopy%precis = canopy%precis - met%precip_sn
      ssoil%ssdn(:,1) = MAX(120.0, ssoil%ssdn(:,1) &
                      & * ssoil%osnowd / MAX(0.01, ssoil%snowd) &
                      & + 120.0 * met%precip_sn / MAX(0.01, ssoil%snowd))
      ssoil%ssdnn = ssoil%ssdn(:,1)
!      WHERE (canopy%precis > 0.0 .and. ssoil%tgg(:,1) < tfrz)
!        ssoil%snowd = MAX(ssoil%snowd + canopy%precis, 0.0)
!!        ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis * hlf &
!!                       & / REAL(ssoil%gammzz(:,1),r_1) + &
!        ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis *  &
!              & ( hlf + cswat*(tfrz-ssoil%tgg(:,1)) ) /  &
!              & ( REAL(ssoil%gammzz(:,1),r_1) + cswat*canopy%precis )  
!        ! change density due to water being added 
!        ssoil%ssdn(:,1) = MIN(750.0, MAX(120.0, ssoil%ssdn(:,1) &
!                        & * ssoil%osnowd / MAX(0.01, ssoil%snowd) &
!                        & + rhowat * canopy%precis / MAX(0.01, ssoil%snowd)))
!        WHERE (soil%isoilm /= 9) ssoil%ssdn(:,1) = MIN(450.0,ssoil%ssdn(:,1))
!!        ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis * hlf  &
!!                         / (ssoil%gammzz(:,1) + cswat * canopy%precis)
!!        ssoil%qasrf = ssoil%qasrf + canopy%precis * hlf/ dels
!        canopy%precis = 0.0
!        ssoil%ssdnn = ssoil%ssdn(:,1)
!      END WHERE
    END WHERE ! (canopy%precis > 0. .and. ssoil%isflag == 0) 
!    PRINT * , 'surfb2',ssoil%snowd,ssoil%osnowd
    WHERE (canopy%precis > 0.0 .and.  ssoil%isflag > 0)
      ! add solid precip
!      ssoil%qasrf = ssoil%qasrf + canopy%precis * hlf/dels
      ssoil%snowd = MAX(ssoil%snowd + met%precip_sn, 0.0)
      canopy%precis = canopy%precis - met%precip_sn  ! remaining liquid precip
      ! update top snow layer with fresh snow
      osm = ssoil%smass(:,1)
      ssoil%smass(:,1) = ssoil%smass(:,1) + met%precip_sn
      ssoil%ssdn(:,1) = MAX(120.0,ssoil%ssdn(:,1)*osm/ssoil%smass(:,1) &
                      & + 120.0 * met%precip_sn/ssoil%smass(:,1))
      ssoil%sdepth(:,1) = MAX(0.02,ssoil%smass(:,1) / ssoil%ssdn(:,1))
      ! add liquid precip
      !WHERE (canopy%precis > 0.0)
      !  ssoil%snowd = MAX(ssoil%snowd + canopy%precis, 0.0)
      !  sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
      !  osm = ssoil%smass(:,1)
!     !   ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + canopy%precis * hlf &
!     !                    & * osm / (sgamm * ssoil%osnowd )
!!    !                    * ssoil%smass(:,1) / (sgamm * ssoil%osnowd )
      !  ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + canopy%precis *  &
      !                   & ( hlf + cswat*(tfrz-ssoil%tggsn(:,1)) )* osm/ssoil%osnowd  &
      !                   & / (sgamm + cswat*canopy%precis*osm/ssoil%osnowd)
      !  ssoil%smass(:,1) = ssoil%smass(:,1) + canopy%precis &
      !                   & * osm/ssoil%osnowd
!     !                   * ssoil%smass(:,1)/ssoil%osnowd
      !  ssoil%ssdn(:,1) = MAX(120.0,MIN(ssoil%ssdn(:,1) * osm/ssoil%smass(:,1) &
      !                  & +  rhowat*(1.0-osm/ssoil%smass(:,1)), 750.0))
      !  WHERE (soil%isoilm /= 9) ssoil%ssdn(:,1) = MIN(450.0,ssoil%ssdn(:,1))
      !  ssoil%sdepth(:,1) = ssoil%smass(:,1)/ssoil%ssdn(:,1)
!
!        sgamm = ssoil%ssdn(:,2) * cgsnow * ssoil%sdepth(:,2)
!        osm = ssoil%smass(:,2)
!!        ssoil%tggsn(:,2) = ssoil%tggsn(:,2) + canopy%precis * hlf &
!!                         & * osm / (sgamm * ssoil%osnowd )
!!!                        * ssoil%smass(:,2) / (sgamm * ssoil%osnowd )
!        ssoil%tggsn(:,2) = ssoil%tggsn(:,2) + canopy%precis * &
!                         & ( hlf + cswat*(tfrz-ssoil%tggsn(:,2)) )* osm/ssoil%osnowd  &
!                         & / (sgamm + cswat*canopy%precis*osm/ssoil%osnowd)
!        ssoil%smass(:,2) = ssoil%smass(:,2) + canopy%precis &
!                         & * osm/ssoil%osnowd
!!                        ssoil%smass(:,2)/ssoil%osnowd
!        ssoil%ssdn(:,2) = MAX(120.0,MIN(ssoil%ssdn(:,2) * osm/ssoil%smass(:,2) &
!                        & + rhowat*(1.0-osm/ssoil%smass(:,2)), 750.0))
!        WHERE (soil%isoilm /= 9) ssoil%ssdn(:,2) = MIN(450.0,ssoil%ssdn(:,2))
!        ssoil%sdepth(:,2) = ssoil%smass(:,2) / ssoil%ssdn(:,2)
!
!        sgamm = ssoil%ssdn(:,3) * cgsnow * ssoil%sdepth(:,3)
!        osm = ssoil%smass(:,3)
!!        ssoil%tggsn(:,3) = ssoil%tggsn(:,3) + canopy%precis * hlf &
!!                         & * osm / (sgamm * ssoil%osnowd )
!        ssoil%tggsn(:,3) = ssoil%tggsn(:,3) + canopy%precis *  &
!                        & ( hlf + cswat*(tfrz-ssoil%tggsn(:,3)) ) * osm/ssoil%osnowd  &
!                        & / (sgamm + cswat*canopy%precis*osm/ssoil%osnowd)
!        ssoil%smass(:,3) = ssoil%smass(:,3) + canopy%precis &
!                         & * osm/ssoil%osnowd
!!                        * ssoil%smass(:,3)/ssoil%osnowd
!        ssoil%ssdn(:,3) = MAX(120.0,MIN(ssoil%ssdn(:,3) * osm/ssoil%smass(:,3) &
!                        & + rhowat*(1.0-osm/ssoil%smass(:,3)), 750.0))
!        WHERE (soil%isoilm /= 9)  ssoil%ssdn(:,3) = MIN(450.0,ssoil%ssdn(:,3))
!        ssoil%sdepth(:,3) = ssoil%smass(:,3)/ssoil%ssdn(:,3)
!
!        canopy%precis = 0.0
!
!      END WHERE
    END WHERE

!    WHERE (ssoil%snowd < 0.1 .and. canopy%fess .gt. 0.0)
!      canopy%fess = MIN(canopy%fess, &
!                 & MAX(0.0,REAL(ssoil%wb(:,1),r_1)-soil%swilt/3.0)* soil%zse(1) & ! allow the top soil
!                 & * 1000.0 * hl / dels)                                        ! moisture below wilting
!      canopy%fess = MIN(canopy%fess, &
!                 & REAL((ssoil%wb(:,1)-ssoil%wbice(:,1)),r_1) * soil%zse(1) &
!                 & * 1000.0 * hl / dels)
!    END WHERE
    ! Calculate snow evaporation total in mm (from W/m2):
    ! EAK 'fess' is for soil evap and 'fes' is for soil evap plus soil puddle evap
    canopy%segg = canopy%fess / hl
    WHERE (ssoil%snowd > 0.1) canopy%segg = canopy%fess / (hl + hlf) ! EAK aug08
    ! Initialise snow evaporation:
    ssoil%evapsn = 0
!    ssoil%qssrf = 0.0
    WHERE (ssoil%snowd > 0.1 ) ssoil%evapsn = dels * canopy%fess / ( hl + hlf )
    ! Snow evaporation from 1 snow layer
    WHERE (ssoil%snowd > 0.1 .and. ssoil%isflag == 0 .and. canopy%fess > 0.0)
      ssoil%evapsn = MIN(ssoil%snowd, ssoil%evapsn ) 
      ssoil%snowd = ssoil%snowd - ssoil%evapsn
      canopy%segg = 0.0
    END WHERE
    ! dew on snow
    WHERE (ssoil%snowd > 0.1 .and. ssoil%isflag == 0 .and. canopy%fess < 0.0)
      ssoil%snowd = ssoil%snowd - ssoil%evapsn
      canopy%segg = 0.0
      !canopy%fess = ssoil%evapsn * (hl + hlf) / dels ! return for hyd. balance
    END WHERE
    ! Snow evaporation from 3 snow layers
    WHERE (ssoil%snowd > 0.1 .and. ssoil%isflag > 0 .and. canopy%fess > 0.0)
      ssoil%evapsn = MIN(0.9*ssoil%smass(:,1), ssoil%evapsn )
      ssoil%snowd = ssoil%snowd - ssoil%evapsn
      ssoil%smass(:,1) = ssoil%smass(:,1)  - ssoil%evapsn
      ssoil%sdepth(:,1) = MAX(0.02,ssoil%smass(:,1) / ssoil%ssdn(:,1))
      canopy%segg = 0.0
    END WHERE
    ! dew on snow
    WHERE (ssoil%snowd > 0.1 .and. ssoil%isflag > 0 .and. canopy%fess < 0.0)
      ssoil%snowd = ssoil%snowd - ssoil%evapsn
      ssoil%smass(:,1) = ssoil%smass(:,1)  - ssoil%evapsn
      ssoil%sdepth(:,1) = MAX(0.02,ssoil%smass(:,1) / ssoil%ssdn(:,1))
      canopy%segg = 0.0
    END WHERE


  END SUBROUTINE snow_accum 


  !-------------------------------------------------------------------------
  SUBROUTINE surfbv (dels, ktau, met, ssoil, soil, veg, canopy )
    REAL(r_1), INTENT(IN)                    :: dels   ! integration time step (s)
    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
    TYPE(met_type), INTENT(INOUT)            :: met    ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)    :: veg
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
!    INTEGER(i_d), PARAMETER      :: nglacier = 0 ! 0 original, 1 off, 2 new Eva
    INTEGER(i_d), PARAMETER      :: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    INTEGER(i_d)                 :: k,i
    REAL(r_1), DIMENSION(mp)     :: rnof5
    REAL(r_1), DIMENSION(mp)     :: sfact
    REAL(r_1), DIMENSION(mp)     :: sgamm
    REAL(r_1), DIMENSION(mp)     :: smasstot
!    REAL(r_1), DIMENSION(mp)     :: snr
!    REAL(r_1), DIMENSION(mp)     :: snrat
    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
    REAL(r_1), DIMENSION(mp)     :: talb ! snow albedo
    REAL(r_1), DIMENSION(mp)     :: tmp ! temporary value
    REAL(r_1), DIMENSION(mp)     :: xxx



    CALL smoisturev ( dels, ktau, ssoil, soil)

    ! Diagnostic block below:
    IF (ntest > 0) THEN
      PRINT * , 'in surfbv after smoisturev '
      PRINT * , 'osnowd,snowd,isflag,ssat,runoff', ssoil%osnowd(idjd), &
       & ssoil%snowd(idjd),ssoil%isflag(idjd),soil%ssat(idjd),ssoil%runoff(idjd)
      PRINT * , 'tggsn_d ', (ssoil%tggsn(idjd,k), k = 1, 3)
    END IF
    DO k = 1, ms
      xxx = REAL(soil%ssat,r_2)
      ssoil%rnof1 = ssoil%rnof1 + REAL((MAX(ssoil%wb(:,k) - xxx, 0.0_r_2) &
                                & * 1000.0),r_1) * soil%zse(k)
      ssoil%wb(:,k) = MIN( ssoil%wb(:,k), xxx )
!      ssoil%wb(:,k) = MIN( ssoil%wb(:,k), soil%ssat )
    END DO
    !  for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
    !  account for soil/ice cracking
!    fracm = MIN(0.2, 1. - MIN(1., ssoil%wb(:,ms) / soil%sfc ) )
!    ssoil%wb(:,ms) = ssoil%wb(:,ms) + fracm*ssoil%rnof1/(1000.0*soil%zse(ms))
!    ssoil%rnof1 = (1. - fracm) * ssoil%rnof1 
!    the code below is replaced, see subroutine smoistv 
!    tmp = MAX(MIN(ssoil%wb(:,ms) - soil%sfc, .99 * ssoil%wb(:,ms) &
!         - ssoil%wbice(:,ms) ) * soil%c3 / 86400., 0.)
!    ssoil%rnof2 = soil%zse(ms) * 1000. * tmp * dels
!    ssoil%wb(:,ms) = ssoil%wb(:,ms) - tmp * dels
!   Scaling  runoff to kg/m^2/s to match rest of the model   

! correct for water imbalance in the lakes
    ssoil%sinfil = 0.0
    !where  ( veg%iveg == 16 )
    where  ( veg%iveg == 17 ) ! MJT fix
      ssoil%sinfil = min( ssoil%rnof1, ssoil%wb_lake + max(0.,canopy%segg)) 
      ssoil%rnof1 = max( 0.0, ssoil%rnof1 - ssoil%sinfil ) 
      ssoil%wb_lake = ssoil%wb_lake - ssoil%sinfil
      ssoil%rnof2 = max( 0.0, ssoil%rnof2 - ssoil%wb_lake )
    endwhere        
!
    ssoil%wbtot = 0.0
    DO k = 1, ms
      ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k),r_1) * 1000.0 * soil%zse(k)
    END DO
    !
    !---  glacier formation
    rnof5 = 0.0
    IF (nglacier == 2) THEN
      smelt1 = 0.0
      WHERE (ssoil%snowd > 1100.0)
      rnof5 = min(0.1,ssoil%snowd - 1100.0)
        !---- change local tg to account for energy - clearly not best method
        WHERE (ssoil%isflag == 0)
          smasstot = 0.0
          ssoil%tgg(:,1) = ssoil%tgg(:,1) - rnof5 * hlf &
                         & / REAL(ssoil%gammzz(:,1),r_1)
          ssoil%snowd = ssoil%snowd - rnof5
        ELSEWHERE
          smasstot = ssoil%smass(:,1) + ssoil%smass(:,2) + ssoil%smass(:,3)
        END WHERE
      END WHERE
!      xxx = ssoil%snowd
!      do i=1,mp
!      if( ssoil%snowd(i) > 1100.0) print 12,i,ssoil%snowd(i), &
!      rnof5(i),ssoil%isflag(i),(ssoil%smass(i,k),k=1,3),smasstot(i), &
!      (ssoil%sdepth(i,k),k=1,3),(ssoil%ssdn(i,k),k=1,3),(ssoil%tggsn(i,k),k=1,3), &
!      ssoil%rnof1(i),ssoil%rnof2(i),canopy%precis(i),met%precip_sn(i), &
!      ssoil%evapsn(i),canopy%fess(i),canopy%fes(i),ssoil%qstss(i),met%qv(i), &
!      ssoil%qstss(i)-met%qv(i)
!      enddo
!12    format(1x,'soilsnowpr_1',i6,f9.3,f8.5,i3,3f7.2,f8.2,3f5.2,3f5.0,1x,3f6.1, &
!      f7.4,f7.4,1x,2f7.4,f7.4,1x,2f7.2,3e12.3)
      DO k = 1, 3
        WHERE (ssoil%snowd > 1100.0  .and.  ssoil%isflag > 0)
          sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
          smelt1(:,k) = MIN(rnof5 * ssoil%smass(:,k) / smasstot, &
                          & 0.2 * ssoil%smass(:,k) )
          ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
          ssoil%snowd = ssoil%snowd - smelt1(:,k)
          !ssoil%tggsn(:,k) = ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm
          !ssoil%dtmlt(:,k) = ssoil%dtmlt(:,k) + smelt1(:,k)*hlf/sgamm
        END WHERE
      END DO
      WHERE (ssoil%isflag > 0 ) rnof5 = smelt1(:,1)+smelt1(:,2)+smelt1(:,3)
    END IF
    !38 format(1x,'Les ssoil',i6,8(f9.3,2x))

    ssoil%rnof1 = ssoil%rnof1 / dels + rnof5/dels
    ssoil%rnof2 = ssoil%rnof2 / dels
    ssoil%runoff = ssoil%rnof1 + ssoil%rnof2 

!    do i=1,mp
!    if( xxx(i) > 1100.0) print 14,i,ssoil%snowd(i), &
!    rnof5(i),ssoil%isflag(i),(ssoil%smass(i,k),k=1,3),smasstot(i), &
!        (ssoil%sdepth(i,k),k=1,3),(ssoil%ssdn(i,k),k=1,3),(ssoil%tggsn(i,k),k=1,3), &
!        ssoil%rnof1(i)*dels,ssoil%rnof2(i)*dels,smelt1(i,1),smelt1(i,2),smelt1(i,3), &
!        smelt1(i,1)+smelt1(i,2)+smelt1(i,3)
!    enddo
!14    format(1x,'soilsnowpr_2',i6,f9.3,f8.5,i3,3f7.2,f8.2,3f5.2,3f5.0,1x,3f6.1, &
!      f8.5,f8.5,1x,4f8.5)

    ! Diagnostic block below:
    IF (ntest > 0) THEN
      PRINT * , 'end surfbv  rnof1,runoff ',ssoil%rnof1(idjd),ssoil%runoff(idjd)
      sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
      PRINT * , 'snowd,isflag,sgamm ', ssoil%snowd(idjd), ssoil%isflag(idjd), &
                                     & sgamm(idjd)
      PRINT * , 'tggsn_d ', (ssoil%tggsn(idjd,k), k = 1, 3)
      PRINT * , 'tgg ', (ssoil%tgg(idjd,k), k = 1, ms)
      PRINT * , 'wb ', (ssoil%wb(idjd,k), k = 1, ms)
    END IF
  END SUBROUTINE surfbv


  !-------------------------------------------------------------------------
!  SUBROUTINE snow_albedo (dels, ktau, met, ssoil, soil )
!    REAL(r_1), INTENT(IN)                    :: dels   ! integration time step (s)
!    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
!    TYPE(met_type), INTENT(INOUT)            :: met    ! all met forcing
!    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil  ! soil+snow variables
!    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
!    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
!    REAL(r_1), DIMENSION(mp)     :: alv ! Snow albedo for visible
!    REAL(r_1), DIMENSION(mp)     :: alir ! Snow albedo for near infra-red (NIR)
!    REAL(r_1), PARAMETER         :: alvo  = 0.95 ! albedo for vis. on new snow
!    REAL(r_1), PARAMETER         :: aliro = 0.65 ! albedo for NIR on new snow
!    REAL(r_1), DIMENSION(mp)     :: ar1 ! crystal growth  (-ve)
!    REAL(r_1), DIMENSION(mp)     :: ar2 ! freezing of melt water
!    REAL(r_1), DIMENSION(mp)     :: ar3
!    REAL(r_1), DIMENSION(mp)     :: dnsnow ! new snow albedo
!    REAL(r_1), DIMENSION(mp)     :: dtau
!    REAL(r_1), DIMENSION(mp)     :: fage !age factor
!    REAL(r_1), DIMENSION(mp)     :: fzenm
!    INTEGER(i_d)                 :: k
!    REAL(r_1), DIMENSION(mp)     :: sfact
!    REAL(r_1), DIMENSION(mp)     :: sgamm
!    REAL(r_1), DIMENSION(mp)     :: smasstot
!    REAL(r_1), DIMENSION(mp)     :: snr
!    REAL(r_1), DIMENSION(mp)     :: snrat
!    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
!    REAL(r_1), DIMENSION(mp)     :: talb ! snow albedo
!    REAL(r_1), DIMENSION(mp)     :: tmp ! temporary value
!    REAL(r_1), DIMENSION(mp)     :: xxx

!    !	 calculate soil/snow albedo
!    !	  cuvrf(i,1) = albsav(iq) ! use surface albedo from indata
!    sfact = 0.68
!    WHERE (soil%albsoil <= 0.14)
!      sfact = 0.5
!    ELSEWHERE (soil%albsoil > 0.14 .and. soil%albsoil <= 0.20)
!      sfact = 0.62
!    END WHERE
!    ssoil%albsoilsn(:,2) = 2. * soil%albsoil / (1.0 + sfact)
!    ssoil%albsoilsn(:,1) = sfact * ssoil%albsoilsn(:,2)
!    ! new snow albedo (needs osnowd from the previous dels)
!    dnsnow = MIN( 1.0, 0.1 * MAX(0.0,ssoil%snowd-ssoil%osnowd) ) ! new snow (cm)
!    ! Snow age depends on snow crystal growth, freezing of melt water,
!    ! accumulation of dirt and amount of new snow.
!    tmp = ssoil%isflag * ssoil%tggsn(:,1) + (1 - ssoil%isflag ) * ssoil%tgg(:,1)
!    tmp = MIN(tmp, 273.15)
!    ar1 = 5000.0 * (1.0 / 273.15 - 1.0 / tmp) ! crystal growth  (-ve)
!    ar2 = 10.0 * ar1 ! freezing of melt water
!    snr = ssoil%snowd / MAX(ssoil%ssdnn, 100.0)
!    ! fixes for Arctic & Antarctic
!    WHERE (soil%isoilm == 9)
!      ar3 = 0.0005
!      dnsnow = MAX(dnsnow, 0.002) !increase refreshing of snow in Antarctic
!      snrat = MIN(1.0, snr / (snr + 0.001) )
!    ELSEWHERE
!      ! accumulation of dirt
!      ar3 = 0.1
!      snrat = MIN(1.0, snr / (snr + 0.01) )
!    END WHERE
!    dtau = 1.0e-6 * (EXP(ar1) + EXP(ar2) + ar3) * dels
!    WHERE (ssoil%snowd <= 1.0)
!      ssoil%snage = 0.0
!    ELSEWHERE
!      ssoil%snage = MAX(0.0, (ssoil%snage + dtau) * (1.0 - dnsnow) )
!    END WHERE
!    fage = 1.0 - 1.0 / (1.0 + ssoil%snage ) !age factor
!    !
!    ! Snow albedo is dependent on zenith angle and  snow age.
!    ! albedo zenith dependence
!    ! alvd = alvo * (1.0-cs*fage); alird = aliro * (1.-cn*fage)
!    ! where cs = 0.2, cn = 0.5, b = 2.0
!    tmp = MAX(0.17365, met%coszen )
!    tmp = MAX(0.01, met%coszen )
!    fzenm = MAX(merge(0.0, (1.0+0.5)/(1.0+4.0*tmp)-0.5, tmp > 0.5), 0.0)
!    tmp = alvo * (1.0 - 0.2 * fage)
!    alv = 0.4 * fzenm * (1.0 - tmp) + tmp
!    tmp = aliro * (1.0 - 0.5 * fage)
!   alir = 0.4 * fzenm * (1.0 - tmp) + tmp
!    talb = 0.5 * (alv + alir) ! snow albedo
!    alss = (1. - snrat) * soil%albsoil + snrat * talb ! canopy free surf albedo
!    ssoil%albsoilsn(:,2) = (1.0 - snrat) * ssoil%albsoilsn(:,2) + snrat * alir
!    ssoil%albsoilsn(:,1) = (1.0 - snrat) * ssoil%albsoilsn(:,1) + snrat * alv

!  END SUBROUTINE snow_albedo 


  !-------------------------------------------------------------------------
  ! SUBROUTINE stempv
  !	 calculates temperatures of the soil
  !	 tgg - new soil/snow temperature
  !	 ga - heat flux from the atmosphere (ground heat flux)
  !	 ccnsw - soil thermal conductivity, including water/ice
  !
  SUBROUTINE stempv(dels, canopy, ssoil, soil)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    INTEGER(i_d), PARAMETER          :: ntest = 0
    REAL(r_2), DIMENSION(mp, -2:ms)  :: at
    REAL(r_2), DIMENSION(mp, -2:ms)  :: bt
    REAL(r_2), DIMENSION(mp, -2:ms)  :: ct
    REAL(r_2), DIMENSION(mp,ms)      :: ccnsw  ! soil thermal conductivity
                                               ! (incl water/ice)
    REAL(r_1), DIMENSION(mp)         :: coefa
    REAL(r_1), DIMENSION(mp)         :: coefb
    REAL(r_2), DIMENSION(mp)         :: dtg
    REAL(r_2), DIMENSION(mp)         :: ew
    REAL(r_2), DIMENSION(mp,-2:ms+1) :: coeff
    INTEGER(i_d)                     :: k
    REAL(r_1), DIMENSION(mp,-2:ms)   :: rhs
!    REAL(r_2), DIMENSION(mp,3)       :: sconds
    REAL(r_1), DIMENSION(mp)         :: sgamm
    REAL(r_2), DIMENSION(mp,ms+3)    :: tmp_mat ! temp. matrix for tggsn & tgg
    REAL(r_2), DIMENSION(mp)         :: xx
    REAL(r_2), DIMENSION(mp)         :: wblfsp 
    !
    at = 0.0
    bt = 1.0
    ct = 0.0
    coeff = 0.0
    DO k = 1, ms
      WHERE (soil%isoilm == 9)
!        ccnsw(:,k) = 1.5
        ccnsw(:,k) = 2.0
      ELSEWHERE
        ew = ssoil%wblf(:,k) * soil%ssat
        ccnsw(:,k) = MIN(soil%cnsd * EXP( ew * LOG(60.0) + ssoil%wbfice(:,k) &
                   & * soil%ssat * LOG(250.0) ), 1.5_r_2) &
                   & * MAX(1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 * soil%ssat &
                   & / MIN(ew, 0.5_r_2 * soil%ssat )) ) )
      END WHERE
    END DO
    WHERE (ssoil%isflag == 0)
      xx = MAX(0., ssoil%snowd / ssoil%ssdnn )
      ccnsw(:,1) = (ccnsw(:,1) - 0.2) * (soil%zse(1) / (soil%zse(1) + xx)) + 0.2
    END WHERE
    DO k = 3, ms
      WHERE (ssoil%isflag == 0)
        coeff(:,k) = 2.0 / ( soil%zse(k-1)/ccnsw(:,k-1)+soil%zse(k)/ccnsw(:,k) )
      END WHERE
    END DO
    k = 1
    WHERE (ssoil%isflag == 0)
      coeff(:,2) = 2.0 / ( (soil%zse(1)+xx)/ccnsw(:,1)+soil%zse(2)/ccnsw(:,2) )
      coefa = 0.0
      coefb = REAL(coeff(:,2),r_1)
      wblfsp = ssoil%wblf(:,k)
      xx=soil%css * soil%rhosoil
      ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
           & + soil%ssat * (wblfsp * cswat * rhowat + ssoil%wbfice(:,k) &
           & * csice * rhowat * 0.9), xx ) * soil%zse(k)
      ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
      dtg = dels / ssoil%gammzz(:,k)
      at(:,k) = - dtg * coeff(:,k)
      ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
      bt(:,k) = 1.0 - at(:,k) - ct(:,k)
    END WHERE
    DO k = 2, ms
      WHERE (ssoil%isflag == 0)
        wblfsp = ssoil%wblf(:,k)
        xx=soil%css * soil%rhosoil
        ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat ) * soil%css * soil%rhosoil &
             & + soil%ssat * (wblfsp * cswat * rhowat + ssoil%wbfice(:,k) &
             & * csice * rhowat * 0.9), xx ) * soil%zse(k)
        dtg = dels / ssoil%gammzz(:,k)
        at(:,k) = - dtg * coeff(:,k)
        ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
        bt(:,k) = 1.0 - at(:,k) - ct(:,k)
      END WHERE
    END DO
    WHERE (ssoil%isflag == 0)
      bt(:,1) = bt(:,1) - canopy%dgdtg * dels / ssoil%gammzz(:,1)
      ssoil%tgg(:,1) = ssoil%tgg(:,1) + (canopy%ga - ssoil%tgg(:,1) &
          & * REAL(canopy%dgdtg,r_1)) * dels / REAL(ssoil%gammzz(:,1),r_1)
    END WHERE
    IF (ntest > 0) THEN
      PRINT * ,'stempvtgg1,ga,gammzz',ssoil%isflag(idjd),ssoil%tgg(idjd,:),canopy%ga(idjd), &
              & ssoil%gammzz(idjd,:)
      PRINT * , 'dgdtg', canopy%dgdtg(idjd),ccnsw(idjd,:)
      PRINT * , 'ssat,css,rhos,cswat,rhowat,csice ', soil%ssat(idjd), &
              & soil%css(idjd), soil%rhosoil(idjd), cswat, rhowat, csice
      PRINT * , 'wblf1,wbfice1,zse1,cgsnow ', ssoil%wblf(idjd,1), &
              & ssoil%wbfice(idjd,1), soil%zse(:), cgsnow
!      PRINT * , 'at ', (at(idjd,k) , k = 1, ms)
!      PRINT * , 'bt ', (bt(idjd,k) , k = 1, ms)
!      PRINT * , 'ct ', (ct(idjd,k) , k = 1, ms)
      PRINT * , 'rhs ', (ssoil%tgg(idjd,k) , k = 1, ms)
    END IF
    coeff(:,1-3) = 0.0
    ! 3-layer snow points done here
    WHERE (ssoil%isflag /= 0)
      ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1)**2 &
                        & + 0.074, 2.51) )
      ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2)**2 &
                        & + 0.074, 2.51) )
      ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3)**2 &
                        & + 0.074, 2.51) )
      coeff(:,-1) = 2.0 / (ssoil%sdepth(:,1) / ssoil%sconds(:,1) &
                       & + ssoil%sdepth(:,2) / ssoil%sconds(:,2) )
      coeff(:,0) = 2.0 / (ssoil%sdepth(:,2) / ssoil%sconds(:,2) &
                      & + ssoil%sdepth(:,3) / ssoil%sconds(:,3) )
      coeff(:,1) = 2.0 / (ssoil%sdepth(:,3) / ssoil%sconds(:,3) &
                      & + soil%zse(1) / ccnsw (:,1) )
    END WHERE
    DO k = 2, ms
      WHERE (ssoil%isflag /= 0)
        coeff(:,k) = 2.0 / (soil%zse(k-1)/ccnsw(:,k-1) + soil%zse(k)/ccnsw(:,k))
      END WHERE
    END DO
    WHERE (ssoil%isflag /= 0)
      coefa = REAL(coeff (:,-1),r_1)
      coefb = REAL(coeff (:,1),r_1)
    END WHERE
!    print *,'stempv coeff isfl=1',coeff(idjd,:),ssoil%sconds(idjd,:)
!    print *,'stempv coeff 2 ', dels,ssoil%sdepth(idjd,:),ssoil%ssdn(idjd,:)
    DO k = 1, 3
      WHERE (ssoil%isflag /= 0)
        sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
        dtg = dels / sgamm
        ! rhs(k-3) = ssoil%tggsn(:,k)  ! A
        at(:,k-3) = - dtg * coeff(:,k-3)
        ct(:,k-3) = - dtg * coeff(:,k-2)
        bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
      END WHERE
    END DO
    DO k = 1, ms
      WHERE (ssoil%isflag /= 0)
        wblfsp = ssoil%wblf(:,k)
        xx=soil%css * soil%rhosoil
        ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat ) * soil%css * soil%rhosoil &
             & + soil%ssat * (wblfsp * cswat * rhowat + ssoil%wbfice(:,k) &
             & * csice * rhowat * 0.9), xx ) * soil%zse(k)
        dtg = dels / ssoil%gammzz(:,k)
        at(:,k) = - dtg * coeff(:,k)
        ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
        bt(:,k) = 1.0 - at(:,k) - ct(:,k)
      END WHERE
    END DO
!    print *,'stempv  at,bt,ct',ssoil%gammzz(idjd,:),at(idjd,:),canopy%dgdtg(idjd)
    WHERE (ssoil%isflag /= 0)
      sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
      ! rhs(1-3) = rhs(1-3)+canopy%ga*dels/sgamm
      ! new code
      bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm
      ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + (canopy%ga - ssoil%tggsn(:,1) &
           & * REAL(canopy%dgdtg,r_1) ) * dels / sgamm
      rhs(:,1-3) = ssoil%tggsn(:,1)
    END WHERE
    IF (ntest > 0) THEN
       IF (ssoil%isflag(idjd) /= 0) THEN
          PRINT * , 'in stempv 3-layer snow code '
          PRINT * , 'ccnsw ', (ccnsw(idjd,k), k = 1, ms)
          PRINT * , 'sdepth d ', (ssoil%sdepth(idjd,k), k = 1, 3)
          PRINT * , 'sconds ', ssoil%sconds(idjd,:)
          !      PRINT * , 'coeff ', coeff
          PRINT * , 'at ', (at(idjd,k), k = - 2, ms)
          PRINT * , 'bt ', (bt(idjd,k), k = - 2, ms)
          PRINT * , 'ct ', (ct(idjd,k), k = - 2, ms)
          PRINT * , 'rhs(tggsn,tgg) ', (ssoil%tggsn(idjd,k), k = 1,3), &
               & (ssoil%tgg(idjd,k), k = 1,ms)
          PRINT * , 'tggsn,ga,sgamm ', ssoil%tggsn(idjd,1), canopy%ga(idjd), &
               & sgamm(idjd)
          PRINT * , 'dgdtg ', canopy%dgdtg(idjd)
       END IF
    END IF
    !
    !     note in the following that tgg and tggsn are processed together
    tmp_mat(:,:3) = REAL(ssoil%tggsn,r_2)
    tmp_mat(:,4:) = REAL(ssoil%tgg,r_2)
    CALL trimb (at, bt, ct, tmp_mat, ms + 3)
    ssoil%tggsn = REAL(tmp_mat(:,:3),r_1)
    ssoil%tgg   = REAL(tmp_mat(:,4:),r_1)

    IF (ntest > 0) THEN
       IF (ssoil%isflag(idjd) /= 0) THEN
          PRINT * , 'afrhs(tggsn,tgg) ', (ssoil%tggsn(idjd,k), k = 1,3), &
              & (ssoil%tgg(idjd,k), k = 1,ms)
          PRINT * , 'aftggsn,ga,sgamm ', ssoil%tggsn(idjd,1), canopy%ga(idjd), &
              & sgamm(idjd)
       END IF
    END IF
 
    canopy%sghflux = coefa * (ssoil%tggsn(:,1) - ssoil%tggsn(:,2) )
!    canopy%sghflux = coefb * (ssoil%tggsn(:,3) - ssoil%tgg(:,1) )
    canopy%ghflux = coefb * (ssoil%tgg(:,1) - ssoil%tgg(:,2) ) ! +ve downwards
  END SUBROUTINE stempv


  !-------------------------------------------------------------------------
  ! SUBROUTINE stempvsn
  !	 calculates temperatures of the three layer snowpack
  !	 tgg - new soil/snow temperature
  !	 ga - heat flux from the atmosphere 
  !	 ccnsw - soil thermal conductivity, including water/ice
  !
  SUBROUTINE stempvsn(dels,soilcond, soilhcp, canopy, ssoil, soil)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    REAL(r_1), DIMENSION(mp), INTENT(IN)     :: soilcond   
    REAL(r_1), DIMENSION(mp), INTENT(IN)     :: soilhcp   
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    INTEGER(i_d), PARAMETER            :: mssn = 1
    INTEGER(i_d), PARAMETER            :: ntest = 0
    REAL(r_2), DIMENSION(mp, -2:mssn)  :: at
    REAL(r_2), DIMENSION(mp, -2:mssn)  :: bt
    REAL(r_2), DIMENSION(mp, -2:mssn)  :: ct
    REAL(r_2), DIMENSION(mp,mssn)      :: ccnsw  ! soil thermal conductivity
                                                 ! (incl water/ice)
    REAL(r_1), DIMENSION(mp)           :: coefa
    REAL(r_1), DIMENSION(mp)           :: coefb
    REAL(r_2), DIMENSION(mp)           :: dtg
    REAL(r_2), DIMENSION(mp)           :: ew
    REAL(r_2), DIMENSION(mp,-2:mssn+1) :: coeff
    INTEGER(i_d)                       :: k
    REAL(r_1), DIMENSION(mp,-2:mssn)   :: rhs
!    REAL(r_2), DIMENSION(mp,3)         :: sconds
    REAL(r_1), DIMENSION(mp)           :: sgamm
    REAL(r_2), DIMENSION(mp,mssn+3)    :: tmp_mat ! temp. matrix for tggsn & tgg
    REAL(r_2), DIMENSION(mp)           :: xx
    REAL(r_2), DIMENSION(mp)           :: wblfsp 
    !
    at = 0.0
    bt = 1.0
    ct = 0.0
    coeff = 0.0
    canopy%dgdtg =0.0
    ccnsw(:,1) = soilcond  ! from Jules
    coeff (:,1-3) = 0.0
    ! 3-layer snow points done here
    WHERE (ssoil%isflag /= 0)
      coeff(:,-1) = 2.0 / (ssoil%sdepth(:,1) / ssoil%sconds(:,1) &
                       & + ssoil%sdepth(:,2) / ssoil%sconds(:,2) )
      coeff(:,0) = 2.0 / (ssoil%sdepth(:,2) / ssoil%sconds(:,2) &
                      & + ssoil%sdepth(:,3) / ssoil%sconds(:,3) )
      coeff(:,1) = 2.0 / (ssoil%sdepth(:,3) / ssoil%sconds(:,3) &
                      & + 0.1 / ccnsw(:,1) )
    END WHERE
    WHERE (ssoil%isflag /= 0)
      coefa = REAL(coeff(:,-1),r_1)
      coefb = REAL(coeff(:,1),r_1)
    END WHERE
    DO k = 1, 3
      WHERE (ssoil%isflag /= 0)
        sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
        dtg = dels / sgamm
        ! rhs(k-3) = ssoil%tggsn(:,k)	  ! A
        at(:,k-3) = - dtg * coeff(:,k-3)
        ct(:,k-3) = - dtg * coeff(:,k-2)
        bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
      END WHERE
    END DO
    DO k = 1, mssn
      WHERE (ssoil%isflag /= 0)
        dtg = dels / soilhcp
        at(:,k) = - dtg * coeff(:,k)
        ct(:,k) = - dtg * coeff(:,k + 1) ! c3(mssn)=0 & not really used
        bt(:,k) = 1.0 - at(:,k) - ct(:,k)
      END WHERE
    END DO
    WHERE (ssoil%isflag /= 0)
      sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
      ! rhs(1-3) = rhs(1-3)+canopy%ga*dels/sgamm
      ! new code
      bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm
      ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + (canopy%ga - ssoil%tggsn(:,1) &
           & * REAL(canopy%dgdtg,r_1) ) * dels / sgamm
      rhs(:,1-3) = ssoil%tggsn(:,1)
    END WHERE
    IF (ntest > 0) THEN
       IF (ssoil%isflag(idjd) /= 0) THEN
          PRINT * , 'in stempv 3-layer snow code '
          PRINT * , 'ccnsw ', (ccnsw(idjd,k), k = 1, mssn)
          PRINT * , 'sdepth d ', (ssoil%sdepth(idjd,k), k = 1, 3)
          PRINT * , 'sconds ', ssoil%sconds
          PRINT * , 'coeff ', coeff
          PRINT * , 'at ', (at(idjd,k) , k = - 2, mssn)
          PRINT * , 'bt ', (bt(idjd,k) , k = - 2, mssn)
          PRINT * , 'ct ', (ct(idjd,k) , k = - 2, mssn)
          PRINT * , 'rhs(tggsn,tgg)',(ssoil%tggsn(idjd,k),k = 1,4)
          PRINT * , 'tggsn,ga,sgamm ', ssoil%tggsn (idjd,1) , canopy%ga(idjd) , &
               & sgamm(idjd)
       END IF
    END IF
    !
    !     note in the following that tgg and tggsn are processed together
    tmp_mat(:,:3) = REAL(ssoil%tggsn,r_2)
    tmp_mat(:,4:) = REAL(ssoil%tgg,r_2)
    CALL trimb (at, bt, ct, tmp_mat, mssn + 3)
    ssoil%tggsn = REAL(tmp_mat(:,:3),r_1)
!    canopy%sghflux = coefb * (ssoil%tggsn(:,3) - ssoil%tgg(:,1) )
!    WHERE ( ssoil%isflag == 1 ) canopy%ga = canopy%sghflux
!    print *,'stempvsn ga',canopy%ga
  END SUBROUTINE stempvsn


  !-------------------------------------------------------------------------
  SUBROUTINE snowcheck(dels, ktau, ssoil, soil, met, ktau_gl )
!  SUBROUTINE snowcheck(dels, ktau_gl, ssoil, soil, met )
    REAL(r_1), INTENT(IN)               :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    INTEGER(i_d), INTENT(IN)            :: ktau_gl 
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    INTEGER(i_d), PARAMETER :: ntest = 0 !  for prints
    INTEGER(i_d)            :: k

    ! Diagnostic block:
    IF (ntest > 0) THEN
      PRINT *, 'in snowcheck ktau= ', ktau, ssoil%isflag(idjd), idjd
      PRINT *, 'dels,ssdn ', dels, (ssoil%ssdn(idjd,k),k=1,3)
      PRINT *, 'osnowd,snowd', ssoil%osnowd(idjd), ssoil%snowd(idjd)
      PRINT *, 'tgg ', (ssoil%tgg(idjd,k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(idjd,k),k=1,3)
    END IF

!!$    IF (ktau <= 1) THEN
!!$      ssoil%ssdn = 120.0
!!$      ssoil%ssdnn = 120.0
!!$      ssoil%tggsn = tfrz
!!$      ssoil%isflag = 0
!!$      ssoil%sdepth(:,1) = ssoil%snowd / ssoil%ssdn(:,1)
!!$      ssoil%sdepth(:,2) = 0.
!!$      ssoil%sdepth(:,3) = 0.
!!$      ssoil%smass(:,1) = ssoil%snowd
!!$      ssoil%smass(:,2) = 0.0     
!!$      ssoil%smass(:,3) = 0.0    
!!$      WHERE( soil%isoilm == 9 )
!!$        ssoil%ssdn(:,1)  = 200.0
!!$        ssoil%ssdn(:,2)  = 300.0
!!$        ssoil%ssdn(:,3)  = 380.0
!!$        ssoil%ssdnn      = 350.0
!!$      ENDWHERE
!!$
!!$    END IF

    WHERE (ssoil%snowd <= 0.0)
      ssoil%isflag = 0
      ssoil%ssdn(:,1) = 120.0
      ssoil%ssdn(:,2) = 120.0
      ssoil%ssdn(:,3) = 120.0
      ssoil%ssdnn = 120.0
      ssoil%tggsn(:,1) = tfrz
      ssoil%tggsn(:,2) = tfrz
      ssoil%tggsn(:,3) = tfrz
      ssoil%sdepth(:,1) = ssoil%snowd / ssoil%ssdn(:,1)
      ssoil%sdepth(:,2) = 0.
      ssoil%sdepth(:,3) = 0.
      ssoil%smass(:,1) = ssoil%snowd
      ssoil%smass(:,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
      ssoil%smass(:,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
    ELSEWHERE (ssoil%snowd < snmin * ssoil%ssdnn) ! original version for CABLE

!    assume land ice points in UM have been initialised with 50,000 cm of snow
!    and 3 layer snow pack can be open only outside land ice points
!    ELSEWHERE ((ssoil%snowd < snmin*ssoil%ssdnn).or.(ssoil%snowd > 40000.0))

      WHERE (ssoil%isflag == 1)
        ssoil%ssdn(:,1) = ssoil%ssdnn
        ssoil%tgg(:,1) = ssoil%tggsn(:,1)
      END WHERE
!      ssoil%ssdnn = MIN( 400.0, MAX(120.0, ssoil%ssdn(:,1)) )
      ssoil%ssdnn = MAX(120.0, ssoil%ssdn(:,1)) 
      ssoil%ssdn(:,1) = ssoil%ssdnn
      ssoil%ssdn(:,2) = ssoil%ssdnn
      ssoil%ssdn(:,3) = ssoil%ssdnn
      ssoil%isflag = 0
      ssoil%tggsn(:,1) = min(tfrz,ssoil%tgg(:,1))
      ssoil%tggsn(:,2) = min(tfrz,ssoil%tgg(:,1))
      ssoil%tggsn(:,3) = min(tfrz,ssoil%tgg(:,1))
      ssoil%sdepth(:,1) = ssoil%snowd / ssoil%ssdn(:,1)
      ssoil%sdepth(:,2) = 0.0     
      ssoil%sdepth(:,3) = 0.0     
      ssoil%smass(:,1) = ssoil%snowd     
      ssoil%smass(:,2) = 0.0     
      ssoil%smass(:,3) = 0.0     
      where( soil%isoilm == 9 .and. ktau_gl <= 2 )
      ssoil%ssdnn = 700.0
      endwhere

    ELSEWHERE ! sufficient snow now for 3 layer snowpack

      WHERE (ssoil%isflag == 0)
        ssoil%tggsn(:,1) = min(tfrz,ssoil%tgg(:,1))
        ssoil%tggsn(:,2) = min(tfrz,ssoil%tgg(:,1))
        ssoil%tggsn(:,3) = min(tfrz,ssoil%tgg(:,1))
        ssoil%ssdn(:,2) = ssoil%ssdn(:,1)
        ssoil%ssdn(:,3) = ssoil%ssdn(:,1)
        WHERE( soil%isoilm == 9 .and. ktau_gl <= 2 )
          ssoil%ssdn(:,1)  = 450.0
          ssoil%ssdn(:,2)  = 580.0
          ssoil%ssdn(:,3)  = 600.0
!          ssoil%tggsn(:,3) = min(210.0,ssoil%tggav)
!          ssoil%tggsn(:,2) = min(213.0,0.5*(ssoil%tggsn(:,3) + met%tk))
!          ssoil%tggsn(:,1) = min(ssoil%tggsn(:,2),met%tk)
        ENDWHERE
!        ssoil%sdepth(:,1) = 0.05
        ssoil%sdepth(:,1) = ssoil%t_snwlr
        ! next 5 lines replaced to fix -ve sdepth (EK 21Dec2007)
!        ssoil%smass(:,1)  = 0.05 * ssoil%ssdn(:,1)
        ssoil%smass(:,1)  =  ssoil%t_snwlr * ssoil%ssdn(:,1)
        ssoil%smass(:,2)  = (ssoil%snowd - ssoil%smass(:,1)) * 0.4
        ssoil%sdepth(:,2) = ssoil%smass(:,2)/ssoil%ssdn(:,2)
        ssoil%smass(:,3)  = (ssoil%snowd - ssoil%smass(:,1)) * 0.6
        ssoil%sdepth(:,3) = ssoil%smass(:,3)/ssoil%ssdn(:,3)
        ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%smass(:,1) + ssoil%ssdn(:,2) &
             & * ssoil%smass(:,2) + ssoil%ssdn(:,3) * ssoil%smass(:,3) ) &
             & / ssoil%snowd

!        ssoil%sdepth(:,2) =  (ssoil%snowd / ssoil%ssdn(:,1) - 0.05) &
!             & * merge(0.4, 0.4, ssoil%snowd > 20.0)
!        ssoil%sdepth(:,3) = (ssoil%snowd / ssoil%ssdn(:,1) - 0.05) &
!             & * merge(0.6, 0.6, ssoil%snowd > 20.0)
!        ssoil%smass(:,1) = 0.05 * ssoil%ssdn(:,1)
!        ssoil%smass(:,2) = ssoil%sdepth(:,2) * ssoil%ssdn(:,2)
!        ssoil%smass(:,3) = ssoil%sdepth(:,3) * ssoil%ssdn(:,3)

      END WHERE
      ssoil%isflag = 1
    END WHERE

  END SUBROUTINE snowcheck 


  !******************************************************************
  SUBROUTINE snowl_adjust(dels, ktau,  ssoil, canopy )
    REAL(r_1), INTENT(IN)               :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE(canopy_type), INTENT(INOUT)    :: canopy
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints
    INTEGER(i_d)             :: k
    REAL(r_2), DIMENSION(mp) :: excd
    REAL(r_2), DIMENSION(mp) :: excm
    REAL(r_2), DIMENSION(mp) :: frac 
    REAL(r_2), DIMENSION(mp) :: xfrac 
    REAL(r_1), DIMENSION(mp) :: osm
   !						~.11 to turn on 3-layer snow
    ! Diagnostic block:
    IF (ntest > 0) THEN
       PRINT *, 'in soilsnowv before stempv,  ktau= ', ktau,ssoil%isflag(idjd)
       PRINT *, 'soilsnowv ', canopy%dgdtg(idjd), canopy%fevc(idjd), &
              & canopy%fess(idjd), canopy%precis(idjd)
       PRINT *, 'ga,dels,ssdn ', canopy%ga(idjd), dels, (ssoil%ssdn(idjd,k),k=1,3)
       PRINT *, 'osnowd,snowd,isflag', ssoil%osnowd(idjd), ssoil%snowd(idjd), &
              & ssoil%isflag(idjd)
       PRINT *, 'tgg ', (ssoil%tgg(idjd,k),k=1,ms)
       PRINT *, 'tggsn ', (ssoil%tggsn(idjd,k),k=1,3)
       PRINT *, 'wb ', (ssoil%wb(idjd,k),k=1,ms)
       PRINT *, 'wbice ', (ssoil%wbice(idjd,k),k=1,ms)
    END IF

    ! adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    WHERE (ssoil%isflag > 0)
!      WHERE ( ssoil%sdepth(:,1) .gt. 0.05 )
      WHERE ( ssoil%sdepth(:,1) .gt.  ssoil%t_snwlr )
!        excd = ssoil%sdepth(:,1) - 0.05
        excd = ssoil%sdepth(:,1) - ssoil%t_snwlr
        excm = excd * ssoil%ssdn(:,1)
        ssoil%sdepth(:,1) = ssoil%sdepth(:,1) - excd
!        ssoil%sdepth(:,1) = 0.05
        osm = ssoil%smass(:,1)
        ssoil%smass(:,1) = ssoil%smass(:,1) - excm

        osm = ssoil%smass(:,2)
        ssoil%smass(:,2) = MAX(0.01_r_2, ssoil%smass(:,2) + excm)
        ssoil%ssdn(:,2) = MAX(120.0_r_2, MIN(750.0_r_2, ssoil%ssdn(:,2)* &
          & osm/ssoil%smass(:,2) + ssoil%ssdn(:,1) * excm / ssoil%smass(:,2)) )
        ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
        ssoil%tggsn(:,2) = ssoil%tggsn(:,2) * osm / ssoil%smass(:,2) + &
                         & ssoil%tggsn(:,1) * excm/ ssoil%smass(:,2)
        ! following line changed to fix -ve sdepth (EK 21Dec2007)
        ssoil%smass(:,3) = MAX(0.01_r_2, &
                          & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
!        ssoil%smass(:,3) = ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2)
!      ELSEWHERE ! ssoil%sdepth(:,1) < 0.05
      ELSEWHERE ! ssoil%sdepth(:,1) <  ssoil%t_snwlr
        ! 1st layer
!        excd = 0.05 - ssoil%sdepth(:,1)
        excd =  ssoil%t_snwlr - ssoil%sdepth(:,1)
        excm = excd * ssoil%ssdn(:,2)
        osm = ssoil%smass(:,1)
        ssoil%smass(:,1) = ssoil%smass(:,1) + excm
!        ssoil%sdepth(:,1) = 0.05
        ssoil%sdepth(:,1) =  ssoil%t_snwlr
        ssoil%ssdn(:,1) = MAX(120.0_r_2, MIN(750.0_r_2, ssoil%ssdn(:,1)* &
        & osm/ssoil%smass(:,1) + ssoil%ssdn(:,2) * excm / ssoil%smass(:,1)) )
        ssoil%tggsn(:,1) = ssoil%tggsn(:,1) * osm / ssoil%smass(:,1) +  &
                         & ssoil%tggsn(:,2) * excm/ ssoil%smass(:,1)
        ! 2nd layer
        ssoil%smass(:,2) = MAX(0.01_r_2, ssoil%smass(:,2) - excm)
        ssoil%sdepth(:,2) = ssoil%smass(:,2)/ssoil%ssdn(:,2)
        ssoil%smass(:,3) = MAX(0.01_r_2, &
                          & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
!        ssoil%smass(:,3) = ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2)
      END WHERE
    END WHERE
!!$    write(0,*) "MRD SSOIL ISFLAG", ssoil%isflag
!!$    write(0,*) "MRD SMASS 2 ", ssoil%smass(:,2)
!!$    write(0,*) "MRD SMASS 3 ", ssoil%smass(:,3)
    xfrac = 0. ! Work around NaN initialisation problem
    WHERE (ssoil%isflag > 0)
      frac = ssoil%smass(:,2) / MAX(0.02, ssoil%smass(:,3))
      ! if frac > 0.6 or frac < 0.74 do nothing HOW TO translate this to xfrac
      xfrac = 2.0/3.0/ frac
      WHERE ( xfrac > 1.0 )
        excm = (xfrac - 1.0) * ssoil%smass(:,2)
        osm = ssoil%smass(:,2)
        ! changed 0.02 to 0.01 to fix -ve sdepth (EK 21Dec2007)
        ssoil%smass(:,2) = MAX(0.01_r_2, ssoil%smass(:,2) + excm)
        ssoil%tggsn(:,2) = ssoil%tggsn(:,2) * osm / ssoil%smass(:,2) +  &
                         & ssoil%tggsn(:,3) * excm/ ssoil%smass(:,2)
        ssoil%ssdn(:,2) = MAX(120.0_r_2, MIN(750.0_r_2, ssoil%ssdn(:,2)* &
           & osm/ssoil%smass(:,2) + ssoil%ssdn(:,3) * excm / ssoil%smass(:,2)) )
        ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
        ssoil%smass(:,3) = MAX(0.01_r_2, &
                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
        ssoil%sdepth(:,3) = MAX(0.02, ssoil%smass(:,3) / ssoil%ssdn(:,3) )
      ELSEWHERE ! xfrac < 1
        excm = (1 - xfrac) * ssoil%smass(:,2)
        ssoil%smass(:,2) = MAX(0.01_r_2, ssoil%smass(:,2) - excm)
        ssoil%sdepth(:,2) = MAX(0.02, ssoil%smass(:,2) / ssoil%ssdn(:,2) )

        osm = ssoil%smass(:,3)
        ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
        ssoil%smass(:,3) = MAX(0.01_r_2, &
                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
        ssoil%tggsn(:,3) = ssoil%tggsn(:,3) * osm / ssoil%smass(:,3) +  &
                         & ssoil%tggsn(:,2) * excm/ ssoil%smass(:,3)
        ssoil%ssdn(:,3) = MAX(120.0_r_2, MIN(750.0_r_2, ssoil%ssdn(:,3)* &
          & osm/ssoil%smass(:,3) + ssoil%ssdn(:,2) * excm / ssoil%smass(:,3)) )
        ssoil%sdepth(:,3) = ssoil%smass(:,3) /  ssoil%ssdn(:,3)
      END WHERE
      ssoil%isflag = 1
      ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%sdepth(:,1) + ssoil%ssdn(:,2) &
             & * ssoil%sdepth(:,2) + ssoil%ssdn(:,3) * ssoil%sdepth(:,3) ) &
             & / (ssoil%sdepth(:,1) + ssoil%sdepth(:,2) + ssoil%sdepth(:,3))
    END WHERE

    ! Diagnostic block:
    IF (ntest > 0) THEN
      PRINT *, 'in soilsnowv before stempv,  ktau= ',ktau
      PRINT *, 'ga,dels,ssdn ', canopy%ga(idjd), dels, (ssoil%ssdn(idjd,k),k=1,3)
      PRINT *, 'osnowd,snowd,isflag', ssoil%osnowd(idjd), ssoil%snowd(idjd), &
             & ssoil%isflag(idjd)
      PRINT *, 'tgg ', (ssoil%tgg(idjd,k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(idjd,k),k=1,3)
      PRINT *, 'wb ', (ssoil%wb(idjd,k),k=1,ms)
      PRINT *, 'wbice ', (ssoil%wbice(idjd,k),k=1,ms)
      PRINT *, 'wblf ', (ssoil%wblf(idjd,k),k=1,ms)
      PRINT *, 'wbfice ', (ssoil%wbfice(idjd,k),k=1,ms)
    END IF
    IF (ntest == 1) THEN
      PRINT *, 'in soilsnow printing wbfice_max'
      PRINT *, 'sdepth c2 ', (ssoil%sdepth(idjd,k),k=1,3)
    END IF

  END SUBROUTINE snowl_adjust


  !******************************************************************
  SUBROUTINE soilfreeze(dels, soil, ssoil)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    REAL(r_1), DIMENSION(mp)           :: sicefreeze
    REAL(r_1), DIMENSION(mp)           :: sicemelt
    REAL(r_1), DIMENSION(mp)           :: xx
    INTEGER(i_d) k
    ! allow for some water (< 5%) to remain unfrozen
    DO k = 1, ms
      WHERE (ssoil%tgg(:,k) < tfrz &
          & .and. 0.80 * ssoil%wb(:,k) - ssoil%wbice(:,k) > .001)
        sicefreeze = MIN( MAX(0.0_r_2,(0.80*ssoil%wb(:,k)-ssoil%wbice(:,k))) &
                        & * soil%zse(k) * 1000.0, &
                        & (tfrz - ssoil%tgg(:,k) ) * ssoil%gammzz(:,k) / hlf )
        ssoil%wbice(:,k) = MIN( ssoil%wbice(:,k) + sicefreeze / (soil%zse(k) &
                         & * 1000.0), 0.80 * ssoil%wb(:,k) )
        xx=soil%css * soil%rhosoil
        ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
             & + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * cswat * rhowat &
             & + ssoil%wbice(:,k) * csice * rhowat * 0.9, xx) &
             & * soil%zse(k)
        WHERE (k == 1 .and. ssoil%isflag == 0)
          ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
        END WHERE
        ssoil%tgg(:,k) = ssoil%tgg(:,k) + REAL(sicefreeze,r_1) &
                       & * hlf / REAL(ssoil%gammzz(:,k),r_1)
      ELSEWHERE (ssoil%tgg(:,k) > tfrz .and. ssoil%wbice(:,k) > 0.)
        sicemelt = MIN(ssoil%wbice(:,k) * soil%zse(k) * 1000.0, &
               (ssoil%tgg(:,k) - tfrz) * ssoil%gammzz(:,k) / hlf)
        ssoil%wbice(:,k) = MAX(0.0_r_2, ssoil%wbice(:,k) - sicemelt &
                         & / (soil%zse(k) * 1000.0) )
        xx = soil%css * soil%rhosoil
        ssoil%gammzz(:,k) = MAX( (1.0-soil%ssat) * soil%css * soil%rhosoil &
             & + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * cswat * rhowat &
             & + ssoil%wbice(:,k) * csice * rhowat * 0.9, xx ) * soil%zse(k)
        WHERE (k == 1 .and. ssoil%isflag == 0)
          ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
        END WHERE
        ssoil%tgg(:,k) = ssoil%tgg(:,k) - REAL(sicemelt,r_1) &
                       & * hlf / REAL(ssoil%gammzz(:,k),r_1)
      END WHERE
    END DO
  END SUBROUTINE soilfreeze


  !******************************************************************
!  SUBROUTINE remove_trans(dels, soil, ssoil, canopy, veg)
    SUBROUTINE remove_trans(dels, soil, ssoil, canopy)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
!    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    REAL(r_1), DIMENSION(mp,ms)   :: evapfb
    REAL(r_1), DIMENSION(mp,0:ms) :: diff 
    REAL(r_1), DIMENSION(mp)      :: xx
    INTEGER(i_d) k

    diff(:,:) = 0.0  
!    diff(:,0) = 0.0
    evapfb(:,:) = 0.0

!    print *,'remove_trans1',canopy%fevc(352),soil%swilt(352)
!    print *,'remove_trans1a',ssoil%wb(352,:)
    DO k = 1,ms
      ! Removing transpiration from soil:
      WHERE (canopy%fevc > 0.0)     ! convert to mm/dels
        ! Calculate the amount (perhaps moisture/ice limited)
        ! which can be removed:
        !xx = ssoil%evapfbl(:,k) + diff(:,k-1)   ! kg/m2
        xx = canopy%fevc * dels / hl * soil%froot(:,k) + diff(:,k-1)   ! kg/m2
        diff(:,k) = ( MAX( 0.0, ssoil%wb(:,k) - soil%swilt) &      ! m3/m3
                  & * soil%zse(k) * 1000.0 - xx) / (soil%zse(k) * 1000.0)
!        diff(:,k) = (MAX( 0.0, MIN( ssoil%wb(:,k) - soil%swilt, &      ! m3/m3
!                  & ssoil%wb(:,k) - ssoil%wbice(:,k) ) ) &
!                  & * soil%zse(k) * 1000.0 - xx) / (soil%zse(k) * 1000.0)
        WHERE ( diff(:,k) .gt. 0.0 )
          ssoil%wb(:,k) = ssoil%wb(:,k) - xx / (soil%zse(k) * 1000.0)
          diff(:,k) = 0.0
          evapfb(:,k) = xx 
        ELSEWHERE
          diff(:,k) = xx
        ENDWHERE
      END WHERE
    END DO

    ! Adjust fevc
!    WHERE (canopy%fevc > 0.)  
!      canopy%fevc = (evapfbl(:,1)+evapfbl(:,2)+evapfbl(:,3)+evapfbl(:,4))*hl/dels
!!        & + evapfbl(:,5)+evapfbl(:,6))*hl/dels
!    END WHERE

!    print *,'remove_trans2',canopy%fevc(352),soil%swilt(352)
!    print *,'remove_trans2a',ssoil%wb(352,:)
!    print *,'remove_trans3a',ssoil%wb(352,:)

  END SUBROUTINE remove_trans 

  !----------------------------------------------------------------------------
  ! SUBROUTINE soil_snow
  !
  ! Replaces following
  ! SUBROUTINE soilsnow (dt_in, ktau_in, ga, dgdtg, condxpr, scondxpr, fev, fes, t, coszen)
  !	    for snow diag prints set ntest to 1 throughout
  !	    or, usefully can edit 'ntest > 0' to 'ktau > nnn'
  !----------------------------------------------------------------------
  ! Inputs:
  !	 dt_in - time step in sec
  !	 ktau_in - time step no.
  !	 ga	 - ground heat flux W/m^2
  !	 dgdtg	 -
  !	 condxpr - total precip reaching the ground (liquid and solid)
  !	 scondxpr - precip (solid only)
  !	 fev   - transpiration (W/m2)
  !	 fess   - soil evaporation (W/m2)
  !	 isoil - soil type
  !	 ivegt - vegetation type
  ! Output
  !	 ssoil
! SUBROUTINE soil_snow(dels, ktau, soil, ssoil, veg,canopy, met, bal)
  SUBROUTINE soil_snow(dels, ktau, soil, ssoil, canopy, met, bal, veg, ktau_gl)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
    INTEGER(i_d), INTENT(IN)                 :: ktau_gl ! global step number
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
!    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    TYPE (balances_type), INTENT(INOUT)      :: bal
    TYPE (veg_parameter_type), INTENT(IN)    :: veg
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints
    INTEGER(i_d)             :: k
!    REAL(r_2), DIMENSION(mp):: sicefreeze
!    REAL(r_2), DIMENSION(mp):: sicemelt
    REAL(r_1), DIMENSION(mp) :: snowmlt
    REAL(r_1), DIMENSION(mp) :: totwet
    REAL(r_1), DIMENSION(mp) :: weting
    REAL(r_1), DIMENSION(mp) :: xxx
    REAL(r_2), DIMENSION(mp) :: sinfil1,sinfil2,sinfil3
    REAL(r_2), DIMENSION(mp) :: xx,deltat
!    REAL(r_2), DIMENSION(mp) :: oldpudsto,newpudsto
    REAL(r_1)                :: zsetot
    REAL(r_1), DIMENSION(mp) :: tgg_old,tggsn_old 
    INTEGER(i_d) :: idjd1,idjd2,idjd3
!%%

    idjd1 = 351
    idjd2 = 8285
    idjd3 = 9985


    ! Diagnostic block:
    IF(ntest>0) THEN
!    IF(mp > 10300) THEN
      print *,'idjd',idjd
      PRINT 10, ktau,mp, ssoil%isflag(idjd),ssoil%osnowd(idjd),ssoil%snowd(idjd), &
                      met%tk(idjd),met%tc(idjd),met%qv(idjd),met%ua(idjd), &
                      met%fsd(idjd,3),met%fld(idjd), met%precip(idjd),canopy%precis(idjd)
10    format(x,'soilsnowv befstempv,ktau=',2i6,2x,i2,2f8.2,f5.0,f5.1,f7.4,f5.1,x,2f6.0,2f6.3)
      PRINT 11, ktau,canopy%dgdtg(idjd),canopy%fevc(idjd),canopy%fevw(idjd), &
                canopy%fess(idjd),canopy%fe(idjd),canopy%fhs(idjd),canopy%fhv(idjd), &
                canopy%fh(idjd),canopy%fnv(idjd),canopy%fns(idjd),canopy%precis(idjd)
11    format(x,'bsoilsnow',i6,10f7.2,2x,f6.3)
      PRINT 12, ktau,canopy%ga(idjd),(ssoil%ssdn(idjd,k),k=1,3), &
                (ssoil%sdepth(idjd,k),k=1,3),(ssoil%smass(idjd,k),k=1,3)
12    format(x,'bga,ssdn,sdepth smass',i6,f7.1,3f6.1,2x,3f8.3,2x,3f8.1)
      PRINT 13, ktau,(ssoil%tgg(idjd,k),k=1,ms),(ssoil%tggsn(idjd,k),k=1,3)
13    format(x,'btgg tggsn ',i6,6f7.2,2x,3f7.2)
      PRINT 14, ktau,(ssoil%wb(idjd,k),k=1,ms),(ssoil%wbice(idjd,k),k=1,ms),soil%froot(idjd,:)
14    format(x,'bwb wbice',i6,6f6.3,2x,6f6.3,6f6.3)
    ENDIF

    zsetot = sum(soil%zse) 
    ssoil%tggav = 0.
    DO k = 1, ms
     ssoil%tggav = ssoil%tggav  + soil%zse(k)*ssoil%tgg(:,k)/zsetot
    END DO
    ssoil%fwtop1 = 0.0
    ssoil%fwtop2 = 0.0
    ssoil%fwtop3 = 0.0
    ssoil%runoff = 0.0
    ssoil%rnof1 = 0.0
    ssoil%rnof2 = 0.0
    ssoil%smelt = 0.0 ! initialise snowmelt
    ssoil%osnowd = ssoil%snowd
    ssoil%dtmlt = 0.0
!    ssoil%qasrf =  0.0
!    ssoil%qfsrf =  0.0
!    ssoil%qssrf =  0.0

    IF (ktau_gl .eq. 1) THEN
!      canopy%ga = 0.0
    print *,'ktau_gl',ktau_gl
      canopy%dgdtg = 0.0
      ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
      !     IF (ntest == 3) snmin = .11 ! to force 3-layer snow for testing
      ssoil%wbtot = 0.0
      DO k = 1, ms
        ssoil%wb(:,k)  = min( soil%ssat,max ( ssoil%wb(:,k), soil%swilt ))
        ssoil%wbtot = ssoil%wbtot + ssoil%wb(:,k) * 1000.0 * soil%zse(k)
      END DO
      ssoil%wb(:,ms-2)  = min( soil%ssat,max ( ssoil%wb(:,ms-2), 0.5*(soil%sfc+soil%swilt) ))
      ssoil%wb(:,ms-1)  = min( soil%ssat,max ( ssoil%wb(:,ms-1), 0.8*soil%sfc ))
      ssoil%wb(:,ms)    = min( soil%ssat,max ( ssoil%wb(:,ms),   soil%sfc) )
!      ssoil%tgg(:,ms-1)  = ssoil%tgg(:,ms-1) - 1.
!      ssoil%tgg(:,ms)    = ssoil%tgg(:,ms) - 2.
!
      DO k = 1, ms
        WHERE (ssoil%tgg(:,k) <= tfrz .and. ssoil%wbice(:,k) <= 0.01)
          ssoil%wbice(:,k) = 0.5 * ssoil%wb(:,k)
        END WHERE
        WHERE (ssoil%tgg(:,k) < tfrz)
          ssoil%wbice(:,k) = 0.80 * ssoil%wb(:,k)
        END WHERE
      END DO
      WHERE (soil%isoilm == 9) 
        ssoil%snowd = 1100.00 
        ssoil%osnowd = 1100.00 
        ssoil%tgg(:,1) = ssoil%tgg(:,1) - 1.0
!        ssoil%tgg(:,2) = ssoil%tgg(:,2) - 1.0
!        ssoil%tgg(:,3) = ssoil%tgg(:,3) - 0.0
        ssoil%wb(:,1) = 0.95 * soil%ssat
        ssoil%wb(:,2) = 0.95 * soil%ssat
        ssoil%wb(:,3) = 0.95 * soil%ssat
        ssoil%wb(:,4) = 0.95 * soil%ssat
        ssoil%wb(:,5) = 0.95 * soil%ssat
        ssoil%wb(:,6) = 0.95 * soil%ssat
        ssoil%wbice(:,1) = 0.90 * ssoil%wb(:,1)
        ssoil%wbice(:,2) = 0.90 * ssoil%wb(:,2)
        ssoil%wbice(:,3) = 0.90 * ssoil%wb(:,3)
        ssoil%wbice(:,4) = 0.90 * ssoil%wb(:,4)
        ssoil%wbice(:,5) = 0.90 * ssoil%wb(:,5)
        ssoil%wbice(:,6) = 0.90 * ssoil%wb(:,6)
      ENDWHERE
      xx=soil%css * soil%rhosoil
      ssoil%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
           & + (ssoil%wb(:,1) - ssoil%wbice(:,1) ) * cswat * rhowat &
           & + ssoil%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)
    END IF
    xx=soil%css * soil%rhosoil
    IF (ktau <= 1) ssoil%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
           & + (ssoil%wb(:,1) - ssoil%wbice(:,1) ) * cswat * rhowat &
           & + ssoil%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)

!   PRINT *, 'aft init',ktau,canopy%ga(idjd),ssoil%gammzz(idjd,1)
    DO k = 1, ms ! for stempv
      ssoil%wblf(:,k) = MAX( 0.01_r_2, (ssoil%wb(:,k) - ssoil%wbice(:,k)) ) &
                      & / soil%ssat
      ssoil%wbfice(:,k) = REAL(ssoil%wbice(:,k),r_1) / soil%ssat
    END DO

    CALL snowcheck (dels, ktau, ssoil, soil, met , ktau_gl)
!    CALL snowcheck (dels, ktau_gl, ssoil, soil, met )

    CALL snowdensity (dels, ssoil, soil)

    CALL snow_accum (dels, ktau, canopy, met, ssoil, soil )

    CALL snow_melting (dels, snowmlt, ktau, ssoil, soil )
    ssoil%smelt = snowmlt

    ! adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...

    CALL snowl_adjust(dels, ktau, ssoil, canopy )

    ! Diagnostic block:
    IF (ntest > 0) THEN
      PRINT *, 'in soilsnowv before stempv,  ktau= ',ktau
      PRINT *, 'ga,dels,ssdn ', canopy%ga(idjd), dels, (ssoil%ssdn(idjd,k),k=1,3)
      PRINT *, 'osnowd,snowd,isflag', ssoil%osnowd(idjd), ssoil%snowd(idjd), &
             & ssoil%isflag(idjd)
      PRINT *, 'tgg ', (ssoil%tgg(idjd,k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(idjd,k),k=1,3)
      PRINT *, 'wb ', (ssoil%wb(idjd,k),k=1,ms)
      PRINT *, 'wbice ', (ssoil%wbice(idjd,k),k=1,ms)
      PRINT *, 'wblf ', (ssoil%wblf(idjd,k),k=1,ms)
      PRINT *, 'wbfice ', (ssoil%wbfice(idjd,k),k=1,ms)
    END IF

    CALL stempv(dels, canopy, ssoil, soil)
    ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)

    CALL snow_melting (dels, snowmlt, ktau, ssoil, soil )
    ssoil%smelt = ssoil%smelt + snowmlt

    CALL remove_trans(dels, soil, ssoil, canopy)

    CALL  soilfreeze(dels, soil, ssoil)

    totwet = canopy%precis + ssoil%smelt
    weting = totwet + max(0.,ssoil%pudsto - canopy%fesp/hl*dels) ! total available liquid including puddle
    xxx=soil%ssat - ssoil%wb(:,1)
    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*rhowat,weting) !soil capacity
    xxx=soil%ssat - ssoil%wb(:,2)
    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*rhowat,weting-sinfil1) !soil capacity
    xxx=soil%ssat - ssoil%wb(:,3)
    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*rhowat,weting-sinfil1-sinfil2) !soil capacity
    ssoil%fwtop1 = sinfil1/ dels - canopy%segg
    ssoil%fwtop2 = sinfil2/ dels
    ssoil%fwtop3 = sinfil3/ dels

!   Puddle for the next time step
    ssoil%pudsto = max( 0., weting - sinfil1 - sinfil2 - sinfil3 )
    ssoil%rnof1 = max(0.,ssoil%pudsto - ssoil%pudsmx)
    ssoil%pudsto = ssoil%pudsto - ssoil%rnof1

    CALL surfbv(dels, ktau,  met, ssoil, soil, veg, canopy )

    where( ssoil%dtmlt(:,1) > 0.0001) canopy%fhs = canopy%fhs+ssoil%dtmlt(:,1)*ssoil%dfh_dtg
    where( ssoil%dtmlt(:,1) > 0.0001) canopy%fes = canopy%fes+ssoil%dtmlt(:,1)* &
                                             (ssoil%cls*ssoil%dfe_ddq * ssoil%ddq_dtg)

    CALL hydraulic_redistribution(dels,ktau,soil,ssoil,canopy,veg, met)

    ! Diagnostic block:
    IF (ntest > 0) THEN
!    IF (mp > 10000) THEN
      PRINT 15, ktau, ssoil%isflag(idjd),ssoil%osnowd(idjd),ssoil%snowd(idjd), &
                      ssoil%rnof1(idjd),ssoil%rnof2(idjd)
15    format(x,'after soilsnowv ktau=',i6,2x,i2,2f8.2,2x,2f6.3)
      PRINT 16, ktau,canopy%dgdtg(idjd),canopy%fevc(idjd),canopy%fevw(idjd), &
                canopy%fess(idjd),canopy%fe(idjd),canopy%fhs(idjd),canopy%fhv(idjd), &
                canopy%fh(idjd),canopy%fnv(idjd),canopy%fns(idjd),canopy%precis(idjd)
16    format(x,'absoilsnow',i6,10f7.2,2x,f6.3)
      PRINT 17, ktau,canopy%ga(idjd),(ssoil%ssdn(idjd,k),k=1,3), &
                (ssoil%sdepth(idjd,k),k=1,3),(ssoil%smass(idjd,k),k=1,3)
17    format(x,'aft ga,ssdn,sdepth smass',i6,f7.1,3f6.1,2x,3f8.3,2x,3f8.1)
      PRINT 18, ktau,(ssoil%tgg(idjd,k),k=1,ms),(ssoil%tggsn(idjd,k),k=1,3)
18    format(x,'aftgg tggsn ',i6,6f7.2,2x,3f7.2)
      PRINT 19, ktau,(ssoil%wb(idjd,k),k=1,ms),(ssoil%wbice(idjd,k),k=1,ms), &
                (ssoil%wblf(idjd,k),k=1,ms),(ssoil%wbfice(idjd,k),k=1,ms)
19    format(x,'awb wbice',i6,6f6.3,x,6f6.3,x,6f6.3,x,6f6.3)
    END IF
    ! Set weighted soil/snow surface temperature:

    ssoil%tss=(1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)
    ssoil%smelt = ssoil%smelt/dels
  
    ssoil%wbtot = 0.0
    DO k = 1, ms
       ssoil%wbtot = ssoil%wbtot + ssoil%wb(:,k) * 1000.0 * soil%zse(k)
    END DO

  END SUBROUTINE soil_snow

  !+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++
  ! Sciences from Ryel et al. Oecologia, 2002; Lee et al., 2005, PNAS
  ! Code by LiLH 16 Feb, 2011
  ! Fixed problem of negative wb in global run by BP Mar 2011
  SUBROUTINE hydraulic_redistribution(dels, ktau, soil, ssoil, canopy, veg, met)
    REAL(r_1),                 INTENT(IN) :: dels ! integration time step (s)
    INTEGER(i_d),              INTENT(IN) :: ktau  ! integration step number
    TYPE(soil_parameter_type), INTENT(IN) :: soil
    TYPE(soil_snow_type),   INTENT(INOUT) :: ssoil
    TYPE(canopy_type),         INTENT(IN) :: canopy
    TYPE(veg_parameter_type),  INTENT(IN) :: veg
    TYPE(met_type), INTENT(INOUT)         :: met ! all met forcing
    INTEGER(i_d) k
    INTEGER(i_d) j
    INTEGER(i_d) ii
    REAL(r_1), DIMENSION(mp,ms)    :: S_VG                ! --
    REAL(r_1), DIMENSION(mp,ms)    :: wpsy                ! MPa
    REAL(r_1), DIMENSION(mp)       :: frootX              ! --
    REAL(r_1), DIMENSION(mp,ms)    :: C_hr                ! --
!    REAL(r_1), DIMENSION(mp,ms)    :: hr                  ! cm/hour
    REAL(r_1), DIMENSION(mp,ms,ms) :: hr_term             ! cm/hour
    REAL(r_1), DIMENSION(mp)       :: Dtran               ! Swith for hr

    REAL(r_1), PARAMETER :: thetas=0.45  ! from Belk et al., 2007, WRR
    REAL(r_1), PARAMETER :: thetar=0.20  ! from Belk et al., 2007, WRR
    ! REAL(r_1), PARAMETER :: alpha_VG = 0.00045  ! from Belk et al., 2007, WRR,
                                                  ! cm^{-1} 1cmH2O=100Pa
    ! REAL(r_1), PARAMETER :: n_VG = 1.40         ! from Belk et al., 2007, WRR
    REAL(r_1), PARAMETER :: n_VG = 2.06           ! -- 2.06
    REAL(r_1), PARAMETER :: m_VG = 1.0-1.0/n_VG   ! --
    REAL(r_1), PARAMETER :: alpha_VG = 0.00423    ! cm^{-1} Note: 1cmH2O=100Pa
    REAL(r_1), PARAMETER :: n_hr = 3.22           ! --
    REAL(r_1), PARAMETER :: wpsy50 = -1.0         ! MPa
    REAL(r_1), PARAMETER :: CRT = 125.0           ! cm MPa^-1 h^-1, default value (0.097) from Ryel et al., 2002
    REAL(r_1), PARAMETER :: wiltParam = 0.5   
    REAL(r_1), PARAMETER :: satuParam = 0.8   
    ! 125.0  for Tumbarumbar 
    ! 125.0  for Dinghushan
    ! 125.0  for Howard Springs
    ! 125.0  for KM83 
    REAL(r_1), DIMENSION(mp,ms,ms) :: hr_perTime
    REAL(r_1), DIMENSION(mp)       :: temp
    REAL(r_1), DIMENSION(mp)       :: available
    REAL(r_1), DIMENSION(mp)       :: accommodate
    REAL(r_1), DIMENSION(mp)       :: totalmoist,totalice
    REAL(r_1), DIMENSION(mp)       :: total2,zsetot
!    REAL(r_1), PARAMETER :: wiltParam = 0.5
!    REAL(r_1), PARAMETER :: satuParam = 0.8
    INTEGER(i_d) :: idjd1

    ! Find initial moisture total for checking purpose
    idjd1 = 351  ! Amazon

    zsetot = sum(soil%zse)
    totalmoist(:) = 0.0
    totalice(:) = 0.0
    DO k=1, ms
      totalmoist(:) = totalmoist(:) + ssoil%wb(:,k)*soil%zse(k)/zsetot
      totalice(:) = totalice(:) + ssoil%wbice(:,k)*soil%zse(k)/zsetot
    ENDDO

    WHERE ( canopy%fevc > 10.0 )  ! restrict to nighttime and cloudy day
      Dtran=0.0
    ELSEWHERE
      Dtran=1.0
    ENDWHERE
    where ( totalice  > 1.e-2 ) Dtran=0.
    
    DO k=1, ms
      S_VG(:,k) = MIN( 1.0, MAX(1.0E-4, ssoil%wb(:,k) - soil%swilt) &
                             / (soil%ssat - soil%swilt) )
      wpsy(:,k) = -1.0/alpha_VG*(S_VG(:,k)**(-1.0/m_VG)-1.0)**(1/n_VG) &
                   *100*1.0E-6     ! VG model, convert from cm to Pa by (*100),
                                   !                           to MPa (*1.0E-6)
      C_hr(:,k) = 1./(1+(wpsy(:,k)/wpsy50)**n_hr)
    ENDDO
!    print *,'zse',soil%zse,ms,mp,soil%swilt(idjd1),soil%ssat(idjd1),soil%sfc(idjd1), &
!                  met%tk(idjd1),met%precip(idjd1)
!    print 14,(ssoil%wb(idjd1,k),k=1,4),totalmoist(idjd1),(wpsy(idjd1,k),k=1,4), &
!                       (C_hr(idjd1,k),k=1,4),canopy%fevc(idjd1),Dtran(idjd1)
!14  format(1x,'hydrre0',4f8.5,2x,f8.5,2(4f8.3),2x,f8.3,f4.1)

!    hr(:,:) = 0.0
    temp(:)        = 0.0
    hr_term(:,:,:) = 0.0    ! unit: cm h^{-1}
    hr_perTime(:,:,:) = 0.0
    ! setting hr_term=0 for top layer, follows Lee et al., 2005, PNAS
    DO k = ms, 3, -1
      DO j = k-1, 2, -1
        temp(:)        = 0.0
        available(:)   = 0.0
        accommodate(:) = 0.0
        WHERE ( soil%froot(:,k) > soil%froot(:,j) )
          frootX=soil%froot(:,k)
        ELSEWHERE
          frootX=soil%froot(:,j)
        ENDWHERE
        hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
                       *(soil%froot(:,k)*soil%froot(:,j))/(1-frootX) * Dtran
        hr_perTime(:,k,j) = hr_term(:,k,j)*1.0E-2/3600.0*dels ! m per timestep
        hr_perTime(:,j,k) = -1.0 * hr_perTime(:,k,j)
        hr_perTime(:,k,j) = hr_perTime(:,k,j)/soil%zse(k)
        hr_perTime(:,j,k) = hr_perTime(:,j,k)/soil%zse(j)
        ! Restricting changes to all broadleaf forests, and
        ! other forests and woody savannas in the tropics
        ! Note that veg types here are based on IGBP classification (BP mar2011)
!        WHERE (.NOT.(veg%iveg==2.OR.veg%iveg==7.OR.veg%iveg==8.OR.veg%iveg==9))
        WHERE (.NOT.veg%iveg==2)
!        WHERE (.NOT.(veg%iveg == 2 .OR. veg%iveg == 2 ))
!                     .OR.((veg%iveg <=5 .OR. veg%iveg ==8) )))
!                     .OR.((veg%iveg <=5 .OR. veg%iveg ==8) &
!                          .AND. patch(:)%latitude > -24.0 &
!                          .AND. patch(:)%latitude <  24.0)  ))
          hr_perTime(:,k,j) = 0.0
          hr_perTime(:,j,k) = 0.0
        ENDWHERE
        WHERE (hr_perTime(:,k,j) < 0.0)
!          available(:)   = MAX(0.0, ssoil%wb(:,k)-soil%swilt(:))
          available(:)   = MAX(0.0, ssoil%wb(:,k)-  &
                          ( soil%swilt(:) + (soil%sfc(:)-soil%swilt(:))/3.) )
          accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,j))
          temp(:) = MAX(hr_perTime(:,k,j), &
                        -1.0*wiltParam*available(:), &
                        -1.0*satuParam*accommodate(:)*soil%zse(j)/soil%zse(k)) 
          hr_perTime(:,k,j) = temp(:)
          hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)
        ELSEWHERE (hr_perTime(:,j,k) < 0.0)
!          available(:)   = MAX(0.0, ssoil%wb(:,j)-soil%swilt(:))
          available(:)   = MAX(0.0, ssoil%wb(:,j)-  &
                           ( soil%swilt(:) + (soil%sfc(:)-soil%swilt(:))/3.) )
          accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,k))
          temp(:) = MAX(hr_perTime(:,j,k), &
                        -1.0*wiltParam*available(:), &
                        -1.0*satuParam*accommodate(:)*soil%zse(k)/soil%zse(j))
          hr_perTime(:,j,k) = temp(:)
          hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)
        ENDWHERE
        ssoil%wb(:,k) = ssoil%wb(:,k) + hr_perTime(:,k,j)
        ssoil%wb(:,j) = ssoil%wb(:,j) + hr_perTime(:,j,k)
!      print 18,k,j,veg%iveg(idjd1),available(idjd1),accommodate(idjd1),hr_perTime(idjd1,j,k)
!18  format(1x,'hydrre1',3i3,3f8.4)
      ENDDO 
    ENDDO

    ! Redistribution should not change the total. Check if it is true.
!    total2(:) = 0.0
!    DO k=1, ms
!      total2(:) = total2(:) + ssoil%wb(:,k)*soil%zse(k)/zsetot
!    ENDDO
!    print 20,(ssoil%wb(idjd1,k),k=1,4),totalmoist(idjd1),total2(idjd1)
!20  format(1x,'hydrre2',4f8.5,2x,2f8.5)


!    IF (MINVAL(totalmoist(:) - total2(:)) < -1.0e-6 .OR. &
!        MAXVAL(totalmoist(:) - total2(:)) >  1.0e-6) THEN
!      PRINT *, 'Problem with hydraulic_redistribution at timestep ', ktau
!      PRINT *, 'Maximum of ', MAXVAL(totalmoist(:) - total2(:)), ' at ', &
!                MAXLOC(totalmoist(:) - total2(:))
!      PRINT *, 'Minimum of ', MINVAL(totalmoist(:) - total2(:)), ' at ', &
!                MINLOC(totalmoist(:) - total2(:))
!      STOP
!    ENDIF

  END SUBROUTINE hydraulic_redistribution
  !+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++


END MODULE soil_snow_module
