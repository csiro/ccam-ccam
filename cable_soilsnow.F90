
MODULE soil_snow_module
  USE physical_constants
  USE define_dimensions       
  USE define_types       
  USE io_variables, ONLY: landpt
  IMPLICIT NONE
  PRIVATE
  REAL(r_1), PARAMETER :: cgsnow = 2090.0 ! specific heat capacity for snow
  REAL(r_1), PARAMETER :: csice = 2.100e3 ! specific heat capacity for ice
  REAL(r_1), PARAMETER :: cswat = 4.218e3 ! specific heat capacity for water
  REAL(r_1), PARAMETER :: cp = capp       ! specific heat capacity for air
  REAL(r_1), PARAMETER :: rhowat = 1000.0 ! density of water
   !jhan:Eva uses .11 -> 1.
  REAL(r_1), PARAMETER :: snmin = 1.    ! for 3-layer;
  REAL(r_1), PARAMETER :: max_ssdn = 750.0
  REAL(r_1), PARAMETER :: max_sconds = 2.51
   !jhn::r2411 uses .85
  REAL(r_1), PARAMETER :: frozen_limit = 0.85 ! EAK Feb2011 (could be 0.95)
    real :: max_glacier_snowd

  ! This module contains the following subroutines:
  PUBLIC soil_snow ! must be available outside this module
  PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust 
   PRIVATE trimb, smoisturev, snow_accum, stempv
!   PRIVATE soilfreeze, remove_trans, snow_albedo
   PRIVATE soilfreeze, remove_trans

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
    INTEGER(i_d), INTENT(IN)                  :: kmax ! no. of discrete layers    
    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: a ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: b ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: c ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(INOUT)  :: rhs ! right hand side of eq
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
  ! SUBROUTINE smoisturev (fwtop,dels,ssoil,soil)
  !      Solves implicit soil moisture equation
  !      Science development by Eva Kowalczyk and John McGregor, CMAR
  !
  SUBROUTINE smoisturev (dels,ssoil,soil,veg)
   use cable_common_module
    REAL(r_1), INTENT(IN)                     :: dels    ! time step size (s)
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    TYPE(veg_parameter_type), INTENT(IN)  :: veg

    INTEGER(i_d), PARAMETER                   :: ntest = 0 ! 2 for funny pre-set
    ! nmeth selects the solution method
    INTEGER(i_d), PARAMETER                   :: nmeth = -1 ! preferred method
    !                                  Values as follows:
    !                                   -1 for simple implicit D
    !                                    1 for fully implicit solution
    !                                    2 for simpler implicit
    !                                    3 for simple implicit D, explicit K 
    !                                    4 for simple implicit D, implicit K
    !                                    0 for simple implicit D, new jlm TVD K
! change dimension of at,bt,ct from 3*ms to ms (BP Jun2010)
    REAL(r_2), DIMENSION(mp,ms)   :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms)   :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms)   :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp)      :: fact
    REAL(r_2), DIMENSION(mp)      :: fact2
    REAL(r_2), DIMENSION(mp)      :: fluxhi
    REAL(r_2), DIMENSION(mp)      :: fluxlo
    REAL(r_2), DIMENSION(mp)      :: hydss  ! hydraulic conductivity
                                            ! adjusted for ice
    INTEGER(i_d)                  :: k
    REAL(r_2), DIMENSION(mp)      :: phi
    REAL(r_2), DIMENSION(mp)      :: pwb
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
    INTEGER(i_d) :: u           ! I/O unit

!jhan:peculiar to offline
!jh   if( cable_runtime%offline .or. cable_runtime%mk3l ) then
!jh      soil%pwb_min = (soil%swilt / soil%ssat ) **soil%ibp2
!jh   endif

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
          where(ssoil%wbice(:,k) > 0.05.or.ssoil%wbice(:,k+1) > 0.01) &
             wh = 0.9*wbl_k + 0.1*wbl_kp
          ! with 50% wbice, reduce hyds by 1.e-5
          ! Calculate hyd conductivity adjusted for ice:
          hydss = soil%hyds
!          hydss = soil%hyds * (1.0 - MIN(2.0_r_2 * ssoil%wbice(:,k) &
!              & / MAX(0.01_r_2,  ssoil%wb(:,k) ), 0.99999_r_2) )
          speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
          ! update wb by TVD method
          rat = delt(:,k - 1) / (delt(:,k)+SIGN(REAL(1.0e-20,r_2), delt(:,k)))
          phi = MAX(0.0_r_2, MIN(1.0_r_2, 2.0_r_2 * rat), &
               MIN(2.0_r_2, rat) ) ! 0 for -ve rat
          fluxhi = wh
          fluxlo = wbl_k
          ! scale speed to grid lengths per dels & limit speed for stability
          ! 1. OK too for stability
          speed_k = MIN(speed_k, REAL(0.5 * soil%zse(k) / dels , r_2))
          fluxh(:,k) = speed_k * (fluxlo + phi * (fluxhi - fluxlo) )
      END DO
      ! calculate drainage (this code replaces the code in the surfb)
      k = ms 
!jhan:ifdef ONLINE_UM
      WHERE( ssoil%wb(:,ms) > soil%sfc(:) )
      !jh:WHERE( soil%albsoil(:,1) .lt. 0.30 )
        wbl_k = MAX(0.001_r_2, ssoil%wb(:,ms) - ssoil%wbice(:,ms) )
        wbl_kp = MAX(0.001_r_2, soil%ssat(:) - ssoil%wbice(:,ms) )
        !wh = 0.9*wbl_k + 0.1*wbl_kp
        wh = MIN(wbl_k, wbl_kp)
        where(ssoil%wbice(:,ms).gt. 0.05) wh = 0.9*wbl_k + 0.1*wbl_kp
       
        ! Calculate hyd conductivity adjusted for ice:
        hydss = soil%hyds
        !hydss = soil%hyds * ( 1. - MIN( 2.0_r_2 * ssoil%wbice(:,ms) &
        !      & / MAX(0.01_r_2, ssoil%wb(:,ms)), 0.99999_r_2 ) )
        speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
        speed_k =  0.5*speed_k / (1. - min(0.5,10.*ssoil%wbice(:,ms)))
        fluxlo = wbl_k
        ! scale speed to grid lengths per dt & limit speed for stability
        speed_k = MIN(0.5*speed_k, 0.5 * soil%zse(ms) / dels)
        fluxh(:,ms) = MAX(0.0,speed_k * fluxlo )
     END WHERE

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
        pwb_wbh = (soil%hsbh * (1. - MIN(2. * MIN(0.1_r_2, MAX( &
                & ssoil%wbice(:,k-1) / MAX(0.01_r_2, ssoil%wb(:,k-1)), &
                & ssoil%wbice(:,k)   / MAX(0.01_r_2, ssoil%wb(:,k)  ) )) &
                & , 0.1_r_2) )) &
                & * MAX(soil%pwb_min, wbh_k * fact)
        ! moisture diffusivity (D) is  wbh*pwb; hsbh includes b
        ! i.e. D(k-.5)/soil%zshh(k)
        z3_k = pwb_wbh / soil%zshh (k)

        ! where dtt=dels/(soil%zse(k)*ssatcurr_k)
        at (:,k) = - dtt(:,k) * z3_k
        ct (:,k-1) = - dtt(:,k-1) * z3_k
      END DO
      bt = 1. - at - ct
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
        !jhan: this has already been done ?
        ssoil%wblf(:,k) = (ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr(:,k)
        !jhan: this has already been done?
        ssoil%wbfice(:,k) = REAL(ssoil%wbice(:,k),r_1) / soil%ssat
        wbficemx = MAX(wbficemx, ssoil%wbfice(:,k) )
        dtt(:,k) = dels / (soil%zse(k) * ssatcurr(:,k) )
     END DO

     IF (nmeth == 1) THEN ! full implicit method
        DO k = 2, ms
           wbh(:,k) = (soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1) &
                * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k-1) )
           fact = wbh(:,k) ** (soil%ibp2 - 1) ! i.e. wbh**(bch+1)
           fact2 = fact * fact
           pwb = soil%hsbh * fact
           ! moisture diffusivity (D) is  wbh*pwb
           ! other term (K) is wbh*soil%hyds*fact2
           z1(:,k) = wbh(:,k) * ( (soil%i2bp3 - 1) * soil%hyds * fact2 &
                - soil%ibp2 * pwb * (ssoil%wblf(:,k) - ssoil%wblf(:,k-1) ) &
                / soil%zshh (k) )
           z2(:,k) = - soil%i2bp3 * soil%hyds * fact2 + soil%ibp2 * pwb &
                * (ssoil%wblf(:,k) - ssoil%wblf(:,k-1) ) / soil%zshh (k)
           z3(:,k) = pwb * wbh(:,k) / soil%zshh (k)
           at(:,k) = dtt(:,k) * (z2(:,k) * 0.5 * soil%zse(k) / soil%zshh (k) &
                - z3(:,k) )
        END DO
        DO k = 1, ms - 1
           ct(:,k) = dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k) &
                / soil%zshh (k+1) - z3(:,k+1) )
           bt(:,k) = 1.0 + dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k+1) &
                / soil%zshh (k+1) + z2(:,k) * 0.5 * soil%zse(MAX(k-1,1)) &
                / soil%zshh (k) + z3(:,k+1) + z3(:,k) )
        END DO
        bt(:,ms) = 1.0 + dtt(:,ms) * (z2(:,ms) * 0.5 * soil%zse(ms) &
                  & / soil%zshh (ms) + z3(:,ms) )
        DO k = 1, ms
          ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) * (z1(:,k+1) - z1(:,k) )
        END DO
     END IF ! (nmeth == 1)
     
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
           bt(:,1) = bt(:,1) + dtt(:,1) * z1mult(:,1+1) &
                * z1(:,1+1) * soil%zse(1+1) / (soil%zse(1) + soil%zse(1+1) )
           DO k = 2, ms
              at(:,k)   = at(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k) &
                   * soil%zse(k) / (soil%zse(k) + soil%zse(k-1) )
              ct(:,k-1) = ct(:,k-1) + dtt(:,k-1) * z1mult(:,k) * z1(:,k) &
                   * soil%zse(k-1) / (soil%zse(k) + soil%zse(k-1) )
              bt(:,k)   = bt(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k) &
                   * soil%zse(k-1) / (soil%zse(k) + soil%zse(k-1) ) &
                   + dtt(:,k) * z1mult(:,k+1) * z1(:,k+1) &
                   * soil%zse(k+1) / (soil%zse(k) + soil%zse(k+1) )
           END DO
        END IF ! (nmeth == 4)
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
                      * ( (z1mult(:,k+1) - 1.0) * z1(:,k+1) &
                      - (z1mult(:,k) - 1.0) * z1(:,k) )
              END WHERE
           ELSE
              WHERE (wbficemx < 0.75)
                 ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) &
                      * (z1(:,k) - z1(:,k+1) )
              END WHERE
           END IF
        END DO
      END IF

!jhan: peculiar to UM      
!jh   if( cable_runtime%um) then
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
      END IF
!jh   END IF

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
      ssoil%wbice(:,k) = MIN(ssoil%wbice(:,k), frozen_limit * ssoil%wb(:,k) )
    END DO

!jhan: peculiar to UM      
!jh   if( cable_runtime%um) then
    IF (ntest > 0) THEN
      totwbc = 0.
      totwblc = 0.
      DO k = 1, ms
        totwbc = totwbc + soil%zse(k) * REAL(ssoil%wb(:,k),r_1)
        totwblc = totwblc + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1)
      END DO
    END IF
!jh    END IF

  END SUBROUTINE smoisturev


  !-------------------------------------------------------------------------

 SUBROUTINE snowdensity (dels, ssoil, soil)
    REAL(r_1), INTENT(IN)   :: dels   ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    WHERE (ssoil%snowd > 0.1 .and. ssoil%isflag == 0)
      ssoil%ssdn(:,1) = MIN(max_ssdn,MAX(120.0, ssoil%ssdn(:,1) + dels &
          & * ssoil%ssdn(:,1) * 3.1e-6 * EXP( -0.03 &
          & * (273.15 - MIN(tfrz, ssoil%tgg(:,1) )) &
          & - merge(0.046, 0.0, ssoil%ssdn(:,1) >= 150.0) &
          & * (ssoil%ssdn(:,1) - 150.0) ) ))
      ssoil%ssdn(:,1) = MIN(max_ssdn,ssoil%ssdn(:,1) + dels * 9.806 &
          & * ssoil%ssdn(:,1) * 0.75 * ssoil%snowd &
          & / (3.0e7 * EXP(0.021 * ssoil%ssdn(:,1) + 0.081 &
          & * (273.15 - MIN(tfrz, ssoil%tgg(:,1))))))
!jhan: peculiar to UM      
!jh   if( cable_runtime%um) then
      WHERE (soil%isoilm /= 9) ssoil%ssdn(:,1) = MIN(450.0,ssoil%ssdn(:,1))
!jh   endif
      ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
                                                        & + 0.074, max_sconds) )
!                                                        & + 0.074, 1.0) )
      ssoil%sconds(:,2) = ssoil%sconds(:,1) 
      ssoil%sconds(:,3) = ssoil%sconds(:,1) 
      ssoil%ssdnn = ssoil%ssdn(:,1)
      ssoil%ssdn(:,2) = ssoil%ssdn(:,1)
      ssoil%ssdn(:,3) = ssoil%ssdn(:,1)
    END WHERE
    WHERE (ssoil%isflag == 1)
       ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * ssoil%ssdn(:,1) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,1))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,1) >= 150.0) &
            * (ssoil%ssdn(:,1) - 150.0) )
       ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * ssoil%ssdn(:,2) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,2))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,2) >= 150.0) &
            * (ssoil%ssdn(:,2) - 150.0) )
       ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * ssoil%ssdn(:,3) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,3))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,3) >= 150.0) &
            * (ssoil%ssdn(:,3) - 150.0) )
       ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * 9.806 * ssoil%ssdn(:,1) &
            * ssoil%t_snwlr*ssoil%ssdn(:,1) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,1) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,1)))))
       ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * 9.806 * ssoil%ssdn(:,2) &
            * (ssoil%t_snwlr * ssoil%ssdn(:,1) + 0.5 * ssoil%smass(:,2) ) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,2) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,2)))))
       ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * 9.806 * ssoil%ssdn(:,3) &
            * (ssoil%t_snwlr*ssoil%ssdn(:,1) + ssoil%smass(:,2) &
            + 0.5*ssoil%smass(:,3)) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,3) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,3)))))
      ssoil%sdepth(:,1) =  ssoil%smass(:,1) / ssoil%ssdn(:,1) 
      ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
      ssoil%sdepth(:,3) =  ssoil%smass(:,3) / ssoil%ssdn(:,3) 
      ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%smass(:,1) + ssoil%ssdn(:,2) &
            * ssoil%smass(:,2) + ssoil%ssdn(:,3) * ssoil%smass(:,3) ) &
            / ssoil%snowd
      ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
                                                        & + 0.074, max_sconds) )
      ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2) ** 2 &
                                                        & + 0.074, max_sconds) )
      ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3) ** 2 &
                                                        & + 0.074, max_sconds) )
    END WHERE
  END SUBROUTINE snowdensity

  !-------------------------------------------------------------------------

  SUBROUTINE snow_melting (dels, snowmlt, ssoil, soil )
      use cable_common_module
    REAL(r_1), INTENT(IN)                 :: dels   ! integration time step (s)
    REAL(r_1), DIMENSION(mp), INTENT(OUT) :: snowmlt ! snow melt   
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)   :: ssoil  ! soil+snow variables
    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
    INTEGER(i_d)                 :: k,j 
    REAL(r_1), DIMENSION(mp)     :: osm
    REAL(r_1), DIMENSION(mp)     :: sgamm
    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
    REAL(r_1), DIMENSION(mp)     :: snowflx

    snowmlt= 0.0
    smelt1 = 0.0
    
    do j=1,mp  
      if (ssoil%snowd(j) > 0.0 .AND. ssoil%isflag(j) == 0 &
          .AND. ssoil%tgg(j,1) >= tfrz ) then
        ! snow covered land
        ! following done in sflux  via  ga= ... +cls*egg + ...
        ! ** land,snow,melting
        snowflx(j) = REAL((ssoil%tgg(j,1) - tfrz) * ssoil%gammzz(j,1),r_1)
        ! prevent snow depth going negative
        snowmlt(j) = MIN(snowflx(j) / hlf, ssoil%snowd(j) )
       !jhan:Eva adds this
       ssoil%dtmlt(j,1) = ssoil%dtmlt(j,1) + snowmlt(j) * hlf / ssoil%gammzz(j,1)
        ssoil%snowd(j) = ssoil%snowd(j) - snowmlt(j)
        ssoil%tgg(j,1) = &
             & REAL(ssoil%tgg(j,1) - snowmlt(j) * hlf / ssoil%gammzz(j,1),r_1)
      ENdif
    end do

    smelt1(:,0) = 0.0
    DO k = 1, 3
      WHERE (ssoil%snowd > 0.0 .and. ssoil%isflag > 0)
        sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
        ! snow melt refreezing
        snowflx = smelt1(:,k-1) * hlf / dels
!jhan: diff. calc. runtime IF not allowed in WHERE loop 
!jh   if( cable_runtime%um) then
        ssoil%tggsn(:,k) = ssoil%tggsn(:,k) + ( snowflx * dels +  &
                         & smelt1(:,k-1)*cswat*(tfrz-ssoil%tggsn(:,k)) ) / &
                         & ( sgamm + cswat*smelt1(:,k-1) )
!jh   else       
!jh          ssoil%tggsn(:,k) = ssoil%tggsn(:,k) + snowflx * dels / sgamm
!jh   endif
        ! increase density due to snowmelt
        osm = ssoil%smass(:,k)
        ssoil%smass(:,k) = ssoil%smass(:,k) + smelt1(:,k-1)
        ssoil%ssdn(:,k) = MAX(120.0,MIN(ssoil%ssdn(:,k) * osm/ssoil%smass(:,k) &
                        & + rhowat*(1.0-osm/ssoil%smass(:,k)), max_ssdn))
!jhan: peculiar to UM      
        WHERE (soil%isoilm /= 9) ssoil%ssdn(:,k) = MIN(450.0,ssoil%ssdn(:,k))
        ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
        sgamm = ssoil%smass(:,k) * cgsnow
        !jh:smelt1(:,k-1) = 0.0
        smelt1(:,k) = 0.0
        ! snow melting
        WHERE (ssoil%tggsn(:,k) > tfrz)
          snowflx = (ssoil%tggsn(:,k) - tfrz) * sgamm
          !jhan:Eva uses
          smelt1(:,k) = MIN(snowflx / hlf, 0.6 * ssoil%smass(:,k) )
          !jh:smelt1(:,k) = MIN(snowflx / hlf, 0.9 * ssoil%smass(:,k) )
          ssoil%dtmlt(:,k) = ssoil%dtmlt(:,k) + smelt1(:,k) * hlf / sgamm
          osm = ssoil%smass(:,k)
          ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
          !jhan:Eva uses
          ssoil%tggsn(:,k) = ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm
          !jh:ssoil%tggsn(:,k) = MIN(ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm, &
           !jh:                    & tfrz)
          ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
        END WHERE
      END WHERE
    END DO
    WHERE (ssoil%snowd > 0.0 .and. ssoil%isflag > 0)
      snowmlt = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
      ssoil%snowd = ssoil%snowd - snowmlt
    END WHERE

  END SUBROUTINE snow_melting


  !-------------------------------------------------------------------------
  SUBROUTINE snow_accum (dels,  canopy, met, ssoil, soil)
   use cable_common_module
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy ! vegetation variables
    TYPE(met_type), INTENT(INOUT)            :: met   ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    INTEGER(i_d), PARAMETER  :: ntest = 0 ! for snow diag prints
    INTEGER(i_d)             :: k
    REAL(r_1), DIMENSION(mp) :: osm
    REAL(r_1), DIMENSION(mp) :: sgamm
    REAL(r_1), DIMENSION(mp) :: snowmlt
    REAL(r_1), DIMENSION(mp) :: xxx
    real, dimension(mp) :: swilt_eff

    WHERE (canopy%precis > 0.0 .and. ssoil%isflag == 0)
      ssoil%snowd = MAX(ssoil%snowd + met%precip_sn, 0.0) ! accumulate solid part
      canopy%precis = canopy%precis - met%precip_sn
      ssoil%ssdn(:,1) = MAX(120.0, ssoil%ssdn(:,1) &
                      & * ssoil%osnowd / MAX(0.01, ssoil%snowd) &
                      & + 120.0 * met%precip_sn / MAX(0.01, ssoil%snowd))
      ssoil%ssdnn = ssoil%ssdn(:,1)
      WHERE (canopy%precis > 0.0 .and. ssoil%tgg(:,1) < tfrz )
        ssoil%snowd = MAX(ssoil%snowd + canopy%precis, 0.0)
        ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis * hlf &
                   & / ( REAL(ssoil%gammzz(:,1),r_1) + cswat*canopy%precis )  
!                       & / REAL(ssoil%gammzz(:,1),r_1) + &
!jh:        ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis *  &
!jh:              & ( hlf + cswat*(tfrz-ssoil%tgg(:,1)) ) /  &
!jh:              & ( REAL(ssoil%gammzz(:,1),r_1) + cswat*canopy%precis )  
        ! change density due to water being added 
        ssoil%ssdn(:,1) = MIN(max_ssdn, MAX(120.0, ssoil%ssdn(:,1) &
                        & * ssoil%osnowd / MAX(0.01, ssoil%snowd) &
                        & + rhowat * canopy%precis / MAX(0.01, ssoil%snowd)))
        WHERE (soil%isoilm /= 9) ssoil%ssdn(:,1) = MIN(450.0,ssoil%ssdn(:,1))
          canopy%precis = 0.0
          ssoil%ssdnn = ssoil%ssdn(:,1)
      END WHERE
    END WHERE ! (canopy%precis > 0. .and. ssoil%isflag == 0) 

    WHERE (canopy%precis > 0.0 .and.  ssoil%isflag > 0)
      ! add solid precip
      ssoil%snowd = MAX(ssoil%snowd + met%precip_sn, 0.0)
      canopy%precis = canopy%precis - met%precip_sn  ! remaining liquid precip
      ! update top snow layer with fresh snow
      osm = ssoil%smass(:,1)
      ssoil%smass(:,1) = ssoil%smass(:,1) + met%precip_sn
      ssoil%ssdn(:,1) = MAX(120.0,ssoil%ssdn(:,1)*osm/ssoil%smass(:,1) &
                      & + 120.0 * met%precip_sn/ssoil%smass(:,1))
      ssoil%sdepth(:,1) = MAX(0.02,ssoil%smass(:,1) / ssoil%ssdn(:,1))
      ! add liquid precip
       WHERE (canopy%precis > 0.0 )
         ssoil%snowd = MAX(ssoil%snowd + canopy%precis, 0.0)
         sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
         osm = ssoil%smass(:,1)
!        ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + canopy%precis *  &
!                         & ( hlf + cswat*(tfrz-ssoil%tggsn(:,1)) )* osm/ssoil%osnowd  &
!                         & / (sgamm + cswat*canopy%precis*osm/ssoil%osnowd)
         ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + canopy%precis * hlf &
                            * osm / (sgamm * ssoil%osnowd )
         ssoil%smass(:,1) = ssoil%smass(:,1) + canopy%precis &
                            * osm/ssoil%osnowd
         ssoil%ssdn(:,1) = MAX(120.0,MIN(ssoil%ssdn(:,1) * osm/ssoil%smass(:,1) &
                        & +  rhowat*(1.0-osm/ssoil%smass(:,1)), max_ssdn))
         WHERE (soil%isoilm /= 9) ssoil%ssdn(:,1) = MIN(450.0,ssoil%ssdn(:,1))
         ssoil%sdepth(:,1) = ssoil%smass(:,1)/ssoil%ssdn(:,1)
! layer 2
         sgamm = ssoil%ssdn(:,2) * cgsnow * ssoil%sdepth(:,2)
         osm = ssoil%smass(:,2)
!        ssoil%tggsn(:,2) = ssoil%tggsn(:,2) + canopy%precis * &
!                         & ( hlf + cswat*(tfrz-ssoil%tggsn(:,2)) )* osm/ssoil%osnowd  &
!                         & / (sgamm + cswat*canopy%precis*osm/ssoil%osnowd)
         ssoil%tggsn(:,2) = ssoil%tggsn(:,2) + canopy%precis * hlf &
                            * osm / (sgamm * ssoil%osnowd )
         ssoil%smass(:,2) = ssoil%smass(:,2) + canopy%precis &
                            & * osm/ssoil%osnowd
         ssoil%ssdn(:,2) = MAX(120.0,MIN(ssoil%ssdn(:,2) * osm/ssoil%smass(:,2) &
                        & + rhowat*(1.0-osm/ssoil%smass(:,2)), max_ssdn))
         WHERE (soil%isoilm /= 9) ssoil%ssdn(:,2) = MIN(450.0,ssoil%ssdn(:,2))
         ssoil%sdepth(:,2) = ssoil%smass(:,2) / ssoil%ssdn(:,2)
! layer 3        
         sgamm = ssoil%ssdn(:,3) * cgsnow * ssoil%sdepth(:,3)
         osm = ssoil%smass(:,3)
!        ssoil%tggsn(:,3) = ssoil%tggsn(:,3) + canopy%precis *  &
!                        & ( hlf + cswat*(tfrz-ssoil%tggsn(:,3)) ) * osm/ssoil%osnowd  &
!                        & / (sgamm + cswat*canopy%precis*osm/ssoil%osnowd)
         ssoil%tggsn(:,3) = ssoil%tggsn(:,3) + canopy%precis * hlf &
                           * osm / (sgamm * ssoil%osnowd )
         ssoil%smass(:,3) = ssoil%smass(:,3) + canopy%precis &
                           & * osm/ssoil%osnowd
         ssoil%ssdn(:,3) = MAX(120.0,MIN(ssoil%ssdn(:,3) * osm/ssoil%smass(:,3) &
                        & + rhowat*(1.0-osm/ssoil%smass(:,3)), max_ssdn))
         WHERE (soil%isoilm /= 9)  ssoil%ssdn(:,3) = MIN(450.0,ssoil%ssdn(:,3))
         ssoil%sdepth(:,3) = ssoil%smass(:,3)/ssoil%ssdn(:,3)

         canopy%precis = 0.0
      END WHERE
    END WHERE

!jhan:decrease wilting point in UM
   if( cable_runtime%um) then
      swilt_eff = soil%swilt/3
   else
      swilt_eff = soil%swilt
   endif

!    WHERE (ssoil%snowd < 0.1 .AND. canopy%fess > 0.0)
!       canopy%fess = MIN(canopy%fess, &
!            MAX(0.0,(REAL(ssoil%wb(:,1),r_1)-swilt_eff))* soil%zse(1) &
!            * 1000.0 * hl / dels)
!       canopy%fess = MIN(canopy%fess, &
!            REAL((ssoil%wb(:,1)-ssoil%wbice(:,1)),r_1) * soil%zse(1) &
!            * 1000.0 * hl / dels)
!    END WHERE
    ! Calculate snow evaporation total in mm (from W/m2):
   !jhan: Eva uses
    ! EAK 'fess' is for soil evap and 'fes' is for soil evap plus soil puddle evap
    canopy%segg = canopy%fess / hl
    canopy%segg = (canopy%fess + canopy%fes_cor) / hl
    !WHERE (ssoil%snowd > 0.1) canopy%segg = canopy%fess / (hl + hlf) ! EAK aug08
    !WHERE (ssoil%snowd > 0.1) canopy%segg = (canopy%fess + canopy%fes_cor) / (hl + hlf) ! EAK aug08
    !jh:canopy%segg = canopy%fes / hl
    !jh:WHERE (ssoil%snowd > 0.1) canopy%segg = canopy%fes / (hl + hlf) ! EAK aug08

    ! Initialise snow evaporation:
    ssoil%evapsn = 0
    ! Snow evaporation and dew on snow
    !jhan:r2411 uses different scheme here CHECK w Eva
   !jhan: Eva uses fes -> fess
    WHERE (ssoil%snowd > 0.1)
      ssoil%evapsn = dels * (canopy%fess + canopy%fes_cor) / ( hl + hlf )
      xxx = ssoil%evapsn
      WHERE (ssoil%isflag == 0 .and. canopy%fess+ canopy%fes_cor.gt.0.0) &
            & ssoil%evapsn = MIN(ssoil%snowd, xxx ) 
      WHERE ( ssoil%isflag  > 0 .and. canopy%fess+ canopy%fes_cor.gt.0.0) &
            & ssoil%evapsn = MIN(0.9*ssoil%smass(:,1), xxx )
      ssoil%snowd = ssoil%snowd - ssoil%evapsn
      WHERE ( ssoil%isflag > 0 )
        ssoil%smass(:,1) = ssoil%smass(:,1)  - ssoil%evapsn
        ssoil%sdepth(:,1) = MAX(0.02,ssoil%smass(:,1) / ssoil%ssdn(:,1))
      END WHERE
      canopy%segg = 0.0
!      canopy%fess = ssoil%evapsn * (hl + hlf) / dels ! return for hyd. balance
    END WHERE

  END SUBROUTINE snow_accum 


  !-------------------------------------------------------------------------
  SUBROUTINE surfbv (dels, met, ssoil, soil, veg, canopy )
   use cable_common_module
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(met_type), INTENT(INOUT)            :: met    ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE(veg_parameter_type), INTENT(IN)     :: veg
    TYPE(canopy_type), INTENT(IN)         :: canopy

    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
!jhan:default to UM nglacier =0
    INTEGER(i_d), PARAMETER      :: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    INTEGER(i_d)                 :: k
    REAL(r_1), DIMENSION(mp)     :: rnof5
    REAL(r_1), DIMENSION(mp)     :: sfact
    REAL(r_1), DIMENSION(mp)     :: sgamm
    REAL(r_1), DIMENSION(mp)     :: smasstot
    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
    REAL(r_1), DIMENSION(mp)     :: talb ! snow albedo
    REAL(r_1), DIMENSION(mp)     :: tmp ! temporary value
    REAL(r_1), DIMENSION(mp)     :: xxx

    CALL smoisturev ( dels, ssoil, soil, veg)

    DO k = 1, ms
      xxx = REAL(soil%ssat,r_2)
      ssoil%rnof1 = ssoil%rnof1 + REAL((MAX(ssoil%wb(:,k) - xxx, 0.0_r_2) &
                                & * 1000.0),r_1) * soil%zse(k)
      ssoil%wb(:,k) = MAX(0.01,MIN( ssoil%wb(:,k), xxx ))
!      ssoil%wb(:,k) = MIN( ssoil%wb(:,k), soil%ssat )
    END DO
    ! for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
    ! account for soil/ice cracking
    !   fracm = MIN(0.2, 1. - MIN(1., ssoil%wb(:,ms) / soil%sfc ) )
    !   ssoil%wb(:,ms) = ssoil%wb(:,ms) &
    !                  + fracm*ssoil%rnof1/(1000.0*soil%zse(ms))
    !   ssoil%rnof1 = (1. - fracm) * ssoil%rnof1 
    ! the code below is replaced, see subroutine smoistv 
    !   tmp = MAX(MIN(ssoil%wb(:,ms) - soil%sfc, .99 * ssoil%wb(:,ms) &
    !       - ssoil%wbice(:,ms) ) * soil%c3 / 86400., 0.)
    !   ssoil%rnof2 = soil%zse(ms) * 1000. * tmp * dels
    !   ssoil%wb(:,ms) = ssoil%wb(:,ms) - tmp * dels

    !   Scaling  runoff to kg/m^2/s to match rest of the model
    ssoil%sinfil = 0.0
    ! MJT - CCAM has own lake scheme
!    where  ( veg%iveg == 16 ) ! CSIRO - MJT
    where  ( veg%iveg == 17 ) ! IGBP - MJT
      ssoil%sinfil = min( ssoil%rnof1, ssoil%wb_lake + max(0.,canopy%segg))
      ssoil%rnof1 = max( 0.0, ssoil%rnof1 - ssoil%sinfil )
      ssoil%wb_lake = ssoil%wb_lake - ssoil%sinfil
      ssoil%rnof2 = max( 0.0, ssoil%rnof2 - ssoil%wb_lake )
    endwhere

!    ssoil%rnof1 = ssoil%rnof1 / dels
!    ssoil%rnof2 = ssoil%rnof2 / dels
!    ssoil%runoff = ssoil%rnof1 + ssoil%rnof2
!    ssoil%wbtot = 0.0
!    DO k = 1, ms
!      ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k),r_1) * 1000.0 * soil%zse(k)
!    END DO
    !
!jhan: diff min
    !---  glacier formation
    rnof5= 0.

    IF (nglacier == 2) THEN
      smelt1=0.
      WHERE (ssoil%snowd > max_glacier_snowd)
        rnof5 = min(0.1,ssoil%snowd - max_glacier_snowd)
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
      DO k = 1, 3
        WHERE (ssoil%snowd > max_glacier_snowd  .and.  ssoil%isflag > 0)
          sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
          smelt1(:,k) = MIN(rnof5 * ssoil%smass(:,k) / smasstot, &
                          & 0.2 * ssoil%smass(:,k) )
          ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
          ssoil%snowd = ssoil%snowd - smelt1(:,k)
        END WHERE
      END DO
      WHERE (ssoil%isflag > 0 ) rnof5 = smelt1(:,1)+smelt1(:,2)+smelt1(:,3)
    END IF

    ssoil%rnof1 = ssoil%rnof1 / dels + rnof5/dels
    ssoil%rnof2 = ssoil%rnof2 / dels
    ssoil%runoff = ssoil%rnof1 + ssoil%rnof2 

  END SUBROUTINE surfbv

  !-------------------------------------------------------------------------
  ! SUBROUTINE stempv
  !	 calculates temperatures of the soil
  !	 tgg - new soil/snow temperature
  !	 ga - heat flux from the atmosphere (ground heat flux)
  !	 ccnsw - soil thermal conductivity, including water/ice
  !
  SUBROUTINE stempv(dels, canopy, ssoil, soil)
!   use arraydiag_m
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
    INTEGER(i_d)                     :: j,k
    REAL(r_1), DIMENSION(mp,-2:ms)   :: rhs
    REAL(r_1), DIMENSION(mp)         :: sgamm
    REAL(r_2), DIMENSION(mp,ms+3)    :: tmp_mat ! temp. matrix for tggsn & tgg
    REAL(r_2), DIMENSION(mp)         :: xx
    REAL(r_2), DIMENSION(mp)         :: wblfsp 
    real :: snow_ccnsw, exp_arg
    logical :: direct2min = .false.

!   logical stats

    at = 0.0
    bt = 1.0
    ct = 0.0
    coeff = 0.0

    snow_ccnsw = 2.0

   DO k = 1, ms
      do j = 1, mp
         if(soil%isoilm(j) == 9) then
            ccnsw(j,k) = snow_ccnsw
         ELSE
            ew(j) = ssoil%wblf(j,k) * soil%ssat(j)
            
            exp_arg =  (ew(j) * LOG(60.0) ) + (ssoil%wbfice(j,k) * soil%ssat(j) * LOG(250.0) )
            if(exp_arg > 30) direct2min = .true.
            if(direct2min) then
               ccnsw(j,k) = 1.5 * MAX(1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 * soil%ssat(j) &
                      & / MIN(ew(j), 0.5_r_2 * soil%ssat(j) )) ) )
            else         
               ccnsw(j,k) = MIN(soil%cnsd(j) * EXP( exp_arg ), 1.5_r_2) &
                      & * MAX(1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 * soil%ssat(j) &
                      & / MIN(ew(j), 0.5_r_2 * soil%ssat(j) )) ) )
            endif          
            direct2min = .false.
         ENDif 
      END DO
    END DO
    xx = 0. 
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
          ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
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
   
    coeff(:,1-3) = 0.0  ! SO DOES THIS MEAN coeff(:,-2) ??
    ! 3-layer snow points done here
    WHERE (ssoil%isflag /= 0)
      ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1)**2 &
                        & + 0.074, max_sconds) )
      ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2)**2 &
                        & + 0.074, max_sconds) )
      ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3)**2 &
                        & + 0.074, max_sconds) )
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
    DO k = 1, 3
      WHERE (ssoil%isflag /= 0)
        sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
        dtg = dels / sgamm
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

    WHERE (ssoil%isflag /= 0)
      sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
      ! new code
      bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm
      ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + (canopy%ga - ssoil%tggsn(:,1) &
           & * REAL(canopy%dgdtg,r_1) ) * dels / sgamm
      rhs(:,1-3) = ssoil%tggsn(:,1)
    END WHERE
 
    !     note in the following that tgg and tggsn are processed together
    tmp_mat(:,1:3) = REAL(ssoil%tggsn,r_2) 
    tmp_mat(:,4:(ms+3)) = REAL(ssoil%tgg,r_2)

    CALL trimb (at, bt, ct, tmp_mat, ms + 3)
   
    ssoil%tggsn = REAL(tmp_mat(:,:3),r_1)
    ssoil%tgg   = REAL(tmp_mat(:,4:(ms+3)),r_1)
    canopy%sghflux = coefa * (ssoil%tggsn(:,1) - ssoil%tggsn(:,2) )
    canopy%ghflux = coefb * (ssoil%tgg(:,1) - ssoil%tgg(:,2) ) ! +ve downwards

    
  END SUBROUTINE stempv



  !-------------------------------------------------------------------------
 SUBROUTINE snowcheck(dels, ssoil, soil, met )
   use cable_common_module
    REAL(r_1), INTENT(IN)               :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    INTEGER(i_d), PARAMETER :: ntest = 0 !  for prints
    INTEGER(i_d)            :: k,j
      logical :: cable_runtime_coupled = .false.

      do j=1,mp
         if(ssoil%snowd(j) <= 0.0) then
            ssoil%isflag(j) = 0
            ssoil%ssdn(j,:) = 120.0
            ssoil%ssdnn(j) = 120.0
            ssoil%tggsn(j,:) = tfrz
            ssoil%sdepth(j,1) = ssoil%snowd(j) / ssoil%ssdn(j,1)
             
            ssoil%sdepth(j,2) = 0.
            ssoil%sdepth(j,3) = 0.

            ssoil%smass(j,1) = ssoil%snowd(j)
            ssoil%smass(j,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
            ssoil%smass(j,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
      
         ELSEif (ssoil%snowd(j) < snmin * ssoil%ssdnn(j)) then

            if(ssoil%isflag(j) == 1) then
               ssoil%ssdn(j,1) = ssoil%ssdnn(j)
               ssoil%tgg(j,1) = ssoil%tggsn(j,1)
            ENDif 

            ssoil%isflag(j) = 0
            ssoil%ssdnn(j) = min(400.0, MAX(120.0, ssoil%ssdn(j,1)) ) 
     
            ssoil%tggsn(j,:) = min(tfrz,ssoil%tgg(j,1))

            ssoil%sdepth(j,1) = ssoil%snowd(j) / ssoil%ssdn(j,1)
            ssoil%sdepth(j,2) = 0.0     
            ssoil%sdepth(j,3) = 0.0     

            ssoil%smass(j,1) = ssoil%snowd(j)     
            ssoil%smass(j,2) = 0.0     
            ssoil%smass(j,3) = 0.0     

            ssoil%ssdn(j,:) = ssoil%ssdnn(j)
            
            if(.NOT.cable_runtime_coupled) then
               if( soil%isoilm(j) == 9 .and. ktau_gl <= 2 ) &
                  ssoil%ssdnn(j) = 700.0
            endif

         ELSE ! sufficient snow now for 3 layer snowpack

            if (ssoil%isflag(j) == 0) then
               ssoil%tggsn(j,:) = min(tfrz,ssoil%tgg(j,1))

               ssoil%ssdn(j,2) = ssoil%ssdn(j,1)
               ssoil%ssdn(j,3) = ssoil%ssdn(j,1)

               if(.NOT.cable_runtime_coupled) then
                  if( soil%isoilm(j) == 9 .and. ktau_gl <= 2 ) then
                     ssoil%ssdn(j,1)  = 450.0
                     ssoil%ssdn(j,2)  = 580.0
                     ssoil%ssdn(j,3)  = 600.0
                  endif
               endif
               
               ssoil%sdepth(j,1) = ssoil%t_snwlr(j)
               
               ssoil%smass(j,1)  =  ssoil%t_snwlr(j) * ssoil%ssdn(j,1)
               
               ssoil%smass(j,2)  = (ssoil%snowd(j) - ssoil%smass(j,1)) * 0.4
               ssoil%smass(j,3)  = (ssoil%snowd(j) - ssoil%smass(j,1)) * 0.6
               
               ssoil%sdepth(j,2) = ssoil%smass(j,2)/ssoil%ssdn(j,2)
               ssoil%sdepth(j,3) = ssoil%smass(j,3)/ssoil%ssdn(j,3)
               
               ssoil%ssdnn(j) = (ssoil%ssdn(j,1) * ssoil%smass(j,1) + ssoil%ssdn(j,2) &
                    & * ssoil%smass(j,2) + ssoil%ssdn(j,3) * ssoil%smass(j,3) ) &
                    & / ssoil%snowd(j)
            ENDif 
            ssoil%isflag(j) = 1
         END if
      enddo

  END SUBROUTINE snowcheck 


  !******************************************************************
  SUBROUTINE snowl_adjust(dels, ssoil, canopy )
    REAL(r_1), INTENT(IN)               :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE(canopy_type), INTENT(INOUT)    :: canopy
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints
    INTEGER(i_d)             :: k
    REAL(r_2), DIMENSION(mp) :: excd
    REAL(r_2), DIMENSION(mp) :: excm
    REAL(r_2), DIMENSION(mp) :: frac 
    REAL(r_2), DIMENSION(mp) :: xfrac 
    REAL(r_1), DIMENSION(mp) :: osm
! AJA ##############################################
    INTEGER(i_d) :: api ! active patch counter
! ##################################################
    
    ! adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    WHERE (ssoil%isflag > 0)
       WHERE ( ssoil%sdepth(:,1) > ssoil%t_snwlr )
          excd = ssoil%sdepth(:,1) - ssoil%t_snwlr
          excm = excd * ssoil%ssdn(:,1)
          ssoil%sdepth(:,1) = ssoil%sdepth(:,1) - REAL(excd,r_1)
          osm = ssoil%smass(:,1)
          ssoil%smass(:,1) = ssoil%smass(:,1) - REAL(excm,r_1)
             
          osm = ssoil%smass(:,2)
          ssoil%smass(:,2) = MAX(0.01, ssoil%smass(:,2) + REAL(excm,r_1))
          ssoil%ssdn(:,2) = REAL(MAX(120.0_r_2, MIN(REAL(max_ssdn,r_2), ssoil%ssdn(:,2) &
                        & * osm/ssoil%smass(:,2) + ssoil%ssdn(:,1) * excm &
                        & / ssoil%smass(:,2))),r_1)
          ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
          ssoil%tggsn(:,2) = REAL(ssoil%tggsn(:,2) * osm / ssoil%smass(:,2) &
                         & + ssoil%tggsn(:,1) * excm/ ssoil%smass(:,2),r_1)
          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssoil%smass(:,3) = MAX(0.01, &
                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
          ! ssoil%smass(:,3) = ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2)
       ELSEWHERE ! ssoil%sdepth(:,1) < ssoil%t_snwlr
          ! 1st layer
          excd = ssoil%t_snwlr - ssoil%sdepth(:,1)
          excm = excd * ssoil%ssdn(:,2)
          osm = ssoil%smass(:,1)
          ssoil%smass(:,1) = ssoil%smass(:,1) + REAL(excm,r_1)
          ssoil%sdepth(:,1) = ssoil%t_snwlr
          ssoil%ssdn(:,1) = REAL(MAX(120.0_r_2, MIN(REAL(max_ssdn,r_2), ssoil%ssdn(:,1) &
                        & * osm/ssoil%smass(:,1) + ssoil%ssdn(:,2) * excm &
                        & / ssoil%smass(:,1))),r_1)
          ssoil%tggsn(:,1) = REAL(ssoil%tggsn(:,1) * osm / ssoil%smass(:,1) &
                         & + ssoil%tggsn(:,2) * excm/ ssoil%smass(:,1),r_1)
          ! 2nd layer
          ssoil%smass(:,2) = MAX(0.01, ssoil%smass(:,2) - REAL(excm,r_1))
          ssoil%sdepth(:,2) = ssoil%smass(:,2)/ssoil%ssdn(:,2)
          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssoil%smass(:,3) = MAX(0.01, &
                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
          ! ssoil%smass(:,3) = ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2)
        END WHERE
   ! AJA REPLACING THE WHERE LOOP #########################  
      END WHERE 
    ! AJA END WHERE

    DO  api=1,mp
      IF(ssoil%isflag(api).gt.0)THEN
        frac(api) = ssoil%smass(api,2) / MAX(0.02, ssoil%smass(api,3))
        ! if frac > 0.6 or frac < 0.74 do nothing HOW TO translate this to xfrac
        xfrac(api) = 2.0/3.0/ frac(api)
        IF(xfrac(api) > 1.0 )THEN
           excm(api) = (xfrac(api) - 1.0) * ssoil%smass(api,2)
           osm(api) = ssoil%smass(api,2)
           ! changed 0.02 to 0.01 to fix -ve sdepth (EK 21Dec2007)
           ssoil%smass(api,2) = MAX(0.01, ssoil%smass(api,2) + REAL(excm(api),r_1))
           ssoil%tggsn(api,2) = ssoil%tggsn(api,2) * osm(api) / ssoil%smass(api,2) +  &
                         & ssoil%tggsn(api,3) * REAL(excm(api),r_1)/ ssoil%smass(api,2)
           ssoil%ssdn(api,2) = MAX(120.0, MIN(max_ssdn, ssoil%ssdn(api,2)* &
              & osm(api)/ssoil%smass(api,2) + ssoil%ssdn(api,3) &
              * REAL(excm(api),r_1) /ssoil%smass(api,2)) )
           ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
           ssoil%smass(api,3) = MAX(0.01, &
                         & ssoil%snowd(api) - ssoil%smass(api,1) - ssoil%smass(api,2))
           ssoil%sdepth(api,3) = MAX(0.02, ssoil%smass(api,3) / ssoil%ssdn(api,3) )
        ELSE! xfrac < 1
           excm(api) = (1 - xfrac(api)) * ssoil%smass(api,2)
           ssoil%smass(api,2) = MAX(0.01, ssoil%smass(api,2) - REAL(excm(api),r_1))
           ssoil%sdepth(api,2) = MAX(0.02, ssoil%smass(api,2) / ssoil%ssdn(api,2) )

           osm(api) = ssoil%smass(api,3)
           ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
           ssoil%smass(api,3) = MAX(0.01, &
                         & ssoil%snowd(api) - ssoil%smass(api,1) - ssoil%smass(api,2))
!          if (ssoil%smass(api,3).eq.0) write(6,*) 'rmlcheck ',api,ssoil%snowd(api), &
! ssoil%smass(api,1),ssoil%smass(api,2),xfrac(api),osm(api),ssoil%sdepth(api,2), &
! ssoil%ssdn(api,2),excm(api)
           ssoil%tggsn(api,3) = ssoil%tggsn(api,3) * osm(api) / ssoil%smass(api,3) +  &
                         & ssoil%tggsn(api,2) * REAL(excm(api),r_1)/ ssoil%smass(api,3)
           ssoil%ssdn(api,3) = MAX(120.0, MIN(max_ssdn, ssoil%ssdn(api,3)* &
                osm(api)/ssoil%smass(api,3) + ssoil%ssdn(api,2) &
                * REAL(excm(api),r_1) / ssoil%smass(api,3)) )
           ssoil%sdepth(api,3) = ssoil%smass(api,3) /  ssoil%ssdn(api,3)
        END IF
        ssoil%isflag(api) = 1
        ssoil%ssdnn(api) = (ssoil%ssdn(api,1) * ssoil%sdepth(api,1) + ssoil%ssdn(api,2) &
             & * ssoil%sdepth(api,2) + ssoil%ssdn(api,3) * ssoil%sdepth(api,3) ) &
             & / (ssoil%sdepth(api,1) + ssoil%sdepth(api,2) + ssoil%sdepth(api,3))
      END IF
   END DO

  END SUBROUTINE snowl_adjust


  !******************************************************************
  SUBROUTINE soilfreeze(dels, soil, ssoil)
   use cable_common_module
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    REAL(r_2), DIMENSION(mp)           :: sicefreeze
    REAL(r_2), DIMENSION(mp)           :: sicemelt
    REAL(r_1), DIMENSION(mp)           :: xx
    INTEGER(i_d) k

!jhan:changing to reals() changed answer    
    xx = 0.
    DO k = 1, ms
       WHERE (ssoil%tgg(:,k) < tfrz &
            & .AND. frozen_limit * ssoil%wb(:,k) - ssoil%wbice(:,k) > .001)
          sicefreeze = MIN( MAX(0.0_r_2,(frozen_limit*ssoil%wb(:,k)-ssoil%wbice(:,k))) &
               & * soil%zse(k) * 1000.0, &
               & (tfrz - ssoil%tgg(:,k) ) * ssoil%gammzz(:,k) / hlf )
          ssoil%wbice(:,k) = MIN( ssoil%wbice(:,k) + sicefreeze / (soil%zse(k) &
               * 1000.0), frozen_limit * ssoil%wb(:,k) )
          xx=soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( &
               REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2) &
               + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * REAL(cswat * rhowat,r_2) &
               + ssoil%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2), REAL(xx,r_2)) &
               * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) + REAL(sicefreeze,r_1) &
               * hlf / REAL(ssoil%gammzz(:,k),r_1)
       ELSEWHERE (ssoil%tgg(:,k) > tfrz .AND. ssoil%wbice(:,k) > 0.)
          sicemelt = MIN(ssoil%wbice(:,k) * soil%zse(k) * 1000.0, &
               & (ssoil%tgg(:,k) - tfrz) * ssoil%gammzz(:,k) / hlf)
          ssoil%wbice(:,k) = MAX(0.0_r_2, ssoil%wbice(:,k) - sicemelt &
               / (soil%zse(k) * 1000.0) )
          xx = soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( &
               REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2) &
               + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * REAL(cswat*rhowat,r_2) &
               + ssoil%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2), &
               REAL(xx,r_2) ) * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) - REAL(sicemelt,r_1) &
               * hlf / REAL(ssoil%gammzz(:,k),r_1)
       END WHERE
    END DO
  END SUBROUTINE soilfreeze


  !******************************************************************
  SUBROUTINE remove_trans(dels, soil, ssoil, canopy, veg)
    ! Removes transpiration water from soil.
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
!    REAL(r_2), DIMENSION(mp,ms)   :: evapfb
    REAL(r_2), DIMENSION(mp,0:ms) :: diff 
    REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
    INTEGER(i_d) k

!    diff(:,:) = 0.0  
!    evapfb(:,:) = 0.0

!    DO k = 1,ms
!      WHERE (canopy%fevc > 0.0)
!        ssoil%wb(:,k) = ssoil%wb(:,k)-canopy%evapfbl(:,k)/(soil%zse(k)*1000.0)
!      END WHERE
!    END DO
!    diff(:,:) = 0.0  
!    DO k = 1,ms
!       ! Removing transpiration from soil:
!       WHERE (canopy%fevc > 0.0)     ! convert to mm/dels
!          ! Calculate the amount (perhaps moisture/ice limited)
!          ! which can be removed:
!          xx = canopy%fevc * dels / hl * veg%froot(:,k) + diff(:,k-1)
!          !jhan:Eva's r2411
!          xx = ssoil%evapfbl(:,k) + diff(:,k-1)   ! kg/m2
!          diff(:,k) = ( MAX( 0.0, ssoil%wb(:,k) - soil%swilt) &      ! m3/m3
!                  & * soil%zse(k) * 1000.0 - xx) / (soil%zse(k) * 1000.0)
!
!          !diff(:,k) = (MAX( 0.0, MIN( ssoil%wb(:,k) - soil%swilt, &      ! m3/m3
!          !   & ssoil%wb(:,k) - ssoil%wbice(:,k) ) ) &
!          !   & * soil%zse(k) * 1000.0 - xx) / (soil%zse(k) * 1000.0)
!
!!          diff(:,k) = (MAX( 0.0_r_2, MIN( ssoil%wb(:,k) - &
!!               REAL(soil%swilt,r_2), ssoil%wb(:,k) - ssoil%wbice(:,k) ) ) &
!!               * REAL(soil%zse(k),r_2) * 1000.0_r_2 - xx)  &
!!               / REAL(soil%zse(k) * 1000.0,r_2)
!          WHERE ( diff(:,k) > 0.0 )
!            ssoil%wb(:,k) = ssoil%wb(:,k) - xx / (soil%zse(k) * 1000.0)
!            diff(:,k) = 0.0
!            evapfb(:,k) = xx / (soil%zse(k) * 1000.0)
!          ELSEWHERE
!            diff(:,k) = xx
!          ENDWHERE
!!!          !    evapfbl(:,k) = (MIN(canopy%fevc * dels / hl * veg%froot(:,k), &
!!!          !         MAX(0._r_2, MIN(ssoil%wb(:,k) &
!!!          !         - soil%swilt,ssoil%wb(:,k)-ssoil%wbice(:,k))) &
!!!          !         * soil%zse(k) * 1000.)) / (soil%zse(k) * 1000.)
!!!          !      ! Remove this amount from  the soil:
!!!          !      ssoil%wb(:,k) = ssoil%wb(:,k) -   evapfbl(:,k)
!        END WHERE
!     END DO
     xx = 0.; xxd = 0.; diff(:,:) = 0.
     DO k = 1,ms
     ! Removing transpiration from soil:
       WHERE (canopy%fevc > 0.0 )     ! convert to mm/dels
         ! Calculate the amount (perhaps moisture/ice limited)
         ! which can be removed:
         xx = canopy%fevc * dels / hl * veg%froot(:,k) + diff(:,k-1)   ! kg/m2
         !xx =  ssoil%evapfbl(:,k) +  diff(:,k-1)   ! kg/m2
         !diff(:,k) = (MAX( 0.0, ssoil%wb(:,k) - soil%swilt) &      ! m3/m3
         !          & * soil%zse(k)*1000.0 - xx) / (soil%zse(k) * 1000.0)
         diff(:,k) = MAX( 0.0, ssoil%wb(:,k) - soil%swilt) &      ! m3/m3
                   & * soil%zse(k)*1000.0
        !        diff(:,k) = (MAX( 0.0, MIN( ssoil%wb(:,k) - soil%swilt, &      ! m3/m3
        !                  & ssoil%wb(:,k) - ssoil%wbice(:,k) ) ) &
        !                  & * soil%zse(k) * 1000.0 - xx) / (soil%zse(k) * 1000.0)
         xxd = xx - diff(:,k)
         !WHERE ( diff(:,k) .gt. 0.0 )
         WHERE ( xxd .gt. 0.0 )
           ssoil%wb(:,k) = ssoil%wb(:,k) - diff(:,k) / (soil%zse(k)*1000.0)
           diff(:,k) = xxd
!          evapfbl(:,k) = xx
         ELSEWHERE
           ssoil%wb(:,k) = ssoil%wb(:,k) - xx / (soil%zse(k)*1000.0)
           diff(:,k) = 0.0
         ENDWHERE
       END WHERE
     END DO
!
!!    PRINT *, 'Before adjusting fevc'
!     ! Adjust fevc
!     WHERE (canopy%fevc > 0.) ! convert to mm/dels
!        canopy%fevc = (evapfbl(:,1)*soil%zse(1)+evapfbl(:,2)*soil%zse(2) &
!            & +evapfbl(:,3)*soil%zse(3)+evapfbl(:,4)*soil%zse(4)+evapfbl(:,5) &
!            & *soil%zse(5)+evapfbl(:,6)*soil%zse(6))*1000.*hl/dels
!     END WHERE

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
  !	 fes   - soil evaporation (W/m2)
  !	 isoil - soil type
  !	 ivegt - vegetation type
  ! Output
  !	 ssoil
  SUBROUTINE soil_snow(dels, soil, ssoil, canopy, met, bal, veg)
   use cable_common_module
!  use arraydiag_m
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    TYPE (balances_type), INTENT(INOUT)      :: bal
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints
    INTEGER(i_d)             :: k
    REAL(r_1), DIMENSION(mp) :: snowmlt
    REAL(r_1), DIMENSION(mp) :: totwet
    REAL(r_1), DIMENSION(mp) :: weting
    REAL(r_1), DIMENSION(mp) :: xxx, tgg_old, tggsn_old
    REAL(r_2), DIMENSION(mp) :: xx,deltat,sinfil1,sinfil2,sinfil3 
    REAL(r_1)                :: zsetot
    logical :: cable_runtime_coupled = .false.
    integer, save :: ktau =0 
    INTEGER(i_d) :: idjd1,idjd2,idjd3,idjd4,idjd5

    ! Find initial moisture total for checking purpose
    idjd1 = 2356 
    idjd2 = 340
    idjd3 = 14470 ! Sahara
    !idjd3 = 381 Cabauw
    idjd4 = 351  ! Amazon
    idjd5 = 351  ! Amazon
      ktau = ktau +1 
!   logical stats
!   integer m

   if( cable_runtime%um) then
      max_glacier_snowd = 50000.0
   else
      max_glacier_snowd = 1100.0
   endif


!jhan: peculiar to UM        
    zsetot = sum(soil%zse) 
    ssoil%tggav = 0.
    DO k = 1, ms
     ssoil%tggav = ssoil%tggav  + soil%zse(k)*ssoil%tgg(:,k)/zsetot
    END DO

!    IF(ntest>0) THEN
    !    IF(mp > 10300) THEN
!      print *,'idjd1',idjd1,mp
!      PRINT *,'f', ktau,mp, ssoil%isflag(idjd1),ssoil%osnowd(idjd1),ssoil%snowd(idjd1)
!      PRINT *,'f', ktau,mp, ssoil%isflag(idjd1),ssoil%osnowd(idjd1),ssoil%snowd(idjd1)
! PRINT *,'s_1', met%tk(idjd1),met%tc(idjd1)
! PRINT *,'s_11', met%tk(idjd1),met%tc(idjd1),met%qv(idjd1)
! PRINT *,'s_12', met%tk(idjd1),met%tc(idjd1),met%qv(idjd1),met%ua(idjd1)
! PRINT *,'s_2', met%tk(idjd1),met%tc(idjd1),met%qv(idjd1),met%ua(idjd1),met%fsd(idjd1,:)
! PRINT *,'3',  met%fsd(idjd1,2),met%fld(idjd1), met%precip(idjd1),canopy%precis(idjd1)

!      PRINT 10, ktau,mp, ssoil%isflag(idjd1),ssoil%osnowd(idjd1),ssoil%snowd(idjd1), &
!                met%tk(idjd1),met%qv(idjd1),met%ua(idjd1),met%fsd(idjd1,1), &
!               met%fsd(idjd1,2),met%fld(idjd1), met%precip(idjd1),canopy%precis(idjd1)
! 10    format(x,'soilsnowv befstempv,ktau=',2i6,2x,i2,2f8.2,f6.1,f7.4,f5.1,x,3f6.0,2f6.3)
!       PRINT 11, ktau,canopy%dgdtg(idjd1),canopy%fevc(idjd1),canopy%fevw(idjd1), &
!       canopy%fess(idjd1),canopy%fe(idjd1),canopy%fhs(idjd1),canopy%fhv(idjd1), &
!       canopy%fh(idjd1),canopy%fnv(idjd1),canopy%fns(idjd1),canopy%precis(idjd1)
! 11    format(x,'bsoilsnow',i6,10f7.2,2x,f6.3)
!       PRINT 12, ktau,canopy%ga(idjd1),(ssoil%ssdn(idjd1,k),k=1,3), &
!                (ssoil%sdepth(idjd1,k),k=1,3),(ssoil%smass(idjd1,k),k=1,3)
! 12    format(x,'bga,ssdn,sdepth smass',i6,f7.1,3f6.1,2x,3f8.3,2x,3f8.1)
!       PRINT 13, ktau,(ssoil%tgg(idjd1,k),k=1,ms),(ssoil%tggsn(idjd1,k),k=1,3)
! 13    format(x,'btgg tggsn ',i6,6f7.2,2x,3f7.2)
!       PRINT 14, ktau,(ssoil%wb(idjd1,k),k=1,ms),(ssoil%wbice(idjd1,k),k=1,ms),veg%froot(idjd1,:)
! 14    format(x,'bwb wbice',i6,6f6.3,2x,6f6.3,6f6.3)
!    ssoil%wbtot = 0.0
!    DO k = 1, ms
!       ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k)*1000.0*soil%zse(k),r_1)
!    END DO
!      PRINT 16, idjd1, ktau_gl,met%tk(idjd1),met%qv(idjd1),met%ua(idjd1),met%fsd(idjd1,1)+ &
!               met%fsd(idjd1,2),met%fld(idjd1), met%precip(idjd1),canopy%precis(idjd1), &
!               canopy%ga(idjd1),canopy%dgdtg(idjd1),canopy%fevc(idjd1),canopy%fevw(idjd1), &
!               canopy%fess(idjd1),canopy%fe(idjd1),canopy%fhs(idjd1),canopy%fhv(idjd1), &
!               canopy%fh(idjd1),canopy%fnv(idjd1),canopy%fns(idjd1),(ssoil%tgg(idjd1,k),k=1,ms), &
!               (ssoil%wb(idjd1,k),k=1,ms),ssoil%wbtot(idjd1), &
!               canopy%precis(idjd1)-(canopy%fevc(idjd1)+canopy%fess(idjd1))*dels/hl, &
!               ssoil%rnof1(idjd1)+ssoil%rnof2(idjd1),ssoil%rtsoil(idjd1),soil%swilt(idjd1),soil%sfc(idjd1)
!16    format('idjd1bef',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,f6.0,2f5.2)
!      PRINT 161, idjd2, ktau_gl,met%tk(idjd2),met%qv(idjd2),met%ua(idjd2),met%fsd(idjd2,1)+ &
!               met%fsd(idjd2,2),met%fld(idjd2), met%precip(idjd2),canopy%precis(idjd2), &
!               canopy%ga(idjd2),canopy%dgdtg(idjd2),canopy%fevc(idjd2),canopy%fevw(idjd2), &
!               canopy%fess(idjd2),canopy%fe(idjd2),canopy%fhs(idjd2),canopy%fhv(idjd2), &
!               canopy%fh(idjd2),canopy%fnv(idjd2),canopy%fns(idjd2),(ssoil%tgg(idjd2,k),k=1,ms), &
!               (ssoil%wb(idjd2,k),k=1,ms),ssoil%wbtot(idjd2), &
!               canopy%precis(idjd2)-(canopy%fevc(idjd2)+canopy%fess(idjd2))*dels/hl, &
!               ssoil%rnof1(idjd2)+ssoil%rnof2(idjd2),ssoil%rtsoil(idjd2),soil%swilt(idjd2),soil%sfc(idjd2)
!161    format('idjd2bef',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,f6.0,2f5.2)
!      PRINT 162, idjd3, ktau_gl,met%tk(idjd3),met%qv(idjd3),met%ua(idjd3),met%fsd(idjd3,1)+ &
!               met%fsd(idjd3,2),met%fld(idjd3), met%precip(idjd3),canopy%precis(idjd3), &
!               canopy%ga(idjd3),canopy%dgdtg(idjd3),canopy%fevc(idjd3),canopy%fevw(idjd3), &
!               canopy%fess(idjd3),canopy%fe(idjd3),canopy%fhs(idjd3),canopy%fhv(idjd3), &
!               canopy%fh(idjd3),canopy%fnv(idjd3),canopy%fns(idjd3),(ssoil%tgg(idjd3,k),k=1,ms), &
!               (ssoil%wb(idjd3,k),k=1,ms),ssoil%wbtot(idjd3), &
!               canopy%precis(idjd3)-(canopy%fevc(idjd3)+canopy%fess(idjd3))*dels/hl, &
!               ssoil%rnof1(idjd3)+ssoil%rnof2(idjd3),ssoil%rtsoil(idjd3),soil%swilt(idjd3),soil%sfc(idjd3)
!162    format('idjd3bef',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,f6.0,2f5.2)
!     ENDIF


!jhan: UM uses spatially explicit var, parametrized in Mk3L          
   if( cable_runtime%offline .or. cable_runtime%mk3l ) then
        ssoil%t_snwlr = 0.05
   endif

    ssoil%fwtop1 = 0.0
    ssoil%fwtop2 = 0.0
    ssoil%fwtop3 = 0.0
    ssoil%runoff = 0.0 ! initialise total runoff
    ssoil%rnof1 = 0.0 ! initialise surface runoff
    ssoil%rnof2 = 0.0 ! initialise deep drainage
    ssoil%smelt = 0.0 ! initialise snowmelt
    ssoil%dtmlt = 0.0 
    ssoil%osnowd = ssoil%snowd


if(.NOT.cable_runtime_coupled) then
!jhan: IF block peculiar to UM        
    IF (ktau_gl <= 1) THEN
      canopy%dgdtg = 0.0
      ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
      ssoil%wbtot = 0.0
      DO k = 1, ms
        !jh:ssoil%wbtot = ssoil%wbtot + ssoil%wb(:,k) * 1000.0 * soil%zse(k)
         !jhan:Eva uses
        ssoil%wb(:,k)  = min( soil%ssat,max ( ssoil%wb(:,k), soil%swilt ))
      END DO

      !jhan:Eva uses
      ssoil%wb(:,ms-2)  = min( soil%ssat,max ( ssoil%wb(:,ms-2), 0.5*(soil%sfc+soil%swilt) ))
      ssoil%wb(:,ms-1)  = min( soil%ssat,max ( ssoil%wb(:,ms-1), 0.8*soil%sfc ))
      ssoil%wb(:,ms)    = min( soil%ssat,max ( ssoil%wb(:,ms),   soil%sfc) )
      
      DO k = 1, ms
        WHERE (ssoil%tgg(:,k) <= tfrz .and. ssoil%wbice(:,k) <= 0.01)
         !jhan:Eva uses 0.1 -> 0.5
          ssoil%wbice(:,k) = 0.5 * ssoil%wb(:,k)
        END WHERE
        WHERE (ssoil%tgg(:,k) < tfrz)
          ssoil%wbice(:,k) = frozen_limit * ssoil%wb(:,k)
        END WHERE
      END DO

      !jhan:this whereblock added in r2411
      WHERE (soil%isoilm == 9) 
        ssoil%snowd = max_glacier_snowd
        ssoil%osnowd = max_glacier_snowd
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
      
!! RML DON'T NEED THESE 2 LINES NOW ADDED BELOW IF BLOCK ??
!jhan:have to check where%gammzz is initialized and whether it needs it ktau/ktau_gl etc

      xx=soil%css * soil%rhosoil
      ssoil%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
           & + (ssoil%wb(:,1) - ssoil%wbice(:,1) ) * cswat * rhowat &
           & + ssoil%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)
    END IF
endif  ! if(.NOT.cable_runtime_coupled)

    xx=soil%css * soil%rhosoil
    IF (ktau <= 1) ssoil%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
             & + (ssoil%wb(:,1) - ssoil%wbice(:,1) ) * cswat * rhowat &
             & + ssoil%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1) + &
             & (1. - ssoil%isflag) * cgsnow * ssoil%snowd



    DO k = 1, ms ! for stempv
       ! Set liquid soil water fraction (fraction of saturation value):
       ssoil%wblf(:,k) = MAX( 0.01_r_2, (ssoil%wb(:,k) - ssoil%wbice(:,k)) ) &
            & / REAL(soil%ssat,r_2)
       ! Set ice soil water fraction (fraction of saturation value):
       ssoil%wbfice(:,k) = REAL(ssoil%wbice(:,k),r_1) / soil%ssat
    END DO
  
    CALL snowcheck (dels, ssoil, soil, met )

    CALL snowdensity (dels, ssoil, soil)

    CALL snow_accum (dels, canopy, met, ssoil, soil )

    CALL snow_melting (dels, snowmlt, ssoil, soil )

    ! Add snow melt to global snow melt variable:
    ssoil%smelt = snowmlt

    ! Adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    CALL snowl_adjust(dels, ssoil, canopy )

    CALL stempv(dels, canopy, ssoil, soil)
    !jhan: Eva uses
    ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)

    CALL snow_melting (dels, snowmlt, ssoil, soil )
    
    ! Add new snow melt to global snow melt variable: 
    ssoil%smelt = ssoil%smelt + snowmlt

    !jhan: Eva r2411 comments out
    !where(  ssoil%dtmlt(:,1) > 0.0001) canopy%fhs = canopy%fhs+ssoil%dtmlt(:,1)*ssoil%dfh_dtg
    !where( ssoil%dtmlt(:,1) > 0.0001) canopy%fes = canopy%fes+ssoil%dtmlt(:,1)* &
    !                                    (ssoil%cls*ssoil%dfe_ddq * ssoil%ddq_dtg)

!jhan: looks like UM does this diff, not investigated yet
!jhan:#ifdef ONLINE_UM 
    CALL remove_trans(dels, soil, ssoil, canopy, veg)

!    CALL  soilfreeze(dels, soil, ssoil)
!
!
!   !jhan:Eva uses this block
!    totwet = canopy%precis + ssoil%smelt
!    weting = totwet + max(0.,ssoil%pudsto - canopy%fesp/hl*dels) ! total available liquid including puddle
!    xxx=soil%ssat - ssoil%wb(:,1)
!   
!      !jhan:Eva's r2411 uses 0.95
!    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*rhowat, weting) !soil capacity
!    xxx=soil%ssat - ssoil%wb(:,2)
!    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*rhowat, weting - sinfil1) !soil capacity
!    xxx=soil%ssat - ssoil%wb(:,3)
!    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*rhowat,weting-sinfil1-sinfil2)
!    ssoil%fwtop1 = sinfil1 / dels - canopy%segg          ! net water flux to the soil
!    ssoil%fwtop2 = sinfil2 / dels           ! net water flux to the soil
!    ssoil%fwtop3 = sinfil3 / dels           ! net water flux to the soil
!!   Puddle for the next time step
!    ssoil%pudsto = max( 0., weting - sinfil1 - sinfil2 - sinfil3 )
!    ssoil%rnof1 = max(0.,ssoil%pudsto - ssoil%pudsmx)
!    ssoil%pudsto = ssoil%pudsto - ssoil%rnof1
!
!
!    CALL surfbv(dels, met, ssoil, soil, veg, canopy )

!    ssoil%wbtot = 0.0
!    DO k = 1, ms
!       ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k)*1000.0*soil%zse(k),r_1)
!    END DO
!      PRINT 17, idjd1, ktau,met%tk(idjd1),met%qv(idjd1),met%ua(idjd1),met%fsd(idjd1,1)+ &
!               met%fsd(idjd1,2),met%fld(idjd1), met%precip(idjd1),canopy%precis(idjd1), &
!               canopy%ga(idjd1),canopy%dgdtg(idjd1),canopy%fevc(idjd1),canopy%fevw(idjd1), &
!               canopy%fess(idjd1),canopy%fe(idjd1),canopy%fhs(idjd1),canopy%fhv(idjd1), &
!               canopy%fh(idjd1),canopy%fnv(idjd1),canopy%fns(idjd1),(ssoil%tgg(idjd1,k),k=1,ms), &
!               (ssoil%wb(idjd1,k),k=1,ms),ssoil%wbtot(idjd1), &
!               canopy%precis(idjd1)-(canopy%fevc(idjd1))*dels/hl, &
!!               canopy%precis(idjd1)-(canopy%fevc(idjd1)+canopy%fess(idjd1))*dels/hl, &
!               ssoil%rnof1(idjd1)+ssoil%rnof2(idjd1),cable_runtime_coupled
!17    format('idjd1af1',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,a3)
!      PRINT 171, idjd2, ktau,met%tk(idjd2),met%qv(idjd2),met%ua(idjd2),met%fsd(idjd2,1)+ &
!               met%fsd(idjd2,2),met%fld(idjd2), met%precip(idjd2),canopy%precis(idjd2), &
!               canopy%ga(idjd2),canopy%dgdtg(idjd2),canopy%fevc(idjd2),canopy%fevw(idjd2), &
!               canopy%fess(idjd2),canopy%fe(idjd2),canopy%fhs(idjd2),canopy%fhv(idjd2), &
!               canopy%fh(idjd2),canopy%fnv(idjd2),canopy%fns(idjd2),(ssoil%tgg(idjd2,k),k=1,ms), &
!               (ssoil%wb(idjd2,k),k=1,ms),ssoil%wbtot(idjd2), &
!               canopy%precis(idjd2)-(canopy%fevc(idjd2))*dels/hl, &
!               ssoil%rnof1(idjd2)+ssoil%rnof2(idjd2),cable_runtime_coupled
!171    format('idjd2af1',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,a3)
!      PRINT 172, idjd3, ktau,met%tk(idjd3),met%qv(idjd3),met%ua(idjd3),met%fsd(idjd3,1)+ &
!               met%fsd(idjd3,2),met%fld(idjd3), met%precip(idjd3),canopy%precis(idjd3), &
!               canopy%ga(idjd3),canopy%dgdtg(idjd3),canopy%fevc(idjd3),canopy%fevw(idjd3), &
!               canopy%fess(idjd3),canopy%fe(idjd3),canopy%fhs(idjd3),canopy%fhv(idjd3), &
!               canopy%fh(idjd3),canopy%fnv(idjd3),canopy%fns(idjd3),(ssoil%tgg(idjd3,k),k=1,ms), &
!               (ssoil%wb(idjd3,k),k=1,ms),ssoil%wbtot(idjd3), &
!               canopy%precis(idjd3)-(canopy%fevc(idjd3))*dels/hl, &
!               ssoil%rnof1(idjd3)+ssoil%rnof2(idjd3),cable_runtime_coupled
!172    format('idjd3af1',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,a3)

    CALL  soilfreeze(dels, soil, ssoil)


   !jhan:Eva uses this block
    totwet = canopy%precis + ssoil%smelt 
    weting = totwet + max(0.,ssoil%pudsto - canopy%fesp/hl*dels) ! total available liquid including puddle
    xxx=soil%ssat - ssoil%wb(:,1)
   
      !jhan:Eva's r2411 uses 0.95
    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*rhowat, weting) !soil capacity
    xxx=soil%ssat - ssoil%wb(:,2)
    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*rhowat, weting - sinfil1) !soil capacity
    xxx=soil%ssat - ssoil%wb(:,3)
    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*rhowat,weting-sinfil1-sinfil2)
    ssoil%fwtop1 = sinfil1 / dels - canopy%segg          ! net water flux to the soil
    ssoil%fwtop2 = sinfil2 / dels           ! net water flux to the soil
    ssoil%fwtop3 = sinfil3 / dels           ! net water flux to the soil
!   Puddle for the next time step
    ssoil%pudsto = max( 0., weting - sinfil1 - sinfil2 - sinfil3 )
    ssoil%rnof1 = max(0.,ssoil%pudsto - ssoil%pudsmx)
    ssoil%pudsto = ssoil%pudsto - ssoil%rnof1

    CALL surfbv(dels, met, ssoil, soil, veg, canopy )

!___
    canopy%fhs_cor = ssoil%dtmlt(:,1)*ssoil%dfh_dtg
    canopy%fes_cor = ssoil%dtmlt(:,1)*(ssoil%cls*ssoil%dfe_ddq * ssoil%ddq_dtg)

    canopy%fhs = canopy%fhs+canopy%fhs_cor
    canopy%fes = canopy%fes+canopy%fes_cor

    CALL hydraulic_redistribution(dels,soil,ssoil,canopy,veg, met)

    ssoil%smelt = ssoil%smelt/dels

    ! Set weighted soil/snow surface temperature
    ssoil%tss=(1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)

    ssoil%wbtot = 0.0
    DO k = 1, ms
       ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k)*1000.0*soil%zse(k),r_2)
    END DO

!    IF(ntest>0) THEN
    !    IF(mp > 10300) THEN
!      print *,'idjd1after',idjd1
!      PRINT 10, ktau,mp, ssoil%isflag(idjd1),ssoil%osnowd(idjd1),ssoil%snowd(idjd1), &
!                met%tk(idjd1),met%qv(idjd1),met%ua(idjd1), met%fsd(idjd1,1), &
!               met%fsd(idjd1,2),met%fld(idjd1), met%precip(idjd1),canopy%precis(idjd1)
!! 10    format(x,'soilsnowv befstempv,ktau=',2i6,2x,i2,2f8.2,f5.0,f5.1,f7.4,f5.1,x,2f6.0,2f6.3)
!       PRINT 11, ktau,canopy%dgdtg(idjd1),canopy%fevc(idjd1),canopy%fevw(idjd1), &
!       canopy%fess(idjd1),canopy%fe(idjd1),canopy%fhs(idjd1),canopy%fhv(idjd1), &
!       canopy%fh(idjd1),canopy%fnv(idjd1),canopy%fns(idjd1),canopy%precis(idjd1)
!! 11    format(x,'bsoilsnow',i6,10f7.2,2x,f6.3)
!       PRINT 12, ktau,canopy%ga(idjd1),(ssoil%ssdn(idjd1,k),k=1,3), &
!                (ssoil%sdepth(idjd1,k),k=1,3),(ssoil%smass(idjd1,k),k=1,3)
!! 12    format(x,'bga,ssdn,sdepth smass',i6,f7.1,3f6.1,2x,3f8.3,2x,3f8.1)
!       PRINT 13, ktau,(ssoil%tgg(idjd1,k),k=1,ms),(ssoil%tggsn(idjd1,k),k=1,3)
!! 13    format(x,'btgg tggsn ',i6,6f7.2,2x,3f7.2)
!       PRINT 14, ktau,(ssoil%wb(idjd1,k),k=1,ms),(ssoil%wbice(idjd1,k),k=1,ms),veg%froot(idjd1,:)
!! 14    format(x,'bwb wbice',i6,6f6.3,2x,6f6.3,6f6.3)
!      PRINT 18, idjd1, ktau,met%tk(idjd1),met%qv(idjd1),met%ua(idjd1),met%fsd(idjd1,1)+ &
!               met%fsd(idjd1,2),met%fld(idjd1), met%precip(idjd1),canopy%precis(idjd1), &
!               canopy%ga(idjd1),canopy%dgdtg(idjd1),canopy%fevc(idjd1),canopy%fevw(idjd1), &
!               canopy%fess(idjd1),canopy%fe(idjd1),canopy%fhs(idjd1),canopy%fhv(idjd1), &
!               canopy%fh(idjd1),canopy%fnv(idjd1),canopy%fns(idjd1),(ssoil%tgg(idjd1,k),k=1,ms), &
!               (ssoil%wb(idjd1,k),k=1,ms),ssoil%wbtot(idjd1), &
!               canopy%precis(idjd1)-(canopy%fevc(idjd1)+canopy%fess(idjd1))*dels/hl, &
!               ssoil%rnof1(idjd1)+ssoil%rnof2(idjd1),canopy%fesp(idjd1)
!18    format('idjd1af2',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,f6.1)
!      PRINT 181, idjd2, ktau,met%tk(idjd2),met%qv(idjd2),met%ua(idjd2),met%fsd(idjd2,1)+ &
!               met%fsd(idjd2,2),met%fld(idjd2), met%precip(idjd2),canopy%precis(idjd2), &
!               canopy%ga(idjd2),canopy%dgdtg(idjd2),canopy%fevc(idjd2),canopy%fevw(idjd2), &
!               canopy%fess(idjd2),canopy%fe(idjd2),canopy%fhs(idjd2),canopy%fhv(idjd2), &
!               canopy%fh(idjd2),canopy%fnv(idjd2),canopy%fns(idjd2),(ssoil%tgg(idjd2,k),k=1,ms), &
!               (ssoil%wb(idjd2,k),k=1,ms),ssoil%wbtot(idjd2), &
!               canopy%precis(idjd2)-(canopy%fevc(idjd2)+canopy%fess(idjd2))*dels/hl, &
!               ssoil%rnof1(idjd2)+ssoil%rnof2(idjd2),canopy%fesp(idjd2)
!181    format('idjd2af2',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,f6.1)
!      PRINT 182, idjd3, ktau,met%tk(idjd3),met%qv(idjd3),met%ua(idjd3),met%fsd(idjd3,1)+ &
!               met%fsd(idjd3,2),met%fld(idjd3), met%precip(idjd3),canopy%precis(idjd3), &
!               canopy%ga(idjd3),canopy%dgdtg(idjd3),canopy%fevc(idjd3),canopy%fevw(idjd3), &
!               canopy%fess(idjd3),canopy%fe(idjd3),canopy%fhs(idjd3),canopy%fhv(idjd3), &
!               canopy%fh(idjd3),canopy%fnv(idjd3),canopy%fns(idjd3),(ssoil%tgg(idjd3,k),k=1,ms), &
!               (ssoil%wb(idjd3,k),k=1,ms),ssoil%wbtot(idjd3), &
!               canopy%precis(idjd3)-(canopy%fevc(idjd3)+canopy%fess(idjd3))*dels/hl, &
!               ssoil%rnof1(idjd3)+ssoil%rnof2(idjd3),canopy%fesp(idjd3)
!182    format('idjd3af2',i5,i3,f6.1,f7.4,f5.1,2f5.0,x,2f7.4,f6.1,f5.0,x,4f6.1,x,3f6.1,x,2f5.0, &
!              6f6.1,6(f5.4),1x,f9.4,2f7.4,f6.1)
!    ENDIF

  END SUBROUTINE soil_snow

  !+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++
  ! Sciences from Ryel et al. Oecologia, 2002; Lee et al., 2005, PNAS
  ! Code by LiLH 16 Feb, 2011
  ! Fixed problem of negative wb in global run by BP Mar 2011
  SUBROUTINE hydraulic_redistribution(dels, soil, ssoil, canopy, veg, met)
    REAL(r_1),                 INTENT(IN) :: dels ! integration time step (s)
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
    INTEGER(i_d) :: idjd1,idjd2,idjd3,idjd4,idjd5

    ! Find initial moisture total for checking purpose
    idjd1 = 738  
    idjd2 = 13940
    idjd3 = 14470  ! Sahara
    !idjd3 = 351  ! Cabau
    idjd4 = 351  ! Amazon
    idjd5 = 351  ! Amazon

    zsetot = sum(soil%zse)
    totalmoist(:) = 0.0
    totalice(:) = 0.0
    DO k=1, ms
      totalmoist(:) = totalmoist(:) + ssoil%wb(:,k)*soil%zse(k)/zsetot
      totalice(:) = totalice(:) + ssoil%wbice(:,k)*soil%zse(k)/zsetot
    ENDDO

    Dtran=0.0
    WHERE ( canopy%fevc < 10.0 .and.  totalice  < 1.e-2)   Dtran=1.0
    
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
        frootX= max(0.01,max( veg%froot(:,k),veg%froot(:,j)))
        hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
                       *(veg%froot(:,k)*veg%froot(:,j))/(1-frootX) * Dtran
        hr_perTime(:,k,j) = hr_term(:,k,j)*1.0E-2/3600.0*dels ! m per timestep
        hr_perTime(:,j,k) = -1.0 * hr_perTime(:,k,j)
        hr_perTime(:,k,j) = hr_perTime(:,k,j)/soil%zse(k)
        hr_perTime(:,j,k) = hr_perTime(:,j,k)/soil%zse(j)
        ! Restricting changes to all broadleaf forests, and
        ! other forests and woody savannas in the tropics
        ! Note that veg types here are based on IGBP classification (BP mar2011)
!        WHERE (.NOT.(veg%iveg==2.OR.veg%iveg==7.OR.veg%iveg==8.OR.veg%iveg==9))
!        WHERE (.NOT.veg%iveg==2)
        ! MJT - use IGBP instead of CSIRO types
        !WHERE (.NOT.(veg%iveg == 2 .OR. veg%iveg == 7 ))
        WHERE (.NOT.veg%iveg==2)
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

    WHERE ( met%tk < tfrz + 5.  ) Dtran=0.0
      DO k=1, ms
        S_VG(:,k) = MIN( 1.0, MAX(1.0E-4, ssoil%wb(:,k) - soil%swilt) &
                    / (soil%ssat - soil%swilt) )
        wpsy(:,k) = -1.0/alpha_VG*(S_VG(:,k)**(-1.0/m_VG)-1.0)**(1/n_VG) &
                     *100*1.0E-6     ! VG model, convert from cm to Pa by (*100),
           !                           to MPa (*1.0E-6)
        C_hr(:,k) = 1./(1+(wpsy(:,k)/wpsy50)**n_hr)
      ENDDO                                                                                                   

      hr_term(:,:,:) = 0.0    ! unit: cm h^{-1}
      hr_perTime(:,:,:) = 0.0  
      DO k = 1,ms-2
        DO j = k+1,ms-1
         temp(:)        = 0.0
         available(:)   = 0.0
         accommodate(:) = 0.0
         frootX= max(0.01,max( veg%froot(:,k),veg%froot(:,j)))
         hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
            *(max(0.01,veg%froot(:,k))*max(0.01,veg%froot(:,j)))/(1-frootX)*Dtran
         hr_perTime(:,k,j) = hr_term(:,k,j)*1.0E-2/3600.0*dels ! m per timestep
         hr_perTime(:,j,k) = -1.0 * hr_perTime(:,k,j)
         hr_perTime(:,k,j) = hr_perTime(:,k,j)/soil%zse(k)
         hr_perTime(:,j,k) = hr_perTime(:,j,k)/soil%zse(j)
        ! Restricting changes to all broadleaf forests, and
        ! other forests and woody savannas in the tropics
        ! Note that veg types here are based on IGBP classification (BP mar2011)
        !        WHERE (.NOT.(veg%iveg == 1 .OR. veg%iveg == 6 ))
        ! MJT - use IGBP instead of CSIRO types
        !WHERE (.NOT.(veg%iveg == 2 .OR. veg%iveg == 7 ))
        WHERE (.NOT.(veg%iveg == 2))
           hr_perTime(:,k,j) = 0.0
           hr_perTime(:,j,k) = 0.0
        ENDWHERE
        WHERE (hr_perTime(:,k,j) < 0.0)
           available(:)   = MAX(0.0, ssoil%wb(:,k)- soil%sfc(:))
           accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,j))
           temp(:) = MAX(hr_perTime(:,k,j), &
                         -1.0*wiltParam*available(:), &
                         -1.0*satuParam*accommodate(:)*soil%zse(j)/soil%zse(k))
           hr_perTime(:,k,j) = temp(:)
           hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)
        ELSEWHERE (hr_perTime(:,j,k) < 0.0)
           available(:)   = MAX(0.0, ssoil%wb(:,j)- soil%sfc(:))
           accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,k))
           temp(:) = MAX(hr_perTime(:,j,k), &
                    -1.0*wiltParam*available(:), &
                    -1.0*satuParam*accommodate(:)*soil%zse(k)/soil%zse(j))
           hr_perTime(:,j,k) = temp(:)
           hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)
        ENDWHERE
        ssoil%wb(:,k) = ssoil%wb(:,k) + hr_perTime(:,k,j)
        ssoil%wb(:,j) = ssoil%wb(:,j) + hr_perTime(:,j,k)
      ENDDO
     ENDDO
                           

!    DO k=1, ms
!      total2(:) = total2(:) + ssoil%wb(:,k)*soil%zse(k)/zsetot
!    ENDDO
!    print 20,(ssoil%wb(idjd1,k),k=1,4),totalmoist(idjd1),total2(idjd1)
!20  format(1x,'hydrre2',4f8.5,2x,2f8.5)


!    IF (MINVAL(totalmoist(:) - total2(:)) < -1.0e-6 .OR. &
!        MAXVAL(totalmoist(:) - total2(:)) >  1.0e-6) THEN
!      PRINT *, 'Maximum of ', MAXVAL(totalmoist(:) - total2(:)), ' at ', &
!                MAXLOC(totalmoist(:) - total2(:))
!      PRINT *, 'Minimum of ', MINVAL(totalmoist(:) - total2(:)), ' at ', &
!                MINLOC(totalmoist(:) - total2(:))
!      STOP
!    ENDIF

  END SUBROUTINE hydraulic_redistribution
  !+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++


END MODULE soil_snow_module
