! soil_snow.f90
!
! Soil and snow routines source file for CABLE, CSIRO land surface model
!
! Science development by Eva Kowalczyk, CSIRO Marine and Atmospheric Research
! 
! Fortran-95 coding by Harvey Davies, Gab Abramowitz and Martin Dix
! bugs to gabsun@gmail.com.

MODULE soil_snow_module
  USE physical_constants ! from cable_variables.f90
  USE define_types       ! from cable_variables.f90
  IMPLICIT NONE
  PRIVATE
  REAL(r_1), PARAMETER :: cgsnow = 2090.  ! specific heat capacity for snow
  REAL(r_1), PARAMETER :: csice = 2.100e3 ! specific heat capacity for ice
  REAL(r_1), PARAMETER :: cswat = 4.218e3 ! specific heat capacity for water
  REAL(r_1), PARAMETER :: hl = 2.5104e6   ! latent heat of evaporation
  REAL(r_1), PARAMETER :: hlf = 0.335e6   ! latent heat of fusion
  REAL(r_1), PARAMETER :: cp = 1004.64
  REAL(r_1), PARAMETER :: rhowat = 1000.  ! density of water
  ! This module contains the following subroutines:
  PUBLIC soil_snow ! must be available outside this module
  PRIVATE trimb, smoisturev, surfbv, stempv ! internal subroutines
  INTEGER(i_d), PARAMETER, PRIVATE	      :: idjd = 3179
CONTAINS

  !----------------------------------------------------------------------
  ! SUBROUTINE trimb
  !	 like trim, but work arrays in work3f
  !	 rhs initially contains rhs; leaves with answer (jlm)
  !	 n.b. this one does not assume b = 1-a-c
  !
  !	 this routine solves the system
  !	   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kmax-1
  !	   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)	       for k=1
  !	   and	 a(k)*u(k-1)+b(k)*u(k)=rhs(k)	       for k=kmax
  !
  !	 the Thomas algorithm is used for solving sets of linear equation

  SUBROUTINE trimb (a, b, c, rhs, kmax)
    REAL(r_2), DIMENSION(:,:), INTENT(IN) :: a ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(IN) :: b ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(IN) :: c ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(:,:), INTENT(INOUT) :: rhs ! right hand side of eq
    INTEGER(i_d), INTENT(IN)		     :: kmax ! number of discrete layers
    INTEGER(i_d)			     :: k ! do lloop counter
    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: e !
    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: g !
    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: temp !
    !
    e (:,1) = c (:,1) / b (:,1)
    DO k = 2, kmax - 1
       temp(:,k) = 1. / (b(:,k) - a(:,k) * e(:,k-1) )
       e(:,k) = c(:,k) * temp(:,k)
    END DO
    g (:,1) = rhs (:,1) / b (:,1)
    DO k = 2, kmax - 1
       g (:,k) = (rhs(:,k) - a(:,k) * g(:,k-1) ) * temp(:,k)
    END DO
    !
    !	 do back substitution to give answer now
    rhs(:,kmax) = (rhs(:,kmax) - &
         a(:,kmax) * g(:,kmax-1)) / (b(:,kmax) - a(:,kmax) * e(:,kmax-1))
    DO k = kmax - 1, 1, - 1
       rhs (:,k) = g (:,k) - e (:,k) * rhs (:,k + 1)
    END DO
  END SUBROUTINE trimb
  !-------------------------------------------------------------------------
  SUBROUTINE smoisturev (fwtop,dt,ktau,ssoil,soil)
    ! Solves implicit soil moisture equation
    REAL(r_1), DIMENSION(mp), INTENT(INOUT)   :: fwtop ! water flux into the surface (precip-evap)
    REAL(r_1), INTENT(IN)		      :: dt    ! time step size (s)
    INTEGER(i_d), INTENT(IN)		      :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    INTEGER(i_d), PARAMETER		      :: ntest = 0 ! 2 for funny pre-set for idjd
    INTEGER(i_d), PARAMETER		      :: nmeth = - 1 ! Values as follows:
    !					 -1 for simple implicit D
    !					 1 for fully implicit solution
    !					 2 for simpler implicit
    !					 3 for simple implicit D, explicit K jlm pref
    !					 4 for simple implicit D, implicit K
    !					 0 for simple implicit D, new jlm TVD K
    REAL(r_2), DIMENSION(mp,3*ms) :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,3*ms) :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,3*ms) :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp)      :: fact   ! 
    REAL(r_2), DIMENSION(mp)	  :: fact2  !
    REAL(r_1), DIMENSION(mp)      :: fluxhi ! 
    REAL(r_1), DIMENSION(mp)	  :: fluxlo ! 
    REAL(r_1), DIMENSION(mp)      :: hydss  ! hydrolic conductivity adjusted for ice
    INTEGER(i_d), DIMENSION(mp)	  :: iqmn   ! for testing purposes only
    INTEGER(i_d), DIMENSION(mp)	  :: iqmx   ! for testing purposes only
    INTEGER(i_d)	          :: k      !
    REAL(r_1), DIMENSION(mp)	  :: phi    ! 
    REAL(r_2), DIMENSION(mp)	  :: pwb    ! 
    REAL(r_2), DIMENSION(:),ALLOCATABLE,SAVE :: pwb_min !
    REAL(r_1), DIMENSION(mp)	  :: rat    !   
    REAL(r_1), DIMENSION(mp)	  :: speed_k    !
    REAL(r_2), DIMENSION(mp)      :: ssatcurr_k !
    REAL(r_1), DIMENSION(mp)	  :: wblfmn !
    REAL(r_1), DIMENSION(mp)	  :: wblfmx !
    REAL(r_2), DIMENSION(mp,ms+1) :: wbh    !
    REAL(r_2), DIMENSION(mp,ms+1) :: z1     !
    REAL(r_2), DIMENSION(mp,ms+1) :: z2
    REAL(r_2), DIMENSION(mp,ms+1) :: z3
    REAL(r_1), DIMENSION(mp,ms+1) :: z1mult
    REAL(r_1), DIMENSION(mp,0:ms) :: fluxh
    REAL(r_1), DIMENSION(mp,0:ms) :: delt
    REAL(r_2), DIMENSION(mp,0:ms) :: dtt
    REAL(r_2), DIMENSION(mp)	  :: pwb_wbh
    REAL(r_2), DIMENSION(mp,ms)	  :: ssatcurr
    REAL(r_1), DIMENSION(mp)	  :: totwba ! diagnostic
    REAL(r_1), DIMENSION(mp)      :: totwbb
    REAL(r_1), DIMENSION(mp)      :: totwbc
    REAL(r_1), DIMENSION(mp)	  :: totwblb
    REAL(r_1), DIMENSION(mp)	  :: totwblc
    REAL(r_1), DIMENSION(mp)      :: wbficemx
    REAL(r_2), DIMENSION(mp)      :: wbh_k
    REAL(r_2), DIMENSION(mp)	  :: wbl_k
    REAL(r_2), DIMENSION(mp)      :: wbl_kp
    REAL(r_2), DIMENSION(mp)	  :: wh
    REAL(r_2), DIMENSION(mp)	  :: z3_k
    REAL(r_1), DIMENSION(mp)	  :: zsetot
    LOGICAL :: is_open     ! Is file open?
    INTEGER :: u	       ! I/O unit
    !
    IF (ktau == 1) THEN
       ALLOCATE(pwb_min(mp))
       ! Set working variable:
       pwb_min = (soil%swilt / soil%ssat ) **soil%ibp2
    END IF
    ! Block below for testing purposes only: - - - - - - - - - - - -
    IF (ntest > 0) THEN
         PRINT * , 'entering smoisturev fwtop,i2bp3,swilt,sfc,ssat: ', &
         fwtop(idjd), soil%i2bp3(idjd) , soil%swilt(idjd) , soil%sfc(idjd) , soil%ssat(idjd)
       u = 97
       inquire (u, opened=is_open)
       if (.not. is_open) then
          open (u, file='f97.txt', status='replace')
          write (u, *) 'ktau', ' fwtop', ' i2bp3', ' swilt', ' sfc', ' ssat'
       END IF
       write (u, *) ktau, fwtop, soil%i2bp3 , soil%swilt , soil%sfc , soil%ssat
       IF (ntest == 2) THEN ! just to test conservation
          IF (ktau == 1) ssoil%wb(:,ms) = soil%swilt
          fwtop = 0.
       END IF
       WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k) , k = 1, ms)
       WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k) , k = 1, ms)
       totwba = 0.
       DO k = 1, ms
          totwba = totwba + soil%zse(k) * ssoil%wb(:,k)
       END DO
    END IF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
    at = 0.
    bt = 1.
    ct = 0.
    z1mult(:,1) = 0.	! corresponds to 2b+3
    z1mult(:,ms + 1) = 0.  ! corresponds to 2b+3
    z1(:,1) = 0.  !  i.e. K(.5),	value at surface
    z1(:,ms + 1) = 0.  !  i.e. K(ms+.5), value at bottom
    ! nmeth: equation solution technique
    IF (nmeth <= 0) THEN
       !	  jlm split TVD version
       ! all land points
       delt (:,0) = 0.
       fluxh (:,0) = 0.
       fluxh (:,ms) = 0.
       DO k = 1, ms - 1
          ! Calculate amount of liquid soil water:
          wbl_k = max (0.01, ssoil%wb(:,k) - ssoil%wbice(:,k) )
          wbl_kp = max (0.01, ssoil%wb(:,k + 1) - ssoil%wbice(:,k + 1) )
          ! Calculate difference in liq soil water b/w consecutive layers:
          delt (:,k) = wbl_kp - wbl_k
          !	    especially to allow for isolated frozen layers, use min speed
          wh = min (wbl_k, wbl_kp)
          !	    with 50% wbice, reduce hyds by 1.e-5
          ! Calculate hyd conductivity adjusted for ice:
          hydss = soil%hyds * (1. - min (2. * ssoil%wbice(:,k) / max (0.01, &
               ssoil%wb(:,k) ), .99999) )
          speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
          !	    update wb by TVD method
          rat = delt (:,k - 1) / (delt (:,k) + sign (1.e-20, delt (:,k) ) )
          phi = max (0., min (1., 2. * rat), min (2., rat) ) ! 0 for -ve rat
          fluxhi = wh
          fluxlo = wbl_k
          ! Block below for testing purposes only:
          IF (ntest > 0) THEN
              PRINT * , 'in TVD for k= ', k
              PRINT * , 'wbl,wh,hydss ', wbl_k(idjd), wh(idjd), hydss(idjd)
              PRINT * , 'speeda,speedb,fluxhi,fluxlo,delt,rat,phi ', &
                   speed_k(idjd), .5 * soil%zse(k) / dt, fluxhi(idjd), fluxlo(idjd), delt (idjd,k) , rat(idjd), phi(idjd)
          END IF
          !  scale speed to grid lengths per dt & limit speed for stability
          !  1. OK too for stability
          speed_k = min (speed_k, .5 * soil%zse(k) / dt)
          fluxh (:,k) = speed_k * (fluxlo + phi * (fluxhi - fluxlo) )
       END DO
       !
       !	 update wb by TVD method
       DO k = ms, 1, - 1
          IF (nmeth ==  - 1) THEN ! each new wb constrained by ssat
             fluxh (:,k - 1) = min (fluxh (:,k - 1), (soil%ssat &
                  - ssoil%wb(:,k) ) * soil%zse(k) / dt + fluxh (:,k) )
          END IF
          ssoil%wb(:,k) = ssoil%wb(:,k) + dt * (fluxh (:,k-1) - fluxh (:,k) ) / soil%zse(k)
          !	re-calculate wblf
          ssatcurr_k = soil%ssat - ssoil%wbice(:,k)
          dtt(:,k) = dt / (soil%zse(k) * ssatcurr_k)
          !	this defn of wblf has different meaning from previous one in surfbv
          !	N.B. are imposing wbice<wb, so wblf <1
          ssoil%wblf(:,k) = (ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr_k
       END DO
       !
       ! wbh_k represents wblf(k-.5)
       DO k = 2, ms
          ssatcurr_k = soil%ssat - ssoil%wbice(:,k)
          wbh_k = (soil%zse(k) * ssoil%wblf(:,k - 1) + soil%zse(k - 1) &
               * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k - 1) )
          ! i.e. wbh**(bch+1)
          fact = wbh_k** (soil%ibp2 - 1)
          !	with 50% wbice, reduce hbsh by 1.e-5
          pwb_wbh = (soil%hsbh * (1. - min (2. * min (0.99, max ( &
               ssoil%wbice(:,k-1) / max (0.01, ssoil%wb(:,k-1) ), &
               ssoil%wbice(:,k)   / max (0.01, ssoil%wb(:,k) )) ) &
               , .99999) )) &
               * max (pwb_min, wbh_k * fact)
          !	moisture diffusivity (D) is  wbh*pwb; hsbh includes b
          !  i.e. D(k-.5)/soil%zshh(k)
          z3_k = pwb_wbh / soil%zshh (k)
          !	     PRINT * , 'z3_k ', z3_k
          ! where dtt=dt/(soil%zse(k)*ssatcurr_k)
          at (:,k) = - dtt(:,k) * z3_k
          ct (:,k - 1) = - dtt(:,k - 1) * z3_k
       END DO
       bt = 1. - at - ct
       ! Block below for testing purposes only:
       IF (ntest > 0) THEN
          PRINT * , 'midway through nmeth<=0'
          PRINT * , 'fluxh ', (fluxh (idjd,k) , k = 1, ms)
          WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k) , k = 1, ms)
          WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k) , k = 1, ms)
          totwbb = 0.
          totwblb = 0.
          DO k = 1, ms
             totwbb = totwbb + soil%zse(k) * ssoil%wb(:,k) ! diagnostic
             totwblb = totwblb + soil%zse(k) * ssoil%wblf(:,k) ! diagnostic
          END DO
          PRINT * , 'nmeth, b+2, 2b+3: ', nmeth, soil%ibp2(idjd) , soil%i2bp3(idjd)
          WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k) , k = 1, ms)
          WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k) , k = 1, ms)
          WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k) , k = 1, ms)
          PRINT * , 'zse ', soil%zse
          PRINT * , 'zshh ', soil%zshh
          PRINT * , 'dtt ', (dtt (idjd, k) , k = 1, ms)
          PRINT * , 'at ', (at (idjd, k) , k = 1, ms)
          PRINT * , 'bt ', (bt (idjd, k) , k = 1, ms)
          PRINT * , 'ct ', (ct (idjd, k) , k = 1, ms)
       END IF
       ssoil%wblf(:,1) = ssoil%wblf(:,1) + dtt(:,1) * fwtop / rhowat
    END IF
    IF (nmeth > 0) THEN
       wbficemx = 0.
       DO k = 1, ms
          ssatcurr(:,k) = soil%ssat - ssoil%wbice(:,k)
          !	    this defn of wblf has different meaning from previous one in surfbv
          !	    N.B. are imposing wbice<wb, so wblf <1
          ssoil%wblf(:,k) = (ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr(:,k)
          ssoil%wbfice(:,k) = ssoil%wbice(:,k) / soil%ssat
          wbficemx = max (wbficemx, ssoil%wbfice(:,k) )
          dtt(:,k) = dt / (soil%zse(k) * ssatcurr(:,k) )
       END DO
       IF (nmeth == 1) THEN ! full implicit method
          DO k = 2, ms
             !	  wbh(k)=min(1.,ww(k)*wblf(:,k-1)+(1.-ww(k))*wblf(:,k))
             !	  jlm: this is same as:
             wbh(:,k) = (soil%zse(k) * ssoil%wblf(:,k - 1) + soil%zse(k - 1) &
                  * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k - 1) )
             fact = wbh(:,k) ** (soil%ibp2 - 1) ! i.e. wbh**(bch+1)
             fact2 = fact * fact
             pwb = soil%hsbh * fact
             !	  moisture diffusivity (D) is  wbh*pwb
             !	  other term (K) is wbh*soil%hyds*fact2
             z1(:,k) = wbh(:,k) * ( (soil%i2bp3 - 1) * soil%hyds * fact2 - &
                  soil%ibp2 * pwb * (ssoil%wblf(:,k) - ssoil%wblf(:,k - 1) ) &
                  / soil%zshh (k) )
             z2(:,k) = - soil%i2bp3 * soil%hyds * fact2 + soil%ibp2 &
                  * pwb * (ssoil%wblf(:,k) - ssoil%wblf(:,k - 1) ) / soil%zshh (k)
             z3(:,k) = pwb * wbh(:,k) / soil%zshh (k)
             at (:,k) = dtt(:,k) * (z2(:,k) * .5 * soil%zse(k) / soil%zshh (k) - z3(:,k) )
          END DO
          DO k = 1, ms - 1
             ct (:,k) = dtt(:,k) * ( - z2(:,k + 1) * .5 * soil%zse(k) &
                  / soil%zshh (k + 1) - z3(:,k + 1) )
             bt (:,k) = 1. + dtt(:,k) * ( - z2(:,k + 1) * .5 * soil%zse( &
                  k + 1) / soil%zshh (k + 1) + z2(:,k) * .5 * soil%zse(max(k-1,1)) &
                  / soil%zshh (k) + z3(:,k + 1) + z3(:,k) )
          END DO
          bt (:,ms) = 1. + dtt(:,ms) * (z2(:,ms) * .5 * soil%zse(ms) &
               / soil%zshh (ms) + z3(:,ms) )
          DO k = 1, ms
             ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) * (z1(:,k + 1) - z1(:,k) )
          END DO
       END IF
       IF (nmeth >= 2) THEN ! part implicit method
          DO k = 2, ms
             z1mult(:,k) = soil%i2bp3 ! corresponds to 2b+3
          END DO
          DO k = 2, ms ! wbh(k) represents wblf(k-.5)
             wbh(:,k) = (soil%zse(k) * ssoil%wblf(:,k - 1) + soil%zse(k - 1) &
                  * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k - 1) )
             fact = wbh(:,k) ** (soil%ibp2 - 1) ! i.e. wbh**(bch+1)
             IF (nmeth == 2) pwb_wbh = soil%hsbh * wbh(:,k) * fact
             IF (nmeth >= 3) pwb_wbh = soil%hsbh * max (pwb_min, wbh(:,k) * fact)
             fact2 = fact * fact
             !	  moisture diffusivity (D) is  wbh*pwb
             !	  other term (K) is wbh*soil%hyds*fact2
             z1(:,k) = soil%hyds * fact2 !  i.e. K(k-.5)/wbh(:,k)
             z3(:,k) = pwb_wbh / soil%zshh (k) !  i.e. D(k-.5)/soil%zshh(k)
             at (:,k) = - dtt(:,k) * z3(:,k)
             ct (:,k - 1) = - dtt(:,k - 1) * z3(:,k)
          END DO
          bt = 1. - at - ct
          IF (nmeth == 4) THEN ! for simple implicit D, implicit K
             bt (:,1) = bt (:,1) + dtt(:,1) * z1mult(:,1 + 1) &
                  * z1(:,1 + 1) * soil%zse(1 + 1) / (soil%zse(1) + soil%zse(1 + 1) )
             DO k = 2, ms
                at (:,k) = at (:,k) - dtt(:,k) * z1mult(:,k) &
                     * z1(:,k) * soil%zse(k) / (soil%zse(k) + soil%zse(k - 1) )
                ct (:,k - 1) = ct (:,k - 1) + dtt(:,k - 1) * z1mult(:,k) &
                     * z1(:,k) * soil%zse(k - 1) / (soil%zse(k) + soil%zse(k-1) )
                bt (:,k) = bt (:,k) - dtt(:,k) * z1mult(:,k) &
                     * z1(:,k) * soil%zse(k - 1) / (soil%zse(k) + soil%zse(k - 1) ) &
                     + dtt(:,k) * z1mult(:,k + 1) * z1(:,k + 1) * soil%zse(k + 1) &
                     / (soil%zse(k) + soil%zse(k + 1) )
             END DO
             ! (nmeth == 4)
          END IF
          DO k = 2, ms
             !  i.e. now K(k-.5)
             z1(:,k) = wbh(:,k) * z1(:,k)
          END DO
          !
          ! the following top & bottom b.c.'s will preserve a uniform column
          !	     z1(1) =z1(2)   ! simple dk/dz=0
          !	     z1(ms+1)=z1(ms) ! simple dk/dz=0
          !	N.B. z1 are here +ve
          z1(:,1) = min (z1(:,2), z1(:,ms) )
          z1(:,ms + 1) = z1(:,1)
          ! no gravit. term if too much ice 11/12/00
          DO k = 1, ms
             IF (nmeth == 4) THEN
                WHERE (wbficemx < .75)
                   ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) * ( (z1mult(:,k+1) - 1.) &
                        * z1(:,k + 1) - (z1mult(:,k) - 1.) * z1(:,k) )
                END WHERE
             ELSE
                WHERE (wbficemx < .75)
                   ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) * (z1(:,k) - z1(:,k + 1) )
                END WHERE
             END IF
          END DO
       END IF
       IF (ntest > 0) THEN
          wblfmx = maxval(ssoil%wblf, 2)
          wblfmn = minval(ssoil%wblf, 2)
          iqmx   = maxloc(ssoil%wblf, 2)
          iqmn   = minloc(ssoil%wblf, 2)
       END IF
       ! Block below for testing purposes only:
       IF (ntest > 0) THEN
          totwbb = 0.
          totwblb = 0.
          DO k = 1, ms
             ! diagnostic
             totwbb = totwbb + soil%zse(k) * ssoil%wb(:,k)
             ! diagnostic
             totwblb = totwblb + soil%zse(k) * ssoil%wblf(:,k)
          END DO
          PRINT * , 'nmeth, b+2, 2b+3: ', nmeth, soil%ibp2(idjd) , soil%i2bp3(idjd)
          WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k) , k = 1, ms)
          WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k) , k = 1, ms)
          WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k) , k = 1, ms)
          WRITE (6, " ('wbh  ', 7f8.3) ") wbh(idjd,:)
          WRITE (6, " ('ssatcurr', 6f8.3) ") ssatcurr(idjd,:)
          PRINT * , 'pwb_wbh,pwb_min* for ms ', pwb_wbh(idjd), soil%hsbh(idjd) * pwb_min(idjd)
          PRINT * , 'wblfmx,wblfmn,iqmx,iqmn ', wblfmx(idjd), wblfmn(idjd), iqmx(idjd), iqmn(idjd)
          PRINT * , 'zse ', soil%zse
          PRINT * , 'zshh ', soil%zshh
          PRINT * , 'at ', (at (idjd,k) , k = 1, ms)
          PRINT * , 'bt ', (bt (idjd,k) , k = 1, ms)
          PRINT * , 'ct ', (ct (idjd,k) , k = 1, ms)
       END IF
       IF (nmeth == 3) THEN
          ! artificial fix applied here for safety (explicit nmeth only)
          DO k = 1, ms
             ssoil%wblf(:,k) = max (0., min (ssoil%wblf(:,k), 1.) )
          END DO
          ! (nmeth == 3)
       END IF
       ssoil%wblf(:,1) = ssoil%wblf(:,1) + dtt(:,1) * fwtop / rhowat
    END IF
    CALL trimb (at, bt, ct, ssoil%wblf, ms)
    DO k = 1, ms
       ssatcurr(:,k) = soil%ssat - ssoil%wbice(:,k)
       ssoil%wb(:,k) = ssoil%wblf(:,k) * ssatcurr(:,k) + ssoil%wbice(:,k)
       ssoil%wbice(:,k) = min (ssoil%wbice(:,k), .99 * ssoil%wb(:,k) )
    END DO
    ! Block below for testing purposes only:
    IF (ntest > 0) THEN
       PRINT * , 'at end of smoisturev,fwtop ', fwtop(idjd)
       PRINT * , 'tgg ', (ssoil%tgg(idjd, k) , k = 1, ms)
       WRITE (6, " ('wb   ', 6f8.3) ")  (ssoil%wb(idjd,k) , k = 1, ms)
       WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(idjd,k) , k = 1, ms)
       WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(idjd,k) , k = 1, ms)
       totwbc = 0.
       totwblc = 0.
       zsetot = 0.
       DO k = 1, ms
          totwbc = totwbc + soil%zse(k) * ssoil%wb(:,k)
          totwblc = totwblc + soil%zse(k) * ssoil%wblf(:,k)
          zsetot = zsetot + soil%zse(k)
       END DO
       PRINT * , 'totwba,totwbb,totwbc ', totwba(idjd), totwbb(idjd), totwbc(idjd)
       PRINT * , 'totwblb,totwblc ', totwblb(idjd), totwblc(idjd)
       PRINT * , 'with totwbc/zsetot: ', totwbc(idjd) / zsetot(idjd)
    END IF
  END SUBROUTINE smoisturev

  SUBROUTINE surfbv (dt, ktau, canopy, met, ssoil, soil)
    ! Calculates amount of 
    REAL(r_1), INTENT(IN)	:: dt   ! integration time step (s)
    INTEGER(i_d), INTENT(IN)	:: ktau ! integration step number
    TYPE(canopy_type), INTENT(INOUT)	:: canopy ! vegetation variables
    TYPE(met_type), INTENT(INOUT)	:: met    ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    INTEGER(i_d), PARAMETER	:: ntest = 0 ! for snow diag prints
    INTEGER(i_d), PARAMETER	:: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    REAL(r_1), DIMENSION(mp)	:: alv ! Snow albedo for visible
    REAL(r_1), DIMENSION(mp)	:: alir ! Snow albedo for near infra-red
    REAL(r_1), PARAMETER	:: alvo  = 0.95 ! albedo for vis. on a new snow
    REAL(r_1), PARAMETER	:: aliro = 0.65 ! albedo for near-infr. on a new snow
    REAL(r_1), DIMENSION(mp)	:: ar1 ! crystal growth  (-ve)
    REAL(r_1), DIMENSION(mp)	:: ar2 ! freezing of melt water
    REAL(r_1), DIMENSION(mp)	:: ar3
    REAL(r_1), DIMENSION(mp)	:: dnsnow ! new snow albedo
    REAL(r_1), DIMENSION(mp)	:: dtau
    REAL(r_1), DIMENSION(mp)	:: evapsn
    REAL(r_1), DIMENSION(mp)	:: fage !age factor
    REAL(r_1), DIMENSION(mp)	:: fwtop
    REAL(r_1), DIMENSION(mp)	:: fracm
    REAL(r_1), DIMENSION(mp)	:: fzenm
    INTEGER(i_d)		:: k
    REAL(r_1), DIMENSION(mp)	:: rnof5
    REAL(r_1), DIMENSION(mp)	:: segg
    REAL(r_1), DIMENSION(mp)	:: sfact
    REAL(r_1), DIMENSION(mp)	:: sgamm
    REAL(r_1), DIMENSION(mp)	:: smasstot
    REAL(r_1), DIMENSION(mp)	:: smelt
    REAL(r_1), DIMENSION(mp)	:: snowflx
    REAL(r_1), DIMENSION(mp)	:: snr
    REAL(r_1), DIMENSION(mp)	:: snrat
    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
    REAL(r_1), DIMENSION(mp)	:: talb ! snow albedo
    REAL(r_1), DIMENSION(mp)	:: tmp ! temporary value
    REAL(r_1), DIMENSION(mp)	:: totwet
    REAL(r_1), DIMENSION(mp)	:: weting
    !
    IF (ntest > 0) THEN
       PRINT * , 'entering surfbv	condxpr', canopy%precis(idjd)
       PRINT * , 'osnowd,snowd,isflag', ssoil%osnowd(idjd),ssoil%snowd(idjd),ssoil%isflag(idjd)
       PRINT * , 'tggsn ', (ssoil%tggsn (idjd, k) , k = 1, 3)
       PRINT * , 'tgg ', (ssoil%tgg (idjd, k) , k = 1, ms)
       PRINT * , 'wb ', (ssoil%wb (idjd, k) , k = 1, ms)
       PRINT * , 'wbice ', (ssoil%wbice(idjd,k) , k = 1, ms)
       PRINT * , 'gammzz ', (ssoil%gammzz(idjd,k) , k = 1, ms)
       PRINT * , 'albnew ', (ssoil%albsoilsn(idjd,k) , k = 1, 2)
    END IF
    ! Initialise runoff variables:
    ssoil%runoff = 0.
    ssoil%rnof1 = 0.
    ssoil%rnof2 = 0.
    smelt = 0. ! initialise snowmelt
    ssoil%osnowd = ssoil%snowd ! 
    !
    ! just using ncondxpr=1 treatment now
    WHERE (canopy%precis > 0.)
       ssoil%snowd = max (ssoil%snowd + met%precips, 0.)
       canopy%precis = canopy%precis - met%precips
       WHERE (ssoil%isflag == 0)
          WHERE (canopy%precis > 0. .and. ssoil%tgg(:,1) < tfrz)
             ssoil%snowd = max (ssoil%snowd + canopy%precis, 0.)
             ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis * hlf / ssoil%gammzz(:,1)
             canopy%precis = 0.
          END WHERE
       ELSEWHERE ! i.e. ssoil%isflag=1
          WHERE (canopy%precis > 0.)
             ssoil%snowd = max (ssoil%snowd + canopy%precis, 0.)
             sgamm = ssoil%ssdn(:,1) * 2105. * ssoil%sdepth(:,1)
             ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + canopy%precis * hlf &
                  * ssoil%smass(:,1) / (sgamm * ssoil%osnowd )
             sgamm = ssoil%ssdn(:,2) * 2105. * ssoil%sdepth(:,2)
             ssoil%tggsn(:,2) = ssoil%tggsn(:,2) + canopy%precis * hlf &
                  * ssoil%smass(:,2) / (sgamm * ssoil%osnowd )
             sgamm = ssoil%ssdn(:,3) * 2105. * ssoil%sdepth(:,3)
             ssoil%tggsn(:,3) = ssoil%tggsn(:,3) + canopy%precis * hlf &
                  * ssoil%smass(:,3) / (sgamm * ssoil%osnowd )
             canopy%precis = 0.
          END WHERE
       END WHERE
    END WHERE
    ! (Potentially) reduce soil evaporation for reduction in available water 
    ! due to transpiration:
    WHERE (ssoil%snowd < 0.1 .and. canopy%fes .gt. 0. )
       canopy%fes= min(canopy%fes,max(0.,(ssoil%wb(:,1)-soil%swilt))* soil%zse(1) &
            * 1000. * hl / dt)
       canopy%fes = min(canopy%fes,(ssoil%wb(:,1)-ssoil%wbice(:,1)) * soil%zse(1) &
            * 1000. * hl / dt)
    END WHERE
    ! Calculate soil evaporation total in mm (from W/m2):
    segg = canopy%fes / hl
    ! Initialise snow evaporation:
    evapsn = 0
    ! Snow evaporation and melting
    WHERE (ssoil%snowd > .1)
       evapsn = min (ssoil%snowd, dt * canopy%fes / (hl + hlf) )
       ssoil%snowd = ssoil%snowd - evapsn
       segg = 0.
       canopy%fes = evapsn * (hl + hlf) / dt ! return for hyd. balance
    END WHERE
    WHERE (ssoil%snowd > 0.)
       WHERE (ssoil%isflag == 0)
          !	      snow covered land
          !	      following done in sflux  via  ga= ... +cls*egg + ...
          WHERE (ssoil%tgg(:,1) >= tfrz)
             !**		land,snow,melting
             snowflx = (ssoil%tgg(:,1) - tfrz) * ssoil%gammzz(:,1)
             !		prevent snow depth going negative
             smelt = min (snowflx / hlf, ssoil%snowd )
             ssoil%snowd = ssoil%snowd - smelt
             ssoil%tgg(:,1) = ssoil%tgg(:,1) - smelt * hlf / ssoil%gammzz(:,1)
          END WHERE
       ELSEWHERE ! 3-layer scheme,  isflag=1
          smelt1(:,0) = 0.0
       END WHERE
    END WHERE
    DO k = 1, 3
       WHERE (ssoil%snowd > 0.0 .and. ssoil%isflag > 0)
          sgamm = ssoil%ssdn(:,k) * 2105. * ssoil%sdepth(:,k)
          !		WHERE (soil%isoilm == 9)
          ! snow melt refreezing
          snowflx = smelt1(:,k - 1) * hlf / dt
          ssoil%tggsn(:,k) = ssoil%tggsn(:,k) + snowflx * dt / sgamm
          ! increase density due to snowmelt
          ssoil%ssdn(:,k) = min ( (ssoil%smass(:,k) + smelt1(:,k - 1) ) &
               / (ssoil%smass(:,k) / ssoil%ssdn(:,k) + smelt1(:,k - 1) / 1000.), 450.)
          ssoil%smass(:,k) = ssoil%smass(:,k) + smelt1(:,k - 1)
          sgamm = ssoil%smass(:,k) * 2105.
          smelt1(:,k - 1) = 0.
          !		END WHERE
          smelt1(:,k) = 0.
          WHERE (ssoil%tggsn(:,k) > tfrz)
             snowflx = (ssoil%tggsn(:,k) - tfrz) * sgamm
             smelt1(:,k) = min (snowflx / hlf, 0.9 * ssoil%smass(:,k) )
             ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
             ssoil%tggsn(:,k) = min (ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm, tfrz)
          END WHERE
          ssoil%snowd = ssoil%snowd - smelt
       END WHERE
    END DO
    WHERE (ssoil%snowd > 0.0 .and. ssoil%isflag > 0)
       smelt = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
       ssoil%snowd = ssoil%snowd - smelt
    END WHERE
    totwet = canopy%precis + smelt
    ssoil%rnof1 = max (0., totwet - dt / 172.8) ! 86400/500 = 172.8
    weting = totwet - ssoil%rnof1
    ssoil%rnof1 = ssoil%rnof1 + max (0., &
         weting - .99 * min ( (soil%ssat - ssoil%wb(:,1) ) * soil%zse(1) * 1000., weting))
    weting = totwet - ssoil%rnof1
    fwtop = weting / dt - segg
    ! Diagnostic block below:
    IF (ntest > 0) THEN
       PRINT * , 'in surfbv before smoisturev  condxpr', canopy%precis(idjd)
       PRINT * , 'osnowd,snowd,isflag',ssoil%osnowd(idjd),ssoil%snowd(idjd),ssoil%isflag(idjd)
       PRINT * , 'tggsn_c ', (ssoil%tggsn (idjd, k) , k = 1, 3)
       PRINT * , 'tgg ', (ssoil%tgg (idjd, k) , k = 1, ms)
       PRINT * , 'smass ', (ssoil%smass(idjd,k) , k = 1, 3)
       PRINT * , 'ssdn ', (ssoil%ssdn(idjd,k) , k = 1, 3)
       PRINT * , 'fwtop ', fwtop(idjd)
    END IF
    CALL smoisturev (fwtop, dt, ktau, ssoil, soil)
    ! Diagnostic block below:
    IF (ntest > 0) THEN
       PRINT * , 'in surfbv after smoisturev '
       PRINT * , 'osnowd,snowd,isflag,ssat,runoff', ssoil%osnowd(idjd) , &
            ssoil%snowd(idjd),ssoil%isflag(idjd),soil%ssat(idjd),ssoil%runoff(idjd)
       PRINT * , 'tggsn_d ', (ssoil%tggsn (idjd, k) , k = 1, 3)
    END IF
    DO k = 1, ms
       ssoil%rnof1 = ssoil%rnof1 + max (ssoil%wb(:,k) - soil%ssat, 0.) * 1000. * soil%zse(k)
       ssoil%wb(:,k) = min (ssoil%wb(:,k), soil%ssat )
    END DO
    !  for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
    !  account for soil/ice cracking
    fracm = min (0.2, 1. - min (1., ssoil%wb(:,ms) / soil%sfc ) )
    ssoil%wb(:,ms) = ssoil%wb(:,ms) + fracm * ssoil%rnof1 / (1000. * soil%zse(ms))
    ssoil%rnof1 = (1. - fracm) * ssoil%rnof1
    tmp = max (min (ssoil%wb(:,ms) - soil%sfc, .99 * ssoil%wb(:,ms) &
         - ssoil%wbice(:,ms) ) * soil%c3 / 86400., 0.)
    ssoil%rnof2 = soil%zse(ms) * 1000. * tmp * dt
    ssoil%wb(:,ms) = ssoil%wb(:,ms) - tmp * dt
    ssoil%runoff = ssoil%rnof1 + ssoil%rnof2
    ssoil%wbtot = 0.
    DO k = 1, ms
       ssoil%wbtot = ssoil%wbtot + ssoil%wb(:,k) * 1000. * soil%zse(k)
    END DO
    !
    !---  glacier formation
    IF (nglacier == 2) THEN
       WHERE (ssoil%snowd > 1000.)
          rnof5 = ssoil%snowd - 1000.
          ssoil%runoff = ssoil%runoff + rnof5
          !----	   change local tg to account for energy - clearly not best method
          WHERE (ssoil%isflag == 0)
             smasstot = 0.0
             ssoil%tgg(:,1) = ssoil%tgg(:,1) - rnof5 * hlf / ssoil%gammzz(:,1)
             ssoil%snowd = 1000.
          ELSEWHERE
             smasstot = ssoil%smass(:,1) + ssoil%smass(:,2) + ssoil%smass(:,3)
          END WHERE
       END WHERE
       DO k = 1, 3
          WHERE (ssoil%snowd > 1000.  .and.  ssoil%isflag > 0)
             sgamm = ssoil%ssdn(:,k) * 2105. * ssoil%sdepth(:,k)
             smelt1(:,k) = min (rnof5 * ssoil%smass(:,k) / smasstot, 0.9 * ssoil%smass(:,k) )
             ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
             ssoil%snowd = ssoil%snowd - smelt1(:,k)
             ssoil%tggsn(:,k) = ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm
          END WHERE
       END DO
    END IF
    ! Diagnostic block below:
    IF (ntest > 0) THEN
       PRINT * , 'end surfbv  rnof1,runoff ', ssoil%rnof1(idjd), ssoil%runoff(idjd)
       sgamm = ssoil%ssdn(:,1) * 2105. * ssoil%sdepth(:,1)
       PRINT * , 'snowd,isflag,sgamm ', ssoil%snowd(idjd), ssoil%isflag(idjd),sgamm(idjd)
       PRINT * , 'tggsn_d ', (ssoil%tggsn (idjd, k) , k = 1, 3)
       PRINT * , 'tgg ', (ssoil%tgg (idjd, k) , k = 1, ms)
       PRINT * , 'wb ', (ssoil%wb(idjd,k) , k = 1, ms)
    END IF
    !	 calulate soil/snow albedo
    !	  cuvrf(i,1) = albsav(iq) ! use surface albedo from indata
    sfact = 0.68
    WHERE (soil%albsoil <= 0.14)
       sfact = 0.5
    ELSEWHERE (soil%albsoil > 0.14 .and. soil%albsoil <= 0.20)
       sfact = 0.62
    END WHERE
    ssoil%albsoilsn(:,2) = 2. * soil%albsoil / (1. + sfact)
    ssoil%albsoilsn(:,1) = sfact * ssoil%albsoilsn(:,2)
    !
    !	       new snow albedo (needs osnowd from the previous dt)
    dnsnow = min (1., .1 * max (0., ssoil%snowd - ssoil%osnowd ) ) ! new snow (cm H2O)
    !		Snow age depends on snow crystal growth, freezing of melt water,
    !		accumulation of dirt and amount of new snow.
    tmp = ssoil%isflag * ssoil%tggsn(:,1) + (1 - ssoil%isflag ) * ssoil%tgg(:,1)
    tmp = min (tmp, 273.15)
    ar1 = 5000. * (1. / 273.15 - 1. / tmp) ! crystal growth  (-ve)
    ar2 = 10. * ar1 ! freezing of melt water
    snr = ssoil%snowd / max (ssoil%ssdnn, 100.)
    ! fixes for Arctic & Antarctic
    WHERE (soil%isoilm == 9)
       ar3 = .0005
       dnsnow = max (dnsnow, .002) !increase refreshing of snow in Antarctic
       snrat = min (1., snr / (snr + .001) )
    ELSEWHERE
       ! accumulation of dirt
       ar3 = .1
       snrat = min (1., snr / (snr + .01) )
    END WHERE
    dtau = 1.e-6 * (exp(ar1) + exp(ar2) + ar3) * dt
    WHERE (ssoil%snowd <= 1.0)
       ssoil%snage = 0.
    ELSEWHERE
       ssoil%snage = max (0., (ssoil%snage + dtau) * (1. - dnsnow) )
    END WHERE
    fage = 1. - 1. / (1. + ssoil%snage ) !age factor
    !
    !	    Snow albedo is dependent on zenith angle and  snow age.
    !	    albedo zenith dependence
    !	    alvd = alvo * (1.0-cs*fage); alird = aliro * (1.-cn*fage)
    !		    where cs = 0.2, cn = 0.5, b = 2.0
    tmp = max (.17365, met%coszen )
    fzenm = max(merge(0.0, (1. + 1. / 2.) / (1. + 2. * 2. * tmp) - 1. / 2., tmp > 0.5), 0.)
    tmp = alvo * (1.0 - 0.2 * fage)
    alv = .4 * fzenm * (1. - tmp) + tmp
    tmp = aliro * (1. - .5 * fage)
    alir = .4 * fzenm * (1.0 - tmp) + tmp
    talb = .5 * (alv + alir) ! snow albedo
    ! alss = (1. - snrat) * soil%albsoil + snrat * talb ! canopy free surface albedo
    ssoil%albsoilsn(:,2) = (1. - snrat) * ssoil%albsoilsn(:,2) + snrat * alir
    ssoil%albsoilsn(:,1) = (1. - snrat) * ssoil%albsoilsn(:,1) + snrat * alv

    !  print *,'soilsn',soil%albsoil,ssoil%albsoilsn,talb
  END SUBROUTINE surfbv

  ! SUBROUTINE stempv
  !	 calculates temperatures of the soil
  !	 tgg - new soil/ice temperature
  !	 ga - heat flux from the atmosphere (ground heat flux)
  !	 ccnsw - soil thermal conductivity, including water/ice

  SUBROUTINE stempv(dt, canopy, ssoil, soil)
    REAL(r_1), INTENT(IN)		:: dt ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)	:: canopy
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    INTEGER(i_d), PARAMETER	        :: ntest = 0
    REAL(r_2), DIMENSION(mp, -2:ms)	:: at
    REAL(r_2), DIMENSION(mp, -2:ms)	:: bt
    REAL(r_2), DIMENSION(mp, -2:ms)	:: ct
    REAL(r_2), DIMENSION(mp,ms)		:: ccnsw  ! soil thermal conductivity (incl water/ice)
    REAL(r_1), DIMENSION(mp)		:: coefa
    REAL(r_1), DIMENSION(mp)		:: coefb
    REAL(r_2), DIMENSION(mp)		:: dtg
    REAL(r_2), DIMENSION(mp)		:: ew
    REAL(r_2), DIMENSION(mp,-2:ms+1)	:: coeff
    INTEGER(i_d)			:: k
    REAL(r_1), DIMENSION(mp,-2:ms)	:: rhs
    REAL(r_2), DIMENSION(mp,3)		:: sconds
    REAL(r_1), DIMENSION(mp)		:: sgamm
    REAL(r_2), DIMENSION(mp,ms+3)	:: tmp_mat ! temp. matrix for tggsn & tgg
    REAL(r_2), DIMENSION(mp)		:: xx
    !
    at = 0.
    bt = 1.
    ct = 0.
    coeff = 0.
    DO k = 1, ms
       WHERE (soil%isoilm == 9)
          ccnsw(:,k) = 1.5
       ELSEWHERE
          ew = ssoil%wblf(:,k) * soil%ssat
          ccnsw(:,k) = min (soil%cnsd * exp (ew * log (60.) + ssoil%wbfice(:,k) &
               * soil%ssat * log (250.) ), 1.5) * &
               max (1., sqrt (min (2., .5 * soil%ssat / min (ew, .5 * soil%ssat )) ) )
       END WHERE
    END DO
    WHERE (ssoil%isflag == 0)
       xx = max (0., ssoil%snowd / ssoil%ssdnn )
       ccnsw(:,1) = (ccnsw(:,1) - 0.2) * (soil%zse(1) / (soil%zse(1) + xx)) + 0.2
    END WHERE
    DO k = 3, ms
       WHERE (ssoil%isflag == 0)
          coeff (:,k) = 2.0 / (soil%zse(k-1) / ccnsw (:,k-1) + soil%zse(k) / ccnsw (:,k) )
       END WHERE
    END DO
    k = 1
    WHERE (ssoil%isflag == 0)
       coeff (:,2) = 2. / ( (soil%zse(1) + xx) / ccnsw (:,1) + soil%zse(2) / ccnsw (:,2) )
       coefa = 0.
       coefb = coeff (:,2)
       ssoil%gammzz(:,k) = max ( (1. - soil%ssat ) * soil%css * soil%rhosoil &
            + soil%ssat * (ssoil%wblf(:,k) * cswat * rhowat &
            + ssoil%wbfice(:,k) * csice * rhowat * .9) &
            , soil%css * soil%rhosoil ) * soil%zse(k)
       ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
       dtg = dt / ssoil%gammzz(:,k)
       at (:,k) = - dtg * coeff (:,k)
       ct (:,k) = - dtg * coeff (:,k + 1) ! c3(ms)=0 & not really used
       bt (:,k) = 1. - at (:,k) - ct (:,k)
    END WHERE
    DO k = 2, ms
       WHERE (ssoil%isflag == 0)
          ssoil%gammzz(:,k) = max ( (1. - soil%ssat ) * soil%css * soil%rhosoil &
               + soil%ssat * (ssoil%wblf(:,k) * cswat * rhowat + ssoil%wbfice(:,k) &
               * csice * rhowat * .9), soil%css * soil%rhosoil ) * soil%zse(k)
          dtg = dt / ssoil%gammzz(:,k)
          at (:,k) = - dtg * coeff (:,k)
          ct (:,k) = - dtg * coeff (:,k + 1) ! c3(ms)=0 & not really used
          bt (:,k) = 1. - at (:,k) - ct (:,k)
       END WHERE
    END DO
    WHERE (ssoil%isflag == 0)
       bt (:,1) = bt (:,1) - canopy%dgdtg * dt / ssoil%gammzz(:,1)
       ssoil%tgg(:,1) = ssoil%tgg(:,1) + (canopy%ga - ssoil%tgg(:,1) * canopy%dgdtg) &
            * dt / ssoil%gammzz(:,1)
    END WHERE
    IF (ntest > 0) THEN
       PRINT * , 'tgg1,ga,gammzz ', ssoil%tgg(idjd,1) , canopy%ga(idjd) , ssoil%gammzz(idjd,1)
       PRINT * , 'dgdtg,degdt,dfgdt ', canopy%dgdtg(idjd)
       PRINT * , 'ssat,css,rhos,cswat,rhowat,csice ', soil%ssat(idjd) , &
            soil%css(idjd) , soil%rhosoil(idjd) , cswat, rhowat, csice
       PRINT * , 'wblf1,wbfice1,zse1,cgsnow ', ssoil%wblf(idjd,1) , &
            ssoil%wbfice(idjd,1) , soil%zse(1) , cgsnow
       PRINT * , 'at ', (at (idjd,k) , k = 1, ms)
       PRINT * , 'bt ', (bt (idjd,k) , k = 1, ms)
       PRINT * , 'ct ', (ct (idjd,k) , k = 1, ms)
       PRINT * , 'rhs ', (ssoil%tgg(idjd,k) , k = 1, ms)
    END IF
    coeff (:,1 - 3) = 0.
    ! 3-layer snow points done here
    WHERE (ssoil%isflag /= 0)
       sconds(:,1) = max (0.2, min (2.576e-6 * ssoil%ssdn(:,1) ** 2 + .074, 1.) )
       sconds(:,2) = max (0.2, min (2.576e-6 * ssoil%ssdn(:,2) ** 2 + .074, 1.) )
       sconds(:,3) = max (0.2, min (2.576e-6 * ssoil%ssdn(:,3) ** 2 + .074, 1.) )
       coeff(:,-1) = 2. / (ssoil%sdepth(:,1) / sconds(:,1) + ssoil%sdepth(:,2) / sconds(:,2) )
       coeff(:,0)	= 2. / (ssoil%sdepth(:,2) / sconds(:,2) + ssoil%sdepth(:,3) / sconds(:,3) )
       coeff(:,1)	= 2. / (ssoil%sdepth(:,3) / sconds(:,3) + soil%zse(1) / ccnsw (:,1) )
    END WHERE
    DO k = 2, ms
       WHERE (ssoil%isflag /= 0)
          coeff (:,k) = 2. / (soil%zse(k-1) / ccnsw (:,k-1) + soil%zse(k) / ccnsw (:,k) )
       END WHERE
    END DO
    WHERE (ssoil%isflag /= 0)
       coefa = coeff (:,-1)
       coefb = coeff (:,1)
    END WHERE
    DO k = 1, 3
       WHERE (ssoil%isflag /= 0)
          sgamm = ssoil%ssdn(:,k) * 2105. * ssoil%sdepth(:,k)
          dtg = dt / sgamm
          !	       rhs(k-3) = ssoil%tggsn(:,k)	  ! A
          at (:,k - 3) = - dtg * coeff (:,k - 3)
          ct (:,k - 3) = - dtg * coeff (:,k - 2)
          bt (:,k - 3) = 1. - at (:,k - 3) - ct (:,k - 3)
       END WHERE
    END DO
    DO k = 1, ms
       WHERE (ssoil%isflag /= 0)
          ssoil%gammzz(:,k) = max ( (1. - soil%ssat ) * soil%css * soil%rhosoil &
               + soil%ssat * (ssoil%wblf(:,k) * cswat * rhowat + ssoil%wbfice(:,k) &
               * csice * rhowat * .9), soil%css * soil%rhosoil ) * soil%zse(k)
          dtg = dt / ssoil%gammzz(:,k)
          at (:,k) = - dtg * coeff (:,k)
          ct (:,k) = - dtg * coeff (:,k + 1) ! c3(ms)=0 & not really used
          bt (:,k) = 1. - at (:,k) - ct (:,k)
       END WHERE
    END DO
    WHERE (ssoil%isflag /= 0)
       sgamm = ssoil%ssdn(:,1) * 2105. * ssoil%sdepth(:,1)
       !	    rhs(1-3) = rhs(1-3)+canopy%ga*dt/sgamm
       !	    new code
       bt (:,- 2) = bt (:,- 2) - canopy%dgdtg * dt / sgamm
       ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + (canopy%ga - ssoil%tggsn(:,1) &
            * canopy%dgdtg ) * dt / sgamm
       rhs(:,1 - 3) = ssoil%tggsn(:,1)
    END WHERE
    IF (ssoil%isflag(idjd) /= 0 .and.  ntest > 0) THEN
       PRINT * , 'in stempv 3-layer snow code '
       PRINT * , 'ccnsw ', (ccnsw (idjd,k) , k = 1, ms)
       PRINT * , 'sdepth d ', (ssoil%sdepth(idjd,k) , k = 1, 3)
       PRINT * , 'sconds ', sconds
       PRINT * , 'coeff ', coeff
       PRINT * , 'at ', (at (idjd,k) , k = - 2, ms)
       PRINT * , 'bt ', (bt (idjd,k) , k = - 2, ms)
       PRINT * , 'ct ', (ct (idjd,k) , k = - 2, ms)
       PRINT * , 'rhs(tggsn,tgg) ', (ssoil%tggsn(idjd,k), k = 1,3), (ssoil%tgg(idjd,k), k = 1,ms)
       PRINT * , 'tggsn,ga,sgamm ', ssoil%tggsn (idjd,1) , canopy%ga(idjd) ,sgamm(idjd)
       PRINT * , 'dgdtg ', canopy%dgdtg(idjd)
    END IF
    !
    !     note in the following that tgg and tggsn are processed together
    tmp_mat(:,:3) = ssoil%tggsn
    tmp_mat(:,4:) = ssoil%tgg
    CALL trimb (at, bt, ct, tmp_mat, ms + 3)
    ssoil%tggsn = tmp_mat(:,:3)
    ssoil%tgg   = tmp_mat(:,4:)
    canopy%sghflux = coefa * (ssoil%tggsn(:,1) - ssoil%tggsn(:,2) )
    canopy%ghflux = coefb * (ssoil%tgg(:,1) - ssoil%tgg(:,2) ) ! +ve downwards
  END SUBROUTINE stempv

  !--------------------------------------------------------------------------------
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

  SUBROUTINE soil_snow(dt, ktau, soil, ssoil, canopy, met)
    REAL(r_1), INTENT(IN)	:: dt ! integration time step (s)
    INTEGER(i_d), INTENT(IN)	:: ktau ! integration step number
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(canopy_type), INTENT(INOUT)	     :: canopy
    TYPE(met_type), INTENT(INOUT)	     :: met ! all met forcing
    INTEGER(i_d), PARAMETER	:: ntest = 0 !  for prints
    INTEGER(i_d)		:: k
    REAL(r_2), DIMENSION(mp)	:: excd
    REAL(r_2), DIMENSION(mp)	:: excm
    REAL(r_1), DIMENSION(mp)	:: osm
    REAL(r_2), DIMENSION(mp)	:: sicefreeze
    REAL(r_2), DIMENSION(mp)	:: sicemelt
    REAL(r_1), DIMENSION(mp)	:: sd1
    REAL(r_1), DIMENSION(mp)	:: sm1
    REAL(r_1), DIMENSION(mp)	:: tr1
    REAL(r_1), DIMENSION(mp,ms) :: evapfbl
    REAL(r_1)			:: snmin = 0.11 ! 1000. for 1-layer;
    !						~.11 to turn on 3-layer snow
    ! Diagnostic block:
    IF(ntest>0) THEN
       PRINT *,'in soilsnowv before stempv,  ktau= ',ktau,ssoil%isflag(idjd)
       PRINT *,'soilsnowv ',canopy%dgdtg(idjd),canopy%fevc(idjd),canopy%fes(idjd),canopy%precis(idjd)
       PRINT *,'ga,dt,ssdn ',canopy%ga(idjd),dt,(ssoil%ssdn(idjd,k),k=1,3)
       PRINT *,'osnowd,snowd,isflag', ssoil%osnowd(idjd),ssoil%snowd(idjd),ssoil%isflag(idjd)
       PRINT *,'tgg ',(ssoil%tgg(idjd,k),k=1,ms)
       PRINT *,'tggsn ',(ssoil%tggsn(idjd,k),k=1,3)
       PRINT *,'wb ',(ssoil%wb(idjd,k),k=1,ms)
       PRINT *,'wbice ',(ssoil%wbice(idjd,k),k=1,ms)
    END IF

    IF (ktau <= 1) THEN
       ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
    !   IF (ntest == 3) snmin = .11 ! to force 3-layer snow for testing
       DO k = 1, ms
          WHERE (ssoil%tgg(:,k) <= tfrz .and. ssoil%wbice(:,k) <= 0.01)
             ssoil%wbice(:,k) = 0.1 * ssoil%wb(:,k)
          END WHERE
          WHERE (ssoil%tgg(:,k) < tfrz)
             ssoil%wbice(:,k) = 0.99 * ssoil%wb(:,k)
          END WHERE
       END DO
    END IF
    DO k = 1, ms ! for stempv
       ssoil%wblf(:,k) = max (0.01, (ssoil%wb(:,k) - ssoil%wbice(:,k) ) ) / soil%ssat
       ssoil%wbfice(:,k) = ssoil%wbice(:,k) / soil%ssat
    END DO
    tr1 = 0. ! Requires this initialisation with Intel compiler
    WHERE (ssoil%snowd <= 0.)
       ssoil%isflag = 0
       ssoil%ssdn(:,1) = 140.
       ssoil%ssdnn = 140.
       ssoil%tggsn(:,1) = tfrz
       ssoil%tggsn(:,2) = tfrz
       ssoil%tggsn(:,3) = tfrz
       ssoil%sdepth(:,1) = 0.
       ssoil%sdepth(:,2) = 0.
       ssoil%sdepth(:,3) = 0.
    ELSEWHERE (ssoil%snowd < snmin * ssoil%ssdnn)
       WHERE (ssoil%isflag == 1)
          ssoil%ssdn(:,1) = ssoil%ssdnn
       END WHERE
       ssoil%ssdn(:,1) = max (140., ssoil%ssdn(:,1) + dt * ssoil%ssdn(:,1) &
            * 2.8e-6 * exp ( - .03 * &
            (273.15 - min(tfrz, merge(ssoil%tggsn(:,1), ssoil%tgg(:,1), &
            ssoil%isflag == 1))) &
            - merge(0.046, 0.0, ssoil%ssdn(:,1) >= 150.) &
            * (ssoil%ssdn(:,1) - 150.) ) )
       ssoil%ssdn(:,1) = (140.0 - ssoil%ssdn(:,1)) * max(0., 1. &
            - ssoil%osnowd / ssoil%snowd ) + ssoil%ssdn(:,1)
       ssoil%ssdnn = ssoil%ssdn(:,1)
       ssoil%isflag = 0
       ssoil%tggsn(:,1) = tfrz
       ssoil%tggsn(:,2) = tfrz
       ssoil%tggsn(:,3) = tfrz
       ssoil%sdepth(:,1) = ssoil%snowd / ssoil%ssdn(:,1)
    ELSEWHERE ! sufficient snow now
       WHERE (ssoil%isflag == 0)
          ssoil%tggsn(:,1) = ssoil%tgg(:,1)
          ssoil%tggsn(:,2) = ssoil%tgg(:,2)
          ssoil%tggsn(:,3) = ssoil%tgg(:,3)
          ssoil%ssdn(:,2) = ssoil%ssdn(:,1)
          ssoil%ssdn(:,3) = ssoil%ssdn(:,1)
          ssoil%sdepth(:,1) = .07
          ssoil%sdepth(:,2) = max (.02, (ssoil%snowd / ssoil%ssdn(:,1) - .07) &
               * merge(0.3, 0.45, ssoil%snowd > 20.0))
          ssoil%sdepth(:,3) = max (.02, (ssoil%snowd / ssoil%ssdn(:,1) - .07) &
               * merge(0.7, 0.55, ssoil%snowd > 20.0))
          ssoil%smass(:,1) = .07 * ssoil%ssdn(:,1)
          ssoil%smass(:,2) = ssoil%sdepth(:,2) * ssoil%ssdn(:,2)
          ssoil%smass(:,3) = ssoil%sdepth(:,3) * ssoil%ssdn(:,3)
       END WHERE
       ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dt * ssoil%ssdn(:,1) * 3.1e-6 * &
            exp( - .03 * (273.15 - min(tfrz, ssoil%tggsn(:,1))) &
            - merge(0.046, 0.0, ssoil%ssdn(:,1) >= 150.) &
            * (ssoil%ssdn(:,1) - 150.) )
       ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dt * ssoil%ssdn(:,2) * 3.1e-6 * &
            exp( - .03 * (273.15 - min(tfrz, ssoil%tggsn(:,2))) &
            - merge(0.046, 0.0, ssoil%ssdn(:,1) >= 150.) &
            * (ssoil%ssdn(:,2) - 150.) )
       ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dt * ssoil%ssdn(:,3) * 3.1e-6 * &
            exp( - .03 * (273.15 - min(tfrz, ssoil%tggsn(:,3))) &
            - merge(0.046, 0.0, ssoil%ssdn(:,1) >= 150.) &
            * (ssoil%ssdn(:,3) - 150.) )
       ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dt * 9.806 * .5 * .07 * ssoil%ssdn(:,1) &
            * ssoil%ssdn(:,1) &
            / (3.e7 * exp (.021 * ssoil%ssdn(:,1) + .081 &
            * (273.15 - min(tfrz, ssoil%tggsn(:,1)))))
       ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dt * 9.806 * ssoil%ssdn(:,2) &
            * (.07 * ssoil%ssdn(:,1) + .5 * ssoil%smass(:,2) ) &
            / (3.e7 * exp (.021 * ssoil%ssdn(:,2) + .081 &
            * (273.15 - min(tfrz, ssoil%tggsn(:,2)))))
       ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dt * 9.806 * ssoil%ssdn(:,3) &
            * (.07 * ssoil%ssdn(:,1) + ssoil%smass(:,2) + .5 * ssoil%smass(:,3) ) &
            / (3.e7 * exp (.021 * ssoil%ssdn(:,3) + .081 &
            * (273.15 - min(tfrz, ssoil%tggsn(:,3)))))
       tr1 = ssoil%snowd - ssoil%osnowd
       WHERE (tr1 >= 0.)
          ssoil%ssdn(:,1) = max ( (ssoil%smass(:,1) + tr1) / (ssoil%smass(:,1) &
               / ssoil%ssdn(:,1) + tr1 / 140.), 140.)
          osm = ssoil%smass(:,1)
          ssoil%smass(:,1) = .07 * ssoil%ssdn(:,1)
          ssoil%sdepth(:,1) = .07
          excm = osm + tr1 - ssoil%smass(:,1)
          excd = excm / ssoil%ssdn(:,1)
          osm = ssoil%smass(:,2)
          ssoil%smass(:,2) = max (.01, ssoil%smass(:,2) + .4 * excm)
          ssoil%ssdn(:,2) = max (140., min (500., ssoil%smass(:,2) / &
               (osm / ssoil%ssdn(:,2) + .4 * excd) ) )
          ssoil%sdepth(:,2) = max (.02, ssoil%smass(:,2) / ssoil%ssdn(:,2) )
          osm = ssoil%smass(:,3)
          ssoil%smass(:,3) = max (.01, ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2) )
          ssoil%sdepth(:,3) = max (.02, osm / ssoil%ssdn(:,3) &
               + .6 * excm / ssoil%ssdn(:,2) )
          ssoil%ssdn(:,3) = max (140., min (500., ssoil%smass(:,3) / ssoil%sdepth(:,3) ) )
          WHERE (ssoil%ssdn(:,3) < ssoil%ssdn(:,2) )
             ssoil%ssdn(:,3) = ssoil%ssdn(:,2)
             ssoil%sdepth(:,3) = max (.02, ssoil%smass(:,3) / ssoil%ssdn(:,3) )
          END WHERE
       ELSEWHERE ! snow melting
          ssoil%sdepth(:,1) = .07
          sm1 = max (.01, ssoil%smass(:,1) ) !current mass of the 1st layer
          ! after snow melt
          ! current depth of the 1st layer after density update
          sd1 = max (.005, ssoil%smass(:,1) / ssoil%ssdn(:,1) )
          excd = .07 - sd1 ! required extra depth
          !
          ! add mass from the layer below
          ssoil%smass(:,1) = max (140. * .07, &
               sd1 * ssoil%ssdn(:,1) + excd * ssoil%ssdn(:,2) )
          ssoil%ssdn(:,1) = ssoil%smass(:,1) / .07 ! new density
          excm = ssoil%smass(:,1) - sm1 ! mass differnce
          ! substract only pr
          ssoil%smass(:,2) = max (.01, ssoil%smass(:,2) - min(ssoil%smass(:,2) &
               / (ssoil%smass(:,3) + ssoil%smass(:,2) ), .9) * excm)
          ! fraction of excm from 2nd layer
          ssoil%sdepth(:,2) = max (.02, ssoil%smass(:,2) / ssoil%ssdn(:,2) )
          ! adjust temp.
          ssoil%tggsn(:,1) = ssoil%tggsn(:,1) * sm1 / ssoil%smass(:,1) + &
               (1. - sm1 / ssoil%smass(:,1) ) * ssoil%tggsn(:,2)
          ssoil%smass(:,3) = max (.01, ssoil%snowd - ssoil%smass(:,1) &
               - ssoil%smass(:,2) )
          ssoil%sdepth(:,3) = max (.02, ssoil%smass(:,3) / ssoil%ssdn(:,3) )
          WHERE (ssoil%smass(:,3) < ssoil%smass(:,2) )
             ssoil%smass(:,2) = .45 * (ssoil%snowd - ssoil%smass(:,1) )
             ssoil%sdepth(:,2) = max (.02, ssoil%smass(:,2) / ssoil%ssdn(:,2) )
             ssoil%smass(:,3) = max (.01, ssoil%snowd - ssoil%smass(:,1) &
                  - ssoil%smass(:,2) )
             ssoil%sdepth(:,3) = max (.02, ssoil%smass(:,3) / ssoil%ssdn(:,3) )
          END WHERE
       END WHERE
       ssoil%isflag = 1
       ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%sdepth(:,1) + ssoil%ssdn(:,2) &
            * ssoil%sdepth(:,2) + ssoil%ssdn(:,3) * ssoil%sdepth(:,3) ) &
            / (ssoil%sdepth(:,1) + ssoil%sdepth(:,2) + ssoil%sdepth(:,3))
    END WHERE
    ! Diagnostic block:
    IF (ntest > 0) THEN
       PRINT *,'in soilsnowv before stempv,  ktau= ',ktau
       PRINT *,'ga,dt,ssdn ',canopy%ga(idjd),dt,(ssoil%ssdn(idjd,k),k=1,3)
       PRINT *,'osnowd,snowd,isflag', ssoil%osnowd(idjd),ssoil%snowd(idjd),ssoil%isflag(idjd)
       PRINT *,'tgg ',(ssoil%tgg(idjd,k),k=1,ms)
       PRINT *,'tggsn ',(ssoil%tggsn(idjd,k),k=1,3)
       PRINT *,'wb ',(ssoil%wb(idjd,k),k=1,ms)
       PRINT *,'wbice ',(ssoil%wbice(idjd,k),k=1,ms)
       PRINT *,'wblf ',(ssoil%wblf(idjd,k),k=1,ms)
       PRINT *,'wbfice ',(ssoil%wbfice(idjd,k),k=1,ms)
    END IF
    IF (ntest == 1) THEN
       PRINT *,'in soilsnow printing wbfice_max'
       PRINT *,'sdepth c2 ',(ssoil%sdepth(idjd,k),k=1,3)
    END IF
    CALL stempv(dt, canopy, ssoil, soil)
    DO k = 1, ms
       WHERE (ssoil%tgg(:,k) < tfrz .and. .99 * ssoil%wb(:,k) - ssoil%wbice(:,k) > .001)
          sicefreeze = min (max (0., (.99 * ssoil%wb(:,k) - ssoil%wbice(:,k))) &
               * soil%zse(k) * 1000., (tfrz - ssoil%tgg(:,k) ) * ssoil%gammzz(:,k) / hlf)
          ssoil%wbice(:,k) = min (ssoil%wbice(:,k) + sicefreeze / (soil%zse(k) &
               * 1000.), .99 * ssoil%wb(:,k) )
          ssoil%gammzz(:,k) = max ( (1. - soil%ssat ) * soil%css * soil%rhosoil &
               + (ssoil%wb(:,k) - ssoil%wbice(:,k) ) * cswat * rhowat &
               + ssoil%wbice(:,k) * csice * rhowat * .9, soil%css * soil%rhosoil) &
               * soil%zse(k)
          WHERE (k == 1 .and. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) + sicefreeze * hlf / ssoil%gammzz(:,k)
       ELSEWHERE (ssoil%tgg(:,k) > tfrz .and. ssoil%wbice(:,k) > 0.)
          sicemelt = min (ssoil%wbice(:,k) * soil%zse(k) * 1000., &
               (ssoil%tgg(:,k) - tfrz) * ssoil%gammzz(:,k) / hlf)
          ssoil%wbice(:,k) = max (0., ssoil%wbice(:,k) - sicemelt / (soil%zse(k) * 1000.) )
          ssoil%gammzz(:,k) = max ( (1. - soil%ssat ) * soil%css * soil%rhosoil + &
               (ssoil%wb(:,k) - ssoil%wbice(:,k) ) * cswat * rhowat + ssoil%wbice(:,k) &
               * csice * rhowat * .9, soil%css * soil%rhosoil ) * soil%zse(k)
          WHERE (k == 1 .and. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) - sicemelt * hlf / ssoil%gammzz(:,k)
       END WHERE
    END DO
    DO k = 1,ms
       ! Removing transpiration from soil:
       WHERE (canopy%fevc > 0.) ! convert to mm/dt 
          ! Calculate the (perhaps moisture/ice limited) amount which can be removed:
          evapfbl(:,k) = (min (canopy%fevc * dt / hl * soil%froot(:,k), &
               max (0., min(ssoil%wb(:,k) - soil%swilt,ssoil%wb(:,k)-ssoil%wbice(:,k))) &
               * soil%zse(k) * 1000.)) / (soil%zse(k) * 1000.)
          ! Remove this amount from  the soil:
          ssoil%wb(:,k) = ssoil%wb(:,k) -   evapfbl(:,k)
       END WHERE
    END DO
    ! Adjust fevc
    WHERE (canopy%fevc > 0.) ! convert to mm/dt                                           
       canopy%fevc = (evapfbl(:,1)*soil%zse(1)+evapfbl(:,2)*soil%zse(2) &
            +evapfbl(:,3)*soil%zse(3)+evapfbl(:,4)*soil%zse(4)+evapfbl(:,5) &
            *soil%zse(5)+evapfbl(:,6)*soil%zse(6))*1000.*hl/dt
       
    END WHERE

    CALL surfbv (dt, ktau, canopy, met, ssoil, soil)

    ! Diagnostic block:
    IF (ntest > 0) THEN
       PRINT *,'after surfbv,isflag ',ssoil%isflag(idjd)
       PRINT *,'tgg ',(ssoil%tgg(idjd,k),k=1,ms)
       PRINT *,'wb ',(ssoil%wb(idjd,k),k=1,ms)
       PRINT *,'wblf ',(ssoil%wblf(idjd,k),k=1,ms)
       PRINT *,'wbfice ',(ssoil%wbfice(idjd,k),k=1,ms)
    END IF
    ! Set weighted soil/snow surface temperature:
    ssoil%tss=(1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)

  END SUBROUTINE soil_snow

END MODULE soil_snow_module
