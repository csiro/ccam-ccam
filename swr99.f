      module swr99_m
      implicit none
      public :: swr99
      private :: swr99p, swr99_work, pack_index
      contains
      subroutine swr99 ( fsw_out, hsw_out, grdflx_out, ufsw_out, 
     &                   dfsw_out, press_in, press2_in, coszro_in,
     &                   taudar_in, rh2o_in, rrco2, ssolar,
     &                   qo3_in, nclds_in, ktopsw_in, kbtmsw_in,
     &                   cirab_in, cirrf_in, cuvrf_in, camt_in,
     &                   swr_out, cldoff ) ! MJT cable


!     This driver routine just calculates the size and passes it to 
!     swr99p. This is necessary so that swr99p can use automatic rather
!     than allocatable arrays.

      use newmpar_m
      
      include 'hcon.h'
      include 'rdparm.h'

      real, intent(in), dimension(imax,lp1)  :: press_in, press2_in
      real, intent(in), dimension(imax)      :: coszro_in, taudar_in
      real, intent(in), dimension(imax,l)    :: rh2o_in, qo3_in
      integer, intent(in), dimension(imax)   :: nclds_in
      integer, intent(in), dimension(imax,lp1)  :: ktopsw_in, kbtmsw_in
      real, intent(in), dimension(imax,lp1)  :: cirab_in, cirrf_in, 
     &                                          cuvrf_in, camt_in
      real, intent(in) :: rrco2, ssolar

      real, intent(out), dimension(imax,lp1) :: fsw_out,  
     &                                          ufsw_out, dfsw_out
      real, intent(out), dimension(imax,l)   :: hsw_out 
      real, intent(out), dimension(imax)     :: grdflx_out,swr_out ! MJT cable
      
      integer  :: ipts
      logical, intent(in) :: cldoff
      
      if (cldoff) then
        ipts = count ( coszro_in > zero )
      else
        ipts=count(coszro_in>zero.and.nclds_in>0) ! MJT rad
      end if

      call swr99p      ( fsw_out, hsw_out, grdflx_out, ufsw_out, 
     &                   dfsw_out, press_in, press2_in, coszro_in,
     &                   taudar_in, rh2o_in, rrco2, ssolar,
     &                   qo3_in, nclds_in, ktopsw_in, kbtmsw_in,
     &                   cirab_in, cirrf_in, cuvrf_in, camt_in, ipts,
     &                   swr_out, cldoff ) ! MJT cable

      end subroutine swr99

      subroutine swr99p ( fsw_out, hsw_out, grdflx_out, ufsw_out, 
     &                   dfsw_out, press_in, press2_in, coszro_in,
     &                   taudar_in, rh2o_in, rrco2, ssolar,
     &                   qo3_in, nclds_in, ktopsw_in, kbtmsw_in,
     &                   cirab_in, cirrf_in, cuvrf_in, camt_in, ipts,
     &                   swr_out, cldoff ) ! MJT cable


cfpp$ noconcur r
c===>    *********************************************************
c     -sw- radiation code............................
c        inputs:common block radin (except temp);
c        output:fsw,hsw,grdflx(common block swocom).
c===>    *********************************************************
c
c
      use newmpar_m

      include 'hcon.h'
      include 'rdparm.h'

      real, intent(in), dimension(imax,lp1)  :: press_in, press2_in
      real, intent(in), dimension(imax)      :: coszro_in, taudar_in
      real, intent(in), dimension(imax,l)    :: rh2o_in, qo3_in
      integer, intent(in), dimension(imax)   :: nclds_in
      integer, intent(in), dimension(imax,lp1)  :: ktopsw_in, kbtmsw_in
      real, intent(in), dimension(imax,lp1)  :: cirab_in, cirrf_in, 
     &                                          cuvrf_in, camt_in
      real, intent(in) :: rrco2, ssolar
      integer, intent(in) :: ipts

      real, intent(out), dimension(imax,lp1) :: fsw_out,  
     &                                          ufsw_out, dfsw_out
      real, intent(out), dimension(imax,l)   :: hsw_out 
      real, intent(out), dimension(imax)     :: grdflx_out,swr_out ! MJT cable

      real, dimension(:,:),allocatable,save :: fsw_sav,  
     &                                  ufsw_sav, dfsw_sav ! MJT rad
      real, dimension(:,:),allocatable,save :: hsw_sav     ! MJT rad
      real, dimension(:),allocatable,save   :: grdflx_sav,swr_sav ! MJT rad

      real, dimension(ipts,lp1)      :: press, press2
      real,  dimension(ipts)         :: coszro, taudar
      real,  dimension(ipts,l)       :: rh2o, qo3
      integer,  dimension(ipts)      :: nclds
      integer,  dimension(ipts,lp1)  :: ktopsw, kbtmsw
      real,  dimension(ipts,lp1)     :: cirab, cirrf, cuvrf, camt
      real,  dimension(ipts,lp1)     :: fsw, ufsw, dfsw
      real,  dimension(ipts,l)       :: hsw
      real, dimension(ipts)          :: grdflx,swr ! MJT cable

      integer :: k
      logical, intent(in) :: cldoff
      logical, dimension(imax) :: lcoszn
      integer, dimension(ipts) :: pindex
      
      if (.not.allocated(fsw_sav)) then
        allocate(fsw_sav(imax,lp1),ufsw_sav(imax,lp1))
        allocate(dfsw_sav(imax,lp1),hsw_sav(imax,l))
        allocate(grdflx_sav(imax),swr_sav(imax))
      end if

      lcoszn = coszro_in > zero
      if (.not.cldoff) then          ! MJT rad
        lcoszn=lcoszn.and.nclds_in>0 ! MJT rad
        grdflx_out=grdflx_sav        ! MJT rad
        ufsw_out=ufsw_sav            ! MJT rad
        dfsw_out=dfsw_sav            ! MJT rad
        fsw_out=fsw_sav              ! MJT rad
        hsw_out=hsw_sav              ! MJT rad
      end if                         ! MJT rad

!     If there are no sunlit points set everything to zero and return 
!     immediately
      if ( ipts == 0 ) then
       if (cldoff) then ! MJT rad
        grdflx_out = zero
        ufsw_out   = zero
        dfsw_out   = zero
        fsw_out    = zero
        hsw_out    = zero
        grdflx_sav = zero ! MJT rad
        ufsw_sav   = zero ! MJT rad
        dfsw_sav   = zero ! MJT rad
        fsw_sav    = zero ! MJT rad
        hsw_sav    = zero ! MJT rad
       end if ! MJT rad
       return
      endif

      call pack_index ( lcoszn, pindex )
!     pindex contains indices of the sunlit points

      coszro = coszro_in(pindex)
      taudar = taudar_in(pindex)
      nclds  = nclds_in(pindex)
      do k=1,l
       rh2o(:,k) = rh2o_in(pindex,k)
       qo3(:,k)  = qo3_in(pindex,k)
      end do
      do k=1,lp1
       press(:,k)  = press_in(pindex,k)
       press2(:,k) = press2_in(pindex,k)
       ktopsw(:,k) = ktopsw_in(pindex,k)
       kbtmsw(:,k) = kbtmsw_in(pindex,k)
       cirab(:,k)  = cirab_in(pindex,k)
       cirrf(:,k)  = cirrf_in(pindex,k)
       cuvrf(:,k)  = cuvrf_in(pindex,k)
       camt(:,k)   = camt_in(pindex,k)
      end do

      call swr99_work (  fsw, hsw, grdflx, ufsw, dfsw, press, press2,
     &                   coszro, taudar, rh2o, rrco2, ssolar,
     &                   qo3, nclds, ktopsw, kbtmsw, cirab,
     &                   cirrf, cuvrf, camt, ipts,
     &                   swr ) ! MJT cable

!     Unpack results  ( Use my version to get around NEC compiler error )
      if (cldoff) then ! MJT rad
       grdflx_out = 0.0
       swr_out = 0.5 ! MJT cable
       fsw_out = 0.0
       hsw_out = 0.0
       ufsw_out = 0.0
       dfsw_out = 0.0
      end if ! MJT rad
      grdflx_out(pindex) = grdflx
      swr_out(pindex) = swr ! MJT cable
      do k=1,l
       hsw_out(pindex,k)  = hsw(:,k)
      end do
      do k=1,lp1
       fsw_out(pindex,k)  = fsw(:,k)
       ufsw_out(pindex,k) = ufsw(:,k)
       dfsw_out(pindex,k) = dfsw(:,k)
      end do
      
      if (cldoff) then             ! MJT rad
        grdflx_sav=grdflx_out        ! MJT rad
        ufsw_sav=ufsw_out            ! MJT rad
        dfsw_sav=dfsw_out            ! MJT rad
        fsw_sav=fsw_out              ! MJT rad
        hsw_sav=hsw_out              ! MJT rad
      end if                         ! MJT rad

      end subroutine swr99p

!-------------------------------------------
      subroutine pack_index ( mask, pindex ) 
       logical, dimension(:), intent(in)  :: mask
       integer, dimension(:), intent(out) :: pindex
       integer :: i, j
!      integer, dimension(size(mask)) :: ival

       if ( count(mask) /= size(pindex) ) then
           print*, " Size error in unpack_index ", 
     &              count(mask), size(pindex)
            stop
         end if

!      ival = (/ (i,i=1,size(mask)) /)
!      pindex = pack ( ival, mask )

!        Do this explicitly for better vectorisation on bragg.
         j = 0
         do i=1,size(mask)
            if ( mask(i) ) then
               j = j + 1
               pindex(j) = i
            end if
         end do

      end subroutine pack_index
!--------------------------------------------
      subroutine swr99_work(fsw, hsw, grdflx, ufsw, dfsw, press,press2,
     &                   coszro, taudar, rh2o, rrco2, ssolar,
     &                   qo3, nclds, ktopsw, kbtmsw, cirab,
     &                   cirrf, cuvrf, camt, ipts,
     &                   swr ) ! MJT cable

      use newmpar_m
      use parm_m
      
      include 'hcon.h'
      include 'rdparm.h'

      integer, intent(in) :: ipts
      real, intent(in), dimension(ipts,lp1)  :: press, press2
      real, intent(in), dimension(ipts)      :: coszro, taudar
      real, intent(in), dimension(ipts,l)    :: rh2o, qo3
      integer, intent(in), dimension(ipts)   :: nclds
      integer, intent(in), dimension(ipts,lp1)  :: ktopsw, kbtmsw
      real, intent(in), dimension(ipts,lp1)  :: cirab, cirrf, 
     &                                          cuvrf, camt
      real, intent(in) :: rrco2, ssolar

      real, intent(out), dimension(ipts,lp1) :: fsw, ufsw, dfsw
      real, intent(out), dimension(ipts,l)   :: hsw
      real, intent(out), dimension(ipts)     :: grdflx,swr ! MJT cable

      real dfn(ipts,lp1),ufn(ipts,lp1),
     &  ttd(ipts,lp1),ttu(ipts,lp1),ttdb1(ipts,lp1),ttub1(ipts,lp1),
     &  tco2(ipts,llp2),ud(ipts,lp1),ur(ipts,lp1),
     &  pp(ipts,lp1),dp(ipts,l),pptop(ipts,lp1),dpcld(ipts,lp1),
     &  dfntop(ipts),ct(ipts,lp1),
     &  tucl1(ipts,lp1),tucl1i(ipts,lp1),tdcl1i(ipts,lp1),
     &  ufntrn(ipts,lp1),dfntrn(ipts,lp1),
     &  temp1(ipts),
     &  refl(ipts),secz(ipts),rray(ipts)

      double precision alfa(ipts,lp1),alfau(ipts,lp1),uo3(ipts,llp2),
     &                 ff(lp1,ipts),ffco2(lp1,ipts),ffo3(lp1,ipts),
     &                 pr2(ipts,l),cr(ipts,lp1),
     &  tclu(ipts,lp1),tcld(ipts,lp1),
     &  tdcl1(ipts,lp1),tdcl2(ipts,lp1),
     &  ufnclu(ipts,lp1),dfnclu(ipts,lp1)
      double precision alfa1(lp1),alfau1(lp1) ! MJT CHANGE - mr
      double precision lrd,ltdc               ! MJT CHANGE - mr
      real ctemp                              ! MJT CHANGE - mr
      logical ctest(lp1)                      ! MJT CHANGE - mr
      integer pos(1),k1,k2                    ! MJT CHANGE - mr

      real, dimension(ipts,l)    :: du, duco2, duo3
      real, dimension(ipts,llp2) :: uco2, ao3

      integer :: i, k, kk, nx, kclds, ip, jtop, j2min, j1max
      integer :: j1min, j3max
      real :: htemp, tempf, tempg

!     From GFDL set up
!---specification of data statements:
!         abcff=absorption coefficients for bands in k-distri-
!     bution. originally given by lacis and hansen, revised by
!     ramaswamy
!         pwts=corresponding weights assigned to bands in the
!     k-distribution
!         reflo3,rrayav= reflection coefficients given by
!     lacis and hansen to account for effects of rayleigh scattering
!     in the visible frequencies (band 1)
!         cfco2,cfo3=conversion factors from gm/cm**2 to cm-atm(stp)
!
!---the following are the coefficients for the 12-band shortwave
!   radiation code, specified by ramaswamy.
      real, parameter, dimension(nb) :: abcff = 
     &  (/ 4.0e-5, 4.0e-5, 0.002, 0.035, 0.377, 1.95, 9.40, 44.6, 190.,
     &     989., 2706., 39011. /)
      real, parameter, dimension(nb) :: pwts = 
     &  (/ 0.5000, 0.121416, 0.0698, 0.1558, 0.0631, 0.0362, 0.0243,
     &     0.0158, 0.0087, 0.001467, 0.002342, 0.001075 /)

      real, parameter :: cfco2 = 508.96, cfo3 = 466.64
      real, parameter :: reflo3 = 1.9
      real, parameter :: rrayav = 0.144
      !
!
!---find the maximum number of clouds in the latitude row
      kclds = maxval ( nclds )
!
!   calculate secant of zenith angle (secz),flux pressures(pp),layer
!   width (dp) and pressure scaling factor (pr2).

      secz = h35e1/sqrt(h1224e3*coszro*coszro+one)
!     secz = 1.0/coszro
      temp1 = 1.0/press(:,lp1)

      do k=1,lp1
       pp(:,k) = press2(:,k)
      end do

      do k=1,l
       dp(:,k) = press2(:,k+1) - press2(:,k)
       pr2(:,k) = haf * ( press2(:,k) + press2(:,k+1) )
      end do
      do k=1,l
       pr2(:,k) = pr2(:,k)*temp1
      end do

!---set up angular factor to be multiplied to the optical factor.
!   above the highest cloud, this is the secant of the zenith
!   angle (modified slightly for refractive effects) (=secz).
!   below the highest cloud, this is diffctr (o3difctr for ozone,
!   in accordance with lacis-hansen parameterization).this factor
!   is used regardless of cloud amount-and for direct and diffuse
!   radiation (this is not a 2 1/2 stream model!)

      ff    = diffctr
      ffco2 = diffctr
      ffo3  = o3difctr

      do ip=1,ipts
       jtop = ktopsw(ip,nclds(ip)+1)
       do k=1,jtop
          ffo3 (k,ip) = secz(ip)
          ffco2(k,ip) = secz(ip)
          ff   (k,ip) = secz(ip)
       end do
      end do

!     calculate pressure-weighted optical paths for each layer
!     in units of cm-atm. pressure weighting is using pr2.
!     du= value for h2o;duco2 for co2;duo3 for o3.

      duo3  = (ginv*cfo3)*qo3*dp
      duco2 = (rrco2*ginv*cfco2)*dp*pr2
      du    = ginv*rh2o*dp*pr2

!     obtain the optical path from the top of the atmosphere to the
!     flux pressure. angular factors are now included. ud=downward
!     path for h2o,with ur the upward path for h2o. corresponding
!     quantities for co2,o3 are udco2/urco2 and udo3/uro3.

      ud(:,1)   = zero
      uco2(:,1) = zero
      uo3(:,1)  = zero

      do k=2,lp1
       ud(:,k)   = ud(:,k-1)   + du(:,k-1)*ff(k,:)
       uco2(:,k) = uco2(:,k-1) + duco2(:,k-1)*ffco2(k,:)
       uo3(:,k)  = uo3(:,k-1)  + duo3(:,k-1)*ffo3(k,:)
      end do
      ur(:,lp1)    = ud(:,lp1)
      uco2(:,llp2) = uco2(:,lp1)
      uo3(:,llp2)  = uo3(:,lp1)

      do k=l,1,-1
       ur(:,k) = ur(:,k+1) + du(:,k)*diffctr
       uco2(:,lp1+k) = uco2(:,lp1+k+1) + duco2(:,k)*diffctr
       uo3(:,lp1+k) = uo3(:,lp1+k+1) + duo3(:,k)*reflo3
      end do

!---for the skyhi model only, obtain the oxygen optical path,using
!   the o3 angular integration.because of the equivalencing,this
!   (and all other) optical path calculations must be complete before
!   transmissivities are computed.
!     do 321 k=2,ko2
!     do 321 i=1,ipts
!     uo2(ib(k)+i)=pp(k,i)*ffo3(k,i)
!321  continue
!---end skyhi model only
!
!     calculate co2 absorptions . they will be used in near infrared
!     bands.since the absorption amount is given (in the formula used
!     below, derived from sasamori) in terms of the total solar flux,
!     and the absorption is only included in the near ir (50 percent
!     of the solar spectrum), the absorptions are multiplied by 2.
!       since code actually requires transmissions, these are the
!     values actually stored in tco2.

      tco2(:,2:) = 
     &       1.0-two*(h235m3*exp(hp26*log(uco2(:,2:)+h129m2))-h75826m4)


!     now calculate ozone absorptions. these will be used in
!     the visible band.just as in the co2 case, since this band is
!     50 percent of the solar spectrum,the absorptions are multiplied
!     by 2.
      htemp=h1036e2*h1036e2*h1036e2
      ao3(:,2:)=two*uo3(:,2:)*
     &    (h1p082*exp(hmp805*log(one+h1386e2*uo3(:,2:)))+
     &     h658m2/(one+htemp*uo3(:,2:)*uo3(:,2:)*uo3(:,2:))+
     &     h2118m2/(one+uo3(:,2:)*(h42m2+h323m4*uo3(:,2:))))

!---for skyhi model only
!     calculate o2 absorptions (in visible band)
!     formula is: abs=1.02e-5*uo2(k)**.3795 for uo2<673.9057
!                 abs=4.51e-6*uo2(k)**.5048 for uo2>673.9057
!     the absorption is constant below 35 km (skyhi level 12)
!     do 844 i=ipts1,ipts*ko2
!     if ((uo2(i)-h67390e2).le.zero) then
!       ao2(i)=1.02e-5*uo2(i)**.3795
!     else
!       ao2(i)=4.51e-6*uo2(i)**.5048
!     endif
!844  continue
!     do 624 k=ko21,llp1
!     do 624 i=1,ipts
!     ao2(ib(k)+i)=ao2(ib(ko2)+i)
!624  continue
!---add o2 absorption to o3 absorption (needed only for skyhi).
!     do 33 i=ipts1,ipllp1
!     ao3(i)=ao3(i)+ao2(i)
!33   continue
!---end skyhi model only
!

!---execute the lacis-hansen reflectivity parameterization
!---for the visible band
      rray = hp219/(one+hp816*coszro)
      refl = rray + (one-rray)*(one-rrayav)*cuvrf(:,1)/
     &              (one-cuvrf(:,1)*rrayav)

!   MJT CHANGE - mr (move this outside nx loop)
!---obtain the pressure at the top,bottom and the thickness of
!   "thick" clouds (those at least 2 layers thick). this is used
!   later is obtaining fluxes inside the thick clouds, for all
!   frequency bands.
      do kk=1,kclds
        do i=1,ipts
          if ( (kbtmsw(i,kk+1)-1) > ktopsw(i,kk+1) ) then
             pptop(i,kk) = pp(i,ktopsw(i,kk+1))
             dpcld(i,kk) = one/(pptop(i,kk)-pp(i,kbtmsw(i,kk+1)))
          endif
        end do
      end do
!
      dfsw = 0.0
      ufsw = 0.0
!     start frequency loop (on nx) here
      BAND: do nx=1,nb

!     calculate entering flux at the top for each band(in cgs units)
      dfntop = ssolar*h69766e5*coszro*taudar*pwts(nx)
!
!
      if ( nx == 1 ) then
!        Band 1 (visible) includes O3, O2 and (negligible) H2O absorption
        do k=1,l
          ttdb1(:,k+1) = exp(hm1ez*min(fifty,abcff(1)*ud(:,k+1)))
          ttub1(:,k) = exp(hm1ez*min(fifty,abcff(1)*ur(:,k)))
          ttd(:,k+1) = ttdb1(:,k+1)*(1.-ao3(:,k+1))
          ttu(:,k) = ttub1(:,k)*(1.-ao3(:,lp1+k))
        end do
      else if ( nx == 2 ) then
!       The water vapor transmission function for band 2 is equal to
!       that of band 1 (saved as ttdb1,ttub1)
       do k=1,l
          ttd(:,k+1) = ttdb1(:,k+1)*tco2(:,k+1)
          ttu(:,k)   = ttub1(:,k)*tco2(:,lp1+k)
       end do
      else
!       Calculate water vapor transmission functions for near infrared
!       bands. Include CO2 transmission (tdco2/tuco2), which
!       is the same for all infrared bands.
       do k=1,l
          ttd(:,k+1) = exp(hm1ez*min(fifty,abcff(nx)*ud(:,k+1)))*
     &                   tco2(:,k+1)
            ttu(:,k)   = exp(hm1ez*min(fifty,abcff(nx)*ur(:,k)))*
     &                   tco2(:,lp1+k)
         end do

      end if

!---at this point,include ttd(1),ttu(lp1), noting that ttd(1)=1 for

      ttd(:,1)   = 1.
      ttu(:,lp1) = ttd(:,lp1)
!
!
!***if no clouds at all exist in the row, the calculations are
!   drastically simplified.
      if ( kclds == 0 ) then
       dfn = ttd
       if ( nx == 1 ) then  ! Visible
          temp1 = refl*dfn(:,lp1)/ttu(:,lp1)
       else  ! IR
          temp1 = cirrf(:,1)*dfn(:,lp1)/ttu(:,lp1)
       end if
       do k=1,lp1
          ufn(:,k) = temp1*ttu(:,k)
       end do
      endif
!
!***compute normal case: at least 1 pt has a cloud.
!

!
!---the following calculation is done for nl=1 case only
!
!   MJT CHANGE - mr (moved outside nx loop)
!---obtain the pressure at the top,bottom and the thickness of
!   "thick" clouds (those at least 2 layers thick). this is used
!   later is obtaining fluxes inside the thick clouds, for all
!   frequency bands.
!      do kk=1,kclds
!      do i=1,ipts
!         if ( (kbtmsw(i,kk+1)-1) > ktopsw(i,kk+1) ) then
!            pptop(i,kk) = pp(i,ktopsw(i,kk+1))
!            dpcld(i,kk) = one/(pptop(i,kk)-pp(i,kbtmsw(i,kk+1)))
!         endif
!      end do
!      end do
!
!***for execution of the cloud loop, it is necessary to separate out
!   transmission fctns at the top and bottom of the clouds, for
!   each band n. the required quantities are:
!      ttd(i,ktopsw(i,k),nl)  k runs from 1 to nclds(i)+1:
!      ttd(i,kbtmsw(i,k),nl)  k runs from 1 to nclds(i)+1:
!      ttu(i,ktopsw(i,k),nl)  k runs from 1 to nclds(i)+1:
!      and inverses of the above. the above quantities are stored
!      in tdcl1,tdcl2,tucl1,tdcl1i,tdcl2i,tucli,respectively, as
!      they have multiple use in the pgm.
!
!---for first cloud layer (ground) tdcl1,tucl1 are known:
      tdcl1(:,1) = ttd(:,lp1)
      tucl1(:,1) = ttu(:,lp1)
      tdcl1i(:,1) = one/tdcl1(:,1)
      tucl1i(:,1) = tdcl1i(:,1)
      tdcl2(:,1) = tdcl1(:,1)

!---for remaining ones,use gathers:
      do kk=2,kclds+1
       do i=1,ipts
          tdcl1(i,kk) = ttd(i,ktopsw(i,kk))
          tucl1(i,kk) = ttu(i,ktopsw(i,kk))
          tdcl2(i,kk) = ttd(i,kbtmsw(i,kk))
       end do
      end do
!---compute inverses
      do k=2,kclds+1
         tdcl1i(:,k) = one/tdcl1(:,k)
         tucl1i(:,k) = one/tucl1(:,k)
      end do


      if ( kclds /= 0 ) then
      
      if ( nmr==2 ) then ! MJT CHANGE - mr
      ! modified M/R overlap

      if ( nx == 1 ) then
!---the first cloud is the ground; its properties are given
!   by refl (the transmission (0) is irrelevant for now!).
       cr(:,1) = refl

!***obtain cloud reflection and transmission coefficients for
!   remaining clouds (if any) in the visible band
!---the maximum no of clouds in the row (kclds) is used. this creates
!   extra work (may be removed in a subsequent update).

!   note 100% cloud cover (removed camt term)
       cr(:,2:kclds+1)=cuvrf(:,2:kclds+1)
       ct(:,2:kclds+1)=one-cr(:,2:kclds+1)
  
      else  !  IR, no Rayleigh scattering

!   note 100% cloud cover (removed camt term)
       cr(:,1:kclds+1)=cirrf(:,1:kclds+1)
       ct(:,1:kclds+1)=one-(cirrf(:,1:kclds+1)+cirab(:,1:kclds+1))

      end if

      ! calculate transmission through clouds
      alfa(:,1) = cr(:,1)
      alfau(:,1) = zero
      do k=2,kclds+1
        alfa(:,k) = cr(:,1)*(tdcl1(:,1)*tdcl1i(:,k))**2 ! clear sky
        alfau(:,k) = tdcl1(:,k-1)*tdcl1i(:,k)
      end do
      do i=1,ipts
        k=2 
        do while (k.le.nclds(i)+1)
          kk=k ! k=bottom of cloud block
          do while (kbtmsw(i,kk+1).eq.ktopsw(i,kk).and.kk.lt.lp1
     &              .and.kbtmsw(i,kk+1).gt.1) ! find top of cloud block
            kk=kk+1 ! kk=top of cloud block
          end do
          ctest=.false.
          ctest(k:kk)=.true. ! cloud layers in current path through cloud block
          ctemp=zero
          alfa(i,k:kk)=zero  ! initalise weighted sum at zero
          alfau(i,k:kk)=zero
          alfa1(k-1)=alfa(i,k-1)
          ! loop over possible paths through cloud block (maximum overlap levels within cloud block)
          do while (any(ctest(k:kk)))
            pos=minloc(camt(i,k:kk),ctest(k:kk))
            k1=pos(1)-1+k ! level of smallest cloud fraction from remaning levels
            ctemp=camt(i,k1)-ctemp ! overlap fraction for current path
            lrd=alfa1(k-1)    ! last non-trivial reflection
            ltdc=tdcl1(i,k-1) ! tdcl1 at last non-trivial reflection
            do k2=k,kk
              if (ctest(k2)) then !cloud
                ! transmissivity from top of cloud layer k2 to top of cloud layer k2-1
                tclu(i,k2-1)=ltdc*tdcl1i(i,k2)*ct(i,k2)
                ! transmissivity from bottom of cloud layer k2 to top of cloud layer k2-1
                tcld(i,k2-1)=ltdc/tdcl2(i,k2)
                ! recursion relation for reflection of the current cloud layer
                ! and all clouds below
                alfa1(k2)=tclu(i,k2-1)*tclu(i,k2-1)*lrd
     &            /(1.-tcld(i,k2-1)*tcld(i,k2-1)*lrd*cr(i,k2))
     &            +cr(i,k2)
                ! ratio of downwards flux between level k2 and k2-1 (for calculating dfncld)
                alfau1(k2)=tdcl1(i,k2-1)*tdcl1i(i,k2)*ct(i,k2)
     &            /(1.-tcld(i,k2-1)*tcld(i,k2-1)*lrd*cr(i,k2))
                lrd=alfa1(k2) 
                ltdc=tdcl1(i,k2)
              else !clear sky
                tclu(i,k2-1)=tdcl1(i,k2-1)*tdcl1i(i,k2)
                alfa1(k2)=tclu(i,k2-1)*tclu(i,k2-1)*alfa1(k2-1)
                alfau1(k2)=tclu(i,k2-1)
              end if
            end do
            alfa(i,k:kk)=alfa(i,k:kk)+ctemp*alfa1(k:kk) ! sum of reflections over all paths
            alfau(i,k:kk)=alfau(i,k:kk)+ctemp*alfau1(k:kk)
            ctemp=camt(i,k1)
            ctest(k1)=.false.
          end do
!dir$ ivdep
          alfa(i,k:kk)=alfa(i,k:kk)
     &      +(1.-ctemp)*alfa(i,k-1)*(tdcl1(i,k-1)*tdcl1i(i,k:kk))**2
          alfau(i,k:kk)=alfau(i,k:kk)
     &      +(1.-ctemp)*tdcl1(i,k-1)*tdcl1i(i,k:kk)
          k=kk+1
        end do
        alfa(i,k:kclds+1)=alfa(i,k-1)
     &    *(tdcl1(i,k-1)*tdcl1i(i,k:kclds+1))**2 ! clear sky
        ! alfau is already calculated for the clear sky
      end do

!     calculate ufn at cloud tops and dfn at cloud bottoms
!---note that ufnclu(i,kclds+1) gives the upward flux at the top
!   of the highest real cloud (if nclds(i)=kclds). it gives the flux
!   at the top of the atmosphere if nclds(i) < kclds. in the first
!   case, tdcl1 equals the transmission fctn to the top of the
!   highest cloud, as we want. in the second case, tdcl1=1, so ufnclu
!   equals alfa. this is also correct.
      ufnclu(:,kclds+1) = alfa(:,kclds+1)*tdcl1(:,kclds+1)
      dfnclu(:,kclds+1) = tdcl1(:,kclds+1)

!---this calculation is the reverse of the recursion relation used
!  above
      do kk=kclds,1,-1
       dfnclu(:,kk) = dfnclu(:,kk+1)*alfau(:,kk+1)
       ufnclu(:,kk) = dfnclu(:,kk)*alfa(:,kk)
      end do
      
      else ! usual

      if ( nx == 1 ) then
!---the first cloud is the ground; its properties are given
!   by refl (the transmission (0) is irrelevant for now!).
       cr(:,1) = refl

!***obtain cloud reflection and transmission coefficients for
!   remaining clouds (if any) in the visible band
!---the maximum no of clouds in the row (kclds) is used. this creates
!   extra work (may be removed in a subsequent update).
       do k=2,kclds+1
        cr(:,k) = cuvrf(:,k)*camt(:,k)
        ct(:,k) = one-cr(:,k)
       end do

      else  !  IR, no Rayleigh scattering

       do k=1,kclds+1
        cr(:,k) = cirrf(:,k)*camt(:,k)
        ct(:,k) = one-camt(:,k)*(cirrf(:,k)+cirab(:,k))
       end do

      end if

!---compute the transmissivity from the top of cloud (k+1) to the
!   top of cloud (k). the cloud transmission (ct) is included. this
!   quantity is called tclu (index k). also, obtain the transmissivity
!   from the bottom of cloud (k+1) to the top of cloud (k)(a path
!   entirely outside clouds). this quantity is called tcld (index k).
      do k=1,kclds
         tclu(:,k) = tdcl1(:,k)*tdcl1i(:,k+1)*ct(:,k+1)
         tcld(:,k) = tdcl1(:,k)/tdcl2(:,k+1)
      end do
!
!***the following is the recursion relation for alfa: the reflection
!   coefficient for a system including the cloud in question and the
!   flux coming out of the cloud system including all clouds below
!   the cloud in question.
!---alfau is alfa without the reflection of the cloud in question
      alfa(:,1) = cr(:,1)
      alfau(:,1) = zero
      do kk=2,kclds+1
!---again,excessive calculations-may be changed later!
       alfau(:,kk) = tclu(:,kk-1)*tclu(:,kk-1)*alfa(:,kk-1) / 
     &            (1.-tcld(:,kk-1)*tcld(:,kk-1)*alfa(:,kk-1)*cr(:,kk))
         alfa(:,kk) = alfau(:,kk) + cr(:,kk)
      end do
!
!     calculate ufn at cloud tops and dfn at cloud bottoms
!---note that ufnclu(i,kclds+1) gives the upward flux at the top
!   of the highest real cloud (if nclds(i)=kclds). it gives the flux
!   at the top of the atmosphere if nclds(i) < kclds. in the first
!   case, tdcl1 equals the transmission fctn to the top of the
!   highest cloud, as we want. in the second case, tdcl1=1, so ufnclu
!   equals alfa. this is also correct.
      ufnclu(:,kclds+1) = alfa(:,kclds+1)*tdcl1(:,kclds+1)
      dfnclu(:,kclds+1) = tdcl1(:,kclds+1)

!---this calculation is the reverse of the recursion relation used
!  above
      do kk=kclds,1,-1
       ufnclu(:,kk) = ufnclu(:,kk+1)*alfau(:,kk+1)/
     &                  (alfa(:,kk+1)*tclu(:,kk))
         dfnclu(:,kk) = ufnclu(:,kk)/alfa(:,kk)
      end do
      
      end if
      
!     now obtain dfn and ufn for levels between the clouds
      do k=1,kclds+1
       ufntrn(:,k) = ufnclu(:,k)*tucl1i(:,k)
       dfntrn(:,k) = dfnclu(:,k)*tdcl1i(:,k)
      end do
!---case of kk=1 (from the ground to the bottom of the lowest cloud)
      j2min = minval(kbtmsw(:,2))
      do k=j2min,lp1
        where ( k >= kbtmsw(:,2) )
          ufn(:,k) = ufntrn(:,1)*ttu(:,k)
          dfn(:,k) = dfntrn(:,1)*ttd(:,k)
        end where
      end do
      do kk=2,kclds+1
       j1max = maxval(ktopsw(:,kk))
       j2min = minval(kbtmsw(:,kk+1))
       do k=j2min,j1max
          where ( k >= kbtmsw(:,kk+1) .and. k <= ktopsw(:,kk) .and. 
     &              ktopsw(:,kk) > 1 )
             ufn(:,k) = ufntrn(:,kk)*ttu(:,k)
             dfn(:,k) = dfntrn(:,kk)*ttd(:,k)
          endwhere
       end do
         j3max = maxval(kbtmsw(:,kk))
         j1min = minval(ktopsw(:,kk))
         if ( j3max - j1min > 1 ) then
            do i=1,ipts
               if ( kbtmsw(i,kk) - ktopsw(i,kk) > 1 ) then
                  tempf = (ufnclu(i,kk)-ufn(i,kbtmsw(i,kk)))*
     $                     dpcld(i,kk-1)
                  tempg = (dfnclu(i,kk)-dfn(i,kbtmsw(i,kk)))*
     $                     dpcld(i,kk-1)
                  do k = ktopsw(i,kk)+1, kbtmsw(i,kk)-1
                     ufn(i,k) = ufnclu(i,kk) +
     $                          tempf*(pp(i,k)-pptop(i,kk-1))
                     dfn(i,k) = dfnclu(i,kk) +
     $                          tempg*(pp(i,k)-pptop(i,kk-1))
                  end do
               end if 
            end do
         end if
      end do

      endif  ! kclds /= 0
!
!     scale visible band fluxes by solar flux at the top of the
!     atmosphere (dfntop(i))
!     dfsw/ufsw will be the fluxes, summed over all bands
      do k=1,lp1
       dfsw(:,k) = dfsw(:,k) + dfn(:,k)*dfntop
       ufsw(:,k) = ufsw(:,k) + ufn(:,k)*dfntop
      end do

      !--------------------------------------------------------------
      ! MJT cable
      if (nx.eq.1) swr=dfsw(:,lp1) ! store VIS
      !--------------------------------------------------------------


      end do BAND ! end of frequency loop (over nx)
!
      grdflx = dfsw(:,lp1) - ufsw(:,lp1)
!
!---obtain the net flux and the shortwave heating rate
      fsw = ufsw - dfsw

      do k=1,l
	 hsw(:,k) = fsw(:,k+1) - fsw(:,k)
      end do

      !--------------------------------------------------------------
      ! MJT cable
      swr=swr/dfsw(:,lp1) ! calculate ratio of VIS to TOTAL
      !--------------------------------------------------------------

      return
      end subroutine swr99_work

      end module swr99_m
