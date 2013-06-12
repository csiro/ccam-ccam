      ! CCAM nudging/assimilation routines
      
      ! These routines preturb the regional model with the large scale circulation of the host model.
      ! Currently, relaxiation, far-field and scale-selective filter options are supported for both
      ! the atmosphere and ocean.
      
      ! We support both 1D and 2D versions of the scale-selective filter.  2D is exact, but expensive.
      ! Current tests suggest the 1D is a good approximation of the 2D filter where the grid stretching
      ! is not too large.

      ! nbd/=0     Far-field or relaxation nudging
      ! mbd/=0     Spectral filter (1D and 2D versions, see nud_uv)
      ! nud_uv =1    Nudge winds (=9 for 2D filter)
      ! nud_t  =1    Nudge air temperature
      ! nud_qg =1    Nudge mixing ratio
      ! nud_p  =1    Nudge surface pressure
      ! nud_sst=1    Nudge water temperature (numbers greater than mbd control strength)
      ! nud_sss=1    Nudge salinity
      ! nud_ouv=1    Nudge water currents
      ! nud_sfh=1    Nudge water free surface height
      ! kbotdav      Lowest atmospheric level to nudge
      ! ktopdav      Highest atmospheric level to nudge
      ! ktopmlo      Deepest water level to nudge
      ! kbotmlo      Shallowest water level to nudge
      ! mloalpha     Weight of water nudging strength

      !--------------------------------------------------------------
      ! FAR-FIELD NUDGING AND RELAXIATION ROUTINES
      ! Called for nbd/=0
      subroutine nestin(iaero)
      
      use arrays_m                     ! Atmosphere dyamics prognostic arrays
      use cc_mpi                       ! CC MPI routines
      use davb_m                       ! Far-field nudging (host store)
      use diag_m                       ! Diagnostic routines
      use indices_m                    ! Grid index arrays
      use latlong_m                    ! Lat/lon coordinates
      use mlo                          ! Ocean physics and prognostic arrays
      use pbl_m                        ! Boundary layer arrays
      use soil_m                       ! Soil and surface data
      use soilsnow_m                   ! Soil, snow and surface data
      
      implicit none
      
      include 'newmpar.h'              ! Grid parameters
      include 'dates.h'                ! Date data
      include 'parm.h'                 ! Model configuration
      include 'stime.h'                ! File date data

      integer, dimension(ifull) :: dumm
      integer, save :: num,mtimea,mtimeb
      integer iq,iaero,k,i,wl,ierr
      integer kdate_r,ktime_r,kdhour,kdmin,iabsdate
      real timerm,cona,conb,rduma
      real, save :: rdumg = -999.
      real, dimension(2) :: dumbb
      real, dimension(:,:), allocatable, save :: ta,ua,va,qa
      real, dimension(:,:), allocatable, save :: tb,ub,vb,qb,ocndep
      real, dimension(:), allocatable, save :: psla,pslb,tssa,tssb
      real, dimension(:), allocatable, save :: sicedepb,fraciceb
      real, dimension(:,:,:), allocatable, save :: sssa,sssb
      real, dimension(ifull) :: zsb,duma,timelt
      real, dimension(ifull,wlev,4) :: dumaa
      real, dimension(ifull,ms) :: dumg
      real, dimension(ifull,kl) :: dumv
      real, dimension(ifull,3) :: dums
      character(len=12) dimnam
      data num/0/,mtimea/0/,mtimeb/-1/
      
!     mtimer, mtimeb are in minutes
      if(ktau<100.and.myid==0)then
        write(6,*) 'in nestin ktau,mtimer,mtimea,mtimeb ',
     &                        ktau,mtimer,mtimea,mtimeb
        write(6,*) 'with kdate_s,ktime_s >= ',kdate_s,ktime_s
      end if

      ! Load next host model timestep for nudging
      if(mtimer>mtimeb) then  ! allows for dt<1 minute
      
        ! Intialise nudging
        if (.not.allocated(ta)) then
          ! Allocate host data arrays
          allocate(ta(ifull,kl),ua(ifull,kl),va(ifull,kl),qa(ifull,kl))
          allocate(tb(ifull,kl),ub(ifull,kl),vb(ifull,kl),qb(ifull,kl))
          allocate(psla(ifull),pslb(ifull),tssa(ifull),tssb(ifull))
          allocate(sicedepb(ifull),fraciceb(ifull),ocndep(ifull,2))
          allocate(sssa(ifull,wlev,4),sssb(ifull,wlev,4))

          ! Save host atmospheric data
          if ( myid==0 ) write(6,*)
     &      'set nesting fields to those already read in via indata'
          pslb(:)=psl(:)
          tssb(:)=tss(:)
          sicedepb(:)=sicedep(:)
          fraciceb(:)=fracice(:)
          tb(:,:)=t(1:ifull,:)
          qb(:,:)=qg(1:ifull,:)
          ub(:,:)=u(1:ifull,:)
          vb(:,:)=v(1:ifull,:)

          ! Save host ocean data
          if (nmlo/=0) then
            ocndep=0.
            sssb(:,:,1)=293.16
            sssb(:,:,2)=34.72
            sssb(:,:,3)=0.
            sssb(:,:,4)=0.
            do i=1,4
              do k=1,wlev
                call mloexport(i-1,sssb(:,k,i),k,0)
              end do
            end do
          end if
        
          ! record time of saved data
          mtimeb=mtimer
        endif       ! (.not.allocated(ta))
      
!       transfer mtimeb fields to mtimea and update sice variables
        mtimea=mtimeb
        psla(:)=pslb(:)
        tssa(:)=tssb(:)
        ta(1:ifull,:)=tb(1:ifull,:)
        qa(1:ifull,:)=qb(1:ifull,:)
        ua(1:ifull,:)=ub(1:ifull,:)
        va(1:ifull,:)=vb(1:ifull,:)
        if (nmlo/=0) then
          sssa(:,:,:)=sssb(:,:,:)
        end if

        ! Read sea-ice data from host when not using
        ! AMIP SSTs or Mixed-Layer-Ocean sea-ice      
        if(namip==0.and.nmlo==0)then
!         check whether present ice points should change to/from sice points
          sicedep(:)=sicedepb(:)
          fracice(:)=fraciceb(:)
!         ensure that sice is only over sea
          do iq=1,ifull
            if(fraciceb(iq)>0..and.fracice(iq)==0.)then
!             N.B. if already a sice point, keep present tice (in tggsn)
              tggsn(iq,1)=min(271.2,tssb(iq),tb(iq,1)+.04*6.5) ! for 40 m lev1
            endif  ! (fraciceb(iq)==0..and.fracice(iq))
            if(fracice(iq)<.02)fracice(iq)=0.
            if(land(iq))then
              sicedep(iq)=0.
              fracice(iq)=0.
            else
              if(fracice(iq)>0..and.sicedep(iq)==0.)then
!               assign to 2. in NH and 1. in SH (according to spo)
!               do this in indata, amipdata and nestin because of onthefly
                if(rlatt(iq)>0.)then
                  sicedep(iq)=2.
                else
                  sicedep(iq)=1.
                endif ! (rlatt(iq)>0.)
              elseif(fracice(iq)==0..and.sicedep(iq)>0.)then  ! e.g. from Mk3  
                fracice(iq)=1.
              endif  ! (fracice(iq)>0..and.sicedep(iq)==0.) .. elseif ..
            endif    ! (land(iq))
          enddo     ! iq loop
        endif ! (namip==0)

        ! Read host atmospheric and ocean data for nudging      
        if(abs(io_in)==1)then
          call onthefly(1,kdate_r,ktime_r,
     &                 pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb, 
     &                 dumg,dumg,dumg,duma,dumv,dumv,dumv,dums,dums,
     &                 dums,duma,duma,dumm,iaero,sssb,ocndep)
        else
          write(6,*) 'ERROR: Nudging requires abs(io_in)=1'
          call ccmpi_abort(-1)
        endif   ! (io_in==1)
        tssb(:) = abs(tssb(:))
        if (mydiag) then
          write (6,"('zsb# nestin  ',9f7.1)") diagvals(zsb)
          write (6,"('tssb# nestin ',9f7.1)") diagvals(tssb) 
        end if
   
        ! determine time corrosponding to new host nudging data
        kdhour=ktime_r/100-ktime/100
        kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
        mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     &               +60*kdhour+kdmin

        if ( myid == 0 ) then
          write(6,*) 'nesting file has: kdate_r,ktime_r,kdhour,kdmin ',
     &                             kdate_r,ktime_r,kdhour,kdmin
          write(6,*) 'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
          write(6,*) 'giving mtimeb = ',mtimeb
!         print additional information
          write(6,*) ' kdate ',kdate,' ktime ',ktime
          write(6,*) 'timeg,mtimer,mtimea,mtimeb: ',
     &                timeg,mtimer,mtimea,mtimeb
          write(6,*) 'ds ',ds
        end if

!       ensure qb big enough, but not too big in top levels (from Sept '04)
        qb(1:ifull,:)=max(qb(1:ifull,:),0.)

!       following is useful if troublesome data is read in
        if(mod(ktau,nmaxpr)==0.or.ktau==2.or.diag)then
          if ( myid == 0 ) then
            write(6,*) 'following max/min values printed from nestin'
          end if
          call maxmin(ub,'ub',ktau,1.,kl)
          call maxmin(vb,'vb',ktau,1.,kl)
          call maxmin(tb,'tb',ktau,1.,kl)
          call maxmin(qb,'qb',ktau,1.e3,kl)
          if ( myid == 0 ) then
            write(6,*) 'following are really psl not ps'
          end if
          call maxmin(pslb,'ps',ktau,100.,1)
        endif

!       in these cases redefine pslb, tb and (effectively) zsb using zs
!       this keeps fine-mesh land mask & zs
!       presently simplest to do whole pslb, tb (& qb) arrays
        if(nmaxpr==1.and.mydiag)then
          write(6,*) 'zs (idjd) :',zs(idjd)
          write(6,*) 'zsb (idjd) :',zsb(idjd)
          write (6,"('100*psl.wesn ',2p5f8.3)") psl(idjd),psl(iw(idjd)),
     &              psl(ie(idjd)),psl(is(idjd)),psl(in(idjd))
          write (6,"('ps.wesn ',-2p5f9.3)") ps(idjd),
     &           ps(iw(idjd)),ps(ie(idjd)),ps(is(idjd)),ps(in(idjd))
          write(6,*) 'pslb in(idjd) :',pslb(idjd)
          write(6,*) 'now call retopo from nestin'
        endif
        call retopo(pslb,zsb,zs(1:ifull),tb,qb)
        if(nmaxpr==1.and.mydiag)then
          write (6,"('100*pslb.wesn ',2p5f8.3)") pslb(idjd),
     &       pslb(iw(idjd)),pslb(ie(idjd)),pslb(is(idjd)),pslb(in(idjd))
          write(6,*) 'pslb out(idjd) :',pslb(idjd)
          write(6,*) 'after pslb print; num= ',num
        endif

        ! display diagnostics      
        if(num==0)then
          num=1
          call printa('zs  ',zs        ,ktau,0  ,ia,ib,ja,jb,0.,.01)
          call printa('zsb ',zsb       ,ktau,0  ,ia,ib,ja,jb,0.,.01)
          call printa('psl ',psl       ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
          call printa('pslb',pslb      ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
          call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
          call printa('tb  ',tb,ktau,nlv,ia,ib,ja,jb,200.,1.)
          call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('ub  ',ub,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vb  ',vb,ktau,nlv,ia,ib,ja,jb,0.,1.)
          return
        endif   !  num==0
      
      end if ! (mtimer>mtimeb)

!     now use tt, uu, vv arrays for time interpolated values
      timerm=ktau*dt/60.   ! real value in minutes (in case dt < 60 seconds)
      cona=(mtimeb-timerm)/real(mtimeb-mtimea)
      conb=(timerm-mtimea)/real(mtimeb-mtimea)
      psls(:)=cona*psla(:)+conb*pslb(:)
      tt (:,:)=cona*ta(:,:)+conb*tb(:,:)
      qgg(:,:)=cona*qa(:,:)+conb*qb(:,:)
      uu (:,:)=cona*ua(:,:)+conb*ub(:,:)
      vv (:,:)=cona*va(:,:)+conb*vb(:,:)

!     calculate time interpolated tss 
      if(namip==0)then     ! namip SSTs/sea-ice take precedence
        if (nmlo==0) then
          ! SSTs read from host model
          where (.not.land)
            tss=cona*tssa+conb*tssb
            tgg(:,1)=tss
          end where  ! (.not.land)
        else
          if (nud_sst/=0.or.nud_sss/=0.or.nud_ouv/=0.or.
     &        nud_sfh/=0) then
            ! nudge mlo
            dumaa=cona*sssa+conb*sssb
            wl=wlev
            if (rdumg<-1.) then
              rduma=maxval(ocndep(:,1)) ! check if 3D data exists
              dumbb(1)=rduma
              call ccmpi_allreduce(dumbb(1:1),dumbb(2:2),"max",
     &                             comm_world)
              rdumg=dumbb(2)
            end if
            if (rdumg<0.5) wl=1
            if (wl==1) then ! switch to 2D if 3D data is missing
              !call mloexpmelt(timelt)
              call mloexport(0,timelt,1,0)
              dumaa(:,1,1)=cona*tssa+conb*tssb
              where (fraciceb>0.)
                dumaa(:,1,1)=timelt
              end where
            end if
            call mlonudge(dumaa(:,:,1),dumaa(:,:,2),
     &                    dumaa(:,:,3:4),ocndep(:,2),wl)
          end if
        endif ! nmlo==0 ..else..
      endif   ! namip==0
      
      return
      end subroutine nestin


      !--------------------------------------------------------------
      ! SCALE SELECTIVE FILTER ASSIMILATION
      ! Called for mbd/=0
      subroutine nestinb(iaero)

      use arrays_m                     ! Atmosphere dyamics prognostic arrays
      use cc_mpi                       ! CC MPI routines
      use diag_m                       ! Diagnostic routines
      use indices_m                    ! Grid index arrays
      use latlong_m                    ! Lat/lon coordinates
      use mlo                          ! Ocean physics and prognostic arrays
      use pbl_m                        ! Boundary layer arrays
      use soil_m                       ! Soil and surface data
      use soilsnow_m                   ! Soil, snow and surface data
 
      implicit none
 
      include 'newmpar.h'              ! Grid parameters
      include 'dates.h'                ! Date data
      include 'parm.h'                 ! Model configuration
      include 'stime.h'                ! File date data
 
      integer, dimension(ifull) :: dumm
      integer, save :: mtimeb = -1
      integer kdate_r,ktime_r,iaero,wl,ierr
      integer iabsdate,iq,k,kdhour,kdmin
      real ds_r,timeg_b,rduma
      real, save :: rdumg = -999.
      real, dimension(2) :: dumbb
      real, dimension(:,:), allocatable, save :: tb,ub,vb,qb,ocndep
      real, dimension(:), allocatable, save :: pslb,tssb,fraciceb
      real, dimension(:), allocatable, save :: sicedepb
      real, dimension(:,:,:), allocatable, save :: sssb
      real, dimension(ifull) :: zsb,pslc,duma,timelt
      real, dimension(ifull,ms) :: dumg
      real, dimension(ifull,kl) :: dumv
      real, dimension(ifull,kl) :: uc,vc,tc,qc
      real, dimension(ifull,3) :: dums
 
      ! allocate arrays on first call     
      if (.not.allocated(tb)) then
        allocate(tb(ifull,kl),ub(ifull,kl),vb(ifull,kl),qb(ifull,kl))
        allocate(pslb(ifull),tssb(ifull),fraciceb(ifull))
        allocate(sicedepb(ifull),ocndep(ifull,2))
        allocate(sssb(ifull,wlev,4))
      end if

!     mtimer, mtimeb are in minutes
      if(ktau<100.and.myid==0)then
        write(6,*) 'in nestinb ktau,mtimer,mtimeb,io_in ',
     &                      ktau,mtimer,mtimeb,io_in
        write(6,*) 'with kdate_s,ktime_s >= ',kdate_s,ktime_s
      end if

      ! Load new host data to be ready for next call to filter
      if (mtimer>mtimeb) then

!       following (till end of subr) reads in next bunch of data in readiness
!       read tb etc  - for globpea, straight into tb etc
        if (abs(io_in)==1) then
          call onthefly(1,kdate_r,ktime_r,
     &                 pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb, 
     &                 dumg,dumg,dumg,duma,dumv,dumv,dumv,dums,dums,
     &                 dums,duma,duma,dumm,iaero,sssb,ocndep)
        else
          write(6,*) 'ERROR: Scale-selective filter requires ',
     &               'abs(io_in)=1'
          call ccmpi_abort(-1)
        endif   ! (abs(io_in)==1)
        tssb(:) = abs(tssb(:))  ! moved here Mar '03
        if (mydiag) then
          write (6,"('zsb# nestinb  ',9f7.1)") diagvals(zsb)
          write (6,"('tssb# nestinb ',9f7.1)") diagvals(tssb) 
        end if

        ! calculate time for next filter call   
        kdhour=ktime_r/100-ktime/100   ! integer hour diff from Oct '05
        kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
        mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     &                +60*kdhour+kdmin
        if ( myid == 0 ) then
          write(6,*) 'nestinb file has: kdate_r,ktime_r,kdhour,kdmin ',
     &                                  kdate_r,ktime_r,kdhour,kdmin
          write(6,*) 'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
          write(6,*) 'giving mtimeb = ',mtimeb
!         print additional information
          write(6,*) ' kdate ',kdate,' ktime ',ktime
          write(6,*) 'timeg,mtimer,mtimeb: ',
     &                timeg,mtimer,mtimeb
          write(6,*) 'ds ',ds
        end if

        if(mod(ktau,nmaxpr)==0.or.ktau==2.or.diag)then
!         following is useful if troublesome data is read in
          if ( myid == 0 ) then
            write(6,*) 'following max/min values printed from nestinb'
          end if
          call maxmin(ub,'ub',ktau,1.,kl)
          call maxmin(vb,'vb',ktau,1.,kl)
          call maxmin(tb,'tb',ktau,1.,kl)
          call maxmin(qb,'qb',ktau,1.e3,kl)
          if ( myid == 0 ) then
            write(6,*) 
     &      'following in nestinb after read pslb are psl not ps'
          end if
          call maxmin(pslb,'pB',ktau,100.,1)
        endif

!       in these cases redefine pslb, tb and (effectively) zsb using zs
!       this keeps fine-mesh land mask & zs
!       presently simplest to do whole pslb, tb (& qb) arrays
        if(nmaxpr==1.and.mydiag)then
          write(6,*) 'zs (idjd) :',zs(idjd)
          write(6,*) 'zsb (idjd) :',zsb(idjd)
          write(6,*) 'pslb in(idjd) :',pslb(idjd)
          write(6,*) 
     &     'call retopo from nestin; psl# prints refer to pslb'
        endif
        call retopo(pslb,zsb,zs(1:ifull),tb,qb)

        ! Removed by MJT for conservation
!       ensure qb big enough, but not too big in top levels (from Sept '04)
        !qb(1:ifull,:)=max(qb(1:ifull,:),qgmin)
        !do k=kl-2,kl
        !  qb(1:ifull,k)=min(qb(1:ifull,k),10.*qgmin)
        !enddo

      end if ! ((mtimer>mtimeb).or.firstcall)

      ! Apply filter to model data using previously loaded host data
      if ((mtimer==mtimeb).and.(mod(nint(ktau*dt),60)==0)) then

        ! atmospheric nudging if required
        if (nud_p/=0.or.nud_t/=0.or.nud_uv/=0.or.nud_q/=0) then
          pslc(:)=pslb(:)-psl(1:ifull)
          uc(:,:)=ub(:,:)-u(1:ifull,:)
          vc(:,:)=vb(:,:)-v(1:ifull,:)
          tc(:,:)=tb(:,:)-t(1:ifull,:)
          qc(:,:)=qb(:,:)-qg(1:ifull,:)
          call getspecdata(pslc,uc,vc,tc,qc)
        end if

        ! specify sea-ice if not AMIP or Mixed-Layer-Ocean
        if(namip==0) then  ! namip SSTs/sea-ice take precedence
          if (nmlo==0) then
!           following sice updating code copied from nestin June '08      
!           check whether present ice points should change to/from sice points
            sicedep(:)=sicedepb(:)  ! from Jan 06
            fracice(:)=fraciceb(:)
!           because of new zs etc, ensure that sice is only over sea
            do iq=1,ifull
              if(fraciceb(iq)>0.)then
!               N.B. if already a sice point, keep present tice (in tggsn)
                if(fracice(iq)==0.)then
                  tggsn(iq,1)=min(271.2,tssb(iq),tb(iq,1)+.04*6.5) ! for 40 m lev1
                endif  ! (fracice(iq)==0.)
              endif  ! (fraciceb(iq)==0.)
              if(fracice(iq)<.02)fracice(iq)=0.
              if(land(iq))then
                sicedep(iq)=0.
                fracice(iq)=0.
              else
                if(fracice(iq)>0..and.sicedep(iq)==0.)then
!                 assign to 2. in NH and 1. in SH (according to spo)
!                 do this in indata, amipdata and nestin because of onthefly
                  if(rlatt(iq)>0.)then
                    sicedep(iq)=2.
                  else
                    sicedep(iq)=1.
                  endif ! (rlatt(iq)>0.)
                elseif(fracice(iq)==0..and.sicedep(iq)>0.)then  ! e.g. from Mk3  
                  fracice(iq)=1.
                endif  ! (fracice(iq)>0..and.sicedep(iq)==0.) .. elseif ..
              endif    ! (land(iq))
            enddo     ! iq loop

!           update tss 
            where (.not.land)
              tss=tssb
              tgg(:,1)=tss
            end where  ! (.not.land(iq))
          else
            ! nudge Mixed-Layer-Ocean
            if (nud_sst/=0.or.nud_sss/=0.or.nud_ouv/=0.or.
     &          nud_sfh/=0) then
              wl=wlev
              if (rdumg<-1.) then
                rduma=maxval(ocndep(:,1)) ! check for 3D data
                dumbb(1)=rduma
                call ccmpi_allreduce(dumbb(1:1),dumbb(2:2),"max",
     &                               comm_world)
                rdumg=dumbb(2)
              end if
              if (rdumg<0.5) wl=1
              if (wl==1) then ! switch to 2D data if 3D is missing
                call mloexport(0,timelt,1,0)
                sssb(:,1,1)=tssb
                where (fraciceb>0.)
                  sssb(:,1,1)=timelt
                end where
              end if
              call mlofilterhub(sssb(:,:,1),sssb(:,:,2),
     &                          sssb(:,:,3:4),ocndep(:,2),wl)
            end if
          end if ! (nmlo==0)
        end if ! (namip==0)
      end if ! (mod(nint(ktau*dt),60)==0)

      return
      end subroutine nestinb

      !--------------------------------------------------------------
      ! This subroutine gathers and distributes data for the
      ! scale-selective filter
      subroutine getspecdata(pslb,ub,vb,tb,qb)

      use arrays_m                     ! Atmosphere dyamics prognostic arrays
      use cc_mpi                       ! CC MPI routines
      use nharrs_m                     ! Non-hydrostatic atmosphere arrays
      use savuvt_m                     ! Saved dynamic arrays
      use savuv1_m                     ! Saved dynamic arrays
      use sigs_m                       ! Atmosphere sigma levels
      use vecsuv_m                     ! Map to cartesian coordinates
      use work3sav_m                   ! Water and tracer saved arrays
      use xyzinfo_m, only : x,y,z,wts  ! Grid coordinate arrays
      
      implicit none

      include 'newmpar.h'              ! Grid parameters
      include 'const_phys.h'           ! Physical constants
      include 'parm.h'                 ! Model configuration
      include 'parmgeom.h'             ! Coordinate data

      integer iq,k,ierr,kb,kln,klx,klt,klc
      real, dimension(ifull), intent(in) :: pslb
      real, dimension(ifull) :: costh,sinth
      real, dimension(ifull,kl), intent(in) :: ub,vb,tb,qb
      real, dimension(ifull,kblock) :: dum
      real, dimension(ifull) :: psld
      real, dimension(ifull,kblock) :: ud,vd,wd,td,qd
      real den,polenx,poleny,polenz,zonx,zony,zonz
      logical lblock

      ! nud_uv=0 (no preturbing of winds)
      ! nud_uv=1 (1D scale-selective filter)
      ! nud_uv=3 (JLM preturb zonal winds with 1D filter)
      ! nud_uv=9 (2D scale-selective filter)

      ! zonal wind option
      if (nud_uv==3) then
        polenx=-cos(rlat0*pi/180.)
        poleny=0.
        polenz=sin(rlat0*pi/180.)
        do iq=1,ifull
         zonx=            -polenz*y(iq)
         zony=polenz*x(iq)-polenx*z(iq)
         zonz=polenx*y(iq)
         den=sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) ) 
         costh(iq)= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
         sinth(iq)=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
        enddo
      end if
      
      ! store pressure for nudging
      psld=pslb
     
      ! Loop over maximum block size
      ! kblock can be reduced to save memory
      do kb=kbotdav,ktopdav,kblock
        if (myid == 0) then     
          write(6,*) "Gather data for spectral filter      ",kb
        end if      
        kln=kb                       ! lower limit of block
        klx=min(kb+kblock-1,ktopdav) ! upper limit of block
        klt=klx-kln+1                ! number of levels in block
        lblock=(kb==kbotdav)         ! flag for first loop (include psl)
        
        ! store winds, air temperature and water vapor for nudging
        if (nud_uv==3) then
          do k=kln,klx
            klc=k-kln+1
            ud(:,klc)=costh(:)*ub(:,k)  ! uzon
     &               -sinth(:)*vb(:,k)
          end do
        else
          ud(:,1:klt)=ub(:,kln:klx)
          vd(:,1:klt)=vb(:,kln:klx)
        end if
        td(:,1:klt)=tb(:,kln:klx)
        qd(:,1:klt)=qb(:,kln:klx)

        !-----------------------------------------------------------------------
        ! select nudging option
        if (nud_uv==9) then 
          if (myid==0) then
            write(6,*) "Two dimensional spectral filter      ",kb
          end if
          call slowspecmpi(.1*real(mbd)/(pi*schmidt)
     &                  ,psld(:),ud(:,1:klt),vd(:,1:klt)
     &                  ,td(:,1:klt),qd(:,1:klt),lblock,klt)
        else if (mod(6,nproc)==0.or.mod(nproc,6)==0) then
          if (myid==0) then
            write(6,*) "Separable 1D filter (MPI optimised)  ",kb
          end if
          call specfastmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                  ,psld(:),ud(:,1:klt),vd(:,1:klt)
     &                  ,td(:,1:klt),qd(:,1:klt),lblock,klt)
        else
          if (myid==0) then
            write(6,*) "Separable 1D filter (MPI)            ",kb
          end if
          call fourspecmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                  ,psld(:),ud(:,1:klt),vd(:,1:klt)
     &                  ,td(:,1:klt),qd(:,1:klt),lblock,klt)
        endif  ! (nud_uv<0) .. else ..
        !-----------------------------------------------------------------------

        if (myid==0) then
          write(6,*) "Distribute data from spectral filter ",kb
        end if
        
        if (nud_uv==3) then
          dum(:,1:klt)=ud(:,1:klt)
          do k=1,klt
            ud(1:ifull,k)=costh(:)*dum(:,k)
            vd(1:ifull,k)=-sinth(:)*dum(:,k)	  
          end do
        end if
        if (nud_uv/=0) then
          u(1:ifull,kln:klx)=u(1:ifull,kln:klx)+ud(:,1:klt)
          v(1:ifull,kln:klx)=v(1:ifull,kln:klx)+vd(:,1:klt)
          savu(1:ifull,kln:klx)=savu(1:ifull,kln:klx)+ud(:,1:klt)
          savu1(1:ifull,kln:klx)=savu1(1:ifull,kln:klx)+ud(:,1:klt)
          savu2(1:ifull,kln:klx)=savu2(1:ifull,kln:klx)+ud(:,1:klt)
          savv(1:ifull,kln:klx)=savv(1:ifull,kln:klx)+vd(:,1:klt)
          savv1(1:ifull,kln:klx)=savv1(1:ifull,kln:klx)+vd(:,1:klt)
          savv2(1:ifull,kln:klx)=savv2(1:ifull,kln:klx)+vd(:,1:klt)
        end if
        if (nud_t>0) then
          t(1:ifull,kln:klx)=t(1:ifull,kln:klx)+td(:,1:klt)
        end if
        if (nud_q>0) then
          qg(1:ifull,kln:klx)=max(qg(1:ifull,kln:klx)
     &     +qd(:,1:klt),0.)
          qgsav(:,kln:klx)=max(qgsav(:,kln:klx)
     &     +qd(:,1:klt),0.)
        end if
      
      end do
      
      if (nud_p>0) then
        psl(1:ifull)=psl(1:ifull)+psld(:)
        ps=1.e5*exp(psl(1:ifull))
      end if
      if (nud_t>0) then
        phi(:,1)=bet(1)*t(1:ifull,1)
        do k=2,kl
          phi(:,k)=phi(:,k-1)+bet(k)*t(1:ifull,k)
     &                      +betm(k)*t(1:ifull,k-1)
        end do
        phi=phi+phi_nh
      end if

      return
      end subroutine getspecdata

      !---------------------------------------------------------------------------------
      ! Slow 2D spectral downscaling - MPI version
      ! This option is an exact treatment of the filter
      subroutine slowspecmpi(cin,pslb,ub,vb,tb,qb,lblock,klt)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      use vecsuv_m          ! Map to cartesian coordinates
      use xyzinfo_m         ! Grid coordinate arrays
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration

      integer, intent(in) :: klt
      integer i,j,n,iq,iqg,k,ierr
      real, intent(in) :: cin
      real, dimension(ifull), intent(inout) :: pslb
      real, dimension(ifull,klt), intent(inout) :: ub,vb
      real, dimension(ifull,klt), intent(inout) :: tb,qb
      real, dimension(ifull_g) :: psls
      real, dimension(ifull_g,klt) :: uu,vv,ww
      real, dimension(ifull_g,klt) :: tt,qgg
      real, dimension(ifull_g) :: r,da,db
      real, dimension(klt) :: ud,vd,wd
      real cq,psum
      logical, intent(in) :: lblock

      cq=sqrt(4.5)*cin

      if (nud_p>0.and.lblock) then
        call ccmpi_gatherall(pslb(:), psls(:))
      end if
      if (nud_uv==3) then
        call ccmpi_gatherall(ub(:,1:klt),ww(:,1:klt))
      else if (nud_uv>0) then
        call ccmpi_gatherall(ub(:,1:klt),uu(:,1:klt))
        call ccmpi_gatherall(vb(:,1:klt),vv(:,1:klt))
        do k=1,klt
          da=ub(:,k)
          db=vb(:,k)
          uu(:,k)=ax_g(:)*da+bx_g(:)*db
          vv(:,k)=ay_g(:)*da+by_g(:)*db
          ww(:,k)=az_g(:)*da+bz_g(:)*db
        end do
      endif
      if (nud_t>0) then
        call ccmpi_gatherall(tb(:,1:klt),tt(:,1:klt))
      end if
      if (nud_q>0) then
        call ccmpi_gatherall(qb(:,1:klt),qgg(:,1:klt))
      end if

      if (myid==0.and.nmaxpr==1) then
        write(6,*) "Start 2D filter"
      end if
   
      do n=1,npan
        do j=1,jpan
          do i=1,ipan
            iqg=indg(i,j,n)
            iq =indp(i,j,n)
            r(:)=x_g(iqg)*x_g(:)+y_g(iqg)*y_g(:)+z_g(iqg)*z_g(:)
            r(:)=acos(max(min(r(:),1.),-1.))
            r(:)=exp(-(cq*r(:))**2)/(em_g(:)*em_g(:))
            psum=sum(r(:))
            if (nud_p>0.and.lblock) then
              pslb(iq)=sum(r(:)*psls(:))/psum
            end if
            if (nud_uv>0) then
              do k=1,klt
                ud(k)=sum(r(:)*uu(:,k))/psum
                vd(k)=sum(r(:)*vv(:,k))/psum
                wd(k)=sum(r(:)*ww(:,k))/psum
                ub(iq,k)=ax(iq)*ud(k)+ay(iq)*vd(k)+az(iq)*wd(k)
                vb(iq,k)=bx(iq)*ud(k)+by(iq)*vd(k)+bz(iq)*wd(k)
              end do
            end if
            if (nud_t>0) then
              do k=1,klt
                tb(iq,k)=sum(r(:)*tt(:,k))/psum
              end do
            end if
            if (nud_q>0) then
              do k=1,klt
                qb(iq,k)=sum(r(:)*qgg(:,k))/psum
              end do
            end if
          end do
        end do
      end do
 
      if (myid == 0.and.nmaxpr==1) write(6,*) "End 2D filter"

      return
      end subroutine slowspecmpi
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      ! Four pass spectral downscaling (symmetric)
      ! Used when mod(6,nproc)/=0 and mod(nproc,6)/=0 since it is slower than specfastmpi
      subroutine fourspecmpi(myid,cin,psls,uu,vv,tt,qgg,lblock,klt)
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: myid,klt
      integer pn,px,hproc,mproc,ns,ne,npta
      real, intent(in) :: cin
      real, dimension(ifull), intent(inout) :: psls
      real, dimension(ifull,klt), intent(inout) :: uu,vv
      real, dimension(ifull,klt), intent(inout) :: tt,qgg
      logical, intent(in) :: lblock
      
      npta=1                              ! number of panels per processor
      mproc=nproc                         ! number of processors per panel
      pn=0                                ! start panel
      px=5                                ! end panel
      hproc=0                             ! host processor for panel
      call procdiv(ns,ne,il_g,nproc,myid) ! number of rows per processor

      if (myid==0) then
        call spechost(mproc,hproc,npta,pn,px,ns,ne,cin,psls,uu,vv,
     &       tt,qgg,lblock,klt,ifull_g)
      else
        call spechost(mproc,hproc,npta,pn,px,ns,ne,cin,psls,uu,vv,
     &       tt,qgg,lblock,klt,0)
      end if

          
      return
      end subroutine fourspecmpi
      !---------------------------------------------------------------------------------


      !---------------------------------------------------------------------------------
      ! Four pass spectral downscaling (symmetric)
      ! MPI optimised for magic processor numbers 1,2,3,6,12,18,24,30,36,...
      ! (only works for mod(6,nproc)==0 or mod(nproc,6)==0)
      subroutine specfastmpi(myid,cin,psls,uu,vv,tt,qgg,lblock,klt)
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: myid,klt
      integer pn,px,hproc,mproc,ns,ne,npta,ifg
      real, intent(in) :: cin
      real, dimension(ifull), intent(inout) :: psls
      real, dimension(ifull,klt), intent(inout) :: uu,vv
      real, dimension(ifull,klt), intent(inout) :: tt,qgg
      logical, intent(in) :: lblock
      
      npta=max(6/nproc,1)                       ! number of panels per processor
      mproc=max(nproc/6,1)                      ! number of processors per panel
      pn=myid*npta/mproc                        ! start panel
      px=(myid+mproc)*npta/mproc-1              ! end panel
      hproc=pn*mproc/npta                       ! host processor for panel
      call procdiv(ns,ne,il_g,mproc,myid-hproc) ! number of rows per processor

      if (myid==hproc) then
        ifg=ifull_g
      else
        ifg=0
      end if
      if (npta==1) then
        call spechost_n(mproc,hproc,pn,ns,ne,cin,psls,uu,vv,
     &         tt,qgg,lblock,klt,ifg)
      else
        call spechost(mproc,hproc,npta,pn,px,ns,ne,cin,psls,uu,vv,
     &         tt,qgg,lblock,klt,ifg)
      end if

      return
      end subroutine specfastmpi
      !---------------------------------------------------------------------------------


      !---------------------------------------------------------------------------------
      ! This is the main routine for the scale-selective filter
      ! (see spechost_n for a reduced memory version)
      subroutine spechost(mproc,hproc,npta,pn,px,ns,ne,cin,pslb,
     &                    ub,vb,tb,qb,lblock,klt,ifg)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      use vecsuv_m          ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: mproc,hproc,npta,pn,px,ns,ne,klt,ifg
      integer k,ppass,iy,ppn,ppx,nne,nns,iproc,ierr
      integer n,a,b,c,til,colour,rank
      integer :: itag=0
      integer, save :: comm_host
      real, intent(in) :: cin
      real, dimension(ifull), intent(inout) :: pslb
      real, dimension(ifull,klt), intent(inout) :: ub,vb
      real, dimension(ifull,klt), intent(inout) :: tb,qb
      real, dimension(ifg) :: psls
      real, dimension(ifg,klt) :: uu,vv,ww
      real, dimension(ifg,klt) :: tt,qgg
      real, dimension(ifg) :: qp,qsum,zp,x_g,xx_g
      real, dimension(ifg,klt) :: qu,qv,qw,qt,qq
      real, dimension(ifg,klt) :: zu,zv,zw,zt,zq
      real, dimension(ifg*klt) :: dd
      real cq
      logical, intent(in) :: lblock
      logical, save :: first = .true.

      if (first) then
        if (myid==hproc) then
          colour=0
          rank=hproc/mproc
        else
          colour=1
          rank=myid-hproc
        end if
        call ccmpi_commsplit(comm_host,comm_world,colour,rank)
        first=.false.
      end if
      
      til=il_g*il_g
      cq=sqrt(4.5)*cin ! filter length scale

      ! gather data onto myid==0
      if (nud_p>0.and.lblock) then
        if (myid==0) then
          call ccmpi_gather(pslb(:), psls(:))
        else
          call ccmpi_gather(pslb(:))
        end if
      end if
      if (nud_uv==3) then
        if (myid==0) then
          call ccmpi_gather(ub(:,1:klt),ww(:,1:klt))
        else
          call ccmpi_gather(ub(:,1:klt))
        end if
      else if (nud_uv>0) then
        if (myid==0) then
          call ccmpi_gather(ub(:,1:klt),uu(:,1:klt))
          call ccmpi_gather(vb(:,1:klt),vv(:,1:klt))
          ! we assume that the coordinate transform on a single 
          ! processor is faster than sending three components 
          ! of the wind vector with gather
          do k=1,klt
            x_g =uu(:,k)
            xx_g=vv(:,k)
            uu(:,k)=ax_g(:)*x_g+bx_g(:)*xx_g
            vv(:,k)=ay_g(:)*x_g+by_g(:)*xx_g
            ww(:,k)=az_g(:)*x_g+bz_g(:)*xx_g
          end do
        else
          call ccmpi_gather(ub(:,1:klt))
          call ccmpi_gather(vb(:,1:klt))
        end if
      end if
      if (nud_t>0) then
        if (myid==0) then
          call ccmpi_gather(tb(:,1:klt),tt(:,1:klt))
        else
          call ccmpi_gather(tb(:,1:klt))
        end if
      end if
      if (nud_q>0) then
        if (myid==0) then
          call ccmpi_gather(qb(:,1:klt),qgg(:,1:klt))
        else
          call ccmpi_gather(qb(:,1:klt))
        end if
      end if

      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "Start 1D filter"

      ! distribute data over host processors
      if (myid==hproc) then
        qsum(:)=1./(em_g(:)*em_g(:))
        if (nud_p>0.and.lblock) then
          call ccmpi_bcast(psls,0,comm_host)
          psls(:)=psls(:)*qsum(:)
        end if
        if (nud_uv>0) then
          call ccmpi_bcast(uu,0,comm_host)
          call ccmpi_bcast(vv,0,comm_host)
          call ccmpi_bcast(ww,0,comm_host)
          do k=1,klt
            uu(:,k)=uu(:,k)*qsum(:)
            vv(:,k)=vv(:,k)*qsum(:)
            ww(:,k)=ww(:,k)*qsum(:)
          end do
        end if
        if (nud_t>0) then
          call ccmpi_bcast(tt,0,comm_host)
          do k=1,klt
            tt(:,k)=tt(:,k)*qsum(:)
          end do
        end if
        if (nud_q>0) then
          call ccmpi_bcast(qgg,0,comm_host)
          do k=1,klt
            qgg(:,k)=qgg(:,k)*qsum(:)
          end do
        end if
      end if

      do ppass=pn,px

        ! reset nudging fields for the next panel
        if (myid==hproc) then
          qsum(:)=1./(em_g(:)*em_g(:))
          if (nud_p>0.and.lblock) then
            qp=psls
          end if
          if (nud_uv>0) then
            do k=1,klt
              qu(:,k)=uu(:,k)
              qv(:,k)=vv(:,k)
              qw(:,k)=ww(:,k)
            end do
          end if
          if (nud_t>0) then
            do k=1,klt
              qt(:,k)=tt(:,k)
            end do
          end if
          if (nud_q>0) then
            do k=1,klt
              qq(:,k)=qgg(:,k)
            end do
          end if
        end if

        ! computations for the local processor group
        call speclocal(mproc,hproc,ns,ne,cq,ppass,qsum,qp,
     &         qu,qv,qw,qt,qq,lblock,klt,ifg)
        
        ! store the filtered results for this panel
        if (myid==hproc) then
          nns=ppass*til+1
          nne=ppass*til+til
          if (nud_p>0.and.lblock) then
            zp(nns:nne)=qp(nns:nne)/qsum(nns:nne)
          end if
          if (nud_uv>0) then
            do k=1,klt
              zu(nns:nne,k)=qu(nns:nne,k)/qsum(nns:nne)
              zv(nns:nne,k)=qv(nns:nne,k)/qsum(nns:nne)
              zw(nns:nne,k)=qw(nns:nne,k)/qsum(nns:nne)
            end do
          end if
          if (nud_t>0) then
            do k=1,klt
              zt(nns:nne,k)=qt(nns:nne,k)/qsum(nns:nne)
            end do
          end if
          if (nud_q>0) then
            do k=1,klt
              zq(nns:nne,k)=qq(nns:nne,k)/qsum(nns:nne)
            end do
          end if
        end if
        
      end do

      if (myid==0.and.nmaxpr==1) write(6,*) "End 1D filter"

      itag=itag+1
      if (myid == 0) then
        if (nmaxpr==1) write(6,*) 
     &    "Receive arrays from all host processors"
        do iproc=mproc,nproc-1,mproc
          ppn=iproc*npta/mproc
          ppx=(iproc+mproc)*npta/mproc-1
          iy=npta*til
          a=til
          c=-til*ppn
          if(nud_p>0.and.lblock)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do ppass=ppn,ppx
              do n=1,til
                zp(n+ppass*til)=dd(n+a*ppass+c)
              end do
            end do
          end if
          iy=npta*til*klt
          b=npta*til
          c=-til*(ppn+npta)
          if(nud_uv>0)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do ppass=ppn,ppx
                do n=1,til
                  zu(n+ppass*til,k)=dd(n+a*ppass+b*k+c)
                end do
              end do
            end do
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do ppass=ppn,ppx
                do n=1,til
                  zv(n+ppass*til,k)=dd(n+a*ppass+b*k+c)
                end do
              end do
            end do
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do ppass=ppn,ppx
                do n=1,til
                  zw(n+ppass*til,k)=dd(n+a*ppass+b*k+c)
                end do
              end do
            end do
          end if
          if(nud_t>0)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do ppass=ppn,ppx
                do n=1,til
                  zt(n+ppass*til,k)=dd(n+a*ppass+b*k+c)
                end do
              end do
            end do
          end if
          if(nud_q>0)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do ppass=ppn,ppx
                do n=1,til
                  zq(n+ppass*til,k)=dd(n+a*ppass+b*k+c)
                end do
              end do
            end do
          end if
        end do
      elseif (myid==hproc) then
        iy=npta*til
        a=til
        c=-til*pn
        if(nud_p>0.and.lblock)then
          do ppass=pn,px
            do n=1,til
              dd(n+a*ppass+c)=zp(n+ppass*til)
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
        iy=npta*til*klt
        b=npta*til
        c=-til*(pn+npta)
        if(nud_uv>0)then
          do k=1,klt
            do ppass=pn,px
              do n=1,til
                dd(n+a*ppass+b*k+c)=zu(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
          do k=1,klt
            do ppass=pn,px
              do n=1,til
                dd(n+a*ppass+b*k+c)=zv(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
          do k=1,klt
            do ppass=pn,px
              do n=1,til
                dd(n+a*ppass+b*k+c)=zw(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
        if(nud_t>0)then
          do k=1,klt
            do ppass=pn,px
              do n=1,til
                dd(n+a*ppass+b*k+c)=zt(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
        if(nud_q>0)then
          do k=1,klt
            do ppass=pn,px
              do n=1,til
                dd(n+a*ppass+b*k+c)=zq(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
      end if      
      
      ! distribute data to processors
      if (nud_p>0.and.lblock) then
        if (myid==0) then
          call ccmpi_distribute(pslb,zp)
        else
          call ccmpi_distribute(pslb)
        end if
      end if
      if (nud_uv==3) then
        if (myid==0) then
          call ccmpi_distribute(ub(:,1:klt),zw(:,1:klt))
        else
          call ccmpi_distribute(ub(:,1:klt))
        end if
      else if (nud_uv>0) then
        if (myid==0) then
          do k=1,klt
            x_g=ax_g(:)*zu(:,k)+ay_g(:)*zv(:,k)+az_g(:)*zw(:,k)
            xx_g=bx_g(:)*zu(:,k)+by_g(:)*zv(:,k)+bz_g(:)*zw(:,k)
            zu(:,k)=x_g
            zv(:,k)=xx_g
          end do
          call ccmpi_distribute(ub(:,1:klt),zu(:,1:klt))
          call ccmpi_distribute(vb(:,1:klt),zv(:,1:klt))
        else
          call ccmpi_distribute(ub(:,1:klt))
          call ccmpi_distribute(vb(:,1:klt))
        end if
      end if
      if (nud_t>0) then
        if (myid==0) then
          call ccmpi_distribute(tb(:,1:klt),zt(:,1:klt))
        else
          call ccmpi_distribute(tb(:,1:klt))
        end if
      end if
      if (nud_q>0) then
        if (myid==0) then
          call ccmpi_distribute(qb(:,1:klt),zq(:,1:klt))
        else
          call ccmpi_distribute(qb(:,1:klt))
        end if
      end if

      return
      end subroutine spechost
      !---------------------------------------------------------------------------------

      ! This version is for one panel per processor (reduced memory)
      subroutine spechost_n(mproc,hproc,pn,ns,ne,cin,pslb,
     &                    ub,vb,tb,qb,lblock,klt,ifg)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      use vecsuv_m          ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: mproc,hproc,pn,ns,ne,klt,ifg
      integer k,iy,ppn,nne,nns,iproc,ierr
      integer n,a,b,c,til,colour,rank
      integer :: itag=0
      integer, save :: comm_host
      real, intent(in) :: cin
      real cq
      real, dimension(ifull), intent(inout) :: pslb
      real, dimension(ifull,klt), intent(inout) :: ub,vb
      real, dimension(ifull,klt), intent(inout) :: tb,qb
      real, dimension(ifg) :: psls
      real, dimension(ifg,klt) :: uu,vv,ww
      real, dimension(ifg,klt) :: tt,qgg
      real, dimension(ifg) :: qsum,x_g,xx_g
      real, dimension(ifg*klt) :: dd
      logical, intent(in) :: lblock
      logical, save :: first = .true.

      if (first) then
        if (myid==hproc) then
          colour=0
          rank=hproc/mproc
        else
          colour=1
          rank=myid-hproc
        end if
        call ccmpi_commsplit(comm_host,comm_world,colour,rank)
        first=.false.
      end if
      
      til=il_g*il_g
      cq=sqrt(4.5)*cin ! filter length scale

      ! gather data onto myid==0
      if (nud_p>0.and.lblock) then
        if (myid==0) then
          call ccmpi_gather(pslb(:), psls(:))
        else
          call ccmpi_gather(pslb(:))
        end if
      end if
      if (nud_uv==3) then
        if (myid==0) then
          call ccmpi_gather(ub(:,1:klt),ww(:,1:klt))
        else
          call ccmpi_gather(ub(:,1:klt))
        end if
      else if (nud_uv>0) then
        if (myid==0) then
          call ccmpi_gather(ub(:,1:klt),uu(:,1:klt))
          call ccmpi_gather(vb(:,1:klt),vv(:,1:klt))
          do k=1,klt
            x_g =uu(:,k)
            xx_g=vv(:,k)
            uu(:,k)=ax_g(:)*x_g+bx_g(:)*xx_g
            vv(:,k)=ay_g(:)*x_g+by_g(:)*xx_g
            ww(:,k)=az_g(:)*x_g+bz_g(:)*xx_g
          end do
        else
          call ccmpi_gather(ub(:,1:klt))
          call ccmpi_gather(vb(:,1:klt))
        end if
      end if
      if (nud_t>0) then
        if (myid==0) then
          call ccmpi_gather(tb(:,1:klt),tt(:,1:klt))
        else
          call ccmpi_gather(tb(:,1:klt))
        end if
      end if
      if (nud_q>0) then
        if (myid==0) then
          call ccmpi_gather(qb(:,1:klt),qgg(:,1:klt))
        else
          call ccmpi_gather(qb(:,1:klt))
        end if
      end if

      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "Start 1D filter"

      ! distribute data over host processors
      if (myid==hproc) then
        qsum(:)=1./(em_g(:)*em_g(:))
        if (nud_p>0.and.lblock) then
          call ccmpi_bcast(psls,0,comm_host)
          psls(:)=psls(:)*qsum(:)
        end if
        if (nud_uv>0) then
          call ccmpi_bcast(uu,0,comm_host)
          call ccmpi_bcast(vv,0,comm_host)
          call ccmpi_bcast(ww,0,comm_host)
          do k=1,klt
            uu(:,k)=uu(:,k)*qsum(:)
            vv(:,k)=vv(:,k)*qsum(:)
            ww(:,k)=ww(:,k)*qsum(:)
          end do
        end if
        if (nud_t>0) then
          call ccmpi_bcast(tt,0,comm_host)
          do k=1,klt
            tt(:,k)=tt(:,k)*qsum(:)
          end do
        end if
        if (nud_q>0) then
          call ccmpi_bcast(qgg,0,comm_host)
          do k=1,klt
            qgg(:,k)=qgg(:,k)*qsum(:)
          end do
        end if
      end if

      ! computations for the local processor group
      call speclocal(mproc,hproc,ns,ne,cq,pn,qsum,psls,
     &         uu,vv,ww,tt,qgg,lblock,klt,ifg)

      ! store results on host processor
      if (myid==hproc) then        
        nns=pn*til+1
        nne=pn*til+til
        if (nud_p>0.and.lblock) then
          psls(nns:nne)=psls(nns:nne)/qsum(nns:nne)
        end if
        if (nud_uv>0) then
          do k=1,klt
            uu(nns:nne,k)=uu(nns:nne,k)/qsum(nns:nne)
            vv(nns:nne,k)=vv(nns:nne,k)/qsum(nns:nne)
            ww(nns:nne,k)=ww(nns:nne,k)/qsum(nns:nne)
          end do
        end if
        if (nud_t>0) then
          do k=1,klt
            tt(nns:nne,k)=tt(nns:nne,k)/qsum(nns:nne)
          end do
        end if
        if (nud_q>0) then
          do k=1,klt
            qgg(nns:nne,k)=qgg(nns:nne,k)/qsum(nns:nne)
          end do
        end if
      end if
        
      if (myid==0.and.nmaxpr==1) write(6,*) "End 1D filter"

      ! collect data from host processors
      itag=itag+1
      if (myid == 0) then
        if (nmaxpr==1) write(6,*) 
     &    "Receive arrays from all host processors"
        do iproc=mproc,nproc-1,mproc
          ppn=iproc/mproc
          iy=til
          if(nud_p>0.and.lblock)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do n=1,til
              psls(n+ppn*til)=dd(n)
            end do
          end if
          iy=til*klt
          b=til
          c=-til
          if(nud_uv>0)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do n=1,til
                uu(n+ppn*til,k)=dd(n+b*k+c)
              end do
            end do
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do n=1,til
                vv(n+ppn*til,k)=dd(n+b*k+c)
              end do
            end do
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do n=1,til
                ww(n+ppn*til,k)=dd(n+b*k+c)
              end do
            end do
          end if
          if(nud_t>0)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do n=1,til
                tt(n+ppn*til,k)=dd(n+b*k+c)
              end do
            end do
          end if
          if(nud_q>0)then
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do k=1,klt
              do n=1,til
                qgg(n+ppn*til,k)=dd(n+b*k+c)
              end do
            end do
          end if
        end do
      elseif (myid==hproc) then
        iy=til
        if(nud_p>0.and.lblock)then
          do n=1,til
            dd(n)=psls(n+pn*til)
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
        iy=til*klt
        b=til
        c=-til
        if(nud_uv>0)then
          do k=1,klt
            do n=1,til
              dd(n+b*k+c)=uu(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
          do k=1,klt
            do n=1,til
              dd(n+b*k+c)=vv(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
          do k=1,klt
            do n=1,til
              dd(n+b*k+c)=ww(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
        if(nud_t>0)then
          do k=1,klt
            do n=1,til
              dd(n+b*k+c)=tt(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
        if(nud_q>0)then
          do k=1,klt
            do n=1,til
              dd(n+b*k+c)=qgg(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(dd(1:iy),0,itag,comm_world)
        end if
      end if      
      
      ! distribute data to processors
      if (nud_p>0.and.lblock) then
        if (myid==0) then
          call ccmpi_distribute(pslb,psls)
        else
          call ccmpi_distribute(pslb)
        end if
      end if
      if (nud_uv==3) then
        if (myid==0) then
          call ccmpi_distribute(ub(:,1:klt),ww(:,1:klt))
        else
          call ccmpi_distribute(ub(:,1:klt))
        end if
      else if (nud_uv>0) then
        if (myid==0) then
          do k=1,klt
            x_g=ax_g(:)*uu(:,k)+ay_g(:)*vv(:,k)+az_g(:)*ww(:,k)
            xx_g=bx_g(:)*uu(:,k)+by_g(:)*vv(:,k)+bz_g(:)*ww(:,k)
            uu(:,k)=x_g
            vv(:,k)=xx_g
          end do
          call ccmpi_distribute(ub(:,1:klt),uu(:,1:klt))
          call ccmpi_distribute(vb(:,1:klt),vv(:,1:klt))
        else
          call ccmpi_distribute(ub(:,1:klt))
          call ccmpi_distribute(vb(:,1:klt))
        end if
      end if
      if (nud_t>0) then
        if (myid==0) then
          call ccmpi_distribute(tb(:,1:klt),tt(:,1:klt))
        else
          call ccmpi_distribute(tb(:,1:klt))
        end if
      end if
      if (nud_q>0) then
        if (myid==0) then
          call ccmpi_distribute(qb(:,1:klt),qgg(:,1:klt))
        else
          call ccmpi_distribute(qb(:,1:klt))
        end if
      end if
      
      return
      end subroutine spechost_n
      
      !---------------------------------------------------------------------------------
      ! This code runs between the local processors
      ! Code was moved to this subroutine to help the compiler vectorise the code
      ! This version handles uneven distribution of rows across processors
      subroutine speclocal(mproc,hproc,ns,ne,cq,ppass,qsum,
     &             qp,qu,qv,qw,qt,qq,lblock,klt,ifg)

      use cc_mpi            ! CC MPI routines
      use xyzinfo_m         ! Grid coordinate arrays

      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: mproc,hproc,ns,ne,ppass,klt,ifg
      integer j,k,n,ipass,iy
      integer iproc,ierr
      integer nne,nns,me
      integer a,b,d
      integer :: itag=0
      integer, dimension(4*il_g,il_g) :: igrd
      integer, dimension(0:3) :: maps
      real, intent(in) :: cq
      real, dimension(ifg), intent(inout) :: qp,qsum
      real, dimension(ifg,klt), intent(inout) :: qu,qv
      real, dimension(ifg,klt), intent(inout) :: qw
      real, dimension(ifg,klt), intent(inout) :: qt,qq
      real, dimension(4*il_g) :: ra,xa,ya,za
      real, dimension(4*il_g,ns:ne,klt) :: pu,pv,pw,pt,pq
      real, dimension(4*il_g,ns:ne,klt) :: au,av,aw,at,aq
      real, dimension(4*il_g,ns:ne) :: pp,psum
      real, dimension(4*il_g,ns:ne) :: ap,asum
      real, dimension(4*il_g*klt*(ne-ns+1)) :: dd
      logical, intent(in) :: lblock
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      
      do ipass=0,3
        me=maps(ipass)
        call getiqa(igrd(:,:),me,ipass,ppass,il_g)

        itag=itag+1
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "Receive arrays from local host"
        end if
        if (myid==hproc) then
          do iproc=hproc+1,mproc+hproc-1
            call procdiv(nns,nne,il_g,mproc,iproc-hproc)
            if (nns>nne) exit
            iy=me*(nne-nns+1)
            a=me
            d=-me*nns
            do j=nns,nne
              do n=1,me
                dd(n+a*j+d)=qsum(igrd(n,j))
              end do
            end do
            call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
            if(nud_p>0.and.lblock)then
              do j=nns,nne
                do n=1,me
                  dd(n+a*j+d)=qp(igrd(n,j))
                end do
              end do
              call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
            end if
            iy=me*(nne-nns+1)*klt
            b=me*(nne-nns+1)
            d=-me*(nne+1)
            if(nud_uv>0)then
              do k=1,klt
                do j=nns,nne
                  do n=1,me
                    dd(n+a*j+b*k+d)=qu(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
              do k=1,klt
                do j=nns,nne
                  do n=1,me
                    dd(n+a*j+b*k+d)=qv(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
              do k=1,klt
                do j=nns,nne
                  do n=1,me
                    dd(n+a*j+b*k+d)=qw(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
            end if
            if(nud_t>0)then
              do k=1,klt
                do j=nns,nne
                  do n=1,me
                    dd(n+a*j+b*k+d)=qt(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
            end if
            if(nud_q>0)then
              do k=1,klt
                do j=nns,nne
                  do n=1,me
                    dd(n+a*j+b*k+d)=qq(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
            end if
          end do
          do j=ns,ne
            asum(1:me,j)=qsum(igrd(1:me,j))
          end do
          if (nud_p>0.and.lblock) then
            do j=ns,ne
              ap(1:me,j)=qp(igrd(1:me,j))
            end do
          end if
          if (nud_uv>0) then
            do j=ns,ne
              au(1:me,j,:)=qu(igrd(1:me,j),:)
              av(1:me,j,:)=qv(igrd(1:me,j),:)
              aw(1:me,j,:)=qw(igrd(1:me,j),:)
            end do
          end if
          if (nud_t>0) then
            do j=ns,ne
              at(1:me,j,:)=qt(igrd(1:me,j),:)
            end do
          end if
          if (nud_q>0) then
            do j=ns,ne
              aq(1:me,j,:)=qq(igrd(1:me,j),:)
            end do
          end if
        else
          iy=me*(ne-ns+1)
          a=me
          d=-me*ns
          call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
          do j=ns,ne
            do n=1,me
              asum(n,j)=dd(n+a*j+d)
            end do
          end do
          if(nud_p>0.and.lblock)then
            call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
            do j=ns,ne
              do n=1,me
                ap(n,j)=dd(n+a*j+d)
              end do
            end do
          endif
          iy=me*(ne-ns+1)*klt
          b=me*(ne-ns+1)
          d=-me*(ne+1)
          if(nud_uv>0)then
            call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
            do k=1,klt
              do j=ns,ne
                do n=1,me
                  au(n,j,k)=dd(n+a*j+b*k+d)
                end do
              end do
            end do
            call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
            do k=1,klt
              do j=ns,ne
                do n=1,me
                  av(n,j,k)=dd(n+a*j+b*k+d)
                end do
              end do
            end do
            call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
            do k=1,klt
              do j=ns,ne
                do n=1,me
                  aw(n,j,k)=dd(n+a*j+b*k+d)
                end do
              end do
            end do
          end if
          if(nud_t>0)then
            call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
            do k=1,klt
              do j=ns,ne
                do n=1,me
                  at(n,j,k)=dd(n+a*j+b*k+d)
                end do
              end do
            end do
          end if
          if(nud_q>0)then
            call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
            do k=1,klt
              do j=ns,ne
                do n=1,me
                  aq(n,j,k)=dd(n+a*j+b*k+d)
                end do
              end do
            end do
          end if
        end if

        if (myid==0.and.nmaxpr==1) write(6,*) "Start convolution"

        do j=ns,ne
          xa(1:me)=x_g(igrd(1:me,j))
          ya(1:me)=y_g(igrd(1:me,j))
          za(1:me)=z_g(igrd(1:me,j))
          do n=1,il_g
            ra(1:me)=xa(n)*xa(1:me)+ya(n)*ya(1:me)+za(n)*za(1:me)
            ra(1:me)=acos(max(min(ra(1:me),1.),-1.))
            ra(1:me)=exp(-(cq*ra(1:me))**2)
            ! can also use the lines below which integrate the gaussian
            ! analytically over the length element (but slower)
            !ra(1)=2.*erf(cq*0.5*(ds/rearth)
            !ra(2:me)=erf(cq*(ra(2:me)+0.5*(ds/rearth)))  ! redefine ra(:) as wgt(:)
     &      !        -erf(cq*(ra(2:me)-0.5*(ds/rearth)))  ! (correct units are 1/cq)
            psum(n,j)=sum(ra(1:me)*asum(1:me,j))
            if (nud_p>0.and.lblock) then
              pp(n,j)=sum(ra(1:me)*ap(1:me,j))
            end if
            if (nud_uv>0) then
              do k=1,klt
                pu(n,j,k)=sum(ra(1:me)*au(1:me,j,k))
                pv(n,j,k)=sum(ra(1:me)*av(1:me,j,k))
                pw(n,j,k)=sum(ra(1:me)*aw(1:me,j,k))
              end do
            end if
            if (nud_t>0) then
              do k=1,klt
                pt(n,j,k)=sum(ra(1:me)*at(1:me,j,k))
              end do
            end if
            if (nud_q>0) then
              do k=1,klt
                pq(n,j,k)=sum(ra(1:me)*aq(1:me,j,k))
              end do
            end if
          end do
        end do

        if (myid==0.and.nmaxpr==1) write(6,*) "End convolution"

        itag=itag+1
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "Send arrays to local host"
        end if
        if (myid==hproc) then
          ! process points on host proc
          do j=ns,ne
            qsum(igrd(1:il_g,j))=psum(1:il_g,j)
          end do
          if (nud_p>0.and.lblock) then
            do j=ns,ne
              qp(igrd(1:il_g,j))=pp(1:il_g,j)
            end do
          end if
          if (nud_uv>0) then
            do j=ns,ne
              qu(igrd(1:il_g,j),:)=pu(1:il_g,j,:)
              qv(igrd(1:il_g,j),:)=pv(1:il_g,j,:)
              qw(igrd(1:il_g,j),:)=pw(1:il_g,j,:)
            end do
          end if
          if (nud_t>0) then
            do j=ns,ne
              qt(igrd(1:il_g,j),:)=pt(1:il_g,j,:)
            end do
          end if
          if (nud_q>0) then
            do j=ns,ne
              qq(igrd(1:il_g,j),:)=pq(1:il_g,j,:)
            end do
          end if
          ! process points from other procs
          do iproc=hproc+1,mproc+hproc-1
            call procdiv(nns,nne,il_g,mproc,iproc-hproc)
            if (nns>nne) exit
            iy=il_g*(nne-nns+1)
            a=il_g
            d=-il_g*nns
            call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
            do j=nns,nne
              do n=1,il_g
                qsum(igrd(n,j))=dd(n+a*j+d)
              end do
            end do
            if(nud_p>0.and.lblock)then
              call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
              do j=nns,nne
                do n=1,il_g
                  qp(igrd(n,j))=dd(n+a*j+d)
                end do
              end do
            end if
            iy=il_g*(nne-nns+1)*klt
            b=il_g*(nne-nns+1)
            d=-il_g*(nne+1)
            if(nud_uv>0)then
              call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
              do k=1,klt
                do j=nns,nne
                  do n=1,il_g
                    qu(igrd(n,j),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
              do k=1,klt
                do j=nns,nne
                  do n=1,il_g
                    qv(igrd(n,j),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
              do k=1,klt
                do j=nns,nne
                  do n=1,il_g
                    qw(igrd(n,j),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
            if(nud_t>0)then
              call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
              do k=1,klt
                do j=nns,nne
                  do n=1,il_g
                    qt(igrd(n,j),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
            if(nud_q>0)then
              call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
              do k=1,klt
                do j=nns,nne
                  do n=1,il_g
                    qq(igrd(n,j),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
          end do
       else
          iy=il_g*(ne-ns+1)
          a=il_g
          d=-il_g*ns
          do j=ns,ne
            do n=1,il_g
              dd(n+a*j+d)=psum(n,j)
            end do
          end do
          call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
          if(nud_p>0.and.lblock)then
            do j=ns,ne
              do n=1,il_g
                dd(n+a*j+d)=pp(n,j)
              end do
            end do
            call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
          end if
          iy=il_g*(ne-ns+1)*klt
          b=il_g*(ne-ns+1)
          d=-il_g*(ne+1)
          if(nud_uv>0)then
            do k=1,klt
              do j=ns,ne
                do n=1,il_g
                  dd(n+a*j+b*k+d)=pu(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
            do k=1,klt
              do j=ns,ne
                do n=1,il_g
                  dd(n+a*j+b*k+d)=pv(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
            do k=1,klt
              do j=ns,ne
                do n=1,il_g
                  dd(n+a*j+b*k+d)=pw(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
          end if
          if(nud_t>0)then
            do k=1,klt
              do j=ns,ne
                do n=1,il_g
                  dd(n+a*j+b*k+d)=pt(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
          end if
          if(nud_q>0)then
            do k=1,klt
              do j=ns,ne
                do n=1,il_g
                  dd(n+a*j+b*k+d)=pq(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
          end if
        end if  
      end do
      
      return  
      end subroutine speclocal

      !---------------------------------------------------------------------------------
      ! Map from 1D convolution to global index
      subroutine getiqa(iq,ne,ipass,ppass,il_g)
      
      implicit none
      
      integer, intent(in) :: ne,ipass,ppass,il_g
      integer, dimension(4*il_g,il_g), intent(out) :: iq
      integer sn,n,j,a,b,c
      
      do sn=1,ne,il_g

        select case(ppass*100+ipass*10+(sn-1)/il_g)
          case(0,300,520)                        ! panel 5   - x pass
            a=il_g
            b=-1
            c=5*il_g*il_g+1
          case(10,220,310)                       ! panel 2   - x pass
            a=1
            b=il_g
            c=il_g*(2*il_g-1)
          case(20,21,200,321,500)                ! panel 0,1 - y pass
            a=il_g
            b=1
            c=-il_g
          case(22)                               ! panel 3   - y pass
            a=1
            b=-il_g
            c=il_g*(4*il_g-2)
          case(23,323)                           ! panel 4   - y pass
            a=1
            b=-il_g
            c=il_g*(5*il_g-3)
          case(30,100,400)                       ! panel 0   - z pass
            a=1
            b=il_g
            c=-il_g
          case(31,223,232,331,523,532)           ! panel 2   - z pass ! panel 4   - x pass ! panel 3   - z pass
            a=il_g
            b=-1
            c=il_g*il_g+1
          case(32,332)                           ! panel 5   - z pass
            a=1
            b=il_g
            c=il_g*(5*il_g-3)
          case(110,222,330,410)                  ! panel 3   - z pass ! panel 5   - x pass
            a=il_g
            b=-1
            c=3*il_g*il_g+1
          case(120)                              ! panel 1   - x pass
            a=1
            b=il_g
            c=il_g*(il_g-1)
          case(121,421)                          ! panel 2   - x pass
            a=1
            b=il_g
            c=2*il_g*(il_g-1)
          case(122,123,230,423)                  ! panel 4,5 - x pass ! panel 2   - z pass
            a=il_g
            b=-1
            c=2*il_g*il_g+1
          case(130)                              ! panel 1   - y pass
            a=il_g
            b=1
            c=il_g*(il_g-1)
          case(131,431)                          ! panel 3   - y pass
            a=1
            b=-il_g
            c=il_g*(4*il_g-1)
          case(132,322,432)                      ! panel 0,1 - y pass
            a=il_g
            b=1
            c=-il_g*(2*il_g+1)
          case(210,320,510)                      ! panel 3   - y pass
            a=1
            b=-il_g
            c=4*il_g*il_g
          case(221,521)                          ! panel 1   - x pass
            a=1
            b=il_g
            c=il_g*(il_g-2)
          case(231,531)                          ! panel 0   - z pass
            a=1
            b=il_g
            c=-2*il_g
          case(420)                              ! panel 4   - x pass
            a=il_g
            b=-1
            c=4*il_g*il_g+1
          case(422)                              ! panel 1   - x pass
            a=1
            b=il_g
            c=il_g*(il_g-3)
          case(430)                              ! panel 4   - y pass
            a=1
            b=-il_g
            c=5*il_g*il_g
          case(522)                              ! panel 2   - x pass
            a=1
            b=il_g
            c=il_g*(2*il_g-3)
          case(530)                              ! panel 5   - z pass
            a=1
            b=il_g
            c=il_g*(5*il_g-1)            
          case DEFAULT
            write(6,*) "Invalid index ",ppass,ipass,sn,
     &              ppass*100+ipass*10+(sn-1)/il_g
            stop
        end select

        do j=1,il_g  
          do n=sn,sn+il_g-1
            iq(n,j)=a*n+b*j+c
          end do
        end do
  
      end do

      return
      end subroutine getiqa

      !---------------------------------------------------------------------------------
      ! Partition grid for each processor
      subroutine procdiv(ns,ne,ifull_g,nproc,myid)
      
      implicit none
      
      integer, intent(in) :: ifull_g,nproc,myid
      integer, intent(out) :: ns,ne
      integer npt,resid
      
      npt=ifull_g/nproc
      resid=mod(ifull_g,nproc)
      if ((myid+1)<=resid) then
        ns=myid*(npt+1)+1
        ne=(myid+1)*(npt+1)
      else
        ns=resid+myid*npt+1
        ne=resid+(myid+1)*npt
      end if
      
      return
      end subroutine procdiv

      !---------------------------------------------------------------------------------
      ! This subroutine gathers and distributes data for the
      ! MLO scale-selective filter
      subroutine mlofilterhub(sstb,sssb,suvb,sfh,wl)

      use cc_mpi                  ! CC MPI routines
      use mlo, only : mloimport,  ! Ocean physics and prognostic arrays
     &  mloexport,mloexpdep,wlev
      use mlodynamics             ! Ocean dynamics routines
      use soil_m                  ! Soil and surface data
      use vecsuv_m                ! Map to cartesian coordinates
      
      implicit none

      include 'newmpar.h'         ! Grid parameters      
      include 'parm.h'            ! Model configuration
      
      integer, intent(in) :: wl
      integer k,ka,kb,kc,ke,kln,klx,klt,kbb
      real, dimension(ifull), intent(in) :: sfh
      real, dimension(ifull,wlev), intent(in) :: sstb,sssb
      real, dimension(ifull,wlev,2), intent(in) :: suvb
      real, dimension(ifull,1) :: diffh_l
      real, dimension(ifull,kblock) :: diff_l,diffs_l
      real, dimension(ifull,kblock) :: diffu_l,diffv_l
      real, dimension(ifull) :: old,oldt
      real, dimension(ifull,ktopmlo:kbotmlo) :: rho,nsq
      logical lblock
      integer, parameter :: tempfix=1 ! delta temp (0=linear, 1=buoyancy)
      real, parameter :: rho0=1030.   ! linear density offset
      real, parameter :: a0=-0.3      ! linear density temp gradient
      real, parameter :: miss = 999999.
      
      
      kc=min(kbotmlo,ktopmlo+wl-1)

      if (nud_sfh/=0) then
        old=sfh
        call mloexport(4,old,0,0)
        where (.not.land)
          diffh_l(:,1)=sfh-old
        elsewhere
          diffh_l(:,1)=miss
        end where
      end if
      
      do kbb=ktopmlo,kc,kblock
      
        if (myid==0) then
          write(6,*) "Gather data for MLO filter ",kbb
        end if
            
        kln=kbb
        klx=min(kbb+kblock-1,kc)
        klt=klx-kln+1
        lblock=(kbb==ktopmlo)
      
        if (nud_sst/=0) then
          do k=kln,klx
            kb=k-kln+1
            old=sstb(:,k)
            call mloexport(0,old,k,0)
            where (.not.land)
              diff_l(:,kb)=sstb(:,k)-old
            elsewhere
              diff_l(:,kb)=miss
            end where
          end do
        end if

        if (nud_sss/=0) then
          do k=kln,klx
            kb=k-kln+1
            old=sssb(:,k)
            call mloexport(1,old,k,0)
            where (.not.land)
              diffs_l(:,kb)=sssb(:,k)-old
            elsewhere
              diffs_l(:,kb)=miss
            end where
          end do
        end if

        if (nud_ouv/=0) then
          do k=kln,klx
            kb=k-kln+1
            old=suvb(:,k,1)
            call mloexport(2,old,k,0)
            where (.not.land)
              diffu_l(:,kb)=suvb(:,k,1)-old
            elsewhere
              diffu_l(:,kb)=miss
            end where
          end do
          do k=kln,klx
            kb=k-kln+1
            old=suvb(:,k,2)
            call mloexport(3,old,k,0)
            where (.not.land)
              diffv_l(:,kb)=suvb(:,k,2)-old
            elsewhere
              diffv_l(:,kb)=miss
            end where
          end do
        end if

        if ((nud_uv/=9.and.abs(nmlo)/=1).or.namip/=0) then
          call mlofilterfast(diff_l(:,1:klt),diffs_l(:,1:klt),
     &                       diffu_l(:,1:klt),diffv_l(:,1:klt),
     &                       diffh_l(:,1),miss,lblock,klt)
        else
          call mlofilter(diff_l(:,1:klt),diffs_l(:,1:klt),
     &                   diffu_l(:,1:klt),diffv_l(:,1:klt),
     &                   diffh_l(:,1),miss,lblock,klt)
        end if

        if (myid==0) then
          write(6,*) "Distribute data for MLO filter ",kbb
        end if

        if (nud_sst/=0) then
          ! correct temp pertubation to minimise change in buoyancy
          if (tempfix==1.and.kc==1) then
            if (ktopmlo/=1) then
              write(6,*) "ERROR: nud_sst with SST input"
              write(6,*) "requires ktopmlo=1"
              stop
            end if
            old=293.
            do k=1,kbotmlo
              call mloexport(0,old,k,0)
              rho(:,k)=rho0+a0*old ! linear approximation to density
            end do
            do k=1,kbotmlo-1
              !nsq=-2.*grav*(rho(:,k)-rho(:,k+1))/((dep(:,k+1)-dep(:,k))*(rho(:,k)+rho(:,k+1)))
              nsq(:,k)=-(rho(:,k)-rho(:,k+1))/(rho(:,k)+rho(:,k+1))
              !nsq(:,k)=max(nsq(:,k),0.)
            end do
            ! nsq(:,k-1)=-(a0*(oldt-old))/(2.*rho0+a0*(oldt+old))
            ! old*(1.-nsq(:,k-1))=oldt*(1.+nsq(:,k-1))+2.*rho0/a0*nsq(:,k-1)
            call mloexport(0,old,1,0)
            old=old+diff_l(:,1)*10./real(mloalpha)
            old=max(old,271.)
            call mloimport(0,old,1,0)
            oldt=old
            do k=2,kbotmlo
              old=(oldt*(1.+nsq(:,k-1))+2.*nsq(:,k-1)*rho0/a0)
     &          /(1.-nsq(:,k-1))
              old=max(old,271.)  
              call mloimport(0,old,k,0)
              oldt=old
            end do
          else
            do k=kln,klx
              ka=min(wl,k)
              kb=k-kln+1
              old=sstb(:,ka)
              call mloexport(0,old,k,0)
              old=old+diff_l(:,kb)*10./real(mloalpha)
              old=max(old,271.)
              call mloimport(0,old,k,0)
            end do
            if (klx==kc) then
              do k=kc+1,kbotmlo
                old=sstb(:,ka)
                call mloexport(0,old,k,0)
                old=old+diff_l(:,kb)*10./real(mloalpha) ! kb saved from above loop
                old=max(old,271.)	  
                call mloimport(0,old,k,0)
              end do
            end if
          end if
        end if

        if (nud_sss/=0) then
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=sssb(:,ka)
            call mloexport(1,old,k,0)
            old=old+diffs_l(:,kb)*10./real(mloalpha)
            old=max(old,0.)
            call mloimport(1,old,k,0)
          end do
          if (klx==kc) then
            do k=kc+1,kbotmlo
              old=sssb(:,ka)
              call mloexport(1,old,k,0)
              old=old+diffs_l(:,kb)*10./real(mloalpha) ! kb saved from above loop
              old=max(old,0.)
              call mloimport(1,old,k,0)
            end do
          end if
        end if

        if (nud_ouv/=0) then
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=suvb(:,ka,1)
            call mloexport(2,old,k,0)
            old=old+diffu_l(:,kb)*10./real(mloalpha)
            call mloimport(2,old,k,0)
            if (allocated(oldu1)) then
              oldu1(:,k)=oldu1(:,k)+diffu_l(:,kb)*10./real(mloalpha)
              oldu2(:,k)=oldu2(:,k)+diffu_l(:,kb)*10./real(mloalpha)
            end if
          end do
          if (klx==kc) then
            do k=kc+1,kbotmlo
              old=suvb(:,ka,1)
              call mloexport(2,old,k,0)
              old=old+diffu_l(:,kb)*10./real(mloalpha) ! kb saved from above loop
              call mloimport(2,old,k,0)
              if (allocated(oldu1)) then
                oldu1(:,k)=oldu1(:,k)+diffu_l(:,kb)*10./real(mloalpha)
                oldu2(:,k)=oldu2(:,k)+diffu_l(:,kb)*10./real(mloalpha)
              end if
            end do
          end if
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=suvb(:,ka,2)
            call mloexport(3,old,k,0)
            old=old+diffv_l(:,kb)*10./real(mloalpha)
            call mloimport(3,old,k,0)
            if (allocated(oldv1)) then
              oldv1(:,k)=oldv1(:,k)+diffv_l(:,kb)*10./real(mloalpha)
              oldv2(:,k)=oldv2(:,k)+diffv_l(:,kb)*10./real(mloalpha)
            end if
          end do
          if (klx==kc) then
            do k=kc+1,kbotmlo
              old=suvb(:,ka,2)
              call mloexport(3,old,k,0)
              old=old+diffv_l(:,kb)*10./real(mloalpha)
              call mloimport(3,old,k,0)
              if (allocated(oldv1)) then
                oldv1(:,k)=oldv1(:,k)+diffv_l(:,kb)*10./real(mloalpha)
                oldv2(:,k)=oldv2(:,k)+diffv_l(:,kb)*10./real(mloalpha)
              end if
            end do
          end if
        end if
      
      end do
     
      if (nud_sfh/=0) then
        old=sfh
        call mloexport(4,old,0,0)
        old=old+diffh_l(:,1)*10./real(mloalpha)
        call mloimport(4,old,0,0)
      end if

      return
      end subroutine mlofilterhub
      
      !---------------------------------------------------------------------------------
      ! 2D Filter for MLO 
      subroutine mlofilter(diff_l,diffs_l,diffu_l,diffv_l,
     &                     diffh_l,miss,lblock,kd)

      use cc_mpi                  ! CC MPI routines
      use vecsuv_m                ! Map to cartesian coordinates

      implicit none

      include 'newmpar.h'         ! Grid parameters
      include 'parm.h'            ! Model configuration

      integer, intent(in) :: kd
      integer k,ierr
      real, intent(in) :: miss
      real, dimension(ifull,1), intent(inout) :: diffh_l
      real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
      real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
      real, dimension(ifull_g,1) :: diffh_g
      real, dimension(ifull_g,kd) :: diff_g,diffs_g
      real, dimension(ifull_g,kd) :: diffu_g,diffv_g
      real, dimension(ifull_g,kd) :: diffw_g
      real, dimension(ifull_g) :: x_g, xx_g
      logical, intent(in) :: lblock
      logical, dimension(ifull_g) :: landg

      if (nud_sst/=0) then
        call ccmpi_gatherall(diff_l(:,1:kd),diff_g(:,1:kd))
        landg=abs(diff_g(:,1)-miss)<0.1
      end if
      if (nud_sss/=0) then
        call ccmpi_gatherall(diffs_l(:,1:kd),diffs_g(:,1:kd))
        landg=abs(diffs_g(:,1)-miss)<0.1
      end if
      if (nud_ouv/=0) then
        call ccmpi_gatherall(diffu_l(:,1:kd),diffu_g(:,1:kd))
        call ccmpi_gatherall(diffv_l(:,1:kd),diffv_g(:,1:kd))
        do k=1,kd
          x_g =diffu_g(:,k)
          xx_g=diffv_g(:,k)
          diffu_g(:,k)=ax_g*x_g+bx_g*xx_g
          diffv_g(:,k)=ay_g*x_g+by_g*xx_g
          diffw_g(:,k)=az_g*x_g+bz_g*xx_g
          where (abs(x_g-miss)<0.1)
            diffu_g(:,k)=miss
            diffv_g(:,k)=miss
            diffw_g(:,k)=miss
          end where
        end do        
        landg=abs(diffw_g(:,1)-miss)<0.1
      end if
      if (nud_sfh/=0.and.lblock) then
        call ccmpi_gatherall(diffh_l(:,1),diffh_g(:,1))
        landg=abs(diffh_g(:,1)-miss)<0.1
      end if

      if (myid==0) then
        write(6,*) "MLO 2D scale-selective filter"
        if (kd==1) then
          write(6,*) "Single level filter"
        else
          write(6,*) "Multiple level filter"
        end if
      end if

      call mlofilterhost(diff_g,diffs_g,
     &                   diffu_g,diffv_g,diffw_g,
     &                   diffh_g,kd,miss,landg,lblock)

      if (myid==0.and.nmaxpr==1) then
        write(6,*) "MLO end 2D filter"
      end if
      
      if (nud_sst/=0) then
        diff_l=diff_g(1:ifull,:)
      end if
      if (nud_sss/=0) then
        diffs_l=diffs_g(1:ifull,:)
      end if
      if (nud_ouv/=0) then
        diffu_l=diffu_g(1:ifull,:)
        diffv_l=diffv_g(1:ifull,:)
      end if
      if (nud_sfh/=0.and.lblock) then
        diffh_l=diffh_g(1:ifull,:)
      end if

      return
      end subroutine mlofilter

      subroutine mlofilterhost(diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                         diffh_g,kd,miss,landg,lblock)

      use cc_mpi             ! CC MPI routines
      use map_m              ! Grid map arrays
      use vecsuv_m           ! Map to cartesian coordinates
      use xyzinfo_m          ! Grid coordinate arrays

      implicit none

      include 'newmpar.h'    ! Grid parameters
      include 'const_phys.h' ! Physical constants
      include 'parm.h'       ! Model configuration
      include 'parmgeom.h'   ! Coordinate data

      integer, intent(in) :: kd
      integer i,j,n,iqq,iqqg,k
      real, intent(in) :: miss
      real nsum,cq
      real, dimension(ifull_g,1), intent(inout) :: diffh_g
      real, dimension(ifull_g,kd), intent(inout) :: diff_g,diffs_g
      real, dimension(ifull_g,kd), intent(inout) :: diffu_g,diffv_g
      real, dimension(ifull_g,kd), intent(inout) :: diffw_g
      real, dimension(ifull_g) :: rr,mm,nn
      real, dimension(ifull) :: ddh
      real, dimension(ifull,kd) :: dd,dds,ddu,ddv
      real, dimension(kd) :: tu,tv,tw
      logical, intent(in) :: lblock
      logical, dimension(ifull_g), intent(in) :: landg

      ! eventually will be replaced with mbd once full ocean coupling is complete
      cq=sqrt(4.5)*.1*real(max(nud_sst,nud_sss,nud_ouv,nud_sfh,mbd))
     &   /(pi*schmidt)

      dd=0.
      dds=0.
      ddu=0.
      ddv=0.
      ddh=0.
      mm=1./(em_g*em_g)
      where(.not.landg)
        nn=mm
      elsewhere
        nn=0.
      end where
      if (nud_sst/=0) then
        do k=1,kd
          diff_g(:,k)=diff_g(:,k)*nn
        end do
      end if
      if (nud_sss/=0) then
        do k=1,kd
          diffs_g(:,k)=diffs_g(:,k)*nn
        end do
      end if
      if (nud_ouv/=0) then
        do k=1,kd
          diffu_g(:,k)=diffu_g(:,k)*nn
          diffv_g(:,k)=diffv_g(:,k)*nn
          diffw_g(:,k)=diffw_g(:,k)*nn
        end do
      end if
      if (nud_sfh/=0.and.lblock) then
        diffh_g(:,1)=diffh_g(:,1)*nn
      end if
      do n=1,npan
        do j=1,jpan
          do i=1,ipan
            iqqg=indg(i,j,n)
            iqq=indp(i,j,n)
            if (.not.landg(iqqg)) then
              rr(:)=x_g(iqqg)*x_g(:)+y_g(iqqg)*y_g(:)+z_g(iqqg)*z_g(:)
              rr(:)=acos(max(min(rr(:),1.),-1.))
              rr(:)=exp(-(cq*rr(:))**2)
              nsum=sum(rr(:)*mm(:))
              if (nud_sst/=0) then
                do k=1,kd
                  dd(iqq,k)=sum(rr(:)*diff_g(:,k))/nsum
                end do
              end if
              if (nud_sss/=0) then
                do k=1,kd
                  dds(iqq,k)=sum(rr(:)*diffs_g(:,k))/nsum
                end do
              end if
              if (nud_ouv/=0) then
                do k=1,kd
                  tu(k)=sum(rr(:)*diffu_g(:,k))/nsum
                  tv(k)=sum(rr(:)*diffv_g(:,k))/nsum
                  tw(k)=sum(rr(:)*diffw_g(:,k))/nsum
                  ddu(iqq,k)=ax(iqq)*tu(k)+ay(iqq)*tv(k)+az(iqq)*tw(k)
                  ddv(iqq,k)=bx(iqq)*tu(k)+by(iqq)*tv(k)+bz(iqq)*tw(k)
                end do
              end if
              if (nud_sfh/=0.and.lblock) then
                ddh(iqq)=sum(rr(:)*diffh_g(:,1))/nsum
              end if
            end if
          end do
        end do
      end do

      if (nud_sst/=0) then
        diff_g(1:ifull,:)=dd(:,:)
      end if
      if (nud_sss/=0) then
        diffs_g(1:ifull,:)=dds(:,:)
      end if
      if (nud_ouv/=0) then
        diffu_g(1:ifull,:)=ddu(:,:)
        diffv_g(1:ifull,:)=ddv(:,:)
      end if
      if (nud_sfh/=0.and.lblock) then
        diffh_g(1:ifull,1)=ddh(:)
      end if

      return
      end subroutine mlofilterhost

      ! 1D filer for mlo
      subroutine mlofilterfast(diff_l,diffs_l,diffu_l,diffv_l,
     &                         diffh_l,miss,lblock,kd)

      use cc_mpi                  ! CC MPI routines

      implicit none

      include 'newmpar.h'         ! Grid parameters
      include 'const_phys.h'      ! Physical constants
      include 'parm.h'            ! Model configuration
      include 'parmgeom.h'        ! Coordinate data

      integer, intent(in) :: kd
      integer pn,px,hproc,mproc,ns,ne,npta,ifg
      real, intent(in) :: miss      
      real, dimension(ifull,1), intent(inout) :: diffh_l
      real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
      real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
      real cq
      logical, intent(in) :: lblock
      
      ! eventually will be replaced with mbd once full ocean coupling is complete
      cq=sqrt(4.5)*.1*real(max(nud_sst,nud_sss,nud_ouv,nud_sfh,mbd))
     &   /(pi*schmidt)
      
      if(mod(6,nproc)==0.or.mod(nproc,6)==0)then
        if (myid==0) then
          write(6,*) "MLO 1D scale-selective filter (MPI optimised)"
          if (kd==1) then
            write(6,*) "Single level filter"
          else
            write(6,*) "Multiple level filter"
          end if
        end if
        npta=max(6/nproc,1)                       ! number of panels per processor
        mproc=max(nproc/6,1)                      ! number of processors per panel
        pn=myid*npta/mproc                        ! start panel
        px=pn+npta-1                              ! end panel
        hproc=pn*mproc/npta                       ! host processor for panel
        call procdiv(ns,ne,il_g,mproc,myid-hproc) ! number of rows per processor
      else
        if (myid==0) then
          write(6,*) "MLO 1D scale-selective filter (MPI)"
          if (kd==1) then
            write(6,*) "Single level filter"
          else
            write(6,*) "Multiple level filter"
          end if
        end if        
        npta=1                              ! number of panels per processor
        mproc=nproc                         ! number of processors per panel
        pn=0                                ! start panel
        px=5                                ! end panel
        hproc=0                             ! host processor for panel
        call procdiv(ns,ne,il_g,nproc,myid) ! number of rows per processor
      end if

      if (myid==hproc) then
        ifg=ifull_g
      else
        ifg=0
      end if
      if (pn==px) then
        call mlospechost_n(mproc,hproc,pn,ns,ne,cq,
     &                   diff_l,diffs_l,diffu_l,diffv_l,
     &                   diffh_l,miss,lblock,kd,ifg)
      else
        call mlospechost(mproc,hproc,npta,pn,px,ns,ne,cq,
     &                   diff_l,diffs_l,diffu_l,diffv_l,
     &                   diffh_l,miss,lblock,kd,ifg)
      end if

      if (myid==0.and.nmaxpr==1) then
        write(6,*) "MLO end 1D filter"
      end if

      return
      end subroutine mlofilterfast

      subroutine mlospechost(mproc,hproc,npta,pn,px,ns,ne,cq,
     &                       diff_l,diffs_l,diffu_l,diffv_l,
     &                       diffh_l,miss,lblock,kd,ifg)

      use cc_mpi             ! CC MPI routines
      use map_m              ! Grid map arrays
      use vecsuv_m           ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: mproc,hproc,npta,pn,px,ns,ne,kd,ifg
      integer ppass,iy,ppn,ppx,nne,nns,iproc,ierr
      integer n,a,b,c,k,til,colour,rank
      integer :: itag=0
      integer, save :: comm_host
      real, intent(in) :: cq,miss
      real, dimension(ifull,1), intent(inout) :: diffh_l
      real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
      real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
      real, dimension(ifg,1) :: diffh_g
      real, dimension(ifg,kd) :: diff_g,diffs_g
      real, dimension(ifg,kd) :: diffu_g,diffv_g
      real, dimension(ifg,kd) :: diffw_g
      real, dimension(ifg) :: qsum,rsum,zph,qph,x_g,xx_g
      real, dimension(ifg,kd) :: zp,zps,zpu,zpv,zpw
      real, dimension(ifg,kd) :: qp,qps,qpu,qpv,qpw
      real, dimension(ifg*kd) :: zz
      logical, intent(in) :: lblock
      logical, dimension(ifg) :: landg
      logical, save :: first = .true.

      if (first) then
        if (myid==hproc) then
          colour=0
          rank=hproc/mproc
        else
          colour=1
          rank=myid-hproc
        end if
        call ccmpi_commsplit(comm_host,comm_world,colour,rank)
        first=.false.
      end if
      
      til=il_g*il_g 

      ! gather data onto myid==0
      if (nud_sst/=0) then
        if (myid==0) then
          call ccmpi_gather(diff_l(:,1:kd),diff_g(:,1:kd))
        else
          call ccmpi_gather(diff_l(:,1:kd))
        end if
      end if
      if (nud_sss/=0) then
        if (myid==0) then
          call ccmpi_gather(diffs_l(:,1:kd),diffs_g(:,1:kd))
        else
          call ccmpi_gather(diffs_l(:,1:kd))
        end if
      end if
      if (nud_ouv/=0) then
        if (myid==0) then
          call ccmpi_gather(diffu_l(:,1:kd),diffu_g(:,1:kd))
          call ccmpi_gather(diffv_l(:,1:kd),diffv_g(:,1:kd))
          do k=1,kd
            x_g =diffu_g(:,k)
            xx_g=diffv_g(:,k)
            diffu_g(:,k)=ax_g*x_g+bx_g*xx_g
            diffv_g(:,k)=ay_g*x_g+by_g*xx_g
            diffw_g(:,k)=az_g*x_g+bz_g*xx_g
            where (abs(x_g-miss)<0.1)
              diffu_g(:,k)=miss
              diffv_g(:,k)=miss
              diffw_g(:,k)=miss
            end where
          end do
        else
          call ccmpi_gather(diffu_l(:,1:kd),diffu_g(:,1:kd))
          call ccmpi_gather(diffv_l(:,1:kd),diffv_g(:,1:kd))
        end if        
      end if
      if (nud_sfh/=0.and.lblock) then
        if (myid==0) then
          call ccmpi_gather(diffh_l(:,1),diffh_g(:,1))
        else
          call ccmpi_gather(diffh_l(:,1))
        end if
      end if
     
      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "MLO Start 1D filter"

      ! distribute data over host processors
      if (myid==hproc) then
        qsum(:)=1./(em_g(:)*em_g(:))
        if (nud_sst/=0) then
          call ccmpi_bcast(diff_g,0,comm_host)
          landg=abs(diff_g(:,1)-miss)<0.1
        end if
        if (nud_sss/=0) then
          call ccmpi_bcast(diffs_g,0,comm_host)
          landg=abs(diffs_g(:,1)-miss)<0.1
        end if
        if (nud_ouv/=0) then
          call ccmpi_bcast(diffu_g,0,comm_host)
          call ccmpi_bcast(diffv_g,0,comm_host)
          call ccmpi_bcast(diffw_g,0,comm_host)
          landg=abs(diffw_g(:,1)-miss)<0.1
        end if
        if (nud_sfh/=0.and.lblock) then
          call ccmpi_bcast(diffh_g,0,comm_host)
          landg=abs(diffh_g(:,1)-miss)<0.1
        end if        
        
        where(.not.landg) ! land/sea mask
          rsum(:)=qsum(:)
        elsewhere
          rsum(:)=0.
        end where

        if (nud_sst/=0) then
          do k=1,kd
            diff_g(:,k)=diff_g(:,k)*rsum
          end do
        end if
        if (nud_sss/=0) then
          do k=1,kd
            diffs_g(:,k)=diffs_g(:,k)*rsum
          end do
        end if 
        if (nud_ouv/=0) then
          do k=1,kd
            diffu_g(:,k)=diffu_g(:,k)*rsum
            diffv_g(:,k)=diffv_g(:,k)*rsum
            diffw_g(:,k)=diffw_g(:,k)*rsum
          end do
        end if
        if (nud_sfh/=0.and.lblock) then
          diffh_g(:,1)=diffh_g(:,1)*rsum
        end if

        zp=0.
        zps=0.
        zpu=0.
        zpv=0.
        zpw=0.
        zph=0.
      end if

      do ppass=pn,px

        if (myid==hproc) then
          qsum(:)=1./(em_g(:)*em_g(:))

          if (nud_sst/=0) then
            do k=1,kd
              qp(:,k)=diff_g(:,k)
            end do
          end if
          if (nud_sss/=0) then
            do k=1,kd
              qps(:,k)=diffs_g(:,k)
            end do
          end if 
          if (nud_ouv/=0) then
            do k=1,kd
              qpu(:,k)=diffu_g(:,k)
              qpv(:,k)=diffv_g(:,k)
              qpw(:,k)=diffw_g(:,k)
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            qph=diffh_g(:,1)
          end if
        end if

        ! computations for the local processor group
        call mlospeclocal(mproc,hproc,ns,ne,cq,ppass,qsum,
     &                      qp,qps,qpu,qpv,qpw,qph,kd,lblock,
     &                      ifg)
        
        ! store results on host processors
        if (myid==hproc) then
          nns=ppass*til+1
          nne=ppass*til+til
          if (nud_sst/=0) then
            do k=1,kd
              where (qsum(nns:nne)>1.E-8)
                zp(nns:nne,k)=qp(nns:nne,k)/qsum(nns:nne)
              end where
            end do
          end if
          if (nud_sss/=0) then
            do k=1,kd
              where (qsum(nns:nne)>1.E-8)
                zps(nns:nne,k)=qps(nns:nne,k)/qsum(nns:nne)
              end where
            end do
          end if
          if (nud_ouv/=0) then
            do k=1,kd
              where (qsum(nns:nne)>1.E-8)
                zpu(nns:nne,k)=qpu(nns:nne,k)/qsum(nns:nne)
                zpv(nns:nne,k)=qpv(nns:nne,k)/qsum(nns:nne)
                zpw(nns:nne,k)=qpw(nns:nne,k)/qsum(nns:nne)
              end where
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            where (qsum(nns:nne)>1.E-8)
              zph(nns:nne)=qph(nns:nne)/qsum(nns:nne)
            end where
          end if
        end if
        
      end do

      if (myid==0.and.nmaxpr==1) write(6,*) "MLO End 1D filter"

      itag=itag+1
      if (myid == 0) then
        if (nmaxpr==1) write(6,*) "MLO Receive arrays from all proc"
        do iproc=mproc,nproc-1,mproc
          ppn=iproc*npta/mproc
          ppx=(iproc+mproc)*npta/mproc-1
          iy=npta*til*kd
          a=til
          b=npta*til
          c=-til*(ppn+npta)
          if (nud_sst/=0) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do ppass=ppn,ppx
                do n=1,til
                  zp(n+ppass*til,k)=zz(n+a*ppass+b*k+c)
                end do
              end do
            end do
          end if
          if (nud_sss/=0) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do ppass=ppn,ppx
                do n=1,til
                  zps(n+ppass*til,k)=zz(n+a*ppass+b*k+c)
                end do
              end do
            end do
          end if
          if (nud_ouv/=0) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do ppass=ppn,ppx
                do n=1,til
                  zpu(n+ppass*til,k)=zz(n+a*ppass+b*k+c)
                end do
              end do
            end do
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do ppass=ppn,ppx
                do n=1,til
                  zpv(n+ppass*til,k)=zz(n+a*ppass+b*k+c)
                end do
              end do
            end do
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do ppass=ppn,ppx
                do n=1,til
                  zpw(n+ppass*til,k)=zz(n+a*ppass+b*k+c)
                end do
              end do
            end do
          end if
          iy=npta*til
          a=til
          c=-til*ppn
          if (nud_sfh/=0.and.lblock) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do ppass=ppn,ppx
              do n=1,til
                zph(n+ppass*til)=zz(n+a*ppass+c)
              end do
            end do
          end if
        end do
      elseif (myid==hproc) then
        iy=npta*til*kd
        a=til
        b=npta*til
        c=-til*(pn+npta)
        if (nud_sst/=0) then
          do k=1,kd
            do ppass=pn,px
              do n=1,til
                zz(n+a*ppass+b*k+c)=zp(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
        if (nud_sss/=0) then
          do k=1,kd
            do ppass=pn,px
              do n=1,til
                zz(n+a*ppass+b*k+c)=zps(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
        if (nud_ouv/=0) then
          do k=1,kd
            do ppass=pn,px
              do n=1,til
                zz(n+a*ppass+b*k+c)=zpu(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
          do k=1,kd
            do ppass=pn,px
              do n=1,til
                zz(n+a*ppass+b*k+c)=zpv(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
          do k=1,kd
            do ppass=pn,px
              do n=1,til
                zz(n+a*ppass+b*k+c)=zpw(n+ppass*til,k)
              end do
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
        iy=npta*til
        a=til
        c=-til*pn
        if (nud_sfh/=0.and.lblock) then
          do ppass=pn,px
            do n=1,til
              zz(n+a*ppass+c)=zph(n+ppass*til)
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
      end if      
      
      ! distribute data to processors
      if (nud_sst/=0) then
        if (myid == 0) then
          call ccmpi_distribute(diff_l(:,1:kd), zp(:,1:kd))
        else
          call ccmpi_distribute(diff_l(:,1:kd))
        end if
      end if
      if (nud_sss/=0) then
        if (myid == 0) then
          call ccmpi_distribute(diffs_l(:,1:kd), zps(:,1:kd))
        else
          call ccmpi_distribute(diffs_l(:,1:kd))
        end if
      end if
      if (nud_ouv/=0) then
        if (myid == 0) then
          do k=1,kd
            x_g =ax_g*zpu(:,k)+ay_g*zpv(:,k)
     &        +az_g*zpw(:,k)
            xx_g=bx_g*zpu(:,k)+by_g*zpv(:,k)
     &        +bz_g*zpw(:,k)
            zpu(:,k)=x_g
            zpv(:,k)=xx_g
          end do
          call ccmpi_distribute(diffu_l(:,1:kd), zpu(:,1:kd))
          call ccmpi_distribute(diffv_l(:,1:kd), zpv(:,1:kd))
        else
          call ccmpi_distribute(diffu_l(:,1:kd))
          call ccmpi_distribute(diffv_l(:,1:kd))
        end if
      end if
      if (nud_sfh/=0.and.lblock) then
        if (myid == 0) then
          call ccmpi_distribute(diffh_l(:,1), zph(:))
        else
          call ccmpi_distribute(diffh_l(:,1))
        end if
      end if

      return
      end subroutine mlospechost
      !---------------------------------------------------------------------------------
      ! memory reduced version
      subroutine mlospechost_n(mproc,hproc,pn,ns,ne,cq,
     &                       diff_l,diffs_l,diffu_l,diffv_l,
     &                       diffh_l,miss,lblock,kd,ifg)

      use cc_mpi             ! CC MPI routines
      use map_m              ! Grid map arrays
      use vecsuv_m          ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: mproc,hproc,pn,ns,ne,kd,ifg
      integer iy,ppn,nne,nns,iproc,ierr
      integer n,a,b,c,k,til,colour,rank
      integer :: itag=0
      integer, save :: comm_host
      real, intent(in) :: cq,miss
      real, dimension(ifull,1), intent(inout) :: diffh_l
      real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
      real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
      real, dimension(ifg,1) :: diffh_g
      real, dimension(ifg,kd) :: diff_g,diffs_g
      real, dimension(ifg,kd) :: diffu_g,diffv_g,diffw_g
      real, dimension(ifg) :: qsum,rsum,x_g,xx_g
      real, dimension(ifg*kd) :: zz
      logical, intent(in) :: lblock
      logical, dimension(ifg) :: landg
      logical, save :: first = .true.

      if (first) then
        if (myid==hproc) then
          colour=0
          rank=hproc/mproc
        else
          colour=1
          rank=myid-hproc
        end if
        call ccmpi_commsplit(comm_host,comm_world,colour,rank)
        first=.false.
      end if
      
      til=il_g*il_g 

      ! gather data onto myid==0
      if (nud_sst/=0) then
        if (myid==0) then
          call ccmpi_gather(diff_l(:,1:kd),diff_g(:,1:kd))
        else
          call ccmpi_gather(diff_l(:,1:kd))
        end if
      end if
      if (nud_sss/=0) then
        if (myid==0) then
          call ccmpi_gather(diffs_l(:,1:kd),diffs_g(:,1:kd))
        else
          call ccmpi_gather(diffs_l(:,1:kd))
        end if
      end if
      if (nud_ouv/=0) then
        if (myid==0) then
          call ccmpi_gather(diffu_l(:,1:kd),diffu_g(:,1:kd))
          call ccmpi_gather(diffv_l(:,1:kd),diffv_g(:,1:kd))
          do k=1,kd
            x_g =diffu_g(:,k)
            xx_g=diffv_g(:,k)
            diffu_g(:,k)=ax_g*x_g+bx_g*xx_g
            diffv_g(:,k)=ay_g*x_g+by_g*xx_g
            diffw_g(:,k)=az_g*x_g+bz_g*xx_g
            where (abs(x_g-miss)<0.1)
              diffu_g(:,k)=miss
              diffv_g(:,k)=miss
              diffw_g(:,k)=miss
            end where
          end do
        else
          call ccmpi_gather(diffu_l(:,1:kd),diffu_g(:,1:kd))
          call ccmpi_gather(diffv_l(:,1:kd),diffv_g(:,1:kd))
        end if        
      end if
      if (nud_sfh/=0.and.lblock) then
        if (myid==0) then
          call ccmpi_gather(diffh_l(:,1),diffh_g(:,1))
        else
          call ccmpi_gather(diffh_l(:,1))
        end if
      end if
      
      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "MLO Start 1D filter"

      ! distribute data over host processors
      if (myid==hproc) then
        qsum(:)=1./(em_g(:)*em_g(:))
        if (nud_sst/=0) then
          call ccmpi_bcast(diff_g,0,comm_host)
          landg=abs(diff_g(:,1)-miss)<0.1
        end if
        if (nud_sss/=0) then
          call ccmpi_bcast(diffs_g,0,comm_host)
          landg=abs(diffs_g(:,1)-miss)<0.1
        end if
        if (nud_ouv/=0) then
          call ccmpi_bcast(diffu_g,0,comm_host)
          call ccmpi_bcast(diffv_g,0,comm_host)
          call ccmpi_bcast(diffw_g,0,comm_host)
          landg=abs(diffw_g(:,1)-miss)<0.1
        end if
        if (nud_sfh/=0.and.lblock) then
          call ccmpi_bcast(diffh_g,0,comm_host)
          landg=abs(diffh_g(:,1)-miss)<0.1
        end if        
        
        where(.not.landg) ! land/sea mask
          rsum(:)=qsum(:)
        elsewhere
          rsum(:)=0.
        end where

        if (nud_sst/=0) then
          do k=1,kd
            diff_g(:,k)=diff_g(:,k)*rsum
          end do
        end if
        if (nud_sss/=0) then
          do k=1,kd
            diffs_g(:,k)=diffs_g(:,k)*rsum
          end do
        end if 
        if (nud_ouv/=0) then
          do k=1,kd
            diffu_g(:,k)=diffu_g(:,k)*rsum
            diffv_g(:,k)=diffv_g(:,k)*rsum
            diffw_g(:,k)=diffw_g(:,k)*rsum
          end do
        end if
        if (nud_sfh/=0.and.lblock) then
          diffh_g(:,1)=diffh_g(:,1)*rsum
        end if
      end if

      ! computations for the local processor group
      call mlospeclocal(mproc,hproc,ns,ne,cq,pn,qsum,
     &                    diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                    diffh_g,kd,lblock,ifg)
      
      ! store results on host processor  
      if (myid==hproc) then
        nns=pn*til+1
        nne=pn*til+til
        if (nud_sst/=0) then
          do k=1,kd
            where (qsum(nns:nne)>1.E-8)
              diff_g(nns:nne,k)=diff_g(nns:nne,k)/qsum(nns:nne)
            elsewhere
              diff_g(nns:nne,k)=0.
            end where
          end do
        end if
        if (nud_sss/=0) then
          do k=1,kd
            where (qsum(nns:nne)>1.E-8)
              diffs_g(nns:nne,k)=diffs_g(nns:nne,k)/qsum(nns:nne)
            elsewhere
              diffs_g(nns:nne,k)=0.
            end where
          end do
        end if
        if (nud_ouv/=0) then
          do k=1,kd
            where (qsum(nns:nne)>1.E-8)
              diffu_g(nns:nne,k)=diffu_g(nns:nne,k)/qsum(nns:nne)
              diffv_g(nns:nne,k)=diffv_g(nns:nne,k)/qsum(nns:nne)
              diffw_g(nns:nne,k)=diffw_g(nns:nne,k)/qsum(nns:nne)
            elsewhere
              diffu_g(nns:nne,k)=0.
              diffv_g(nns:nne,k)=0.
              diffw_g(nns:nne,k)=0.
            end where
          end do
        end if
        if (nud_sfh/=0.and.lblock) then
          where (qsum(nns:nne)>1.E-8)
            diffh_g(nns:nne,1)=diffh_g(nns:nne,1)/qsum(nns:nne)
          elsewhere
            diffh_g(nns:nne,1)=0.
          end where
        end if
      end if

      if (myid==0.and.nmaxpr==1) write(6,*) "MLO End 1D filter"

      itag=itag+1
      if (myid == 0) then
        if (nmaxpr==1) write(6,*) "MLO Receive arrays from all proc"
        do iproc=mproc,nproc-1,mproc
          ppn=iproc/mproc
          iy=til*kd
          b=til
          c=-til
          if (nud_sst/=0) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do n=1,til
                diff_g(n+ppn*til,k)=zz(n+b*k+c)
              end do
            end do
          end if
          if (nud_sss/=0) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do n=1,til
                diffs_g(n+ppn*til,k)=zz(n+b*k+c)
              end do
             end do
          end if
          if (nud_ouv/=0) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do n=1,til
                diffu_g(n+ppn*til,k)=zz(n+b*k+c)
              end do
            end do
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do n=1,til
                diffv_g(n+ppn*til,k)=zz(n+b*k+c)
              end do
            end do
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do k=1,kd
              do n=1,til
                diffw_g(n+ppn*til,k)=zz(n+b*k+c)
              end do
            end do
          end if
          iy=til
          if (nud_sfh/=0.and.lblock) then
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do n=1,til
              diffh_g(n+ppn*til,1)=zz(n)
            end do
          end if
        end do
      elseif (myid==hproc) then
        iy=til*kd
        b=til
        c=-til
        if (nud_sst/=0) then
          do k=1,kd
            do n=1,til
              zz(n+b*k+c)=diff_g(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
        if (nud_sss/=0) then
          do k=1,kd
            do n=1,til
              zz(n+b*k+c)=diffs_g(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
        if (nud_ouv/=0) then
          do k=1,kd
            do n=1,til
              zz(n+b*k+c)=diffu_g(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
          do k=1,kd
            do n=1,til
              zz(n+b*k+c)=diffv_g(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
          do k=1,kd
            do n=1,til
              zz(n+b*k+c)=diffw_g(n+pn*til,k)
            end do
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
        iy=til
        if (nud_sfh/=0.and.lblock) then
          do n=1,til
            zz(n)=diffh_g(n+pn*til,1)
          end do
          call ccmpi_ssend(zz(1:iy),0,itag,comm_world)
        end if
      end if      

      ! distribute data to processors
      if (nud_sst/=0) then
        if (myid == 0) then
          call ccmpi_distribute(diff_l(:,1:kd), diff_g(:,1:kd))
        else
          call ccmpi_distribute(diff_l(:,1:kd))
        end if
      end if
      if (nud_sss/=0) then
        if (myid == 0) then
          call ccmpi_distribute(diffs_l(:,1:kd), diffs_g(:,1:kd))
        else
          call ccmpi_distribute(diffs_l(:,1:kd))
        end if
      end if
      if (nud_ouv/=0) then
        if (myid == 0) then
          do k=1,kd
            x_g =ax_g*diffu_g(:,k)+ay_g*diffv_g(:,k)
     &        +az_g*diffw_g(:,k)
            xx_g=bx_g*diffu_g(:,k)+by_g*diffv_g(:,k)
     &        +bz_g*diffw_g(:,k)
            diffu_g(:,k)=x_g
            diffv_g(:,k)=xx_g
          end do
          call ccmpi_distribute(diffu_l(:,1:kd), diffu_g(:,1:kd))
          call ccmpi_distribute(diffv_l(:,1:kd), diffv_g(:,1:kd))
        else
          call ccmpi_distribute(diffu_l(:,1:kd))
          call ccmpi_distribute(diffv_l(:,1:kd))
        end if
      end if
      if (nud_sfh/=0.and.lblock) then
        if (myid == 0) then
          call ccmpi_distribute(diffh_l(:,1:1), diffh_g(:,1:1))
        else
          call ccmpi_distribute(diffh_l(:,1:1))
        end if
      end if

      return
      end subroutine mlospechost_n

      
      !---------------------------------------------------------------------------------
      ! This version is for asymmetric decomposition
      subroutine mlospeclocal(mproc,hproc,ns,ne,cq,ppass,
     &             qsum,qp,qps,qpu,qpv,qpw,qph,kd,lblock,
     &             ifg)

      use cc_mpi             ! CC MPI routines
      use xyzinfo_m          ! Grid coordinate arrays
     
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: mproc,hproc,ns,ne,ppass,kd,ifg
      integer j,n,ipass,iy
      integer iproc,ierr
      integer nne,nns,me
      integer a,b,d,k
      integer :: itag=0
      integer, dimension(4*il_g,il_g) :: igrd
      integer, dimension(0:3) :: maps
      real, intent(in) :: cq
      real, dimension(ifg), intent(inout) :: qph,qsum
      real, dimension(ifg,kd), intent(inout) :: qp,qps
      real, dimension(ifg,kd), intent(inout) :: qpu,qpv
      real, dimension(ifg,kd), intent(inout) :: qpw
      real, dimension(4*il_g) :: rr,xa,ya,za
      real, dimension(4*il_g,ns:ne) :: asum,psum
      real, dimension(4*il_g,ns:ne) :: aph,pph
      real, dimension(4*il_g,ns:ne,kd) :: ap,aps,pp,pps
      real, dimension(4*il_g,ns:ne,kd) :: apu,apv,apw,ppu,ppv,ppw
      real, dimension(4*il_g*kd*(ne-ns+1)) :: zz
      logical, intent(in) :: lblock
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      
      do ipass=0,3
        me=maps(ipass)
        call getiqa(igrd(:,:),me,ipass,ppass,il_g)

        itag=itag+1
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "MLO Recieve arrays from local host"
        end if
        if (myid==hproc) then
          do iproc=hproc+1,mproc+hproc-1
            call procdiv(nns,nne,il_g,mproc,iproc-hproc)
            if (nns>nne) exit
            iy=me*(nne-nns+1)
            a=me
            d=-me*nns
            do j=nns,nne
              do n=1,me
                zz(n+a*j+d)=qsum(igrd(n,j))
              end do
            end do
            call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
            if (nud_sfh/=0.and.lblock) then
              do j=nns,nne
                do n=1,me
                  zz(n+a*j+d)=qph(igrd(n,j))
                end do
              end do
              call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
            end if
            iy=me*(nne-nns+1)*kd
            b=me*(nne-nns+1)
            d=-me*(nne+1)
            if (nud_sst/=0) then
              do k=1,kd
                do j=nns,nne
                  do n=1,me
                    zz(n+a*j+b*k+d)=qp(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
            end if
            if (nud_sss/=0) then
              do k=1,kd
                do j=nns,nne
                  do n=1,me
                    zz(n+a*j+b*k+d)=qps(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
            end if
            if (nud_ouv/=0) then
              do k=1,kd
                do j=nns,nne
                  do n=1,me
                    zz(n+a*j+b*k+d)=qpu(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
              do k=1,kd
                do j=nns,nne
                  do n=1,me
                    zz(n+a*j+b*k+d)=qpv(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
              do k=1,kd
                do j=nns,nne
                  do n=1,me
                    zz(n+a*j+b*k+d)=qpw(igrd(n,j),k)
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
            end if
          end do
          ! prepare data on host processor
          do j=ns,ne
            asum(1:me,j)=qsum(igrd(1:me,j))
          end do
          if (nud_sst/=0) then
            do j=ns,ne
              ap(1:me,j,:)=qp(igrd(1:me,j),:)
            end do
          end if
          if (nud_sss/=0) then
            do j=ns,ne
              aps(1:me,j,:)=qps(igrd(1:me,j),:)
            end do
          end if
          if (nud_ouv/=0) then
            do j=ns,ne
              apu(1:me,j,:)=qpu(igrd(1:me,j),:)
              apv(1:me,j,:)=qpv(igrd(1:me,j),:)
              apw(1:me,j,:)=qpw(igrd(1:me,j),:)
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            do j=ns,ne
              aph(1:me,j)=qph(igrd(1:me,j))
            end do
          end if
        else
          iy=me*(ne-ns+1)
          a=me
          d=-me*ns
          call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
          do j=ns,ne
            do n=1,me
              asum(n,j)=zz(n+a*j+d)
            end do
          end do
          if (nud_sfh/=0.and.lblock) then
            call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
            do j=ns,ne
              do n=1,me
                aph(n,j)=zz(n+a*j+d)
              end do
            end do
          end if
          iy=me*(ne-ns+1)*kd
          b=me*(ne-ns+1)
          d=-me*(ne+1)
          if (nud_sst/=0) then
            call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
            do k=1,kd
              do j=ns,ne
                do n=1,me
                  ap(n,j,k)=zz(n+a*j+b*k+d)
                end do
              end do
            end do
          end if
          if (nud_sss/=0) then
            call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
            do k=1,kd
              do j=ns,ne
                do n=1,me
                  aps(n,j,k)=zz(n+a*j+b*k+d)
                end do
              end do
            end do  
          end if
          if (nud_ouv/=0) then
            call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
            do k=1,kd
              do j=ns,ne
                do n=1,me
                  apu(n,j,k)=zz(n+a*j+b*k+d)
                end do
              end do
            end do
            call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
            do k=1,kd
              do j=ns,ne
                do n=1,me
                  apv(n,j,k)=zz(n+a*j+b*k+d)
                end do
              end do
            end do
            call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
            do k=1,kd
              do j=ns,ne
                do n=1,me
                  apw(n,j,k)=zz(n+a*j+b*k+d)
                end do
              end do
            end do  
          end if          
        end if

        if (myid==0.and.nmaxpr==1) write(6,*) "MLO start conv"

        do j=ns,ne
          xa(1:me)=x_g(igrd(1:me,j))
          ya(1:me)=y_g(igrd(1:me,j))
          za(1:me)=z_g(igrd(1:me,j))
          do n=1,il_g
            rr(1:me)=xa(n)*xa(1:me)+ya(n)*ya(1:me)+za(n)*za(1:me)
            rr(1:me)=acos(max(min(rr(1:me),1.),-1.))
            rr(1:me)=exp(-(cq*rr(1:me))**2)
            psum(n,j)=sum(rr(1:me)*asum(1:me,j))
            if (nud_sst/=0) then
              do k=1,kd
                pp(n,j,k)=sum(rr(1:me)*ap(1:me,j,k))
              end do
            end if
            if (nud_sss/=0) then
              do k=1,kd
                pps(n,j,k)=sum(rr(1:me)*aps(1:me,j,k))
              end do
            end if
            if (nud_ouv/=0) then
              do k=1,kd
                ppu(n,j,k)=sum(rr(1:me)*apu(1:me,j,k))
                ppv(n,j,k)=sum(rr(1:me)*apv(1:me,j,k))
                ppw(n,j,k)=sum(rr(1:me)*apw(1:me,j,k))
              end do
            end if
            if (nud_sfh/=0.and.lblock) then
              pph(n,j)=sum(rr(1:me)*aph(1:me,j))
            end if
          end do
        end do

        if (myid==0.and.nmaxpr==1) write(6,*) "MLO end conv"

        itag=itag+1
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "MLO Send arrays to local host"
        end if
        if (myid==hproc) then
          ! process points on host proc
          do j=ns,ne
            qsum(igrd(1:il_g,j))=psum(1:il_g,j)
          end do
          if (nud_sst/=0) then
            do j=ns,ne
              qp(igrd(1:il_g,j),:)=pp(1:il_g,j,:)
            end do
          end if
          if (nud_sss/=0) then
            do j=ns,ne
              qps(igrd(1:il_g,j),:)=pps(1:il_g,j,:)
            end do
          end if
          if (nud_ouv/=0) then
            do j=ns,ne
              qpu(igrd(1:il_g,j),:)=ppu(1:il_g,j,:)
              qpv(igrd(1:il_g,j),:)=ppv(1:il_g,j,:)
              qpw(igrd(1:il_g,j),:)=ppw(1:il_g,j,:)
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            do j=ns,ne
              qph(igrd(1:il_g,j))=pph(1:il_g,j)
            end do
          end if
          ! process points from other procs
          do iproc=hproc+1,mproc+hproc-1
            call procdiv(nns,nne,il_g,mproc,iproc-hproc)
            if (nns>nne) exit
            iy=il_g*(nne-nns+1)
            a=il_g
            d=-il_g*nns
            call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
            do j=nns,nne
              do n=1,il_g
                qsum(igrd(n,j))=zz(n+a*j+d)
              end do
            end do
            if (nud_sfh/=0.and.lblock) then
              call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
              do j=nns,nne
                do n=1,il_g
                  qph(igrd(n,j))=zz(n+a*j+d)
                end do
              end do
            end if
            iy=il_g*(nne-nns+1)*kd
            b=il_g*(nne-nns+1)
            d=-il_g*(nne+1)
            if (nud_sst/=0) then
              call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
              do k=1,kd
                do j=nns,nne
                  do n=1,il_g
                    qp(igrd(n,j),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
            if (nud_sss/=0) then
              call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
              do k=1,kd
                do j=nns,nne
                  do n=1,il_g
                    qps(igrd(n,j),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
            if (nud_ouv/=0) then
              call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
              do k=1,kd
                do j=nns,nne
                  do n=1,il_g
                    qpu(igrd(n,j),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
              do k=1,kd
                do j=nns,nne
                  do n=1,il_g
                    qpv(igrd(n,j),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
              do k=1,kd
                do j=nns,nne
                  do n=1,il_g
                    qpw(igrd(n,j),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
          end do
        else
          iy=il_g*(ne-ns+1)
          a=il_g
          d=-il_g*ns
          do j=ns,ne
            do n=1,il_g
              zz(n+a*j+d)=psum(n,j)
            end do
          end do
          call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
          if (nud_sfh/=0.and.lblock) then
            do j=ns,ne
              do n=1,il_g
                zz(n+a*j+d)=pph(n,j)
              end do
            end do
            call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
          end if
          iy=il_g*(ne-ns+1)*kd
          b=il_g*(ne-ns+1)
          d=-il_g*(ne+1)
          if (nud_sst/=0) then
            do k=1,kd
              do j=ns,ne
                do n=1,il_g
                  zz(n+a*j+b*k+d)=pp(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
          end if
          if (nud_sss/=0) then
            do k=1,kd
              do j=ns,ne
                do n=1,il_g
                  zz(n+a*j+b*k+d)=pps(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
          end if
          if (nud_ouv/=0) then
            do k=1,kd
              do j=ns,ne
                do n=1,il_g
                  zz(n+a*j+b*k+d)=ppu(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
            do k=1,kd
              do j=ns,ne
                do n=1,il_g
                  zz(n+a*j+b*k+d)=ppv(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
            do k=1,kd
              do j=ns,ne
                do n=1,il_g
                  zz(n+a*j+b*k+d)=ppw(n,j,k)
                end do
              end do
            end do
            call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
          end if
        end if
          
      end do
      
      return  
      end subroutine mlospeclocal

      ! Relaxtion method for mlo
      subroutine mlonudge(new,sssb,suvb,sfh,wl)

      use mlo, only : mloimport, ! Ocean physics and prognostic arrays
     &  mloexport,wlev
      
      implicit none

      include 'newmpar.h'        ! Grid parameters
      include 'parm.h'           ! Model configuration

      integer, intent(in) :: wl
      integer k,ka,i
      real, dimension(ifull), intent(in) :: sfh
      real, dimension(ifull,wlev), intent(in) :: new,sssb
      real, dimension(ifull,wlev,2), intent(in) :: suvb
      real, dimension(ifull) :: old
      real wgt
      
      wgt=dt/real(nud_hrs*3600)
      if (nud_sst/=0) then
        do k=ktopmlo,kbotmlo
          ka=min(k,wl)
          old=new(:,ka)
          call mloexport(0,old,k,0)
          old=old*(1.-wgt)+new(:,ka)*wgt
          old=max(old,271.)
          call mloimport(0,old,k,0)
        end do
      end if
      
      if (nud_sss/=0) then
        do k=ktopmlo,kbotmlo
          ka=min(k,wl)
          old=sssb(:,ka)
          call mloexport(1,old,k,0)
          old=old*(1.-wgt)+sssb(:,ka)*wgt
          old=max(old,0.)	  
          call mloimport(1,old,k,0)
        end do
      end if
      
      if (nud_ouv/=0) then
        do i=2,3
          do k=ktopmlo,kbotmlo
            ka=min(k,wl)
            old=suvb(:,ka,i-1)
            call mloexport(i,old,k,0)
            old=old*(1.-wgt)+suvb(:,ka,i-1)*wgt
            call mloimport(i,old,k,0)
          end do
        end do
      end if

      if (nud_sfh/=0) then
        old=sfh
        call mloexport(4,old,0,0)
        old=old*(1.-wgt)+sfh*wgt
        call mloimport(4,old,0,0)
      end if
      
      return
      end subroutine mlonudge
      
