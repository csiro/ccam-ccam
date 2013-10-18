      ! CCAM nudging/assimilation routines
      
      ! These routines preturb the regional model with the large scale circulation of the host model.
      ! Currently, relaxiation, far-field and scale-selective filter options are supported for both
      ! the atmosphere and ocean.
      
      ! We support both 1D and 2D versions of the scale-selective filter.  2D is exact, but expensive.
      ! Current tests suggest the 1D is a good approximation of the 2D filter where the grid stretching
      ! is not too large.

      ! nbd/=0       Far-field or relaxation nudging
      ! mbd/=0       Spectral filter (1D and 2D versions, see nud_uv)
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
#ifdef debug
        if (mydiag) then
          write (6,"('zsb# nestin  ',9f7.1)") diagvals(zsb)
          write (6,"('tssb# nestin ',9f7.1)") diagvals(tssb) 
        end if
#endif
   
        ! determine time corrosponding to new host nudging data
        kdhour=ktime_r/100-ktime/100
        kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
        mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     &               +60*kdhour+kdmin

#ifdef debug
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
#endif

!       ensure qb big enough, but not too big in top levels (from Sept '04)
        qb(1:ifull,:)=max(qb(1:ifull,:),0.)

!       following is useful if troublesome data is read in
#ifdef debug
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
#endif
        call retopo(pslb,zsb,zs(1:ifull),tb,qb)
#ifdef debug
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
#endif
      
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
#ifdef debug
        if (mydiag) then
          write (6,"('zsb# nestinb  ',9f7.1)") diagvals(zsb)
          write (6,"('tssb# nestinb ',9f7.1)") diagvals(tssb) 
        end if
#endif

        ! calculate time for next filter call   
        kdhour=ktime_r/100-ktime/100   ! integer hour diff from Oct '05
        kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
        mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     &                +60*kdhour+kdmin
#ifdef debug
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
#endif
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
        else
          if (myid==0) then
#ifdef uniform_decomp
            write(6,*) "Separable 1D filter (MPI)            ",kb
#else
            write(6,*) "Separable 1D filter (MPI optimised)  ",kb
#endif
          end if
          call specfastmpi(.1*real(mbd)/(pi*schmidt)
     &                  ,psld(:),ud(:,1:klt),vd(:,1:klt)
     &                  ,td(:,1:klt),qd(:,1:klt),lblock,klt)
        endif  ! (nud_uv==9) .. else ..
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
      use vecsuv_m          ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration

      integer, intent(in) :: klt
      integer i,j,n,iq,iqg,k
      real, intent(in) :: cin
      real, dimension(ifull), intent(inout) :: pslb
      real, dimension(ifull,klt), intent(inout) :: ub,vb
      real, dimension(ifull,klt), intent(inout) :: tb,qb
      real, dimension(ifull,klt) :: wb
      real, dimension(ifull_g,klt) :: tt
      real, dimension(ifull) :: da,db
      real, dimension(klt) :: ud,vd,wd
      real cq
      logical, intent(in) :: lblock

      cq=sqrt(4.5)*cin

      ! Create global map on each processor.  This can require a lot of memory
      if (nud_p>0.and.lblock) then
        call ccmpi_gatherall(pslb(:), tt(:,1))
        call slowspecmpi_work(cin,tt(:,1),pslb,1)
      end if
      if (nud_uv==3) then
        call ccmpi_gatherall(ub(:,1:klt),tt(:,1:klt))
        call slowspecmpi_work(cin,tt,ub,klt)
      else if (nud_uv>0) then
        do k=1,klt
          da=ub(:,k)
          db=vb(:,k)
          ub(:,k)=ax(1:ifull)*da+bx(1:ifull)*db
          vb(:,k)=ay(1:ifull)*da+by(1:ifull)*db
          wb(:,k)=az(1:ifull)*da+bz(1:ifull)*db
        end do
        call ccmpi_gatherall(ub(:,1:klt),tt(:,1:klt))
        call slowspecmpi_work(cin,tt,ub,klt)
        call ccmpi_gatherall(vb(:,1:klt),tt(:,1:klt))
        call slowspecmpi_work(cin,tt,vb,klt)
        call ccmpi_gatherall(wb(:,1:klt),tt(:,1:klt))
        call slowspecmpi_work(cin,tt,wb,klt)
      endif
      if (nud_t>0) then
        call ccmpi_gatherall(tb(:,1:klt),tt(:,1:klt))
        call slowspecmpi_work(cin,tt,tb,klt)
      end if
      if (nud_q>0) then
        call ccmpi_gatherall(qb(:,1:klt),tt(:,1:klt))
        call slowspecmpi_work(cin,tt,qb,klt)
      end if

      return
      end subroutine slowspecmpi

      subroutine slowspecmpi_work(cin,tt,tb,klt)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      use xyzinfo_m         ! Grid coordinate arrays
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters

      integer, intent(in) :: klt
      integer i,j,n,iq,iqg,k
      real, intent(in) :: cin
      real, dimension(ifull,klt), intent(out) :: tb
      real, dimension(ifull_g,klt), intent(in) :: tt
      real, dimension(ifull_g) :: r
      real cq,psum

#ifdef debug
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "Start 2D filter"
      end if
#endif
   
      do n=1,npan
        do j=1,jpan
          do i=1,ipan
            iqg=indg(i,j,n)
            iq =indp(i,j,n)
            r(:)=x_g(iqg)*x_g(:)+y_g(iqg)*y_g(:)+z_g(iqg)*z_g(:)
            r(:)=acos(max(min(r(:),1.),-1.))
            r(:)=exp(-(cq*r(:))**2)/(em_g(:)*em_g(:))
            psum=sum(r(:))
            do k=1,klt
              tb(iq,k)=sum(r(:)*tt(:,k))/psum
            end do
          end do
        end do
      end do
 
#ifdef debug
      if (myid == 0.and.nmaxpr==1) write(6,*) "End 2D filter"
#endif

      return
      end subroutine slowspecmpi_work
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      ! Four pass spectral downscaling
      subroutine specfastmpi(cin,psls,uu,vv,tt,qgg,lblock,klt)
      
      use cc_mpi             ! CC MPI routines
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: klt
      integer ifg,xpan
      real, intent(in) :: cin
      real, dimension(ifull), intent(inout) :: psls
      real, dimension(ifull,klt), intent(inout) :: uu,vv
      real, dimension(ifull,klt), intent(inout) :: tt,qgg
      logical, intent(in) :: lblock
      
      if (joff==0) then
        ifg=ifull_g
      else
        ifg=0
      end if
      xpan=max(ipan,jpan)
      if (npta==1) then
        ! reduced memory version
        call spechost_n(cin,psls,uu,vv,tt,qgg,lblock,klt,ifg,xpan)
      else
        ! normal version
        call spechost(cin,psls,uu,vv,tt,qgg,lblock,klt,ifg,xpan)
      end if

      return
      end subroutine specfastmpi
      !---------------------------------------------------------------------------------


      !---------------------------------------------------------------------------------
      ! This is the main routine for the scale-selective filter
      ! (see spechost_n for a reduced memory version)
      subroutine spechost(cin,pslb,ub,vb,tb,qb,lblock,klt,ifg,xpan)

      use cc_mpi            ! CC MPI routines
      use vecsuv_m          ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: klt,ifg,xpan
      integer k,ppass,nne,nns,oos,ooe,iproc
      integer j,n,til,ipoff,jpoff,npoff
      real, intent(in) :: cin
      real, dimension(ifull), intent(inout) :: pslb
      real, dimension(ifull,klt), intent(inout) :: ub,vb
      real, dimension(ifull,klt), intent(inout) :: tb,qb
      real, dimension(ifull,klt) :: wb
      real, dimension(ifg,klt) :: uu,vv,ww
      real, dimension(ifg,klt) :: tt,qgg
      real, dimension(ifg) :: psls
      real, dimension(ifg) :: qp
      real, dimension(ifull) :: xa_l,xb_l
      real, dimension(ifg,klt) :: qu,qv,qw,qt,qq
      real, dimension(ifg*klt) :: dd,gg
      real, dimension(ifull*klt) :: ff
      real cq
      logical, intent(in) :: lblock

      til=il_g*il_g
      cq=sqrt(4.5)*cin ! filter length scale

      ! gather data onto myid==0
      if (myid==0) then
        if (nud_p>0.and.lblock) then
          call ccmpi_gather(pslb(:), psls(:))
        end if
        if (nud_uv==3) then
          call ccmpi_gather(ub(:,1:klt),ww(:,1:klt))
        else if (nud_uv>0) then
          do k=1,klt
            xa_l=ub(:,k)
            xb_l=vb(:,k)
            ub(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
            vb(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
            wb(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
          end do
          call ccmpi_gather(ub(:,1:klt),uu(:,1:klt))
          call ccmpi_gather(vb(:,1:klt),vv(:,1:klt))
          call ccmpi_gather(wb(:,1:klt),ww(:,1:klt))
        end if
        if (nud_t>0) then
          call ccmpi_gather(tb(:,1:klt),tt(:,1:klt))
        end if
        if (nud_q>0) then
          call ccmpi_gather(qb(:,1:klt),qgg(:,1:klt))
        end if
      else
        if (nud_p>0.and.lblock) then
          call ccmpi_gather(pslb(:))
        end if
        if (nud_uv==3) then
          call ccmpi_gather(ub(:,1:klt))
        else if (nud_uv>0) then
          do k=1,klt
            xa_l=ub(:,k)
            xb_l=vb(:,k)
            ub(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
            vb(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
            wb(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
          end do
          call ccmpi_gather(ub(:,1:klt))
          call ccmpi_gather(vb(:,1:klt))
          call ccmpi_gather(wb(:,1:klt))
        end if
        if (nud_t>0) then
          call ccmpi_gather(tb(:,1:klt))
        end if          
        if (nud_q>0) then
          call ccmpi_gather(qb(:,1:klt))
        end if
      end if
      
#ifdef debug
      if (myid==0) write(6,*) "Start 1D filter"
#endif

      ! distribute data over host processors
      if (myid==hproc) then
        if (nud_p>0.and.lblock) then
          call ccmpi_bcast(psls,0,comm_host)
        end if
        if (nud_uv>0) then
          call ccmpi_bcast(uu,0,comm_host)
          call ccmpi_bcast(vv,0,comm_host)
          call ccmpi_bcast(ww,0,comm_host)
        end if
        if (nud_t>0) then
          call ccmpi_bcast(tt,0,comm_host)
        end if
        if (nud_q>0) then
          call ccmpi_bcast(qgg,0,comm_host)
        end if
      end if
      
      do ppass=pprocn,pprocx

        ! reset nudging fields for the next panel
        if (myid==hproc) then
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
        call speclocal(cq,ppass,qp,qu,qv,qw,qt,qq,lblock,klt,
     &         ifg,xpan)
      
        ! distribute results  
        if(nud_p>0.and.lblock)then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do j=nns,nne
                do n=oos,ooe
                  dd(n-oos+1+ipan*(j-nns)+ipan*jpan*(iproc-hproc))=
     &              qp(n+il_g*(j-1)+ppass*til)
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:til),ff(1:ipan*jpan),0,comm_proc)
          do n=1,ipan*jpan
            pslb(n+ipan*jpan*(ppass-pprocn))=ff(n)
          end do
        end if
        if (nud_uv==3) then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,klt
                do j=nns,nne
                  do n=oos,ooe
                    dd(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*klt*(iproc-hproc))=qw(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &           comm_proc)
          do k=1,klt
            do n=1,ipan*jpan
              ub(n+ipan*jpan*(ppass-pprocn),k)=ff(n+ipan*jpan*(k-1))
            end do
          end do
        else if(nud_uv>0)then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,klt
                do j=nns,nne
                  do n=oos,ooe
                    dd(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*klt*(iproc-hproc))=qu(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &           comm_proc)
          do k=1,klt
            do n=1,ipan*jpan
              ub(n+ipan*jpan*(ppass-pprocn),k)=ff(n+ipan*jpan*(k-1))
            end do
          end do
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,klt
                do j=nns,nne
                  do n=oos,ooe
                    dd(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*klt*(iproc-hproc))=qv(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &           comm_proc)
          do k=1,klt
            do n=1,ipan*jpan
              vb(n+ipan*jpan*(ppass-pprocn),k)=ff(n+ipan*jpan*(k-1))
            end do
          end do
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,klt
                do j=nns,nne
                  do n=oos,ooe
                    dd(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*klt*(iproc-hproc))=qw(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &           comm_proc)
          do k=1,klt
            do n=1,ipan*jpan
              wb(n+ipan*jpan*(ppass-pprocn),k)=ff(n+ipan*jpan*(k-1))
            end do
          end do
          do k=1,klt
            do n=1+ipan*jpan*(ppass-pprocn),ipan*jpan*(ppass-pprocn+1)
              xa_l(n)=ax(n)*ub(n,k)+ay(n)*vb(n,k)
     &               +az(n)*wb(n,k)
              xb_l(n)=bx(n)*ub(n,k)+by(n)*vb(n,k)
     &               +bz(n)*wb(n,k)
              ub(n,k)=xa_l(n)
              vb(n,k)=xb_l(n)
            end do
          end do
        end if
        if(nud_t>0)then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,klt
                do j=nns,nne
                  do n=oos,ooe
                    dd(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*klt*(iproc-hproc))=qt(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &           comm_proc)
          do k=1,klt
            do n=1,ipan*jpan
              tb(n+ipan*jpan*(ppass-pprocn),k)=ff(n+ipan*jpan*(k-1))
            end do
          end do
        end if
        if(nud_q>0)then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,klt
                do j=nns,nne
                  do n=oos,ooe
                    dd(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*klt*(iproc-hproc))=qq(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &           comm_proc)
          do k=1,klt
            do n=1,ipan*jpan
              qb(n+ipan*jpan*(ppass-pprocn),k)=ff(n+ipan*jpan*(k-1))
            end do
          end do
        end if
        
      end do

#ifdef debug
      if (myid==0) then
        write(6,*) "End 1D filter"
      end if
#endif

      return
      end subroutine spechost
      !---------------------------------------------------------------------------------

      ! This version is for one panel per processor (reduced memory)
      subroutine spechost_n(cin,pslb,ub,vb,tb,qb,lblock,klt,ifg,xpan)

      use cc_mpi            ! CC MPI routines
      use vecsuv_m          ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: klt,ifg,xpan
      integer k,iproc,nns,nne,oos,ooe
      integer ipoff,jpoff,npoff
      integer n,j,til
      real, intent(in) :: cin
      real cq
      real, dimension(ifull), intent(inout) :: pslb
      real, dimension(ifull,klt), intent(inout) :: ub,vb
      real, dimension(ifull,klt), intent(inout) :: tb,qb
      real, dimension(ifull,klt) :: wb
      real, dimension(ifg,klt) :: uu,vv,ww
      real, dimension(ifg,klt) :: tt,qgg
      real, dimension(ifg*klt) :: dd,gg
      real, dimension(ifull*klt) :: ff
      real, dimension(ifg) :: psls
      real, dimension(ifull) :: xa_l,xb_l
      logical, intent(in) :: lblock

      til=il_g*il_g
      cq=sqrt(4.5)*cin ! filter length scale

      ! gather data
      if (nud_p>0.and.lblock) then
        call ccmpi_gatherx(dd(1:til),pslb(1:ifull),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do j=nns,nne
              do n=oos,ooe
                gg(n+il_g*(j-1))=dd(n-oos+1+ipan*(j-nns)
     &            +ifull*(iproc-hproc))
              end do
            end do
          end do
          call ccmpi_allgatherx(psls(1:ifull_g),gg(1:til),comm_host)
        end if
      end if
      if (nud_uv==3) then
        do k=1,klt
          do n=1,ifull
            ff(n+ifull*(k-1))=ub(n,k)
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  gg(n+il_g*(j-1)+til*(k-1))=dd(n-oos+1+ipan*(j-nns)
     &              +ifull*(k-1)+ifull*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(dd(1:ifull_g*klt),gg(1:til*klt),
     &           comm_host)
          do iproc=0,5
            do k=1,klt
              do n=1,til
                ww(n+til*iproc,k)=dd(n+til*(k-1)+til*klt*iproc)
              end do
            end do
          end do
        end if
      else if (nud_uv>0) then
        do k=1,klt
          xa_l=ub(:,k)
          xb_l=vb(:,k)
          ub(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
          vb(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
          wb(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
        end do
        do k=1,klt
          do n=1,ifull
            ff(n+ifull*(k-1))=ub(n,k)
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  gg(n+il_g*(j-1)+til*(k-1))=dd(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(dd(1:ifull_g*klt),gg(1:til*klt),
     &           comm_host)
          do iproc=0,5
            do k=1,klt
              do n=1,til
                uu(n+til*iproc,k)=dd(n+til*(k-1)+til*klt*iproc)
              end do
            end do
          end do
        end if
        do k=1,klt
          do n=1,ifull
            ff(n+ifull*(k-1))=vb(n,k)
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  gg(n+il_g*(j-1)+til*(k-1))=dd(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(dd(1:ifull_g*klt),gg(1:til*klt),
     &           comm_host)
          do iproc=0,5
            do k=1,klt
              do n=1,til
                vv(n+til*iproc,k)=dd(n+til*(k-1)+til*klt*iproc)
              end do
            end do
          end do
        end if
        do k=1,klt
          do n=1,ifull
            ff(n+ifull*(k-1))=wb(n,k)
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  gg(n+il_g*(j-1)+til*(k-1))=dd(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(dd(1:ifull_g*klt),gg(1:til*klt),
     &           comm_host)
          do iproc=0,5
            do k=1,klt
              do n=1,til
                ww(n+til*iproc,k)=dd(n+til*(k-1)+til*klt*iproc)
              end do
            end do
          end do
        end if
      end if
      if (nud_t>0) then
        do k=1,klt
          do n=1,ifull
            ff(n+ifull*(k-1))=tb(n,k)
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  gg(n+il_g*(j-1)+til*(k-1))=dd(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(dd(1:ifull_g*klt),gg(1:til*klt),
     &           comm_host)
          do iproc=0,5
            do k=1,klt
              do n=1,til
                tt(n+til*iproc,k)=dd(n+til*(k-1)+til*klt*iproc)
              end do
            end do
          end do
        end if
      end if
      if (nud_q>0) then
        do k=1,klt
          do n=1,ifull
            ff(n+ifull*(k-1))=qb(n,k)
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  gg(n+il_g*(j-1)+til*(k-1))=dd(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(dd(1:ifull_g*klt),gg(1:til*klt),
     &           comm_host)
          do iproc=0,5
            do k=1,klt
              do n=1,til
                qgg(n+til*iproc,k)=dd(n+til*(k-1)+til*klt*iproc)
              end do
            end do
          end do
        end if
      end if

#ifdef debug
      if (myid==0) write(6,*) "Start 1D filter"
#endif

      ! computations for the local processor group
      call speclocal(cq,pprocn,psls,uu,vv,ww,tt,qgg,lblock,klt,ifg,
     &               xpan)

      ! distribute data to processors
      if(nud_p>0.and.lblock)then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do j=nns,nne
              do n=oos,ooe
                dd(n-oos+1+ipan*(j-nns)+ifull*(iproc-hproc))=
     &            psls(n+il_g*(j-1)+pprocn*til)
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(dd(1:til),pslb(1:ifull),0,comm_proc)
      end if
      if (nud_uv==3) then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  dd(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))=ww(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        do k=1,klt
          do n=1,ifull
            ub(n,k)=ff(n+ifull*(k-1))
          end do
        end do
      else if(nud_uv>0)then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  dd(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))=uu(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        do k=1,klt
          do n=1,ifull
            ub(n,k)=ff(n+ifull*(k-1))
          end do
        end do
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  dd(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))=vv(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        do k=1,klt
          do n=1,ifull
            vb(n,k)=ff(n+ifull*(k-1))
          end do
        end do
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  dd(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))=ww(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        do k=1,klt
          do n=1,ifull
            wb(n,k)=ff(n+ifull*(k-1))
          end do
        end do
        do k=1,klt
          xa_l=ax(1:ifull)*ub(:,k)+ay(1:ifull)*vb(:,k)
     &        +az(1:ifull)*wb(:,k)
          xb_l=bx(1:ifull)*ub(:,k)+by(1:ifull)*vb(:,k)
     &        +bz(1:ifull)*wb(:,k)
          ub(:,k)=xa_l
          vb(:,k)=xb_l
        end do
      end if
      if(nud_t>0)then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  dd(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))=tt(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        do k=1,klt
          do n=1,ifull
            tb(n,k)=ff(n+ifull*(k-1))
          end do
        end do
      end if
      if(nud_q>0)then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  dd(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*klt*(iproc-hproc))=qgg(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(dd(1:til*klt),ff(1:ifull*klt),0,comm_proc)
        do k=1,klt
          do n=1,ifull
            qb(n,k)=ff(n+ifull*(k-1))
          end do
        end do
      end if
      
      return
      end subroutine spechost_n
      
      !---------------------------------------------------------------------------------
      ! This code runs between the local processors
      
      ! hproc is the captian processor that farms out the convolution to the comm_host
      ! group of processors.
      subroutine speclocal(cq,ppass,qp,qu,qv,qw,qt,
     &                     qq,lblock,klt,ifg,xpan)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      use xyzinfo_m         ! Grid coordinate arrays

      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: ppass,klt,ifg,xpan
      integer j,k,n,ipass,iproc
      integer ipoff,jpoff,npoff
      integer nne,nns,ooe,oos,me,ns,ne,os,oe
      integer til
      integer a,b,c,sn,sy,jj,nn
      integer, dimension(0:3) :: astr,bstr,cstr
      integer, dimension(0:3) :: maps
      real, intent(in) :: cq
      real, dimension(ifg,klt), intent(inout) :: qu,qv,qw
      real, dimension(ifg,klt), intent(inout) :: qt,qq
      real, dimension(ifg), intent(inout) :: qp
      real, dimension(ifg) :: qsum
      real, dimension(4*il_g,xpan,klt) :: au,av,aw,at,aq
      real, dimension(4*il_g,xpan) :: ap,asum
      real, dimension(xpan,xpan,klt) :: pu,pv
      real, dimension(xpan,xpan,klt) :: pw
      real, dimension(xpan,xpan,klt) :: pt,pq
      real, dimension(xpan,xpan) :: pp,psum
      real, dimension(4*il_g) :: ra,xa,ya,za
      real, dimension(4*il_g*il_g*klt) :: dd
      real, dimension(4*il_g*xpan*klt) :: ff
      logical, intent(in) :: lblock
      
      ! panels 1,2 and 3 are matched so that the output is the same
      ! as the ax array contents on ipass==3.  panels 0, 4 and 5
      ! need ipan and jpan to be switched in the convolutions.
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      til=il_g*il_g
      astr=0
      bstr=0
      cstr=0

      ns=joff+1
      ne=joff+jpan
      os=ioff+1
      oe=ioff+ipan
      
      do ipass=0,2
        me=maps(ipass)
        call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug        
        if (myid==0) then
          write(6,*) "Receive arrays from local host"
        end if
#endif

        do j=1,jpan
          jj=j+ns-1
          do sn=1,me,il_g
            sy=(sn-1)/il_g
            a=astr(sy)
            b=bstr(sy)
            c=cstr(sy)
            do n=sn,sn+il_g-1
              asum(n,j)=1./em_g(a*n+b*jj+c)**2
            end do
          end do
        end do

        ! hproc should be ioff=0 and joff=0
        if(nud_p>0.and.lblock)then
          if (ioff==0) then
            if (myid==hproc) then
              do j=1,il_g
                do sn=1,me,il_g
                  sy=(sn-1)/il_g
                  a=astr(sy)
                  b=bstr(sy)
                  c=cstr(sy)
                  do n=sn,sn+il_g-1
                    dd(n+me*(j-1))=qp(a*n+b*j+c)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:me*il_g),ff(1:me*jpan),0,
     &             comm_cols)
          end if
          call ccmpi_bcast(ff(1:me*jpan),0,comm_rows)
          do j=1,jpan
            do n=1,me
              ap(n,j)=ff(n+me*(j-1))*asum(n,j)
            end do
          end do
        end if
        if(nud_uv>0)then
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,klt
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      dd(n+me*(k-1)+me*klt*(j-1))=qu(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:me*il_g*klt),ff(1:me*jpan*klt),
     &             0,comm_cols)
          end if
          call ccmpi_bcast(ff(1:me*jpan*klt),0,comm_rows)
          do k=1,klt
            do j=1,jpan
              do n=1,me
                au(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))*asum(n,j)
              end do
            end do
          end do
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,klt
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      dd(n+me*(k-1)+me*klt*(j-1))=qv(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:me*il_g*klt),ff(1:me*jpan*klt),
     &             0,comm_cols)
          end if
          call ccmpi_bcast(ff(1:me*jpan*klt),0,comm_rows)
          do k=1,klt
            do j=1,jpan
              do n=1,me
                av(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))*asum(n,j)
              end do
            end do
          end do
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,klt
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      dd(n+me*(k-1)+me*klt*(j-1))=qw(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:me*il_g*klt),ff(1:me*jpan*klt),
     &             0,comm_cols)
          end if
          call ccmpi_bcast(ff(1:me*jpan*klt),0,comm_rows)
          do k=1,klt
            do j=1,jpan
              do n=1,me
                aw(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))*asum(n,j)
              end do
            end do
          end do
        end if
        if(nud_t>0)then
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,klt
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      dd(n+me*(k-1)+me*klt*(j-1))=qt(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:me*il_g*klt),ff(1:me*jpan*klt),
     &             0,comm_cols)
          end if
          call ccmpi_bcast(ff(1:me*jpan*klt),0,comm_rows)
          do k=1,klt
            do j=1,jpan
              do n=1,me
                at(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))*asum(n,j)
              end do
            end do
          end do
        end if
        if(nud_q>0)then
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,klt
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      dd(n+me*(k-1)+me*klt*(j-1))=qq(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:me*il_g*klt),ff(1:me*jpan*klt),
     &             0,comm_cols)
          end if
          call ccmpi_bcast(ff(1:me*jpan*klt),0,comm_rows)
          do k=1,klt
            do j=1,jpan
              do n=1,me
                aq(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))*asum(n,j)
              end do
            end do
          end do
        end if

#ifdef debug
        if (myid==0) write(6,*) "Start convolution"
#endif

        do j=1,jpan
          jj=j+ns-1
          do sn=1,me,il_g
            sy=(sn-1)/il_g
            a=astr(sy)
            b=bstr(sy)
            c=cstr(sy)
            do n=sn,sn+il_g-1
              xa(n)=x_g(a*n+b*jj+c)
              ya(n)=y_g(a*n+b*jj+c)
              za(n)=z_g(a*n+b*jj+c)
            end do
          end do
          do n=1,ipan
            nn=n+os-1
            ra(1:me)=xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me)
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

#ifdef debug
        if (myid==0) then
          write(6,*) "End convolution"
          write(6,*) "Send arrays to local host"
        end if
#endif

        ! unpacking grid
        a=astr(0)
        b=bstr(0)
        c=cstr(0)

        ! gather data on host processors
        ! gather psum
        do j=1,jpan
          do n=1,ipan
            ff(n+ipan*(j-1))=psum(n,j)
          end do
        end do
        call ccmpi_gatherx(dd(1:il_g*ipan),ff(1:ipan*jpan),0,comm_cols)
        if (joff==0) then
          do jpoff=0,il_g-1,jpan
            sy=jpoff/jpan
            nns=jpoff+1
            nne=jpoff+jpan
            do j=nns,nne
              do n=os,oe
                qsum(a*n+b*j+c)=dd(n-os+1+ipan*(j-nns)
     &                            +ipan*jpan*sy)
              end do
            end do
          end do
        end if
        if (nud_p>0.and.lblock) then
          do j=1,jpan
            do n=1,ipan
              ff(n+jpan*(j-1))=pp(n,j)
            end do
          end do
          call ccmpi_gatherx(dd(1:il_g*ipan),ff(1:ipan*jpan),0,
     &           comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do j=nns,nne
                do n=os,oe
                  qp(a*n+b*j+c)=dd(n-os+1+jpan*(j-nns)
     &                            +ipan*jpan*sy)
                end do
              end do
            end do
          end if
        end if
        if (nud_uv>0) then
          do k=1,klt
            do j=1,jpan
              do n=1,ipan
                ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pu(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(dd(1:il_g*ipan*klt),ff(1:ipan*jpan*klt),
     &           0,comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,klt
                do j=nns,nne
                  do n=os,oe
                    qu(a*n+b*j+c,k)=dd(n-os+1+jpan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*klt*sy)
                  end do
                end do
              end do
            end do
          end if
          do k=1,klt
            do j=1,jpan
              do n=1,ipan
                ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pv(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(dd(1:il_g*ipan*klt),ff(1:ipan*jpan*klt),
     &           0,comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,klt
                do j=nns,nne
                  do n=os,oe
                    qv(a*n+b*j+c,k)=dd(n-os+1+jpan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*klt*sy)
                  end do
                end do
              end do
            end do
          end if
          do k=1,klt
            do j=1,jpan
              do n=1,ipan
                ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pw(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(dd(1:il_g*ipan*klt),ff(1:ipan*jpan*klt),
     &           0,comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,klt
                do j=nns,nne
                  do n=os,oe
                    qw(a*n+b*j+c,k)=dd(n-os+1+jpan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*klt*sy)
                  end do
                end do
              end do
            end do
          end if
        end if
        if (nud_t>0) then
          do k=1,klt
            do j=1,jpan
              do n=1,ipan
                ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pt(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(dd(1:il_g*ipan*klt),ff(1:ipan*jpan*klt),
     &           0,comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,klt
                do j=nns,nne
                  do n=os,oe
                    qt(a*n+b*j+c,k)=dd(n-os+1+jpan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*klt*sy)
                  end do
                end do
              end do
            end do
          end if
        end if
        if (nud_q>0) then
          do k=1,klt
            do j=1,jpan
              do n=1,ipan
                ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pq(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(dd(1:il_g*ipan*klt),ff(1:ipan*jpan*klt),
     &           0,comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,klt
                do j=nns,nne
                  do n=os,oe
                    qq(a*n+b*j+c,k)=dd(n-os+1+jpan*(j-nns)
     &                 +ipan*jpan*(k-1)+ipan*jpan*klt*sy)
                  end do
                end do
              end do
            end do
          end if
        end if

      end do

      ns=ioff+1
      ne=ioff+ipan
      os=joff+1
      oe=joff+jpan

      ipass=3
      me=maps(ipass)
      call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug        
      if (myid==0) then
        write(6,*) "Receive arrays from local host"
      end if
#endif

      if (joff==0) then
        do j=ns,ne
          do sn=1,me,il_g
            sy=(sn-1)/il_g
            a=astr(sy)
            b=bstr(sy)
            c=cstr(sy)
            do n=sn,sn+il_g-1
              ff(n+me*(j-ns))=qsum(a*n+b*j+c)
            end do
          end do
        end do
      end if
      call ccmpi_bcast(ff(1:me*ipan),0,comm_cols)
      do j=1,ipan
        do n=1,me
          asum(n,j)=ff(n+me*(j-1))
        end do
      end do
      if(nud_p>0.and.lblock)then
        if (joff==0) then
          do j=ns,ne
            do sn=1,me,il_g
              sy=(sn-1)/il_g
              a=astr(sy)
              b=bstr(sy)
              c=cstr(sy)
              do n=sn,sn+il_g-1
                ff(n+me*(j-ns))=qp(a*n+b*j+c)
              end do
            end do
          end do
        end if
        call ccmpi_bcast(ff(1:me*ipan),0,comm_cols)
        do j=1,ipan
          do n=1,me
            ap(n,j)=ff(n+me*(j-1))
          end do
        end do
      end if
      if(nud_uv>0)then
        if (joff==0) then
          do k=1,klt
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  ff(n+me*(k-1)+me*klt*(j-ns))=qu(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(ff(1:me*ipan*klt),0,comm_cols)
        do k=1,klt
          do j=1,ipan
            do n=1,me
              au(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))
            end do
          end do
        end do
        if (joff==0) then
          do k=1,klt
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  ff(n+me*(k-1)+me*klt*(j-ns))=qv(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(ff(1:me*ipan*klt),0,comm_cols)
        do k=1,klt
          do j=1,ipan
            do n=1,me
              av(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))
            end do
          end do
        end do
        if (joff==0) then
          do k=1,klt
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  ff(n+me*(k-1)+me*klt*(j-ns))=qw(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(ff(1:me*ipan*klt),0,comm_cols)
        do k=1,klt
          do j=1,ipan
            do n=1,me
              aw(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))
            end do
          end do
        end do
      end if
      if(nud_t>0)then
        if (joff==0) then
          do k=1,klt
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  ff(n+me*(k-1)+me*klt*(j-ns))=qt(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(ff(1:me*ipan*klt),0,comm_cols)
        do k=1,klt
          do j=1,ipan
            do n=1,me
              at(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))
            end do
          end do
        end do
      end if
      if(nud_q>0)then
        if (joff==0) then
          do k=1,klt
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  ff(n+me*(k-1)+me*klt*(j-ns))=qq(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(ff(1:me*ipan*klt),0,comm_cols)
        do k=1,klt
          do j=1,ipan
            do n=1,me
              aq(n,j,k)=ff(n+me*(k-1)+me*klt*(j-1))
            end do
          end do
        end do
      end if
 
#ifdef debug
      if (myid==0) write(6,*) "Start convolution"
#endif

      do j=1,ipan
        jj=j+ns-1
        do sn=1,me,il_g
          sy=(sn-1)/il_g
          a=astr(sy)
          b=bstr(sy)
          c=cstr(sy)
          do n=sn,sn+il_g-1
            xa(n)=x_g(a*n+b*jj+c)
            ya(n)=y_g(a*n+b*jj+c)
            za(n)=z_g(a*n+b*jj+c)
          end do
        end do
        do n=1,jpan
          nn=n+os-1
          ra(1:me)=xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me)
          ra(1:me)=acos(max(min(ra(1:me),1.),-1.))
          ra(1:me)=exp(-(cq*ra(1:me))**2)
          ! can also use the lines below which integrate the gaussian
          ! analytically over the length element (but slower)
          !ra(1)=2.*erf(cq*0.5*(ds/rearth)
          !ra(2:me)=erf(cq*(ra(2:me)+0.5*(ds/rearth)))  ! redefine ra(:) as wgt(:)
     &    !        -erf(cq*(ra(2:me)-0.5*(ds/rearth)))  ! (correct units are 1/cq)
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

#ifdef debug
      if (myid==0) then
        write(6,*) "End convolution"
        write(6,*) "Send arrays to local host"
      end if
#endif

      ! unpacking grid
      a=astr(0)
      b=bstr(0)
      c=cstr(0)

      ! gather data on host processors
      if (nud_p>0.and.lblock) then
        do j=1,ipan
          do n=1,jpan
            ff(n+jpan*(j-1))=pp(n,j)/psum(n,j)
          end do
        end do
        call ccmpi_gatherx(dd(1:til),ff(1:ipan*jpan),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do j=nns,nne
              do n=oos,ooe
                qp(a*n+b*j+c)=dd(n-oos+1+jpan*(j-nns)
     &                          +ipan*jpan*(iproc-hproc))
              end do
            end do
          end do
        end if
      end if
      if (nud_uv>0) then
        do k=1,klt
          do j=1,ipan
            do n=1,jpan
              ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pu(n,j,k)/psum(n,j)
            end do
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  qu(a*n+b*j+c,k)=dd(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
        do k=1,klt
          do j=1,ipan
            do n=1,jpan
              ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pv(n,j,k)/psum(n,j)
            end do
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  qv(a*n+b*j+c,k)=dd(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
        do k=1,klt
          do j=1,ipan
            do n=1,jpan
              ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pw(n,j,k)/psum(n,j)
            end do
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  qw(a*n+b*j+c,k)=dd(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
      end if
      if (nud_t>0) then
        do k=1,klt
          do j=1,ipan
            do n=1,jpan
              ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pt(n,j,k)/psum(n,j)
            end do
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  qt(a*n+b*j+c,k)=dd(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
      end if
      if (nud_q>0) then
        do k=1,klt
          do j=1,ipan
            do n=1,jpan
              ff(n+jpan*(j-1)+ipan*jpan*(k-1))=pq(n,j,k)/psum(n,j)
            end do
          end do
        end do
        call ccmpi_gatherx(dd(1:til*klt),ff(1:ipan*jpan*klt),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,klt
              do j=nns,nne
                do n=oos,ooe
                  qq(a*n+b*j+c,k)=dd(n-oos+1+jpan*(j-nns)
     &               +ipan*jpan*(k-1)+ipan*jpan*klt*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
      end if

      return  
      end subroutine speclocal

      !---------------------------------------------------------------------------------
      ! Map from 1D convolution to global index
      subroutine getiqa(a,b,c,ne,ipass,ppass,il_g)
      
      implicit none
      
      integer, intent(in) :: ne,ipass,ppass,il_g
      integer, dimension(0:3), intent(out) :: a,b,c
      integer sn,sy
      
      do sn=1,ne,il_g
        sy=(sn-1)/il_g

        select case(ppass*100+ipass*10+sy)
          case(0)                                ! panel 5   - x pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=il_g*(5*il_g-1)
          case(10)                               ! panel 2   - x pass
            a(sy)=-1
            b(sy)=il_g
            c(sy)=2*il_g*il_g+1
          case(20,21,321)                        ! panel 0,1 - y pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=-il_g
          case(22)                               ! panel 3   - y pass
            a(sy)=1
            b(sy)=-il_g
            c(sy)=il_g*(4*il_g-2)
          case(23,323)                           ! panel 4   - y pass
            a(sy)=1
            b(sy)=-il_g
            c(sy)=il_g*(5*il_g-3)
          case(30,100)                           ! panel 0   - z pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=-il_g
          case(31,223,523,532)                   ! panel 2   - z pass ! panel 4   - x pass ! panel 3   - z pass
            a(sy)=il_g
            b(sy)=-1
            c(sy)=il_g*il_g+1
          case(32)                               ! panel 5   - z pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(5*il_g-3)
          case(110)                              ! panel 3   - z pass
            a(sy)=-il_g
            b(sy)=1
            c(sy)=4*il_g*il_g
          case(120)                              ! panel 1   - x pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(il_g-1)
          case(121,421)                          ! panel 2   - x pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=2*il_g*(il_g-1)
          case(122,123,423)                      ! panel 4,5 - x pass ! panel 2   - z pass
            a(sy)=il_g
            b(sy)=-1
            c(sy)=2*il_g*il_g+1
          case(130)                              ! panel 1   - y pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=il_g*(il_g-1)
          case(131)                              ! panel 3   - y pass
            a(sy)=1
            b(sy)=-il_g
            c(sy)=il_g*(4*il_g-1)
          case(132,322)                          ! panel 0,1 - y pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=-il_g*(2*il_g+1)
          case(200)                              ! panel 0   - y pass
            a(sy)=-il_g
            b(sy)=1
            c(sy)=il_g*il_g
          case(210)                              ! panel 3   - y pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(3*il_g-1)
          case(220)                              ! panel 2   - x pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(2*il_g-1)
          case(221,521)                          ! panel 1   - x pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(il_g-2)
          case(222,410)                          ! panel 3   - z pass ! panel 5   - x pass
            a(sy)=il_g
            b(sy)=-1
            c(sy)=3*il_g*il_g+1
          case(230)                              ! panel 2   - z pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=il_g*(2*il_g-1)
          case(231)                              ! panel 0   - z pass
            a(sy)=1
            b(sy)=-il_g
            c(sy)=il_g*(il_g-1)
          case(232)                              ! panel 3   - z pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=il_g*(il_g-1)
          case(300)                              ! panel 5   - x pass
            a(sy)=-il_g
            b(sy)=-1
            c(sy)=il_g*(6*il_g+1)+1
          case(310)                              ! panel 2   - x pass
            a(sy)=1
            b(sy)=-il_g
            c(sy)=3*il_g*il_g
          case(320)                              ! panel 3   - y pass
            a(sy)=1
            b(sy)=-il_g
            c(sy)=4*il_g*il_g
          case(330)                              ! panel 3   - z pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=il_g*(3*il_g-1)
          case(331)                              ! panel 2   - z pass
            a(sy)=il_g
            b(sy)=1
            c(sy)=il_g*(il_g-1)
          case(332)                              ! panel 5   - z pass
            a(sy)=1
            b(sy)=-il_g
            c(sy)=il_g*(6*il_g-2)
          case(400)                              ! panel 0   - z pass
            a(sy)=-1
            b(sy)=-il_g
            c(sy)=il_g*(il_g+1)+1
          case(420)                              ! panel 4   - x pass
            a(sy)=il_g
            b(sy)=-1
            c(sy)=4*il_g*il_g+1
          case(422)                              ! panel 1   - x pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(il_g-3)
          case(430)                              ! panel 4   - y pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(4*il_g-1)
          case(431)                              ! panel 3   - y pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(3*il_g-2)
          case(432)                              ! panel 0   - y pass
            a(sy)=il_g
            b(sy)=-1
            c(sy)=-2*il_g*il_g+1
          case(500)                              ! panel 0   - y pass
            a(sy)=il_g
            b(sy)=-1
            c(sy)=1
          case(510)                              ! panel 3   - y pass
            a(sy)=-1
            b(sy)=-il_g
            c(sy)=il_g*(4*il_g+1)+1
          case(520)                              ! panel 5   - x pass
            a(sy)=il_g
            b(sy)=-1
            c(sy)=5*il_g*il_g+1
          case(522)                              ! panel 2   - x pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(2*il_g-3)
          case(530)                              ! panel 5   - z pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=il_g*(5*il_g-1)            
          case(531)                              ! panel 0   - z pass
            a(sy)=1
            b(sy)=il_g
            c(sy)=-2*il_g
          case DEFAULT
            write(6,*) "Invalid index ",ppass,ipass,sn,
     &              ppass*100+ipass*10+sy
            stop
        end select
 
      end do

      return
      end subroutine getiqa

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
          write(6,*) "Gather data for MLO filter     ",kbb
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
      real, dimension(ifull,kd) :: diffw_l
      real, dimension(ifull_g,kd) :: diff_g
      real, dimension(ifull) :: xa_l, xb_l
      logical, intent(in) :: lblock
      logical, dimension(ifull_g) :: landg

      if (myid==0) then
        write(6,*) "MLO 2D scale-selective filter"
        if (kd==1) then
          write(6,*) "Single level filter"
        else
          write(6,*) "Multiple level filter"
        end if
      end if

      if (nud_sst/=0) then
        call ccmpi_gatherall(diff_l(:,1:kd),diff_g(:,1:kd))
        landg=abs(diff_g(:,1)-miss)<0.1
        call mlofilterhost(diff_g,diff_l,kd,miss,landg)
      end if
      if (nud_sss/=0) then
        call ccmpi_gatherall(diffs_l(:,1:kd),diff_g(:,1:kd))
        landg=abs(diff_g(:,1)-miss)<0.1
        call mlofilterhost(diff_g,diffs_l,kd,miss,landg)
      end if
      if (nud_ouv/=0) then
        do k=1,kd
          xa_l=diffu_l(:,k)
          xb_l=diffv_l(:,k)
          diffu_l(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
          diffv_l(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
          diffw_l(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
          where (abs(xa_l-miss)<0.1)
            diffu_l(:,k)=miss
            diffv_l(:,k)=miss
            diffw_l(:,k)=miss
          end where
        end do        
        call ccmpi_gatherall(diffu_l(:,1:kd),diff_g(:,1:kd))
        landg=abs(diff_g(:,1)-miss)<0.1
        call mlofilterhost(diff_g,diffu_l,kd,miss,landg)
        call ccmpi_gatherall(diffv_l(:,1:kd),diff_g(:,1:kd))
        call mlofilterhost(diff_g,diffv_l,kd,miss,landg)
        call ccmpi_gatherall(diffw_l(:,1:kd),diff_g(:,1:kd))
        call mlofilterhost(diff_g,diffw_l,kd,miss,landg)
      end if
      if (nud_sfh/=0.and.lblock) then
        call ccmpi_gatherall(diffh_l(:,1),diff_g(:,1))
        landg=abs(diff_g(:,1)-miss)<0.1
        call mlofilterhost(diff_g,diffh_l,1,miss,landg)
      end if

#ifdef debug
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "MLO end 2D filter"
      end if
#endif
      
      return
      end subroutine mlofilter

      subroutine mlofilterhost(diff_g,dd,kd,miss,landg)

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
      real, dimension(ifull_g,kd), intent(inout) :: diff_g
      real, dimension(ifull_g) :: rr,mm,nn
      real, dimension(ifull) :: ddh
      real, dimension(ifull,kd), intent(out) :: dd
      logical, dimension(ifull_g), intent(in) :: landg

      ! eventually will be replaced with mbd once full ocean coupling is complete
      cq=sqrt(4.5)*.1*real(max(nud_sst,nud_sss,nud_ouv,nud_sfh,mbd))
     &   /(pi*schmidt)

      dd=0.
      mm=1./(em_g*em_g)
      where(.not.landg)
        nn=mm
      elsewhere
        nn=0.
      end where
      do k=1,kd
        diff_g(:,k)=diff_g(:,k)*nn
      end do

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
              do k=1,kd
                dd(iqq,k)=sum(rr(:)*diff_g(:,k))/nsum
              end do
            end if
          end do
        end do
      end do

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
      integer ifg,xpan
      real, intent(in) :: miss      
      real, dimension(ifull,1), intent(inout) :: diffh_l
      real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
      real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
      real cq
      logical, intent(in) :: lblock
      
      ! eventually will be replaced with mbd once full ocean coupling is complete
      cq=sqrt(4.5)*.1*real(max(nud_sst,nud_sss,nud_ouv,nud_sfh,mbd))
     &   /(pi*schmidt)
      
#ifdef uniform_decomp
      if (myid==0) then
        write(6,*) "MLO 1D scale-selective filter (MPI)"
        if (kd==1) then
          write(6,*) "Single level filter"
        else
          write(6,*) "Multiple level filter"
        end if
      end if        
#else
      if (myid==0) then
        write(6,*) "MLO 1D scale-selective filter (MPI optimised)"
        if (kd==1) then
          write(6,*) "Single level filter"
        else
          write(6,*) "Multiple level filter"
        end if
      end if
#endif

      if (joff==0) then
        ifg=ifull_g
      else
        ifg=0
      end if
      xpan=max(ipan,jpan)
      if (pprocn==pprocx) then
        call mlospechost_n(cq,diff_l,diffs_l,diffu_l,diffv_l,
     &                   diffh_l,miss,lblock,kd,ifg,xpan)
      else
        call mlospechost(cq,diff_l,diffs_l,diffu_l,diffv_l,
     &                   diffh_l,miss,lblock,kd,ifg,xpan)
      end if

#ifdef debug
      if (myid==0.and.nmaxpr==1) then
        write(6,*) "MLO end 1D filter"
      end if
#endif

      return
      end subroutine mlofilterfast

      subroutine mlospechost(cq,diff_l,diffs_l,diffu_l,diffv_l,
     &                       diffh_l,miss,lblock,kd,ifg,xpan)

      use cc_mpi             ! CC MPI routines
      use vecsuv_m           ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: kd,ifg,xpan
      integer ppass,nne,nns,oos,ooe,iproc
      integer n,j,k,til,ipoff,jpoff,npoff
      real, intent(in) :: cq,miss
      real, dimension(ifull,1), intent(inout) :: diffh_l
      real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
      real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
      real, dimension(ifull,kd) :: diffw_l
      real, dimension(ifg,1) :: diffh_g
      real, dimension(ifg,kd) :: diff_g,diffs_g
      real, dimension(ifg,kd) :: diffu_g,diffv_g
      real, dimension(ifg,kd) :: diffw_g
      real, dimension(ifg) :: qph
      real, dimension(ifg,kd) :: qp,qps,qpu,qpv,qpw
      real, dimension(ifg*kd) :: zz,xx
      real, dimension(ifull*kd) :: yy
      real, dimension(ifull) :: xa_l,xb_l
      logical, intent(in) :: lblock

      til=il_g*il_g 

      ! gather data onto myid==0
      if (myid==0) then
        if (nud_sst/=0) then
          call ccmpi_gather(diff_l(:,1:kd),diff_g(:,1:kd))
        end if
        if (nud_sss/=0) then
          call ccmpi_gather(diffs_l(:,1:kd),diffs_g(:,1:kd))
        end if
        if (nud_ouv/=0) then
          do k=1,kd
            xa_l=diffu_l(:,k)
            xb_l=diffv_l(:,k)
            diffu_l(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
            diffv_l(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
            diffw_l(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
            where (abs(xa_l-miss)<0.1)
              diffu_l(:,k)=miss
              diffv_l(:,k)=miss
              diffw_l(:,k)=miss
            end where
          end do
          call ccmpi_gather(diffu_l(:,1:kd),diffu_g(:,1:kd))
          call ccmpi_gather(diffv_l(:,1:kd),diffv_g(:,1:kd))
          call ccmpi_gather(diffw_l(:,1:kd),diffw_g(:,1:kd))
        end if
        if (nud_sfh/=0.and.lblock) then
          call ccmpi_gather(diffh_l(:,1),diffh_g(:,1))
        end if
      else
        if (nud_sst/=0) then
          call ccmpi_gather(diff_l(:,1:kd))
        end if
        if (nud_sss/=0) then
          call ccmpi_gather(diffs_l(:,1:kd))
        end if
        if (nud_ouv/=0) then
          do k=1,kd
            xa_l=diffu_l(:,k)
            xb_l=diffv_l(:,k)
            diffu_l(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
            diffv_l(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
            diffw_l(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
            where (abs(xa_l-miss)<0.1)
              diffu_l(:,k)=miss
              diffv_l(:,k)=miss
              diffw_l(:,k)=miss
            end where
          end do
          call ccmpi_gather(diffu_l(:,1:kd))
          call ccmpi_gather(diffv_l(:,1:kd))
          call ccmpi_gather(diffw_l(:,1:kd))
        end if        
        if (nud_sfh/=0.and.lblock) then
          call ccmpi_gather(diffh_l(:,1))
        end if
      end if
     
#ifdef debug
      if (myid==0) write(6,*) "MLO Start 1D filter"
#endif

      do ppass=pprocn,pprocx

        if (myid==hproc) then
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
        call mlospeclocal(cq,ppass,qp,qps,qpu,qpv,qpw,qph,kd,
     &                    lblock,ifg,xpan,miss)
        
        ! distribute data to processors
        if (nud_sst/=0) then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,kd
                do j=nns,nne
                  do n=oos,ooe
                    zz(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*kd*(iproc-hproc))=qp(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ipan*jpan*kd),0,
     &           comm_proc)
          do k=1,kd
            do n=1,ipan*jpan
              diff_l(n+ipan*jpan*(ppass-pprocn),k)=yy(n+ipan*jpan*(k-1))
            end do
          end do
        end if
        if (nud_sss/=0) then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,kd
                do j=nns,nne
                  do n=oos,ooe
                    zz(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ifull*kd*(iproc-hproc))=qps(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ipan*jpan*kd),0,
     &           comm_proc)
          do k=1,kd
            do n=1,ipan*jpan
              diffs_l(n+ipan*jpan*(ppass-pprocn),k)=
     &          yy(n+ipan*jpan*(k-1))
            end do
          end do
        end if
        if (nud_ouv/=0) then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,kd
                do j=nns,nne
                  do n=oos,ooe
                    zz(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*kd*(iproc-hproc))=qpu(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ipan*jpan*kd),0,
     &           comm_proc)
          do k=1,kd
            do n=1,ipan*jpan
              diffu_l(n+ipan*jpan*(ppass-pprocn),k)=
     &          yy(n+ipan*jpan*(k-1))
            end do
          end do
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,kd
                do j=nns,nne
                  do n=oos,ooe
                    zz(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*kd*(iproc-hproc))=qpv(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ipan*jpan*kd),0,
     &           comm_proc)
          do k=1,kd
            do n=1,ipan*jpan
              diffv_l(n+ipan*jpan*(ppass-pprocn),k)=
     &          yy(n+ipan*jpan*(k-1))
            end do
          end do
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do k=1,kd
                do j=nns,nne
                  do n=oos,ooe
                    zz(n-oos+1+ipan*(j-nns)+ipan*jpan*(k-1)
     &                +ipan*jpan*kd*(iproc-hproc))=qpw(n+il_g*(j-1)
     &                +ppass*til,k)
                  end do
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ipan*jpan*kd),0,
     &           comm_proc)
          do k=1,kd
            do n=1,ipan*jpan
              diffw_l(n+ipan*jpan*(ppass-pprocn),k)=
     &          yy(n+ipan*jpan*(k-1))
            end do
          end do
          do k=1,kd
            do n=1+ipan*jpan*(ppass-pprocn),ipan*jpan*(ppass-pprocn+1)
              xa_l(n)=ax(n)*diffu_l(n,k)+ay(n)*diffv_l(n,k)
     &               +az(n)*diffw_l(n,k)
              xb_l(n)=bx(n)*diffu_l(n,k)+by(n)*diffv_l(n,k)
     &               +bz(n)*diffw_l(n,k)
              diffu_l(n,k)=xa_l(n)
              diffv_l(n,k)=xb_l(n)
            end do
          end do
        end if
        if (nud_sfh/=0.and.lblock) then
          if (myid==hproc) then
            do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
              call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                             ipan,jpan)
#else
              call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                              nyproc,ipan,jpan,npan)
#endif
              nns=jpoff+1
              nne=jpoff+jpan
              oos=ipoff+1
              ooe=ipoff+ipan
              do j=nns,nne
                do n=oos,ooe
                  zz(n-oos+1+ipan*(j-nns)
     &              +ipan*jpan*(iproc-hproc))=qph(n+il_g*(j-1)
     &              +ppass*til)
                end do
              end do
            end do
          end if
          call ccmpi_scatterx(zz(1:ifg),yy(1:ipan*jpan),0,
     &           comm_proc)
          do n=1,ipan*jpan
            diffh_l(n+ipan*jpan*(ppass-pprocn),1)=yy(n)
          end do
        end if
        
      end do

#ifdef debug
      if (myid==0) then
        write(6,*) "MLO End 1D filter"
      end if
#endif

      return
      end subroutine mlospechost
      !---------------------------------------------------------------------------------
      ! memory reduced version
      subroutine mlospechost_n(cq,diff_l,diffs_l,diffu_l,diffv_l,
     &                       diffh_l,miss,lblock,kd,ifg,xpan)

      use cc_mpi             ! CC MPI routines
      use vecsuv_m           ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: kd,ifg,xpan
      integer nne,nns,oos,ooe,iproc
      integer ipoff,jpoff,npoff
      integer n,j,k,til
      real, intent(in) :: cq,miss
      real, dimension(ifull,1), intent(inout) :: diffh_l
      real, dimension(ifull,kd), intent(inout) :: diff_l,diffs_l
      real, dimension(ifull,kd), intent(inout) :: diffu_l,diffv_l
      real, dimension(ifull,kd) :: diffw_l
      real, dimension(ifg,1) :: diffh_g
      real, dimension(ifg,kd) :: diff_g,diffs_g
      real, dimension(ifg,kd) :: diffu_g,diffv_g,diffw_g
      real, dimension(ifg*kd) :: zz,xx
      real, dimension(ifull*kd) :: yy
      real, dimension(ifull) :: xa_l,xb_l
      logical, intent(in) :: lblock
      
      til=il_g*il_g 
      
      ! gather data
      if (nud_sst/=0) then
        do k=1,kd
          do n=1,ifull
            yy(n+ifull*(k-1))=diff_l(n,k)
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ifull*kd),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  xx(n+il_g*(j-1)+til*(k-1))=zz(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(zz(1:ifull_g*kd),xx(1:til*kd),
     &           comm_host)
          do iproc=0,5
            do k=1,kd
              do n=1,til
                diff_g(n+til*iproc,k)=zz(n+til*(k-1)+til*kd*iproc)
              end do
            end do
          end do
        end if
      end if
      if (nud_sss/=0) then
        do k=1,kd
          do n=1,ifull
            yy(n+ifull*(k-1))=diffs_l(n,k)
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ifull*kd),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  xx(n+il_g*(j-1)+til*(k-1))=zz(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(zz(1:ifull_g*kd),xx(1:til*kd),
     &           comm_host)
          do iproc=0,5
            do k=1,kd
              do n=1,til
                diffs_g(n+til*iproc,k)=zz(n+til*(k-1)+til*kd*iproc)
              end do
            end do
          end do
        end if
      end if
      if (nud_ouv/=0) then
        do k=1,kd
          xa_l=diffu_l(:,k)
          xb_l=diffv_l(:,k)
          diffu_l(:,k)=ax(1:ifull)*xa_l+bx(1:ifull)*xb_l
          diffv_l(:,k)=ay(1:ifull)*xa_l+by(1:ifull)*xb_l
          diffw_l(:,k)=az(1:ifull)*xa_l+bz(1:ifull)*xb_l
          where (abs(xa_l-miss)<0.1)
            diffu_l(:,k)=miss
            diffv_l(:,k)=miss
            diffw_l(:,k)=miss
          end where
        end do
        do k=1,kd
          do n=1,ifull
            yy(n+ifull*(k-1))=diffu_l(n,k)
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ifull*kd),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  xx(n+il_g*(j-1)+til*(k-1))=zz(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(zz(1:ifull_g*kd),xx(1:til*kd),
     &           comm_host)
          do iproc=0,5
            do k=1,kd
              do n=1,til
                diffu_g(n+til*iproc,k)=zz(n+til*(k-1)+til*kd*iproc)
              end do
            end do
          end do
        end if
        do k=1,kd
          do n=1,ifull
            yy(n+ifull*(k-1))=diffv_l(n,k)
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ifull*kd),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  xx(n+il_g*(j-1)+til*(k-1))=zz(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(zz(1:ifull_g*kd),xx(1:til*kd),
     &           comm_host)
          do iproc=0,5
            do k=1,kd
              do n=1,til
                diffv_g(n+til*iproc,k)=zz(n+til*(k-1)+til*kd*iproc)
              end do
            end do
          end do
        end if
        do k=1,kd
          do n=1,ifull
            yy(n+ifull*(k-1))=diffw_l(n,k)
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ifull*kd),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  xx(n+il_g*(j-1)+til*(k-1))=zz(n-oos+1
     &              +ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
          call ccmpi_allgatherx(zz(1:ifull_g*kd),xx(1:til*kd),
     &           comm_host)
          do iproc=0,5
            do k=1,kd
              do n=1,til
                diffw_g(n+til*iproc,k)=zz(n+til*(k-1)+til*kd*iproc)
              end do
            end do
          end do
        end if
      end if        
      if (nud_sfh/=0.and.lblock) then
        do n=1,ifull
          yy(n)=diffh_l(n,1)
        end do
        call ccmpi_gatherx(zz(1:til),yy(1:ifull),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do j=nns,nne
              do n=oos,ooe
                xx(n+il_g*(j-1))=zz(n-oos+1
     &            +ipan*(j-nns)
     &            +ifull*(iproc-hproc))
              end do
            end do
          end do
          call ccmpi_allgatherx(zz(1:ifull_g),xx(1:til),
     &           comm_host)
          do iproc=0,5
            do n=1,til
              diffh_g(n+til*iproc,1)=zz(n+til*iproc)
            end do
          end do
        end if
      end if
      
#ifdef debug
      if (myid==0) write(6,*) "MLO Start 1D filter"
#endif

        ! computations for the local processor group
      call mlospeclocal(cq,pprocn,diff_g,diffs_g,diffu_g,
     &                  diffv_g,diffw_g,diffh_g,kd,lblock,ifg,
     &                  xpan,miss)

#ifdef debug
      if (myid==0) then
        write(6,*) "MLO End 1D filter"
        write(6,*) "MLO Receive arrays from hproc proc"
      end if
#endif

      ! distribute data to processors
      if (nud_sst/=0) then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  zz(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))=diff_g(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ifull*kd),0,comm_proc)
        do k=1,kd
          do n=1,ifull
            diff_l(n,k)=yy(n+ifull*(k-1))
          end do
        end do
      end if
      if (nud_sss/=0) then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  zz(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))=diffs_g(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ifull*kd),0,comm_proc)
        do k=1,kd
          do n=1,ifull
            diffs_l(n,k)=yy(n+ifull*(k-1))
          end do
        end do
      end if
      if (nud_ouv/=0) then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  zz(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))=diffu_g(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ifull*kd),0,comm_proc)
        do k=1,kd
          do n=1,ifull
            diffu_l(n,k)=yy(n+ifull*(k-1))
          end do
        end do
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  zz(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))=diffv_g(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ifull*kd),0,comm_proc)
        do k=1,kd
          do n=1,ifull
            diffv_l(n,k)=yy(n+ifull*(k-1))
          end do
        end do
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  zz(n-oos+1+ipan*(j-nns)+ifull*(k-1)
     &              +ifull*kd*(iproc-hproc))=diffw_g(n+il_g*(j-1)
     &              +pprocn*til,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(zz(1:ifg*kd),yy(1:ifull*kd),0,comm_proc)
        do k=1,kd
          do n=1,ifull
            diffw_l(n,k)=yy(n+ifull*(k-1))
          end do
        end do
        do k=1,kd
          xa_l=ax(1:ifull)*diffu_l(:,k)+ay(1:ifull)*diffv_l(:,k)
     &        +az(1:ifull)*diffw_l(:,k)
          xb_l=bx(1:ifull)*diffu_l(:,k)+by(1:ifull)*diffv_l(:,k)
     &        +bz(1:ifull)*diffw_l(:,k)
          diffu_l(:,k)=xa_l
          diffv_l(:,k)=xb_l
        end do
      end if
      if (nud_sfh/=0.and.lblock) then
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
            nns=jpoff+1
            nne=jpoff+jpan
            oos=ipoff+1
            ooe=ipoff+ipan
            do j=nns,nne
              do n=oos,ooe
                zz(n-oos+1+ipan*(j-nns)
     &            +ifull*(iproc-hproc))=diffh_g(n+il_g*(j-1)
     &            +pprocn*til,1)
              end do
            end do
          end do
        end if
        call ccmpi_scatterx(zz(1:ifg),diffh_l(1:ifull,1),0,comm_proc)
      end if

      return
      end subroutine mlospechost_n

      
      !---------------------------------------------------------------------------------
      ! This version is for asymmetric decomposition
      subroutine mlospeclocal(cq,ppass,qp,qps,qpu,qpv,qpw,qph,
     &             kd,lblock,ifg,xpan,miss)

      use cc_mpi             ! CC MPI routines
      use map_m              ! Grid map arrays
      use xyzinfo_m          ! Grid coordinate arrays
     
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: ppass,kd,ifg,xpan
      integer j,n,ipass,ns,ne,os,oe
      integer iproc,istep
      integer ipoff,jpoff,npoff
      integer nne,nns,me,ooe,oos
      integer k,til,sn,sy,a,b,c,jj,nn
      integer, dimension(0:3) :: astr,bstr,cstr
      integer, dimension(0:3) :: maps
      real, intent(in) :: cq,miss
      real, dimension(ifg), intent(inout) :: qph
      real, dimension(ifg,kd), intent(inout) :: qp,qps
      real, dimension(ifg,kd), intent(inout) :: qpu,qpv
      real, dimension(ifg,kd), intent(inout) :: qpw
      real, dimension(ifg) :: qsum
      real, dimension(4*il_g) :: rr,xa,ya,za
      real, dimension(4*il_g,xpan) :: aph,asum,rsum
      real, dimension(4*il_g,xpan,kd) :: ap,aps,apu,apv,apw
      real, dimension(xpan,xpan) :: pph,psum
      real, dimension(xpan,xpan,kd) :: pp,pps
      real, dimension(xpan,xpan,kd) :: ppu,ppv
      real, dimension(xpan,xpan,kd) :: ppw
      real, dimension(4*il_g*il_g*kd) :: zz
      real, dimension(4*il_g*xpan*kd) :: yy
      logical, intent(in) :: lblock
      logical, dimension(4*il_g,xpan) :: landl
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      til=il_g*il_g
      astr=0
      bstr=0
      cstr=0

      ns=joff+1
      ne=joff+jpan
      os=ioff+1
      oe=ioff+ipan
      
      do ipass=0,2
        me=maps(ipass)
        call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
        if (myid==0) then
          write(6,*) "MLO Recieve arrays from local host"
        end if
#endif
        ! myid==hproc should have ioff=0 and joff=0
        if (nud_sfh/=0.and.lblock) then
          if (ioff==0) then
            if (myid==hproc) then
              do j=1,il_g
                do sn=1,me,il_g
                  sy=(sn-1)/il_g
                  a=astr(sy)
                  b=bstr(sy)
                  c=cstr(sy)
                  do n=sn,sn+il_g-1
                    zz(n+me*(j-1))=qph(a*n+b*j+c)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:me*il_g),yy(1:me*jpan),0,
     &             comm_cols)
          end if
          call ccmpi_bcast(yy(1:me*jpan),0,comm_rows)
          do j=1,jpan
            do n=1,me
              aph(n,j)=yy(n+me*(j-1))
            end do
          end do
        end if
        if (nud_sst/=0) then
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,kd
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      zz(n+me*(k-1)+me*kd*(j-1))=qp(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:me*il_g),yy(1:me*jpan*kd),0,
     &             comm_cols)
          end if
          call ccmpi_bcast(yy(1:me*jpan*kd),0,comm_rows)
          do k=1,kd
            do j=1,jpan
              do n=1,me
                ap(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
              end do
            end do
          end do
        end if
        if (nud_sss/=0) then
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,kd
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      zz(n+me*(k-1)+me*kd*(j-1))=qps(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:me*il_g),yy(1:me*jpan*kd),0,
     &             comm_cols)
          end if
          call ccmpi_bcast(yy(1:me*jpan*kd),0,comm_rows)
          do k=1,kd
            do j=1,jpan
              do n=1,me
                aps(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
              end do
            end do
          end do
        end if
        if (nud_ouv/=0) then
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,kd
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      zz(n+me*(k-1)+me*kd*(j-1))=qpu(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:me*il_g),yy(1:me*jpan*kd),0,
     &             comm_cols)
          end if
          call ccmpi_bcast(yy(1:me*jpan*kd),0,comm_rows)
          do k=1,kd
            do j=1,jpan
              do n=1,me
                apu(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
              end do
            end do
          end do
          if (ioff==0) then          
            if (myid==hproc) then
              do k=1,kd
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      zz(n+me*(k-1)+me*kd*(j-1))=qpv(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:me*il_g),yy(1:me*jpan*kd),0,
     &             comm_cols)
          end if
          call ccmpi_bcast(yy(1:me*jpan*kd),0,comm_rows)
          do k=1,kd
            do j=1,jpan
              do n=1,me
                apv(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
              end do
            end do
          end do
          if (ioff==0) then
            if (myid==hproc) then
              do k=1,kd
                do j=1,il_g
                  do sn=1,me,il_g
                    sy=(sn-1)/il_g
                    a=astr(sy)
                    b=bstr(sy)
                    c=cstr(sy)
                    do n=sn,sn+il_g-1
                      zz(n+me*(k-1)+me*kd*(j-1))=qpw(a*n+b*j+c,k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:me*il_g),yy(1:me*jpan*kd),0,
     &             comm_cols)
          end if
          call ccmpi_bcast(yy(1:me*jpan*kd),0,comm_rows)
          do k=1,kd
            do j=1,jpan
              do n=1,me
                apw(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
              end do
            end do
          end do
        end if

        if (nud_sfh/=0.and.lblock) then
          landl(1:me,1:jpan)=abs(aph(1:me,1:jpan)-miss)<0.1
        end if
        if (nud_sst/=0) then
          landl(1:me,1:jpan)=abs(ap(1:me,1:jpan,1)-miss)<0.1
        end if
        if (nud_sss/=0) then
          landl(1:me,1:jpan)=abs(aps(1:me,1:jpan,1)-miss)<0.1
        end if
        if (nud_ouv/=0) then
          landl(1:me,1:jpan)=abs(apw(1:me,1:jpan,1)-miss)<0.1
        end if
        do j=1,jpan
          jj=j+ns-1
          do sn=1,me,il_g
            sy=(sn-1)/il_g
            a=astr(sy)
            b=bstr(sy)
            c=cstr(sy)
            do n=sn,sn+il_g-1
              asum(n,j)=1./em_g(a*n+b*jj+c)**2
            end do
          end do
        end do
        where (landl(1:me,1:jpan))
          rsum(1:me,1:jpan)=0.
        elsewhere
          rsum(1:me,1:jpan)=asum(1:me,1:jpan)
        end where
        if (nud_sfh/=0.and.lblock) then
          do j=1,jpan
            do n=1,me
              aph(n,j)=aph(n,j)*rsum(n,j)
            end do
          end do          
        end if
        if (nud_sst/=0) then
          do k=1,kd
            do j=1,jpan
              do n=1,me
                ap(n,j,k)=ap(n,j,k)*rsum(n,j)
              end do
            end do
          end do
        end if
        if (nud_sss/=0) then
          do k=1,kd
            do j=1,jpan
              do n=1,me
                aps(n,j,k)=aps(n,j,k)*rsum(n,j)
              end do
            end do
          end do
        end if
        if (nud_ouv/=0) then
          do k=1,kd
            do j=1,jpan
              do n=1,me
                apu(n,j,k)=apu(n,j,k)*rsum(n,j)
                apv(n,j,k)=apv(n,j,k)*rsum(n,j)
                apw(n,j,k)=apw(n,j,k)*rsum(n,j)
              end do
            end do
          end do          
        end if

#ifdef debug
        if (myid==0) write(6,*) "MLO start convolution"
#endif

        do j=1,jpan
          jj=j+ns-1
          do sn=1,me,il_g
            sy=(sn-1)/il_g
            a=astr(sy)
            b=bstr(sy)
            c=cstr(sy)
            do n=sn,sn+il_g-1
              xa(n)=x_g(a*n+b*jj+c)
              ya(n)=y_g(a*n+b*jj+c)
              za(n)=z_g(a*n+b*jj+c)
            end do
          end do
          do n=1,ipan
            nn=n+os-1
            rr(1:me)=xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me)
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

#ifdef debug
        if (myid==0) then
          write(6,*) "MLO end conv"
          write(6,*) "MLO Send arrays to local host"
        end if
#endif

        ! unpack grid
        a=astr(0)
        b=bstr(0)
        c=cstr(0)

        ! gather data on host processors
        do j=1,jpan
          do n=1,ipan
            yy(n+ipan*(j-1))=psum(n,j)
          end do
        end do
        call ccmpi_gatherx(zz(1:il_g*ipan),yy(1:ipan*jpan),0,
     &         comm_cols)
        if (joff==0) then
          do jpoff=0,il_g-1,jpan
            sy=jpoff/jpan
            nns=jpoff+1
            nne=jpoff+jpan
            do j=nns,nne
              do n=os,oe
                qsum(a*n+b*j+c)=zz(n-os+1+ipan*(j-nns)
     &                            +ipan*jpan*sy)
              end do
            end do
          end do
        end if
        if (nud_sfh/=0.and.lblock) then
          do j=1,jpan
            do n=1,ipan
              yy(n+ipan*(j-1))=pph(n,j)
            end do
          end do
          call ccmpi_gatherx(zz(1:il_g*ipan),yy(1:ipan*jpan),0,
     &           comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do j=nns,nne
                do n=os,oe
                  qph(a*n+b*j+c)=zz(n-os+1+ipan*(j-nns)
     &                             +ipan*jpan*sy)
                end do
              end do
            end do
          end if
        end if
        if (nud_sst/=0) then
          do k=1,kd
            do j=1,jpan
              do n=1,ipan
                yy(n+ipan*(j-1)+ipan*jpan*(k-1))=pp(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(zz(1:il_g*ipan*kd),yy(1:ipan*jpan*kd),0,
     &           comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,kd
                do j=nns,nne
                  do n=os,oe
                    qp(a*n+b*j+c,k)=zz(n-os+1+ipan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*kd*sy)
                  end do
                end do
              end do
            end do
          end if
        end if
        if (nud_sss/=0) then
          do k=1,kd
            do j=1,jpan
              do n=1,ipan
                yy(n+ipan*(j-1)+ipan*jpan*(k-1))=pps(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(zz(1:il_g*ipan*kd),yy(1:ipan*jpan*kd),0,
     &           comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,kd
                do j=nns,nne
                  do n=os,oe
                    qps(a*n+b*j+c,k)=zz(n-os+1+ipan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*kd*sy)
                  end do
                end do
              end do
            end do
          end if
        end if
        if (nud_ouv/=0) then
          do k=1,kd
            do j=1,jpan
              do n=1,ipan
                yy(n+ipan*(j-1)+ipan*jpan*(k-1))=ppu(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(zz(1:il_g*ipan*kd),yy(1:ipan*jpan*kd),0,
     &           comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,kd
                do j=nns,nne
                  do n=os,oe
                    qpu(a*n+b*j+c,k)=zz(n-os+1+ipan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*kd*sy)
                  end do
                end do
              end do
            end do
          end if
          do k=1,kd
            do j=1,jpan
              do n=1,ipan
                yy(n+ipan*(j-1)+ipan*jpan*(k-1))=ppv(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(zz(1:il_g*ipan*kd),yy(1:ipan*jpan*kd),0,
     &           comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,kd
                do j=nns,nne
                  do n=os,oe
                    qpv(a*n+b*j+c,k)=zz(n-os+1+ipan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*kd*sy)
                  end do
                end do
              end do
            end do
          end if
          do k=1,kd
            do j=1,jpan
              do n=1,ipan
                yy(n+ipan*(j-1)+ipan*jpan*(k-1))=ppw(n,j,k)
              end do
            end do
          end do
          call ccmpi_gatherx(zz(1:il_g*ipan*kd),yy(1:ipan*jpan*kd),0,
     &           comm_cols)
          if (joff==0) then
            do jpoff=0,il_g-1,jpan
              sy=jpoff/jpan
              nns=jpoff+1
              nne=jpoff+jpan
              do k=1,kd
                do j=nns,nne
                  do n=os,oe
                    qpw(a*n+b*j+c,k)=zz(n-os+1+ipan*(j-nns)
     &                +ipan*jpan*(k-1)+ipan*jpan*kd*sy)
                  end do
                end do
              end do
            end do
          end if
        end if
          
      end do

      ns=ioff+1
      ne=ioff+ipan
      os=joff+1
      oe=joff+jpan

      ipass=3
      me=maps(ipass)
      call getiqa(astr,bstr,cstr,me,ipass,ppass,il_g)

#ifdef debug
      if (myid==0) then
        write(6,*) "MLO Recieve arrays from local host"
      end if
#endif

      if (joff==0) then
        do j=ns,ne
          do sn=1,me,il_g
            sy=(sn-1)/il_g
            a=astr(sy)
            b=bstr(sy)
            c=cstr(sy)
            do n=sn,sn+il_g-1
              yy(n+me*(j-ns))=qsum(a*n+b*j+c)
            end do
          end do
        end do
      end if
      call ccmpi_bcast(yy(1:me*ipan),0,comm_cols)
      do j=1,ipan
        do n=1,me
          asum(n,j)=yy(n+me*(j-1))
        end do
      end do
      if (nud_sfh/=0.and.lblock) then
        if (joff==0) then
          do j=ns,ne
            do sn=1,me,il_g
              sy=(sn-1)/il_g
              a=astr(sy)
              b=bstr(sy)
              c=cstr(sy)
              do n=sn,sn+il_g-1
                yy(n+me*(j-ns))=qph(a*n+b*j+c)
              end do
            end do
          end do
        end if
        call ccmpi_bcast(yy(1:me*ipan),0,comm_cols)
        do j=1,ipan
          do n=1,me
            aph(n,j)=yy(n+me*(j-1))
          end do
        end do
      end if
      if (nud_sst/=0) then
        if (joff==0) then
          do k=1,kd
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  yy(n+me*(k-1)+me*kd*(j-ns))=qp(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(yy(1:me*ipan*kd),0,comm_cols)
        do k=1,kd
          do j=1,ipan
            do n=1,me
              ap(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
            end do
          end do
        end do
      end if
      if (nud_sss/=0) then
        if (joff==0) then
          do k=1,kd
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  yy(n+me*(k-1)+me*kd*(j-ns))=qps(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(yy(1:me*ipan*kd),0,comm_cols)
        do k=1,kd
          do j=1,ipan
            do n=1,me
              aps(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
            end do
          end do
        end do
      end if
      if (nud_ouv/=0) then
        if (joff==0) then
          do k=1,kd
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  yy(n+me*(k-1)+me*kd*(j-ns))=qpu(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(yy(1:me*ipan*kd),0,comm_cols)
        do k=1,kd
          do j=1,ipan
            do n=1,me
              apu(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
            end do
          end do
        end do
        if (joff==0) then
          do k=1,kd
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  yy(n+me*(k-1)+me*kd*(j-ns))=qpv(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(yy(1:me*ipan*kd),0,comm_cols)
        do k=1,kd
          do j=1,ipan
            do n=1,me
              apv(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
            end do
          end do
        end do
        if (joff==0) then
          do k=1,kd
            do j=ns,ne
              do sn=1,me,il_g
                sy=(sn-1)/il_g
                a=astr(sy)
                b=bstr(sy)
                c=cstr(sy)
                do n=sn,sn+il_g-1
                  yy(n+me*(k-1)+me*kd*(j-ns))=qpw(a*n+b*j+c,k)
                end do
              end do
            end do
          end do
        end if
        call ccmpi_bcast(yy(1:me*ipan*kd),0,comm_cols)
        do k=1,kd
          do j=1,ipan
            do n=1,me
              apw(n,j,k)=yy(n+me*(k-1)+me*kd*(j-1))
            end do
          end do
        end do
      end if

#ifdef debug
      if (myid==0) write(6,*) "MLO start convolution"
#endif

      do j=1,ipan
        jj=j+ns-1
        do sn=1,me,il_g
          sy=(sn-1)/il_g
          a=astr(sy)
          b=bstr(sy)
          c=cstr(sy)
          do n=sn,sn+il_g-1
            xa(n)=x_g(a*n+b*jj+c)
            ya(n)=y_g(a*n+b*jj+c)
            za(n)=z_g(a*n+b*jj+c)
          end do
        end do
        do n=1,jpan
          nn=n+os-1
          rr(1:me)=xa(nn)*xa(1:me)+ya(nn)*ya(1:me)+za(nn)*za(1:me)
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

#ifdef debug
      if (myid==0) then
        write(6,*) "MLO end conv"
        write(6,*) "MLO Send arrays to local host"
      end if
#endif

      ! unpack grid
      a=astr(0)
      b=bstr(0)
      c=cstr(0)

      ! gather data on host processors
      if (nud_sfh/=0.and.lblock) then
        do j=1,ipan
          do n=1,jpan
            if (psum(n,j)>1.E-8) then
              yy(n+jpan*(j-1))=pph(n,j)/psum(n,j)
            else
              yy(n+jpan*(j-1))=0.
            end if
          end do
        end do
        call ccmpi_gatherx(zz(1:til),yy(1:ipan*jpan),0,comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do j=nns,nne
              do n=oos,ooe
                qph(a*n+b*j+c)=zz(n-oos+1+jpan*(j-nns)
     &                           +ipan*jpan*(iproc-hproc))
              end do
            end do
          end do
        end if
      end if
      if (nud_sst/=0) then
        do k=1,kd
          do j=1,ipan
            do n=1,jpan
              if (psum(n,j)>1.E-8) then
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=pp(n,j,k)/psum(n,j)
              else
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=0.
              end if
            end do
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ipan*jpan*kd),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  qp(a*n+b*j+c,k)=zz(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
      end if
      if (nud_sss/=0) then
        do k=1,kd
          do j=1,ipan
            do n=1,jpan
              if (psum(n,j)>1.E-8) then
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=pps(n,j,k)/psum(n,j)
              else
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=0.
              end if
            end do
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ipan*jpan*kd),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  qps(a*n+b*j+c,k)=zz(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
      end if
      if (nud_ouv/=0) then
        do k=1,kd
          do j=1,ipan
            do n=1,jpan
              if (psum(n,j)>1.E-8) then
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=ppu(n,j,k)/psum(n,j)
              else
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=0.
              end if
            end do
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ipan*jpan*kd),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  qpu(a*n+b*j+c,k)=zz(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
        do k=1,kd
          do j=1,ipan
            do n=1,jpan
              if (psum(n,j)>1.E-8) then
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=ppv(n,j,k)/psum(n,j)
              else
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=0.
              end if
            end do
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ipan*jpan*kd),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  qpv(a*n+b*j+c,k)=zz(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
        do k=1,kd
          do j=1,ipan
            do n=1,jpan
              if (psum(n,j)>1.E-8) then
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=ppw(n,j,k)/psum(n,j)
              else
                yy(n+jpan*(j-1)+ipan*jpan*(k-1))=0.
              end if
            end do
          end do
        end do
        call ccmpi_gatherx(zz(1:til*kd),yy(1:ipan*jpan*kd),0,
     &         comm_proc)
        if (myid==hproc) then
          do iproc=hproc,hproc+mproc-1
#ifdef uniform_decomp
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,
     &                           ipan,jpan)
#else
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,
     &                            nyproc,ipan,jpan,npan)
#endif
            nns=ipoff+1
            nne=ipoff+ipan
            oos=jpoff+1
            ooe=jpoff+jpan
            do k=1,kd
              do j=nns,nne
                do n=oos,ooe
                  qpw(a*n+b*j+c,k)=zz(n-oos+1+jpan*(j-nns)
     &              +ipan*jpan*(k-1)+ipan*jpan*kd*(iproc-hproc))
                end do
              end do
            end do
          end do
        end if
      end if
          
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
      
