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
              tggsn(iq,1)=min(271.2,tssb(iq),tb(iq,1)+.04*6.5) ! for 40 m lev1 ! MJT seaice
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
!          write (6,"('100*psl.wesn ',2p5f8.3)") psl(idjd),
!     &           psl(iw(idjd)),psl(ie(idjd)),psl(is(idjd)),psl(in(idjd))
!          write (6,"('ps.wesn ',-2p5f9.3)") ps(idjd),
!     &          ps(iw(idjd)),ps(ie(idjd)),ps(is(idjd)),ps(in(idjd))
          write(6,*) 'pslb in(idjd) :',pslb(idjd)
          write(6,*) 
     &     'call retopo from nestin; psl# prints refer to pslb'
        endif
        call retopo(pslb,zsb,zs(1:ifull),tb,qb)
!        if(nmaxpr==1.and.mydiag)then
!           write (6,"('100*pslb.wesn ',2p5f8.3)") pslb(idjd),
!     &       pslb(iw(idjd)),pslb(ie(idjd)),pslb(is(idjd)),pslb(in(idjd))
!        endif

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
                !call mloexpmelt(timelt)
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
      real, dimension(ifull,kblock) :: diffu,diffv
      real, dimension(ifull_g) :: psld,x_g,xx_g
      real, dimension(ifull_g,kblock) :: ud,vd,wd,td,qd
      real den,polenx,poleny,polenz,zonx,zony,zonz
      logical disflag,lblock

      ! nud_uv<0 (JLM 3-pass filter)
      ! nud_uv=0 (no preturbing of winds)
      ! nud_uv=1 (1D scale-selective filter)
      ! nud_uv=3 (JLM preturb zonal winds with 1D filter)
      ! nud_uv=9 (2D scale-selective filter)

      if (nud_p>0) then
        call ccmpi_gatherall(pslb(:), psld(:))
      end if
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
      
      ! Loop over maximum block size
      ! kblock can be reduced to save memory
      do kb=kbotdav,ktopdav,kblock
        if (myid == 0) then     
          write(6,*) "Gather data for spectral filter ",kb
        end if      
        kln=kb
        klx=min(kb+kblock-1,ktopdav)
        klt=klx-kln+1
        lblock=(kb==kbotdav)
        if (nud_uv==3) then
          do k=kln,klx
            klc=k-kln+1
            diffu(:,klc)=costh(:)*ub(:,k)  ! uzon
     &                  -sinth(:)*vb(:,k)
          end do
          call ccmpi_gatherall(diffu(:,1:klt),wd(:,1:klt))
        else if (nud_uv/=0) then
          call ccmpi_gatherall(ub(:,kln:klx),ud(:,1:klt))
          call ccmpi_gatherall(vb(:,kln:klx),vd(:,1:klt))
          do k=1,klt
            x_g=ud(:,k)
            xx_g=vd(:,k)
            ud(:,k)=ax_g(:)*x_g+bx_g(:)*xx_g
            vd(:,k)=ay_g(:)*x_g+by_g(:)*xx_g
            wd(:,k)=az_g(:)*x_g+bz_g(:)*xx_g
          end do
        end if
        if (nud_t>0) then
          call ccmpi_gatherall(tb(:,kln:klx),td(:,1:klt))
        end if
        if (nud_q>0) then
          call ccmpi_gatherall(qb(:,kln:klx),qd(:,1:klt))
        end if

        !-----------------------------------------------------------------------
        if (nud_uv<0) then 
          if (myid==0) then
            write(6,*) "JLM Fast spectral filter ",kb
            if (nproc>1) then
              write(6,*) "Option is currently disabled"
              write(6,*) "for multi-processor configurations"
              call ccmpi_abort(-1)
            end if
            call fastspec((.1*real(mbd)/(pi*schmidt))**2
     &         ,psld(:),ud(:,1:klt),vd(:,1:klt),wd(:,1:klt)
     &         ,td(:,1:klt),qd(:,1:klt),lblock,klt)
          end if
          disflag=.true.
        else if (nud_uv==9) then 
          if (myid==0) write(6,*) "Two dimensional spectral filter ",kb
          call slowspecmpi(.1*real(mbd)/(pi*schmidt)
     &                  ,psld(:),ud(:,1:klt),vd(:,1:klt),wd(:,1:klt)
     &                  ,td(:,1:klt),qd(:,1:klt),lblock,klt)
          disflag=.false.
        else if (mod(6,nproc)==0.or.mod(nproc,6)==0) then
          if (myid==0) write(6,*) 
     &      "Separable 1D filter (MPI optimised) ",kb
          call specfastmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                  ,psld(:),ud(:,1:klt),vd(:,1:klt),wd(:,1:klt)
     &                  ,td(:,1:klt),qd(:,1:klt),lblock,klt)
          disflag=.true.
        else
          if (myid==0) write(6,*) "Separable 1D filter (MPI) ",kb
          call fourspecmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                  ,psld(:),ud(:,1:klt),vd(:,1:klt),wd(:,1:klt)
     &                  ,td(:,1:klt),qd(:,1:klt),lblock,klt)
          disflag=.true.
        endif  ! (nud_uv<0) .. else ..
        !-----------------------------------------------------------------------

        if (myid==0) then
          write(6,*) "Distribute data from spectral filter ",kb
        end if
        if (nud_uv==3) then
          if (disflag) then
            if (myid==0) then
              call ccmpi_distribute(diffu(:,1:klt),wd(:,1:klt))
            else
              call ccmpi_distribute(diffu(:,1:klt))
            end if
          else
            diffu(:,1:klt)=wd(:,1:klt)
          end if
          do k=1,klt
            ud(1:ifull,k)=costh(:)*diffu(:,k)
            vd(1:ifull,k)=-sinth(:)*diffu(:,k)	  
          end do
        else if (nud_uv/=0) then
          if (disflag) then
            if (myid==0) then
              do k=1,klt
                x_g=ax_g(:)*ud(:,k)+ay_g(:)*vd(:,k)+az_g(:)*wd(:,k)
                xx_g=bx_g(:)*ud(:,k)+by_g(:)*vd(:,k)+bz_g(:)*wd(:,k)
                ud(:,k)=x_g
                vd(:,k)=xx_g
              end do
              call ccmpi_distribute(diffu(:,1:klt),ud(:,1:klt))
              call ccmpi_distribute(diffv(:,1:klt),vd(:,1:klt))
            else
              call ccmpi_distribute(diffu(:,1:klt))
              call ccmpi_distribute(diffv(:,1:klt))
            end if
          else
            do k=1,klt
              diffu(:,k)=ax*ud(1:ifull,k)+ay*vd(1:ifull,k)
     &          +az*wd(1:ifull,k)
              diffv(:,k)=bx*ud(1:ifull,k)+by*vd(1:ifull,k)
     &          +bz*wd(1:ifull,k)
            end do
          end if
        end if
        if (nud_uv/=0) then
          u(1:ifull,kln:klx)=u(1:ifull,kln:klx)+diffu(:,1:klt)
          v(1:ifull,kln:klx)=v(1:ifull,kln:klx)+diffv(:,1:klt)
          savu(1:ifull,kln:klx)=savu(1:ifull,kln:klx)+diffu(:,1:klt)
          savu1(1:ifull,kln:klx)=savu1(1:ifull,kln:klx)+diffu(:,1:klt)
          savu2(1:ifull,kln:klx)=savu2(1:ifull,kln:klx)+diffu(:,1:klt)
          savv(1:ifull,kln:klx)=savv(1:ifull,kln:klx)+diffv(:,1:klt)
          savv1(1:ifull,kln:klx)=savv1(1:ifull,kln:klx)+diffv(:,1:klt)
          savv2(1:ifull,kln:klx)=savv2(1:ifull,kln:klx)+diffv(:,1:klt)
        end if
        if (nud_t>0) then
          if (disflag) then
            if (myid==0) then
              call ccmpi_distribute(diffu(:,1:klt),td(:,1:klt))
            else
              call ccmpi_distribute(diffu(:,1:klt))
            end if
          else
            diffu(:,1:klt)=td(1:ifull,1:klt)
          end if
          t(1:ifull,kln:klx)=t(1:ifull,kln:klx)+diffu(:,1:klt)
        end if
        if (nud_q>0) then
          if (disflag) then
            if (myid==0) then
              call ccmpi_distribute(diffu(:,1:klt),qd(:,1:klt))
            else
              call ccmpi_distribute(diffu(:,1:klt))
            end if
          else
            diffu(:,1:klt)=qd(1:ifull,1:klt)
          end if
          qg(1:ifull,kln:klx)=max(qg(1:ifull,kln:klx)
     &     +diffu(:,1:klt),0.)
          qgsav(:,kln:klx)=max(qgsav(:,kln:klx)
     &     +diffu(:,1:klt),0.)
        end if
      
      end do
      
      if (nud_p>0) then
        if (disflag) then
          if (myid==0) then
            call ccmpi_distribute(diffu(:,1), psld(:))
          else
            call ccmpi_distribute(diffu(:,1))
          end if
        else
          diffu(:,1)=psld(1:ifull)
        end if
        psl(1:ifull)=psl(1:ifull)+diffu(:,1)
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

      !--------------------------------------------------------------
      ! Fast spectral downscaling (JLM version)
      ! Only works with single processor
      subroutine fastspec(cutoff2,psla,ua,va,wa,ta,qa,lblock,klt)
   
      use indices_m          ! Grid index arrays
      use map_m              ! Grid map arrays
      use xyzinfo_m          ! Grid coordinate arrays

      implicit none
!      
      include 'newmpar.h'    ! Grid parameters
      include 'const_phys.h' ! Physical constants
      include 'parm.h'       ! Model configuration
      include 'parmgeom.h'   ! Coordinate data  

      integer, parameter :: ntest=0
      integer, intent(in) :: klt 
      integer i,j,k,n,n1,iq,iq1
      real, intent(in) :: cutoff2
      real, dimension(ifull_g), intent(inout) :: psla
      real, dimension(ifull_g,klt), intent(inout) :: ua,va
      real, dimension(ifull_g,klt), intent(inout) :: wa
      real, dimension(ifull_g,klt), intent(inout) :: ta,qa
      real, dimension(ifull_g) :: psls,sumwt
      real, dimension(ifull_g) :: psls2
      real, dimension(:), allocatable, save :: xx,yy,zz
      real, dimension(ifull_g,klt) :: uu,vv,ww,tt,qgg
      real, dimension(ifull_g,klt) :: uu2,vv2,ww2,tt2,qgg2
      real emmin,dist,dist1,wt,wt1,xxmax,yymax,zzmax
      logical, intent(in) :: lblock
      
      ! myid must = 0 to get here.  So there is no need to check.
 
      if (.not.allocated(xx)) then
        allocate(xx(ifull_g),yy(ifull_g),zz(ifull_g))

!       set up geometry for filtering through panel 1
!       x pass on panels 1, 2, 4, 5
!       y pass on panels 0, 1, 3, 4
!       z pass on panels 0, 2, 3, 5
        xx=0.
        yy=0.
        zz=0.
        do iq=1+il_g*il_g,3*il_g*il_g
          xx(iq)=xx(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
        enddo
         do iq=1+4*il_g*il_g,6*il_g*il_g
          xx(iq)=xx(is_g(iq))+sqrt((x_g(iq)-x_g(is_g(iq)))**2+
     &           (y_g(iq)-y_g(is_g(iq)))**2+(z_g(iq)-z_g(is_g(iq)))**2)
        enddo
        do iq=1,2*il_g*il_g
          yy(iq)=yy(is_g(iq))+sqrt((x_g(iq)-x_g(is_g(iq)))**2+
     &           (y_g(iq)-y_g(is_g(iq)))**2+(z_g(iq)-z_g(is_g(iq)))**2)
        enddo
        do iq=1+3*il_g*il_g,5*il_g*il_g
          yy(iq)=yy(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
        enddo
        if(mbd>0)then
         do iq=1,il_g*il_g
          zz(iq)=zz(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
         enddo
         do iq=1+2*il_g*il_g,4*il_g*il_g
          zz(iq)=zz(is_g(iq))+sqrt((x_g(iq)-x_g(is_g(iq)))**2+
     &           (y_g(iq)-y_g(is_g(iq)))**2+(z_g(iq)-z_g(is_g(iq)))**2)
         enddo
         do iq=1+5*il_g*il_g,6*il_g*il_g
          zz(iq)=zz(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
         enddo
        endif  ! (mbd>0)
        if(ntest>0)then
          do iq=1,144
           write(6,*) 'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
          do iq=il_g*il_g,il_g*il_g+il_g
           write(6,*) 'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
         do iq=4*il_g*il_g-il_g,4*il_g*il_g
           write(6,*) 'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
         do iq=5*il_g*il_g-il_g,5*il_g*il_g
           write(6,*) 'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
          write(6,*) 'xx mid:'   
          do i=1,48
           write(6,*) 'i xx',i,xx(il_g*il_g*1.5+i)
          enddo
          do i=1,48
           write(6,*) 'i xx',i+il_g,xx(il_g*il_g*2.5+i)
          enddo
          do i=1,48
           write(6,*) 'i xx',i+2*il_g,xx(il_g*il_g*4-il_g/2+i*il_g)
          enddo
          do i=1,48
           write(6,*) 'i xx',i+3*il_g,xx(il_g*il_g*5-il_g/2+i*il_g)
          enddo
          write(6,*) 'yy mid:'   
          do j=1,96
           write(6,*) 'j yy',j,yy(-il_g/2+j*il_g)
          enddo
          do j=1,48
           write(6,*) 'j yy',j+2*il_g,yy(il_g*il_g*3.5+j)
          enddo
          do j=1,48
           write(6,*) 'j yy',j+3*il_g,yy(il_g*il_g*4.5+j)
          enddo
!         wrap-around values defined by xx(il_g,5*il_g+j),j=1,il_g; yy(i,5*il_g),i=1,il_g
          write(6,*) 'wrap-round values'
          do j=1,il_g
           write(6,*) 'j,xx ',j,xx(6*il_g*il_g+1-j)       ! xx(il_g+1-j,il_g,5)
          enddo
          do i=1,il_g
           write(6,*) 'i,yy ',i,yy(5*il_g*il_g+il_g-il_g*i)   ! yy(il_g,il_g+1-i,4)
          enddo
          do j=1,il_g
           write(6,*) 'j,zz ',j,zz(5*il_g*il_g+il_g*j)      ! zz(il_g,j,5)
          enddo
        endif  ! ntest>0
      endif    ! .not.allocated

      qgg(1:ifull_g,1:klt)=0.
      tt(1:ifull_g,1:klt)=0.
      uu(1:ifull_g,1:klt)=0.
      vv(1:ifull_g,1:klt)=0.
      ww(1:ifull_g,1:klt)=0.
      psls(1:ifull_g)=0.
      qgg2(1:ifull_g,1:klt)=0.
      tt2(1:ifull_g,1:klt)=0.
      uu2(1:ifull_g,1:klt)=0.
      vv2(1:ifull_g,1:klt)=0.
      ww2(1:ifull_g,1:klt)=0.
      psls2(1:ifull_g)=0.
      sumwt(1:ifull_g)=1.e-20   ! for undefined panels
      emmin=sqrt(cutoff2)*ds/rearth
      write(6,*) 'schmidt,cutoff,kbotdav ',schmidt,sqrt(cutoff2),kbotdav 
      
      do j=1,il_g                ! doing x-filter on panels 1,2,4,5
       xxmax=xx(il_g*(6*il_g-1)+il_g+1-j)
       write(6,*) 'j,xxmax ',j,xxmax
       do n=1,4*il_g
        if(n<=il_g)iq=il_g*(il_g+j-1)+n                   ! panel 1
        if(n>il_g.and.n<=2*il_g)iq=il_g*(2*il_g+j-2)+n      ! panel 2
        if(n>2*il_g)iq=il_g*(2*il_g+n-1)+il_g+1-j           ! panel 4,5
        
        if (em_g(iq)>emmin) then
        
        do n1=n,4*il_g
!        following test shows on sx6 don't use "do n1=m+1,4*il_g"
!        if(n==4*il_g)write(6,*) 'problem for i,n,n1 ',i,n,n1
         if(n1<=il_g)iq1=il_g*(il_g+j-1)+n1               ! panel 1
         if(n1>il_g.and.n1<=2*il_g)iq1=il_g*(2*il_g+j-2)+n1 ! panel 2
         if(n1>2*il_g)iq1=il_g*(2*il_g+n1-1)+il_g+1-j       ! panel 4,5
         dist1=abs(xx(iq)-xx(iq1))
         dist=min(dist1,xxmax-dist1)
         wt=exp(-4.5*dist*dist*cutoff2)
         wt1=wt/em_g(iq1)
         wt=wt/em_g(iq)
         if(n==n1)wt1=0.  ! so as not to add in twice
c        if(iq==10345.or.iq1==10345)
c    &     write(6,*) 'iq,iq1,n,n1,xx,xx1,dist1,dist,wt,wt1 ',         
c    &              iq,iq1,n,n1,xx(iq),xx(iq1),dist1,dist,wt,wt1 
         sumwt(iq)=sumwt(iq)+wt1
         sumwt(iq1)=sumwt(iq1)+wt
!        producing "x-filtered" version of pslb-psl etc
c        psls(iq)=psls(iq)+wt1*(pslb(iq1)-psl(iq1))
c        psls(iq1)=psls(iq1)+wt*(pslb(iq)-psl(iq))
         psls(iq)=psls(iq)+wt1*psla(iq1)
         psls(iq1)=psls(iq1)+wt*psla(iq)
         do k=1,klt ! MJT nestin
          qgg(iq,k)=qgg(iq,k)+wt1*qa(iq1,k)
          qgg(iq1,k)=qgg(iq1,k)+wt*qa(iq,k)
          tt(iq,k)=tt(iq,k)+wt1*ta(iq1,k)
          tt(iq1,k)=tt(iq1,k)+wt*ta(iq,k)
          uu(iq,k)=uu(iq,k)+wt1*ua(iq1,k)
          uu(iq1,k)=uu(iq1,k)+wt*ua(iq,k)
          vv(iq,k)=vv(iq,k)+wt1*va(iq1,k)
          vv(iq1,k)=vv(iq1,k)+wt*va(iq,k)
          ww(iq,k)=ww(iq,k)+wt1*wa(iq1,k)
          ww(iq1,k)=ww(iq1,k)+wt*wa(iq,k)
         enddo  ! k loop
c        write(6,*) 'n,n1,dist,wt,wt1 ',n,n1,dist,wt,wt1
        enddo   ! n1 loop
        else
          sumwt(iq)=1.
        end if
       enddo    ! n loop
      enddo     ! j loop      
      if(nud_uv==-1)then
        do iq=1,ifull_g
         psls2(iq)=psls(iq)/sumwt(iq)
         do k=1,klt
          qgg2(iq,k)=qgg(iq,k)/sumwt(iq)
          tt2(iq,k)=tt(iq,k)/sumwt(iq)
          uu2(iq,k)=uu(iq,k)/sumwt(iq)
          vv2(iq,k)=vv(iq,k)/sumwt(iq)
          ww2(iq,k)=ww(iq,k)/sumwt(iq)
         enddo
        enddo
      else  ! original fast scheme
        do iq=1,ifull_g
         if(sumwt(iq)/=1.e-20)then
           psla(iq)=psls(iq)/sumwt(iq)
           do k=1,klt
            qa(iq,k)=qgg(iq,k)/sumwt(iq)
            ta(iq,k)=tt(iq,k)/sumwt(iq)
            ua(iq,k)=uu(iq,k)/sumwt(iq)
            va(iq,k)=vv(iq,k)/sumwt(iq)
            wa(iq,k)=ww(iq,k)/sumwt(iq)
           enddo
         endif  ! (sumwt(iq)/=1.e-20)
        enddo
      endif  ! (nud_uv==-1) .. else ..
      
      qgg(1:ifull_g,1:klt)=0.
      tt(1:ifull_g,1:klt)=0.
      uu(1:ifull_g,1:klt)=0.
      vv(1:ifull_g,1:klt)=0.
      ww(1:ifull_g,1:klt)=0.
      psls(1:ifull_g)=0.
      sumwt(1:ifull_g)=1.e-20   ! for undefined panels
      
      do i=1,il_g                ! doing y-filter on panels 0,1,3,4
       yymax=yy(il_g*(5*il_g-i+1))  
       do n=1,4*il_g
        if(n<=2*il_g)iq=il_g*(n-1)+i                      ! panel 0,1
        if(n>2*il_g.and.n<=3*il_g)iq=il_g*(4*il_g-i-2)+n      ! panel 3
        if(n>3*il_g)iq=il_g*(5*il_g-i-3)+n                  ! panel 4       
        if (em_g(iq)>emmin) then       
        do n1=n,4*il_g
         if(n1<=2*il_g)iq1=il_g*(n1-1)+i                  ! panel 0,1
         if(n1>2*il_g.and.n1<=3*il_g)iq1=il_g*(4*il_g-i-2)+n1 ! panel 3
         if(n1>3*il_g)iq1=il_g*(5*il_g-i-3)+n1              ! panel 4
         dist1=abs(yy(iq)-yy(iq1))
         dist=min(dist1,yymax-dist1)
         wt=exp(-4.5*dist*dist*cutoff2)
         wt1=wt/em_g(iq1)
         wt=wt/em_g(iq)
         if(n==n1)wt1=0.  ! so as not to add in twice
         sumwt(iq)=sumwt(iq)+wt1
         sumwt(iq1)=sumwt(iq1)+wt
!        producing "y-filtered" version of pslb-psl etc
         psls(iq)=psls(iq)+wt1*psla(iq1)
         psls(iq1)=psls(iq1)+wt*psla(iq)
         do k=1,klt
          qgg(iq,k)=qgg(iq,k)+wt1*qa(iq1,k)
          qgg(iq1,k)=qgg(iq1,k)+wt*qa(iq,k)
          tt(iq,k)=tt(iq,k)+wt1*ta(iq1,k)
          tt(iq1,k)=tt(iq1,k)+wt*ta(iq,k)
          uu(iq,k)=uu(iq,k)+wt1*ua(iq1,k)
          uu(iq1,k)=uu(iq1,k)+wt*ua(iq,k)
          vv(iq,k)=vv(iq,k)+wt1*va(iq1,k)
          vv(iq1,k)=vv(iq1,k)+wt*va(iq,k)
          ww(iq,k)=ww(iq,k)+wt1*wa(iq1,k)
          ww(iq1,k)=ww(iq1,k)+wt*wa(iq,k)
         enddo  ! k loop
        enddo   ! n1 loop
        else
          sumwt(iq)=1.
        end if
       enddo    ! n loop
      enddo     ! i loop
      if(nud_uv==-1)then
        do iq=1,ifull_g
         psls2(iq)=psls2(iq)+psls(iq)/sumwt(iq)
         do k=1,klt
          qgg2(iq,k)=qgg2(iq,k)+qgg(iq,k)/sumwt(iq)
          tt2(iq,k)=tt2(iq,k)+tt(iq,k)/sumwt(iq)
          uu2(iq,k)=uu2(iq,k)+uu(iq,k)/sumwt(iq)
          vv2(iq,k)=vv2(iq,k)+vv(iq,k)/sumwt(iq)
          ww2(iq,k)=ww2(iq,k)+ww(iq,k)/sumwt(iq)
         enddo
        enddo
      else  ! original fast scheme
        do iq=1,ifull_g
         if(sumwt(iq)/=1.e-20)then
           psla(iq)=psls(iq)/sumwt(iq)
           do k=1,klt
            qa(iq,k)=qgg(iq,k)/sumwt(iq)
            ta(iq,k)=tt(iq,k)/sumwt(iq)
            ua(iq,k)=uu(iq,k)/sumwt(iq)
            va(iq,k)=vv(iq,k)/sumwt(iq)
            wa(iq,k)=ww(iq,k)/sumwt(iq)
           enddo
         endif  ! (sumwt(iq)/=1.e-20)
        enddo
      endif  ! (nud_uv==-1) .. else ..

      if(mbd>=0) then
       qgg(1:ifull_g,1:klt)=0.
       tt(1:ifull_g,1:klt)=0.
       uu(1:ifull_g,1:klt)=0.
       vv(1:ifull_g,1:klt)=0.
       ww(1:ifull_g,1:klt)=0.
       psls(1:ifull_g)=0.
       sumwt(1:ifull_g)=1.e-20   ! for undefined panels
    
       do j=1,il_g                ! doing "z"-filter on panels 0,2,3,5
        zzmax=zz(5*il_g*il_g+il_g*j)
        write(6,*) 'j,zzmax ',j,zzmax
        do n=1,4*il_g
         if(n<=il_g)iq=il_g*(j-1)+n                     ! panel 0
         if(n>il_g.and.n<=3*il_g)iq=il_g*(il_g+n-1)+il_g+1-j  ! panel 2,3
         if(n>3*il_g)iq=il_g*(5*il_g+j-4)+n               ! panel 5
        
         if (em_g(iq)>emmin) then
        
         do n1=n,4*il_g
          if(n1<=il_g)iq1=il_g*(j-1)+n1                     ! panel 0
          if(n1>il_g.and.n1<=3*il_g)iq1=il_g*(il_g+n1-1)+il_g+1-j ! panel 2,3
          if(n1>3*il_g)iq1=il_g*(5*il_g+j-4)+n1               ! panel 5
          dist1=abs(zz(iq)-zz(iq1))
          dist=min(dist1,zzmax-dist1)
          wt=exp(-4.5*dist*dist*cutoff2)
          wt1=wt/em_g(iq1)
          wt=wt/em_g(iq)
          if(n==n1)wt1=0.  ! so as not to add in twice
          sumwt(iq)=sumwt(iq)+wt1
          sumwt(iq1)=sumwt(iq1)+wt
!         producing "z"-filtered version of pslb-psl etc
          psls(iq)=psls(iq)+wt1*psla(iq1)
          psls(iq1)=psls(iq1)+wt*psla(iq)
          do k=1,klt ! MJT nestin
           qgg(iq,k)=qgg(iq,k)+wt1*qa(iq1,k)
           qgg(iq1,k)=qgg(iq1,k)+wt*qa(iq,k)
           tt(iq,k)=tt(iq,k)+wt1*ta(iq1,k)
           tt(iq1,k)=tt(iq1,k)+wt*ta(iq,k)
           uu(iq,k)=uu(iq,k)+wt1*ua(iq1,k)
           uu(iq1,k)=uu(iq1,k)+wt*ua(iq,k)
           vv(iq,k)=vv(iq,k)+wt1*va(iq1,k)
           vv(iq1,k)=vv(iq1,k)+wt*va(iq,k)
           ww(iq,k)=ww(iq,k)+wt1*wa(iq1,k)
           ww(iq1,k)=ww(iq1,k)+wt*wa(iq,k)
          enddo  ! k loop
         enddo   ! n1 loop
         else
           sumwt(iq)=1.
         end if
        enddo    ! n loop
       enddo     ! j loop      
      if(nud_uv==-1)then
        write(6,*) 'in nestinb nud_uv ',nud_uv
        do iq=1,ifull_g
         psls2(iq)=psls2(iq)+psls(iq)/sumwt(iq)
         do k=1,klt ! MJT nestin
          qgg2(iq,k)=qgg2(iq,k)+qgg(iq,k)/sumwt(iq)
          tt2(iq,k)=tt2(iq,k)+tt(iq,k)/sumwt(iq)
          uu2(iq,k)=uu2(iq,k)+uu(iq,k)/sumwt(iq)
          vv2(iq,k)=vv2(iq,k)+vv(iq,k)/sumwt(iq)
          ww2(iq,k)=ww2(iq,k)+ww(iq,k)/sumwt(iq)
         enddo
        enddo
        psla(1:ifull_g)=.5*psls2(1:ifull_g)
        qa(1:ifull_g,1:klt)=.5*qgg2(1:ifull_g,1:klt)
        ta(1:ifull_g,1:klt)=.5*tt2(1:ifull_g,1:klt)
        ua(1:ifull_g,1:klt)=.5*uu2(1:ifull_g,1:klt)
        va(1:ifull_g,1:klt)=.5*vv2(1:ifull_g,1:klt)
        wa(1:ifull_g,1:klt)=.5*ww2(1:ifull_g,1:klt)
      else  ! original fast scheme
        write(6,*) 'in nestinb  nud_uv ',nud_uv
        do iq=1,ifull_g
         if(sumwt(iq)/=1.e-20)then
           psla(iq)=psls(iq)/sumwt(iq)
           do k=1,klt ! MJT nestin
            qa(iq,k)=qgg(iq,k)/sumwt(iq)
            ta(iq,k)=tt(iq,k)/sumwt(iq)
            ua(iq,k)=uu(iq,k)/sumwt(iq)
            va(iq,k)=vv(iq,k)/sumwt(iq)
            wa(iq,k)=ww(iq,k)/sumwt(iq)
           enddo
         endif  ! (sumwt(iq)/=1.e-20)
        enddo
      endif  ! (nud_uv==-1) .. else ..
      end if ! (mbd>=0)

      return
      end subroutine fastspec

      !---------------------------------------------------------------------------------
      ! Slow 2D spectral downscaling - MPI version
      subroutine slowspecmpi(cin,psls,uu,vv,ww,tt,qgg,lblock,klt)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      use xyzinfo_m         ! Grid coordinate arrays
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration

      integer, intent(in) :: klt
      integer i,j,n,iq,iqg,k,ierr
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,klt), intent(inout) :: uu,vv
      real, dimension(ifull_g,klt), intent(inout) :: ww
      real, dimension(ifull_g,klt), intent(inout) :: tt,qgg
      real, dimension(ifull_g) :: r
      real, dimension(ifull) :: pd
      real, dimension(ifull,klt) :: ud,vd,wd,td,qd
      real cq,psum
      logical, intent(in) :: lblock

      cq=sqrt(4.5)*cin

      if (myid==0.and.nmaxpr==1) then
        write(6,*) "Start 2D filter"
      end if
   
      do n=1,npan
        do j=1,jpan
          do i=1,ipan
            iqg=indg(i,j,n)
            iq=indp(i,j,n)
            r(:)=x_g(iqg)*x_g(:)+y_g(iqg)*y_g(:)+z_g(iqg)*z_g(:)
            r(:)=acos(max(min(r(:),1.),-1.))
            r(:)=exp(-(cq*r(:))**2)/(em_g(:)**2)
            psum=sum(r(:))
            if (nud_p>0.and.lblock) then
              pd(iq)=sum(r(:)*psls(:))/psum
            end if
            if (nud_uv>0) then
              do k=1,klt
                ud(iq,k)=sum(r(:)*uu(:,k))/psum
                vd(iq,k)=sum(r(:)*vv(:,k))/psum
                wd(iq,k)=sum(r(:)*ww(:,k))/psum
              end do        
            end if
            if (nud_t>0) then
              do k=1,klt
                td(iq,k)=sum(r(:)*tt(:,k))/psum
              end do
            end if
            if (nud_q>0) then
              do k=1,klt
                qd(iq,k)=sum(r(:)*qgg(:,k))/psum
              end do
            end if
          end do
        end do
      end do
 
      if (nud_p>0.and.lblock) then
        psls(1:ifull)=pd
      end if
      if (nud_uv>0) then
        uu(1:ifull,:)=ud(:,:)
        vv(1:ifull,:)=vd(:,:)
        ww(1:ifull,:)=wd(:,:)
      end if
      if (nud_t>0) then
        tt(1:ifull,:)=td(:,:)
      end if
      if (nud_q>0) then
        qgg(1:ifull,:)=qd(:,:)
      end if
 
      if (myid == 0.and.nmaxpr==1) write(6,*) "End 2D filter"

      return
      end subroutine slowspecmpi
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      ! Four pass spectral downscaling (symmetric)
      ! Used when mod(6,nproc)/=0 and mod(nproc,6)/=0 since it is slower than specfastmpi
      subroutine fourspecmpi(myid,cin,psls,uu,vv,ww,tt,qgg,lblock,klt)
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: myid,klt
      integer pn,px,hproc,mproc,ns,ne,npta
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,klt), intent(inout) :: uu,vv
      real, dimension(ifull_g,klt), intent(inout) :: ww
      real, dimension(ifull_g,klt), intent(inout) :: tt,qgg
      logical, intent(in) :: lblock
      
      npta=1                              ! number of panels per processor
      mproc=nproc                         ! number of processors per panel
      pn=0                                ! start panel
      px=5                                ! end panel
      hproc=0                             ! host processor for panel
      call procdiv(ns,ne,il_g,nproc,myid) ! number of rows per processor

      call spechost(mproc,hproc,npta,pn,px,ns,ne,cin,psls,uu,vv,
     &       ww,tt,qgg,lblock,klt)
          
      return
      end subroutine fourspecmpi
      !---------------------------------------------------------------------------------


      !---------------------------------------------------------------------------------
      ! Four pass spectral downscaling (symmetric)
      ! MPI optimised for magic processor numbers 1,2,3,6,12,18,24,30,36,...
      ! (only works for mod(6,nproc)==0 or mod(nproc,6)==0)
      subroutine specfastmpi(myid,cin,psls,uu,vv,ww,tt,qgg,lblock,klt)
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: myid,klt
      integer pn,px,hproc,mproc,ns,ne,npta
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,klt), intent(inout) :: uu,vv
      real, dimension(ifull_g,klt), intent(inout) :: ww
      real, dimension(ifull_g,klt), intent(inout) :: tt,qgg
      logical, intent(in) :: lblock
      
      npta=max(6/nproc,1)                       ! number of panels per processor
      mproc=max(nproc/6,1)                      ! number of processors per panel
      pn=myid*npta/mproc                        ! start panel
      px=(myid+mproc)*npta/mproc-1              ! end panel
      hproc=pn*mproc/npta                       ! host processor for panel
      call procdiv(ns,ne,il_g,mproc,myid-hproc) ! number of rows per processor

      if (npta==1) then
        call spechost_n(mproc,hproc,pn,ns,ne,cin,psls,uu,vv,
     &         ww,tt,qgg,lblock,klt)
      else
        call spechost(mproc,hproc,npta,pn,px,ns,ne,cin,psls,uu,vv,
     &         ww,tt,qgg,lblock,klt)
      end if

      return
      end subroutine specfastmpi
      !---------------------------------------------------------------------------------


      !---------------------------------------------------------------------------------
      ! This is the main routine for the scale-selective filter
      subroutine spechost(mproc,hproc,npta,pn,px,ns,ne,cin,psls,
     &                    uu,vv,ww,tt,qgg,lblock,klt)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: mproc,hproc,npta,pn,px,ns,ne,klt
      integer k,ppass,iy,ppn,ppx,nne,nns,iproc,ierr
      integer n,a,b,c,til,lsize
      integer :: itag=0
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,klt), intent(inout) :: uu,vv
      real, dimension(ifull_g,klt), intent(inout) :: ww
      real, dimension(ifull_g,klt), intent(inout) :: tt,qgg
      real, dimension(ifull_g) :: qp,qsum,zp
      real, dimension(ifull_g,klt) :: qu,qv,qw,qt,qq
      real, dimension(ifull_g,klt) :: zu,zv,zw,zt,zq
      real, dimension(ifull_g*klt) :: dd
      real cq
      logical, intent(in) :: lblock
      
      til=il_g*il_g
      cq=sqrt(4.5)*cin ! filter length scale

      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "Start 1D filter"

      qsum(:)=1./(em_g(:)*em_g(:))
      if (nud_p>0.and.lblock) then
        psls(:)=psls(:)*qsum(:)
      end if
      if (nud_uv>0) then
        do k=1,klt
          uu(:,k)=uu(:,k)*qsum(:)
          vv(:,k)=vv(:,k)*qsum(:)
          ww(:,k)=ww(:,k)*qsum(:)
        end do
      end if
      if (nud_t>0) then
        do k=1,klt
          tt(:,k)=tt(:,k)*qsum(:)
        end do
      end if
      if (nud_q>0) then
        do k=1,klt
          qgg(:,k)=qgg(:,k)*qsum(:)
        end do
      end if

      do ppass=pn,px

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

        ! computations for the local processor group
        if (mod(il_g,mproc)==0) then
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "Symmetric decomposition"
          end if
          lsize=ifull_g/mproc
          call speclocal_s(mproc,hproc,ns,ne,cq,ppass,qsum,qp,
     &           qu,qv,qw,qt,qq,lsize,lblock,klt)        
        else
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "Asymmetric decomposition"
          end if
          call speclocal(mproc,hproc,ns,ne,cq,ppass,qsum,qp,
     &           qu,qv,qw,qt,qq,lblock,klt)
        end if
        
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
      
      if (nud_p>0.and.lblock) then
        psls(:)=zp(:)
      end if
      if (nud_uv>0) then
        uu(:,:)=zu(:,:)
        vv(:,:)=zv(:,:)
        ww(:,:)=zw(:,:)
      end if
      if (nud_t>0) then
        tt(:,:)=zt(:,:)
      end if
      if (nud_q>0) then
        qgg(:,:)=zq(:,:)
      end if

      return
      end subroutine spechost
      !---------------------------------------------------------------------------------

      ! This version is for one panel per processor (reduced memory)
      subroutine spechost_n(mproc,hproc,pn,ns,ne,cin,psls,
     &                    uu,vv,ww,tt,qgg,lblock,klt)

      use cc_mpi            ! CC MPI routines
      use map_m             ! Grid map arrays
      
      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: mproc,hproc,pn,ns,ne,klt
      integer k,iy,ppn,nne,nns,iproc,ierr
      integer n,a,b,c,til,lsize
      integer :: itag=0
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,klt), intent(inout) :: uu,vv
      real, dimension(ifull_g,klt), intent(inout) :: ww
      real, dimension(ifull_g,klt), intent(inout) :: tt,qgg
      real, dimension(ifull_g) :: qsum
      real, dimension(ifull_g*klt) :: dd
      real cq
      logical, intent(in) :: lblock
      
      til=il_g*il_g
      cq=sqrt(4.5)*cin ! filter length scale

      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "Start 1D filter"

      qsum(:)=1./(em_g(:)*em_g(:))
      if (nud_p>0.and.lblock) then
        psls(:)=psls(:)*qsum(:)
      end if
      if (nud_uv>0) then
        do k=1,klt
          uu(:,k)=uu(:,k)*qsum(:)
          vv(:,k)=vv(:,k)*qsum(:)
          ww(:,k)=ww(:,k)*qsum(:)
        end do
      end if
      if (nud_t>0) then
        do k=1,klt
          tt(:,k)=tt(:,k)*qsum(:)
        end do
      end if
      if (nud_q>0) then
        do k=1,klt
          qgg(:,k)=qgg(:,k)*qsum(:)
        end do
      end if

      ! computations for the local processor group
      if (mod(il_g,mproc)==0) then
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "Symmetric decomposition"
        end if
        lsize=ifull_g/mproc
        call speclocal_s(mproc,hproc,ns,ne,cq,pn,qsum,psls,
     &         uu,vv,ww,tt,qgg,lsize,lblock,klt)        
      else
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "Asymmetric decomposition"
        end if
        call speclocal(mproc,hproc,ns,ne,cq,pn,qsum,psls,
     &         uu,vv,ww,tt,qgg,lblock,klt)
      end if
        
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
        
      if (myid==0.and.nmaxpr==1) write(6,*) "End 1D filter"

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
      
      return
      end subroutine spechost_n
      
      !---------------------------------------------------------------------------------
      ! This code runs between the local processors
      ! Code was moved to this subroutine to help the compiler vectorise the code
      ! This version handles uneven distribution of rows across processors
      subroutine speclocal(mproc,hproc,ns,ne,cq,ppass,qsum,
     &             qp,qu,qv,qw,qt,qq,lblock,klt)

      use cc_mpi            ! CC MPI routines
      use xyzinfo_m         ! Grid coordinate arrays

      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: mproc,hproc,ns,ne,ppass,klt
      integer j,k,n,ipass,kpass,iy
      integer iproc,ierr
      integer nne,nns,me
      integer a,b,c,d
      integer :: itag=0
      integer, dimension(4*il_g,il_g,0:3) :: igrd
      integer, dimension(0:3) :: maps
      integer, parameter, dimension(2:3) :: kn=(/0,3/)
      integer, parameter, dimension(2:3) :: kx=(/2,3/)
      real, intent(in) :: cq
      real, dimension(ifull_g), intent(inout) :: qp,qsum
      real, dimension(ifull_g,klt), intent(inout) :: qu,qv
      real, dimension(ifull_g,klt), intent(inout) :: qw
      real, dimension(ifull_g,klt), intent(inout) :: qt,qq
      real, dimension(4*il_g,klt) :: pu,pv,pw,pt,pq
      real, dimension(4*il_g,klt) :: au,av,aw,at,aq
      real, dimension(4*il_g) :: pp,ap,psum,asum,ra,xa,ya,za
      real, dimension(ifull_g*klt) :: dd
      logical, intent(in) :: lblock
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      
      do ipass=0,3
        me=maps(ipass)
        call getiqa(igrd(:,:,ipass),me,ipass,ppass,il_g)

        if (ipass==3) then
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
                  dd(n+a*j+d)=qsum(igrd(n,j,ipass))
                end do
              end do
              call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
              if(nud_p>0.and.lblock)then
                do j=nns,nne
                  do n=1,me
                    dd(n+a*j+d)=qp(igrd(n,j,ipass))
                  end do
                end do
                call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
              end if
              iy=me*(nne-nns+1)*klt
              b=me*(nne-nns+1)
              d=-me*(nns+(nne-nns+1))
              if(nud_uv>0)then
                do k=1,klt
                  do j=nns,nne
                    do n=1,me
                      dd(n+a*j+b*k+d)=qu(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
                do k=1,klt
                  do j=nns,nne
                    do n=1,me
                      dd(n+a*j+b*k+d)=qv(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
                do k=1,klt
                  do j=nns,nne
                    do n=1,me
                      dd(n+a*j+b*k+d)=qw(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
              end if
              if(nud_t>0)then
                do k=1,klt
                  do j=nns,nne
                    do n=1,me
                      dd(n+a*j+b*k+d)=qt(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
              end if
              if(nud_q>0)then
                do k=1,klt
                  do j=nns,nne
                    do n=1,me
                      dd(n+a*j+b*k+d)=qq(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(dd(1:iy),iproc,itag,comm_world)
              end if
            end do
          else
            iy=me*(ne-ns+1)
            a=me
            d=-me*ns
            call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
            do j=ns,ne
              do n=1,me
                qsum(igrd(n,j,ipass))=dd(n+a*j+d)
              end do
            end do
            if(nud_p>0.and.lblock)then
              call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
              do j=ns,ne
                do n=1,me
                  qp(igrd(n,j,ipass))=dd(n+a*j+d)
                end do
              end do
            endif
            iy=me*(ne-ns+1)*klt
            b=me*(ne-ns+1)
            d=-me*(ns+(ne-ns+1))
            if(nud_uv>0)then
              call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
              do k=1,klt
                do j=ns,ne
                  do n=1,me
                    qu(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
              do k=1,klt
                do j=ns,ne
                  do n=1,me
                    qv(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
              do k=1,klt
                do j=ns,ne
                  do n=1,me
                    qw(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
            if(nud_t>0)then
              call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
              do k=1,klt
                do j=ns,ne
                  do n=1,me
                    qt(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
            if(nud_q>0)then
              call ccmpi_recv(dd(1:iy),hproc,itag,comm_world)
              do k=1,klt
                do j=ns,ne
                  do n=1,me
                    qq(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
          end if
        end if

        if (myid==0.and.nmaxpr==1) write(6,*) "Start convolution"

        do j=ns,ne
          asum(1:me)=qsum(igrd(1:me,j,ipass))
          xa(1:me)=x_g(igrd(1:me,j,ipass))
          ya(1:me)=y_g(igrd(1:me,j,ipass))
          za(1:me)=z_g(igrd(1:me,j,ipass))
          if (nud_p>0.and.lblock) then
            ap(1:me)=qp(igrd(1:me,j,ipass))
          end if
          if (nud_uv>0) then
            au(1:me,:)=qu(igrd(1:me,j,ipass),:)
            av(1:me,:)=qv(igrd(1:me,j,ipass),:)
            aw(1:me,:)=qw(igrd(1:me,j,ipass),:)
          end if
          if (nud_t>0) then
            at(1:me,:)=qt(igrd(1:me,j,ipass),:)
          end if
          if (nud_q>0) then
            aq(1:me,:)=qq(igrd(1:me,j,ipass),:)
          end if
          do n=1,il_g
            ra(1:me)=xa(n)*xa(1:me)+ya(n)*ya(1:me)+za(n)*za(1:me)
            ra(1:me)=acos(max(min(ra(1:me),1.),-1.))
            ra(1:me)=exp(-(cq*ra(1:me))**2)
            ! can also use the lines below which integrate the gaussian
            ! analytically over the length element (but slower)
            !ra(1)=2.*erf(cq*0.5*(ds/rearth)
            !ra(2:me)=erf(cq*(ra(2:me)+0.5*(ds/rearth)))  ! redefine ra(:) as wgt(:)
     &      !        -erf(cq*(ra(2:me)-0.5*(ds/rearth)))  ! (correct units are 1/cq)
            psum(n)=sum(ra(1:me)*asum(1:me))
            if (nud_p>0.and.lblock) then
              pp(n)=sum(ra(1:me)*ap(1:me))
            end if
            if (nud_uv>0) then
              do k=1,klt
                pu(n,k)=sum(ra(1:me)*au(1:me,k))
                pv(n,k)=sum(ra(1:me)*av(1:me,k))
                pw(n,k)=sum(ra(1:me)*aw(1:me,k))
              end do
            end if
            if (nud_t>0) then
              do k=1,klt
                pt(n,k)=sum(ra(1:me)*at(1:me,k))
              end do
            end if
            if (nud_q>0) then
              do k=1,klt
                pq(n,k)=sum(ra(1:me)*aq(1:me,k))
              end do
            end if
          end do
          qsum(igrd(1:il_g,j,ipass))=psum(1:il_g)
          if (nud_p>0.and.lblock) then
            qp(igrd(1:il_g,j,ipass))=pp(1:il_g)
          end if
          if (nud_uv>0) then
            qu(igrd(1:il_g,j,ipass),:)=pu(1:il_g,:)
            qv(igrd(1:il_g,j,ipass),:)=pv(1:il_g,:)
            qw(igrd(1:il_g,j,ipass),:)=pw(1:il_g,:)
          end if
          if (nud_t>0) then
            qt(igrd(1:il_g,j,ipass),:)=pt(1:il_g,:)
          end if
          if (nud_q>0) then
            qq(igrd(1:il_g,j,ipass),:)=pq(1:il_g,:)
          end if
        end do

        if (myid==0.and.nmaxpr==1) write(6,*) "End convolution"

        if (ipass==2.or.ipass==3) then
          itag=itag+1
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "Send arrays to local host"
          end if
          if (myid==hproc) then
            do iproc=hproc+1,mproc+hproc-1
              call procdiv(nns,nne,il_g,mproc,iproc-hproc)
              if (nns>nne) exit
              iy=il_g*(nne-nns+1)*(kx(ipass)-kn(ipass)+1)
              a=il_g
              c=il_g*(nne-nns+1)
              d=-il_g*(nns+(nne-nns+1)*kn(ipass))
              call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
              do kpass=kn(ipass),kx(ipass)
                do j=nns,nne
                  do n=1,il_g
                    qsum(igrd(n,j,kpass))=dd(n+a*j+c*kpass+d)
                  end do
                end do
              end do
              if(nud_p>0.and.lblock)then
                call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
                do kpass=kn(ipass),kx(ipass)
                  do j=nns,nne
                    do n=1,il_g
                      qp(igrd(n,j,kpass))=dd(n+a*j+c*kpass+d)
                    end do
                  end do
                end do
              end if
              iy=il_g*(nne-nns+1)*klt
     &           *(kx(ipass)-kn(ipass)+1)
              b=il_g*(nne-nns+1)
              c=il_g*(nne-nns+1)*klt
              d=-il_g*(nns+(nne-nns+1)
     &          *(1+klt*kn(ipass)))
              if(nud_uv>0)then
                call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qu(igrd(n,j,kpass),k)
     &                    =dd(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
                call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qv(igrd(n,j,kpass),k)
     &                    =dd(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
                call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qw(igrd(n,j,kpass),k)
     &                    =dd(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
              end if
              if(nud_t>0)then
                call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qt(igrd(n,j,kpass),k)
     &                    =dd(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
              end if
              if(nud_q>0)then
                call ccmpi_recv(dd(1:iy),iproc,itag,comm_world)
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qq(igrd(n,j,kpass),k)
     &                    =dd(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
              end if
            end do
          else
            iy=il_g*(ne-ns+1)*(kx(ipass)-kn(ipass)+1)
            a=il_g
            c=il_g*(ne-ns+1)
            d=-il_g*(ns+(ne-ns+1)*kn(ipass))
            do kpass=kn(ipass),kx(ipass)
              do j=ns,ne
                do n=1,il_g
                  dd(n+a*j+c*kpass+d)=qsum(igrd(n,j,kpass))
                end do
              end do
            end do
            call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
            if(nud_p>0.and.lblock)then
              do kpass=kn(ipass),kx(ipass)
                do j=ns,ne
                  do n=1,il_g
                    dd(n+a*j+c*kpass+d)=qp(igrd(n,j,kpass))
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
            end if
            iy=il_g*(ne-ns+1)*klt
     &         *(kx(ipass)-kn(ipass)+1)
            b=il_g*(ne-ns+1)
            c=il_g*(ne-ns+1)*klt
            d=-il_g*(ns+(ne-ns+1)
     &        *(1+klt*kn(ipass)))
            if(nud_uv>0)then
              do k=1,klt
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      dd(n+a*j+b*k+c*kpass+d)=qu(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
              do k=1,klt
                do kpass=kn(ipass),kx(ipass)
                   do j=ns,ne
                    do n=1,il_g
                      dd(n+a*j+b*k+c*kpass+d)=qv(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
              do k=1,klt
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      dd(n+a*j+b*k+c*kpass+d)=qw(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
            end if
            if(nud_t>0)then
              do k=1,klt
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      dd(n+a*j+b*k+c*kpass+d)=qt(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
            end if
            if(nud_q>0)then
              do k=1,klt
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      dd(n+a*j+b*k+c*kpass+d)=qq(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(dd(1:iy),hproc,itag,comm_world)
            end if
          end if
        end if
          
      end do
      
      return  
      end subroutine speclocal

      !---------------------------------------------------------------------------------
      ! This code runs between the local processors
      ! Code was moved to this subroutine to help the compiler vectorise the code
      ! This version only handles an even distribution of rows across processors
      subroutine speclocal_s(mproc,hproc,ns,ne,cq,ppass,qsum,
     &             qp,qu,qv,qw,qt,qq,lsize,lblock,klt)

      use cc_mpi            ! CC MPI routines
      use xyzinfo_m         ! Grid coordinate arrays

      implicit none
      
      include 'newmpar.h'   ! Grid parameters
      include 'parm.h'      ! Model configuration
      
      integer, intent(in) :: mproc,hproc,ns,ne,ppass,lsize,klt
      integer j,k,n,ipass,kpass,iy,ix,me,a,b,c
      integer iproc,ierr
      integer :: itag=0
      integer, save :: comm_face
      integer, dimension(4*il_g,il_g,0:3) :: igrd
      integer, dimension(0:3) :: maps
      integer, parameter, dimension(2:3) :: kn=(/0,3/)
      integer, parameter, dimension(2:3) :: kx=(/2,3/)
      real, intent(in) :: cq
      real, dimension(ifull_g), intent(inout) :: qp,qsum
      real, dimension(ifull_g,klt), intent(inout) :: qu,qv
      real, dimension(ifull_g,klt), intent(inout) :: qw
      real, dimension(ifull_g,klt), intent(inout) :: qt,qq
      real, dimension(4*il_g,klt) :: pu,pv,pw,pt,pq
      real, dimension(4*il_g,klt) :: au,av,aw,at,aq
      real, dimension(4*il_g) :: pp,ap,psum,asum,ra,xa,ya,za
      real, dimension(ifull_g*klt) :: dd
      real, dimension(lsize*klt) :: ee
      logical, intent(in) :: lblock
      logical, save :: first=.true.
      
      if (first) then
        call ccmpi_commsplit(comm_face,comm_world,hproc,myid-hproc)
        first=.false.
      end if
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      dd=0.
      ee=0.
      
      do ipass=0,3
        me=maps(ipass)
        call getiqa(igrd(:,:,ipass),me,ipass,ppass,il_g)

        if (ipass==3) then
          itag=itag+1
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "Receive arrays from local host"
          end if
          iy=me*il_g
          ix=iy/mproc
          if (myid==hproc) then
            do j=ne+1,il_g
              do n=1,me
                dd(n+me*(j-1))=qsum(igrd(n,j,ipass))
              end do
            end do
          end if
          call ccmpi_scatterx(dd(1:iy),ee(1:ix),0,comm_face)
          if (myid/=hproc) then
            do j=ns,ne
              do n=1,me
                qsum(igrd(n,j,ipass))=ee(n+me*(j-ns))
              end do
            end do
          end if
          if(nud_p>0.and.lblock)then
            if (myid==hproc) then
              do j=ne+1,il_g
                do n=1,me
                  dd(n+me*(j-1))=qp(igrd(n,j,ipass))
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do n=1,me
                  qp(igrd(n,j,ipass))=ee(n+me*(j-ns))
                end do
              end do
            end if
          end if
          iy=me*il_g*klt
          ix=iy/mproc
          a=me*klt
          if(nud_uv>0)then
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt
                  do n=1,me
                    dd(n+me*(k-1)+a*(j-1))
     &                =qu(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do n=1,me
                    qu(igrd(n,j,ipass),k)
     &                =ee(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt
                  do n=1,me
                    dd(n+me*(k-1)+a*(j-1))
     &                =qv(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do n=1,me
                    qv(igrd(n,j,ipass),k)
     &                =ee(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt
                  do n=1,me
                    dd(n+me*(k-1)+a*(j-1))
     &                =qw(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do n=1,me
                    qw(igrd(n,j,ipass),k)
     &                =ee(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
          end if
          if(nud_t>0)then
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt
                  do n=1,me
                    dd(n+me*(k-1)+a*(j-1))
     &                =qt(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do n=1,me
                    qt(igrd(n,j,ipass),k)
     &                =ee(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
          end if
          if(nud_q>0)then
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt
                  do n=1,me
                    dd(n+me*(k-1)+a*(j-1))
     &                =qq(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do n=1,me
                    qq(igrd(n,j,ipass),k)
     &                =ee(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
          end if
        end if

        if (myid==0.and.nmaxpr==1) write(6,*) "Start convolution"

        do j=ns,ne
          asum(1:me)=qsum(igrd(1:me,j,ipass))
          xa(1:me)=x_g(igrd(1:me,j,ipass))
          ya(1:me)=y_g(igrd(1:me,j,ipass))
          za(1:me)=z_g(igrd(1:me,j,ipass))
          if (nud_p>0.and.lblock) then
            ap(1:me)=qp(igrd(1:me,j,ipass))
          end if
          if (nud_uv>0) then
            au(1:me,:)=qu(igrd(1:me,j,ipass),:)
            av(1:me,:)=qv(igrd(1:me,j,ipass),:)
            aw(1:me,:)=qw(igrd(1:me,j,ipass),:)
          end if
          if (nud_t>0) then
            at(1:me,:)=qt(igrd(1:me,j,ipass),:)
          end if
          if (nud_q>0) then
            aq(1:me,:)=qq(igrd(1:me,j,ipass),:)
          end if
          do n=1,il_g
            ra(1:me)=xa(n)*xa(1:me)+ya(n)*ya(1:me)+za(n)*za(1:me)
            ra(1:me)=acos(max(min(ra(1:me),1.),-1.))
            ra(1:me)=exp(-(cq*ra(1:me))**2)
            ! can also use the lines below which integrate the gaussian
            ! analytically over the length element (but slower)
            !ra(1)=2.*erf(cq*0.5*(ds/rearth)
            !ra(2:me)=erf(cq*(ra(2:me)+0.5*(ds/rearth)))  ! redefine ra(:) as wgt(:)
     &      !        -erf(cq*(ra(2:me)-0.5*(ds/rearth)))  ! (correct units are 1/cq)
            psum(n)=sum(ra(1:me)*asum(1:me))
            if (nud_p>0.and.lblock) then
              pp(n)=sum(ra(1:me)*ap(1:me))
            end if
            if (nud_uv>0) then
              do k=1,klt
                pu(n,k)=sum(ra(1:me)*au(1:me,k))
                pv(n,k)=sum(ra(1:me)*av(1:me,k))
                pw(n,k)=sum(ra(1:me)*aw(1:me,k))
              end do
            end if
            if (nud_t>0) then
              do k=1,klt
                pt(n,k)=sum(ra(1:me)*at(1:me,k))
              end do
            end if
            if (nud_q>0) then
              do k=1,klt
                pq(n,k)=sum(ra(1:me)*aq(1:me,k))
              end do
            end if
          end do
          qsum(igrd(1:il_g,j,ipass))=psum(1:il_g)
          if (nud_p>0.and.lblock) then
            qp(igrd(1:il_g,j,ipass))=pp(1:il_g)
          end if
          if (nud_uv>0) then
            qu(igrd(1:il_g,j,ipass),:)=pu(1:il_g,:)
            qv(igrd(1:il_g,j,ipass),:)=pv(1:il_g,:)
            qw(igrd(1:il_g,j,ipass),:)=pw(1:il_g,:)
          end if
          if (nud_t>0) then
            qt(igrd(1:il_g,j,ipass),:)=pt(1:il_g,:)
          end if
          if (nud_q>0) then
            qq(igrd(1:il_g,j,ipass),:)=pq(1:il_g,:)
          end if
        end do

        if (myid==0.and.nmaxpr==1) write(6,*) "End convolution"

        if (ipass==2.or.ipass==3) then
          itag=itag+1
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "Send arrays to local host"
          end if
          iy=il_g*il_g*(kx(ipass)-kn(ipass)+1)
          ix=iy/mproc
          a=il_g*(kx(ipass)-kn(ipass)+1)
          c=kn(ipass)
          if (myid/=hproc) then
            do j=ns,ne
              do kpass=kn(ipass),kx(ipass)
                do n=1,il_g
                  ee(n+il_g*(kpass-c)+a*(j-ns))
     &              =qsum(igrd(n,j,kpass))
                end do
              end do
            end do
          end if
          call ccmpi_gatherx(dd(1:iy),ee(1:ix),0,comm_face)
          if (myid==hproc) then
            do j=ne+1,il_g
              do kpass=kn(ipass),kx(ipass)
                do n=1,il_g
                  qsum(igrd(n,j,kpass))
     &              =dd(n+il_g*(kpass-c)+a*(j-1))
                end do
              end do
            end do
          end if
          if(nud_p>0.and.lblock)then
            if (myid/=hproc) then
              do j=ns,ne
                do kpass=kn(ipass),kx(ipass)
                  do n=1,il_g
                    ee(n+il_g*(kpass-c)+a*(j-ns))
     &                =qp(igrd(n,j,kpass))
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do kpass=kn(ipass),kx(ipass)
                  do n=1,il_g
                    qp(igrd(n,j,kpass))
     &                =dd(n+il_g*(kpass-c)+a*(j-1))
                  end do
                end do
              end do
            end if
          end if
          iy=il_g*il_g*klt
     &       *(kx(ipass)-kn(ipass)+1)
          ix=iy/mproc
          a=il_g*(kx(ipass)-kn(ipass)+1)
          b=il_g*(kx(ipass)-kn(ipass)+1)*klt
          c=kn(ipass)
          if(nud_uv>0)then
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      ee(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qu(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt           
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qu(igrd(n,j,kpass),k)
     &                  =dd(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      ee(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qv(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt          
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qv(igrd(n,j,kpass),k)
     &                  =dd(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      ee(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qw(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt   
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qw(igrd(n,j,kpass),k)
     &                  =dd(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
          end if
          if(nud_t>0)then
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      ee(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qt(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt           
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qt(igrd(n,j,kpass),k)
     &                  =dd(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
          end if
          if(nud_q>0)then
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,klt
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      ee(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qq(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(dd(1:iy),ee(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,klt        
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qq(igrd(n,j,kpass),k)
     &                  =dd(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
          end if
        end if
          
      end do
      
      return  
      end subroutine speclocal_s

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
      real, dimension(ifull_g,1) :: diffh_g
      real, dimension(ifull_g,kblock) :: diff_g,diffs_g
      real, dimension(ifull_g,kblock) :: diffu_g,diffv_g,diffw_g
      real, dimension(ifull,kblock) :: diff
      real, dimension(ifull_g) :: x_g,xx_g
      real, dimension(ifull) :: old,oldt
      real, dimension(ifull,ktopmlo:kbotmlo) :: rho,nsq
      logical disflag,lblock
      integer, parameter :: tempfix=1 ! delta temp (0=linear, 1=buoyancy)
      real, parameter :: rho0=1030.   ! linear density offset
      real, parameter :: a0=-0.3      ! linear density temp gradient
      real, parameter :: miss = 999999.
      
      
      diff_g=miss
      diffs_g=miss
      diffu_g=miss
      diffv_g=miss
      diffw_g=miss
      diffh_g=miss

      kc=min(kbotmlo,ktopmlo+wl-1)

      if (nud_sfh/=0) then
        diff(:,1)=miss
        old=sfh
        call mloexport(4,old,0,0)
        where (.not.land)
          diff(:,1)=sfh-old
        end where
        call ccmpi_gatherall(diff(:,1),diffh_g(:,1))
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
          diff=miss
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=sstb(:,ka)
            call mloexport(0,old,k,0)
            where (.not.land)
              diff(:,kb)=sstb(:,ka)-old
            end where
          end do
          call ccmpi_gatherall(diff(:,1:klt),diff_g(:,1:klt))
        end if

        if (nud_sss/=0) then
          diff=miss
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=sssb(:,ka)
            call mloexport(1,old,k,0)
            where (.not.land)
              diff(:,kb)=sssb(:,ka)-old
            end where
          end do
          call ccmpi_gatherall(diff(:,1:klt),diffs_g(:,1:klt))
        end if

        if (nud_ouv/=0) then
          diff=miss
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=suvb(:,ka,1)
            call mloexport(2,old,k,0)
            where (.not.land)
              diff(:,kb)=suvb(:,ka,1)-old
            end where
          end do
          call ccmpi_gatherall(diff(:,1:klt),diffu_g(:,1:klt))
          diff=miss
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=suvb(:,ka,2)
            call mloexport(3,old,k,0)
            where (.not.land)
              diff(:,kb)=suvb(:,ka,2)-old
            end where
          end do
          call ccmpi_gatherall(diff(:,1:klt),diffv_g(:,1:klt))
          do k=1,klt
            x_g=diffu_g(:,k)
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
        end if

        if ((nud_uv/=9.and.abs(nmlo)/=1).or.namip/=0) then
          call mlofilterfast(diff_g(:,1:klt),diffs_g(:,1:klt),
     &                       diffu_g(:,1:klt),diffv_g(:,1:klt),
     &                       diffw_g(:,1:klt),diffh_g(:,1),
     &                       miss,lblock,klt)
          disflag=.true.
        else
          call mlofilter(diff_g(:,1:klt),diffs_g(:,1:klt),
     &                   diffu_g(:,1:klt),diffv_g(:,1:klt),
     &                   diffw_g(:,1:klt),diffh_g(:,1),
     &                   miss,lblock,klt)
          disflag=.false.
        end if

        if (myid==0) then
          write(6,*) "Distribute data for MLO filter ",kbb
        end if

        if (nud_sst/=0) then
          if (disflag) then
            if (myid == 0) then
              call ccmpi_distribute(diff(:,1:klt), diff_g(:,1:klt))
            else
              call ccmpi_distribute(diff(:,1:klt))
            end if
          else
            diff(:,1:klt)=diff_g(1:ifull,1:klt)
          end if
          ! correct temp pertubation to minimise change in buoyancy
          if (tempfix==1.and.kc==1) then
            if (ktopmlo/=1) then
              write(6,*) "ERROR: nud_sst with SST input"
              write(6,*) "requires ktopmlo=1"
              stop
            end if
            old=293.
            do k=ktopmlo,kbotmlo
              call mloexport(0,old,k,0)
              rho(:,k)=rho0+a0*old ! linear approximation to density
            end do
            do k=ktopmlo,kbotmlo-1
              !nsq=-2.*grav*(rho(:,k)-rho(:,k+1))/((dep(:,k+1)-dep(:,k))*(rho(:,k)+rho(:,k+1)))
              nsq(:,k)=-(rho(:,k)-rho(:,k+1))/(rho(:,k)+rho(:,k+1))
              !nsq(:,k)=max(nsq(:,k),0.)
            end do
            ! nsq(:,k-1)=-(a0*(oldt-old))/(2.*rho0+a0*(oldt+old))
            ! old*(1.-nsq(:,k-1))=oldt*(1.+nsq(:,k-1))+2.*rho0/a0*nsq(:,k-1)
            call mloexport(0,old,1,0)
            old=old+diff(:,1)*10./real(mloalpha)
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
              old=old+diff(:,kb)*10./real(mloalpha)
              old=max(old,271.)
              call mloimport(0,old,k,0)
            end do
            if (klx==kc) then
              do k=kc+1,kbotmlo
                old=sstb(:,ka)
                call mloexport(0,old,k,0)
                old=old+diff(:,kb)*10./real(mloalpha) ! kb saved from above loop
                old=max(old,271.)	  
                call mloimport(0,old,k,0)
              end do
            end if
          end if
        end if

        if (nud_sss/=0) then
          if (disflag) then
            if (myid == 0) then
              call ccmpi_distribute(diff(:,1:klt), diffs_g(:,1:klt))
            else
              call ccmpi_distribute(diff(:,1:klt))
            end if
          else
            diff(:,1:klt)=diffs_g(1:ifull,1:klt)
          end if
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=sssb(:,ka)
            call mloexport(1,old,k,0)
            old=old+diff(:,kb)*10./real(mloalpha)
            old=max(old,0.)
            call mloimport(1,old,k,0)
          end do
          if (klx==kc) then
            do k=kc+1,kbotmlo
              old=sssb(:,ka)
              call mloexport(1,old,k,0)
              old=old+diff(:,kb)*10./real(mloalpha) ! kb saved from above loop
              old=max(old,0.)
              call mloimport(1,old,k,0)
            end do
          end if
        end if

        if (nud_ouv/=0) then
          if (disflag) then
            if (myid == 0) then
              do k=1,klt
                x_g=ax_g*diffu_g(:,k)+ay_g*diffv_g(:,k)
     &            +az_g*diffw_g(:,k)
                xx_g=bx_g*diffu_g(:,k)+by_g*diffv_g(:,k)
     &            +bz_g*diffw_g(:,k)
                diffu_g(:,k)=x_g
                diffv_g(:,k)=xx_g
              end do
              call ccmpi_distribute(diff(:,1:klt), diffu_g(:,1:klt))
            else
              call ccmpi_distribute(diff(:,1:klt))
            end if
          else
            do k=1,klt
              x_g(1:ifull)=ax*diffu_g(1:ifull,k)+ay*diffv_g(1:ifull,k)
     &          +az*diffw_g(1:ifull,k)
              xx_g(1:ifull)=bx*diffu_g(1:ifull,k)+by*diffv_g(1:ifull,k)
     &          +bz*diffw_g(1:ifull,k)
              diffu_g(1:ifull,k)=x_g(1:ifull)
              diffv_g(1:ifull,k)=xx_g(1:ifull)
            end do
            diff(:,1:klt)=diffu_g(1:ifull,1:klt)
          end if
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=suvb(:,ka,1)
            call mloexport(2,old,k,0)
            old=old+diff(:,kb)*10./real(mloalpha)
            call mloimport(2,old,k,0)
            if (allocated(oldu1)) then
              oldu1(:,k)=oldu1(:,k)+diff(:,kb)*10./real(mloalpha)
              oldu2(:,k)=oldu2(:,k)+diff(:,kb)*10./real(mloalpha)
            end if
          end do
          if (klx==kc) then
            do k=kc+1,kbotmlo
              old=suvb(:,ka,1)
              call mloexport(2,old,k,0)
              old=old+diff(:,kb)*10./real(mloalpha) ! kb saved from above loop
              call mloimport(2,old,k,0)
              if (allocated(oldu1)) then
                oldu1(:,k)=oldu1(:,k)+diff(:,kb)*10./real(mloalpha)
                oldu2(:,k)=oldu2(:,k)+diff(:,kb)*10./real(mloalpha)
              end if
            end do
          end if
          if (disflag) then
            if (myid == 0) then
              call ccmpi_distribute(diff(:,1:klt), diffv_g(:,1:klt))
            else
              call ccmpi_distribute(diff(:,1:klt))
            end if
          else
            diff(:,1:klt)=diffv_g(1:ifull,1:klt)
          end if
          do k=kln,klx
            ka=min(wl,k)
            kb=k-kln+1
            old=suvb(:,ka,2)
            call mloexport(3,old,k,0)
            old=old+diff(:,kb)*10./real(mloalpha)
            call mloimport(3,old,k,0)
            if (allocated(oldv1)) then
              oldv1(:,k)=oldv1(:,k)+diff(:,kb)*10./real(mloalpha)
              oldv2(:,k)=oldv2(:,k)+diff(:,kb)*10./real(mloalpha)
            end if
          end do
          if (klx==kc) then
            do k=kc+1,kbotmlo
              old=suvb(:,ka,2)
              call mloexport(3,old,k,0)
              old=old+diff(:,kb)*10./real(mloalpha)
              call mloimport(3,old,k,0)
              if (allocated(oldv1)) then
                oldv1(:,k)=oldv1(:,k)+diff(:,kb)*10./real(mloalpha)
                oldv2(:,k)=oldv2(:,k)+diff(:,kb)*10./real(mloalpha)
              end if
            end do
          end if
        end if
      
      end do
     
      if (nud_sfh/=0) then
        if (disflag) then
          if (myid == 0) then
            call ccmpi_distribute(diff(:,1:1), diffh_g(:,1:1))
          else
            call ccmpi_distribute(diff(:,1:1))
          end if
        else
          diff(:,1)=diffh_g(1:ifull,1)
        end if
        old=sfh
        call mloexport(4,old,0,0)
        old=old+diff(:,1)*10./real(mloalpha)
        call mloimport(4,old,0,0)
      end if

      return
      end subroutine mlofilterhub
      
      !---------------------------------------------------------------------------------
      ! 2D Filter for MLO 
      subroutine mlofilter(diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                     diffh_g,miss,lblock,kd)

      use cc_mpi                  ! CC MPI routines

      implicit none

      include 'newmpar.h'         ! Grid parameters
      include 'parm.h'            ! Model configuration

      integer, intent(in) :: kd
      integer ierr
      real, intent(in) :: miss
      real, dimension(ifull_g,1), intent(inout) :: diffh_g
      real, dimension(ifull_g,kd), intent(inout) :: diff_g,diffs_g
      real, dimension(ifull_g,kd), intent(inout) :: diffu_g,diffv_g
      real, dimension(ifull_g,kd), intent(inout) :: diffw_g
      logical, intent(in) :: lblock
      logical, dimension(ifull_g) :: landg

      if (nud_sst/=0) then
        landg=abs(diff_g(:,1)-miss)<0.1
      end if

      if (nud_sss/=0) then
        landg=abs(diffs_g(:,1)-miss)<0.1
      end if

      if (nud_ouv/=0) then
        landg=abs(diffw_g(:,1)-miss)<0.1
      end if

      if (nud_sfh/=0.and.lblock) then
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

      return
      end subroutine mlofilter

      subroutine mlofilterhost(diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                         diffh_g,kd,miss,landg,lblock)

      use cc_mpi             ! CC MPI routines
      use map_m              ! Grid map arrays
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
      real, dimension(ifull,kd) :: dd,dds,ddu,ddv,ddw
      logical, intent(in) :: lblock
      logical, dimension(ifull_g), intent(in) :: landg

      ! eventually will be replaced with mbd once full ocean coupling is complete
      cq=sqrt(4.5)*.1*real(max(nud_sst,nud_sss,nud_ouv,nud_sfh,mbd))
     &   /(pi*schmidt)

      mm=1./(em_g*em_g)
      nn=0.
      dd=0.
      dds=0.
      ddu=0.
      ddv=0.
      ddw=0.
      ddh=0.
      where(.not.landg)
        nn=mm
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
                  ddu(iqq,k)=sum(rr(:)*diffu_g(:,k))/nsum
                  ddv(iqq,k)=sum(rr(:)*diffv_g(:,k))/nsum
                  ddw(iqq,k)=sum(rr(:)*diffw_g(:,k))/nsum
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
      subroutine mlofilterfast(diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                         diffh_g,miss,lblock,kd)

      use cc_mpi                  ! CC MPI routines

      implicit none

      include 'newmpar.h'         ! Grid parameters
      include 'const_phys.h'      ! Physical constants
      include 'parm.h'            ! Model configuration
      include 'parmgeom.h'        ! Coordinate data

      integer, intent(in) :: kd
      integer pn,px,hproc,mproc,ns,ne,npta
      real, intent(in) :: miss      
      real, dimension(ifull_g,1), intent(inout) :: diffh_g
      real, dimension(ifull_g,kd), intent(inout) :: diff_g,diffs_g
      real, dimension(ifull_g,kd), intent(inout) :: diffu_g,diffv_g
      real, dimension(ifull_g,kd), intent(inout) :: diffw_g
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
        px=(myid+mproc)*npta/mproc-1              ! end panel
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

      if (pn==px) then
        call mlospechost_n(mproc,hproc,pn,ns,ne,cq,
     &                   diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                   diffh_g,miss,lblock,kd)
      else
        call mlospechost(mproc,hproc,npta,pn,px,ns,ne,cq,
     &                   diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                   diffh_g,miss,lblock,kd)
      end if

      if (myid==0.and.nmaxpr==1) then
        write(6,*) "MLO end 1D filter"
      end if

      return
      end subroutine mlofilterfast

      subroutine mlospechost(mproc,hproc,npta,pn,px,ns,ne,cq,
     &                       diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                       diffh_g,miss,lblock,kd)

      use cc_mpi             ! CC MPI routines
      use map_m              ! Grid map arrays
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: mproc,hproc,npta,pn,px,ns,ne,kd
      integer ppass,iy,ppn,ppx,nne,nns,iproc,ierr
      integer n,a,b,c,k,til,lsize
      integer :: itag=0
      real, intent(in) :: cq,miss
      real, dimension(ifull_g,1), intent(inout) :: diffh_g
      real, dimension(ifull_g,kd), intent(inout) :: diff_g,diffs_g
      real, dimension(ifull_g,kd), intent(inout) :: diffu_g,diffv_g
      real, dimension(ifull_g,kd), intent(inout) :: diffw_g
      real, dimension(ifull_g) :: qsum,rsum,zph,qph
      real, dimension(ifull_g,kd) :: zp,zps,zpu,zpv,zpw
      real, dimension(ifull_g,kd) :: qp,qps,qpu,qpv,qpw
      real, dimension(ifull_g*kd) :: zz
      logical, intent(in) :: lblock
      logical, dimension(ifull_g) :: landg
      
      til=il_g*il_g 

      if (nud_sst/=0) then
        landg=abs(diff_g(:,1)-miss)<0.1
      end if
      if (nud_sss/=0) then
        landg=abs(diffs_g(:,1)-miss)<0.1
      end if
      if (nud_ouv/=0) then
        landg=abs(diffw_g(:,1)-miss)<0.1
      end if
      if (nud_sfh/=0.and.lblock) then
        landg=abs(diffh_g(:,1)-miss)<0.1
      end if
      
      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "MLO Start 1D filter"

      zp=0.
      zps=0.
      zpu=0.
      zpv=0.
      zpw=0.
      zph=0.

      do ppass=pn,px

        rsum=0.
        qsum(:)=1./(em_g(:)*em_g(:))
        where(.not.landg) ! land/sea mask
          rsum(:)=qsum(:)
        end where

        if (nud_sst/=0) then
          do k=1,kd
            qp(:,k)=diff_g(:,k)*rsum
          end do
        end if
        if (nud_sss/=0) then
          do k=1,kd
            qps(:,k)=diffs_g(:,k)*rsum
          end do
        end if 
        if (nud_ouv/=0) then
          do k=1,kd
            qpu(:,k)=diffu_g(:,k)*rsum
            qpv(:,k)=diffv_g(:,k)*rsum
            qpw(:,k)=diffw_g(:,k)*rsum
          end do
        end if
        if (nud_sfh/=0.and.lblock) then
          qph=diffh_g(:,1)*rsum
        end if

        ! computations for the local processor group
        if (mod(il_g,mproc)==0) then
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "MLO Symmetric decomposition"
          end if
          lsize=ifull_g/mproc
          call mlospeclocal_s(mproc,hproc,ns,ne,cq,ppass,qsum,
     &                      qp,qps,qpu,qpv,qpw,qph,kd,lsize,lblock)
        else
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "MLO Asymmetric decomposition"
          end if
          call mlospeclocal(mproc,hproc,ns,ne,cq,ppass,qsum,
     &                      qp,qps,qpu,qpv,qpw,qph,kd,lblock)
        end if
        
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
      
      diff_g(:,:)=zp(:,:)
      diffs_g(:,:)=zps(:,:)
      diffu_g(:,:)=zpu(:,:)
      diffv_g(:,:)=zpv(:,:)
      diffw_g(:,:)=zpw(:,:)
      diffh_g(:,1)=zph(:)

      return
      end subroutine mlospechost
      !---------------------------------------------------------------------------------

      subroutine mlospechost_n(mproc,hproc,pn,ns,ne,cq,
     &                       diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                       diffh_g,miss,lblock,kd)

      use cc_mpi             ! CC MPI routines
      use map_m              ! Grid map arrays
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: mproc,hproc,pn,ns,ne,kd
      integer iy,ppn,nne,nns,iproc,ierr
      integer n,a,b,c,k,til,lsize
      integer :: itag=0
      real, intent(in) :: cq,miss
      real, dimension(ifull_g,1), intent(inout) :: diffh_g
      real, dimension(ifull_g,kd), intent(inout) :: diff_g,diffs_g
      real, dimension(ifull_g,kd), intent(inout) :: diffu_g,diffv_g
      real, dimension(ifull_g,kd), intent(inout) :: diffw_g
      real, dimension(ifull_g) :: qsum,rsum
      real, dimension(ifull_g*kd) :: zz
      logical, intent(in) :: lblock
      logical, dimension(ifull_g) :: landg
      
      til=il_g*il_g 

      if (nud_sst/=0) then
        landg=abs(diff_g(:,1)-miss)<0.1
      end if
      if (nud_sss/=0) then
        landg=abs(diffs_g(:,1)-miss)<0.1
      end if
      if (nud_ouv/=0) then
        landg=abs(diffw_g(:,1)-miss)<0.1
      end if
      if (nud_sfh/=0.and.lblock) then
        landg=abs(diffh_g(:,1)-miss)<0.1
      end if
      
      if (ns>ne) return
      if (myid==0.and.nmaxpr==1) write(6,*) "MLO Start 1D filter"

      rsum=0.
      qsum(:)=1./(em_g(:)*em_g(:))
      where(.not.landg) ! land/sea mask
        rsum(:)=qsum(:)
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

      ! computations for the local processor group
      if (mod(il_g,mproc)==0) then
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "MLO Symmetric decomposition"
        end if
        lsize=ifull_g/mproc
        call mlospeclocal_s(mproc,hproc,ns,ne,cq,pn,qsum,
     &                    diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                    diffh_g,kd,lsize,lblock)
      else
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "MLO Asymmetric decomposition"
        end if
        call mlospeclocal(mproc,hproc,ns,ne,cq,pn,qsum,
     &                    diff_g,diffs_g,diffu_g,diffv_g,diffw_g,
     &                    diffh_g,kd,lblock)
      end if
        
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

      if (myid==0.and.nmaxpr==1) write(6,*) "MLO End 1D filter"

      itag=itag+1
      if (myid == 0) then
        if (nmaxpr==1) write(6,*) "MLO Receive arrays from all proc"
        do iproc=mproc,nproc-1,mproc
          ppn=iproc/mproc
          iy=til*kd
          a=til
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
        a=til
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

      return
      end subroutine mlospechost_n

      
      !---------------------------------------------------------------------------------
      ! This version is for asymmetric decomposition
      subroutine mlospeclocal(mproc,hproc,ns,ne,cq,ppass,
     &             qsum,qp,qps,qpu,qpv,qpw,qph,kd,lblock)

      use cc_mpi             ! CC MPI routines
      use xyzinfo_m          ! Grid coordinate arrays
     
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: mproc,hproc,ns,ne,ppass,kd
      integer j,n,ipass,kpass,iy
      integer iproc,ierr
      integer nne,nns,me
      integer a,b,c,d,k
      integer :: itag=0
      integer, dimension(4*il_g,il_g,0:3) :: igrd
      integer, dimension(0:3) :: maps
      integer, parameter, dimension(2:3) :: kn=(/0,3/)
      integer, parameter, dimension(2:3) :: kx=(/2,3/)
      real, intent(in) :: cq
      real, dimension(ifull_g), intent(inout) :: qph
      real, dimension(ifull_g,kd), intent(inout) :: qp,qps
      real, dimension(ifull_g,kd), intent(inout) :: qpu,qpv
      real, dimension(ifull_g,kd), intent(inout) :: qpw
      real, dimension(ifull_g), intent(inout) :: qsum
      real, dimension(4*il_g) :: rr,xa,ya,za
      real, dimension(4*il_g) :: asum,psum
      real, dimension(4*il_g) :: aph,pph
      real, dimension(4*il_g,kd) :: ap,aps,pp,pps
      real, dimension(4*il_g,kd) :: apu,apv,apw,ppu,ppv,ppw
      real, dimension(ifull_g*kd) :: zz
      logical, intent(in) :: lblock
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      
      do ipass=0,3
        me=maps(ipass)
        call getiqa(igrd(:,:,ipass),me,ipass,ppass,il_g)

        if (ipass==3) then
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
                  zz(n+a*j+d)=qsum(igrd(n,j,ipass))
                end do
              end do
              call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
              if (nud_sfh/=0.and.lblock) then
                do j=nns,nne
                  do n=1,me
                    zz(n+a*j+d)=qph(igrd(n,j,ipass))
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
                      zz(n+a*j+b*k+d)=qp(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
              end if
              if (nud_sss/=0) then
                do k=1,kd
                  do j=nns,nne
                    do n=1,me
                      zz(n+a*j+b*k+d)=qps(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
              end if
              if (nud_ouv/=0) then
                do k=1,kd
                  do j=nns,nne
                    do n=1,me
                      zz(n+a*j+b*k+d)=qpu(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
                do k=1,kd
                  do j=nns,nne
                    do n=1,me
                      zz(n+a*j+b*k+d)=qpv(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
                do k=1,kd
                  do j=nns,nne
                    do n=1,me
                      zz(n+a*j+b*k+d)=qpw(igrd(n,j,ipass),k)
                    end do
                  end do
                end do
                call ccmpi_ssend(zz(1:iy),iproc,itag,comm_world)
              end if
            end do
          else
            iy=me*(ne-ns+1)
            a=me
            d=-me*ns
            call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
            do j=ns,ne
              do n=1,me
                qsum(igrd(n,j,ipass))=zz(n+a*j+d)
              end do
            end do
            if (nud_sfh/=0.and.lblock) then
              call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
              do j=ns,ne
                do n=1,me
                  qph(igrd(n,j,ipass))=zz(n+a*j+d)
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
                    qp(igrd(n,j,ipass),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
            end if
            if (nud_sss/=0) then
              call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
              do k=1,kd
                do j=ns,ne
                  do n=1,me
                    qps(igrd(n,j,ipass),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do  
            end if
            if (nud_ouv/=0) then
              call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
              do k=1,kd
                do j=ns,ne
                  do n=1,me
                    qpu(igrd(n,j,ipass),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
              do k=1,kd
                do j=ns,ne
                  do n=1,me
                    qpv(igrd(n,j,ipass),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do
              call ccmpi_recv(zz(1:iy),hproc,itag,comm_world)
              do k=1,kd
                do j=ns,ne
                  do n=1,me
                    qpw(igrd(n,j,ipass),k)=zz(n+a*j+b*k+d)
                  end do
                end do
              end do  
            end if          
          end if
        end if

        if (myid==0.and.nmaxpr==1) write(6,*) "MLO start conv"

        do j=ns,ne
          xa(1:me)=x_g(igrd(1:me,j,ipass))
          ya(1:me)=y_g(igrd(1:me,j,ipass))
          za(1:me)=z_g(igrd(1:me,j,ipass))
          asum(1:me)=qsum(igrd(1:me,j,ipass))
          if (nud_sst/=0) then
            do k=1,kd
              ap(1:me,k)=qp(igrd(1:me,j,ipass),k)
            end do
          end if
          if (nud_sss/=0) then
            do k=1,kd
              aps(1:me,k)=qps(igrd(1:me,j,ipass),k)
            end do
          end if
          if (nud_ouv/=0) then
            do k=1,kd
              apu(1:me,k)=qpu(igrd(1:me,j,ipass),k)
              apv(1:me,k)=qpv(igrd(1:me,j,ipass),k)
              apw(1:me,k)=qpw(igrd(1:me,j,ipass),k)
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            aph(1:me)=qph(igrd(1:me,j,ipass))
          end if
          do n=1,il_g
            rr(1:me)=xa(n)*xa(1:me)+ya(n)*ya(1:me)+za(n)*za(1:me)
            rr(1:me)=acos(max(min(rr(1:me),1.),-1.))
            rr(1:me)=exp(-(cq*rr(1:me))**2)
            psum(n)=sum(rr(1:me)*asum(1:me))
            if (nud_sst/=0) then
              do k=1,kd
                pp(n,k)=sum(rr(1:me)*ap(1:me,k))
              end do
            end if
            if (nud_sss/=0) then
              do k=1,kd
                pps(n,k)=sum(rr(1:me)*aps(1:me,k))
              end do
            end if
            if (nud_ouv/=0) then
              do k=1,kd
                ppu(n,k)=sum(rr(1:me)*apu(1:me,k))
                ppv(n,k)=sum(rr(1:me)*apv(1:me,k))
                ppw(n,k)=sum(rr(1:me)*apw(1:me,k))
              end do
            end if
            if (nud_sfh/=0.and.lblock) then
              pph(n)=sum(rr(1:me)*aph(1:me))
            end if
          end do
          qsum(igrd(1:il_g,j,ipass))=psum(1:il_g)
          if (nud_sst/=0) then
            do k=1,kd
              qp(igrd(1:il_g,j,ipass),k)=pp(1:il_g,k)
            end do
          end if
          if (nud_sss/=0) then
            do k=1,kd
              qps(igrd(1:il_g,j,ipass),k)=pps(1:il_g,k)
            end do
          end if
          if (nud_ouv/=0) then
            do k=1,kd
              qpu(igrd(1:il_g,j,ipass),k)=ppu(1:il_g,k)
              qpv(igrd(1:il_g,j,ipass),k)=ppv(1:il_g,k)
              qpw(igrd(1:il_g,j,ipass),k)=ppw(1:il_g,k)
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            qph(igrd(1:il_g,j,ipass))=pph(1:il_g)
          end if
        end do

        if (myid==0.and.nmaxpr==1) write(6,*) "MLO end conv"

        if (ipass==2.or.ipass==3) then
          itag=itag+1
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "MLO Send arrays to local host"
          end if
          if (myid==hproc) then
            do iproc=hproc+1,mproc+hproc-1
              call procdiv(nns,nne,il_g,mproc,iproc-hproc)
              if (nns>nne) exit
              iy=il_g*(nne-nns+1)*(kx(ipass)-kn(ipass)+1)
              a=il_g
              c=il_g*(nne-nns+1)
              d=-il_g*(nns+(nne-nns+1)*kn(ipass))
              call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
              do kpass=kn(ipass),kx(ipass)
                do j=nns,nne
                  do n=1,il_g
                    qsum(igrd(n,j,kpass))=zz(n+a*j+c*kpass+d)
                  end do
                end do
              end do
              if (nud_sfh/=0.and.lblock) then
                call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
                do kpass=kn(ipass),kx(ipass)
                  do j=nns,nne
                    do n=1,il_g
                      qph(igrd(n,j,kpass))=zz(n+a*j+c*kpass+d)
                    end do
                  end do
                end do
              end if
              iy=il_g*(nne-nns+1)*kd*(kx(ipass)-kn(ipass)+1)
              b=il_g*(nne-nns+1)
              c=il_g*(nne-nns+1)*kd
              d=-il_g*(nns+(nne-nns+1)*(1+kd*kn(ipass)))
              if (nud_sst/=0) then
                call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qp(igrd(n,j,kpass),k)=zz(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
              end if
              if (nud_sss/=0) then
                call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qps(igrd(n,j,kpass),k)=zz(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
              end if
              if (nud_ouv/=0) then
                call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qpu(igrd(n,j,kpass),k)=zz(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
                call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qpv(igrd(n,j,kpass),k)=zz(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
                call ccmpi_recv(zz(1:iy),iproc,itag,comm_world)
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qpw(igrd(n,j,kpass),k)=zz(n+a*j+b*k+c*kpass+d)
                      end do
                    end do
                  end do
                end do
              end if
            end do
          else
            iy=il_g*(ne-ns+1)*(kx(ipass)-kn(ipass)+1)
            a=il_g
            c=il_g*(ne-ns+1)
            d=-il_g*(ns+(ne-ns+1)*kn(ipass))
            do kpass=kn(ipass),kx(ipass)
              do j=ns,ne
                do n=1,il_g
                  zz(n+a*j+c*kpass+d)=qsum(igrd(n,j,kpass))
                end do
              end do
            end do
            call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
            if (nud_sfh/=0.and.lblock) then
              do kpass=kn(ipass),kx(ipass)
                do j=ns,ne
                  do n=1,il_g
                    zz(n+a*j+c*kpass+d)=qph(igrd(n,j,kpass))
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
            end if
            iy=il_g*(ne-ns+1)*kd*(kx(ipass)-kn(ipass)+1)
            b=il_g*(ne-ns+1)
            c=il_g*(ne-ns+1)*kd
            d=-il_g*(ns+(ne-ns+1)*(1+kd*kn(ipass)))
            if (nud_sst/=0) then
              do k=1,kd
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      zz(n+a*j+b*k+c*kpass+d)=qp(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
            end if
            if (nud_sss/=0) then
              do k=1,kd
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      zz(n+a*j+b*k+c*kpass+d)=qps(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
            end if
            if (nud_ouv/=0) then
              do k=1,kd
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      zz(n+a*j+b*k+c*kpass+d)=qpu(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
              do k=1,kd
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      zz(n+a*j+b*k+c*kpass+d)=qpv(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
              do k=1,kd
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      zz(n+a*j+b*k+c*kpass+d)=qpw(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
              call ccmpi_ssend(zz(1:iy),hproc,itag,comm_world)
            end if
          end if
        end if
          
      end do
      
      return  
      end subroutine mlospeclocal

      !---------------------------------------------------------------------------------
      ! This version is for asymmetric decomposition
      subroutine mlospeclocal_s(mproc,hproc,ns,ne,cq,ppass,
     &             qsum,qp,qps,qpu,qpv,qpw,qph,kd,lsize,lblock)

      use cc_mpi             ! CC MPI routines
      use xyzinfo_m          ! Grid coordinate arrays
     
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      include 'parm.h'       ! Model configuration
      
      integer, intent(in) :: mproc,hproc,ns,ne,ppass,kd
      integer, intent(in) :: lsize
      integer j,n,ipass,kpass,iy,ix
      integer iproc,ierr
      integer nne,nns,me
      integer a,b,c,k
      integer :: itag=0
      integer, save :: comm_face
      integer, dimension(4*il_g,il_g,0:3) :: igrd
      integer, dimension(0:3) :: maps
      integer, parameter, dimension(2:3) :: kn=(/0,3/)
      integer, parameter, dimension(2:3) :: kx=(/2,3/)
      real, intent(in) :: cq
      real, dimension(ifull_g), intent(inout) :: qph
      real, dimension(ifull_g,kd), intent(inout) :: qp,qps
      real, dimension(ifull_g,kd), intent(inout) :: qpu,qpv
      real, dimension(ifull_g,kd), intent(inout) :: qpw
      real, dimension(ifull_g), intent(inout) :: qsum
      real, dimension(4*il_g) :: rr,xa,ya,za
      real, dimension(4*il_g) :: asum,psum
      real, dimension(4*il_g) :: aph,pph
      real, dimension(4*il_g,kd) :: ap,aps,pp,pps
      real, dimension(4*il_g,kd) :: apu,apv,apw,ppu,ppv,ppw
      real, dimension(ifull_g*kd) :: zz
      real, dimension(lsize*kd) :: yy
      logical, intent(in) :: lblock
      logical, save :: first=.true.

      if (first) then
        call ccmpi_commsplit(comm_face,comm_world,hproc,myid-hproc)
        first=.false.
      end if
      
      maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      zz=0.
      yy=0.
      
      do ipass=0,3
        me=maps(ipass)
        call getiqa(igrd(:,:,ipass),me,ipass,ppass,il_g)

        if (ipass==3) then
          itag=itag+1
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "MLO Recieve arrays from local host"
          end if
          iy=me*il_g
          ix=iy/mproc
          if (myid==hproc) then
            do j=ne+1,il_g
              do n=1,me
                zz(n+me*(j-1))=qsum(igrd(n,j,ipass))
              end do
            end do
          end if
          call ccmpi_scatterx(zz(1:iy),yy(1:ix),0,comm_face)
          if (myid/=hproc) then
            do j=ns,ne
              do n=1,me
                qsum(igrd(n,j,ipass))=yy(n+me*(j-ns))
              end do
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            if (myid==hproc) then
              do j=ne+1,il_g
                do n=1,me
                  zz(n+me*(j-1))=qph(igrd(n,j,ipass))
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do n=1,me
                  qph(igrd(n,j,ipass))=yy(n+me*(j-ns))
                end do
              end do
            end if
          end if
          iy=me*il_g*kd
          ix=iy/mproc
          a=me*kd
          if (nud_sst/=0) then
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd
                  do n=1,me
                    zz(n+me*(k-1)+a*(j-1))
     &                =qp(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do n=1,me
                    qp(igrd(n,j,ipass),k)
     &                =yy(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
          end if
          if (nud_sss/=0) then
             if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd
                  do n=1,me
                    zz(n+me*(k-1)+a*(j-1))
     &                =qps(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do n=1,me
                    qps(igrd(n,j,ipass),k)
     &                =yy(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
          end if
          if (nud_ouv/=0) then
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd
                  do n=1,me
                    zz(n+me*(k-1)+a*(j-1))
     &                =qpu(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do n=1,me
                    qpu(igrd(n,j,ipass),k)
     &                =yy(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd
                  do n=1,me
                    zz(n+me*(k-1)+a*(j-1))
     &                =qpv(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do n=1,me
                    qpv(igrd(n,j,ipass),k)
     &                =yy(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd
                  do n=1,me
                    zz(n+me*(k-1)+a*(j-1))
     &                =qpw(igrd(n,j,ipass),k)
                  end do
                end do
              end do
            end if
            call ccmpi_scatterx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do n=1,me
                    qpw(igrd(n,j,ipass),k)
     &                =yy(n+me*(k-1)+a*(j-ns))
                  end do
                end do
              end do
            end if
          end if
        end if

        if (myid==0.and.nmaxpr==1) write(6,*) "MLO start conv"

        do j=ns,ne
          xa(1:me)=x_g(igrd(1:me,j,ipass))
          ya(1:me)=y_g(igrd(1:me,j,ipass))
          za(1:me)=z_g(igrd(1:me,j,ipass))
          asum(1:me)=qsum(igrd(1:me,j,ipass))
          if (nud_sst/=0) then
            do k=1,kd
              ap(1:me,k)=qp(igrd(1:me,j,ipass),k)
            end do
          end if
          if (nud_sss/=0) then
            do k=1,kd
              aps(1:me,k)=qps(igrd(1:me,j,ipass),k)
            end do
          end if
          if (nud_ouv/=0) then
            do k=1,kd
              apu(1:me,k)=qpu(igrd(1:me,j,ipass),k)
              apv(1:me,k)=qpv(igrd(1:me,j,ipass),k)
              apw(1:me,k)=qpw(igrd(1:me,j,ipass),k)
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            aph(1:me)=qph(igrd(1:me,j,ipass))
          end if
          do n=1,il_g
            rr(1:me)=xa(n)*xa(1:me)+ya(n)*ya(1:me)+za(n)*za(1:me)
            rr(1:me)=acos(max(min(rr(1:me),1.),-1.))
            rr(1:me)=exp(-(cq*rr(1:me))**2)
            psum(n)=sum(rr(1:me)*asum(1:me))
            if (nud_sst/=0) then
              do k=1,kd
                pp(n,k)=sum(rr(1:me)*ap(1:me,k))
              end do
            end if
            if (nud_sss/=0) then
              do k=1,kd
                pps(n,k)=sum(rr(1:me)*aps(1:me,k))
              end do
            end if
            if (nud_ouv/=0) then
              do k=1,kd
                ppu(n,k)=sum(rr(1:me)*apu(1:me,k))
                ppv(n,k)=sum(rr(1:me)*apv(1:me,k))
                ppw(n,k)=sum(rr(1:me)*apw(1:me,k))
              end do
            end if
            if (nud_sfh/=0.and.lblock) then
              pph(n)=sum(rr(1:me)*aph(1:me))
            end if
          end do
          qsum(igrd(1:il_g,j,ipass))=psum(1:il_g)
          if (nud_sst/=0) then
            do k=1,kd
              qp(igrd(1:il_g,j,ipass),k)=pp(1:il_g,k)
            end do
          end if
          if (nud_sss/=0) then
            do k=1,kd
              qps(igrd(1:il_g,j,ipass),k)=pps(1:il_g,k)
            end do
          end if
          if (nud_ouv/=0) then
            do k=1,kd
              qpu(igrd(1:il_g,j,ipass),k)=ppu(1:il_g,k)
              qpv(igrd(1:il_g,j,ipass),k)=ppv(1:il_g,k)
              qpw(igrd(1:il_g,j,ipass),k)=ppw(1:il_g,k)
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            qph(igrd(1:il_g,j,ipass))=pph(1:il_g)
          end if
        end do

        if (myid==0.and.nmaxpr==1) write(6,*) "MLO end conv"

        if (ipass==2.or.ipass==3) then
          itag=itag+1
          if (myid==0.and.nmaxpr==1) then
            write(6,*) "MLO Send arrays to local host"
          end if
          iy=il_g*il_g*(kx(ipass)-kn(ipass)+1)
          ix=iy/mproc
          a=il_g*(kx(ipass)-kn(ipass)+1)
          c=kn(ipass)
          if (myid/=hproc) then
            do j=ns,ne
              do kpass=kn(ipass),kx(ipass)
                do n=1,il_g
                  yy(n+il_g*(kpass-c)+a*(j-ns))
     &              =qsum(igrd(n,j,kpass))
                end do
              end do
            end do
          end if
          call ccmpi_gatherx(zz(1:iy),yy(1:ix),0,comm_face)
          if (myid==hproc) then
            do j=ne+1,il_g
              do kpass=kn(ipass),kx(ipass)
                do n=1,il_g
                  qsum(igrd(n,j,kpass))
     &              =zz(n+il_g*(kpass-c)+a*(j-1))
                end do
              end do
            end do
          end if
          if (nud_sfh/=0.and.lblock) then
            if (myid/=hproc) then
              do j=ns,ne
                do kpass=kn(ipass),kx(ipass)
                  do n=1,il_g
                    yy(n+il_g*(kpass-c)+a*(j-ns))
     &                =qph(igrd(n,j,kpass))
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do kpass=kn(ipass),kx(ipass)
                  do n=1,il_g
                    qph(igrd(n,j,kpass))
     &                =zz(n+il_g*(kpass-c)+a*(j-1))
                  end do
                end do
              end do
            end if
          end if
          iy=il_g*il_g*kd*(kx(ipass)-kn(ipass)+1)
          ix=iy/mproc
          a=il_g*(kx(ipass)-kn(ipass)+1)
          b=il_g*(kx(ipass)-kn(ipass)+1)*kd
          c=kn(ipass)
          if (nud_sst/=0) then
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      yy(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qp(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd             
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qp(igrd(n,j,kpass),k)
     &                  =zz(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
          end if
          if (nud_sss/=0) then
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      yy(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qps(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd             
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qps(igrd(n,j,kpass),k)
     &                  =zz(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
          end if
          if (nud_ouv/=0) then
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      yy(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qpu(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd             
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qpu(igrd(n,j,kpass),k)
     &                  =zz(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      yy(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qpv(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd             
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qpv(igrd(n,j,kpass),k)
     &                  =zz(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
            if (myid/=hproc) then
              do j=ns,ne
                do k=1,kd
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      yy(n+il_g*(kpass-c)+a*(k-1)+b*(j-ns))
     &                  =qpw(igrd(n,j,kpass),k)
                    end do
                  end do
                end do
              end do
            end if
            call ccmpi_gatherx(zz(1:iy),yy(1:ix),0,comm_face)
            if (myid==hproc) then
              do j=ne+1,il_g
                do k=1,kd             
                  do kpass=kn(ipass),kx(ipass)
                    do n=1,il_g
                      qpw(igrd(n,j,kpass),k)
     &                  =zz(n+il_g*(kpass-c)+a*(k-1)+b*(j-1))
                    end do
                  end do
                end do
              end do
            end if
          end if
        end if
          
      end do
      
      return  
      end subroutine mlospeclocal_s
      
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
      
