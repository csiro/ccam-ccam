      subroutine onthefly(nested,kdate_r,ktime_r,
     .                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
!     following not used or returned if called by nestin (i.e.nested=1)   
     .                    tgg,wb,wbice,snowd,
     .                    tggsn,smass,ssdn,ssdnn,snage,isflag,
     .                    rtsoil, ! MJT cable
     .                    urban) ! MJT urban
!     Target points use values interpolated to their 4 grid "corners";
!     these corner values are then averaged to the grid centres
!     N.B. this means will get different fields with io_in=-1 from io_in=1
!     Called by either indata or nestin
!     nested=0  for calls from indata; 1  for calls from nestin     

      use cc_mpi
      use define_dimensions, only : ncs, ncp ! MJT cable      
      use utilities
      implicit none
      integer, parameter :: ntest=0
      integer, parameter :: nord=3        ! 1 for bilinear, 3 for bicubic
!     related to cctocc4                       

!     Note: 1) The arrays are replaced in place
!           2) kl is assumed to be the same for both grids
      include 'newmpar.h'
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      include 'carbpools.h' ! MJT cable
      include 'const_phys.h'
      include 'latlong.h'  ! rlatt,rlongg,
c     include 'map.h'  ! zs,land & used for giving info after all setxyz
      include 'parm.h'
      include 'sigs.h'
      include 'soil.h'
      include 'stime.h'   ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'vegpar.h' ! MJT cable
      include 'vecsuv_g.h'
      include 'vvel.h'
      include 'xyzinfo_g.h'  ! x,y,z,wts
      include 'mpif.h'
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      real sigin
      integer ik,jk,kk
      common/sigin/ik,jk,kk,sigin(kl)  ! for vertint, infile
      ! rlong4 needs to be shared with setxyz. These are global arrays.
      real, dimension(ifull_g,4) :: rlong4, rlat4
      common/workglob/rlong4,rlat4   ! shared with setxyz
      ! Used in the global interpolation
      real, dimension(ifull_g) :: u_g, v_g, t_g, qg_g, uc, vc, wc,
     &                            ucc, vcc, wcc, uct, vct, wct
      real, dimension(ifull_g) :: ax_t, ay_t, az_t, bx_t, by_t, bz_t
      real, dimension(ifull) :: tss_l, tss_s, pmsl
      logical, dimension(ifull) :: land_t
      real, dimension(ifull_g,4) :: xg4, yg4
      integer, dimension(ifull_g,4) :: nface4
      real rotpoles(3,3),rotpole(3,3)

!     These are local arrays, not the versions in arrays.h
!     Use in call to infile, so are dimensioned ifull rather than ifull_g
      real psl(ifull),zss(ifull),tss(ifull),fracice(ifull),
     & wb(ifull,ms),wbice(ifull,ms),snowd(ifull),sicedep(ifull),
     & t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     & tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     & ssdnn(ifull),snage(ifull),rtsoil(ifull), ! MJT cable
     & urban(ifull,12) ! MJT urban
      integer isoilm_h(ifull) ! MJT lsmask
      ! Dummy variables here replace the aliasing use of aa, bb in infile call
      real, dimension(ifull) :: dum5
      integer isflag(ifull)
      ! Will get odd results unless this is on process 0 ???
!!    integer, parameter :: id1=3, jd1=60
      real :: rlong0_t, rlat0_t, schmidt_t,  rlong0x, rlat0x, schmidtx,
     &         ds_t, timegb, spval
      integer :: nemi, id_t, jd_t, idjd1, kdate_r, ktime_r,
     &           nested, i, j, k, m, iq, id2, jd2, idjd2, ii, jj,
     &           np, numneg, norder, idjd_t
      entry onthefl(nested,kdate_r,ktime_r,
     .                    psl,zss,tss,sicedep,fracice,t,u,v,qg)

      if ( myid==0 ) print *,'entering onthefly for nested = ',nested
!     save cc target file geometry
      rlong0_t=rlong0
      rlat0_t=rlat0
      schmidt_t=schmidt
      ds_t=ds
c     zs_t(:)=zs(1:ifull)
      land_t(:)=land(:)  
!        save "target" ax, bx etc  - used in transforming target u & v
         ax_t(:) = ax_g(:)
         ay_t(:) = ay_g(:)
         az_t(:) = az_g(:)
         bx_t(:) = bx_g(:)
         by_t(:) = by_g(:)
         bz_t(:) = bz_g(:)

c     start of processing loop 
      nemi=3   !  MJT lsmask
      if(ktau<3.and.myid==0)print *,'search for kdate_s,ktime_s >= ',
     &                                          kdate_s,ktime_s
      id_t=id
      jd_t=jd
      idjd_t=idjd ! Only id+il_g*(jd-1) if it's on process 0.
!!    activate id1,jd1 only if you know approx. corresponding "source" id,jd      
!!    id=id1
!!    jd=jd1
!!    idjd=id+il_g*(jd-1)
      idjd1=idjd

      ! infile reads and distributes data to appropriate processors, so
      ! all processors must call it here.
      ! illegal aliasing of arguments removed now
      if(nested==0)then
        call infile(nested,kdate_r,ktime_r,timegb,ds,
     &            psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &            tgg,wb,wbice,dum5,snowd,  ! dum5 is alb
     &            tggsn,smass,ssdn,ssdnn,snage,isflag,rtsoil, ! MJT cable
     &            isoilm_h,urban) ! MJT lsmask ! MJT urban
      else
        call infil(nested,kdate_r,ktime_r,timegb,ds,
     &            psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &            isoilm_h) ! MJT lsmask
      endif   
!     N.B. above infile call returns values for ik,jk,kk of source data
     
!     Purpose of setxyz call is to get rlat4 rlong4 (and so xx4 yy4) 
!     for the source grid. Only process 0 needs to do this here
      if ( myid==0 ) then
         rlong0=rlong0x
         rlat0=rlat0x
         schmidt=schmidtx
!        N.B. -ve ik in call setxyz preservess target rlat4, rlong4         
c        call setxyz(myid,-ik)      ! for source data geometry        ******************
         call setxyz(-ik,xx4,yy4,myid)  ! for source data geometry    ****     
!     rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!     rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!     rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
         rotpoles = calc_rotpole(rlong0,rlat0)
         if(ktau<3)then
            print *,'nfly,nord ',nfly,nord
            print *,'kdate_r,ktime_r,ktau,ds',
     &               kdate_r,ktime_r,ktau,ds
            print *,'ds,ds_t ',ds,ds_t
            if ( nproc==1 ) print *,'a zss(idjd1) ',zss(idjd1)
            print *,'rotpoles:'
            do i=1,3
               print 9,(i,j,j=1,3),(rotpoles(i,j),j=1,3)
            enddo
         endif                  ! (ktau<3)

      id=id_t
      jd=jd_t
      idjd=idjd_t
!     rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!     rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!     rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
         rotpole = calc_rotpole(rlong0_t,rlat0_t)
         if(nmaxpr==1)then   ! already in myid==0 loop
            print *,'in onthefly rotpole:'
            do i=1,3
               print 9,(i,j,j=1,3),(rotpole(i,j),j=1,3)
 9             format(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)
            enddo
            print *,'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
         endif                  ! (nmaxpr==1)

         if(nmaxpr==1)then  ! already in myid==0 loop
            print *,'before latltoij for id,jd: ',id,jd
            if ( nproc==1 ) then
               ! Diagnostics will only be correct if nproc==1
               print *,'rlong4(1-4) ',(rlong4(idjd,m),m=1,4)
               print *,'rlat4(1-4) ',(rlat4(idjd,m),m=1,4)
            end if
            print *,'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,schmidtx 
         endif                  ! (nmaxpr==1)
         do m=1,4
            do iq=1,ifull_g
               call latltoij(rlong4(iq,m),rlat4(iq,m),         !input
     &                       xg4(iq,m),yg4(iq,m),nface4(iq,m), !output (source)
     &                       xx4,yy4,ik)
            enddo
         enddo
         if(nproc==1.and.nmaxpr==1)then
           ! Diagnostics will only be correct if nproc==1
           id2=nint(xg4(idjd,1))
           jd2=il*nface4(idjd,1)+nint(yg4(idjd,1))
           idjd2=id2+il*(jd2-1)
            print *,'after latltoij giving id2,jd2,idjd2: ',
     .                                     id2,jd2,idjd2
            print *,'nface4(1-4) ',(nface4(idjd,m),m=1,4)
            print *,'xg4(1-4) ',(xg4(idjd,m),m=1,4)
            print *,'yg4(1-4) ',(yg4(idjd,m),m=1,4)
            if(nested==0)then
              write(6,"('wb_s(1)#  ',9f7.3)") 
     .            ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
              write(6,"('wb_s(ms)# ',9f7.3)") 
     .            ((wb(ii+(jj-1)*il,ms),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
            endif  ! (nested==0)
         endif

      end if ! myid==0

      if(nfly==2)then         ! needs pmsl in this case (preferred)
         call mslp(pmsl,psl,zss,t(1:ifull,:))  
      endif

      ! All the following processing is done on processor 0
      ! Avoid memory blow out by only having single level global arrays
      do k=1,kk
         if ( myid==0 ) then
            call ccmpi_gather(u(:,k), u_g)
            call ccmpi_gather(v(:,k), v_g)
            call ccmpi_gather(t(:,k), t_g)
            call ccmpi_gather(qg(:,k), qg_g)
            do iq=1,ik*jk
!              first set up winds in Cartesian "source" coords
               uc(iq)=ax_g(iq)*u_g(iq) + bx_g(iq)*v_g(iq)
               vc(iq)=ay_g(iq)*u_g(iq) + by_g(iq)*v_g(iq)
               wc(iq)=az_g(iq)*u_g(iq) + bz_g(iq)*v_g(iq)
!              now convert to winds in "absolute" Cartesian components
               ucc(iq)=uc(iq)*rotpoles(1,1)+vc(iq)*rotpoles(1,2)
     &                                +wc(iq)*rotpoles(1,3)
               vcc(iq)=uc(iq)*rotpoles(2,1)+vc(iq)*rotpoles(2,2)
     &                                +wc(iq)*rotpoles(2,3)
               wcc(iq)=uc(iq)*rotpoles(3,1)+vc(iq)*rotpoles(3,2)
     &                                +wc(iq)*rotpoles(3,3)
            enddo               ! iq loop
            if(ktau<3.and.k==1.and.nproc==1)then
               print *,'uc,vc,wc: ',uc(id),vc(idjd1),wc(idjd1)
               print *,'ucc,vcc,wcc: ',ucc(idjd1),vcc(idjd1),wcc(idjd1)
               print *,'calling ints4 for k= ',k
            endif

!      interpolate all required arrays to new C-C positions
!      don't need to do map factors and Coriolis on target grid
            np=0                ! controls prints in ints4
            call ints4(t_g, nface4,xg4,yg4,nord,ik)  ! ints4 on source grid
            call ints4(qg_g,nface4,xg4,yg4,nord,ik)
            call ints4(ucc, nface4,xg4,yg4,nord,ik)
            call ints4(vcc, nface4,xg4,yg4,nord,ik)
            call ints4(wcc, nface4,xg4,yg4,nord,ik)

c      ********************** N.B. tracers not ready yet
c      if(iltin>1)then
c        do ntr=1,ntracin
c         call ints4(tr(1,k,ntr),nface4,xg4,yg4,nord,ik)
c        enddo
c      endif
 
            do iq=1,ifull_g
!       now convert to "target" Cartesian components (transpose used)
               uct(iq)=ucc(iq)*rotpole(1,1)+vcc(iq)*rotpole(2,1)
     &                           +wcc(iq)*rotpole(3,1)
               vct(iq)=ucc(iq)*rotpole(1,2)+vcc(iq)*rotpole(2,2)
     &                           +wcc(iq)*rotpole(3,2)
               wct(iq)=ucc(iq)*rotpole(1,3)+vcc(iq)*rotpole(2,3)
     &                           +wcc(iq)*rotpole(3,3)
!       then finally to "target" local x-y components
               u_g(iq) = ax_t(iq)*uct(iq) + ay_t(iq)*vct(iq) +
     &                   az_t(iq)*wct(iq)
               v_g(iq) = bx_t(iq)*uct(iq) + by_t(iq)*vct(iq) +
     &                   bz_t(iq)*wct(iq)
            enddo               ! iq loop
            if(ktau<3.and.k==1.and.nproc==1)then
               ! This only works if idjd is on processor 0
               print *,'interp. ucc,vcc,wcc: ',ucc(idjd),vcc(idjd),
     &                  wcc(idjd)
               print *,'uct,vct,wct: ',uct(idjd),vct(idjd),wct(idjd)
               print *,'ax,ay,az ',ax_t(idjd),ay_t(idjd),az_t(idjd)
               print *,'bx,by,bz ',bx_t(idjd),by_t(idjd),bz_t(idjd)
               print *,'final u , v: ',u_g(idjd),v_g(idjd)
            endif
            call ccmpi_distribute(u(:,k), u_g)
            call ccmpi_distribute(v(:,k), v_g)
            call ccmpi_distribute(t(:,k), t_g)
            call ccmpi_distribute(qg(:,k), qg_g)
         else ! myid /= 0
            call ccmpi_gather(u(:,k))
            call ccmpi_gather(v(:,k))
            call ccmpi_gather(t(:,k))
            call ccmpi_gather(qg(:,k))
            call ccmpi_distribute(u(:,k))
            call ccmpi_distribute(v(:,k))
            call ccmpi_distribute(t(:,k))
            call ccmpi_distribute(qg(:,k))
         end if ! myid==0
      enddo  ! k loop

!     below we interpolate quantities which may be affected by land-sea mask

!     set up land-sea mask from either tss or zss
      !-------------------------------------------
      ! MJT lsmask
      if(nemi==3)then 
        land(:)=isoilm_h(:).gt.0
        if (any(isoilm_h(:).lt.0)) nemi=2
      end if
      !-------------------------------------------
      if(nemi==2)then
         numneg=0
         do iq=1,ifull
            if(tss(iq)>0)then ! over land
               land(iq)=.true.
            else                ! over sea
               land(iq)=.false.
               numneg=numneg+1
            endif               ! (tss(iq)>0) .. else ..
         enddo
         if(numneg==0)nemi=1  ! should be using zss in that case
      endif                     !  (nemi==2)
      if ( myid==0 )print *,'using nemi = ',nemi
      if(nemi==1)then
         land(:) = zss(:) > 0.
      endif                     !  (nemi==1)

      spval=999.
      do iq=1,ifull
         if(land(iq))then       ! over land
            tss_l(iq)=tss(iq)
            tss_s(iq)=spval
            sicedep(iq)=spval
            fracice(iq)=spval
         else                   ! over sea
            numneg=numneg+1
            tss_s(iq)=abs(tss(iq))
            tss_l(iq)=spval
         endif  !   (land(iq)) .. else ..
      enddo     ! iq loop
      
      if(nproc==1.and.nmaxpr==1)then
        print *,'before fill tss ',tss(idjd2)
        print *,'before fill tss_l, tss_s ',tss_l(idjd2),tss_s(idjd2)
        print *,'before fill/ints4 sicedep ',sicedep(idjd2)
c        print *,'before fill wb'
c        write(6,"('wb_s(1)#  ',9f7.3)") 
c     .          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
      endif  ! (nproc==1.and.nmaxpr==1)
      call fill_cc(tss_l,spval)
      call fill_cc(tss_s,spval)
      call fill_cc(sicedep,spval)
      call fill_cc(fracice,spval)
      if(nproc==1.and.nmaxpr==1)then
        print *,'after fill tss_l, tss_s ',tss_l(idjd2),tss_s(idjd2)
        print *,'after fill sicedep ',sicedep(idjd2)
c        print *,'after fill wb'
c        write(6,"('wb_s(1)#  ',9f7.3)") 
c     .          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
        print *,'before ints4 psl(idjd2),zss(idjd2) ',
     .                        psl(idjd),zss(idjd2)
      endif  ! (nproc==1.and.nmaxpr==1)

      if(nfly>0)then
        norder=1
      else
        norder=nord
      endif

      ! The routine doints4 does the gather, calls ints4 and redistributes
      call doints4(psl ,     nface4,xg4,yg4,norder,ik)
      call doints4(zss ,nface4,xg4,yg4,norder,ik)  
      if(nfly==2)then
        call doints4(pmsl,   nface4,xg4,yg4,nord,ik)
!       invert pmsl to get psl
        call to_psl(pmsl,psl,zss,t)  
      endif  ! (nfly==2)
      call doints4(tss_l ,   nface4,xg4,yg4,nord,ik)
      call doints4(tss_s ,   nface4,xg4,yg4,nord,ik)
      if(nproc==1.and.nmaxpr==1)then
         print *,'after ints4 idjd,zss(idjd) ',idjd,zss(idjd)
         print *,'after ints4 psl,pmsl ',psl(idjd),pmsl(idjd)
c        print *,'after ints4 wb_t'
c        write(6,"('wb_t(1)#  ',9f7.3)") 
c     .           ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
      endif  ! (nproc==1.and.nmaxpr==1)
      call doints4(sicedep,nface4,xg4,yg4,nord,ik)
      call doints4(fracice,nface4,xg4,yg4,nord,ik)
      if ( nproc==1 ) print *,'after ints4 sicedep ',sicedep(idjd)
      if(nested==0)then
        do iq=1,ifull
         if(.not.land(iq))then       
            snowd(iq)=spval
            do k=1,ms
               tgg(iq,k)=spval
               wb(iq,k)=spval
            enddo
         endif  !   (.not.land(iq)) 
        enddo   ! iq loop
        call fill_cc(snowd,spval)
        do k=1,ms
         call fill_cc(tgg(1,k),spval)
         call fill_cc(wb(1,k),spval)
        enddo
        call doints4(snowd,  nface4,xg4,yg4,nord,ik)
        do k=1,ms
         call doints4(tgg(1,k),nface4,xg4,yg4,nord,ik)
         call doints4(wb(1,k) ,nface4,xg4,yg4,nord,ik)
        enddo
        !--------------------------------------------------
        ! MJT urban
        if (nurban.ne.0) then
          do k=1,12
            where ((.not.land(:)).or.(urban(:,k).ge.399.))
              urban(:,k)=spval
            end where
            call fill_cc(urban(:,k),spval)
            call doints4(urban(:,k),nface4,xg4,yg4,nord,ik)
          end do
        end if
        !--------------------------------------------------
        !--------------------------------------------------
        ! MJT cable
        if ((nsib.eq.4).or.(nsib.eq.6)) then
          call doints4(rtsoil(:),nface4,xg4,yg4,nord,ik)
          do k=1,ncp
            call doints4(cplant(:,k),nface4,xg4,yg4,nord,ik)
          end do
          do k=1,ncs
            call doints4(csoil(:,k),nface4,xg4,yg4,nord,ik)
          end do
          call doints4(cansto(:),nface4,xg4,yg4,nord,ik)
        end if
        !--------------------------------------------------	           	        
c       incorporate target land mask effects for initial fields
        do iq=1,ifull
         if(land_t(iq))then
           tss(iq)=tss_l(iq)
         else
           tgg(iq,1)=tss_s(iq)   ! no sign switch in CCAM
         endif
        enddo  ! iq loop
!       onthefly not yet handling wbice, tggsn, smass etc  
!       so set infile-style defaults
        smass(:,:)=0.
        tggsn(:,:)=280.     ! just a default
        isflag(:) = 0
        snage(:)  = 0.
        do iq=1,ifull
         if(snowd(iq)>100.)then
           ssdn(iq,1)=240.
         else
           ssdn(iq,1) = 140.
         endif  ! (snowd(iq)>100.)
         do k=2,3
          ssdn(iq,k)=ssdn(iq,1)
         enddo
         ssdnn(iq)  = ssdn(iq,1)
         if(snowd(iq)>0.)tgg(iq,1)=min(tgg(iq,1),270.1)
        enddo   ! iq loop
c        do k=1,ms     ! wbice presets done near end of indata from Dec 07
c!         following linearly from 0 to .99 for tgg=tfrz down to tfrz-5
c          wbice(:,k)=
c     &         min(.99,max(0.,.99*(273.1-tgg(:,k))/5.))*wb(:,k) ! jlm
c         enddo ! ms
      endif    ! (nested==0)

c     incorporate other target land mask effects
      do iq=1,ifull
         if(land_t(iq))then
           tss(iq)=tss_l(iq)
           sicedep(iq)=0.
           fracice(iq)=0.
         else
           tss(iq)=tss_s(iq)   ! no sign switch in CCAM
           if(sicedep(iq)<.05)then ! for sflux
             sicedep(iq)=0.
             fracice(iq)=0.
           endif
         endif
      enddo  ! iq loop
      if(nproc==1.and.nmaxpr==1)then
         print *,'after ints tss_l, tss_s ',tss_l(idjd),tss_s(idjd)
         print *,'after ints tss',tss(idjd)
      endif  ! (nproc==1.and.nmaxpr==1)

!     end of processing loop

      rlong0x=rlong0_t  ! for indata cross-check
      rlat0x=rlat0_t
      schmidtx=schmidt_t

!     restore target zs and land arrays
c     zs(1:ifull) = zs_t(:)
      land(:) = land_t(:)  

      rlong0=rlong0_t
      rlat0=rlat0_t
      schmidt=schmidt_t
      ds=ds_t
      id=id_t
      jd=jd_t
      idjd=idjd_t
      if ( myid==0 ) then  
!       call setxyz(myid,il_g)  ! for target  ** produces rlong4,rlat4 **
        call setxyz(il_g,xx4,yy4,myid) ! for target ** produces rlong4,rlat4 **
        print *,'after target setxyz'
      end if

      end subroutine onthefly

      subroutine doints4(s,nface4 ,xg4 ,yg4,nord,ik)  ! does calls to intsb
      use cc_mpi
      implicit none
      include 'newmpar.h'
      real, dimension(ifull), intent(inout) :: s
      integer, intent(in), dimension(ifull_g,4) :: nface4
      real, intent(in), dimension(ifull_g,4) :: xg4, yg4
      integer, intent(in) :: ik, nord
      real, dimension(ifull_g) :: s_g

      if ( myid ==0 ) then
         call ccmpi_gather(s,s_g)
         call ints4(s_g,nface4 ,xg4 ,yg4,nord,ik) 
         call ccmpi_distribute(s,s_g)
      else
         call ccmpi_gather(s)
         call ccmpi_distribute(s)
      endif
      end subroutine doints4

      subroutine ints4(s,nface4 ,xg4 ,yg4,nord,ik)  ! does calls to intsb
      use cc_mpi, only : mydiag
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      integer, parameter :: ntest=0
      real, dimension(ifull_g), intent(inout) :: s
      integer, intent(in), dimension(ifull_g,4) :: nface4
      real, intent(in), dimension(ifull_g,4) :: xg4, yg4
      integer, intent(in) :: ik, nord
      real wrk(ifull_g,4)
      integer :: iq, m
      integer :: idx = 25, jdx = 218, idjdx=10441

      if(nord==1)then
         do m=1,4
            call ints_blb(s,wrk(1,m),nface4(1,m),xg4(1,m),yg4(1,m),ik)
         enddo
      else
         do m=1,4
            call intsb(s,wrk(1,m),nface4(1,m),xg4(1,m),yg4(1,m),ik)
         enddo
      endif   ! (nord==1)  .. else ..
      if(ntest>0.and.mydiag)then
         print *,'in ints4 for id,jd,nord: ',id,jd,nord
         print *,'nface4(1-4) ',(nface4(idjd,m),m=1,4)
         print *,'xg4(1-4) ',(xg4(idjd,m),m=1,4)
         print *,'yg4(1-4) ',(yg4(idjd,m),m=1,4)
         print *,'wrk(1-4) ',(wrk(idjd,m),m=1,4)
      endif
!     average 4 m values to central value
      do iq=1,ifull_g
         s(iq)=.25*(wrk(iq,1)+wrk(iq,2)+wrk(iq,3)+wrk(iq,4))
      enddo

      end subroutine ints4

      subroutine intsb(s,sout,nface,xg,yg,ik)   ! N.B. sout here
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
c     later may wish to save idel etc between array calls
c     this one does linear interp in x on outer y sides
c     doing x-interpolation before y-interpolation
!     This is a global routine 
      include 'newmpar.h'
      include 'parm.h'
      real, dimension(ik*ik), intent(in) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ifull_g) :: nface
      real, intent(in), dimension(ifull_g) :: xg, yg
      real sx(-1:ik+2,-1:ik+2,0:npanels)
      real r(4)

      include 'indices_g.h' ! in,is,iw,ie,inn,iss,iww,iee
      integer :: ind, i, j, n, iq
      ind(i,j,n)=i+(j-1)*ik+n*ik*ik  ! *** for n=0,npanels
c     this is intsb           EW interps done first
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
      do n=0,npanels
       do j=1,ik
        do i=1,ik
         sx(i,j,n)=s(ind(i,j,n))
        enddo  ! i loop
        sx(0,j,n)=s(iw_g(ind(1,j,n)))
        sx(-1,j,n)=s(iww_g(ind(1,j,n)))
        sx(ik+1,j,n)=s(ie_g(ind(ik,j,n)))
        sx(ik+2,j,n)=s(iee_g(ind(ik,j,n)))
       enddo   ! j loop
       do i=1,ik
        sx(i,0,n)=s(is_g(ind(i,1,n)))
        sx(i,-1,n)=s(iss_g(ind(i,1,n)))
        sx(i,ik+1,n)=s(in_g(ind(i,ik,n)))
        sx(i,ik+2,n)=s(inn_g(ind(i,ik,n)))
       enddo  ! i loop
c      for ew interpolation, sometimes need (different from ns):
c          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
       sx(-1,0,n)=s(lwws_g(n))
       sx(0,0,n)=s(lws_g(n))
       sx(0,-1,n)=s(lwss_g(n))
       sx(ik+1,0,n)=s(les_g(n))
       sx(ik+2,0,n)=s(lees_g(n))
       sx(ik+1,-1,n)=s(less_g(n))
       sx(-1,ik+1,n)=s(lwwn_g(n))
       sx(0,ik+2,n)=s(lwnn_g(n))
       sx(ik+2,ik+1,n)=s(leen_g(n))
       sx(ik+1,ik+2,n)=s(lenn_g(n))
       sx(0,ik+1,n)   =s(iwn_g(ind(1,ik,n)))
       sx(ik+1,ik+1,n)=s(ien_g(ind(ik,ik,n)))
      enddo    ! n loop

      do iq=1,ifull_g   ! runs through list of target points
c       if(iq==idjd)print *,'iq,nface,xg,yg ',
c    .                         iq,nface(iq),xg(iq),yg(iq)
        n=nface(iq)
        idel=int(xg(iq))
        xxg=xg(iq)-idel
c       yg here goes from .5 to il +.5
        jdel=int(yg(iq))
        yyg=yg(iq)-jdel
c       if(iq==idjd)then
c         print *,'iq,idel,xxg,jdel,yyg,n ',
c    .             iq,idel,xxg,jdel,yyg,n
c         print *,'sx nn=1',sx(idel  ,jdel-1,n),sx(idel+1,jdel-1,n)
c         print *,'sx nn=2',sx(idel-1,jdel  ,n),sx(idel  ,jdel  ,n)
c    .                     ,sx(idel+1,jdel  ,n),sx(idel+2,jdel  ,n)
c         print *,'sx nn=3',sx(idel-1,jdel+1,n),sx(idel  ,jdel+1,n)
c    .                     ,sx(idel+1,jdel+1,n),sx(idel+2,jdel+1,n)
c         print *,'sx nn=4',sx(idel  ,jdel+2,n),sx(idel+1,jdel+2,n)
c       endif
        do nn=2,3       ! N.B.
         c1=sx(idel-1,jdel+nn-2,n)
         c2=sx(idel  ,jdel+nn-2,n)
         c3=sx(idel+1,jdel+nn-2,n)
         c4=sx(idel+2,jdel+nn-2,n)
         r(nn)=((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     .     -xxg*(1.+xxg)*c4/3.)
     .     +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        enddo    ! nn loop
c       r       ={(1-x     )*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c         -x     *(1+x     )*c4/3}
c         +x    *(1+x     )*(2-x     )*c3}/2
        do nn=1,4,3       ! N.B.
         c2=sx(idel  ,jdel+nn-2,n)
         c3=sx(idel+1,jdel+nn-2,n)
         r(nn)=(1.-xxg)*c2 +xxg*c3
        enddo    ! nn loop
c       array(iq)=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
c    .             -yyg*(1.+yyg)*r(4)/3.)
c    .             +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
c      following does Bermejo Staniforth
        aaa=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
     .      -yyg*(1.+yyg)*r(4)/3.)
     .      +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
c       if(iq==idjd)print *,'before BS aaa sx sx+ sx0+ sx++',aaa,
c    .                         sx(idel,jdel,n),sx(idel+1,jdel,n),
c    .                         sx(idel,jdel+1,n),sx(idel+1,jdel+1,n)
        aaa=min( aaa , max( sx(idel,jdel,n),sx(idel+1,jdel,n),
     .                      sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
        sout(iq)=max( aaa , min( sx(idel,jdel,n),sx(idel+1,jdel,n),
     .                        sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
      enddo    ! iq loop

      end subroutine intsb

      subroutine ints_blb(s,sout,nface,xg,yg,ik) 
c     this one does bi-linear interpolation only
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      real, dimension(ik*ik), intent(in) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ifull_g) :: nface
      real, intent(in), dimension(ifull_g) :: xg, yg
      real sx(-1:ik+2,-1:ik+2,0:npanels)
      include 'indices_g.h' ! in,is,iw,ie,inn,iss,iww,iee
      integer :: ind, i, j, n, iq, idel, jdel, ik
      real :: xxg, yyg
      ind(i,j,n)=i+(j-1)*ik+n*ik*ik  ! *** for n=0,npanels
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
c                    but for bi-linear only need 0:il+1 &  0:il+1
      do n=0,npanels
       do j=1,ik
        do i=1,ik
         sx(i,j,n)=s(ind(i,j,n))
        enddo  ! i loop
       enddo   ! j loop
       do j=1,ik
        sx(0,j,n)=s(iw_g(ind(1,j,n)))
        sx(ik+1,j,n)=s(ie_g(ind(ik,j,n)))
       enddo   ! j loop
       do i=1,ik
        sx(i,0,n)=s(is_g(ind(i,1,n)))
        sx(i,ik+1,n)=s(in_g(ind(i,ik,n)))
       enddo  ! i loop

       sx(0,0,n)=s(lws_g(n))
       sx(ik+1,0,n)=s(les_g(n))
       sx(0,ik+1,n)   =s(iwn_g(ind(1,ik,n)))
       sx(ik+1,ik+1,n)=s(ien_g(ind(ik,ik,n)))
      enddo    ! n loop

      do iq=1,ifull_g  ! runs through list of target points
       n=nface(iq)
       idel=int(xg(iq))
       xxg=xg(iq)-idel
       jdel=int(yg(iq))
       yyg=yg(iq)-jdel
c      print *,'iq,idel,jdel,n ',iq,idel,jdel,n
       sout(iq)=yyg*(xxg*sx(idel+1,jdel+1,n)
     .               +(1.-xxg)*sx(idel,jdel+1,n))
     .    +(1.-yyg)*(xxg*sx(idel+1,jdel,n)
     .               +(1.-xxg)*sx(idel,jdel,n))
      enddo    ! iq loop

      end subroutine ints_blb

      subroutine fill_cc(a_io,value)
c     routine fills in interior of an array which has undefined points
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'indices.h'
      include 'mpif.h'
      real a_io(ifull)         ! input and output array
      real value            ! array value denoting undefined
      real b(ifull), a(ifull+iextra)
      integer :: num, nrem, iq, neighb, nrem_g, nrem_gmin, ierr
      real :: av, avx
      
      a(1:ifull) = a_io(:)
      num=0
      nrem_g = 1    ! Just for first iteration
      nrem_gmin = 1 ! Just for first iteration
!     nrem_gmin used to avoid infinite loops, e.g. for no sice
      do while ( nrem_g > 0 .and. nrem_gmin<ifull )
         ! This has to loop until all are finished otherwise the bounds call
         ! doesn't work.
         call bounds(a)
         nrem=0
         num=num+1
         do iq=1,ifull
            b(iq)=a(iq)
            if(a(iq)==value)then
               neighb=0
               av=0.
c              if(a(in(iq)).ne.value.and.a(in(iq)).ne.2.)then
               if(a(in(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(in(iq))
               endif
c              if(a(ie(iq)).ne.value.and.a(ie(iq)).ne.2.)then
               if(a(ie(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(ie(iq))
               endif
c              if(a(iw(iq)).ne.value.and.a(iw(iq)).ne.2.)then
               if(a(iw(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(iw(iq))
               endif
c              if(a(is(iq)).ne.value.and.a(is(iq)).ne.2.)then
               if(a(is(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(is(iq))
               endif
               if(neighb>0)then
                  b(iq)=av/neighb
                  avx=av
               else
                  nrem=nrem+1   ! current number of points without a neighbour
               endif
            endif
         end do
         do iq=1,ifull
            a(iq)=b(iq)
         enddo
c         call MPI_AllReduce(nrem, nrem_g, 1, MPI_INTEGER, MPI_SUM, 
c    &                      MPI_COMM_WORLD, ierr )
         call MPI_AllReduce(nrem, nrem_g, 1, MPI_INTEGER, MPI_MAX, 
     &                      MPI_COMM_WORLD, ierr )
         call MPI_AllReduce(nrem, nrem_gmin, 1, MPI_INTEGER, MPI_MIN, 
     &                      MPI_COMM_WORLD, ierr )
         if(nrem_g>0.and.myid==0)then
c           if(num<=2) then
c              print *,'after 1/2 time thru fill loop num,nrem,avx = ',
c    &                                           num,nrem_g,avx
c           end if
         endif                  ! (nrem>0)
      end do
      a_io(1:ifull) = a(1:ifull)
      return
      end
