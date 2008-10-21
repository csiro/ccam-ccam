      subroutine onthefly(nested,kdate_r,ktime_r,
     .                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
!     following not used or returned if called by nestin (i.e.nested=1)   
     .                    tgg,wb,wbice,snowd,qfg,qlg,   ! 0808
     .                    tggsn,smass,ssdn,ssdnn,snage,isflag,
     .                    rtsoil,urban) ! MJT cable ! MJT urban
!     Target points use values interpolated to their 4 grid "corners";
!     these corner values are then averaged to the grid centres
!     N.B. this means will get different fields with io_in=-1 from io_in=1
!     Called by either indata or nestin
!     nested=0  for calls from indata; 1  for calls from nestin     
!     ****  qfg and qlg not yet interpolated in ontheflyx
      use cc_mpi
      use utilities
      implicit none
      integer, parameter :: ntest=0
      integer, parameter :: nord=3        ! 1 for bilinear, 3 for bicubic
!     related to cctocc4                       

!     Note: 1) The arrays are replaced in place
!           2) kl is assumed to be the same for both grids
      include 'newmpar.h'
      include 'const_phys.h'
c     include 'map.h'  ! zs,land & used for giving info after all setxyz
      include 'parm.h'
      include 'sigs.h'
      include 'soil.h'
      include 'stime.h'   ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'vecsuv_g.h'
      include 'vvel.h'
      include 'mpif.h'
      real sigin
      integer ik,jk,kk
      common/sigin/ik,jk,kk,sigin(40)  ! for vertint, infile
      logical, dimension(ifull) :: land_t

      include 'darcdf.h'    ! idnc, ncid
      include 'netcdf.inc'
      integer, parameter :: nihead=54,nrhead=14
      integer nahead(nihead),ier,ilen,itype
      integer, save :: ncidold=-1,iarchi

!     These are local arrays, not the versions in arrays.h
!     Use in call to infile, so are dimensioned ifull rather than ifull_g
      real psl(ifull),zss(ifull),tss(ifull),fracice(ifull),
     & wb(ifull,ms),wbice(ifull,ms),snowd(ifull),sicedep(ifull),
     & t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     & tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     & ssdnn(ifull),snage(ifull),qfg(ifull,kl),qlg(ifull,kl),
     & rtsoil(ifull),urban(ifull,12) ! MJT cable ! MJT urban
      ! Dummy variables here replace the aliasing use of aa, bb in infile call
      integer isflag(ifull)
      ! Will get odd results unless this is on process 0 ???
!!    integer, parameter :: id1=3, jd1=60
      integer ::  kdate_r, ktime_r, nested

      if ( myid==0 )then
        print *,'entering onthefly for nested = ',nested
        if(ncid.ne.ncidold)then
          iarchi=1
          ncidold=ncid
          call ncainq(ncid,ncglobal,'int_header',itype,ilen,ier)
          call ncagt(ncid,ncglobal,'int_header',nahead,ier)
          ik=nahead(1)
          jk=nahead(2)
          kk=nahead(3)
          print *,'in onthefly ktau,ncid,iarchi,ik,jk,kk ',
     &                         ktau,ncid,iarchi,ik,jk,kk
        endif
      endif  ! ( myid==0 )
      call MPI_Bcast(ik,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_Bcast(jk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_Bcast(kk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!     save cc land file 
      land_t(:)=land(:)  

c     start of processing loop 
      if(ktau<3.and.myid==0)print *,'search for kdate_s,ktime_s >= ',
     &                                          kdate_s,ktime_s
!!    activate id1,jd1 only if you know approx. corresponding "source" id,jd      
!!    id=id1
!!    jd=jd1
!!    idjd=id+il_g*(jd-1)
      call ontheflyx(nested,land_t,kdate_r,ktime_r,
     .                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
!     following not used or returned if called by nestin (i.e.nested=1)   
     .                    tgg,wb,wbice,snowd,qfg,qlg,  ! 0808
     .                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     .                    rtsoil,urban) ! MJT cable ! MJT urban
      return
      end
      subroutine ontheflyx(nested,land_t,kdate_r,ktime_r,
     .                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
!     following not used or returned if called by nestin (i.e.nested=1)   
     .                    tgg,wb,wbice,snowd,qfg,qlg,
     .                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     .                    rtsoil,urban) ! MJT cable ! MJT urban
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
c     include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      real*8 xx4(1+4*ik,1+4*ik),yy4(1+4*ik,1+4*ik)
      real*8, dimension(ik*ik*6):: z_a,x_a,y_a
c**   xx4 & yy4 only used in indata & here, so no need to redefine after
c**   onthefly; sometime can get rid of common/bigxy4
      include 'carbpools.h' ! MJT cable      
      include 'const_phys.h'
      include 'latlong_g.h'  ! rlatt_g,rlongg_g,
c     include 'map.h'  ! zs,land & used for giving info after all setxyz
      include 'parm.h'
      include 'parmgeom.h'  ! rlong0,rlat0,schmidt  
      include 'sigs.h'
      include 'soil.h'
      include 'stime.h'   ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'vegpar.h' ! MJT cable
      include 'vecsuv_g.h'
      include 'vvel.h'
      include 'mpif.h'

!     These are local arrays, not the versions in arrays.h
!     Use in call to infile, so are dimensioned ifull rather than ifull_g
      real, dimension(ifull) :: psl,zss,tss,fracice,
     &                          snowd,sicedep,ssdnn,snage
      real, dimension(ifull,ms) :: wb,wbice,tgg
      real, dimension(ifull,3) :: tggsn,smass,ssdn
      real, dimension(ifull,kl) :: t,u,v,qg,qfg,qlg
      real, dimension(ifull) :: rtsoil ! MJT cable
      real, dimension(ifull,12) :: urban ! MJT urban
      integer, dimension(ifull) :: isflag
      real, dimension(ik*ik*6) :: psl_a,zss_a,tss_a,fracice_a,dum5,
     &      snowd_a,sicedep_a,ssdnn_a,snage_a,pmsl_a,  tss_l_a,tss_s_a
      real, dimension(ik*ik*6,ms) :: wb_a,wbice_a,tgg_a
      real, dimension(ik*ik*6,3) :: tggsn_a,smass_a,ssdn_a
      real, dimension(ik*ik*6,kk) :: t_a,u_a,v_a,qg_a,qfg_a,qlg_a
      real, dimension(ik*ik*6) :: rtsoil_a,cansto_a ! MJT cable
      real, dimension(ik*ik*6,ncp) :: cplant_a ! MJT cable
      real, dimension(ik*ik*6,ncs) :: csoil_a  ! MJT cable
      real, dimension(ik*ik*6,12) :: urban_a ! MJT urban
      integer, dimension(ik*ik*6) :: isflag_a
      real ::  rlong0x, rlat0x, schmidtx, spval
      integer ::  kdate_r, ktime_r, nemi, id2,jd2,idjd2,
     &            nested, i, j, k, m, iq, ii, jj, np, numneg

      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      integer ik, kk
      ! rlong4 needs to be shared with setxyz. These are global arrays.
      real, dimension(ifull_g,4) :: rlong4, rlat4
      common/workglob/rlong4,rlat4   ! shared with setxyz
      ! Used in the global interpolation
      real, dimension(ik*ik*6) :: ucc, vcc, wcc, uc, vc, wc
      real, dimension(ifull_g) :: u_g, v_g, t_g, qg_g,
     &                            uct_g, vct_g, wct_g
      real  uct_gg, vct_gg, wct_gg
      real, dimension(ifull) :: tss_l, tss_s, pmsl
      real, dimension(ifull_g,4) :: xg4, yg4
      integer, dimension(ifull_g,4) :: nface4
      real rotpoles(3,3),rotpole(3,3)
      real, dimension(ik*ik*6):: axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a
      real, dimension(ik*ik*6):: wts_a  ! not used here or defined in call setxyz
      real ::   timegb
      logical, dimension(ik*ik*6) :: land_a
      logical, dimension(ifull) :: land_t
      integer, dimension(ik*ik*6) :: isoilm_a ! MJT lsmask
   
      nemi=3   !  MJT lsmask
      if(m_fly==1)then
        rlong4(:,1)=rlongg_g(:)*180./pi
        rlat4(:,1)=rlatt_g(:)*180./pi
      endif

      ! infile reads and distributes data (not from onthefly) to appropriate processors, so
      ! all processors must call it here.
      ! illegal aliasing of arguments removed now
   !   if(nested==0)then ! MJT bug
        call infile(nested,kdate_r,ktime_r,timegb,ds,
     &       psl_a,zss_a,tss_a,sicedep_a,fracice_a,t_a,u_a,v_a,qg_a,
     &       tgg_a,wb_a,wbice_a,dum5,snowd_a,qfg_a,qlg_a,                ! dum5 is alb
     &       tggsn_a,smass_a,ssdn_a,ssdnn_a,snage_a,isflag_a,ik*ik*6,
     &       kk,rtsoil_a,isoilm_a,urban_a,cplant_a,csoil_a,cansto_a) ! MJT cable ! MJT lsmask ! MJT urban
   !   else
   !     call infile(nested,kdate_r,ktime_r,timegb,ds,
   !  &       psl_a,zss_a,tss_a,sicedep_a,fracice_a,t_a,u_a,v_a,qg_a,
   !  &       tgg_a,wb_a,wbice_a,dum5,snowd_a,qfg_a,qlg_a,                ! dum5 is alb
   !  &       tggsn_a,smass_a,ssdn_a,ssdnn_a,snage_a,isflag_a,ik*ik*6,
   !  &       kk,rtsoil_a,isoilm_a,urban_a,cplant_a,csoil_a,cansto_a) ! MJT cable ! MJT lsmask ! MJT urban
   !   endif   
!     N.B. above infile call returns values for ik,jk,kk of source data
!     Purpose of setxyz call is to get geometry (and so xx4 yy4) 
!     for the source grid. Only process 0 needs to do this here

ccc      if ( myid==0 ) then
!        N.B. -ve ik in call setxyz preserves TARGET rlat4, rlong4     
!        following setxyz call is for source data geom    ****   
         do iq=1,ik*ik*6
          axs_a(iq)=iq
          ays_a(iq)=iq
          azs_a(iq)=iq
         enddo      
         call setxyz(ik,rlong0x,rlat0x,-schmidtx,
     &    x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4,
     &    myid)
         rotpoles = calc_rotpole(rlong0x,rlat0x)
         if(ktau<3)then
            print *,'m_fly,nord ',m_fly,nord
            print *,'kdate_r,ktime_r,ktau,ds',
     &               kdate_r,ktime_r,ktau,ds
            if ( nproc==1 ) print *,'a zss(idjd) ',zss(idjd)
            print *,'rotpoles:'
            do i=1,3
               print 9,(i,j,j=1,3),(rotpoles(i,j),j=1,3)
            enddo
         endif                  ! (ktau<3)

!        rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!        rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!        rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
         rotpole = calc_rotpole(rlong0,rlat0)
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
         do m=1,m_fly  !  was 4, now may be set to 1 in namelist
            do iq=1,ifull_g
               call latltoij(rlong4(iq,m),rlat4(iq,m),         !input
     &                       rlong0x,rlat0x,schmidtx,          !input
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
     .                                    id2,jd2,idjd2
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
      if ( myid==0 ) then
        call mslpx(pmsl_a,psl_a,zss_a,t_a,ik*ik*6)  ! needs pmsl (preferred) 
      end if ! myid==0

      ! All the following processing is done on processor 0
      ! Avoid memory blow out by only having single level global arrays
      ! - irrelevant for onthefly after 0808
      do k=1,kk
         if ( myid==0 ) then
c            call ccmpi_gather(u_a(:,k), u_g)
c            call ccmpi_gather(v_a(:,k), v_g)
c            call ccmpi_gather(t_a(:,k), t_g)
c            call ccmpi_gather(qg_a(:,k), qg_g)
            do iq=1,ik*ik*6
!              first set up winds in Cartesian "source" coords
               uc(iq)=axs_a(iq)*u_a(iq,k) + bxs_a(iq)*v_a(iq,k)
               vc(iq)=ays_a(iq)*u_a(iq,k) + bys_a(iq)*v_a(iq,k)
               wc(iq)=azs_a(iq)*u_a(iq,k) + bzs_a(iq)*v_a(iq,k)
!              now convert to winds in "absolute" Cartesian components
               ucc(iq)=uc(iq)*rotpoles(1,1)+vc(iq)*rotpoles(1,2)
     &                                +wc(iq)*rotpoles(1,3)
               vcc(iq)=uc(iq)*rotpoles(2,1)+vc(iq)*rotpoles(2,2)
     &                                +wc(iq)*rotpoles(2,3)
               wcc(iq)=uc(iq)*rotpoles(3,1)+vc(iq)*rotpoles(3,2)
     &                                +wc(iq)*rotpoles(3,3)
            enddo               ! iq loop
            if(ktau<3.and.k==1.and.nproc==1)then
               print *,'uc,vc,wc: ',uc(id),vc(idjd),wc(idjd)
               print *,'ucc,vcc,wcc: ',ucc(idjd),vcc(idjd),wcc(idjd)
               print *,'calling ints4 for k= ',k
            endif

!      interpolate all required arrays to new C-C positions
!      don't need to do map factors and Coriolis on target grid
            np=0                ! controls prints in ints4
            call ints4(t_a(1,k),   t_g,nface4,xg4,yg4,nord,ik)  ! ints4 on source grid
            call ints4(qg_a(1,k), qg_g,nface4,xg4,yg4,nord,ik)
            call ints4(ucc,        uct_g, nface4,xg4,yg4,nord,ik)
            call ints4(vcc,        vct_g, nface4,xg4,yg4,nord,ik)
            call ints4(wcc,        wct_g, nface4,xg4,yg4,nord,ik)

c      ********************** N.B. tracers not ready yet
c      if(iltin>1)then
c        do ntr=1,ntracin
c         call ints4(tr(1,k,ntr),nface4,xg4,yg4,nord,ik)
c        enddo
c      endif
 
            do iq=1,ifull_g
!              now convert to "target" Cartesian components (transpose used)
               uct_gg=uct_g(iq)*rotpole(1,1)+vct_g(iq)*rotpole(2,1)
     &                           +wct_g(iq)*rotpole(3,1)
               vct_gg=uct_g(iq)*rotpole(1,2)+vct_g(iq)*rotpole(2,2)
     &                           +wct_g(iq)*rotpole(3,2)
               wct_gg=uct_g(iq)*rotpole(1,3)+vct_g(iq)*rotpole(2,3)
     &                           +wct_g(iq)*rotpole(3,3)
!              then finally to "target" local x-y components
               u_g(iq) = ax_g(iq)*uct_gg + ay_g(iq)*vct_gg +
     &                   az_g(iq)*wct_gg
               v_g(iq) = bx_g(iq)*uct_gg + by_g(iq)*vct_gg +
     &                   bz_g(iq)*wct_gg
            enddo               ! iq loop
            if(ktau<3.and.k==1.and.nproc==1)then
              ! This only works if idjd is on processor 0
              print *,'interp. ucc,vcc,wcc: ',uct_g(idjd),vct_g(idjd),
     &                 wct_g(idjd)
              print *,'uct,vct,wct:',uct_g(idjd),vct_g(idjd),wct_g(idjd)
              print *,'ax,ay,az ',ax_g(idjd),ay_g(idjd),az_g(idjd)
              print *,'bx,by,bz ',bx_g(idjd),by_g(idjd),bz_g(idjd)
              print *,'final u , v: ',u_g(idjd),v_g(idjd)
            endif
            call ccmpi_distribute(u(:,k), u_g)
            call ccmpi_distribute(v(:,k), v_g)
            call ccmpi_distribute(t(:,k), t_g)
            call ccmpi_distribute(qg(:,k), qg_g)
         else ! myid /= 0
c            call ccmpi_gather(u(:,k))
c            call ccmpi_gather(v(:,k))
c            call ccmpi_gather(t(:,k))
c            call ccmpi_gather(qg(:,k))
            call ccmpi_distribute(u(:,k))
            call ccmpi_distribute(v(:,k))
            call ccmpi_distribute(t(:,k))
            call ccmpi_distribute(qg(:,k))
         endif ! myid==0
      enddo  ! k loop

!     below we interpolate quantities which may be affected by land-sea mask
!     set up land-sea mask from either tss or zss
      if(myid==0)then
       !-------------------------------------------
       ! MJT lsmask
       if(nemi==3)then 
         land_a(:)=isoilm_a(:).gt.0
         if (any(isoilm_a(:).lt.0)) nemi=2
       end if
       !-------------------------------------------
       if(nemi==2)then
         numneg=0
         do iq=1,ik*ik*6
            if(tss_a(iq)>0)then ! over land
               land_a(iq)=.true.
            else                ! over sea
               land_a(iq)=.false.
               numneg=numneg+1
            endif               ! (tss(iq)>0) .. else ..
         enddo
         if(numneg==0)nemi=1  ! should be using zss in that case
       endif                     !  (nemi==2)
       print *,'using nemi = ',nemi
      
       if(nemi==1)then
         land_a(:) = zss_a(:) > 0.
       endif                     !  (nemi==1)

       spval=999.
       do iq=1,ik*ik*6
         if(land_a(iq))then       ! over land
            tss_l_a(iq)=tss_a(iq)
            tss_s_a(iq)=spval
            sicedep_a(iq)=spval
            fracice_a(iq)=spval
         else                   ! over sea
            numneg=numneg+1
            tss_s_a(iq)=abs(tss_a(iq))
            tss_l_a(iq)=spval
         endif  !   (land_a(iq)) .. else ..
       enddo     ! iq loop
      
       if(nproc==1.and.nmaxpr==1)then
        print *,'before fill tss ',tss_a(idjd2)
        print *,'before fill tss_l_a, tss_s_a ',
     &                       tss_l_a(idjd2),tss_s_a(idjd2)
        print *,'before fill/ints4 sicedep ',sicedep_a(idjd2)
c        print *,'before fill wb'
c        write(6,"('wb_s(1)#  ',9f7.3)") 
c     .          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
       endif  ! (nproc==1.and.nmaxpr==1)
       call fill_cc(tss_l_a,spval,ik)
       call fill_cc(tss_s_a,spval,ik)
       call fill_cc(sicedep_a,spval,ik)
       call fill_cc(fracice_a,spval,ik)
       if(nproc==1.and.nmaxpr==1)then
        print *,'after fill tss_l, tss_s ',tss_l_a(idjd2),tss_s_a(idjd2)
        print *,'after fill sicedep ',sicedep_a(idjd2)
c        print *,'after fill wb'
c        write(6,"('wb_s(1)#  ',9f7.3)") 
c     .          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
        print *,'before ints4 psl_a(idjd2),zss_a(idjd2) ',
     .                        psl_a(idjd),zss_a(idjd2)
       endif  ! (nproc==1.and.nmaxpr==1)

      endif   ! (myid==0)

!     The routine doints4 does the gather, calls ints4 and redistributes
      call doints4(zss_a ,zss , nface4,xg4,yg4,1,ik)       ! bilinear for zss
      call doints4(pmsl_a,pmsl, nface4,xg4,yg4,nord,ik)
!     invert pmsl to get psl
      call to_pslx(pmsl,psl,zss,t,ifull)  ! on target grid
      call doints4(tss_l_a , tss_l,  nface4,xg4,yg4,nord,ik)
      call doints4(tss_s_a ,tss_s,   nface4,xg4,yg4,nord,ik)
      if(nproc==1.and.nmaxpr==1)then
         print *,'after ints4 idjd,zss(idjd) ',idjd,zss(idjd)
         print *,'zss1-5 ',(zss(iq),iq=1,5)
         print *,'after ints4 psl,pmsl ',psl(idjd),pmsl(idjd)
c        print *,'after ints4 wb_t'
c        write(6,"('wb_t(1)#  ',9f7.3)") 
c     .           ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
      endif  ! (nproc==1.and.nmaxpr==1)
      call doints4(sicedep_a,sicedep,nface4,xg4,yg4,nord,ik)
      call doints4(fracice_a,fracice,nface4,xg4,yg4,nord,ik)
      if ( nproc==1 ) print *,'after ints4 sicedep ',sicedep(idjd)
      if(nested==0)then
       if(myid==0)then
        do iq=1,ik*ik*6
         if(.not.land_a(iq))then       
            snowd_a(iq)=spval
            do k=1,ms
               tgg_a(iq,k)=spval
               wb_a(iq,k)=spval
            enddo
         endif  !   (.not.land_a(iq)) 
        enddo   ! iq loop
        call fill_cc(snowd_a,spval,ik)
        do k=1,ms
         call fill_cc(tgg_a(1,k),spval,ik)
         call fill_cc(wb_a(1,k),spval,ik)
        enddo
       endif  ! (myid==0)
       call doints4(snowd_a,  snowd,nface4,xg4,yg4,nord,ik)
        do k=1,ms
         call doints4(tgg_a(1,k),tgg(1,k),nface4,xg4,yg4,nord,ik)
         call doints4(wb_a(1,k) ,wb(1,k) ,nface4,xg4,yg4,nord,ik)
        enddo
        !--------------------------------------------------
        ! MJT urban
        if (nurban.ne.0) then
          if (myid==0) then
            do k=1,12
              where ((.not.land_a(:)).or.(urban_a(:,k).ge.399.))
                urban_a(:,k)=spval
              end where
              call fill_cc(urban_a(:,k),spval,ik)
            end do
          end if
          do k=1,12
            call doints4(urban_a(:,k),urban(:,k),nface4,xg4,yg4,nord
     &                   ,ik)
          end do
        end if
        !--------------------------------------------------
        !--------------------------------------------------
        ! MJT cable
        if ((nsib.eq.4).or.(nsib.eq.6)) then
          if (myid==0) then
            where (.not.land_a(:))
              rtsoil_a(:)=spval
              cansto_a(:)=spval
            end where
            call fill_cc(rtsoil_a(:),spval,ik)
            call fill_cc(cansto_a(:),spval,ik)
            do k=1,ncp
              where (.not.land_a(:))
                cplant_a(:,k)=spval
              end where
              call fill_cc(cplant_a(:,k),spval,ik)
            end do
            do k=1,ncs
              where (.not.land_a(:))
                csoil_a(:,k)=spval
              end where
              call fill_cc(csoil_a(:,k),spval,ik)
            end do
          end if
          call doints4(rtsoil_a(:),rtsoil(:),nface4,xg4,yg4,nord,ik)
          call doints4(cansto_a(:),cansto(:),nface4,xg4,yg4,nord,ik)
          do k=1,ncp
            call doints4(cplant_a(:,k),cplant(:,k),nface4,xg4,yg4,
     &                   nord,ik)
          end do
          do k=1,ncs
            call doints4(csoil_a(:,k),csoil(:,k),nface4,xg4,yg4,
     &                   nord,ik)
          end do
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

      rlong0x=rlong0  ! just for indata cross-check
      rlat0x=rlat0
      schmidtx=schmidt

!     restore target land array
      land(:) = land_t(:)  
c     write (30+myid,*) 'Z1 myid ',myid
c     close (30+myid)
      end subroutine ontheflyx

      subroutine doints4(s_a,sout,nface4 ,xg4 ,yg4,nord,ik)  ! does calls to intsb
      use cc_mpi
      implicit none
      include 'newmpar.h'
ccc      real, dimension(ik*ik*6), intent(inout) :: s
      real, dimension(ifull), intent(inout) :: sout
      integer, intent(in), dimension(ifull_g,4) :: nface4
      real, intent(in), dimension(ifull_g,4) :: xg4, yg4
      integer, intent(in) :: ik, nord
      real, dimension(ik*ik*6) :: s_a
      real, dimension(ifull_g) :: s_g
c     integer iq

      if ( myid ==0 ) then
ccc         call ccmpi_gather(s,s_a)
         call ints4(s_a,s_g,nface4 ,xg4 ,yg4,nord,ik) 
c        print *,'in doints4 sout1-5',(s_g(iq),iq=1,5)
         call ccmpi_distribute(sout,s_g)
      else
ccc         call ccmpi_gather(s)
         call ccmpi_distribute(sout)
      endif
      end subroutine doints4

      subroutine ints4(s,sout,nface4 ,xg4 ,yg4,nord,ik)  ! does calls to intsb
      use cc_mpi, only : mydiag
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      integer, parameter :: ntest=0
      real, dimension(ik*ik*6), intent(inout) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ifull_g,4) :: nface4
      real, intent(in), dimension(ifull_g,4) :: xg4, yg4
      integer, intent(in) :: ik, nord
      real wrk(ifull_g,4)
      integer :: iq, m
      integer :: idx = 25, jdx = 218, idjdx=10441

      if(nord==1)then
         do m=1,m_fly  !  was 4, now may be 1
            call ints_blb(s,wrk(1,m),nface4(1,m),xg4(1,m),yg4(1,m),ik)
         enddo
      else
         do m=1,m_fly  !  was 4, now may be 1
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
      if(m_fly==1)then
        do iq=1,ifull_g
         sout(iq)=wrk(iq,1)
        enddo
      else
!       average 4 m values to central value
        do iq=1,ifull_g
         sout(iq)=.25*(wrk(iq,1)+wrk(iq,2)+wrk(iq,3)+wrk(iq,4))
        enddo
      endif    ! (m_fly==1)

      end subroutine ints4

      subroutine intsb(s,sout,nface,xg,yg,ik)   ! N.B. sout here
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
c     later may wish to save idel etc between array calls
c     this one does linear interp in x on outer y sides
c     doing x-interpolation before y-interpolation
!     This is a global routine 
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      real, dimension(ik*ik*6), intent(in) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ik*ik*6) :: nface
      integer idel, jdel, ik, nn
      integer :: ind, i, j, n, iq, n_n, n_e, n_w, n_s
      real aaa, c1, c2, c3, c4, xxg, yyg
      real, intent(in), dimension(ik*ik*6) :: xg, yg
      real sx(-1:ik+2,-1:ik+2,0:npanels)
      real r(4)

ch      include 'indices_g.h' ! in,is,iw,ie,inn,iss,iww,iee
      ind(i,j,n)=i+(j-1)*ik+n*ik*ik  ! *** for n=0,npanels
c     this is intsb           EW interps done first
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
      do n=0,npanels
       do j=1,ik
        do i=1,ik
         sx(i,j,n)=s(ind(i,j,n))
        enddo  ! i loop
       enddo   ! j loop
      enddo    ! n loop
      do n=0,npanels,2
       n_w=mod(n+5,6)
       n_e=mod(n+2,6)
       n_n=mod(n+1,6)
       n_s=mod(n+4,6)
       do j=1,ik
        sx(0,j,n)=s(ind(ik,j,n_w))
        sx(-1,j,n)=s(ind(ik-1,j,n_w))
        sx(ik+1,j,n)=s(ind(ik+1-j,1,n_e))
        sx(ik+2,j,n)=s(ind(ik+1-j,2,n_e))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(ind(i,1,n+1))
        sx(i,ik+2,n)=s(ind(i,2,n+1))
        sx(i,0,n)=s(ind(ik,ik+1-i,n_s))
        sx(i,-1,n)=s(ind(ik-1,ik+1-i,n_s))
       enddo
       sx(-1,0,n)=s(ind(ik,2,n_w))       ! wws
       sx(0,-1,n)=s(ind(ik,ik-1,n_s))    ! wss
       sx(0,0,n)=s(ind(ik,1,n_w))        ! ws
       sx(ik+1,0,n)=s(ind(ik,1,n_e))     ! es  
       sx(ik+2,0,n)=s(ind(ik-1,1,n_e))   ! ees 
       sx(-1,ik+1,n)=s(ind(ik,ik-1,n_w)) ! wwn
       sx(0,ik+2,n)=s(ind(ik-1,ik,n_w))  ! wnn
       sx(ik+2,ik+1,n)=s(ind(2,1,n_e))   ! een  
       sx(ik+1,ik+2,n)=s(ind(1,2,n_e))   ! enn  
       sx(0,ik+1,n)   =s(ind(ik,ik,n_w)) ! wn  
       sx(ik+1,ik+1,n)=s(ind(1,1,n_e))   ! en  
       sx(ik+1,-1,n)=s(ind(ik,2,n_e))    ! ess  
      enddo  ! n loop
      do n=1,npanels,2
       n_w=mod(n+4,6)
       n_e=mod(n+1,6)
       n_n=mod(n+2,6)
       n_s=mod(n+5,6)
       do j=1,ik
        sx(0,j,n)=s(ind(ik+1-j,ik,n_w))
        sx(-1,j,n)=s(ind(ik+1-j,ik-1,n_w))
        sx(ik+1,j,n)=s(ind(1,j,n_e))
        sx(ik+2,j,n)=s(ind(2,j,n_e))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(ind(1,ik+1-i,n_n))
        sx(i,ik+2,n)=s(ind(2,ik+1-i,n_n))
        sx(i,0,n)=s(ind(i,ik,n-1))
        sx(i,-1,n)=s(ind(i,ik-1,n-1))
       enddo
       sx(-1,0,n)=s(ind(ik-1,ik,n_w))     ! wws
       sx(0,-1,n)=s(ind(2,ik,n_s))        ! wss
       sx(0,0,n)=s(ind(ik,ik,n_w))        ! ws
       sx(ik+1,0,n)=s(ind(1,1,n_e))       ! es
       sx(ik+2,0,n)=s(ind(1,2,n_e))       ! ees
       sx(-1,ik+1,n)=s(ind(2,ik,n_w))     ! wwn   
       sx(0,ik+2,n)=s(ind(1,ik-1,n_w))    ! wnn  
       sx(ik+2,ik+1,n)=s(ind(1,ik-1,n_e)) ! een  
       sx(ik+1,ik+2,n)=s(ind(2,ik,n_e))   ! enn  
       sx(0,ik+1,n)   =s(ind(1,ik,n_w))   ! wn  
       sx(ik+1,ik+1,n)=s(ind(1,ik,n_e))   ! en  
       sx(ik+1,-1,n)=s(ind(2,1,n_e))     ! ess  
      enddo  ! n loop
c     for ew interpolation, sometimes need (different from ns):
c          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

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
      real, dimension(ik*ik*6), intent(inout) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ik*ik*6) :: nface
      real, intent(in), dimension(ik*ik*6) :: xg, yg
      real sx(-1:ik+2,-1:ik+2,0:npanels)
c     include 'indices_g.h' ! in,is,iw,ie,inn,iss,iww,iee
      integer :: ind, i, j, n, iq, idel, jdel, ik
      integer :: n_n, n_e, n_w, n_s
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
      enddo    ! n loop
      do n=0,npanels,2
       n_w=mod(n+5,6)
       n_e=mod(n+2,6)
       n_n=mod(n+1,6)
       n_s=mod(n+4,6)
       do j=1,ik
        sx(0,j,n)=s(ind(ik,j,n_w))
        sx(ik+1,j,n)=s(ind(ik+1-j,1,n_e))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(ind(i,1,n+1))
        sx(i,0,n)=s(ind(ik,ik+1-i,n_s))
       enddo
       sx(0,0,n)=s(ind(ik,1,n_w))        ! ws
       sx(ik+1,0,n)=s(ind(ik,1,n_e))     ! es  
       sx(0,ik+1,n)   =s(ind(ik,ik,n_w)) ! wn  
       sx(ik+1,ik+1,n)=s(ind(1,1,n_e))   ! en  
      enddo  ! n loop
      do n=1,npanels,2
       n_w=mod(n+4,6)
       n_e=mod(n+1,6)
       n_n=mod(n+2,6)
       n_s=mod(n+5,6)
       do j=1,ik
        sx(0,j,n)=s(ind(ik+1-j,ik,n_w))
        sx(ik+1,j,n)=s(ind(1,j,n_e))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(ind(1,ik+1-i,n_n))
        sx(i,0,n)=s(ind(i,ik,n-1))
       enddo
       sx(0,0,n)=s(ind(ik,ik,n_w))        ! ws
       sx(ik+1,0,n)=s(ind(1,1,n_e))       ! es
       sx(0,ik+1,n)   =s(ind(1,ik,n_w))   ! wn  
       sx(ik+1,ik+1,n)=s(ind(1,ik,n_e))   ! en  
      enddo  ! n loop
      

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

      subroutine fill_cc(a_io,value,ik)
!     this version holds whole array in memory      
c     routine fills in interior of an array which has undefined points
      use cc_mpi
      implicit none
      include 'newmpar.h'
c     include 'indices.h'
      include 'mpif.h'
      real a_io(ik*ik*6)         ! input and output array
      real value            ! array value denoting undefined
c     real b(ik*ik*6), a(ik*ik*6+iextra)
      real b(ik*ik*6), a(ik*ik*6)
      integer :: num, nrem, i, ii, ik, iq, ind, j, n, neighb
      real :: av, avx     
      integer, dimension(ik*ik*6) :: in,ie,iw,is
      integer npann(0:5),npane(0:5),npanw(0:5),npans(0:5)
      data npann/1,103,3,105,5,101/,npane/102,2,104,4,100,0/
      data npanw/5,105,1,101,3,103/,npans/104,0,100,2,102,4/
      ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,npanels
       do iq=1,ik*ik*6
       in(iq)=iq+il
       is(iq)=iq-il
       ie(iq)=iq+1
       iw(iq)=iq-1
      enddo   ! iq loop
      do n=0,npanels
      if(npann(n).lt.100)then
        do ii=1,il
         in(ind(ii,il,n))=ind(ii,1,npann(n))
        enddo    ! ii loop
      else
        do ii=1,il
         in(ind(ii,il,n))=ind(1,il+1-ii,npann(n)-100)
        enddo    ! ii loop
      endif      ! (npann(n).lt.100)
      if(npane(n).lt.100)then
        do ii=1,il
         ie(ind(il,ii,n))=ind(1,ii,npane(n))
        enddo    ! ii loop
      else
        do ii=1,il
         ie(ind(il,ii,n))=ind(il+1-ii,1,npane(n)-100)
        enddo    ! ii loop
      endif      ! (npane(n).lt.100)
      if(npanw(n).lt.100)then
        do ii=1,il
         iw(ind(1,ii,n))=ind(il,ii,npanw(n))
        enddo    ! ii loop
      else
        do ii=1,il
         iw(ind(1,ii,n))=ind(il+1-ii,il,npanw(n)-100)
        enddo    ! ii loop
      endif      ! (npanw(n).lt.100)
      if(npans(n).lt.100)then
        do ii=1,il
         is(ind(ii,1,n))=ind(ii,il,npans(n))
        enddo    ! ii loop
      else
        do ii=1,il
         is(ind(ii,1,n))=ind(il,il+1-ii,npans(n)-100)
        enddo    ! ii loop
      endif      ! (npans(n).lt.100)
      enddo      ! n loop
          
      a(1:ik*ik*6) = a_io(:)
      num=0
      nrem = 1    ! Just for first iteration
c     nrem_gmin = 1 ! Just for first iteration
!     nrem_gmin used to avoid infinite loops, e.g. for no sice
      do while ( nrem > 0)
         ! This has to loop until all are finished otherwise the bounds call
         ! doesn't work.
c808         call bounds(a)
         nrem=0
         num=num+1
         do iq=1,ik*ik*6
            b(iq)=a(iq)
            if(a(iq)==value)then
               neighb=0
               av=0.
               if(a(in(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(in(iq))
               endif
               if(a(ie(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(ie(iq))
               endif
               if(a(iw(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(iw(iq))
               endif
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
         do iq=1,ik*ik*6
            a(iq)=b(iq)
         enddo
c         call MPI_AllReduce(nrem, nrem_g, 1, MPI_INTEGER, MPI_MAX, 
c     &                      MPI_COMM_WORLD, ierr )
c         call MPI_AllReduce(nrem, nrem_gmin, 1, MPI_INTEGER, MPI_MIN, 
c     &                      MPI_COMM_WORLD, ierr )
c        if(nrem_g>0.and.myid==0)then
c         endif                  ! (nrem>0)
      end do
      a_io(1:ik*ik*6) = a(1:ik*ik*6)
      return
      end

      subroutine mslpx(pmsl,psl,zs,t,ifullx)
!     generalized from ifull to allow new onthefly    0808
!     can replace usual mslp sometime
      use cc_mpi, only : mydiag
!     this one will ignore negative zs (i.e. over the ocean)
      implicit none
      integer, parameter :: meth=1 ! 0 for original, 1 for other jlm - always now
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      integer ifullx,iq
      real pmsl(ifullx),psl(ifullx),zs(ifullx),t(ifullx,kl)
      include 'sigs.h'
      integer :: lev
      real c, con, conr, dlnps, phi1, tav, tsurf
      c=grav/stdlapse
      conr=c/rdry
      lev=0
14    lev=lev+1
c     find level just below sig=.9
      if (sig(lev+1).gt..9)go to 14
      con=sig(lev)**(rdry/c)/c
c     if(meth.eq.0)then
c       do iq=1,ifullx
c        pmsl(iq)=ps(iq)*(1.+con*zs(iq)/t(iq,lev))**conr
c       enddo
c     endif  ! (meth.eq.0)
      if(meth.eq.1)then
        do iq=1,ifullx
         phi1=t(iq,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
         tsurf=t(iq,lev)+phi1*stdlapse/grav
         tav=tsurf+max(0.,zs(iq))*.5*stdlapse/grav
         dlnps=max(0.,zs(iq))/(rdry*tav)
         pmsl(iq)=1.e5*exp(psl(iq)+dlnps)
        enddo
      endif  ! (meth.eq.1)
      if(nmaxpr==1.and.mydiag)then
        print *,'meth,lev,sig(lev) ',meth,lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif
      return
      end
      subroutine to_pslx(pmsl,psl,zs,t,ifullx)
!     generalized from ifull to allow new onthefly    0808
!     can replace usual mslp sometime
      use cc_mpi, only : mydiag
!     this one will ignore negative zs (i.e. over the ocean)
      implicit none
      integer, parameter :: meth=1 ! 0 for original, 1 for other jlm - always now
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      integer ifullx,iq
      real pmsl(ifullx),psl(ifullx),zs(ifullx),t(ifullx,kl)
      include 'sigs.h'
      integer :: lev
      real dlnps, phi1, tav, tsurf
      lev=0
14    lev=lev+1
c     find level just below sig=.9
      if (sig(lev+1).gt..9)go to 14
      do iq=1,ifullx
       phi1=t(iq,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
       tsurf=t(iq,lev)+phi1*stdlapse/grav
       tav=tsurf+max(0.,zs(iq))*.5*stdlapse/grav
       dlnps=max(0.,zs(iq))/(rdry*tav)
       psl(iq)=log(1.e-5*pmsl(iq)) -dlnps
      enddo
      if(nmaxpr==1.and.mydiag)then
        print *,'to_psl lev,sig(lev) ',lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif
      return
      end
