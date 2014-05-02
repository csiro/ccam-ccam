      ! Main netcdf input routines.  Host grid is automatically
      ! interpolated to nested model grid.  Three options are
      !   nested=0  Initial conditions
      !   nested=1  Nudging fields
      !   nested=2  Surface data recycling
      
      ! This version supports the parallel file routines contained
      ! in infile.f90.  Hence, restart files do not require any
      ! gathers and scatters.

      ! In the case where the grid needs to be interpolated, a copy
      ! of the input data is sent to all processors and each
      ! processor performs its own interpolation.
      
      subroutine onthefly(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,
     &                    fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg,
     &                    qlg,qrg,tggsn,smass,ssdn,ssdnn,snage,
     &                    isflag,mlodwn,ocndwn)

      use cc_mpi           ! CC MPI routines
      use infile           ! Input file routines
      use mlo              ! Ocean physics and prognostic arrays
      use soil_m           ! Soil and surface data

      implicit none

      include 'newmpar.h'  ! Grid parameters
      include 'darcdf.h'   ! Netcdf data
      include 'parm.h'     ! Model configuration
      include 'stime.h'    ! File date data

      integer, parameter :: nihead=54
      integer, parameter :: nrhead=14

      integer, save :: ik,jk,kk,ok,maxarchi,nsibx
      integer, save :: ncidold = -1
      integer kdate_r,ktime_r,nested,ier,ier2,ilen,itype
      integer idv,mtimer,k,ierx,idvkd,idvkt,idvmt
      integer, dimension(nihead) :: nahead
      integer, dimension(ifull) :: isflag
      real, save :: rlong0x, rlat0x, schmidtx
      real timer
      real, dimension(nrhead) :: ahead
      real, dimension(14) :: rdum
!     These are local arrays, not the versions in arrays.h
!     Use in call to infile, so are dimensioned ifull rather than ifull_g
      real, dimension(ifull) :: psl,zss,tss,fracice,snowd,sicedep
      real, dimension(ifull) :: ssdnn,snage
      real, dimension(ifull,ms) :: wb,wbice,tgg
      real, dimension(ifull,3) :: tggsn,smass,ssdn
      real, dimension(ifull,kl) :: t,u,v,qg,qfg,qlg,qrg
      real, dimension(ifull,wlev,4) :: mlodwn
      real, dimension(ifull,2) :: ocndwn
      logical ltest,newfile,tst

#include "log.h"

      START_LOG(onthefly)
      !--------------------------------------------------------------
      ! pfall indicates all processors have an input file and there
      ! is no need to broadcast metadata (see infile.f90)
      if ( myid==0 .or. pfall )then
        if (myid==0) then
          write(6,*) 'Entering onthefly for nested,ktau = ',
     &                                      nested,ktau
        end if
        if(ncid/=ncidold)then
          if (myid==0) then
            write(6,*) 'Reading new file metadata'
          end if
          iarchi=1
          call ccnf_get_attg(ncid,'int_header',nahead)
          call ccnf_get_attg(ncid,'real_header',ahead)
          ik      =nahead(1)
          jk      =nahead(2)
          kk      =nahead(3)
          nsibx   =nahead(44)
          rlong0x =ahead(5)
          rlat0x  =ahead(6)
          schmidtx=ahead(7)
          if(schmidtx<=0..or.schmidtx>1.)then
            ! backwards compatibility option
            rlong0x =ahead(6)
            rlat0x  =ahead(7)
            schmidtx=ahead(8)
          endif  ! (schmidtx<=0..or.schmidtx>1.)        
          maxarchi=0
          call ccnf_inq_dimlen(ncid,'time',maxarchi)
          ok=0
          call ccnf_inq_dimlen(ncid,'olev',ok,failok=.true.)
          if (myid==0) then
            write(6,*) "Found ik,jk,kk,ok ",ik,jk,kk,ok
            write(6,*) "      maxarchi ",maxarchi
            write(6,*) "      rlong0x,rlat0x,schmidtx ",
     &                        rlong0x,rlat0x,schmidtx
          end if
        end if
        if (myid==0) then
          write(6,*)'Search for kdate_s,ktime_s >= ',
     &                          kdate_s,ktime_s
        end if
        ltest=.true.
        iarchi=iarchi-1
        call ccnf_inq_varid(ncid,'kdate',idvkd,tst)
        call ccnf_inq_varid(ncid,'ktime',idvkt,tst)
        call ccnf_inq_varid(ncid,'mtimer',idvmt,tst)
        if (tst) then
          ierx=1
          call ccnf_inq_varid(ncid,'timer',idvmt,tst)
        else
          ierx=0
        end if
        do while(ltest.and.iarchi<maxarchi)
          ! could read this as one array, but we only usually need to advance 1 step
          iarchi=iarchi+1
          call ccnf_get_var1(ncid,idvkd,iarchi,kdate_r)
          call ccnf_get_var1(ncid,idvkt,iarchi,ktime_r)
          if (ierx==0) then
            call ccnf_get_var1(ncid,idvmt,iarchi,mtimer)
            timer=mtimer/60.
          else
            timer=0.
            call ccnf_get_var1(ncid,idvmt,iarchi,timer)
            mtimer=nint(timer*60.)
          endif
          if (mtimer>0) then
            call datefix(kdate_r,ktime_r,mtimer)
          end if
          ! ltest = .false. when correct date is found
          ltest=2400*(kdate_r-kdate_s)-1200*nsemble
     &              +(ktime_r-ktime_s)<0
        end do
        if (nsemble/=0) then
          kdate_r=kdate_s
          ktime_r=ktime_s
        end if
        if (ltest) then
          ! ran out of file before correct date was located
          ktime_r=-1
        end if
        if (myid==0) then
          write(6,*) 'After search ltest,iarchi =',ltest,iarchi
          write(6,*) '             kdate_r,ktime_r =',kdate_r,ktime_r
        end if
        ! store metadata for possible bcast
        rdum(1)=rlong0x
        rdum(2)=rlat0x
        rdum(3)=schmidtx
        ! kdate_r is too large to represent as a single real, so
        ! we split kdate_r into year, month and day
        rdum(4)=real(kdate_r/10000)
        rdum(5)=real(kdate_r/100-nint(rdum(4))*100)
        rdum(6)=real(kdate_r-nint(rdum(4))*10000-nint(rdum(5))*100)
        rdum(7)=real(ktime_r)
        if (ncid/=ncidold) then
          rdum(8)=1.
        else
          rdum(8)=0.
        end if
        rdum(9) =real(ik)
        rdum(10)=real(jk)
        rdum(11)=real(kk)
        rdum(12)=real(ok)
        rdum(13)=real(iarchi)
        rdum(14)=real(nsibx)
      endif  ! ( myid==0 .or. pfall )

      ! if metadata is not read by all processors, then broadcast
      if (.not.pfall) then
        call ccmpi_bcast(rdum(1:14),0,comm_world)
        rlong0x =rdum(1)
        rlat0x  =rdum(2)
        schmidtx=rdum(3)
        kdate_r =nint(rdum(4))*10000+nint(rdum(5))*100+nint(rdum(6))
        ktime_r =nint(rdum(7))
        newfile =(nint(rdum(8))==1)
        ik      =nint(rdum(9))
        jk      =nint(rdum(10))
        kk      =nint(rdum(11))
        ok      =nint(rdum(12))
        iarchi  =nint(rdum(13))
        nsibx   =nint(rdum(14))
      else
        newfile=(ncid/=ncidold)            
      end if
      
      ! close old file if a new file is opened
      if (newfile) then
        if (ncidold/=-1) then
          if (myid==0) then
            write(6,*) 'Closing old input file'
          end if
          call histclose
        end if
        ncidold=ncid
      end if

      ! trap error if correct date/time is not located      
      if (ktime_r<0) then
        if (nested==2) then
          if (myid==0) then
            write(6,*) "WARN: Cannot locate date/time in input file"
          end if
          return
        else
          write(6,*) "ERROR: Cannot locate date/time in input file"
          call ccmpi_abort(-1)
        end if
      end if
      !--------------------------------------------------------------
      
      ! Here we call ontheflyx with different automatic array
      ! sizes.  This means the arrays are correct for interpolation
      ! and file i/o on myid==0, as well as the arrays are smaller
      ! on myid/=0 when they are not needed.  This way we avoid
      ! having to maintain multiple ontheflyx subroutines.
      
      ! Note that if histrd fails to find a variable, it returns
      ! zero in the output array
      
      if (myid==0) then
        call ontheflyx(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,
     &                 fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg,qlg,
     &                 qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,ik,
     &                 kk,ok,ik,mlodwn,ocndwn,rlong0x,rlat0x,
     &                 schmidtx,nsibx,newfile)
        write(6,*) "Leaving onthefly"
      else
        call ontheflyx(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,
     &                 fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg,qlg,
     &                 qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,ik,
     &                 kk,ok,0,mlodwn,ocndwn,rlong0x,rlat0x,
     &                 schmidtx,nsibx,newfile)
      end if

      END_LOG(onthefly)

      return
      end
      
      ! Read data from netcdf file
      
      ! arrays are typically read as global and then distributed to
      ! processor local arrays.  This allows for more flexibility
      ! with diagnosed fields.  Data is usually read in as 2D
      ! fields which avoids memory problems when the host grid
      ! size is significantly larger than the regional grid size.
      subroutine ontheflyx(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,
     &                    fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg,qlg,
     &                    qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,ik,
     &                    kk,ok,dk,mlodwn,ocndwn,rlong0x,rlat0x,
     &                    schmidtx,nsibx,newfile)
      
      use aerosolldr, only : xtg,ssn,naero           ! LDR aerosol scheme
      use ateb, only : atebdwn                       ! Urban
      use cable_def_types_mod, only : ncs, ncp       ! CABLE dimensions
      use casadimension, only : mplant,mlitter,msoil ! CASA dimensions
      use carbpools_m                                ! Carbon pools
      use cc_mpi                                     ! CC MPI routines
      use cfrac_m                                    ! Cloud fraction
      use cloudmod                                   ! Prognostic strat cloud
      use extraout_m                                 ! Additional diagnostics      
      use infile                                     ! Input file routines
      use latlong_m                                  ! Lat/lon coordinates
      use mlo, only : wlev,micdwn,mloregrid          ! Ocean physics and prognostic arrays
      use mlodynamics                                ! Ocean dynamics
      use morepbl_m                                  ! Additional boundary layer diagnostics
      use nharrs_m, only : phi_nh,lrestart           ! Non-hydrostatic atmosphere arrays
      use nsibd_m, only : isoilm                     ! Land-surface arrays
      use river                                      ! River routing
      use savuvt_m                                   ! Saved dynamic arrays
      use savuv1_m                                   ! Saved dynamic arrays
      use screen_m                                   ! Screen level diagnostics
      use sigs_m                                     ! Atmosphere sigma levels
      use soil_m                                     ! Soil and surface data
      use tkeeps, only : tke,eps,zidry               ! TKE-EPS boundary layer
      use tracers_m                                  ! Tracer data
      use utilities                                  ! Grid utilities
      use vecsuv_m                                   ! Map to cartesian coordinates
      use vvel_m, only : dpsldt,sdot                 ! Additional vertical velocity
      use xarrs_m, only : pslx                       ! Saved dynamic arrays
      use workglob_m                                 ! Additional grid interpolation
      use work2_m                                    ! Diagnostic arrays

      implicit none

      include 'newmpar.h'                            ! Grid parameters
      include 'const_phys.h'                         ! Physical constants
      include 'darcdf.h'                             ! Netcdf data
      include 'kuocom.h'                             ! Convection parameters
      include 'parm.h'                               ! Model configuration
      include 'parmdyn.h'                            ! Dynamics parmaters
      include 'parmgeom.h'                           ! Coordinate data
      include 'soilv.h'                              ! Soil parameters
      include 'stime.h'                              ! File date data

      integer, parameter :: nord=3        ! 1 for bilinear, 3 for bicubic interpolation
      real, parameter :: iotol=1.E-5      ! tolarance for iotest grid matching
      
      integer ik, kk, ok, idv, isoil, nud_test
      integer dk ! controls automatic array size
      integer lev, levkk, ier, ierr, igas
      integer kdate_r, ktime_r, nemi, id2,jd2,idjd2
      integer nested, i, j, k, mm, iq, ii, jj, np, numneg
      integer nsibx
      integer, dimension(:,:), allocatable, save :: nface4
      integer, dimension(:), allocatable, save :: isoilm_a
      integer, dimension(ifull) :: isflag
      integer, dimension(7+3*ms) :: ierc
      integer, dimension(3), save :: iers
      integer, dimension(2) :: dumb
      real(kind=8), dimension(:,:), allocatable, save :: xx4,yy4
      real(kind=8), dimension(dk*dk*6):: z_a,x_a,y_a
      real, dimension(ifull,wlev,4) :: mlodwn
      real, dimension(ifull,ok,2) :: mloin
      real, dimension(ifull,2) :: ocndwn
      real, dimension(ifull,ms) :: wb,wbice,tgg
      real, dimension(ifull,3) :: tggsn,smass,ssdn
      real, dimension(ifull,kl) :: t,u,v,qg,qfg,qlg,qrg
      real, dimension(ifull,kl) :: dum
      real, dimension(ifull,kk) :: u_k,v_k
      real, dimension(3,3), save :: rotpoles,rotpole
      real, dimension(ifull) :: psl,zss,tss,fracice
      real, dimension(ifull) :: snowd,sicedep,ssdnn,snage,dum6
      real, dimension(ifull) :: tss_l, tss_s, pmsl
      real, dimension(ik*ik*6) :: fracice_a,sicedep_a
      real, dimension(ik*ik*6) :: tss_l_a,tss_s_a
      real, dimension(ik*ik*6) :: ucc,vcc,pmsl_a
      real, dimension(dk*dk*6) :: t_a_lev,psl_a,tss_a
      real, dimension(dk*dk*6) :: wts_a  ! not used here or defined in call setxyz
      real, dimension(:,:), allocatable, save :: xg4, yg4
      real, dimension(:), allocatable, save :: sigin,zss_a,ocndep_l
      real, dimension(:), allocatable, save :: axs_a,ays_a,azs_a
      real, dimension(:), allocatable, save :: bxs_a,bys_a,bzs_a
      real, dimension(kk+3) :: dumr
      real, intent(in) :: rlong0x, rlat0x, schmidtx
      real rlongd, rlatd
      character(len=8) vname
      character(len=3) trnum
      logical, dimension(:), allocatable, save :: land_a,sea_a
      logical iotest,newfile,tsstest,tst

      ! land-sea mask method (nemi=3 use soilt, nemi=2 use tgg, nemi=1 use zs)
      nemi=3
      
      ! test if retopo fields are required
      nud_test=1
      if (nud_p==0.and.nud_t==0.and.nud_q==0) nud_test=0
      
      ! Determine if interpolation is required
      iotest=6*ik*ik==ifull_g           .and.
     &       abs(rlong0x-rlong0)<iotol  .and.
     &       abs(rlat0x-rlat0)<iotol    .and.
     &       abs(schmidtx-schmidt)<iotol.and.
     &       nsib==nsibx
      if (iotest) then
        io_in=1   ! no interpolation
      else
        io_in=-1  ! interpolation
      end if
      if (myid==0) then
        write(6,*) "Interpolation iotest,io_in =",iotest,io_in
      end if

      if (.not.allocated(nface4)) then
        allocate(nface4(ifull,4),xg4(ifull,4),yg4(ifull,4))
      end if
      if (newfile) then
        if (allocated(land_a)) deallocate(sigin,land_a,sea_a)
        allocate(land_a(dk*dk*6),sea_a(dk*dk*6))
        allocate(sigin(kk))        
      end if
      
      !--------------------------------------------------------------
      ! Determine input grid coordinates and interpolation arrays
      if (newfile.and..not.iotest) then
       if (allocated(axs_a)) then
         deallocate(axs_a,ays_a,azs_a)
         deallocate(bxs_a,bys_a,bzs_a)
       end if
       allocate(axs_a(dk*dk*6),ays_a(dk*dk*6),azs_a(dk*dk*6))
       allocate(bxs_a(dk*dk*6),bys_a(dk*dk*6),bzs_a(dk*dk*6))
       ! xx4 and yy4 could be replaced with sharded arrays in MPI-3
       allocate(xx4(1+4*ik,1+4*ik),yy4(1+4*ik,1+4*ik))

       if (m_fly==1) then
         rlong4_l(:,1)=rlongg(:)*180./pi
         rlat4_l(:,1)=rlatt(:)*180./pi
       end if    
          
       if ( myid==0 ) then
        write(6,*) "Defining input file grid"
!       N.B. -ve ik in call setxyz preserves TARGET rlat4, rlong4     
!       following setxyz call is for source data geom    ****   
        do iq=1,ik*ik*6
         axs_a(iq)=iq
         ays_a(iq)=iq
         azs_a(iq)=iq
        enddo      
        call setxyz(ik,rlong0x,rlat0x,-schmidtx,x_a,y_a,z_a,wts_a,
     &   axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4,myid)
       end if ! (myid==0)

       call ccmpi_bcastr8(xx4,0,comm_world)
       call ccmpi_bcastr8(yy4,0,comm_world)
     
       rotpoles = calc_rotpole(rlong0x,rlat0x)
       if(ktau<3.and.myid==0)then
          write(6,*)'m_fly,nord ',m_fly,nord
          write(6,*)'kdate_r,ktime_r,ktau,ds',
     &             kdate_r,ktime_r,ktau,ds
          write(6,*)'rotpoles:'
          do i=1,3
             write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)')
     &         (i,j,j=1,3),(rotpoles(i,j),j=1,3)
          enddo
       endif                  ! (ktau<3.and.myid==0)

!      rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!      rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!      rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
       rotpole = calc_rotpole(rlong0,rlat0)
       if(nmaxpr==1.and.myid==0)then
          write(6,*)'in onthefly rotpole:'
          do i=1,3
             write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)')
     &         (i,j,j=1,3),(rotpole(i,j),j=1,3)
          enddo
          write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
          write(6,*)'before latltoij for id,jd: ',id,jd
          write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,
     &               schmidtx
       endif                  ! (nmaxpr==1.and.myid==0)

       do mm=1,m_fly  !  was 4, now may be set to 1 in namelist
         do iq=1,ifull
           call latltoij(rlong4_l(iq,mm),rlat4_l(iq,mm),      !input
     &                   rlong0x,rlat0x,schmidtx,             !input
     &                   xg4(iq,mm),yg4(iq,mm),nface4(iq,mm), !output (source)
     &                   xx4,yy4,ik)
         enddo
       enddo
       deallocate(xx4,yy4)
       
       if(nproc==1.and.nmaxpr==1)then
         ! Diagnostics will only be correct if nproc==1
         id2=nint(xg4(idjd,1))
         jd2=il*nface4(idjd,1)+nint(yg4(idjd,1))
         idjd2=id2+il*(jd2-1)
         write(6,*)'after latltoij giving id2,jd2,idjd2: ',
     &                                  id2,jd2,idjd2
         write(6,*)'nface4(1-4) ',(nface4(idjd,mm),mm=1,4)
         write(6,*)'xg4(1-4) ',(xg4(idjd,mm),mm=1,4)
         write(6,*)'yg4(1-4) ',(yg4(idjd,mm),mm=1,4)
         if(nested==0)then
            write(6,"('wb_s(1)#  ',9f7.3)") 
     &          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
            write(6,"('wb_s(ms)# ',9f7.3)") 
     &          ((wb(ii+(jj-1)*il,ms),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
         endif  ! (nested==0)
       endif
      end if ! newfile .and. .not.iotest
      
       
      ! read time invariant data when file is first opened
      ! need global zss_a for (potentially) landsea mask and psl interpolation
      ! need global isoilm_a for (potentially) landsea mask
      if (newfile) then
       if (myid==0.or.pfall) then
         if (myid==0) then
           write(6,*) "Reading time invariant fields"
         end if
         call ccnf_inq_varid(ncid,'lev',idv,tst)
         if (tst) then
           call ccnf_inq_varid(ncid,'layer',idv,tst)
         end if
         if (tst) then
           call ccnf_get_attg(ncid,'sigma',sigin)
         else
           dumb(1)=1
           dumb(2)=kk
           call ccnf_get_vara(ncid,idv,dumb(1:1),dumb(2:2),sigin)
         end if
         if (myid==0) then
           write(6,'("sigin=",(9f7.4))') (sigin(k),k=1,kk)
         end if
         ! check for missing data
         iers(1:3)=0
         call ccnf_inq_varid(ncid,'mixr',idv,tst)
         if (tst) iers(1)=-1
         call ccnf_inq_varid(ncid,'siced',idv,tst)
         if (tst) iers(2)=-1
         call ccnf_inq_varid(ncid,'fracice',idv,tst)
         if (tst) iers(3)=-1
       end if
       ! bcast data to all processors unless all processes are reading input files
       if (.not.pfall) then
         dumr(1:kk)     =sigin
         dumr(kk+1:kk+3)=real(iers(1:3))
         call ccmpi_bcast(dumr(1:kk+3),0,comm_world)
         sigin    =dumr(1:kk)
         iers(1:3)=nint(dumr(kk+1:kk+3))
       end if
       ! determine whether surface temperature needs to be interpolated (tsstest=.false.)
       tsstest=(iers(2)==0.and.iers(3)==0.and.iotest)
       if (myid==0) then
         write(6,*) "tsstest,iers ",tsstest,iers(1:3)
       end if      
       if (allocated(zss_a)) deallocate(zss_a)
       if (allocated(isoilm_a)) deallocate(isoilm_a)
       if (tsstest) then
         ! load local surface temperature
         allocate(zss_a(ifull))
         call histrd1(ncid,iarchi,ier,'zht',ik,6*ik,zss_a,ifull)
       else
         ! load global surface temperature   
         allocate(zss_a(6*ik*ik)) ! zss_a is allocated for all processors since this memory will
                                  ! be used for a copy of zss_a on all processors when calling
                                  ! doints4.  We could use ucc and copy ucc=zss_a before doints4,
                                  ! if memory becomes an issue, but at the cost of extra copying
         if (myid==0) then
           allocate(isoilm_a(6*ik*ik))
           zss_a=0.
           isoilm_a=-1
         end if
         call histrd1(ncid,iarchi,ier,'zht',ik,6*ik,zss_a,6*ik*ik)
         call histrd1(ncid,iarchi,ier,'soilt',ik,6*ik,ucc,6*ik*ik)
         if (myid==0) then
           if (all(ucc==0.)) ucc=-1.
           isoilm_a=nint(ucc)
         end if
       end if
       if (nmlo/=0.and.abs(nmlo)<=9) then
         if (.not.allocated(ocndep_l)) allocate(ocndep_l(ifull))
         call gethist1(ncid,iarchi,'ocndepth',ik,ocndep_l,iotest,
     &                 nface4,xg4,yg4,nord)
       end if
       if (myid==0) then
         write(6,*) "Finished reading fixed fields"
       end if
      else
       if (myid==0) then
         write(6,*) "Using saved fixed fields"
       end if
       tsstest=(iers(2)==0.and.iers(3)==0.and.iotest)        
      endif ! newfile ..else..
      
      ! detemine the level below sig=0.9 (used to calculate psl)
      lev=0
      do while(sig(lev+1)>0.9) ! nested grid
        lev=lev+1
      end do
      levkk=0
      do while(sigin(levkk+1)>0.9) ! host grid
        levkk=levkk+1
      end do      
      if (myid==0) write(6,*) "Ref height lev,levkk =",
     &                                    lev,levkk

      !--------------------------------------------------------------
      ! Begin reading host data for current time step
      ! psf read when nested=0 or nested=1.and.nud_p/=0
      psl_a=0.
      psl=0.
      if (nested==0.or.(nested==1.and.nud_test/=0)) then
        if (iotest) then
          call histrd1(ncid,iarchi,ier,'psf',ik,6*ik,psl,ifull)
        else
          call histrd1(ncid,iarchi,ier,'psf',ik,6*ik,psl_a,6*ik*ik)
        end if
      endif
      
      ! Read surface temperature 
      ! read global tss to diagnose sea-ice or land-sea mask
      if (tsstest) then
        call histrd1(ncid,iarchi,ier,'tsu',ik,6*ik,tss,ifull)
        zss=zss_a ! used saved zss arrays
      else
        call histrd1(ncid,iarchi,ier,'tsu',ik,6*ik,tss_a,6*ik*ik)
      
        ! set up land-sea mask from either soilt, tss or zss
        if (newfile) then
         if (myid==0) then
          if (nemi==3) then 
           land_a(:)=isoilm_a(:)>0
           numneg=count(.not.land_a)
           if (any(isoilm_a(:)<0)) nemi=2
          end if ! (nemi==3)
          if (nemi==2) then
           numneg=0
           do iq=1,ik*ik*6
             if(tss_a(iq)>0)then ! over land
               land_a(iq)=.true.
             else                ! over sea
               land_a(iq)=.false.
               numneg=numneg+1
             endif               ! (tss(iq)>0) .. else ..
           enddo
           if (numneg==0) nemi=1  ! should be using zss in that case
          endif !  (nemi==2)
          tss_a=abs(tss_a)
          if(nemi==1)then
           land_a(:)=zss_a(:)>0.
           numneg=count(.not.land_a)
          endif ! (nemi==1)
          write(6,*)'Land-sea mask using nemi = ',nemi
          sea_a=.not.land_a
          deallocate(isoilm_a)
         end if ! (myid==0)
        end if ! (newfile)
      end if ! (tsstest) ..else..

      
      !--------------------------------------------------------------
      ! Read ocean data for nudging (sea-ice is read below)
      ! read when nested=0 or nested==1.and.nud/=0 or nested=2
      if (nmlo/=0.and.abs(nmlo)<=9) then
        ! fixed ocean depth
        ocndwn(:,1)=ocndep_l
        ! ocean potential temperature
        ! ocean temperature and soil temperature use the same arrays
        ! as no fractional land or sea cover is allowed in CCAM
        if ((nested/=1.or.nud_sst/=0).and.ok>0) then
          do k=1,ok
            if (k<=ms) then
              write(vname,'("tgg",I1.1)') k
            else
              write(vname,'("tgg",I2.2)') k
            end if
            call filhist1(ncid,iarchi,vname,ik,dk,mloin(:,k,1),iotest,
     &                    nface4,xg4,yg4,nord,land_a)
          end do
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,1),0)
          if (all(mlodwn(:,:,1)==0.)) mlodwn(:,:,1)=293.
        else
          mlodwn(:,:,1)=293.
        end if ! (nestesd/=1.or.nud_sst/=0) ..else..
        ! ocean salinity
        if ((nested/=1.or.nud_sss/=0).and.ok>0) then
          do k=1,ok
            write(vname,'("sal",I2.2)') k
            call filhist1(ncid,iarchi,vname,ik,dk,mloin(:,k,1),iotest,
     &                    nface4,xg4,yg4,nord,land_a)
          end do
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,2),1)
          mlodwn(:,:,2)=max(mlodwn(:,:,2),0.)
        else
          mlodwn(:,:,2)=34.72   
        end if ! (nestesd/=1.or.nud_sss/=0) ..else..
        ! ocean currents
        if ((nested/=1.or.nud_ouv/=0).and.ok>0) then
          do k=1,ok
            if (iotest) then
              write(vname,'("uoc",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     mloin(:,k,1),ifull)
              write(vname,'("voc",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     mloin(:,k,2),ifull)
            else
              ucc=0.
              vcc=0.
              write(vname,'("uoc",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              write(vname,'("voc",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     vcc,6*ik*ik)
              call fill_cc(ucc,dk,0,land_a)
              call fill_cc(vcc,dk,0,land_a)
              call interpwind(ik,mloin(:,k,1),mloin(:,k,2),ucc,vcc,
     &                        axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,
     &                        rotpole,rotpoles,nface4,xg4,yg4,nord,
     &                        dk)
!               interpolate all required arrays to new C-C positions
            end if ! iotest
          end do
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,3),2)
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,2),mlodwn(:,:,4),3)
        else
          mlodwn(:,:,3:4)=0.               
        end if ! (nestesd/=1.or.nud_ouv/=0) ..else..
        ! water surface height
        if (nested/=1.or.nud_sfh/=0) then
          call filhist1(ncid,iarchi,'ocheight',ik,dk,ocndwn(:,2),
     &                  iotest,nface4,xg4,yg4,nord,land_a)
        else
          ocndwn(:,2)=0.
        end if ! (nested/=1.or.nud_sfh/=0) ..else..
      end if
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      ! read sea ice here for prescribed SSTs configuration and for
      ! mixed-layer-ocean
      if (tsstest) then
        call histrd1(ncid,iarchi,ier,'siced',ik,6*ik,sicedep,ifull)
        call histrd1(ncid,iarchi,ier,'fracice',ik,6*ik,fracice,ifull)
      else
        call histrd1(ncid,iarchi,ier,'siced',ik,6*ik,sicedep_a,6*ik*ik)
        call histrd1(ncid,iarchi,ier,'fracice',ik,6*ik,fracice_a,
     &               6*ik*ik)
        if (myid==0) then
         if(iers(2)==0)then  ! i.e. sicedep read in 
          if(iers(3)/=0)then ! i.e. sicedep read in; fracice not read in
            where(sicedep_a>0.)
              fracice_a=1.
            endwhere
          endif  ! (ierr/=0)  fracice
         else     ! sicedep not read in
          if(iers(3)/=0)then  ! neither sicedep nor fracice read in
            sicedep_a(:)=0.  ! Oct 08
            fracice_a(:)=0.
            write(6,*)'pre-setting siced in onthefly from tss'
            where(abs(tss_a) <= 271.6) ! for ERA-Interim
              sicedep_a=1.  ! Oct 08   ! previously 271.2
              fracice_a=1.
            endwhere
          else  ! i.e. only fracice read in;  done in indata, nestin
c***        but needed here for onthefly (different dims) 28/8/08        
            where (fracice_a>.01)
              sicedep_a=2.
            elsewhere
              sicedep_a=0.
              fracice_a=0.
            endwhere
          endif  ! (iers(3)/=0)
         endif    ! (iers(2)/=0) .. else ..    for sicedep
         
         ! interpolate surface temperature and sea-ice
         tss_l_a=abs(tss_a)
         tss_s_a=abs(tss_a)
         call fill_cc(tss_l_a,dk,0,sea_a)
         call fill_cc(tss_s_a,dk,0,land_a)
         call fill_cc(sicedep_a,dk,0,land_a)
         call fill_cc(fracice_a,dk,0,land_a)
        end if ! myid==0

        if (iotest) then
          ! This case occurs for missing sea-ice data
          if (myid==0) then
            call ccmpi_distribute(zss,zss_a)
            call ccmpi_distribute(tss_l,tss_l_a)
            call ccmpi_distribute(tss_s,tss_s_a)
            call ccmpi_distribute(sicedep,sicedep_a)
            call ccmpi_distribute(fracice,fracice_a)
          else
            call ccmpi_distribute(zss)
            call ccmpi_distribute(tss_l)
            call ccmpi_distribute(tss_s)
            call ccmpi_distribute(sicedep)
            call ccmpi_distribute(fracice)
          end if
!         incorporate other target land mask effects
          where (land)
            sicedep=0.
            fracice=0.
            tss=tss_l
          elsewhere
            tss=tss_s
          end where
        else
!         The routine doints4 does the gather, calls ints4 and redistributes
          call doints4(zss_a,    zss,    nface4,xg4,yg4,nord,ik)
          call doints4(tss_l_a,  tss_l,  nface4,xg4,yg4,nord,ik)
          call doints4(tss_s_a,  tss_s,  nface4,xg4,yg4,nord,ik)
          call doints4(fracice_a,fracice,nface4,xg4,yg4,nord,ik)
          call doints4(sicedep_a,sicedep,nface4,xg4,yg4,nord,ik)
!         incorporate other target land mask effects
          do iq=1,ifull
            if(land(iq))then
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
          if (nproc==1.and.nmaxpr==1) then
            write(6,*)'after ints4 idjd,zss(idjd) ',idjd,zss(idjd)
            write(6,*)'zss1-5 ',(zss(iq),iq=1,5)
            write(6,*)'after ints4 psl,pmsl ',psl(idjd),pmsl(idjd)
          end if  ! (nproc==1.and.nmaxpr==1)
        end if ! iotest
      end if ! (tsstest) ..else..

      ! to be depeciated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (nspecial==44.or.nspecial==46) then
       do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if (rlatd>=-43..and.rlatd<=-30.) then
         if (rlongd>=155..and.rlongd<=170.) then
          tss(iq)=tss(iq)+1.
         end if
        end if
       end do
      end if
      if (nspecial==45.or.nspecial==46) then
       do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if (rlatd>=-15..and.rlatd<=-5.) then
         if (rlongd>=150..and.rlongd<=170.) then
          tss(iq)=tss(iq)+1.
         end if
        end if
       end do
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! read atmospheric fields for nested=0 or nested=1.and.nud/=0

      ! air temperature
      ! read for nested=0 or nested=1.and.(nud_t/=0.or.nud_p/=0)
      t=300.
      if (nested==0.or.(nested==1.and.nud_test/=0)) then
        do k=1,kk
          if (iotest) then
            call histrd4s(ncid,iarchi,ier,'temp',ik,6*ik,k,u_k(:,k),
     &                    ifull)                                        !     temperature
          else
            call histrd4s(ncid,iarchi,ier,'temp',ik,6*ik,k,ucc,6*ik*ik) !     temperature
            if (k==levkk.and.myid==0) t_a_lev=ucc ! store for psl calculation below
            call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,ik)  ! ints4 on source grid
          end if ! iotest
        enddo  ! k loop
        call vertint(u_k ,t, 1,kk,sigin)
      end if ! (nested==0.or.(nested==1.and.(nud_t/=0.or.nud_p/=0)))
      ! winds
      ! read for nested=0 or nested=1.and.nud_uv/=0
      u=0.
      v=0.
      if (nested==0.or.(nested==1.and.nud_uv/=0)) then
        do k=1,kk
          ! to reduce memory footprint, we now have to alternatively read
          ! u and v.  This is a bit inefficent for disk accessing,
          ! but makes it possible to downscale large grids
          if (iotest) then
            call histrd4s(ncid,iarchi,ier,'u',ik,6*ik,k,u_k(:,k),
     &                    ifull)                                        !     u wind component
            call histrd4s(ncid,iarchi,ier,'v',ik,6*ik,k,v_k(:,k),
     &                    ifull)                                        !     v wind component
          else
            call histrd4s(ncid,iarchi,ier,'u',ik,6*ik,k,ucc,6*ik*ik)    !     u wind component
            call histrd4s(ncid,iarchi,ier,'v',ik,6*ik,k,vcc,6*ik*ik)    !     v wind component
            call interpwind(ik,u_k(:,k),v_k(:,k),ucc,vcc,axs_a,ays_a,
     &                      azs_a,bxs_a,bys_a,bzs_a,rotpole,
     &                      rotpoles,nface4,xg4,yg4,nord,dk)
!           interpolate all required arrays to new C-C positions
!           don't need to do map factors and Coriolis on target grid
          end if ! iotest
        enddo  ! k loop
        call vertint(u_k ,u, 3,kk,sigin)
        call vertint(v_k ,v, 4,kk,sigin)
      end if ! (nested==0.or.(nested==1.and.nud_uv/=0))
      ! mixing ratio
      ! read for nested=0 or nested=1.and.nud_q/=0
      if (nested==0.or.(nested==1.and.nud_q/=0)) then
        if (iers(1)==0) then
          call gethist4s(ncid,iarchi,'mixr',ik,kk,qg,iotest,
     &                   nface4,xg4,yg4,nord,sigin,2)                   !     mixing ratio
        else
          call gethist4s(ncid,iarchi,'q',ik,kk,qg,iotest,
     &                   nface4,xg4,yg4,nord,sigin,2)                   !     mixing ratio
        end if
      else
        qg=qgmin
      end if ! (nested==0.or.(nested==1.and.nud_q/=0))

      ! re-grid surface pressure by mapping to MSLP, interpolating and then map to surface pressure
      ! requires psl_a, zss, zss_a, t and t_a_lev
      if (nested==0.or.(nested==1.and.nud_test/=0)) then
        if (.not.iotest) then
          if (myid==0) then
            call mslpx(pmsl_a,psl_a,zss_a,t_a_lev,ik*ik*6,sigin(levkk))  ! needs pmsl (preferred)
          end if
          call doints4(pmsl_a,pmsl, nface4,xg4,yg4,nord,ik)
!         invert pmsl to get psl
          call to_pslx(pmsl,psl,zss,t(:,lev),ifull,lev)  ! on target grid
        end if ! .not.iotest
      end if


      !**************************************************************
      ! This is the end of reading the nudging arrays
      !**************************************************************


      !--------------------------------------------------------------
      ! The following data is only read for initial conditions
      if (nested/=1) then

        !--------------------------------------------------
        ! check soil variables
        if (myid==0.or.pfall) then
          ierc(7:7+3*ms)=0
          if (ccycle==0) then
            call ccnf_inq_varid(ncid,'cplant1',idv,tst)
            if (tst) ierc(7)=-1
          else
            call ccnf_inq_varid(ncid,'glai',idv,tst)
            if (tst) ierc(7)=-1
          end if
          do k=1,ms
            write(vname,'("tgg",I1.1)') k
            call ccnf_inq_varid(ncid,vname,idv,tst)
            if (tst) ierc(7+k)=-1
            write(vname,'("wetfrac",I1.1)') k
            call ccnf_inq_varid(ncid,vname,idv,tst)
            if (tst) ierc(7+ms+k)=-1
            write(vname,'("wb",I1.1)') k
            call ccnf_inq_varid(ncid,vname,idv,tst)
            if (tst) ierc(7+2*ms+k)=-1
          end do
        end if
        ! verify if input is a restart file
        if (nested==0) then
          if (myid==0.or.pfall) then
            if (kk==kl.and.iotest) then
              lrestart=.true.
              call ccnf_inq_varid(ncid,'omega',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'zgnhs',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'sdot',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'pslx',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'savu',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'savv',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'savu1',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'savv1',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'savu2',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'savv2',idv,tst)
              if (tst) lrestart=.false.
              call ccnf_inq_varid(ncid,'nstag',idv,tst)
              if (tst) then
                lrestart=.false.
              else 
                call ccnf_get_var1(ncid,idv,iarchi,ierc(3))
              end if
              call ccnf_inq_varid(ncid,'nstagu',idv,tst)
              if (tst) then
                lrestart=.false.
              else 
                call ccnf_get_var1(ncid,idv,iarchi,ierc(4))
              end if
              call ccnf_inq_varid(ncid,'nstagoff',idv,tst)
              if (tst) then
                lrestart=.false.
              else 
                call ccnf_get_var1(ncid,idv,iarchi,ierc(5))
              end if
              if (abs(nmlo)>=3.and.abs(nmlo)<=9) then
                if (ok==wlev) then
                  call ccnf_inq_varid(ncid,'oldu101',idv,tst)
                  if (tst) lrestart=.false.
                  call ccnf_inq_varid(ncid,'oldv101',idv,tst)
                  if (tst) lrestart=.false.
                  call ccnf_inq_varid(ncid,'oldu201',idv,tst)
                  if (tst) lrestart=.false.
                  call ccnf_inq_varid(ncid,'oldv201',idv,tst)
                  if (tst) lrestart=.false.
                  call ccnf_inq_varid(ncid,'ipice',idv,tst)
                  if (tst) lrestart=.false.
                  call ccnf_inq_varid(ncid,'nstagoffmlo',idv,tst)
                  if (tst) then
                    lrestart=.false.
                  else
                    call ccnf_get_var1(ncid,idv,iarchi,ierc(6))
                  end if
                else
                  lrestart=.false.
                end if
              end if
            else
              lrestart=.false.
            end if
            ierc(1:2)=0
            if (lrestart) ierc(1)=1
            call ccnf_inq_varid(ncid,'u10',idv,tst)
            if (tst) ierc(2)=-1
          end if
          if (.not.pfall) then
            call ccmpi_bcast(ierc(1:7+3*ms),0,comm_world)
          end if
          lrestart=(ierc(1)==1)
          if (lrestart) then
            nstag      =ierc(3)
            nstagu     =ierc(4)
            nstagoff   =ierc(5)
            nstagoffmlo=ierc(6)
            if (myid==0) then
              write(6,*) "Continue stagging from"
              write(6,*) "nstag,nstagu,nstagoff ",
     &                    nstag,nstagu,nstagoff
              if (abs(nmlo)>=3.and.abs(nmlo)<=9) then
                write(6,*) "nstagoffmlo ",nstagoffmlo
              end if
            end if
          end if
        else
          if (.not.pfall) then
            call ccmpi_bcast(ierc(7:7+3*ms),0,comm_world)
          end if            
        end if ! nested==0 ..else..
        
        !--------------------------------------------------

        ! SNOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call filhist1(ncid,iarchi,'snd',ik,dk,snowd,iotest,
     &                nface4,xg4,yg4,nord,sea_a)

        ! SOIL TEMPERATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,ms 
          if (ierc(7+k)==0) then
            write(vname,'("tgg",I1.1)') k
          else if (k<=3.and.ierc(7+2)==0) then
            vname="tgg2"
          else if (k<=3) then
            vname="tb3"
          else if (ierc(7+6)==0) then
            vname="tgg6"
          else
            vname="tb2"
          end if
          if (iotest) then
            if (k==1.and.ierc(7+1)/=0) then
              tgg(:,1)=tss
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     tgg(:,k),ifull)
            end if
          else
            if (k==1.and.ierc(7+1)/=0) then
              ucc(1:dk*dk*6)=tss_a(1:dk*dk*6)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
            end if
            call fill_cc(ucc,dk,0,sea_a)
            call doints4(ucc,tgg(:,k),nface4,xg4,yg4,
     &                     nord,ik)
          end if
        end do
        if (.not.iotest) then
          where (snowd>0.)
            tgg(:,1)=min(tgg(:,1),270.1)
          endwhere
        end if

        !--------------------------------------------------
        ! Read MLO sea-ice data
        if (nmlo/=0.and.abs(nmlo)<=9) then
          if (.not.allocated(micdwn)) allocate(micdwn(ifull,11))
          do k=1,7
            select case(k)
              case(1,2,3,4)
                write(vname,'("tggsn",I1.1)') k
                call filhist1(ncid,iarchi,vname,ik,dk,micdwn(:,k),
     &                        iotest,nface4,xg4,yg4,nord,land_a)
                if (all(micdwn(:,k)==0.)) micdwn(:,k)=280.
              case(5)
                micdwn(:,k)=fracice ! read above with nudging arrays
              case(6)
                micdwn(:,k)=sicedep ! read above with nudging arrays
              case(7)
                micdwn(:,k)=snowd*1.E-3
            end select
          end do
          call filhist1(ncid,iarchi,'sto',ik,dk,micdwn(:,8),iotest,
     &                  nface4,xg4,yg4,nord,land_a)
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'uic',ik,6*ik,micdwn(:,9),
     &                   ifull)
            call histrd1(ncid,iarchi,ier,'vic',ik,6*ik,micdwn(:,10),
     &                   ifull)
          else
            call histrd1(ncid,iarchi,ier,'uic',ik,6*ik,ucc,6*ik*ik)
            call histrd1(ncid,iarchi,ier,'vic',ik,6*ik,vcc,6*ik*ik)
            call fill_cc(ucc,dk,0,land_a)
            call fill_cc(vcc,dk,0,land_a)
            call interpwind(ik,micdwn(:,9),micdwn(:,10),ucc,vcc,axs_a,
     &                      ays_a,azs_a,bxs_a,bys_a,bzs_a,rotpole,
     &                      rotpoles,nface4,xg4,yg4,nord,dk)
          end if ! iotest
          call filhist1(ncid,iarchi,'icesal',ik,dk,micdwn(:,11),
     &                  iotest,nface4,xg4,yg4,nord,land_a)
          if (abs(nmlo)>=2) then
            call gethist1(ncid,iarchi,'swater',ik,watbdy,
     &                    iotest,nface4,xg4,yg4,nord)
            call gethist1(ncid,iarchi,'ssalin',ik,salbdy,
     &                    iotest,nface4,xg4,yg4,nord)
          end if
        end if
        !--------------------------------------------------

        ! SOIL MOISTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        wb=20.5
        do k=1,ms
          if (ierc(7+ms+k)==0) then
            write(vname,'("wetfrac",I1.1)') k
          else if (ierc(7+2*ms+k)==0) then
            write(vname,'("wb",I1.1)') k
          else if (k<2.and.ierc(7+2*ms+2)==0) then
            vname="wb2"
          else if (k<2) then
            vname="wfg"
          else if (ierc(7+2*ms+6)==0) then
            vname="wb6"
          else
            vname="wfb"
          end if
          if (iotest) then
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,wb(:,k),ifull)
            if (ierc(7+ms+k)==0) then
              wb(:,k)=wb(:,k)+20. ! flag for fraction of field capacity
            end if
          else
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,ucc,6*ik*ik)
            if (ierc(7+ms+k)==0) then
              ucc=ucc+20.
            end if
            call fill_cc(ucc,dk,0,sea_a)
            call doints4(ucc,wb(:,k),nface4,xg4,yg4,nord,ik)
          end if ! iotest
        end do
        !unpack field capacity into volumetric soil moisture
        if (any(wb(:,:)>10.)) then
          if (mydiag) write(6,*) "Unpacking wetfrac to wb",wb(idjd,1)
          wb(:,:)=wb(:,:)-20.
          do iq=1,ifull
            isoil=isoilm(iq)
            wb(iq,:)=(1.-wb(iq,:))*swilt(isoil)+wb(iq,:)*sfc(isoil)
          end do
          if (mydiag) write(6,*) "giving wb",wb(idjd,1)
        end if

        call filhist1(ncid,iarchi,'wetfac',ik,dk,wetfac,
     &                iotest,nface4,xg4,yg4,nord,sea_a)
        where (.not.land)
          wetfac=1.
        end where

        !--------------------------------------------------
        ! Read 10m wind speeds for special sea roughness length calculations
        if (nested==0) then
          if (ierc(2)==0) then
           call gethist1(ncid,iarchi,'u10',ik,u10,iotest,
     &                   nface4,xg4,yg4,nord)
          else
           u10=sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)*log(10./0.001)
     &                                              /log(zmin/0.001)
          end if
        end if

        !--------------------------------------------------
        ! Read boundary layer height for TKE-eps mixing
        call gethist1(ncid,iarchi,'pblh',ik,pblh,iotest,
     &                nface4,xg4,yg4,nord)
        pblh=max(pblh,1.)
        if (nvmix==6.and.nested==0) then
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'dpblh',ik,6*ik,zidry,ifull)
          else
            zidry=pblh 
          end if ! iotest
          zidry=max(zidry,1.)
        end if

        !--------------------------------------------------
        ! Read CABLE/CASA aggregate C+N+P pools
        if (nsib>=6) then
         if (ccycle==0) then
          if (ierc(7)==0) then
           do k=1,ncp
            write(vname,'("cplant",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,cplant(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
           end do
           do k=1,ncs
            write(vname,'("csoil",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,csoil(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
           end do
          end if
         else
          if (ierc(7)==0) then
           do k=1,mplant
            write(vname,'("cplant",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,cplant(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
            write(vname,'("nplant",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,niplant(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
            write(vname,'("pplant",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,pplant(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
           end do
           do k=1,mlitter
            write(vname,'("clitter",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,clitter(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
            write(vname,'("nlitter",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,nilitter(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
            write(vname,'("plitter",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,plitter(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
           end do         
           do k=1,msoil
            write(vname,'("csoil",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,csoil(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
            write(vname,'("nsoil",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,nisoil(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
            write(vname,'("psoil",I1.1)') k
            call filhist1(ncid,iarchi,vname,ik,dk,psoil(:,k),
     &                    iotest,nface4,xg4,yg4,nord,sea_a)
           end do
           call filhist1(ncid,iarchi,'glai',ik,dk,glai,
     &                   iotest,nface4,xg4,yg4,nord,sea_a)
          end if ! ierc(7)==0
         end if ! ccycle==0 ..else..
        end if ! if nsib==6.or.nsib==7

        !--------------------------------------------------
        ! Read urban data
        if (nurban/=0) then
          if (.not.allocated(atebdwn)) allocate(atebdwn(ifull,24))
          do k=1,24
            select case(k)
              case(1)
                vname='rooftgg1'
              case(2)
                vname='rooftgg2'
              case(3)
                vname='rooftgg3'
              case(4)
                vname='waletgg1'
              case(5)
                vname='waletgg2'
              case(6)
                vname='waletgg3'
              case(7)
                vname='walwtgg1'
              case(8)
                vname='walwtgg2'
              case(9)
                vname='walwtgg3'
              case(10)
                vname='roadtgg1'
              case(11)
                vname='roadtgg2'
              case(12)
                vname='roadtgg3'
              case(13)
                vname='urbnsmc'
              case(14)
                vname='urbnsmr'
              case(15)
                vname='roofwtr'
              case(16)
                vname='roadwtr'
              case(17)
                vname='urbwtrc'
              case(18)
                vname='urbwtrr'
              case(19)
                vname='roofsnd'
              case(20)
                vname='roadsnd'
              case(21)
                vname='roofden'
              case(22)
                vname='roadden'
              case(23)
                vname='roofsna'
              case(24)
                vname='roadsna'
            end select
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     atebdwn(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              where (ucc>=399.)
                ucc=999.
              end where
              call fill_cc(ucc,dk,0,sea_a)
              call doints4(ucc,atebdwn(:,k),nface4,xg4,yg4,nord,ik)
            end if ! iotest
            if (all(atebdwn(:,k)==0.)) then
              select case(k)
                case(1:12)
                  atebdwn(:,k)=300.
                case(21:22)
                  atebdwn(:,k)=100.
                case(23:24)
                  atebdwn(:,k)=0.85
              end select
            end if
          end do
        end if
        !--------------------------------------------------
        
        if (nested==0) then
          ! OMEGA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          dpsldt=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'omega',ik,6*ik,k,
     &                     dpsldt(:,k),ifull)
             dpsldt(:,k)=dpsldt(:,k)/(1.e5*exp(psl))
            enddo  ! k loop
          end if

          ! CLOUD FROZEN WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call gethist4s(ncid,iarchi,'qfg',ik,kk,dum,
     &                   iotest,nface4,xg4,yg4,nord,sigin,5)
          qfg(1:ifull,:)=dum
          ! CLOUD LIQUID WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call gethist4s(ncid,iarchi,'qlg',ik,kk,dum,
     &                   iotest,nface4,xg4,yg4,nord,sigin,5)
          qlg(1:ifull,:)=dum
          ! RAIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call gethist4s(ncid,iarchi,'qrg',ik,kk,dum,
     &                   iotest,nface4,xg4,yg4,nord,sigin,5)
          qrg(1:ifull,:)=dum
          ! CLOUD FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call gethist4s(ncid,iarchi,'cfrac',ik,kk,cfrac,
     &                   iotest,nface4,xg4,yg4,nord,sigin,5)
          ! RAIN FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call gethist4s(ncid,iarchi,'cfrain',ik,kk,dum,
     &                   iotest,nface4,xg4,yg4,nord,sigin,5)
          cffall(1:ifull,:)=dum
          if (ncloud>=3) then
            ! STRAT CLOUD FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call gethist4s(ncid,iarchi,'stratcf',ik,kk,dum,
     &                     iotest,nface4,xg4,yg4,nord,sigin,5)
            stratcloud(1:ifull,:)=dum
            ! STRAT NET TENDENCY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call gethist4s(ncid,iarchi,'strat_nt',ik,kk,nettend,
     &                     iotest,nface4,xg4,yg4,nord,sigin,5)
            ! STRAT NET MASS FLUX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !call gethist4s(ncid,iarchi,'strat_mf',ik,kk,cmflx,
     &      !               iotest,nface4,xg4,yg4,nord,sigin,5)
          end if ! (ncloud>=3)
        end if   ! (nested==0)

        !--------------------------------------------------
        ! TKE-eps data
        if (nvmix==6.and.nested==0) then
          call gethist4s(ncid,iarchi,'tke',ik,kk,dum,iotest,
     &                   nface4,xg4,yg4,nord,sigin,5)
          if (all(dum==0.)) dum=1.5E-4
          tke(1:ifull,:)=dum
          call gethist4s(ncid,iarchi,'eps',ik,kk,dum,iotest,
     &                   nface4,xg4,yg4,nord,sigin,5)
          if (all(dum==0.)) dum=1.E-7
          eps(1:ifull,:)=dum
        end if

        !------------------------------------------------------------
        ! Tracer data
        if (ngas>0) then              
          do igas=1,ngas              
            write(trnum,'(i3.3)') igas
            call gethist4s(ncid,iarchi,'tr'//trnum,ik,kk,dum,
     &                     iotest,nface4,xg4,yg4,nord,sigin,7)
            tr(1:ifull,:,igas)=dum
          enddo                       
        endif                         

        !------------------------------------------------------------
        ! Aerosol data
        if (abs(iaero)>=2) then
          do i=1,naero+2
            select case(i)
              case(1)
                vname='dms'
              case(2)
                vname='so2'
              case(3)
                vname='so4'
              case(4)
                vname='bco'
              case(5)
                vname='bci'
              case(6)
                vname='oco'
              case(7)
                vname='oci'
              case(8)
                vname='dust1'
              case(9)
                vname='dust2'
              case(10)
                vname='dust3'
              case(11)
                vname='dust4'
              case(12)
                vname='seasalt1'
              case(13)
                vname='seasalt2'
              case default
                write(6,*) "ERROR: Unknown aerosol type ",i
                stop
            end select
            call gethist4s(ncid,iarchi,vname,ik,kk,dum,iotest,
     &                     nface4,xg4,yg4,nord,sigin,5)
            if (i<=naero) then
              xtg(1:ifull,:,i)=dum
            else
              ssn(1:ifull,:,i-naero)=dum
            end if
          end do
          if (iaero/=0) then
            ! Factor 1.e3 to convert to g/m2, x 3 to get sulfate from sulfur
            so4t(:)=0.
            do k=1,kl
              so4t(:)=so4t(:)+3.e3*xtg(1:ifull,k,3)*
     &          (-1.e5*exp(psl(1:ifull))*dsig(k))/grav
            enddo
          end if
        end if

        if (nested==0) then
          ! GEOPOTENTIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          phi_nh=0.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'zgnhs',ik,6*ik,k,
     &                     phi_nh(:,k),ifull)
            enddo  ! k loop
          end if

          ! SDOT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          sdot=-999.
          if (lrestart) then
            sdot(:,1)=0.
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'sdot',ik,6*ik,k,
     &                     sdot(:,k+1),ifull)
            enddo  ! k loop
          end if

          ! PSLX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          pslx=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'pslx',ik,6*ik,k,
     &                     pslx(1:ifull,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savu=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savu',ik,6*ik,k,
     &                     savu(:,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savv=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savv',ik,6*ik,k,
     &                     savv(:,k),ifull)
            enddo  ! k loop
          end if

          ! SAVU1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savu1=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savu1',ik,6*ik,k,
     &                     savu1(:,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVV1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savv1=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savv1',ik,6*ik,k,
     &                     savv1(:,k),ifull)
            enddo  ! k loop
          end if

          ! SAVU2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savu2=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savu2',ik,6*ik,k,
     &                     savu2(:,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVV2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savv2=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savv2',ik,6*ik,k,
     &                     savv2(:,k),ifull)
            enddo  ! k loop
          end if
          
          if (abs(nmlo)>=3.and.abs(nmlo)<=9) then
           oldu1=0.
           oldu2=0.
           oldv1=0.
           oldv2=0.
           ipice=0.
           if (lrestart) then
            do k=1,ok
              write(vname,'("oldu1",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     oldu1(:,k),ifull)
              write(vname,'("oldv1",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     oldv1(:,k),ifull)
              write(vname,'("oldu2",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     oldu2(:,k),ifull)
              write(vname,'("oldv2",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     oldv2(:,k),ifull)
            end do
            call histrd1(ncid,iarchi,ier,'ipice',ik,6*ik,
     &                   ipice(1:ifull),ifull)
           end if
          end if
       
        end if ! (nested==0)

        ! SOIL ICE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,ms
          write(vname,'("wbice",I1.1)') k
          call filhist1(ncid,iarchi,vname,ik,dk,wbice(:,k),
     &                  iotest,nface4,xg4,yg4,nord,sea_a)
        end do

        ! SNOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (nmlo==0.or.abs(nmlo)>9) then ! otherwise already read above
         do k=1,3
          write(vname,'("tggsn",I1.1)') k
          call filhist1(ncid,iarchi,vname,ik,dk,tggsn(:,k),
     &                  iotest,nface4,xg4,yg4,nord,sea_a)
          if (all(tggsn(:,k)==0.)) tggsn(:,k)=280.
          where(.not.land)
            tggsn(:,k)=280.
          end where
         end do
        end if
        do k=1,3
          write(vname,'("smass",I1.1)') k
          call filhist1(ncid,iarchi,vname,ik,dk,smass(:,k),
     &                  iotest,nface4,xg4,yg4,nord,sea_a)
        end do
        do k=1,3
          write(vname,'("ssdn",I1.1)') k
          call filhist1(ncid,iarchi,vname,ik,dk,ssdn(:,k),
     &                  iotest,nface4,xg4,yg4,nord,sea_a)
          if (all(ssdn(:,k)==0.)) then
            where (snowd>100.)
              ssdn(:,k)=240.
            elsewhere
              ssdn(:,k)=140.
            end where
          end if
        end do
        ssdnn=ssdn(:,1)
        call filhist1(ncid,iarchi,'snage',ik,dk,snage,
     &                iotest,nface4,xg4,yg4,nord,sea_a)

        call filhist1(ncid,iarchi,'sflag',ik,dk,dum6,
     &                iotest,nface4,xg4,yg4,nord,sea_a)
        isflag=nint(dum6)

        call gethist1(ncid,iarchi,'sgsave',ik,sgsave,
     &                iotest,nface4,xg4,yg4,nord)
        
       endif    ! (nested/=1)

      !**************************************************************
      ! This is the end of reading the initial arrays
      !**************************************************************         
          
      ! tgg holds file surface temperature when no MLO
      if (nmlo==0.or.abs(nmlo)>9) then
        where (.not.land)
          tgg(:,1)=tss
        end where
      end if

      ! set-up for next read of file
      iarchi=iarchi+1
      kdate_s=kdate_r
      ktime_s=ktime_r+1

      if (myid==0.and.nested==0) then
        write(6,*) "Final lrestart ",lrestart
      end if

      return
      end subroutine ontheflyx

      subroutine doints4(s,sout,nface4,xg4,yg4,nord,ik)  ! does calls to intsb
      
      use cc_mpi                 ! CC MPI routines
      
      implicit none
      
      include 'newmpar.h'        ! Grid parameters
      include 'parm.h'           ! Model configuration
      
      integer, parameter :: ntest=0
      integer, intent(in) :: ik, nord
      integer :: iq, mm
      integer :: idx = 25, jdx = 218, idjdx=10441
      integer, intent(in), dimension(ifull,4) :: nface4
      real, dimension(ik*ik*6), intent(inout) :: s
      real, dimension(ifull), intent(inout) :: sout
      real, intent(in), dimension(ifull,4) :: xg4, yg4
      real wrk(ifull,4)

      call ccmpi_bcast(s,0,comm_world)

      if(nord==1)then
         do mm=1,m_fly  !  was 4, now may be 1
            call ints_blb(s,wrk(1,mm),nface4(1,mm),xg4(1,mm),
     &                    yg4(1,mm),ik)
         enddo
      else
         do mm=1,m_fly  !  was 4, now may be 1
            call intsb(s,wrk(1,mm),nface4(1,mm),xg4(1,mm),
     &                 yg4(1,mm),ik)
         enddo
      endif   ! (nord==1)  .. else ..
      if(ntest>0.and.mydiag)then
         write(6,*)'in ints4 for id,jd,nord: ',id,jd,nord
         write(6,*)'nface4(1-4) ',(nface4(idjd,mm),mm=1,4)
         write(6,*)'xg4(1-4) ',(xg4(idjd,mm),mm=1,4)
         write(6,*)'yg4(1-4) ',(yg4(idjd,mm),mm=1,4)
         write(6,*)'wrk(1-4) ',(wrk(idjd,mm),mm=1,4)
      endif
      if(m_fly==1)then
        do iq=1,ifull
         sout(iq)=wrk(iq,1)
        enddo
      else
!       average 4 m values to central value
        do iq=1,ifull
         sout(iq)=.25*(wrk(iq,1)+wrk(iq,2)+wrk(iq,3)+wrk(iq,4))
        enddo
      endif    ! (m_fly==1)

      end subroutine doints4

      subroutine intsb(s,sout,nface,xg,yg,ik)   ! N.B. sout here
      
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
c     later may wish to save idel etc between array calls
c     this one does linear interp in x on outer y sides
c     doing x-interpolation before y-interpolation
!     This is a global routine 

      use cc_mpi, only : indx

      implicit none
      
      include 'newmpar.h'  ! Grid parameters
      include 'parm.h'     ! Model configuration

      integer, dimension(ifull), intent(in) :: nface
      integer :: idel, jdel, ik, nn
      integer :: i, j, n, iq, n_n, n_e, n_w, n_s
      real, dimension(ik*ik*6), intent(in) :: s
      real, dimension(ifull), intent(inout) :: sout
      real aaa, c1, c2, c3, c4, xxg, yyg
      real, intent(in), dimension(ifull) :: xg, yg
      real sx(-1:ik+2,-1:ik+2,0:npanels)
      real r(4)

c     this is intsb           EW interps done first
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
      do n=0,npanels
       do j=1,ik
        do i=1,ik
         sx(i,j,n)=s(indx(i,j,n,ik,ik))
        enddo  ! i loop
       enddo   ! j loop
      enddo    ! n loop
      do n=0,npanels,2
       n_w=mod(n+5,6)
       n_e=mod(n+2,6)
       n_n=mod(n+1,6)
       n_s=mod(n+4,6)
       do j=1,ik
        sx(0,j,n)=s(indx(ik,j,n_w,ik,ik))
        sx(-1,j,n)=s(indx(ik-1,j,n_w,ik,ik))
        sx(ik+1,j,n)=s(indx(ik+1-j,1,n_e,ik,ik))
        sx(ik+2,j,n)=s(indx(ik+1-j,2,n_e,ik,ik))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(indx(i,1,n+1,ik,ik))
        sx(i,ik+2,n)=s(indx(i,2,n+1,ik,ik))
        sx(i,0,n)=s(indx(ik,ik+1-i,n_s,ik,ik))
        sx(i,-1,n)=s(indx(ik-1,ik+1-i,n_s,ik,ik))
       enddo
       sx(-1,0,n)=s(indx(ik,2,n_w,ik,ik))       ! wws
       sx(0,-1,n)=s(indx(ik,ik-1,n_s,ik,ik))    ! wss
       sx(0,0,n)=s(indx(ik,1,n_w,ik,ik))        ! ws
       sx(ik+1,0,n)=s(indx(ik,1,n_e,ik,ik))     ! es  
       sx(ik+2,0,n)=s(indx(ik-1,1,n_e,ik,ik))   ! ees 
       sx(-1,ik+1,n)=s(indx(ik,ik-1,n_w,ik,ik)) ! wwn
       sx(0,ik+2,n)=s(indx(ik-1,ik,n_w,ik,ik))  ! wnn
       sx(ik+2,ik+1,n)=s(indx(2,1,n_e,ik,ik))   ! een  
       sx(ik+1,ik+2,n)=s(indx(1,2,n_e,ik,ik))   ! enn  
       sx(0,ik+1,n)   =s(indx(ik,ik,n_w,ik,ik)) ! wn  
       sx(ik+1,ik+1,n)=s(indx(1,1,n_e,ik,ik))   ! en  
       sx(ik+1,-1,n)=s(indx(ik,2,n_e,ik,ik))    ! ess  
      enddo  ! n loop
      do n=1,npanels,2
       n_w=mod(n+4,6)
       n_e=mod(n+1,6)
       n_n=mod(n+2,6)
       n_s=mod(n+5,6)
       do j=1,ik
        sx(0,j,n)=s(indx(ik+1-j,ik,n_w,ik,ik))
        sx(-1,j,n)=s(indx(ik+1-j,ik-1,n_w,ik,ik))
        sx(ik+1,j,n)=s(indx(1,j,n_e,ik,ik))
        sx(ik+2,j,n)=s(indx(2,j,n_e,ik,ik))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(indx(1,ik+1-i,n_n,ik,ik))
        sx(i,ik+2,n)=s(indx(2,ik+1-i,n_n,ik,ik))
        sx(i,0,n)=s(indx(i,ik,n-1,ik,ik))
        sx(i,-1,n)=s(indx(i,ik-1,n-1,ik,ik))
       enddo
       sx(-1,0,n)=s(indx(ik-1,ik,n_w,ik,ik))     ! wws
       sx(0,-1,n)=s(indx(2,ik,n_s,ik,ik))        ! wss
       sx(0,0,n)=s(indx(ik,ik,n_w,ik,ik))        ! ws
       sx(ik+1,0,n)=s(indx(1,1,n_e,ik,ik))       ! es
       sx(ik+2,0,n)=s(indx(1,2,n_e,ik,ik))       ! ees
       sx(-1,ik+1,n)=s(indx(2,ik,n_w,ik,ik))     ! wwn   
       sx(0,ik+2,n)=s(indx(1,ik-1,n_w,ik,ik))    ! wnn  
       sx(ik+2,ik+1,n)=s(indx(1,ik-1,n_e,ik,ik)) ! een  
       sx(ik+1,ik+2,n)=s(indx(2,ik,n_e,ik,ik))   ! enn  
       sx(0,ik+1,n)   =s(indx(1,ik,n_w,ik,ik))   ! wn  
       sx(ik+1,ik+1,n)=s(indx(1,ik,n_e,ik,ik))   ! en  
       sx(ik+1,-1,n)=s(indx(2,1,n_e,ik,ik))      ! ess  
      enddo  ! n loop
c     for ew interpolation, sometimes need (different from ns):
c          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

      do iq=1,ifull   ! runs through list of target points
        n=nface(iq)
        idel=int(xg(iq))
        xxg=xg(iq)-idel
c       yg here goes from .5 to il +.5
        jdel=int(yg(iq))
        yyg=yg(iq)-jdel
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
        aaa=min( aaa , max( sx(idel,jdel,n),sx(idel+1,jdel,n),
     .                      sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
        sout(iq)=max( aaa , min( sx(idel,jdel,n),sx(idel+1,jdel,n),
     .                        sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
      enddo    ! iq loop

      end subroutine intsb

      subroutine ints_blb(s,sout,nface,xg,yg,ik) 
      
c     this one does bi-linear interpolation only

      use cc_mpi, only : indx

      implicit none
      
      include 'newmpar.h'  ! Grid parameters
      include 'parm.h'     ! Model configuration

      integer :: i, j, n, iq, idel, jdel, ik
      integer :: n_n, n_e, n_w, n_s
      integer, intent(in), dimension(ifull) :: nface
      real, dimension(ik*ik*6), intent(inout) :: s
      real, dimension(ifull), intent(inout) :: sout
      real, intent(in), dimension(ifull) :: xg, yg
      real sx(-1:ik+2,-1:ik+2,0:npanels)
      real :: xxg, yyg

c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
c                    but for bi-linear only need 0:il+1 &  0:il+1
      do n=0,npanels
       do j=1,ik
        do i=1,ik
         sx(i,j,n)=s(indx(i,j,n,ik,ik))
        enddo  ! i loop
       enddo   ! j loop
      enddo    ! n loop
      do n=0,npanels,2
       n_w=mod(n+5,6)
       n_e=mod(n+2,6)
       n_n=mod(n+1,6)
       n_s=mod(n+4,6)
       do j=1,ik
        sx(0,j,n)=s(indx(ik,j,n_w,ik,ik))
        sx(ik+1,j,n)=s(indx(ik+1-j,1,n_e,ik,ik))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(indx(i,1,n+1,ik,ik))
        sx(i,0,n)=s(indx(ik,ik+1-i,n_s,ik,ik))
       enddo
       sx(0,0,n)=s(indx(ik,1,n_w,ik,ik))        ! ws
       sx(ik+1,0,n)=s(indx(ik,1,n_e,ik,ik))     ! es  
       sx(0,ik+1,n)   =s(indx(ik,ik,n_w,ik,ik)) ! wn  
       sx(ik+1,ik+1,n)=s(indx(1,1,n_e,ik,ik))   ! en  
      enddo  ! n loop
      do n=1,npanels,2
       n_w=mod(n+4,6)
       n_e=mod(n+1,6)
       n_n=mod(n+2,6)
       n_s=mod(n+5,6)
       do j=1,ik
        sx(0,j,n)=s(indx(ik+1-j,ik,n_w,ik,ik))
        sx(ik+1,j,n)=s(indx(1,j,n_e,ik,ik))
       enddo
       do i=1,ik
        sx(i,ik+1,n)=s(indx(1,ik+1-i,n_n,ik,ik))
        sx(i,0,n)=s(indx(i,ik,n-1,ik,ik))
       enddo
       sx(0,0,n)=s(indx(ik,ik,n_w,ik,ik))        ! ws
       sx(ik+1,0,n)=s(indx(1,1,n_e,ik,ik))       ! es
       sx(0,ik+1,n)   =s(indx(1,ik,n_w,ik,ik))   ! wn  
       sx(ik+1,ik+1,n)=s(indx(1,ik,n_e,ik,ik))   ! en  
      enddo  ! n loop
      

      do iq=1,ifull  ! runs through list of target points
       n=nface(iq)
       idel=int(xg(iq))
       xxg=xg(iq)-idel
       jdel=int(yg(iq))
       yyg=yg(iq)-jdel
       sout(iq)=yyg*(xxg*sx(idel+1,jdel+1,n)
     .               +(1.-xxg)*sx(idel,jdel+1,n))
     .    +(1.-yyg)*(xxg*sx(idel+1,jdel,n)
     .               +(1.-xxg)*sx(idel,jdel,n))
      enddo    ! iq loop

      end subroutine ints_blb

      subroutine fill_cc(a_io,ik,ndiag,land_a)
      
!     this version holds whole array in memory      
c     routine fills in interior of an array which has undefined points

      use cc_mpi          ! CC MPI routines
      
      implicit none
      
      include 'newmpar.h' ! Grid parameters

      integer :: nrem, i, ii, ik, iq, j, n, neighb, ndiag
      integer :: iminb, imaxb, jminb, jmaxb
      integer, save :: oldik = 0
      integer, dimension(:,:), allocatable, save :: ic
      integer, dimension(0:5) :: imin,imax,jmin,jmax
      integer npann(0:5),npane(0:5),npanw(0:5),npans(0:5)
      real, parameter :: value=999.       ! missing value flag
      real a_io(ik*ik*6)         ! input and output array
      real a(ik*ik*6)
      real :: av     
      logical land_a(ik*ik*6)
      logical, dimension(4) :: mask
      data npann/1,103,3,105,5,101/,npane/102,2,104,4,100,0/
      data npanw/5,105,1,101,3,103/,npans/104,0,100,2,102,4/

      if (myid/=0) return

      where (land_a)
        a_io=value
      end where
      if (all(abs(a_io-value)<1.E-6)) return

      START_LOG(otf_fill)
      
      if (ik/=oldik) then
       oldik=ik
       if (allocated(ic)) then
         deallocate(ic)
       end if
       allocate(ic(4,ik*ik*6))
       do iq=1,ik*ik*6
        ic(1,iq)=iq+ik
        ic(2,iq)=iq-ik
        ic(3,iq)=iq+1
        ic(4,iq)=iq-1
       enddo   ! iq loop
       do n=0,npanels
        if(npann(n)<100)then
         do ii=1,ik
          ic(1,indx(ii,ik,n,ik,ik))=indx(ii,1,npann(n),ik,ik)
         enddo    ! ii loop
        else
         do ii=1,ik
          ic(1,indx(ii,ik,n,ik,ik))=indx(1,ik+1-ii,npann(n)-100,ik,ik)
         enddo    ! ii loop
        endif      ! (npann(n)<100)
        if(npane(n)<100)then
         do ii=1,ik
          ic(3,indx(ik,ii,n,ik,ik))=indx(1,ii,npane(n),ik,ik)
         enddo    ! ii loop
        else
         do ii=1,ik
          ic(3,indx(ik,ii,n,ik,ik))=indx(ik+1-ii,1,npane(n)-100,ik,ik)
         enddo    ! ii loop
        endif      ! (npane(n)<100)
        if(npanw(n)<100)then
         do ii=1,ik
          ic(4,indx(1,ii,n,ik,ik))=indx(ik,ii,npanw(n),ik,ik)
         enddo    ! ii loop
        else
         do ii=1,ik
          ic(4,indx(1,ii,n,ik,ik))=indx(ik+1-ii,ik,npanw(n)-100,ik,ik)
         enddo    ! ii loop
        endif      ! (npanw(n)<100)
        if(npans(n)<100)then
         do ii=1,ik
          ic(2,indx(ii,1,n,ik,ik))=indx(ii,ik,npans(n),ik,ik)
         enddo    ! ii loop
        else
         do ii=1,ik
          ic(2,indx(ii,1,n,ik,ik))=indx(ik,ik+1-ii,npans(n)-100,ik,ik)
         enddo    ! ii loop
        endif      ! (npans(n)<100)
       enddo      ! n loop
      end if ! oldik/=ik

      imin=1
      imax=ik
      jmin=1
      jmax=ik
          
      nrem = 1    ! Just for first iteration
      do while ( nrem > 0)
       nrem=0
       do iq=1,ik*ik*6
        a(iq)=a_io(iq)
       enddo
       ! MJT restricted fill
       do n=0,5
        iminb=ik
        imaxb=1
        jminb=ik
        jmaxb=1
        do j=jmin(n),jmax(n)
         do i=imin(n),imax(n)
          iq=indx(i,j,n,ik,ik)
          if(a(iq)==value)then
           mask=a(ic(:,iq))/=value
           neighb=count(mask)
           if(neighb>0)then
            av=sum(a(ic(:,iq)),mask)
            a_io(iq)=av/real(neighb)
           else
            iminb=min(i,iminb)
            imaxb=max(i,imaxb)
            jminb=min(j,jminb)
            jmaxb=max(j,jmaxb)
            nrem=nrem+1   ! current number of points without a neighbour
           endif
          endif
         end do
        end do
        imin(n)=iminb
        imax(n)=imaxb
        jmin(n)=jminb
        jmax(n)=jmaxb
       end do
      end do
      
      END_LOG(otf_fill)
      
      return
      end subroutine fill_cc

      subroutine mslpx(pmsl,psl,zs,t,ifullx,siglev)
      
!     generalized from ifull to allow new onthefly    0808
!     can replace usual mslp sometime

      use cc_mpi, only : mydiag ! CC MPI routines
      use sigs_m                ! Atmosphere sigma levels
      
!     this one will ignore negative zs (i.e. over the ocean)
      implicit none
      
      include 'newmpar.h'       ! Grid parameters
      include 'const_phys.h'    ! Physical constants
      include 'parm.h'          ! Model configuration

      integer ifullx,iq
      real siglev
      real, dimension(ifullx) :: pmsl,psl,zs,t
      real c, con, conr
      real, dimension(ifullx) :: dlnps, phi1, tav, tsurf
      c=grav/stdlapse
      conr=c/rdry
      con=siglev**(rdry/c)/c
      phi1(:)=t(:)*rdry*(1.-siglev)/siglev ! phi of sig(lev) above sfce
      tsurf(:)=t(:)+phi1(:)*stdlapse/grav
      tav(:)=tsurf(:)+max(0.,zs(:))*.5*stdlapse/grav
      dlnps(:)=max(0.,zs(:))/(rdry*tav(:))
      pmsl(:)=1.e5*exp(psl(:)+dlnps(:))
      return
      end subroutine mslpx
      
      subroutine to_pslx(pmsl,psl,zs,t,ifullx,lev)
!     generalized from ifull to allow new onthefly    0808
!     can replace usual mslp sometime

      use cc_mpi, only : mydiag  ! CC MPI routines
      use sigs_m                 ! Atmosphere sigma levels
      
!     this one will ignore negative zs (i.e. over the ocean)
      implicit none
      
      include 'newmpar.h'        ! Grid parameters
      include 'const_phys.h'     ! Physical constants
      include 'parm.h'           ! Model configuration
      
      integer ifullx,iq,lev
      real, dimension(ifullx) :: pmsl,psl,zs,t
      real, dimension(ifullx) :: dlnps, phi1, tav, tsurf
      phi1(:)=t(:)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
      tsurf(:)=t(:)+phi1(:)*stdlapse/grav
      tav(:)=tsurf(:)+max(0.,zs(:))*.5*stdlapse/grav
      dlnps(:)=max(0.,zs(:))/(rdry*tav(:))
      psl(:)=log(1.e-5*pmsl(:)) -dlnps(:)
      if(nmaxpr==1.and.mydiag)then
        write(6,*)'to_psl lev,sig(lev) ',lev,sig(lev)
        write(6,*)'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd),psl(idjd),pmsl(idjd)
      endif
      return
      end subroutine to_pslx
      
      subroutine interpwind(ik,uct,vct,ucc,vcc,axs_a,ays_a,azs_a,
     &                      bxs_a,bys_a,bzs_a,rotpole,rotpoles,nface4,
     &                      xg4,yg4,nord,dk)
      
      use cc_mpi           ! CC MPI routines
      use vecsuv_m         ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'  ! Grid parameters
      
      integer, intent(in) :: ik,nord,dk
      integer, dimension(ifull), intent(in) :: nface4
      integer iq
      real, dimension(3,3), intent(in) :: rotpole,rotpoles
      real, dimension(ifull), intent(in) :: xg4,yg4
      real, dimension(6*dk*dk), intent(in) :: axs_a,ays_a,azs_a
      real, dimension(6*dk*dk), intent(in) :: bxs_a,bys_a,bzs_a
      real, dimension(6*dk*dk) :: uc,vc,wc
      real, dimension(6*ik*ik), intent(inout) :: ucc,vcc
      real, dimension(6*ik*ik) :: wcc
      real, dimension(ifull), intent(out) :: uct,vct
      real, dimension(ifull) :: wct
      
      do iq=1,dk*dk*6
        ! first set up winds in Cartesian "source" coords            
        uc(iq)=axs_a(iq)*ucc(iq) + bxs_a(iq)*vcc(iq)
        vc(iq)=ays_a(iq)*ucc(iq) + bys_a(iq)*vcc(iq)
        wc(iq)=azs_a(iq)*ucc(iq) + bzs_a(iq)*vcc(iq)
        ! now convert to winds in "absolute" Cartesian components
        ucc(iq)=uc(iq)*rotpoles(1,1)+vc(iq)*rotpoles(1,2)
     &         +wc(iq)*rotpoles(1,3)
        vcc(iq)=uc(iq)*rotpoles(2,1)+vc(iq)*rotpoles(2,2)
     &         +wc(iq)*rotpoles(2,3)
        wcc(iq)=uc(iq)*rotpoles(3,1)+vc(iq)*rotpoles(3,2)
     &         +wc(iq)*rotpoles(3,3)
      end do
      ! interpolate all required arrays to new C-C positions
      ! don't need to do map factors and Coriolis on target grid
      call doints4(ucc,  uct, nface4,xg4,yg4,nord,ik)
      call doints4(vcc,  vct, nface4,xg4,yg4,nord,ik)
      call doints4(wcc,  wct, nface4,xg4,yg4,nord,ik)
      do iq=1,ifull
        ! now convert to "target" Cartesian components (transpose used)
        ucc(iq)=uct(iq)*rotpole(1,1)+vct(iq)*rotpole(2,1)
     &                          +wct(iq)*rotpole(3,1)
        vcc(iq)=uct(iq)*rotpole(1,2)+vct(iq)*rotpole(2,2)
     &                          +wct(iq)*rotpole(3,2)
        wcc(iq)=uct(iq)*rotpole(1,3)+vct(iq)*rotpole(2,3)
     &                          +wct(iq)*rotpole(3,3)
        ! then finally to "target" local x-y components
        uct(iq) = ax(iq)*ucc(iq) + ay(iq)*vcc(iq) +
     &                  az(iq)*wcc(iq)
        vct(iq) = bx(iq)*ucc(iq) + by(iq)*vcc(iq) +
     &                  bz(iq)*wcc(iq)
      enddo               ! iq loop
      
      return
      end subroutine interpwind

      subroutine gethist1(ncid,iarchi,vname,ik,varout,iotest,
     &                    nface4,xg4,yg4,nord)
      
      use infile             ! Input file routines
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      
      integer, intent(in) :: ncid,iarchi,ik,nord
      integer ier
      real, dimension(ifull,4), intent(in) :: nface4,xg4,yg4
      real, dimension(ifull), intent(out) :: varout
      real, dimension(6*ik*ik) :: ucc
      logical, intent(in) :: iotest
      character(len=*), intent(in) :: vname
      
      if (iotest) then
        call histrd1(ncid,iarchi,ier,vname,ik,6*ik,varout,ifull)
      else
        call histrd1(ncid,iarchi,ier,vname,ik,6*ik,ucc,6*ik*ik)
        call doints4(ucc,varout,nface4,xg4,yg4,nord,ik)
      end if ! iotest

      return
      end subroutine gethist1
     
      subroutine filhist1(ncid,iarchi,vname,ik,dk,varout,iotest,
     &                    nface4,xg4,yg4,nord,mask_a)
      
      use infile             ! Input file routines
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      
      integer, intent(in) :: ncid,iarchi,ik,dk,nord
      integer ier
      real, dimension(ifull,4), intent(in) :: nface4,xg4,yg4
      real, dimension(ifull), intent(out) :: varout
      real, dimension(6*ik*ik) :: ucc
      logical, intent(in) :: iotest
      logical, dimension(6*dk*dk), intent(in) :: mask_a
      character(len=*), intent(in) :: vname
      
      if (iotest) then
        call histrd1(ncid,iarchi,ier,vname,ik,6*ik,varout,ifull)
      else
        call histrd1(ncid,iarchi,ier,vname,ik,6*ik,ucc,6*ik*ik)
        call fill_cc(ucc,dk,0,mask_a)
        call doints4(ucc,varout,nface4,xg4,yg4,nord,ik)
      end if ! iotest
      
      return
      end subroutine filhist1
     
      subroutine gethist4s(ncid,iarchi,vname,ik,kk,varout,iotest,
     &                     nface4,xg4,yg4,nord,sigin,vmode)
      
      use infile             ! Input file routines
      
      implicit none
      
      include 'newmpar.h'    ! Grid parameters
      
      integer, intent(in) :: ncid,iarchi,ik,kk,nord,vmode
      integer k,ier
      real, dimension(kk), intent(in) :: sigin
      real, dimension(ifull,4), intent(in) :: nface4,xg4,yg4
      real, dimension(ifull,kl), intent(out) :: varout
      real, dimension(6*ik*ik) :: ucc
      real, dimension(ifull,kk) :: u_k
      logical, intent(in) :: iotest
      character(len=*), intent(in) :: vname
      
      if (iotest) then
        do k=1,kk          
          call histrd4s(ncid,iarchi,ier,vname,ik,6*ik,k,u_k(:,k),ifull)
        end do
      else
        do k=1,kk
          call histrd4s(ncid,iarchi,ier,vname,ik,6*ik,k,ucc,6*ik*ik)
          call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,ik)
        end do
      end if ! iotest
      call vertint(u_k,varout,vmode,kk,sigin)
      
      return
      end subroutine gethist4s   
