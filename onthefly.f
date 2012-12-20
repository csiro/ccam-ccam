      ! Main netcdf input routines.  Host grid is automatically
      ! interpolated to nested model grid.  Three options are
      !   nested=0  Initial conditions
      !   nested=1  Nudging fields
      !   nested=2  Surface data recycling
      
      ! This version supports the parallel file routines contained
      ! in infile.f.  Hence restart files do not require any
      ! gathers and scatters.

      ! In the case where the grid needs to be interpolated, then
      ! the current code gathers the pieces of input from the
      ! IO read processors onto a single processor.  The single
      ! processor then interpolates the input data to the model
      ! grid.  The data is then distributed to all processors.
      ! This approach has smaller message sizes to be sent
      ! (compared to all processors interpolate locally, for
      ! which all processors need a copy of the global grid),
      ! but has the bottle neck of one processor performing
      ! all the interpolation.
      
      subroutine onthefly(nested,kdate_r,ktime_r,
     &                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &                    tgg,wb,wbice,snowd,qfg,qlg,qrg,
     &                    tggsn,smass,ssdn,ssdnn,snage,isflag,
     &                    iaero,mlodwn,ocndwn)

      use cc_mpi           ! CC MPI routines
      use infile           ! Input file routines
      use mlo              ! Ocean physics and prognostic arrays
      use soil_m           ! Soil and surface data

      implicit none

      include 'newmpar.h'  ! Grid parameters
      include 'darcdf.h'   ! Netcdf data
      include 'parm.h'     ! Model configuration
      include 'stime.h'    ! File date data

      integer, parameter :: ntest=0
      integer, parameter :: nihead=54
      integer, parameter :: nrhead=14

      integer, save :: ik,jk,kk,ok,maxarchi
      integer, save :: ncidold=-1
      integer, dimension(nihead) :: nahead
      integer, dimension(ifull) :: isflag
      integer, dimension(8) :: idum
      integer kdate_r,ktime_r,nested,ier,ier2,ilen,itype,iaero
      integer idv,mtimer,k,ierx,idvkd,idvkt,idvmt
      real, dimension(nrhead) :: ahead
      real, dimension(3) :: rdum
      real, save :: rlong0x, rlat0x, schmidtx
      real timer
!     These are local arrays, not the versions in arrays.h
!     Use in call to infile, so are dimensioned ifull rather than ifull_g
      real psl(ifull),zss(ifull),tss(ifull),fracice(ifull),
     & wb(ifull,ms),wbice(ifull,ms),snowd(ifull),sicedep(ifull),
     & t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     & tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     & ssdnn(ifull),snage(ifull),qfg(ifull,kl),
     & qlg(ifull,kl),qrg(ifull,kl),mlodwn(ifull,wlev,4),
     & ocndwn(ifull,2)
      logical ltest,newfile,tst

      call start_log(onthefly_begin)
      !--------------------------------------------------------------
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
          ik=nahead(1)
          jk=nahead(2)
          kk=nahead(3)
          rlong0x =ahead(5)
          rlat0x  =ahead(6)
          schmidtx=ahead(7)
          if(schmidtx<=0..or.schmidtx>1.)then
            rlong0x =ahead(6)
            rlat0x  =ahead(7)
            schmidtx=ahead(8)
          endif  ! (schmidtx<=0..or.schmidtx>1.)        
          maxarchi=0
          call ccnf_inq_dimlen(ncid,'time',maxarchi)
          ok=0
          if (nmlo/=0.and.abs(nmlo)<=9) then
            call ccnf_inq_dimlen(ncid,'olev',ok,failok=.true.)
          end if
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
        ierx=0
        if (tst) then
          ierx=1
          call ccnf_inq_varid(ncid,'timer',idvmt,tst)
        end if
        do while(ltest.and.iarchi<maxarchi)
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
          ltest=2400*(kdate_r-kdate_s)-1200*nsemble
     &              +(ktime_r-ktime_s)<0
        end do
        if (nsemble/=0) then
          kdate_r=kdate_s
          ktime_r=ktime_s
        end if
        if (ltest) then
          ktime_r=-1
        end if
        if (myid==0) then
          write(6,*) 'After search ltest,iarchi =',ltest,iarchi
          write(6,*) '             kdate_r,ktime_r =',kdate_r,ktime_r
        end if
        idum(1)=kdate_r
        idum(2)=ktime_r
        idum(3)=ncid
        idum(4)=ik
        idum(5)=jk
        idum(6)=kk
        idum(7)=ok
        idum(8)=iarchi
        rdum(1)=rlong0x
        rdum(2)=rlat0x
        rdum(3)=schmidtx
      endif  ! ( myid==0 .or. pfall )

      if (.not.pfall) then
        call ccmpi_bcast(idum(1:8),0,comm_world)
        kdate_r=idum(1)
        ktime_r=idum(2)
        ncid=idum(3)
        ik=idum(4)
        jk=idum(5)
        kk=idum(6)
        ok=idum(7)
        iarchi=idum(8)
      end if
      newfile=ncid/=ncidold
      if (newfile) then
        if (.not.pfall) then
          call ccmpi_bcast(rdum(1:3),0,comm_world)
          rlong0x=rdum(1)
          rlat0x=rdum(2)
          schmidtx=rdum(3)
        end if
        if (ncidold/=-1) then
          if (myid==0) then
            write(6,*) 'Closing old input file'
          end if
          call histclose
        end if
        ncidold=ncid
      end if

      if (ktime_r<0) then
        if (nested==2) then
          if (myid==0) then
            write(6,*) "WARN: Cannot locate date/time in input file"
          end if
          return
        end if
        write(6,*) "ERROR: Cannot locate date/time in input file"
        call ccmpi_abort(-1)
      end if
      !--------------------------------------------------------------
      
      ! Here we call ontheflyx with different automatic array
      ! sizes.  This means the arrays are correct for interpolation
      ! and file i/o on myid==0, as well as the arrays are smaller
      ! on myid/=0 when they are not needed.  This way we avoid
      ! having to maintain multiple ontheflyx subroutines.
      if (myid==0) then
        call ontheflyx(nested,kdate_r,ktime_r,
     &                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &                    tgg,wb,wbice,snowd,qfg,qlg,qrg,
     &                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     &                    ok,ik,ifull_g,iaero,mlodwn,ocndwn,rlong0x,
     &                    rlat0x,schmidtx,newfile)
        write(6,*) "Leaving onthefly"
      else
        call ontheflyx(nested,kdate_r,ktime_r,
     &                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &                    tgg,wb,wbice,snowd,qfg,qlg,qrg,
     &                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     &                    ok,0,0,iaero,mlodwn,ocndwn,rlong0x,
     &                    rlat0x,schmidtx,newfile)
      end if

      call end_log(onthefly_end)

      return
      end
      
      ! Read data from netcdf file
      
      ! arrays are typically read as global and then distributed to
      ! processor local arrays.  This allows for more flexibility
      ! with diagnosed fields.  Data is usually read in as 2D
      ! fields which avoids memory problems when the host grid
      ! size is significantly larger than the regional grid size.
      subroutine ontheflyx(nested,kdate_r,ktime_r,
     &                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &                    tgg,wb,wbice,snowd,qfg,qlg,qrg,
     &                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     &                    ok,dk,ifg,iaero,mlodwn,ocndwn,rlong0x,
     &                    rlat0x,schmidtx,newfile)
      
      use aerosolldr, only : xtg,ssn,naero      ! LDR aerosol scheme
      use ateb, only : atebdwn                  ! Urban
      use casadimension, only : mplant,mlitter, ! CASA dimensions
     &      msoil
      use carbpools_m                           ! Carbon pools
      use cc_mpi                                ! CC MPI routines
      use cfrac_m                               ! Cloud fraction
      use define_dimensions, only : ncs, ncp    ! CABLE dimensions
      use infile                                ! Input file routines
      use latlong_m                             ! Lat/lon coordinates
      use mlo, only : wlev,micdwn,mloregrid     ! Ocean physics and prognostic arrays
      use mlodynamics                           ! Ocean dynamics
      use morepbl_m                             ! Additional boundary layer diagnostics
      use nharrs_m, only : phi_nh,lrestart      ! Non-hydrostatic atmosphere arrays
      use nsibd_m, only : isoilm                ! Land-surface arrays
      use savuvt_m                              ! Saved dynamic arrays
      use savuv1_m                              ! Saved dynamic arrays
      use screen_m                              ! Screen level diagnostics
      use sigs_m                                ! Atmosphere sigma levels
      use soil_m                                ! Soil and surface data
      use tkeeps, only : tke,eps                ! TKE-EPS boundary layer
      use tracers_m                             ! Tracer data
      use utilities                             ! Grid utilities
      use vecsuv_m                              ! Map to cartesian coordinates
      use vvel_m, only : dpsldt,sdot            ! Additional vertical velocity
      use xarrs_m, only : pslx                  ! Saved dynamic arrays
      use workglob_m                            ! Additional grid interpolation
      use work2_m                               ! Diagnostic arrays

      implicit none

      include 'newmpar.h'                       ! Grid parameters
      include 'const_phys.h'                    ! Physical constants
      include 'darcdf.h'                        ! Netcdf data
      include 'parm.h'                          ! Model configuration
      include 'parmdyn.h'                       ! Dynamics parmaters
      include 'parmgeom.h'                      ! Coordinate data
      include 'soilv.h'                         ! Soil parameters
      include 'stime.h'                         ! File date data

      integer, parameter :: ntest=0
      integer, parameter :: nord=3        ! 1 for bilinear, 3 for bicubic
      
      integer ik, kk, ok, idv, iaero, isoil, nud_test
      integer dk, ifg ! controls automatic array size
      integer lev,levkk,ier,ierr,igas
      integer kdate_r, ktime_r, nemi, id2,jd2,idjd2
      integer nested, i, j, k, mm, iq, ii, jj, np, numneg
      integer, dimension(:,:), allocatable, save :: nface4
      integer, dimension(ifull) :: isflag
      integer, dimension(:), allocatable, save :: isoilm_a
      integer, dimension(2*ms) :: iera
      integer, dimension(ms) :: ierb
      integer, dimension(6) :: ierc
      integer, dimension(3), save :: iers
      integer, dimension(2) :: dumb
      real*8, dimension(1+4*dk,1+4*dk) :: xx4,yy4
      real*8, dimension(dk*dk*6):: z_a,x_a,y_a
      real, dimension(ifull,wlev,4) :: mlodwn
      real, dimension(ifull,ok,2) :: mloin
      real, dimension(ifull,2) :: ocndwn
      real, dimension(ifull,ms) :: wb,wbice,tgg
      real, dimension(ifull,3) :: tggsn,smass,ssdn
      real, dimension(ifull,kl) :: t,u,v,qg,qfg,qlg,qrg
      real, dimension(ifull,kl) :: dum
      real, dimension(ifull,kk) :: u_k,v_k
      real, dimension(ifull) :: psl,zss,tss,fracice
      real, dimension(ifull) :: snowd,sicedep,ssdnn,snage,dum6
      real, dimension(ifull) :: tss_l, tss_s, pmsl
      real, dimension(:), allocatable, save :: sigin,zss_a,ocndep_l
      real, dimension(:,:), allocatable, save :: xg4, yg4
      real, dimension(ifg) :: uct_g, vct_g
      real, dimension(dk*dk*6) :: tss_a,fracice_a,sicedep_a
      real, dimension(dk*dk*6) :: psl_a,pmsl_a
      real, dimension(dk*dk*6) :: tss_l_a,tss_s_a
      real, dimension(dk*dk*6) :: t_a_lev,ucc,vcc
      real, dimension(dk*dk*6) :: wts_a  ! not used here or defined in call setxyz
      real, dimension(:), allocatable, save :: axs_a,ays_a,azs_a
      real, dimension(:), allocatable, save :: bxs_a,bys_a,bzs_a
      real, dimension(3,3), save :: rotpoles,rotpole
      real, intent(in) :: rlong0x, rlat0x, schmidtx
      real rlongd,rlatd
      character*8 vname
      character*3 trnum
      logical, dimension(dk*dk*6) :: land_a
      logical iotest,newfile,tsstest,tst

      real, parameter :: spval=999.  ! missing value flag
      real, parameter :: iotol=1.E-5 ! tolarance for iotest

      ! internal check (should not occur if code is written correctly)
      if (myid==0.and.ik/=dk) then
        write(6,*) "ERROR: Incorrect automatic array size in onthefly"
        stop
      end if
      
      ! land-sea mask method (nemi=3 use soilt, nemi=2 use tgg, nemi=1 use zs)
      nemi=3
      
      ! test if retopo fields are required
      nud_test=1
      if (nud_p==0.and.nud_t==0.and.nud_q==0) nud_test=0
      
      ! Determine if interpolation is required
      iotest=6*ik*ik==ifull_g.and.abs(rlong0x-rlong0)<iotol.and.
     &       abs(rlat0x-rlat0)<iotol.and.
     &       abs(schmidtx-schmidt)<iotol
      if (iotest) then
        io_in=1   ! no interpolation
      else
        io_in=-1  ! interpolation
      end if
      if (myid==0) write(6,*) "Interpolation iotest,io_in =",
     &                                       iotest,io_in

      !--------------------------------------------------------------
      ! Determine input grid coordinates and interpolation arrays
      if (newfile) then
       if (allocated(nface4)) then
         deallocate(nface4,xg4,yg4)
       end if
       allocate(nface4(ifg,4),xg4(ifg,4),yg4(ifg,4))
       if ( myid==0 ) then
        if (allocated(axs_a)) then
          deallocate(axs_a,ays_a,azs_a)
          deallocate(bxs_a,bys_a,bzs_a)
        end if
        allocate(axs_a(dk*dk*6),ays_a(dk*dk*6),azs_a(dk*dk*6))
        allocate(bxs_a(dk*dk*6),bys_a(dk*dk*6),bzs_a(dk*dk*6))
        write(6,*) "Defining input file grid"
        if (ifg/=ifull_g) then
          write(6,*) "ERROR: Internal array allocation error"
          call ccmpi_abort(-1)
        end if
        if (m_fly==1) then
          rlong4(:,1)=rlongg_g(:)*180./pi
          rlat4(:,1)=rlatt_g(:)*180./pi
        end if
!       N.B. -ve ik in call setxyz preserves TARGET rlat4, rlong4     
!       following setxyz call is for source data geom    ****   
        do iq=1,ik*ik*6
         axs_a(iq)=iq
         ays_a(iq)=iq
         azs_a(iq)=iq
        enddo      
        call setxyz(ik,rlong0x,rlat0x,-schmidtx,
     &   x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4,
     &   myid)
        rotpoles = calc_rotpole(rlong0x,rlat0x)
        if(ktau<3)then
           write(6,*)'m_fly,nord ',m_fly,nord
           write(6,*)'kdate_r,ktime_r,ktau,ds',
     &              kdate_r,ktime_r,ktau,ds
           write(6,*)'rotpoles:'
           do i=1,3
              write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)')
     &          (i,j,j=1,3),(rotpoles(i,j),j=1,3)
           enddo
        endif                  ! (ktau<3)

!       rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!       rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!       rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
        rotpole = calc_rotpole(rlong0,rlat0)
        if(nmaxpr==1)then   ! already in myid==0 loop
           write(6,*)'in onthefly rotpole:'
           do i=1,3
              write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)')
     &          (i,j,j=1,3),(rotpole(i,j),j=1,3)
           enddo
           write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
           write(6,*)'before latltoij for id,jd: ',id,jd
           if ( nproc==1 ) then
              ! Diagnostics will only be correct if nproc==1
              write(6,*)'rlong4(1-4) ',(rlong4(idjd,mm),mm=1,4)
              write(6,*)'rlat4(1-4) ',(rlat4(idjd,mm),mm=1,4)
           end if
           write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,
     &                schmidtx
        endif                  ! (nmaxpr==1)
        do mm=1,m_fly  !  was 4, now may be set to 1 in namelist
           do iq=1,ifull_g
              call latltoij(rlong4(iq,mm),rlat4(iq,mm),          !input
     &                      rlong0x,rlat0x,schmidtx,             !input
     &                      xg4(iq,mm),yg4(iq,mm),nface4(iq,mm), !output (source)
     &                      xx4,yy4,ik)
           enddo
        enddo
        if(nproc==1.and.nmaxpr==1)then
          ! Diagnostics will only be correct if nproc==1
          id2=nint(xg4(idjd,1))
          jd2=il*nface4(idjd,1)+nint(yg4(idjd,1))
          idjd2=id2+il*(jd2-1)
          write(6,*)'after latltoij giving id2,jd2,idjd2: ',
     .                                   id2,jd2,idjd2
          write(6,*)'nface4(1-4) ',(nface4(idjd,mm),mm=1,4)
          write(6,*)'xg4(1-4) ',(xg4(idjd,mm),mm=1,4)
          write(6,*)'yg4(1-4) ',(yg4(idjd,mm),mm=1,4)
          if(nested==0)then
             write(6,"('wb_s(1)#  ',9f7.3)") 
     .           ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
             write(6,"('wb_s(ms)# ',9f7.3)") 
     .           ((wb(ii+(jj-1)*il,ms),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
          endif  ! (nested==0)
        endif
       end if ! (myid==0)
      end if ! newfile

      
      ! special data read for new file
      ! read once when file is first opened
      ! need global zss_a for (potentially) landsea mask and psl interpolation
      ! need global isoilm_a for (potentially) landsea mask
      if (newfile) then
        if (allocated(sigin)) deallocate(sigin)
        allocate(sigin(kk))
        if (myid==0.or.pfall) then
          if (myid==0) then
            write(6,*) "Reading fixed fields"
          end if
          call ccnf_inq_varid(ncid,'lev',idv,tst)
          if (tst) then
            call ccnf_inq_varid(ncid,'layer',idv,tst)
          end if
          if (.not.tst) then
            dumb(1)=1
            dumb(2)=kk
            call ccnf_get_vara(ncid,idv,dumb(1:1),dumb(2:2),sigin)
          else
            call ccnf_get_attg(ncid,'sigma',sigin)
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
        if (.not.pfall) then
          call ccmpi_bcast(sigin,0,comm_world)
          call ccmpi_bcast(iers(1:3),0,comm_world)
        end if
      end if ! newfile
      tsstest=(iers(2)==0.and.iers(3)==0.and.iotest)
      if (newfile) then
        if (myid==0) then
          write(6,*) "tsstest,iers ",tsstest,iers(1:3)
        end if      
        if (allocated(zss_a)) deallocate(zss_a)
        if (allocated(isoilm_a)) deallocate(isoilm_a)
        if (tsstest) then
          allocate(zss_a(ifull))
          call histrd1(ncid,iarchi,ier,'zht',ik,6*ik,zss_a,ifull)
        else
          allocate(zss_a(6*dk*dk))
          if (myid==0) then
            allocate(isoilm_a(6*dk*dk))
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
          ocndep_l=0.
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'ocndepth',ik,6*ik,ocndep_l,
     &                   ifull)
          else
            call histrd1(ncid,iarchi,ier,'ocndepth',ik,6*ik,ucc,
     &                   6*ik*ik)
            call doints4(ucc,ocndep_l,nface4,xg4,yg4,nord,dk,ifg)
          end if ! iotest
        end if
        if (myid==0) then
          write(6,*) "Finished reading fixed fields"
        end if
      else
        if (myid==0) then
          write(6,*) "Using saved fixed fields"
        end if
      endif

#ifdef debug
      ! internal errors which should not occur if code is written correctly
      if (.not.allocated(sigin)) then
        write(6,*) "ERROR: sigin is undefined in onthefly"
        call ccmpi_abort(-1)
      end if
      if (nmlo/=0.and.abs(nmlo)<=9
     &    .and..not.allocated(ocndep_l)) then
        write(6,*) "ERROR: ocndep_l is undefined in onthefly"
        call ccmpi_abort(-1)
      end if
      if (myid==0) then
        if (.not.allocated(zss_a)) then
          write(6,*) "ERROR: zss_a is undefined in onthefly"
          call ccmpi_abort(-1)
        end if
        if (.not.allocated(isoilm_a).and..not.tsstest) then
          write(6,*) "ERROR: isoilm_a is undefined in onthefly"
          call ccmpi_abort(-1)
        end if
      end if
#endif

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
        if(myid==0)then
         if(nemi==3)then 
          land_a(:)=isoilm_a(:)>0
          numneg=count(.not.land_a)
          if (any(isoilm_a(:)<0)) nemi=2
         end if
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
         tss_a=abs(tss_a)
         if(nemi==1)then
          land_a(:) = zss_a(:) > 0.
          numneg=count(.not.land_a)
         endif
         write(6,*)'Land-sea mask using nemi = ',nemi
        end if
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
        mlodwn(:,:,1)=293.
        if ((nested/=1.or.nud_sst/=0).and.ok>0) then
          do k=1,ok
            if (k<=ms) then
              write(vname,'("tgg",I1.1)') k
            else
              write(vname,'("tgg",I2.2)') k
            end if     
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     mloin(:,k,1),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (land_a)
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
!               interpolate all required arrays to new C-C positions
                call ints4(ucc,uct_g, nface4,xg4,yg4,nord,ik)
                call ccmpi_distribute(mloin(:,k,1),uct_g)
              else ! myid /= 0
                call ccmpi_distribute(mloin(:,k,1))
              endif ! myid==0
            end if ! iotest
          end do
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,1),0)
          if (all(mlodwn(:,:,1)==0.)) mlodwn(:,:,1)=293.
        end if ! (nestesd/=1.or.nud_sst/=0) ..else..
        ! ocean salinity
        mlodwn(:,:,2)=34.72
        if ((nested/=1.or.nud_sss/=0).and.ok>0) then
          do k=1,ok
            write(vname,'("sal",I2.2)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     mloin(:,k,1),ifull)
            else
              ucc=34.72
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (land_a)
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
!               interpolate all required arrays to new C-C positions
                call ints4(ucc,uct_g, nface4,xg4,yg4,nord,ik)
                call ccmpi_distribute(mloin(:,k,1),uct_g)              
              else ! myid /= 0
                call ccmpi_distribute(mloin(:,k,1))            
              endif ! myid==0
            end if ! iotest
          end do
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,2),1)
          if (all(mlodwn(:,:,2)==0.)) mlodwn(:,:,2)=34.72
        end if ! (nestesd/=1.or.nud_sss/=0) ..else..
        ! ocean currents
        mlodwn(:,:,3:4)=0.
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
              if (myid==0) then
                where (land_a)
                  ucc=spval
                  vcc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
                call fill_cc(vcc,spval,ik,0)
                call interpwind(ik,uct_g,vct_g,ucc,vcc,axs_a,ays_a,
     &                          azs_a,bxs_a,bys_a,bzs_a,rotpole,
     &                          rotpoles,nface4,xg4,yg4,nord)
!               interpolate all required arrays to new C-C positions
                call ccmpi_distribute(mloin(:,k,1), uct_g)
                call ccmpi_distribute(mloin(:,k,2), vct_g)
              else ! myid /= 0
                call ccmpi_distribute(mloin(:,k,1))
                call ccmpi_distribute(mloin(:,k,2))
              endif ! myid==0
            end if ! iotest
          end do
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,3),2)
          call mloregrid(ok,ocndwn(:,1),mloin(:,:,2),mlodwn(:,:,4),3)
        end if ! (nestesd/=1.or.nud_ouv/=0) ..else..
        ! water surface height
        ocndwn(:,2)=0.
        if (nested/=1.or.nud_sfh/=0) then
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'ocheight',ik,6*ik,
     &                   ocndwn(:,2),ifull)
          else
            ucc=0.
            call histrd1(ncid,iarchi,ier,'ocheight',ik,6*ik,ucc,
     &                   6*ik*ik)
            if (myid==0) then
              where (land_a)
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,ocndwn(:,2),nface4,xg4,yg4,nord,dk,ifg)
          end if ! iotest
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
         do iq=1,ik*ik*6
          if(land_a(iq))then       ! over land
            tss_l_a(iq)=tss_a(iq)
            tss_s_a(iq)=spval
            sicedep_a(iq)=spval
            fracice_a(iq)=spval
          else                   ! over sea
            tss_s_a(iq)=abs(tss_a(iq))
            tss_l_a(iq)=spval
          endif  !   (land_a(iq)) .. else ..
         enddo     ! iq loop
         call fill_cc(tss_l_a,spval,ik,0)
         call fill_cc(tss_s_a,spval,ik,0)
         call fill_cc(sicedep_a,spval,ik,0)
         call fill_cc(fracice_a,spval,ik,0)
        endif   ! (myid==0)

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
          tss=tss_s
          where (land)
            sicedep=0.
            fracice=0.
            tss=tss_l
          end where
        else
!         The routine doints4 does the gather, calls ints4 and redistributes
          call doints4(zss_a , zss,  nface4,xg4,yg4,nord,dk,ifg)
          call doints4(tss_l_a , tss_l,  nface4,xg4,yg4,nord,dk,ifg)
          call doints4(tss_s_a , tss_s,  nface4,xg4,yg4,nord,dk,ifg)
          call doints4(fracice_a , fracice,  nface4,xg4,yg4,nord,dk,ifg)
          call doints4(sicedep_a , sicedep,  nface4,xg4,yg4,nord,dk,ifg)
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
            if (k==levkk) t_a_lev=ucc ! store for psl calculation below
            if (myid==0) then
              call ints4(ucc,uct_g,nface4,xg4,yg4,nord,ik)  ! ints4 on source grid
              call ccmpi_distribute(u_k(:,k),uct_g)
            else ! myid /= 0
              call ccmpi_distribute(u_k(:,k))
            endif ! myid==0
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
          ! but makes it possible to downscale large grids (e.g., C160)
          if (iotest) then
            call histrd4s(ncid,iarchi,ier,'u',ik,6*ik,k,u_k(:,k),
     &                    ifull)                                        !     u wind component
            call histrd4s(ncid,iarchi,ier,'v',ik,6*ik,k,v_k(:,k),
     &                    ifull)                                        !     v wind component
          else
            call histrd4s(ncid,iarchi,ier,'u',ik,6*ik,k,ucc,6*ik*ik)    !     u wind component
            call histrd4s(ncid,iarchi,ier,'v',ik,6*ik,k,vcc,6*ik*ik)    !     v wind component
            if (myid==0) then
              call interpwind(ik,uct_g,vct_g,ucc,vcc,axs_a,ays_a,
     &                        azs_a,bxs_a,bys_a,bzs_a,rotpole,
     &                        rotpoles,nface4,xg4,yg4,nord)
!             interpolate all required arrays to new C-C positions
!             don't need to do map factors and Coriolis on target grid
              call ccmpi_distribute(u_k(:,k), uct_g)
              call ccmpi_distribute(v_k(:,k), vct_g)
            else
              call ccmpi_distribute(u_k(:,k))
              call ccmpi_distribute(v_k(:,k))
            endif ! myid==0
          end if ! iotest
        enddo  ! k loop
        call vertint(u_k ,u, 3,kk,sigin)
        call vertint(v_k ,v, 4,kk,sigin)
      end if ! (nested==0.or.(nested==1.and.nud_uv/=0))
      ! mixing ratio
      ! read for nested=0 or nested=1.and.nud_q/=0
      qg=qgmin
      if (nested==0.or.(nested==1.and.nud_q/=0)) then
        do k=1,kk
          if (iotest) then
            if (iers(1)==0) then
              call histrd4s(ncid,iarchi,ier,'mixr',ik,6*ik,k,u_k(:,k),
     &                      ifull)                                       !     mixing ratio
            else
              call histrd4s(ncid,iarchi,ier,'q',ik,6*ik,k,u_k(:,k),
     &                      ifull)                                       !     mixing ratio
            endif  ! (ier/=0)
          else
            if (iers(1)==0) then
              call histrd4s(ncid,iarchi,ier,'mixr',ik,6*ik,k,ucc,
     &                      6*ik*ik)                                     !     mixing ratio
            else
              call histrd4s(ncid,iarchi,ier,'q',ik,6*ik,k,ucc,
     &                      6*ik*ik)                                     !     mixing ratio
            endif  ! (ier/=0)
            if (myid==0) then
              call ints4(ucc,uct_g, nface4,xg4,yg4,nord,ik)
              call ccmpi_distribute(u_k(:,k),uct_g)
            else
              call ccmpi_distribute(u_k(:,k))         
            endif ! myid==0
          end if ! iotest
        enddo  ! k loop
        call vertint(u_k,qg,2,kk,sigin)
      end if ! (nested==0.or.(nested==1.and.nud_q/=0))

      ! re-grid surface pressure by mapping to MSLP, interpolating and then map to surface pressure
      ! requires psl_a, zss, zss_a, t and t_a_lev
      if (nested==0.or.(nested==1.and.nud_test/=0)) then
        if (.not.iotest) then
!         The routine doints4 does the gather, calls ints4 and redistributes
          if ( myid==0 ) then
            call mslpx(pmsl_a,psl_a,zss_a,t_a_lev,ik*ik*6,sigin(levkk))  ! needs pmsl (preferred)
            call ints4(pmsl_a,uct_g, nface4,xg4,yg4,nord,ik)
            call ccmpi_distribute(pmsl,uct_g)
          else
            call ccmpi_distribute(pmsl)
          end if ! myid==0
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
            call ccnf_inq_varid(ncid,'u10',idv,tst)
            ierc(1:2)=0
            if (lrestart) ierc(1)=1
            if (tst) ierc(2)=-1
          end if
          if (.not.pfall) then
            call ccmpi_bcast(ierc(1:6),0,comm_world)
          end if
          lrestart=(ierc(1)==1)
          if (lrestart) then
            nstag=ierc(3)
            nstagu=ierc(4)
            nstagoff=ierc(5)
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
        end if ! nested==0
        !--------------------------------------------------

        ! SNOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (iotest) then
          call histrd1(ncid,iarchi,ier,'snd',ik,6*ik,snowd,ifull)
        else
          call histrd1(ncid,iarchi,ier,'snd',ik,6*ik,ucc,6*ik*ik)
          if(myid==0)then
            where (.not.land_a)       
                ucc(:)=spval
            end where  !   (.not.land_a(iq)) 
            call fill_cc(ucc,spval,ik,0)
          endif  ! (myid==0)
          call doints4(ucc,  snowd,nface4,xg4,yg4,nord,dk,ifg)
        end if ! iotest

        ! SOIL TEMPERATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (myid==0.or.pfall) then
          iera(1:ms)=0
          do k=1,ms
            write(vname,'("tgg",I1.1)') k
            call ccnf_inq_varid(ncid,vname,idv,tst)
            if (tst) iera(k)=-1
          end do
        end if
        if (.not.pfall) then
          call ccmpi_bcast(iera(1:ms),0,comm_world)
        end if
        do k=1,ms 
          if (iera(k)==0) then
            write(vname,'("tgg",I1.1)') k
          else if (k<=3.and.iera(2)==0) then
            vname="tgg2"
          else if (k<=3) then
            vname="tb3"
          else if (iera(6)==0) then
            vname="tgg6"
          else
            vname="tb2"
          end if
          if (iotest) then
            if (k==1.and.iera(1)/=0) then
              tgg(:,1)=tss
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     tgg(:,k),ifull)
            end if
          else
            if (k==1.and.iera(1)/=0) then
              ucc=tss_a
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
            end if
            if (myid==0) then
              where (.not.land_a(:))
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,tgg(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
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
                if (iotest) then
                  call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                         micdwn(:,k),ifull)
                else
                  call histrd1(ncid,iarchi,ier,vname,ik,6*ik,ucc,
     &                   6*ik*ik)
                  if (myid==0) then
                    where (land_a)
                      ucc=spval
                    end where
                    call fill_cc(ucc,spval,ik,0)
                  end if
                  call doints4(ucc,micdwn(:,k),nface4,xg4,yg4,nord,dk,
     &                         ifg)
                end if
                if (all(micdwn(:,k)==0.)) micdwn(:,k)=280.
              case(5)
                micdwn(:,k)=fracice ! read above with nudging arrays
              case(6)
                micdwn(:,k)=sicedep ! read above with nudging arrays
              case(7)
                micdwn(:,k)=snowd*1.E-3
            end select
          end do
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'sto',ik,6*ik,micdwn(:,k),
     &                 ifull)
          else
            call histrd1(ncid,iarchi,ier,'sto',ik,6*ik,ucc,
     &                 6*ik*ik)
            if (myid==0) then
              where (land_a)
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,micdwn(:,k),nface4,xg4,yg4,nord,dk,
     &               ifg)
          end if
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'uic',ik,6*ik,micdwn(:,9),
     &                   ifull)
            call histrd1(ncid,iarchi,ier,'vic',ik,6*ik,micdwn(:,10),
     &                   ifull)
          else
            call histrd1(ncid,iarchi,ier,'uic',ik,6*ik,ucc,6*ik*ik)
            call histrd1(ncid,iarchi,ier,'vic',ik,6*ik,vcc,6*ik*ik)
            if (myid==0) then
              where (land_a)
                ucc=spval
                vcc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
              call fill_cc(vcc,spval,ik,0)
              call interpwind(ik,uct_g,vct_g,ucc,vcc,axs_a,ays_a,azs_a,
     &                      bxs_a,bys_a,bzs_a,rotpole,rotpoles,nface4,
     &                      xg4,yg4,nord)
              call ccmpi_distribute(micdwn(:,9), uct_g)
              call ccmpi_distribute(micdwn(:,10), vct_g)
            else ! myid /= 0
              call ccmpi_distribute(micdwn(:,9))
              call ccmpi_distribute(micdwn(:,10))
            endif ! myid==0
          end if ! iotest
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'icesal',ik,6*ik,
     &                   micdwn(:,11),ifull)
          else
            call histrd1(ncid,iarchi,ier,'icesal',ik,6*ik,ucc,
     &                   6*ik*ik)
            if (myid==0) then
              where (land_a)
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,micdwn(:,11),nface4,xg4,yg4,nord,dk,
     &             ifg)
          end if
          if (abs(nmlo)>=2.and.abs(nmlo)<=9) then
            if (iotest) then
              call histrd1(ncid,iarchi,ier,'swater',ik,6*ik,
     &                     watbdy(1:ifull),ifull)
            else
              call histrd1(ncid,iarchi,ier,'swater',ik,6*ik,ucc,
     &                     6*ik*ik)
              call doints4(ucc,watbdy(1:ifull),nface4,xg4,yg4,nord,dk,
     &               ifg)
            end if
            if (iotest) then
              call histrd1(ncid,iarchi,ier,'ssalin',ik,6*ik,
     &                     salbdy(1:ifull),ifull)
            else
              call histrd1(ncid,iarchi,ier,'ssalin',ik,6*ik,ucc,
     &                     6*ik*ik)
              call doints4(ucc,salbdy(1:ifull),nface4,xg4,yg4,nord,dk,
     &               ifg)
            end if
          end if
        end if
        !--------------------------------------------------

        ! SOIL MOISTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (myid==0.or.pfall) then
          iera(1:ms)=0
          ierb(1:ms)=0
          do k=1,ms
            write(vname,'("wetfrac",I1.1)') k
            call ccnf_inq_varid(ncid,vname,idv,tst)
            if (tst) iera(k)=-1
            write(vname,'("wb",I1.1)') k
            call ccnf_inq_varid(ncid,vname,idv,tst)
            if (tst) ierb(k)=-1
          end do
        end if
        if (.not.pfall) then
          iera(ms+1:2*ms)=ierb(1:ms)
          call ccmpi_bcast(iera(1:2*ms),0,comm_world)
          ierb(1:ms)=iera(ms+1:2*ms)
        end if
        wb=20.5
        do k=1,ms
          if (iera(k)==0) then
            write(vname,'("wetfrac",I1.1)') k
          else if (ierb(k)==0) then
            write(vname,'("wb",I1.1)') k
          else if (k<2.and.ierb(2)==0) then
            vname="wb2"
          else if (k<2) then
            vname="wfg"
          else if (ierb(6)==0) then
            vname="wb6"
          else
            vname="wfb"
          end if
          if (iotest) then
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   wb(:,k),ifull)
            if (iera(k)==0) then
              wb(:,k)=wb(:,k)+20. ! flag for fraction of field capacity
            end if
          else
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   ucc,6*ik*ik)
            if (iera(k)==0) then
              ucc=ucc+20.
            end if
            if (myid==0) then
              where (.not.land_a(:))
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,wb(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
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

        ! PH - Add wetfac to output for mbase=-19 option
        if (iotest) then
          call histrd1(ncid,iarchi,ier,'wetfac',ik,6*ik,wetfac,
     &                 ifull)
        else
          call histrd1(ncid,iarchi,ier,'wetfac',ik,6*ik,ucc,6*ik*ik)
          if (myid==0) then
            where (.not.land_a(:))
              ucc=spval
            end where
            call fill_cc(ucc,spval,ik,0)
          end if
          call doints4(ucc,wetfac,nface4,xg4,yg4,nord,dk,ifg)
        end if ! iotest
        where (.not.land)
          wetfac=1.
        end where

        !--------------------------------------------------
        ! Read 10m wind speeds for special sea roughness length calculations
        if (nested==0) then
          if (ierc(2)==0) then
           if (iotest) then
            call histrd1(ncid,iarchi,ier,'u10',ik,6*ik,u10,ifull)
           else
            call histrd1(ncid,iarchi,ier,'u10',ik,6*ik,ucc,6*ik*ik)
            call doints4(ucc,u10,nface4,xg4,yg4,nord,dk,ifg)
           end if ! iotest
          else
           u10=sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)*log(10./0.001)
     &                                              /log(zmin/0.001)
          end if
        end if

        !--------------------------------------------------
        ! Read boundary layer height for TKE-eps mixing
        if (nvmix==6.and.nested==0) then
          if (iotest) then
            call histrd1(ncid,iarchi,ier,'pblh',ik,6*ik,pblh,ifull)
          else
            call histrd1(ncid,iarchi,ier,'pblh',ik,6*ik,ucc,6*ik*ik)
            call doints4(ucc,pblh,nface4,xg4,yg4,
     &                     nord,dk,ifg)
          end if ! iotest
          if (all(pblh==0.)) pblh=1000.
        end if

        !--------------------------------------------------
        ! Read CABLE/CASA aggregate C+N+P pools
        if (nsib>=6) then
         if (ccycle==0) then
          if (myid==0.or.pfall) then
            call ccnf_inq_varid(ncid,'cplant1',idv,tst)
            iera(1)=0
            if (tst) iera(1)=-1
          end if
          if (.not.pfall) then
            call ccmpi_bcast(iera(1:1),0,comm_world)
          end if
          if (iera(1)==0) then
           do k=1,ncp
            write(vname,'("cplant",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     cplant(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,cplant(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
           end do
           do k=1,ncs
            write(vname,'("csoil",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     csoil(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,csoil(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
            end if ! iotest
           end do
          end if
         else
          if (myid==0.or.pfall) then
            call ccnf_inq_varid(ncid,'glai',idv,tst)
            iera(1)=0
            if (tst) iera(1)=-1
          end if
          if (.not.pfall) then
            call ccmpi_bcast(iera(1:1),0,comm_world)
          end if
          if (iera(1)==0) then
           do k=1,mplant
            write(vname,'("cplant",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     cplant(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,cplant(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
            write(vname,'("nplant",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     niplant(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,niplant(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
            write(vname,'("pplant",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     pplant(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,pplant(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
           end do
           do k=1,mlitter
            write(vname,'("clitter",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     clitter(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,clitter(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
            write(vname,'("nlitter",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     nilitter(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,nilitter(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
            write(vname,'("plitter",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     plitter(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,plitter(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
            end if ! iotest
           end do         
           do k=1,msoil
            write(vname,'("csoil",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     csoil(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,csoil(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
            write(vname,'("nsoil",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     nisoil(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,nisoil(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
            write(vname,'("psoil",I1.1)') k
            if (iotest) then
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     psoil(:,k),ifull)
            else
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ucc,6*ik*ik)
              if (myid==0) then
                where (.not.land_a(:))
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,psoil(:,k),nface4,xg4,yg4,
     &                   nord,dk,ifg)
            end if ! iotest
           end do
           if (iotest) then
            call histrd1(ncid,iarchi,ier,'glai',ik,6*ik,
     &                   glai,ifull)
           else
            call histrd1(ncid,iarchi,ier,'glai',ik,6*ik,
     &                   ucc,6*ik*ik)
            if (myid==0) then
              where (.not.land_a(:))
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,glai,nface4,xg4,yg4,
     &                 nord,dk,ifg)
           end if ! iotest
          end if ! iera(1)==0
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
              if (myid==0) then
                where (.not.land_a.or.ucc>=399.)
                  ucc=spval
                end where
                call fill_cc(ucc,spval,ik,0)
              end if
              call doints4(ucc,atebdwn(:,k),nface4,xg4,yg4,nord,
     &                     dk,ifg)
            end if ! iotest
            if (all(atebdwn(:,k)==0.)) then
              select case(k)
                case(21,22)
                  atebdwn(:,k)=100.
                case(23,24)
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
        end if

        if (nested==0) then
          ! CLOUD FROZEN WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k=1,kk
           if (iotest) then
             call histrd4s(ncid,iarchi,ier,'qfg',ik,6*ik,k,u_k(:,k),
     &                     ifull) 
           else
             ucc=0. ! dummy for qfg
             call histrd4s(ncid,iarchi,ier,'qfg',ik,6*ik,k,ucc,
     &                     6*ik*ik)
             call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,dk,ifg)
           end if ! iotest
          enddo  ! k loop
          call vertint(u_k,dum,5,kk,sigin)
          qfg(1:ifull,:)=dum
          ! CLOUD LIQUID WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k=1,kk 
           if (iotest) then
             call histrd4s(ncid,iarchi,ier,'qlg',ik,6*ik,k,u_k(:,k),
     &                     ifull)
           else
             ucc=0. ! dummy for qlg
             call histrd4s(ncid,iarchi,ier,'qlg',ik,6*ik,k,ucc,
     &                     6*ik*ik)
             call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,dk,ifg)
           end if ! iotest
          enddo  ! k loop
          call vertint(u_k,dum,5,kk,sigin)
          qlg(1:ifull,:)=dum
          ! RAIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k=1,kk
           if (iotest) then
             call histrd4s(ncid,iarchi,ier,'qrg',ik,6*ik,k,v_k(:,k),
     &                     ifull)
           else 
             vcc=0. ! dummy for qrg
             call histrd4s(ncid,iarchi,ier,'qrg',ik,6*ik,k,vcc,
     &                     6*ik*ik)
             call doints4(vcc,v_k(:,k),nface4,xg4,yg4,nord,dk,ifg)
           end if ! iotest
          enddo  ! k loop
          call vertint(v_k,dum,5,kk,sigin)
          qrg(1:ifull,:)=dum
    !      ! GRAUPLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      do k=1,kk
    !       if (iotest) then
    !         call histrd4s(ncid,iarchi,ier,'qgrau',ik,6*ik,k,
    ! &                     v_k(:,k),ifull)
    !       else 
    !         call histrd4s(ncid,iarchi,ier,'qgrau',ik,6*ik,k,vcc,
    ! &                     6*ik*ik)
    !         call doints4(vcc,v_k(:,k),nface4,xg4,yg4,nord,dk,ifg)
    !       end if ! iotest
    !      enddo  ! k loop
    !      call vertint(v_k,dum,5,kk,sigin)
    !      qgrau(1:ifull,:)=dum
          ! CLOUD FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k=1,kk
           if (iotest) then
             call histrd4s(ncid,iarchi,ier,'cfrac',ik,6*ik,k,
     &                     u_k(:,k),ifull)
           else 
             ucc=0. ! dummy for cfrac
             call histrd4s(ncid,iarchi,ier,'cfrac',ik,6*ik,k,ucc,
     &                     6*ik*ik)
             call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,dk,ifg)
           end if ! iotest
          enddo  ! k loop
          call vertint(u_k,dum,5,kk,sigin)
          cfrac(1:ifull,:)=dum
          ! RAIN FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k=1,kk
           if (iotest) then
             call histrd4s(ncid,iarchi,ier,'cfrain',ik,6*ik,k,
     &                     u_k(:,k),ifull)
           else 
             ucc=0. ! dummy for cffall
             call histrd4s(ncid,iarchi,ier,'cfrain',ik,6*ik,k,ucc,
     &                     6*ik*ik)
             call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,dk,ifg)
           end if ! iotest
          enddo  ! k loop
          call vertint(u_k,dum,5,kk,sigin)
          cffall(1:ifull,:)=dum
        end if ! (nested==0)

        !--------------------------------------------------
        ! TKE-eps data
        if (nvmix==6.and.nested==0) then
          do k=1,kk
            if (iotest) then
              call histrd4s(ncid,iarchi,ier,'tke',ik,6*ik,k,
     &                      u_k(:,k),ifull)
            else
              call histrd4s(ncid,iarchi,ier,'tke',ik,6*ik,k,
     &                    ucc,6*ik*ik)
              call doints4(ucc,u_k(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
            end if ! iotest
          end do
          call vertint(u_k,dum,5,kk,sigin)
          if (all(dum==0.)) dum=1.5E-4
          tke(1:ifull,:)=dum
          do k=1,kk
            if (iotest) then
              call histrd4s(ncid,iarchi,ier,'eps',ik,6*ik,k,
     &                      v_k(:,k),ifull)
            else
              call histrd4s(ncid,iarchi,ier,'eps',ik,6*ik,k,
     &                      vcc,6*ik*ik)
              call doints4(vcc,v_k(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
            end if ! iotest
          end do
          call vertint(v_k,dum,5,kk,sigin)
          if (all(dum==0.)) dum=1.E-7
          eps(1:ifull,:)=dum
        end if

        !------------------------------------------------------------
        ! Tracer data
        if (ngas>0) then              
          do igas=1,ngas              
            write(trnum,'(i3.3)') igas
            do k=1,kk
              if (iotest) then
                call histrd4s(ncid,iarchi,ier,'tr'//trnum,ik,6*ik,k,
     &                   u_k(:,k),ifull)
              else
                ucc=0.
                call histrd4s(ncid,iarchi,ier,'tr'//trnum,ik,6*ik,k,
     &                        ucc,6*ik*ik)
                call doints4(ucc,u_k(:,k),nface4,xg4,yg4,
     &                       nord,dk,ifg)              
              end if ! iotest
            end do
            call vertint(u_k,dum,7,kk,sigin)
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
            do k=1,kk
              if (iotest) then
                call histrd4s(ncid,iarchi,ier,vname,ik,6*ik,k,
     &                        u_k(:,k),ifull)
              else
                call histrd4s(ncid,iarchi,ier,vname,ik,6*ik,k,
     &                          ucc,6*ik*ik)
                call doints4(ucc,u_k(:,k),nface4,xg4,yg4,
     &                       nord,dk,ifg)
              end if ! iotest
            end do
            call vertint(u_k,dum,5,kk,sigin)
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
              so4t(:)=so4t(:)+3.e3*xtg(:,k,3)*(-psl(:)*dsig(k))/grav
            enddo
          end if
        end if

        if (nested==0) then
          ! GEOPOTENTIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          phi_nh=0.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'zgnhs',ik,6*ik,k,
     &                     phi_nh(:,k),ifull)
            enddo  ! k loop
          end if

          ! SDOT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          sdot=-999.
          if (lrestart) then
            sdot(:,1)=0.
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'sdot',ik,6*ik,k,
     &                     sdot(:,k+1),ifull)
            enddo  ! k loop
          end if

          ! PSLX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          pslx=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'pslx',ik,6*ik,k,
     &                     pslx(1:ifull,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savu=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savu',ik,6*ik,k,
     &                     savu(:,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savv=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savv',ik,6*ik,k,
     &                     savv(:,k),ifull)
            enddo  ! k loop
          end if

          ! SAVU1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savu1=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savu1',ik,6*ik,k,
     &                     savu1(:,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVV1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savv1=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savv1',ik,6*ik,k,
     &                     savv1(:,k),ifull)
            enddo  ! k loop
          end if

          ! SAVU2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! only for restart - no interpolation
          savu2=-999.
          if (lrestart) then
            do k=1,kk 
             call histrd4s(ncid,iarchi,ier,'savu2',ik,6*ik,k,
     &                     savu2(:,k),ifull)
            enddo  ! k loop
          end if
          
          ! SAVV2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        ! SOIL ICE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,ms
          write(vname,'("wbice",I1.1)') k
          if (iotest) then
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   wbice(:,k),ifull)
          else
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   ucc,6*ik*ik)
            if (myid==0) then
              where (.not.land_a(:))
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,wbice(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
          end if ! iotest
        end do

        if (nmlo==0.or.abs(nmlo)>9) then ! otherwise already read above
         do k=1,3
          write(vname,'("tggsn",I1.1)') k
          if (iotest) then
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,tggsn(:,k),
     &                   ifull)
          else
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,ucc,
     &                   6*ik*ik)
            if(myid==0)then
              where (.not.land_a)
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            endif  ! (myid==0)
            call doints4(ucc,tggsn(:,k),nface4,xg4,yg4,nord
     &                 ,dk,ifg)
          end if
          if (all(tggsn(:,k)==0.)) tggsn(:,k)=280.
          where(.not.land)
            tggsn(:,k)=280.
          end where
         end do
        end if

        do k=1,3
          write(vname,'("smass",I1.1)') k
          if (iotest) then
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   smass(:,k),ifull)
          else
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   ucc,6*ik*ik)
            if (myid==0) then
              where (.not.land_a(:))
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,smass(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
          end if ! iotest
        end do
        do k=1,3
          write(vname,'("ssdn",I1.1)') k
          if (iotest) then
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     ssdn(:,k),ifull)
          else
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   ucc,6*ik*ik)
            if (myid==0) then
              where (.not.land_a(:))
                ucc=spval
              end where
              call fill_cc(ucc,spval,ik,0)
            end if
            call doints4(ucc,ssdn(:,k),nface4,xg4,yg4,
     &                     nord,dk,ifg)
          end if ! iotest
          if (all(ssdn(:,k)==0.)) then
            where (snowd>100.)
              ssdn(:,k)=240.
            elsewhere
              ssdn(:,k)=140.
            end where
          end if
        end do
        ssdnn=ssdn(:,1)
        
        if (iotest) then
          call histrd1(ncid,iarchi,ier,'snage',ik,6*ik,snage,
     &                 ifull)
        else
          call histrd1(ncid,iarchi,ier,'snage',ik,6*ik,ucc,
     &                 6*ik*ik)
          if (myid==0) then
            where (.not.land_a)
              ucc=spval
            end where
            call fill_cc(ucc,spval,ik,0)
          end if
          call doints4(ucc,snage,nface4,xg4,yg4,
     &                   nord,dk,ifg)
        end if ! iotest

        if (iotest) then
          call histrd1(ncid,iarchi,ier,'sflag',ik,6*ik,dum6,
     &                 ifull)
        else
          call histrd1(ncid,iarchi,ier,'sflag',ik,6*ik,ucc,
     &                 6*ik*ik)
          if (myid==0) then
            where (.not.land_a)
              ucc=spval
            end where
            call fill_cc(ucc,spval,ik,0)
          end if
          call doints4(ucc,dum6,nface4,xg4,yg4,
     &                   nord,dk,ifg)
        end if ! iotest
        isflag=nint(dum6)
        
      endif    ! (nested/=1)

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
      qg(1:ifull,1:kl) = max(qg(1:ifull,1:kl),1.e-6)

      if (myid==0.and.nested==0) then
        write(6,*) "Final lrestart ",lrestart
      end if

      return
      end subroutine ontheflyx

      subroutine doints4(s_a,sout,nface4 ,xg4 ,yg4,nord,ik,ifg)  ! does calls to intsb
      
      use cc_mpi           ! CC MPI routines
      
      implicit none
      
      include 'newmpar.h'  ! Grid parameters

      integer, intent(in) :: ik, nord, ifg      
      integer, intent(in), dimension(ifg,4) :: nface4
      real, dimension(ifull), intent(inout) :: sout
      real, dimension(ik*ik*6), intent(in) :: s_a
      real, intent(in), dimension(ifg,4) :: xg4, yg4
      real, dimension(ifg) :: s_g

      if ( myid ==0 ) then
         if (ifg/=ifull_g) then
           write(6,*) "ERROR: Incorrect ifg in doints4"
           stop
         end if
         if (ik==0) then
           write(6,*) "ERROR: Incorrect ik in doints4"
           stop
         end if
         call ints4(s_a,s_g,nface4 ,xg4 ,yg4,nord,ik)
         call ccmpi_distribute(sout,s_g)
      else
         call ccmpi_distribute(sout)
      endif
      end subroutine doints4

      subroutine ints4(s,sout,nface4 ,xg4 ,yg4,nord,ik)  ! does calls to intsb
      
      use cc_mpi, only : mydiag  ! CC MPI routines
      
      implicit none
      
      include 'newmpar.h'        ! Grid parameters
      include 'parm.h'           ! Model configuration
      
      integer, parameter :: ntest=0
      real, dimension(ik*ik*6), intent(inout) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ifull_g,4) :: nface4
      real, intent(in), dimension(ifull_g,4) :: xg4, yg4
      integer, intent(in) :: ik, nord
      real wrk(ifull_g,4)
      integer :: iq, mm
      integer :: idx = 25, jdx = 218, idjdx=10441

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
      
      include 'newmpar.h'  ! Grid parameters
      include 'parm.h'     ! Model configuration
      
      real, dimension(ik*ik*6), intent(in) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ifull_g) :: nface
      integer idel, jdel, ik, nn
      integer :: ind, i, j, n, iq, n_n, n_e, n_w, n_s
      real aaa, c1, c2, c3, c4, xxg, yyg
      real, intent(in), dimension(ifull_g) :: xg, yg
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
       sx(ik+1,-1,n)=s(ind(2,1,n_e))      ! ess  
      enddo  ! n loop
c     for ew interpolation, sometimes need (different from ns):
c          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

      do iq=1,ifull_g   ! runs through list of target points
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

      implicit none
      
      include 'newmpar.h'  ! Grid parameters
      include 'parm.h'     ! Model configuration
      
      real, dimension(ik*ik*6), intent(inout) :: s
      real, dimension(ifull_g), intent(inout) :: sout
      integer, intent(in), dimension(ifull_g) :: nface
      real, intent(in), dimension(ifull_g) :: xg, yg
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
       sout(iq)=yyg*(xxg*sx(idel+1,jdel+1,n)
     .               +(1.-xxg)*sx(idel,jdel+1,n))
     .    +(1.-yyg)*(xxg*sx(idel+1,jdel,n)
     .               +(1.-xxg)*sx(idel,jdel,n))
      enddo    ! iq loop

      end subroutine ints_blb

      subroutine fill_cc(a_io,value,ik,ndiag)
      
!     this version holds whole array in memory      
c     routine fills in interior of an array which has undefined points

      use cc_mpi          ! CC MPI routines
      
      implicit none
      
      include 'newmpar.h' ! Grid parameters

      real a_io(ik*ik*6)         ! input and output array
      real value            ! array value denoting undefined
c     real b(ik*ik*6), a(ik*ik*6+iextra)
      real b(ik*ik*6), a(ik*ik*6)
      integer :: nrem, i, ii, ik, iq, ind, j, n, neighb, ndiag
      real :: av     
      integer, dimension(ik*ik*6) :: in,ie,iw,is
      integer npann(0:5),npane(0:5),npanw(0:5),npans(0:5)
      data npann/1,103,3,105,5,101/,npane/102,2,104,4,100,0/
      data npanw/5,105,1,101,3,103/,npans/104,0,100,2,102,4/
      ind(i,j,n)=i+(j-1)*ik+n*ik*ik  ! *** for n=0,npanels
      
      if (all(a_io==value)) return
      
       do iq=1,ik*ik*6
       in(iq)=iq+ik
       is(iq)=iq-ik
       ie(iq)=iq+1
       iw(iq)=iq-1
      enddo   ! iq loop
      do n=0,npanels
      if(npann(n)<100)then
        do ii=1,ik
         in(ind(ii,ik,n))=ind(ii,1,npann(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         in(ind(ii,ik,n))=ind(1,ik+1-ii,npann(n)-100)
        enddo    ! ii loop
      endif      ! (npann(n)<100)
      if(npane(n)<100)then
        do ii=1,ik
         ie(ind(ik,ii,n))=ind(1,ii,npane(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         ie(ind(ik,ii,n))=ind(ik+1-ii,1,npane(n)-100)
        enddo    ! ii loop
      endif      ! (npane(n)<100)
      if(npanw(n)<100)then
        do ii=1,ik
         iw(ind(1,ii,n))=ind(ik,ii,npanw(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         iw(ind(1,ii,n))=ind(ik+1-ii,ik,npanw(n)-100)
        enddo    ! ii loop
      endif      ! (npanw(n)<100)
      if(npans(n)<100)then
        do ii=1,ik
         is(ind(ii,1,n))=ind(ii,ik,npans(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         is(ind(ii,1,n))=ind(ik,ik+1-ii,npans(n)-100)
        enddo    ! ii loop
      endif      ! (npans(n)<100)
      enddo      ! n loop
          
      a(1:ik*ik*6) = a_io(:)
      nrem = 1    ! Just for first iteration
c     nrem_gmin = 1 ! Just for first iteration
!     nrem_gmin used to avoid infinite loops, e.g. for no sice
      do while ( nrem > 0)
         ! This has to loop until all are finished otherwise the bounds call
         ! doesn't work.
c808         call bounds(a)
         nrem=0
         do iq=1,ik*ik*6
            b(iq)=a(iq)
            if(a(iq)==value)then
               neighb=0
               av=0.
               if(a(in(iq))/=value)then
                  neighb=neighb+1
                  av=av+a(in(iq))
               endif
               if(a(ie(iq))/=value)then
                  neighb=neighb+1
                  av=av+a(ie(iq))
               endif
               if(a(iw(iq))/=value)then
                  neighb=neighb+1
                  av=av+a(iw(iq))
               endif
               if(a(is(iq))/=value)then
                  neighb=neighb+1
                  av=av+a(is(iq))
               endif
               if(neighb>0)then
                  b(iq)=av/neighb
               else
                  nrem=nrem+1   ! current number of points without a neighbour
               endif
            endif
         end do
         do iq=1,ik*ik*6
            a(iq)=b(iq)
         enddo
      end do
      a_io(1:ik*ik*6) = a(1:ik*ik*6)
      return
      end

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
      real pmsl(ifullx),psl(ifullx),zs(ifullx),t(ifullx)
      real c, con, conr, dlnps, phi1, tav, tsurf
      c=grav/stdlapse
      conr=c/rdry
      con=siglev**(rdry/c)/c
        do iq=1,ifullx
         phi1=t(iq)*rdry*(1.-siglev)/siglev ! phi of sig(lev) above sfce
         tsurf=t(iq)+phi1*stdlapse/grav
         tav=tsurf+max(0.,zs(iq))*.5*stdlapse/grav
         dlnps=max(0.,zs(iq))/(rdry*tav)
         pmsl(iq)=1.e5*exp(psl(iq)+dlnps)
        enddo
      return
      end
      
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
      
      integer ifullx,iq
      real pmsl(ifullx),psl(ifullx),zs(ifullx),t(ifullx)
      integer :: lev
      real dlnps, phi1, tav, tsurf
      do iq=1,ifullx
       phi1=t(iq)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
       tsurf=t(iq)+phi1*stdlapse/grav
       tav=tsurf+max(0.,zs(iq))*.5*stdlapse/grav
       dlnps=max(0.,zs(iq))/(rdry*tav)
       psl(iq)=log(1.e-5*pmsl(iq)) -dlnps
      enddo
      if(nmaxpr==1.and.mydiag)then
        write(6,*)'to_psl lev,sig(lev) ',lev,sig(lev)
        write(6,*)'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd),psl(idjd),pmsl(idjd)
      endif
      return
      end
      
      subroutine interpwind(ik,uct_g,vct_g,ucc,vcc,axs_a,ays_a,azs_a,
     &                      bxs_a,bys_a,bzs_a,rotpole,rotpoles,nface4,
     &                      xg4,yg4,nord)
      
      use vecsuv_m         ! Map to cartesian coordinates
      
      implicit none
      
      include 'newmpar.h'  ! Grid parameters
      
      integer, intent(in) :: ik,nord
      integer, dimension(ifull_g), intent(in) :: nface4
      integer iq
      real, dimension(3,3), intent(in) :: rotpole,rotpoles
      real, dimension(ifull_g), intent(in) :: xg4,yg4
      real, dimension(6*ik*ik), intent(in) :: axs_a,ays_a,azs_a
      real, dimension(6*ik*ik), intent(in) :: bxs_a,bys_a,bzs_a
      real, dimension(6*ik*ik), intent(inout) :: ucc,vcc
      real, dimension(6*ik*ik) :: wcc
      real, dimension(ifull_g), intent(out) :: uct_g,vct_g
      real, dimension(ifull_g) :: wct_g
      real uc,vc,wc,uct_gg,vct_gg,wct_gg
      
      do iq=1,ik*ik*6
        ! first set up winds in Cartesian "source" coords            
        uc=axs_a(iq)*ucc(iq) + bxs_a(iq)*vcc(iq)
        vc=ays_a(iq)*ucc(iq) + bys_a(iq)*vcc(iq)
        wc=azs_a(iq)*ucc(iq) + bzs_a(iq)*vcc(iq)
        ! now convert to winds in "absolute" Cartesian components
        ucc(iq)=uc*rotpoles(1,1)+vc*rotpoles(1,2)+wc*rotpoles(1,3)
        vcc(iq)=uc*rotpoles(2,1)+vc*rotpoles(2,2)+wc*rotpoles(2,3)
        wcc(iq)=uc*rotpoles(3,1)+vc*rotpoles(3,2)+wc*rotpoles(3,3)
      end do
      ! interpolate all required arrays to new C-C positions
      ! don't need to do map factors and Coriolis on target grid
      call ints4(ucc,  uct_g, nface4,xg4,yg4,nord,ik)
      call ints4(vcc,  vct_g, nface4,xg4,yg4,nord,ik)
      call ints4(wcc,  wct_g, nface4,xg4,yg4,nord,ik)
      do iq=1,ifull_g
        ! now convert to "target" Cartesian components (transpose used)
        uct_gg=uct_g(iq)*rotpole(1,1)+vct_g(iq)*rotpole(2,1)
     &                          +wct_g(iq)*rotpole(3,1)
        vct_gg=uct_g(iq)*rotpole(1,2)+vct_g(iq)*rotpole(2,2)
     &                          +wct_g(iq)*rotpole(3,2)
        wct_gg=uct_g(iq)*rotpole(1,3)+vct_g(iq)*rotpole(2,3)
     &                          +wct_g(iq)*rotpole(3,3)
        ! then finally to "target" local x-y components
        uct_g(iq) = ax_g(iq)*uct_gg + ay_g(iq)*vct_gg +
     &                  az_g(iq)*wct_gg
        vct_g(iq) = bx_g(iq)*uct_gg + by_g(iq)*vct_gg +
     &                  bz_g(iq)*wct_gg
      enddo               ! iq loop
      
      return
      end subroutine interpwind
