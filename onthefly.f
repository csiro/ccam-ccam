      ! Main netcdf input routines.  Host grid is automatically
      ! interpolated to nested model grid.  Three options are
      !   nested=0  Initial conditions
      !   nested=1  Nudging fields
      !   nested=2  Surface data recycling
       
      subroutine onthefly(nested,kdate_r,ktime_r,
     .                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     .                    tgg,wb,wbice,snowd,qfg,qlg,
     .                    tggsn,smass,ssdn,ssdnn,snage,isflag,
     .                    iaero,mlodwn,ocndwn)

      use cc_mpi           ! CC MPI routines
      use infile           ! Input file routines
      use mlo              ! Ocean physics and prognostic arrays
      use soil_m           ! Soil and surface data

      implicit none

      include 'newmpar.h'  ! Grid parameters
      include 'darcdf.h'   ! Netcdf data
      include 'mpif.h'     ! MPI parameters
      include 'netcdf.inc' ! Netcdf parameters
      include 'parm.h'     ! Model configuration
      include 'stime.h'    ! File date data

      integer, parameter :: ntest=0
      integer, parameter :: nord=3    ! 1 for bilinear, 3 for bicubic
      integer, parameter :: nihead=54
      integer, parameter :: nrhead=14

      integer, save :: ik,jk,kk,maxarchi
      integer, save :: ncidold=-1
      integer, dimension(nihead) :: nahead
      integer, dimension(ifull) :: isflag
      integer kdate_r,ktime_r,nested,ier,ier2,ilen,itype,iaero
      integer idv,mtimer,k
      real, dimension(nrhead) :: ahead
      real, save :: rlong0x, rlat0x, schmidtx
!     These are local arrays, not the versions in arrays.h
!     Use in call to infile, so are dimensioned ifull rather than ifull_g
      real psl(ifull),zss(ifull),tss(ifull),fracice(ifull),
     & wb(ifull,ms),wbice(ifull,ms),snowd(ifull),sicedep(ifull),
     & t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     & tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     & ssdnn(ifull),snage(ifull),qfg(ifull,kl),qlg(ifull,kl),
     & mlodwn(ifull,wlev,4),ocndwn(ifull,2)
      real timer
      logical ltest,newfile

      !--------------------------------------------------------------
      if ( myid==0 )then
        write(6,*) 'Entering onthefly for nested,ktau = ',
     &                                    nested,ktau
        ! turn OFF fatal netcdf errors; from here on
        call ncpopt(0)
        if(ncid.ne.ncidold)then
          write(6,*) 'Reading new file metadata'
          iarchi=1
          call ncainq(ncid,ncglobal,'int_header',itype,ilen,ier)
          call ncagt(ncid,ncglobal,'int_header',nahead,ier)
          call ncainq(ncid,ncglobal,'real_header',itype,ilen,ier)
          call ncagt(ncid,ncglobal,'real_header',ahead,ier)
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
          idv = ncdid(ncid,'time',ier)
          ier = nf_inq_dimlen(ncid,idv,maxarchi)
          write(6,*) "Found ik,jk,kk,maxarchi = ",ik,jk,kk,maxarchi
          write(6,*) "      rlong0x,rlat0x,schmidtx ",
     &                      rlong0x,rlat0x,schmidtx
        end if
        write(6,*)'Search for kdate_s,ktime_s >= ',
     &                        kdate_s,ktime_s
        ltest=.true.
        iarchi=iarchi-1
        do while(ltest.and.iarchi.lt.maxarchi)
          iarchi=iarchi+1
          idv = ncvid(ncid,'kdate',ier)
          call ncvgt1(ncid,idv,iarchi,kdate_r,ier)
          idv = ncvid(ncid,'ktime',ier)
          call ncvgt1(ncid,idv,iarchi,ktime_r,ier)
          idv = ncvid(ncid,'mtimer',ier)
          if (ier.eq.0) then
            call ncvgt1(ncid,idv,iarchi,mtimer,ier)
            timer=mtimer/60.
          else
            timer=0.
            idv = ncvid(ncid,'timer',ier)
            call ncvgt1(ncid,idv,iarchi,timer,ier)
            mtimer=nint(timer*60.)
          endif
          if (mtimer>0) then
            call datefix(kdate_r,ktime_r,mtimer)
          end if
          ltest=2400*(kdate_r-kdate_s)-1200*nsemble
     &              +(ktime_r-ktime_s)<0
        end do
        if (nsemble.ne.0) then
          kdate_r=kdate_s
          ktime_r=ktime_s
        end if
        if (ltest) then
          ktime_r=-1
        end if
        write(6,*) 'After search ltest,iarchi =',ltest,iarchi
        write(6,*) '             kdate_r,ktime_r =',kdate_r,ktime_r
      endif  ! ( myid==0 )

      call MPI_Bcast(kdate_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_Bcast(ktime_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_Bcast(ncid   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      newfile=ncid.ne.ncidold
      if (newfile) then
        call MPI_Bcast(ik   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        call MPI_Bcast(jk   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        call MPI_Bcast(kk   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        call MPI_Bcast(rlong0x ,1,MPI_REAL,0,MPI_COMM_WORLD,ier)
        call MPI_Bcast(rlat0x  ,1,MPI_REAL,0,MPI_COMM_WORLD,ier)
        call MPI_Bcast(schmidtx,1,MPI_REAL,0,MPI_COMM_WORLD,ier)
        ncidold=ncid
      end if

      if (ktime_r.lt.0) then
        if (nested==2) then
          write(6,*) "WARN: Cannot locate date/time in input file"
          return
        end if
        write(6,*) "ERROR: Cannot locate date/time in input file"
        call MPI_Abort(MPI_COMM_WORLD,-1,ier)
      end if
      !--------------------------------------------------------------
      
      ! Here we call ontheflyx with different automatic array
      ! sizes.  This means the arrays are correct for interpolation
      ! and file i/o on myid==0, as well as the arrays are smaller
      ! on myid.ne.0 when they are not needed.  The code is still
      ! human readable since there is only one ontheflyx subroutine.
      if (myid==0) then
        call ontheflyx(nested,kdate_r,ktime_r,
     &                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &                    tgg,wb,wbice,snowd,qfg,qlg,
     &                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     &                    ik,iaero,mlodwn,ocndwn,rlong0x,rlat0x,
     &                    schmidtx,newfile) ! ik controls automatic array size
        write(6,*) "Leaving onthefly"
      else
        call ontheflyx(nested,kdate_r,ktime_r,
     &                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &                    tgg,wb,wbice,snowd,qfg,qlg,
     &                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     &                    0,iaero,mlodwn,ocndwn,rlong0x,rlat0x,
     &                    schmidtx,newfile) ! 0 controls automatic array size
      end if

      return
      end
      
      ! Read data from netcdf file
      subroutine ontheflyx(nested,kdate_r,ktime_r,
     &                    psl,zss,tss,sicedep,fracice,t,u,v,qg,
     &                    tgg,wb,wbice,snowd,qfg,qlg,
     &                    tggsn,smass,ssdn,ssdnn,snage,isflag,ik,kk,
     &                    dk,iaero,mlodwn,ocndwn,rlong0x,rlat0x,
     &                    schmidtx,newfile)
      
      use aerosolldr, only : xtg,ssn,naero      ! LDR aerosol scheme
      use ateb, only : atebdwn                  ! Urban
      use carbpools_m                           ! Carbon pools
      use cc_mpi                                ! CC MPI routines
      use cfrac_m                               ! Cloud fraction
      use define_dimensions, only : ncs, ncp    ! CABLE dimensions
      use infile                                ! Input file routines
      use latlong_m                             ! Lat/lon coordinates
      use mlo, only : wlev,micdwn,mloregrid     ! Ocean physics and prognostic arrays
      use mlodynamics, only : watbdy            ! Ocean dynamics
      use morepbl_m                             ! Additional boundary layer diagnostics
      use nsibd_m, only : isoilm                ! Land-surface arrays
      use screen_m                              ! Screen level diagnostics
      use sigs_m                                ! Atmosphere sigma levels
      use soil_m                                ! Soil and surface data
      use tkeeps, only : tke,eps,tkesav,epssav  ! TKE-EPS boundary layer
      use tracers_m                             ! Tracer data
      use utilities                             ! Grid utilities
      use vecsuv_m                              ! Map to cartesian coordinates
      use vvel_m                                ! Additional vertical velocity
      use workglob_m                            ! Additional grid interpolation

      implicit none

      include 'newmpar.h'                       ! Grid parameters
      include 'const_phys.h'                    ! Physical constants
      include 'darcdf.h'                        ! Netcdf data
      include 'mpif.h'                          ! MPI parameters
      include 'netcdf.inc'                      ! Netcdf parameters
      include 'parm.h'                          ! Model configuration
      include 'parmgeom.h'                      ! Coordinate data
      include 'soilv.h'                         ! Soil parameters
      include 'stime.h'                         ! File date data

      integer, parameter :: ntest=0
      integer, parameter :: nord=3        ! 1 for bilinear, 3 for bicubic
      
      integer ik, kk, idv, iaero, isoil
      integer dk ! controls automatic array size
      integer lev,levkk,ier,ierr,igas
      integer kdate_r, ktime_r, nemi, id2,jd2,idjd2,
     &        nested, i, j, k, m, iq, ii, jj, np, numneg

c**   xx4 & yy4 only used in indata & here, so no need to redefine after
c**   onthefly; sometime can get rid of common/bigxy4
      real, dimension(:), allocatable, save :: sigin,zss_a,ocndep_l
      real*8 xx4(1+4*dk,1+4*dk),yy4(1+4*dk,1+4*dk)
      real*8, dimension(dk*dk*6):: z_a,x_a,y_a

!     These are local arrays, not the versions in arrays.h
!     Use in call to infile, so are dimensioned ifull rather than ifull_g
      real, dimension(ifull) :: psl,zss,tss,fracice,
     &                          snowd,sicedep,ssdnn,snage,dum6
      real, dimension(ifull,2) :: ocndwn
      real, dimension(ifull,wlev,4) :: mlodwn
      real, dimension(ifull,ms) :: wb,wbice,tgg
      real, dimension(ifull,3) :: tggsn,smass,ssdn
      real, dimension(ifull,kl) :: t,u,v,qg,qfg,qlg
      real, dimension(ifull,kk) :: t_k,u_k,v_k,qg_k
      integer, dimension(ifull) :: isflag
      real, dimension(dk*dk*6) :: psl_a,tss_a,fracice_a,
     &      snowd_a,sicedep_a,pmsl_a,  tss_l_a,tss_s_a
      real, dimension(dk*dk*6,3) :: tggsn_a
      real, dimension(dk*dk*6) :: t_a,qg_a
      real, dimension(dk*dk*6) :: t_a_lev
      real, parameter :: spval=999. ! missing value flag
      real, intent(in) ::  rlong0x, rlat0x, schmidtx

      ! Used in the global interpolation
      real, dimension(ifull_g,4) :: xg4, yg4
      integer, dimension(ifull_g,4) :: nface4
      real, dimension(ifull_g) :: t_g, qg_g,uct_g, vct_g
      real rotpoles(3,3),rotpole(3,3)
      real, dimension(dk*dk*6) :: ucc, vcc
      real, dimension(ifull) :: tss_l, tss_s, pmsl
      real, dimension(dk*dk*6):: axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a
      real, dimension(dk*dk*6):: wts_a  ! not used here or defined in call setxyz
      logical, dimension(dk*dk*6) :: land_a
      integer, dimension(:), allocatable, save :: isoilm_a
      character*8 vname
      character*3 trnum
      logical iotest,newfile

      ! internal check (should not occur if code is written correctly)
      if (myid==0.and.ik.ne.dk) then
        write(6,*) "ERROR: Incorrect automatic array size in onthefly"
        stop
      end if

      ! land-sea mask method (nemi=3 use soilt, nemi=2 use tgg, nemi=1 use zs)
      nemi=3
      
      ! Determine if interpolation is required
      iotest=6*ik*ik.eq.ifull_g.and.abs(rlong0x-rlong0).lt.1.E-5.and.
     &       abs(rlat0x-rlat0).lt.1.E-5.and.
     &       abs(schmidtx-schmidt).lt.1.E-5
      if (iotest) then
        io_in=1   ! no interpolation
      else
        io_in=-1  ! interpolation
      end if
      if (myid==0) write(6,*) "Interpolation iotest,io_in =",
     &                                       iotest,io_in

      !--------------------------------------------------------------
      ! Determine input grid coordinates and interpolation arrays
      if ( myid==0 ) then
        write(6,*) "Defining input file grid"
        if(m_fly==1)then
          rlong4(:,1)=rlongg_g(:)*180./pi
          rlat4(:,1)=rlatt_g(:)*180./pi
        endif
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
              write(6,9)(i,j,j=1,3),(rotpoles(i,j),j=1,3)
           enddo
        endif                  ! (ktau<3)

!       rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!       rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!       rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
        rotpole = calc_rotpole(rlong0,rlat0)
        if(nmaxpr==1)then   ! already in myid==0 loop
           write(6,*)'in onthefly rotpole:'
           do i=1,3
              write(6,9)(i,j,j=1,3),(rotpole(i,j),j=1,3)
 9            format(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)
           enddo
           write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
           write(6,*)'before latltoij for id,jd: ',id,jd
           if ( nproc==1 ) then
              ! Diagnostics will only be correct if nproc==1
              write(6,*)'rlong4(1-4) ',(rlong4(idjd,m),m=1,4)
              write(6,*)'rlat4(1-4) ',(rlat4(idjd,m),m=1,4)
           end if
           write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,
     &                schmidtx
        endif                  ! (nmaxpr==1)
        do m=1,m_fly  !  was 4, now may be set to 1 in namelist
           do iq=1,ifull_g
              call latltoij(rlong4(iq,m),rlat4(iq,m),         !input
     &                      rlong0x,rlat0x,schmidtx,          !input
     &                      xg4(iq,m),yg4(iq,m),nface4(iq,m), !output (source)
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
          write(6,*)'nface4(1-4) ',(nface4(idjd,m),m=1,4)
          write(6,*)'xg4(1-4) ',(xg4(idjd,m),m=1,4)
          write(6,*)'yg4(1-4) ',(yg4(idjd,m),m=1,4)
          if(nested==0)then
             write(6,"('wb_s(1)#  ',9f7.3)") 
     .           ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
             write(6,"('wb_s(ms)# ',9f7.3)") 
     .           ((wb(ii+(jj-1)*il,ms),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
          endif  ! (nested==0)
        endif
      end if ! (myid==0)

      
      ! special data read for new file
      ! read once when file is first opened
      ! need global zss_a for (potentially) landsea mask and psl interpolation
      ! need global isoilm_a for (potentially) landsea mask
      if (newfile) then
        if (allocated(sigin)) deallocate(sigin)
        allocate(sigin(kk))
        if (myid==0) then
          write(6,*) "Reading fixed fields"
          idv = ncvid(ncid,'lev',ier)
          if (ier.eq.0) call ncvgt(ncid,idv,1,kk,sigin,ier)
          if (ier.ne.0) then
            call ncagt(ncid,ncglobal,'sigma',sigin,ier)
          end if
          if (ier.ne.0) then
            call ncagt(ncid,ncglobal,'sigma_lev',sigin,ier)
          end if
          if (ier.ne.0) then
            idv = ncvid(ncid,'layer',ier)
            if (ier.eq.0) call ncvgt(ncid,idv,1,kk,sigin,ier)
          end if
          write(6,'("sigin=",(9f7.4))') (sigin(k),k=1,kk)
        end if
        call MPI_Bcast(sigin  ,kk,MPI_REAL,0,MPI_COMM_WORLD,ier)
        if (myid==0) then
          if (allocated(zss_a)) deallocate(zss_a)
          if (allocated(isoilm_a)) deallocate(isoilm_a)
          allocate(zss_a(6*ik*ik))
          allocate(isoilm_a(6*ik*ik))
          zss_a=0.
          t_a=-1.
        end if
        call histrd1(ncid,iarchi,ier,'zht',ik,6*ik,zss_a,6*ik*ik)
        call histrd1(ncid,iarchi,ier,'soilt',ik,6*ik,t_a,6*ik*ik)
        if (myid==0) then
          isoilm_a=nint(t_a)
        end if
        if (nmlo.ne.0) then
          if (.not.allocated(ocndep_l)) allocate(ocndep_l(ifull))
          t_a=0.
          call histrd1(ncid,iarchi,ier,'ocndepth',ik,6*ik,t_a,
     &                 6*ik*ik)
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(ocndep_l,t_a)
            else
              call ccmpi_distribute(ocndep_l)
            end if
          else
            call doints4(t_a,ocndep_l,nface4,xg4,yg4,nord,ik)
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

      ! interal errors which should not occur if code is written correctly
      if (.not.allocated(sigin)) then
        write(6,*) "ERROR: sigin is undefined in onthefly"
        stop
      end if
      if (.not.allocated(ocndep_l)) then
        write(6,*) "ERROR: ocndep_l is undefined in onthefly"
        stop
      end if
      if (myid==0) then
        if (.not.allocated(zss_a)) then
          write(6,*) "ERROR: zss_a is undefined in onthefly"
          stop
        end if
        if (.not.allocated(isoilm_a)) then
          write(6,*) "ERROR: isoilm_a is undefined in onthefly"
          stop
        end if
      end if

      ! detemine the level below sig=0.9 (used to calculate psl)
      lev=0
      do while(sig(lev+1).gt.0.9) ! nested grid
        lev=lev+1
      end do
      levkk=0
      do while(sigin(levkk+1).gt.0.9) ! host grid
        levkk=levkk+1
      end do      
      if (myid==0) write(6,*) "Ref height lev,levkk =",
     &                                    lev,levkk

      ! could slightly speed up code by also saving zss arrays
      if (iotest) then
        if (myid==0) then
          call ccmpi_distribute(zss,zss_a)
        else
          call ccmpi_distribute(zss)
        end if
      else
        if ( myid==0 ) then
          call ints4(zss_a,  t_g, nface4,xg4,yg4,nord,ik)
          call ccmpi_distribute(zss,t_g)
        else
          call ccmpi_distribute(zss)
        end if ! myid==0
      end if ! iotest

      !--------------------------------------------------------------
      ! Begin reading host data for current time step
      ! psf read when nested=0 or nested=1.and.nud_p.ne.0
      psl_a=0.
      if (nested==0.or.(nested==1.and.nud_p.ne.0)) then
        call histrd1(ncid,iarchi,ier,'psf',ik,6*ik,psl_a,6*ik*ik)
      endif ! nested.ne.2.and.(nested.ne.1.or.nud_p.ne.0)
      ! tsu always read
      tss_a=293.
      call histrd1(ncid,iarchi,ier,'tsu',ik,6*ik,tss_a,6*ik*ik)
      
      ! set up land-sea mask from either soilt, tss or zss
      if(myid==0)then
        if(nemi==3)then 
          land_a(:)=isoilm_a(:).gt.0
          numneg=count(.not.land_a)
          if (any(isoilm_a(:).lt.0)) nemi=2
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

      !--------------------------------------------------------------
      ! Read ocean data for nudging (sea-ice is read below)
      ! read when nested=0 or nested==1.and.nud.ne.0 or nested=2
      if (nmlo.ne.0) then
        ! fixed ocean depth
        ocndwn(:,1)=ocndep_l
        ! ocean potential temperature
        if (nested.ne.1.or.nud_sst.ne.0) then
          do k=1,wlev
            t_a=max(tss_a,271.)
            if (k.le.ms) then
              write(vname,'("tgg",I1.1)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     t_a,6*ik*ik)
            else
              write(vname,'("tgg",I2.2)') k
              call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     t_a,6*ik*ik)
            end if
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(mlodwn(:,k,1),t_a)
              else
                call ccmpi_distribute(mlodwn(:,k,1))
              end if
            else
              if (myid==0) then
                where (land_a)
                  t_a=spval
                end where
                call fill_cc(t_a,spval,ik,0)
!               interpolate all required arrays to new C-C positions
                call ints4(t_a,   t_g, nface4,xg4,yg4,nord,ik)
                call ccmpi_distribute(mlodwn(:,k,1), t_g)
              else ! myid /= 0
                call ccmpi_distribute(mlodwn(:,k,1))
              endif ! myid==0
            end if ! iotest
          end do
        else
          mlodwn(:,:,1)=293.
        end if ! (nestesd.ne.1.or.nud_sst.ne.0) ..else..
        ! ocean salinity
        if (nested.ne.1.or.nud_sss.ne.0) then
          do k=1,wlev
            qg_a=34.72
            write(vname,'("sal",I2.2)') k
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   qg_a,6*ik*ik)
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(mlodwn(:,k,2),qg_a)     
              else
                call ccmpi_distribute(mlodwn(:,k,2))
              end if
            else
              if (myid==0) then
                where (land_a)
                  qg_a=spval
                end where
                call fill_cc(qg_a,spval,ik,0)
!               interpolate all required arrays to new C-C positions
                call ints4(qg_a, qg_g, nface4,xg4,yg4,nord,ik)
                call ccmpi_distribute(mlodwn(:,k,2),qg_g)              
              else ! myid /= 0
                call ccmpi_distribute(mlodwn(:,k,2))            
              endif ! myid==0
            end if ! iotest
          end do
        else
          mlodwn(:,:,2)=34.72
        end if ! (nestesd.ne.1.or.nud_sss.ne.0) ..else..
        ! ocean currents
        if (nested.ne.1.or.nud_ouv.ne.0) then
          do k=1,wlev
            ucc=0.
            vcc=0.
            write(vname,'("uoc",I2.2)') k
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   ucc,6*ik*ik)
            write(vname,'("voc",I2.2)') k
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   vcc,6*ik*ik)
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(mlodwn(:,k,3),ucc)
                call ccmpi_distribute(mlodwn(:,k,4),vcc)
              else
                call ccmpi_distribute(mlodwn(:,k,3))
                call ccmpi_distribute(mlodwn(:,k,4))
              end if
            else
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
                call ccmpi_distribute(mlodwn(:,k,3), uct_g)
                call ccmpi_distribute(mlodwn(:,k,4), vct_g)
              else ! myid /= 0
                call ccmpi_distribute(mlodwn(:,k,3))
                call ccmpi_distribute(mlodwn(:,k,4))
              endif ! myid==0
            end if ! iotest
          end do
        else
          mlodwn(:,:,3:4)=0.
        end if ! (nestesd.ne.1.or.nud_ouv.ne.0) ..else..
        ! water surface height
        if (nested.ne.1.or.nud_sfh.ne.0) then
          t_a=0.
          call histrd1(ncid,iarchi,ier,'ocheight',ik,6*ik,t_a,
     &                 6*ik*ik)
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(ocndwn(:,2),t_a)
            else
              call ccmpi_distribute(ocndwn(:,2))
            end if
          else
            if (myid==0) then
              where (land_a)
                t_a=spval
              end where
              call fill_cc(t_a,spval,ik,0)
            end if
            call doints4(t_a,ocndwn(:,2),nface4,xg4,yg4,nord,ik)
          end if ! iotest
        else
          ocndwn(:,2)=0
        end if ! (nested.ne.1.or.nud_sfh.ne.0) ..else..
        ! vertically interpolate to new grid
        if (nested.ne.1.or.nud_sst.ne.0.or.nud_sss.ne.0.or.
     &      nud_ouv.ne.0) then
          call mloregrid(ocndwn(:,1),mlodwn)
        end if
      end if
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! read sea ice here for prescribed SSTs configuration and for
      ! mixed-layer-ocean
      ! sea-ice is always read
      sicedep_a=0.
      fracice_a=0.
      call histrd1(ncid,iarchi,ier,'siced',ik,6*ik,sicedep_a,6*ik*ik)
      call histrd1(ncid,iarchi,ierr,'fracice',ik,6*ik,fracice_a,6*ik*ik)
      if (myid==0) then
        if(ier==0)then  ! i.e. sicedep read in 
          if(ierr.ne.0)then ! i.e. sicedep read in; fracice not read in
            where(sicedep_a>0.)
              fracice_a=1.
            endwhere
          endif  ! (ierr.ne.0)  fracice
        else     ! sicedep not read in
          if(ierr.ne.0)then  ! neither sicedep nor fracice read in
            sicedep_a(:)=0.  ! Oct 08
            fracice_a(:)=0.
	      write(6,*)'pre-setting siced in onthefly from tss'
            where(abs(tss_a) <= 271.2)
              sicedep_a=1.  ! Oct 08
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
          endif  ! (ierr.ne.0)
        endif    ! (ier.ne.0) .. else ..    for sicedep
  
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
        if(nproc==1.and.nmaxpr==1)then
          write(6,*)'before fill tss ',tss_a(idjd2)
          write(6,*)'before fill tss_l_a, tss_s_a ',
     &                       tss_l_a(idjd2),tss_s_a(idjd2)
          write(6,*)'before fill/ints4 sicedep ',sicedep_a(idjd2)
        endif  ! (nproc==1.and.nmaxpr==1)
        call fill_cc(tss_l_a,spval,ik,0)
        call fill_cc(tss_s_a,spval,ik,0)
        call fill_cc(sicedep_a,spval,ik,0)
        call fill_cc(fracice_a,spval,ik,0)
        if(nproc==1.and.nmaxpr==1)then
          write(6,*)'after fill tss_l, tss_s ',tss_l_a(idjd2),
     &            tss_s_a(idjd2)
          write(6,*)'after fill sicedep ',sicedep_a(idjd2)
          write(6,*)'before ints4 psl_a(idjd2),zss_a(idjd2) ',
     .                        psl_a(idjd),zss_a(idjd2)
        endif  ! (nproc==1.and.nmaxpr==1)
      endif   ! (myid==0)

      if (iotest) then
        if (myid==0) then
          call ccmpi_distribute(tss,tss_a)
          call ccmpi_distribute(sicedep,sicedep_a)
          call ccmpi_distribute(fracice,fracice_a)
        else
          call ccmpi_distribute(tss)
          call ccmpi_distribute(sicedep)
          call ccmpi_distribute(fracice)
        end if
c       incorporate other target land mask effects
        where (land)
          sicedep=0.
          fracice=0.
        end where
      else
!       The routine doints4 does the gather, calls ints4 and redistributes
        call doints4(tss_l_a , tss_l,  nface4,xg4,yg4,nord,ik)
        call doints4(tss_s_a , tss_s,  nface4,xg4,yg4,nord,ik)
        call doints4(fracice_a , fracice,  nface4,xg4,yg4,nord,ik)
        call doints4(sicedep_a , sicedep,  nface4,xg4,yg4,nord,ik)
c       incorporate other target land mask effects
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
        if(nproc==1.and.nmaxpr==1)then
           write(6,*)'after ints4 idjd,zss(idjd) ',idjd,zss(idjd)
           write(6,*)'zss1-5 ',(zss(iq),iq=1,5)
           write(6,*)'after ints4 psl,pmsl ',psl(idjd),pmsl(idjd)
        endif  ! (nproc==1.and.nmaxpr==1)
      end if ! iotest
      if ( nproc==1 ) write(6,*)'after ints4 sicedep ',sicedep(idjd)


      ! read atmospheric fields for nested=0 or nested=1.and.nud.ne.0

      ! air temperature
      ! read for nested=0 or nested=1.and.(nud_t.ne.0.or.nud_p.ne.0)
      if (nested==0.or.(nested==1.and.(nud_t.ne.0.or.nud_p.ne.0)))
     &  then
        do k=1,kk
          call histrd4s(ncid,iarchi,ier,'temp',ik,6*ik,k,t_a,6*ik*ik) !     temperature
          if (k.eq.levkk) t_a_lev=t_a ! store for psl calculation below
          if (myid==0) then
            if (iotest) then
              call ccmpi_distribute(t_k(:,k), t_a)
            else
              call ints4(t_a,    t_g, nface4,xg4,yg4,nord,ik)  ! ints4 on source grid
              call ccmpi_distribute(t_k(:,k), t_g)
            end if ! iotest
          else ! myid /= 0
            call ccmpi_distribute(t_k(:,k))
          endif ! myid==0
        enddo  ! k loop
        call vertint(t_k ,t(1:ifull,:), 1,kk,sigin)
      end if ! (nested==0.or.(nested==1.and.(nud_t.ne.0.or.nud_p.ne.0)))
      ! winds
      ! read for nested=0 or nested=1.and.nud_uv.ne.0
      if (nested==0.or.(nested==1.and.nud_uv.ne.0)) then
        do k=1,kk
          ! to reduce memory footprint, we now have to alternatively read
          ! u and v.  This is a bit inefficent for disk accessing,
          ! but makes it possible to downscale large grids (e.g., C160)
          call histrd4s(ncid,iarchi,ier,'u',ik,6*ik,k,ucc,6*ik*ik)    !     u wind component
          call histrd4s(ncid,iarchi,ier,'v',ik,6*ik,k,vcc,6*ik*ik)    !     v wind component
          if (myid==0) then
            if (iotest) then
              call ccmpi_distribute(u_k(:,k), ucc)
              call ccmpi_distribute(v_k(:,k), vcc)
            else
              call interpwind(ik,uct_g,vct_g,ucc,vcc,axs_a,ays_a,
     &                        azs_a,bxs_a,bys_a,bzs_a,rotpole,
     &                        rotpoles,nface4,xg4,yg4,nord)
!             interpolate all required arrays to new C-C positions
!             don't need to do map factors and Coriolis on target grid
              call ccmpi_distribute(u_k(:,k), uct_g)
              call ccmpi_distribute(v_k(:,k), vct_g)
            end if ! iotest
          else
            call ccmpi_distribute(u_k(:,k))
            call ccmpi_distribute(v_k(:,k))
          endif ! myid==0
        enddo  ! k loop
        call vertint(u_k ,u(1:ifull,:), 3,kk,sigin)
        call vertint(v_k ,v(1:ifull,:), 4,kk,sigin)
      end if ! (nested==0.or.(nested==1.and.nud_uv.ne.0))
      ! mixing ratio
      ! read for nested=0 or nested=1.and.nud_q.ne.0
      if (nested==0.or.(nested==1.and.nud_q.ne.0)) then
        do k=1,kk
          call histrd4s(ncid,iarchi,ier,'mixr',ik,6*ik,k,qg_a,6*ik*ik)!     mixing ratio
          if(ier.ne.0)then                                            !     mixing ratio
            call histrd4s(ncid,iarchi,ier,'q',ik,6*ik,k,qg_a,6*ik*ik) !     mixing ratio
          endif  ! (ier.ne.0)                                         !     mixing ratio
          if (myid==0) then
            if (iotest) then
              call ccmpi_distribute(qg_k(:,k), qg_a)
            else
              call ints4(qg_a,  qg_g, nface4,xg4,yg4,nord,ik)
              call ccmpi_distribute(qg_k(:,k), qg_g)
            end if ! iotest
          else
            call ccmpi_distribute(qg_k(:,k))         
          endif ! myid==0
        enddo  ! k loop
        call vertint(qg_k,qg(1:ifull,:),2,kk,sigin)
      end if ! (nested==0.or.(nested==1.and.nud_q.ne.0))

      ! re-grid surface pressure by mapping to MSLP, interpolating and then map to surface pressure
      ! requires psl_a, zss, zss_a, t and t_a_lev
      if (nested==0.or.(nested==1.and.nud_p.ne.0)) then
        if (iotest) then
          if (myid==0) then
            call ccmpi_distribute(psl,psl_a)
          else
            call ccmpi_distribute(psl)
          end if
        else
!         The routine doints4 does the gather, calls ints4 and redistributes
          if ( myid==0 ) then
            call mslpx(pmsl_a,psl_a,zss_a,t_a_lev,ik*ik*6,sigin(levkk))  ! needs pmsl (preferred)
            call ints4(pmsl_a, t_g, nface4,xg4,yg4,nord,ik)
            call ccmpi_distribute(pmsl,t_g)
          else
            call ccmpi_distribute(pmsl)
          end if ! myid==0
!         invert pmsl to get psl
          call to_pslx(pmsl,psl,zss,t(:,lev),ifull,lev)  ! on target grid
        end if ! iotest
      end if


      !**************************************************************
      ! This is the end of reading the nudging arrays
      !**************************************************************


      !--------------------------------------------------------------
      ! The following data is only read for initial conditions
      if (nested.ne.1) then

        ! SNOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call histrd1(ncid,iarchi,ier,'snd',ik,6*ik,snowd_a,6*ik*ik)
        if (ier.ne.0) then
          where (tss_a.lt.270.)
            snowd_a=min(55.,5.*(271.-abs(tss_a)))
          end where
          where (snowd_a.lt.0.001)
            snowd_a=0.
          end where
        end if

        ! SOIL TEMPERATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,ms 
          t_a=tss_a
          write(vname,'("tgg",I1.1)') k
          call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   t_a,6*ik*ik)
          if (ier.ne.0) then
            if (k.le.3) then
              call histrd1(ncid,iarchi,ier,'tgg2',ik,6*ik,
     &                 t_a,6*ik*ik)
              if (ier.ne.0) then
                call histrd1(ncid,iarchi,ier,'tb3',ik,6*ik,
     &                 t_a,6*ik*ik)
              end if
            else
              call histrd1(ncid,iarchi,ier,'tgg6',ik,6*ik,
     &                 t_a,6*ik*ik)
              if (ier.ne.0) then
                call histrd1(ncid,iarchi,ier,'tb2',ik,6*ik,
     &                 t_a,6*ik*ik)
              end if
            end if
          end if
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(tgg(:,k),t_a)
            else
              call ccmpi_distribute(tgg(:,k))
            end if
          else
            if (myid==0) then
              where (.not.land_a(:))
                t_a=spval
              end where
              call fill_cc(t_a,spval,ik,0)
            end if
            call doints4(t_a,tgg(:,k),nface4,xg4,yg4,
     &                     nord,ik)
          end if
        end do

        !--------------------------------------------------
        ! Read MLO sea-ice data
        if (nmlo.ne.0) then
          if (.not.allocated(micdwn)) allocate(micdwn(ifull,10))
          do k=1,8
            select case(k)
              case(1,2,3)
                t_a=280.
                write(vname,'("tggsn",I1.1)') k
                call histrd1(ncid,iarchi,ier,vname,ik,6*ik,t_a,
     &                 6*ik*ik)
                tggsn_a(:,k)=t_a
              case(4)
                 t_a=272.2
                 call histrd1(ncid,iarchi,ier,'tggsn4',ik,6*ik,t_a,
     &                 6*ik*ik)
              case(5)
                t_a=fracice_a ! read above with nudging arrays
              case(6)
                t_a=sicedep_a ! read above with nudging arrays
              case(7)
                t_a=snowd_a/1000.
              case(8)
                t_a=0.
                call histrd1(ncid,iarchi,ier,'sto',ik,6*ik,t_a,
     &                 6*ik*ik)
            end select
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(micdwn(:,k),t_a)
              else
                call ccmpi_distribute(micdwn(:,k))
              end if
            else
              if (myid==0) then
                where (land_a)
                  t_a=spval
                end where
                call fill_cc(t_a,spval,ik,0)
              end if
              call doints4(t_a,micdwn(:,k),nface4,xg4,yg4,nord,ik)
            end if
          end do
          ucc=0.
          vcc=0.
          call histrd1(ncid,iarchi,ier,'uic',ik,6*ik,ucc,6*ik*ik)
          call histrd1(ncid,iarchi,ier,'vic',ik,6*ik,vcc,6*ik*ik)
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(micdwn(:,9),ucc)
              call ccmpi_distribute(micdwn(:,10),vcc)
            else
              call ccmpi_distribute(micdwn(:,9))
              call ccmpi_distribute(micdwn(:,10))
            end if
          else
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
          if (.not.allocated(watbdy)) then
            allocate(watbdy(ifull+iextra))
            watbdy=0.
          end if
          if (abs(nmlo).ge.2) then
            t_a=0.
            call histrd1(ncid,iarchi,ier,'swater',ik,6*ik,t_a,
     &                   6*ik*ik)
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(watbdy(1:ifull),t_a)
              else
                call ccmpi_distribute(watbdy(1:ifull))
              end if
            else
              call doints4(t_a,watbdy(1:ifull),nface4,xg4,yg4,nord,ik)
            end if
          end if
        end if
        !--------------------------------------------------

        ! SOIL MOISTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,ms
          t_a=20.5 
          write(vname,'("wetfrac",I1.1)') k
          call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   t_a,6*ik*ik)
          if (ier.eq.0) then
            t_a=t_a+20.
          else
            write(vname,'("wb",I1.1)') k
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                     t_a,6*ik*ik)
            if (ier.ne.0) then
              if (k.le.2) then
                call histrd1(ncid,iarchi,ier,'wb2',ik,6*ik,t_a,
     &                 6*ik*ik)
                if (ier.ne.0) then
                  call histrd1(ncid,iarchi,ier,'wfg',ik,6*ik,t_a,
     &                   6*ik*ik)
                end if
              else
                call histrd1(ncid,iarchi,ier,'wb6',ik,6*ik,t_a,
     &                 6*ik*ik)
                if (ier.ne.0) then
                  call histrd1(ncid,iarchi,ier,'wfb',ik,6*ik,t_a,
     &                   6*ik*ik)
                end if
              end if
            end if
          end if
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(wb(:,k),t_a)
            else
              call ccmpi_distribute(wb(:,k))
            end if
          else
            if (myid==0) then
              where (.not.land_a(:))
                t_a=spval
              end where
              call fill_cc(t_a,spval,ik,0)
            end if
            call doints4(t_a,wb(:,k),nface4,xg4,yg4,
     &                     nord,ik)
          end if ! iotest
        end do
        !unpack field capacity into volumetric soil moisture
        if (any(wb(:,:).gt.10.)) then
          if (mydiag) write(6,*) "Unpacking wetfrac to wb",wb(idjd,1)
          wb(:,:)=max(wb(:,:)-20.,0.)
          do iq=1,ifull
            isoil=isoilm(iq)
            wb(iq,:)=(1.-wb(iq,:))*swilt(isoil)+wb(iq,:)*sfc(isoil)
          end do
          if (mydiag) write(6,*) "giving wb",wb(idjd,1)
        end if  

        !--------------------------------------------------
        ! Read 10m wind speeds for special sea roughness length calculations
        call histrd1(ncid,iarchi,ier,'u10',ik,6*ik,t_a,6*ik*ik)
        if (ier==0) then
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(u10,t_a)
            else
              call ccmpi_distribute(u10)
            end if
          else
            call doints4(t_a,u10,nface4,xg4,yg4,nord,ik)
          end if ! iotest
        else
          u10=sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)*log(10./0.001)
     &                                            /log(zmin/0.001)
        end if

        !--------------------------------------------------
        ! Read boundary layer height for TKE-eps mixing
        if (nvmix.eq.6) then
          t_a=1000. ! dummy for pbl
          call histrd1(ncid,iarchi,ier,'pblh',ik,6*ik,t_a,6*ik*ik)
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(pblh,t_a)
            else
              call ccmpi_distribute(pblh)
            end if
          else
            call doints4(t_a,pblh,nface4,xg4,yg4,
     &                     nord,ik)
          end if ! iotest
        end if

        !--------------------------------------------------
        ! Read CABLE aggregate carbon pools
        if (nsib.eq.4.or.nsib.eq.6.or.nsib.eq.7) then
          do k=1,ncp
            t_a=0.
            write(vname,'("cplant",I1.1)') k
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   t_a,6*ik*ik)
            if (ier.eq.0) then
              if (iotest) then
                if (myid==0) then
                  call ccmpi_distribute(cplant(:,k),t_a)
                else
                  call ccmpi_distribute(cplant(:,k))
                end if
              else
                if (myid==0) then
                  where (.not.land_a(:))
                    t_a=spval
                  end where
                  call fill_cc(t_a,spval,ik,0)
                end if
                call doints4(t_a,cplant(:,k),nface4,xg4,yg4,
     &                     nord,ik)
              end if ! iotest
            end if
          end do
          do k=1,ncs
            t_a=0.
            write(vname,'("csoil",I1.1)') k
            call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   t_a,6*ik*ik)
            if (ier.eq.0) then
              if (iotest) then
                if (myid==0) then
                  call ccmpi_distribute(csoil(:,k),t_a)
                else
                  call ccmpi_distribute(csoil(:,k))
                end if
              else
                if (myid==0) then
                  where (.not.land_a(:))
                    t_a=spval
                  end where
                  call fill_cc(t_a,spval,ik,0)
                end if
                call doints4(t_a,csoil(:,k),nface4,xg4,yg4,
     &                       nord,ik)
              end if ! iotest
            end if
          end do
        end if

        !--------------------------------------------------
        ! Read urban data
        if (nurban.ne.0) then
          if (.not.allocated(atebdwn)) allocate(atebdwn(ifull,22))
          do k=1,22
            t_a=999.
            select case(k)
              case(1)
                call histrd1(ncid,iarchi,ier,'rooftgg1',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(2)
                call histrd1(ncid,iarchi,ier,'rooftgg2',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(3)
                call histrd1(ncid,iarchi,ier,'rooftgg3',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(4)
                call histrd1(ncid,iarchi,ier,'waletgg1',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(5)
                call histrd1(ncid,iarchi,ier,'waletgg2',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(6)
                call histrd1(ncid,iarchi,ier,'waletgg3',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(7)
                call histrd1(ncid,iarchi,ier,'walwtgg1',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(8)
                call histrd1(ncid,iarchi,ier,'walwtgg2',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(9)
                call histrd1(ncid,iarchi,ier,'walwtgg3',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(10)
                call histrd1(ncid,iarchi,ier,'roadtgg1',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(11)
                call histrd1(ncid,iarchi,ier,'roadtgg2',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(12)
                call histrd1(ncid,iarchi,ier,'roadtgg3',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(13)
                call histrd1(ncid,iarchi,ier,'urbansm',ik,6*ik,
     &                       t_a,6*ik*ik)
              case(14)
                call histrd1(ncid,iarchi,ier,'roofwtr',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=0.
              case(15)
                call histrd1(ncid,iarchi,ier,'roadwtr',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=0.
              case(16)
                call histrd1(ncid,iarchi,ier,'urblwtr',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=0.
              case(17)
                call histrd1(ncid,iarchi,ier,'roofsnd',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=0.
              case(18)
                call histrd1(ncid,iarchi,ier,'roadsnd',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=0.
              case(19)
                call histrd1(ncid,iarchi,ier,'roofden',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=100.
              case(20)
                call histrd1(ncid,iarchi,ier,'roadden',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=100.
              case(21)
                call histrd1(ncid,iarchi,ier,'roofsna',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=0.85
              case(22)
                call histrd1(ncid,iarchi,ier,'roadsna',ik,6*ik,
     &                       t_a,6*ik*ik)
                if (ier.ne.0) t_a=0.85
            end select
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(atebdwn(:,k),t_a)
              else
                call ccmpi_distribute(atebdwn(:,k))
              end if
            else
              if (myid==0) then
                where (.not.land_a.or.t_a.ge.399.)
                  t_a=spval
                end where
                call fill_cc(t_a,spval,ik,0)
              end if
              call doints4(t_a,atebdwn(:,k),nface4,xg4,yg4,nord
     &                     ,ik)
            end if ! iotest
          end do
        end if
        !--------------------------------------------------
        
        ! CLOUD FROZEN WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,kk 
         ucc=0. ! dummy for qfg
         call histrd4s(ncid,iarchi,ier,'qfg',ik,6*ik,k,ucc,
     &                 6*ik*ik)
         if (iotest) then
           if (myid==0) then
             call ccmpi_distribute(u_k(:,k),ucc)
           else
             call ccmpi_distribute(u_k(:,k))
           end if
         else
           call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,ik)
         end if ! iotest
        enddo  ! k loop
        call vertint(u_k,qfg,5,kk,sigin)
        ! CLOUD LIQUID WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,kk 
         vcc=0. ! dummy for qlg
         call histrd4s(ncid,iarchi,ier,'qlg',ik,6*ik,k,vcc,
     &                 6*ik*ik)
         if (iotest) then
           if (myid==0) then
             call ccmpi_distribute(v_k(:,k),vcc)
           else
             call ccmpi_distribute(v_k(:,k))
           end if
         else
           call doints4(vcc,v_k(:,k),nface4,xg4,yg4,nord,ik)
         end if ! iotest
        enddo  ! k loop
        call vertint(v_k,qlg,5,kk,sigin)

        ! CLOUD FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,kk 
         ucc=0. ! dummy for cfrac
         call histrd4s(ncid,iarchi,ier,'cfrac',ik,6*ik,k,ucc,
     &                 6*ik*ik)
         if (iotest) then
           if (myid==0) then
             call ccmpi_distribute(u_k(:,k),ucc)
           else
             call ccmpi_distribute(u_k(:,k))
           end if
         else
           call doints4(ucc,u_k(:,k),nface4,xg4,yg4,nord,ik)
         end if ! iotest
        enddo  ! k loop
        call vertint(u_k,cfrac,5,kk,sigin)

        !--------------------------------------------------
        ! TKE-eps data
        if (nvmix.eq.6) then
          do k=1,kk
            ucc=1.5E-4 ! dummy for tke
            call histrd4s(ncid,iarchi,ier,'tke',ik,6*ik,k,
     &                    ucc,6*ik*ik)
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(u_k(:,k),ucc)
              else
                call ccmpi_distribute(u_k(:,k))
              end if
            else
              call doints4(ucc,u_k(:,k),nface4,xg4,yg4,
     &                     nord,ik)
            end if ! iotest
          end do
          call vertint(u_k,tke(1:ifull,:),6,kk,sigin)  
          do k=1,kk
            vcc=1.E-7 ! dummy for eps
            call histrd4s(ncid,iarchi,ier,'eps',ik,6*ik,k,
     &                    vcc,6*ik*ik)
            if (iotest) then
              if (myid==0) then
                call ccmpi_distribute(v_k(:,k),vcc)
              else
                call ccmpi_distribute(v_k(:,k))
              end if
            else
              call doints4(vcc,v_k(:,k),nface4,xg4,yg4,
     &                     nord,ik)
            end if ! iotest
          end do
          call vertint(v_k,eps(1:ifull,:),6,kk,sigin)
          tke=max(tke,1.5E-8)
          eps=min(eps,(0.09**0.75)*(tke**1.5)/5.)
          eps=max(eps,(0.09**0.75)*(tke**1.5)/500.)          
          eps=max(eps,1.E-10)       
          tkesav=tke(1:ifull,:)
          epssav=eps(1:ifull,:)
        end if

        !------------------------------------------------------------
        ! Tracer data
        if (ngas>0) then              
          do igas=1,ngas              
            write(trnum,'(i3.3)') igas
            do k=1,kk
              t_a=0.
              call histrd4s(ncid,iarchi,ier,'tr'//trnum,ik,6*ik,k,
     &                 t_a,6*ik*ik)
              if (iotest) then
                if (myid==0) then
                   call ccmpi_distribute(t_k(:,k),t_a)
                else
                  call ccmpi_distribute(t_k(:,k))
                end if
              else
                call doints4(t_a,t_k(:,k),nface4,xg4,yg4,
     &                       nord,ik)              
              end if ! iotest
            end do
            call vertint(t_k,tr(1:ifull,:,igas),7,kk,sigin)
          enddo                       
        endif                         

        !------------------------------------------------------------
        ! Aerosol data
        if (abs(iaero).ge.2) then
          do i=1,naero+2
            do k=1,kk
              ucc=0. ! dummy for aerosol
              select case(i)
                case(1)
                  call histrd4s(ncid,iarchi,ier,'dms',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(2)
                  call histrd4s(ncid,iarchi,ier,'so2',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(3)
                  call histrd4s(ncid,iarchi,ier,'so4',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(4)
                  call histrd4s(ncid,iarchi,ier,'bco',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(5)
                  call histrd4s(ncid,iarchi,ier,'bci',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(6)
                  call histrd4s(ncid,iarchi,ier,'oco',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(7)
                  call histrd4s(ncid,iarchi,ier,'oci',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(8)
                  call histrd4s(ncid,iarchi,ier,'dust1',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(9)
                  call histrd4s(ncid,iarchi,ier,'dust2',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(10)
                  call histrd4s(ncid,iarchi,ier,'dust3',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(11)
                  call histrd4s(ncid,iarchi,ier,'dust4',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(12)
                  call histrd4s(ncid,iarchi,ier,'seasalt1',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case(13)
                  call histrd4s(ncid,iarchi,ier,'seasalt2',ik,6*ik,k,
     &                          ucc,6*ik*ik)
                case default
                  write(6,*) "ERROR: Unknown aerosol type ",i
                  stop
              end select
              if (iotest) then
                if (myid==0) then
                  call ccmpi_distribute(u_k(:,k),ucc)
                else
                  call ccmpi_distribute(u_k(:,k))
                end if
              else
                call doints4(ucc,u_k(:,k),nface4,xg4,yg4,
     &                       nord,ik)
              end if ! iotest
            end do
            if (i.le.naero) then
              call vertint(u_k,xtg(1:ifull,:,i),6,kk,sigin)
            else
              call vertint(u_k,ssn(1:ifull,:,i-naero),6,kk,sigin)
            end if
          end do
          ! Factor 1.e3 to convert to g/m2, x 3 to get sulfate from sulfur
          so4t(:)=0.
          do k=1,kl
            so4t(:)=so4t(:)+3.e3*xtg(:,k,3)*(-psl(:)*dsig(k))/grav
          enddo
        end if

        ! SOIL ICE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,ms
          t_a=0. 
          write(vname,'("wbice",I1.1)') k
          call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   t_a,6*ik*ik)
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(wbice(:,k),t_a)
            else
              call ccmpi_distribute(wbice(:,k))
            end if
          else
            if (myid==0) then
              where (.not.land_a(:))
                t_a=spval
              end where
              call fill_cc(t_a,spval,ik,0)
            end if
            call doints4(t_a,wbice(:,k),nface4,xg4,yg4,
     &                     nord,ik)
          end if ! iotest
        end do

        if (nmlo.eq.0.) then ! otherwise already read above
          tggsn_a(:,:)= 280.
          call histrd1(ncid,iarchi,ier,'tggsn1',ik,6*ik,tggsn_a(:,1),
     &                 6*ik*ik)
          call histrd1(ncid,iarchi,ier,'tggsn2',ik,6*ik,tggsn_a(:,2),
     &                 6*ik*ik)
          call histrd1(ncid,iarchi,ier,'tggsn3',ik,6*ik,tggsn_a(:,3),
     &                 6*ik*ik)
        end if

        if (iotest) then
          if (myid==0) then
            call ccmpi_distribute(snowd,snowd_a)
            do k=1,3
              call ccmpi_distribute(tggsn(:,k),tggsn_a(:,k))
            end do            
          else
            call ccmpi_distribute(snowd)
            do k=1,3
              call ccmpi_distribute(tggsn(:,k))
            end do
          end if
        else
          if(myid==0)then
            do iq=1,ik*ik*6
              if(.not.land_a(iq))then       
                snowd_a(iq)=spval
                tggsn_a(iq,1:3)=spval
              endif  !   (.not.land_a(iq)) 
            enddo   ! iq loop
            call fill_cc(snowd_a,spval,ik,0)
            do k=1,3
              call fill_cc(tggsn_a(:,k),spval,ik,0)
            enddo
          endif  ! (myid==0)
          call doints4(snowd_a,  snowd,nface4,xg4,yg4,nord,ik)
          do k=1,3
            call doints4(tggsn_a(:,k),tggsn(:,k),nface4,xg4,yg4,nord
     &                   ,ik)
          enddo          
          where(.not.land)
            tggsn(:,1)=280.
            tggsn(:,2)=280.
            tggsn(:,3)=280.
          end where
          where (snowd>0.)
            tgg(:,1)=min(tgg(:,1),270.1)
          endwhere
        end if ! iotest

        do k=1,3
          t_a=0.
          write(vname,'("smass",I1.1)') k
          call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   t_a,6*ik*ik)
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(smass(:,k),t_a)
            else
              call ccmpi_distribute(smass(:,k))
            end if
          else
            if (myid==0) then
              where (.not.land_a(:))
                t_a=spval
              end where
              call fill_cc(t_a,spval,ik,0)
            end if
            call doints4(t_a,smass(:,k),nface4,xg4,yg4,
     &                     nord,ik)
          end if ! iotest
        end do
        do k=1,3
          where(snowd_a>100.)
            t_a=240.
          elsewhere
            t_a=140.
          endwhere
          write(vname,'("ssdn",I1.1)') k
          call histrd1(ncid,iarchi,ier,vname,ik,6*ik,
     &                   t_a,6*ik*ik)
          if (iotest) then
            if (myid==0) then
              call ccmpi_distribute(ssdn(:,k),t_a)
            else
              call ccmpi_distribute(ssdn(:,k))
            end if
          else
            if (myid==0) then
              where (.not.land_a(:))
                t_a=spval
              end where
              call fill_cc(t_a,spval,ik,0)
            end if
            call doints4(t_a,ssdn(:,k),nface4,xg4,yg4,
     &                     nord,ik)
          end if ! iotest
        end do
        ssdnn=ssdn(:,1)
        
        t_a=0.
        call histrd1(ncid,iarchi,ier,'snage',ik,6*ik,t_a,
     &               6*ik*ik)
        if (iotest) then
          if (myid==0) then
            call ccmpi_distribute(snage,t_a)
          else
            call ccmpi_distribute(snage)
          end if
        else
          if (myid==0) then
            where (.not.land_a)
              t_a=spval
            end where
            call fill_cc(t_a,spval,ik,0)
          end if
          call doints4(t_a,snage,nface4,xg4,yg4,
     &                   nord,ik)
        end if ! iotest

        t_a=0.
        call histrd1(ncid,iarchi,ier,'sflag',ik,6*ik,t_a,
     &               6*ik*ik)
        if (iotest) then
          if (myid==0) then
            call ccmpi_distribute(dum6,t_a)
          else
            call ccmpi_distribute(dum6)
          end if
        else
          if (myid==0) then
            where (.not.land_a)
              t_a=spval
            end where
            call fill_cc(t_a,spval,ik,0)
          end if
          call doints4(t_a,dum6,nface4,xg4,yg4,
     &                   nord,ik)
        end if ! iotest
        isflag=nint(dum6)
        
      endif    ! (nested.ne.1)

      ! tgg holds file surface temperature when no MLO
      if (nmlo==0) then
        where (.not.land)
          tgg(:,1)=tss
        end where
      end if

      ! set-up for next read of file
      iarchi=iarchi+1
      kdate_s=kdate_r
      ktime_s=ktime_r+1
      qg(1:ifull,1:kl) = max(qg(1:ifull,1:kl),1.e-6)

      return
      end subroutine ontheflyx

      subroutine doints4(s_a,sout,nface4 ,xg4 ,yg4,nord,ik)  ! does calls to intsb
      
      use cc_mpi           ! CC MPI routines
      
      implicit none
      
      include 'newmpar.h'  ! Grid parameters
      
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
         call ccmpi_distribute(sout,s_g)
      else
ccc         call ccmpi_gather(s)
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
         write(6,*)'in ints4 for id,jd,nord: ',id,jd,nord
         write(6,*)'nface4(1-4) ',(nface4(idjd,m),m=1,4)
         write(6,*)'xg4(1-4) ',(xg4(idjd,m),m=1,4)
         write(6,*)'yg4(1-4) ',(yg4(idjd,m),m=1,4)
         write(6,*)'wrk(1-4) ',(wrk(idjd,m),m=1,4)
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
       sx(ik+1,-1,n)=s(ind(2,1,n_e))     ! ess  
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
      include 'mpif.h'    ! MPI parameters

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
      
      if (all(a_io.eq.value)) return
      
       do iq=1,ik*ik*6
       in(iq)=iq+ik
       is(iq)=iq-ik
       ie(iq)=iq+1
       iw(iq)=iq-1
      enddo   ! iq loop
      do n=0,npanels
      if(npann(n).lt.100)then
        do ii=1,ik
         in(ind(ii,ik,n))=ind(ii,1,npann(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         in(ind(ii,ik,n))=ind(1,ik+1-ii,npann(n)-100)
        enddo    ! ii loop
      endif      ! (npann(n).lt.100)
      if(npane(n).lt.100)then
        do ii=1,ik
         ie(ind(ik,ii,n))=ind(1,ii,npane(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         ie(ind(ik,ii,n))=ind(ik+1-ii,1,npane(n)-100)
        enddo    ! ii loop
      endif      ! (npane(n).lt.100)
      if(npanw(n).lt.100)then
        do ii=1,ik
         iw(ind(1,ii,n))=ind(ik,ii,npanw(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         iw(ind(1,ii,n))=ind(ik+1-ii,ik,npanw(n)-100)
        enddo    ! ii loop
      endif      ! (npanw(n).lt.100)
      if(npans(n).lt.100)then
        do ii=1,ik
         is(ind(ii,1,n))=ind(ii,ik,npans(n))
        enddo    ! ii loop
      else
        do ii=1,ik
         is(ind(ii,1,n))=ind(ik,ik+1-ii,npans(n)-100)
        enddo    ! ii loop
      endif      ! (npans(n).lt.100)
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
      if(nmaxpr==1.and.mydiag)then
        write(6,*)'sig(lev) ',siglev
        write(6,*)'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd),psl(idjd),pmsl(idjd)
      endif
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