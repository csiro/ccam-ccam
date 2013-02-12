      ! CCAM netcdf output routines

      ! itype=1     write outfile history file (compressed)
      ! itype=-1    write restart file (uncompressed)
      ! localhist=f single processor output 
      ! localhist=t parallel output for each processor

      !--------------------------------------------------------------
      ! CONFIGURE DIMENSIONS FOR OUTPUT NETCDF FILES
      subroutine outcdf(rundate,nmi,itype,iaero,nstagin)

      use cc_mpi                            ! CC MPI routines
      use infile                            ! Input file routines
      use liqwpar_m                         ! Cloud water mixing ratios
      use mlo, only : wlev                  ! Ocean physics and prognostic arrays
      use parmhdff_m                        ! Horizontal diffusion parameters
      use tracers_m                         ! Tracer data

      implicit none

      include 'newmpar.h'                   ! Grid parameters
      include 'dates.h'                     ! Date data
      include 'filnames.h'                  ! Filenames
      include 'kuocom.h'                    ! Convection parameters
      include 'parm.h'                      ! Model configuration
      include 'parmdyn.h'                   ! Dynamics parameters
      include 'parmgeom.h'                  ! Coordinate data
      include 'parmhor.h'                   ! Horizontal advection parameters
      include 'parmvert.h'                  ! Vertical advection parameters

      integer ixp,iyp,idlev,idnt,idms,idoc
      common/cdfind/ixp,iyp,idlev,idnt,
     &  idms,idoc                           ! Output file dimension data
      integer leap
      common/leap_yr/leap                   ! Leap year (1 to allow leap years)

      integer, parameter :: nihead=54
      integer, parameter :: nrhead=14
      integer, dimension(nihead) :: nahead
      integer, dimension(4), save :: dim,dims,dimo
      integer nmi, itype, iaero, nstagin
      integer xdim,ydim,zdim,tdim,msdim,ocdim
      integer icy, icm, icd, ich, icmi, ics, idv, imode
      integer, save :: idnc=0, iarch=0, idnc0=0
      real, dimension(nrhead) :: ahead
      character(len=180) cdffile
      character(len=33) grdtim
      character(len=20) timorg
      character(len=8) rundate
      character(len=3), dimension(12) :: month
      logical local

      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/


      ! The localhist variable controls whether the local file option
      !  is used.
      !local = localhist .and. itype == 1 ! Only for outfile
      local = localhist                   ! outfile and restfile

      ! Determine file names depending on output
      if(myid==0 .or. local)then
       ! File setup follows
       if(itype==1)then
        ! itype=1 outfile
        iarch=iarch+1
        if(local)then
           write(cdffile,"(a,'.',i6.6)") trim(ofile), myid
        else
           cdffile=ofile
        endif
       else
        ! itype=-1 restfile
        iarch=1
        if(local)then
           write(cdffile,"(a,'.',i6.6)") trim(restfile), myid
        else
           cdffile=restfile
        endif
        idnc=0
       endif ! ( itype==1)then

       if (myid==0.and.iarch>1) then
         write(6,'("outcdf itype,idnc,iarch,cdffile=",3i5," ",a80)')
     &                   itype,idnc,iarch,cdffile
       end if

       ! Open new file
       if(iarch==1)then
        if (myid==0) then
          write(6,'("nccre of itype,cdffile=",i5," ",a80)') 
     &                   itype,cdffile
        end if
        call ccnf_create(cdffile,idnc)
        ! Turn off the data filling
        call ccnf_nofill(idnc)
        ! Create dimensions, lon, runtopo.shlat
        if(local)then
           call ccnf_def_dim(idnc,'longitude',il,xdim)
           call ccnf_def_dim(idnc,'latitude',jl,ydim)
        else
           call ccnf_def_dim(idnc,'longitude',il_g,xdim)
           call ccnf_def_dim(idnc,'latitude',jl_g,ydim)
        endif
        call ccnf_def_dim(idnc,'lev',kl,zdim)
        call ccnf_def_dim(idnc,'zsoil',ms,msdim)
        ocdim=0
        if (abs(nmlo)>0..and.abs(nmlo)<=9) then
          call ccnf_def_dim(idnc,'olev',ol,ocdim)
        end if
        call ccnf_def_dimu(idnc,'time',tdim)
        if (myid==0) then
          write(6,*) "xdim,ydim,zdim,tdim"
          write(6,*)  xdim,ydim,zdim,tdim
        end if

        ! atmosphere dimensions
        dim = (/ xdim, ydim, zdim, tdim /)

        ! soil dimensions
        dims = (/ xdim, ydim, msdim, tdim /)

        ! ocean dimensions
        dimo = (/ xdim, ydim, ocdim, tdim /)

        ! Define coords.
        call ccnf_def_var(idnc,'longitude','float',1,dim(1:1),ixp)
        call ccnf_put_att(idnc,ixp,'point_spacing',4,'even')
        call ccnf_put_att(idnc,ixp,'units',12,'degrees_east')
        call ccnf_def_var(idnc,'latitude','float',1,dim(2:2),iyp)
        call ccnf_put_att(idnc,iyp,'point_spacing',4,'even')
        call ccnf_put_att(idnc,iyp,'units',13,'degrees_north')
        if (myid==0) then
          write(6,*) 'ixp,iyp=',ixp,iyp
        end if

        call ccnf_def_var(idnc,'lev','float',1,dim(3:3),idlev)
        call ccnf_put_att(idnc,idlev,'positive',4,'down')
        call ccnf_put_att(idnc,idlev,'point_spacing',6,'uneven')
        call ccnf_put_att(idnc,idlev,'units',11,'sigma_level')
        call ccnf_put_att(idnc,idlev,'long_name',11,'sigma_level')
        if (myid==0) then
          write(6,*) 'idlev=',idlev
        end if

        call ccnf_def_var(idnc,'zsoil','float',1,dims(3:3),idms)
        call ccnf_put_att(idnc,idms,'point_spacing',6,'uneven')
        call ccnf_put_att(idnc,idms,'units',1,'m')
        if (myid==0) then
          write(6,*) 'idms=',idms
        end if
        
        if (abs(nmlo)>0.and.abs(nmlo)<=9) then
          call ccnf_def_var(idnc,'olev','float',1,dimo(3:3),idoc)
          call ccnf_put_att(idnc,idoc,'point_spacing',6,'uneven')
          call ccnf_put_att(idnc,idoc,'units',11,'sigma_level')
          if (myid==0) then
            write(6,*) 'idoc=',idoc
          end if
        end if

        call ccnf_def_var(idnc,'time','float',1,dim(4:4),idnt)
        call ccnf_put_att(idnc,idnt,'point_spacing',4,'even')
        if (myid==0) then
          write(6,*) 'tdim,idnc=',tdim,idnc
          write(6,*) 'idnt=',idnt
          write(6,*) 'kdate,ktime,ktau=',kdate,ktime,ktau
        end if

        icy=kdate/10000
        icm=max(1,min(12,(kdate-icy*10000)/100))
        icd=max(1,min(31,(kdate-icy*10000-icm*100)))
        if(icy<100)icy=icy+1900
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &               icd,month(icm),icy,ich,icmi,ics
        call ccnf_put_att(idnc,idnt,'time_origin',20,timorg)
        write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",
     &       2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
        call ccnf_put_att(idnc,idnt,'units',33,grdtim)
        if (leap==0) then
          call ccnf_put_att(idnc,idnt,'calendar',6,'noleap')
        end if
        if (myid==0) then
          write(6,*) 'timorg=',timorg
          write(6,*) 'grdtim=',grdtim
        end if

c       create the attributes of the header record of the file
        nahead(1)=il_g       ! needed by cc2hist
        nahead(2)=jl_g       ! needed by cc2hist
        nahead(3)=kl         ! needed by cc2hist
        nahead(4)=m
        nahead(5)=0          ! nsd not used now
        nahead(6)=io_in
        nahead(7)=nbd
        nahead(8)=0          ! not needed now  
        nahead(9)=mex
        nahead(10)=mup
        nahead(11)=nem
        nahead(12)=mtimer
        nahead(13)=nmi
        nahead(14)=nint(dt)  ! needed by cc2hist
        nahead(15)=0         ! not needed now 
        nahead(16)=nhor
        nahead(17)=nkuo
        nahead(18)=khdif
        nahead(19)=kl        ! needed by cc2hist (was kwt)
        nahead(20)=0  !iaa
        nahead(21)=0  !jaa
        nahead(22)=nvad
        nahead(23)=0       ! not needed now      
        nahead(24)=0  !lbd
        nahead(25)=nrun
        nahead(26)=nrunx
        nahead(27)=khor
        nahead(28)=ksc
        nahead(29)=kountr
        nahead(30)=ndiur
        nahead(31)=0  ! spare
        nahead(32)=nhorps
        nahead(33)=nsoil
        nahead(34)=ms        ! needed by cc2hist
        nahead(35)=ntsur
        nahead(36)=nrad
        nahead(37)=kuocb
        nahead(38)=nvmix
        nahead(39)=ntsea
        nahead(40)=0    
        nahead(41)=nextout
        nahead(42)=ilt
        nahead(43)=ntrac     ! needed by cc2hist
        nahead(44)=nsib
        nahead(45)=nrungcm
        nahead(46)=ncvmix
        nahead(47)=ngwd
        nahead(48)=lgwd
        nahead(49)=mup
        nahead(50)=nritch_t
        nahead(51)=ldr
        nahead(52)=nevapls
        nahead(53)=nevapcc
        nahead(54)=nt_adv
        ahead(1)=ds
        ahead(2)=0.  !difknbd
        ahead(3)=0.  ! was rhkuo for kuo scheme
        ahead(4)=0.  !du
        ahead(5)=rlong0     ! needed by cc2hist
        ahead(6)=rlat0      ! needed by cc2hist
        ahead(7)=schmidt    ! needed by cc2hist
        ahead(8)=0.  !stl2
        ahead(9)=0.  !relaxt
        ahead(10)=0.  !hourbd
        ahead(11)=tss_sh
        ahead(12)=vmodmin
        ahead(13)=av_vmod
        ahead(14)=epsp
        if (myid==0) then
          write(6,'("nahead=",(20i4))') nahead
          write(6,*) "ahead=",ahead
        end if
        call ccnf_put_attg(idnc,'int_header',nihead,nahead)
        call ccnf_put_attg(idnc,'real_header',nrhead,ahead)
        call ccnf_put_attg(idnc,'date_header',rundate)

        call ccnf_def_var0(idnc,'ds','float',idv)
        call ccnf_def_var0(idnc,'dt','float',idv)
       endif ! ( iarch=1)then

      endif ! (myid==0.or.local)
      
      ! openhist writes some fields so needs to be called by all processes
      call openhist(iarch,itype,dim,local,idnc,iaero,nstagin)

      if(myid==0.or.local)then
        call ccnf_sync(idnc)
        if(ktau==ntau)then
          call ccnf_close(idnc)
          if (myid==0) then
            write(6,*) "calling ncclos(idnc) ",idnc
          end if
        endif
      endif    ! (myid==0.or.local)

      return   ! outcdf  
      end
      
      !--------------------------------------------------------------
      ! CREATE ATTRIBUTES AND WRITE OUTPUT
      subroutine openhist(iarch,itype,dim,local,idnc,iaero,nstagin)

      use aerosolldr                            ! LDR prognostic aerosols
      use arrays_m                              ! Atmosphere dyamics prognostic arrays
      use ateb, only : atebsave                 ! Urban
      use cable_ccam, only : savetile,          ! CABLE interface
     &      savetiledef
      use casadimension, only : mplant,mlitter, ! CASA dimensions
     &      msoil
      use carbpools_m                           ! Carbon pools
      use cc_mpi                                ! CC MPI routines
      use cfrac_m                               ! Cloud fraction
      use cable_def_types_mod, only : ncs, ncp  ! CABLE dimensions
      use dpsdt_m                               ! Vertical velocity
      use extraout_m                            ! Additional diagnostics
      use gdrag_m                               ! Gravity wave drag
      use histave_m                             ! Time average arrays
      use infile                                ! Input file routines
      use latlong_m                             ! Lat/lon coordinates
      use liqwpar_m                             ! Cloud water mixing ratios
      use map_m                                 ! Grid map arrays
      use mlo, only : wlev,mlosave,mlodiag,     ! Ocean physics and prognostic arrays
     &      mloexpdep
      use mlodynamics                           ! Ocean dynamics
      use morepbl_m                             ! Additional boundary layer diagnostics
      use nharrs_m                              ! Non-hydrostatic atmosphere arrays
      use nsibd_m                               ! Land-surface arrays
      use pbl_m                                 ! Boundary layer arrays
      use prec_m                                ! Precipitation
      use raddiag_m                             ! Radiation diagnostic
      use savuvt_m                              ! Saved dynamic arrays
      use savuv1_m                              ! Saved dynamic arrays
      use screen_m                              ! Screen level diagnostics
      use sigs_m                                ! Atmosphere sigma levels
      use soil_m                                ! Soil and surface data
      use soilsnow_m                            ! Soil, snow and surface data
      use tkeeps, only : tke,eps                ! TKE-EPS boundary layer
      use tracermodule, only : tracmax,tracmin, ! Tracer routines
     &      tracname,writetrpm
      use tracers_m                             ! Tracer data
      use vegpar_m                              ! Vegetation arrays
      use vvel_m                                ! Additional vertical velocity
      use work2_m                               ! Diagnostic arrays
      use xarrs_m, only : pslx                  ! Saved dynamic arrays

      implicit none

      include 'newmpar.h'                       ! Grid parameters
      include 'const_phys.h'                    ! Physical constants
      include 'dates.h'                         ! Date data
      include 'filnames.h'                      ! Filenames
      include 'kuocom.h'                        ! Convection parameters
      include 'parm.h'                          ! Model configuration
      include 'parmdyn.h'                       ! Dynamics parameters
      include 'parmvert.h'                      ! Vertical advection parameters
      include 'soilv.h'                         ! Soil parameters
      include 'trcom2.h'                        ! Station data
      include 'version.h'                       ! Model version data

      integer ixp,iyp,idlev,idnt,idms,idoc
      common/cdfind/ixp,iyp,idlev,idnt,idms,
     &  idoc                                    ! Output file dimension data

      integer i, idkdate, idktau, idktime, idmtimer, idnteg, idnter
      integer idv, iq, isoil, j, k, n, igas, idnc
      integer iarch, itype, iaero, nstagin, idum
      integer, dimension(4), intent(in) :: dim
      integer, dimension(3) :: idim
      integer, dimension(2) :: iduma
      real trmax, trmin
      real, dimension(ms) :: zsoil
      real, dimension(il_g) :: xpnt
      real, dimension(jl_g) :: ypnt
      real, dimension(ifull) :: aa,bb,cc
      real, dimension(ifull) :: ocndep,ocnheight
      real, dimension(ifull,kl) :: tmpry
      real, dimension(ifull,wlev,4) :: mlodwn
      real, dimension(ifull,11) :: micdwn
      real, dimension(ifull,24) :: atebdwn
      character(len=50) expdesc
      character(len=40) lname
      character(len=21) mnam,nnam
      character(len=8) vname
      character(len=3) trnum
      character(len=3), dimension(12) :: mon
      logical, intent(in) :: local
      logical lwrite,lave,lrad,lday,tst

      data mon/'JAN','FEB','MAR','APR','MAY','JUN'
     &        ,'JUL','AUG','SEP','OCT','NOV','DEC'/
     
      lwrite=ktau>0
      lave=mod(ktau,nperavg)==0.or.ktau==ntau
      lave=lave.and.ktau>0
      lrad=mod(ktau,kountr)==0.or.ktau==ntau
      lrad=lrad.and.ktau>0
      lday=mod(ktau,nperday)==0.or.ktau==ntau
      lday=lday.and.ktau>0

      idim(1)=dim(1)
      idim(2)=dim(2)
      idim(3)=dim(4)

      if(myid == 0 .or. local)then
       if (myid==0) then
        write(6,*) 'openhist itype,iarch,idnc=',itype,iarch,idnc
       end if

c      if this is the first archive, set up some global attributes
       if(iarch==1) then

c       Create global attributes
c       Model run number
        if (myid==0) then
         write(6,*) 'dim=',dim
         write(6,*) 'idim=',idim
         write(6,*) 'nrun=',nrun
        end if
        call ccnf_put_attg(idnc,'nrun',nrun)

c       Experiment description
        expdesc = 'CCAM model run'
        call ccnf_put_attg(idnc,'expdesc',expdesc)

c       Model version
        call ccnf_put_attg(idnc,'version',version)

        if(local)then
           call ccnf_put_attg(idnc,'processor_num',myid)
           call ccnf_put_attg(idnc,'nproc',nproc)
#ifdef uniform_decomp
           call ccnf_put_attg(idnc,'decomp','uniform')
#else
           call ccnf_put_attg(idnc,'decomp','face')
#endif
        endif           

c       Sigma levels
        if (myid==0) then
         write(6,*) 'sig=',sig
        end if
        call ccnf_put_attg(idnc,'sigma',kl,sig)

        lname = 'year-month-day at start of run'
        call ccnf_def_var(idnc,'kdate','int',1,dim(4:4),idkdate)
        call ccnf_put_att(idnc,idkdate,'long_name',
     &                         len_trim(lname),lname)

        lname = 'hour-minute at start of run'
        call ccnf_def_var(idnc,'ktime','int',1,dim(4:4),idktime)
        call ccnf_put_att(idnc,idktime,'long_name',
     &                         len_trim(lname),lname)

        lname = 'timer (hrs)'
        call ccnf_def_var(idnc,'timer','float',1,dim(4:4),idnter)
        call ccnf_put_att(idnc,idnter,'long_name',
     &                         len_trim(lname),lname)

        lname = 'mtimer (mins)'
        call ccnf_def_var(idnc,'mtimer','int',1,dim(4:4),idmtimer)
        call ccnf_put_att(idnc,idmtimer,'long_name',
     &                         len_trim(lname),lname)

        lname = 'timeg (UTC)'
        call ccnf_def_var(idnc,'timeg','float',1,dim(4:4),idnteg)
        call ccnf_put_att(idnc,idnteg,'long_name',
     &                         len_trim(lname),lname)

        lname = 'number of time steps from start'
        call ccnf_def_var(idnc,'ktau','int',1,dim(4:4),idktau)
        call ccnf_put_att(idnc,idktau,'long_name',
     &                         len_trim(lname),lname)

        call ccnf_def_var(idnc,'sigma','float',1,dim(3:3),idv)
        call ccnf_put_att(idnc,idv,'positive',
     &                         len_trim('down'),'down')

        lname = 'atm stag direction'
        call ccnf_def_var(idnc,'nstag','int',1,dim(4:4),idv)
        call ccnf_put_att(idnc,idv,'long_name',
     &                         len_trim(lname),lname)

        lname = 'atm unstag direction'
        call ccnf_def_var(idnc,'nstagu','int',1,dim(4:4),idv)
        call ccnf_put_att(idnc,idv,'long_name',
     &                         len_trim(lname),lname)

        lname = 'atm stag offset'
        call ccnf_def_var(idnc,'nstagoff','int',1,dim(4:4),idv)
        call ccnf_put_att(idnc,idv,'long_name',
     &                         len_trim(lname),lname)

        if ((nmlo<0.and.nmlo>=-9).or.
     &      (nmlo>0.and.nmlo<=9.and.itype==-1)) then
          lname = 'ocn stag offset'
          call ccnf_def_var(idnc,'nstagoffmlo','int',1,dim(4:4),idv)
          call ccnf_put_att(idnc,idv,'long_name',
     &                           len_trim(lname),lname)     
        end if

        if (myid==0) then
         write(6,*) 'define attributes of variables'
        end if

c       For time invariant surface fields
        lname = 'Surface geopotential'
        call attrib(idnc,idim(1:2),2,'zht',lname,'m2/s2',-1000.,90.e3,
     &              0,-1)
        lname = 'Std Dev of surface height'
        call attrib(idnc,idim(1:2),2,'he',lname,'m',-1000.,90.e3,
     &              0,-1)
        lname = 'Map factor'
        call attrib(idnc,idim(1:2),2,'map',lname,'none',.001,1500.,
     &              0,itype)
        lname = 'Coriolis factor'
        call attrib(idnc,idim(1:2),2,'cor',lname,'1/sec',-1.5e-4,
     &              1.5e-4,0,itype)
        lname = 'Urban fraction'
        call attrib(idnc,idim(1:2),2,'sigmu',lname,'none',0.,3.25,
     &              0,itype)
        lname = 'Soil type'
        call attrib(idnc,idim(1:2),2,'soilt',lname,'none',0.,65.,
     &              0,itype)
        lname = 'Vegetation type'
        call attrib(idnc,idim(1:2),2,'vegt',lname,'none',0.,65.,
     &              0,itype)

        if ((nmlo<0.and.nmlo>=-9).or.
     &      (nmlo>0.and.nmlo<=9.and.itype==-1)) then
          lname = 'water depth'
          call attrib(idnc,idim(1:2),2,'ocndepth',lname,'m',0.,32500.,
     &                0,itype)
        end if

c       For time varying surface fields
        if (nsib==6.or.nsib==7) then
          lname = 'Stomatal resistance'
          call attrib(idnc,idim,3,'rs',lname,'none',0.,1000.,0,itype)
        else
          lname = 'Minimum stomatal resistance'
          call attrib(idnc,idim(1:2),2,'rsmin',lname,'none',0.,1000.,0,
     &                itype)
        end if
        lname = 'Vegetation fraction'
        call attrib(idnc,idim,3,'sigmf',lname,'none',0.,3.25,0,itype)
        lname ='Scaled Log Surface pressure'
        call attrib(idnc,idim,3,'psf',lname,'none',-1.3,0.2,0,itype)
        lname ='Mean sea level pressure'
        call attrib(idnc,idim,3,'pmsl',lname,'hPa',800.,1200.,0,itype)
        lname = 'Surface roughness'
        call attrib(idnc,idim,3,'zolnd',lname,'m',0.,65.,0,-1) ! -1=long
        lname = 'Leaf area index'
        call attrib(idnc,idim,3,'lai',lname,'none',0.,32.5,0,itype)
        lname = 'Surface temperature'
        call attrib(idnc,idim,3,'tsu',lname,'K',100.,425.,0,itype)
        lname = 'Pan temperature'
        call attrib(idnc,idim,3,'tpan',lname,'K',100.,425.,0,itype)
        lname = 'Precipitation'
        call attrib(idnc,idim,3,'rnd',lname,'mm/day',0.,1300.,0,itype)
        lname = 'Convective precipitation'
        call attrib(idnc,idim,3,'rnc',lname,'mm/day',0.,1300.,0,itype)
        call attrib(idnc,idim,3,'sno','snowfall','mm/day',0.,1300.,0,
     &              itype)
        call attrib(idnc,idim,3,'runoff','Runoff','mm/day',0.,1300.,0,
     &              itype)
        lname = 'Surface albedo'
        call attrib(idnc,idim,3,'alb',lname,'none',0.,1.,0,itype)
        lname = 'Fraction of canopy that is wet'
        call attrib(idnc,idim,3,'fwet',lname,'none',0.,1.,0,itype)

        lname = 'Snow depth (liquid water)'
        call attrib(idnc,idim,3,'snd',lname,'mm',0.,6500.,0,-1)  ! -1=long
        lname = 'Soil temperature lev 1'
        call attrib(idnc,idim,3,'tgg1',lname,'K',100.,425.,0,itype)
        lname = 'Soil temperature lev 2'
        call attrib(idnc,idim,3,'tgg2',lname,'K',100.,425.,0,itype)
        lname = 'Soil temperature lev 3'
        call attrib(idnc,idim,3,'tgg3',lname,'K',100.,425.,0,itype)
        lname = 'Soil temperature lev 4'
        call attrib(idnc,idim,3,'tgg4',lname,'K',100.,425.,0,itype)
        lname = 'Soil temperature lev 5'
        call attrib(idnc,idim,3,'tgg5',lname,'K',100.,425.,0,itype)
        lname = 'Soil temperature lev 6'
        call attrib(idnc,idim,3,'tgg6',lname,'K',100.,425.,0,itype)
 
        if ((nmlo<0.and.nmlo>=-9).or.
     &      (nmlo>0.and.nmlo<=9.and.itype==-1)) then
          do k=ms+1,wlev
           write(lname,'("soil/ocean temperature lev ",I2)') k
           write(vname,'("tgg",I2.2)') k
           call attrib(idnc,idim,3,vname,lname,'K',100.,425.,0,itype)
          end do
          do k=1,wlev
           write(lname,'("ocean salinity lev ",I2)') k
           write(vname,'("sal",I2.2)') k
           call attrib(idnc,idim,3,vname,lname,'PSU',0.,130.,0,itype)
          end do
          do k=1,wlev
           write(lname,'("x-component current lev ",I2)') k
           write(vname,'("uoc",I2.2)') k
           call attrib(idnc,idim,3,vname,lname,'m/s',-65.,65.,0,itype)
           write(lname,'("y-component current lev ",I2)') k
           write(vname,'("voc",I2.2)') k
           call attrib(idnc,idim,3,vname,lname,'m/s',-65.,65.,0,itype)
          end do
          lname = 'water surface height'
          call attrib(idnc,idim,3,'ocheight',lname,'m',-130.,130.,0,
     &                itype)          
          lname = 'Snow temperature lev 1'
          call attrib(idnc,idim,3,'tggsn1',lname,'K',100.,425.,0,
     &                  itype)
          lname = 'Snow temperature lev 2'
          call attrib(idnc,idim,3,'tggsn2',lname,'K',100.,425.,0,
     &                  itype)
          lname = 'Snow temperature lev 3'
          call attrib(idnc,idim,3,'tggsn3',lname,'K',100.,425.,0,
     &                  itype)
          lname = 'Ice temperature lev 4'
          call attrib(idnc,idim,3,'tggsn4',lname,'K',100.,425.,0,itype)
          lname = 'Ice heat store'
          call attrib(idnc,idim,3,'sto',lname,'J/m2',0.,1300.,0,itype)
          lname = 'x-component ice velocity'
          call attrib(idnc,idim,3,'uic',lname,'m/s',-65.,65.,0,itype)
          lname = 'y-component ice velocity'
          call attrib(idnc,idim,3,'vic',lname,'m/s',-65.,65.,0,itype)
          lname = 'Ice salinity'
          call attrib(idnc,idim,3,'icesal',lname,'PSU',0.,130.,0,
     &                itype)
          if (abs(nmlo)>=2) then
            lname = 'Surface water depth'
            call attrib(idnc,idim,3,'swater',lname,'mm',0.,6.5E3,0,
     &                  -1) ! -1 = long
            lname = 'Surface water salinity'
            call attrib(idnc,idim,3,'ssalin',lname,'PSU',0.,130.,0,
     &                  itype)
          end if
        end if

       ! lname = 'Soil moisture lev 1'
       ! call attrib(idnc,idim,3,'wb1',lname,'m3/m3',0.,1.,0)
       ! lname = 'Soil moisture lev 2'
       ! call attrib(idnc,idim,3,'wb2',lname,'m3/m3',0.,1.,0)
       ! lname = 'Soil moisture lev 3'
       ! call attrib(idnc,idim,3,'wb3',lname,'m3/m3',0.,1.,0)
       ! lname = 'Soil moisture lev 4'
       ! call attrib(idnc,idim,3,'wb4',lname,'m3/m3',0.,1.,0)
       ! lname = 'Soil moisture lev 5'
       ! call attrib(idnc,idim,3,'wb5',lname,'m3/m3',0.,1.,0)
       ! lname = 'Soil moisture lev 6'
       ! call attrib(idnc,idim,3,'wb6',lname,'m3/m3',0.,1.,0)
        lname = 'Wetness fraction layer 1' ! 5. for frozen sand
        call attrib(idnc,idim,3,'wetfrac1',lname,'none',-6.5,6.5,0,
     &              itype)
        lname = 'Wetness fraction layer 2'
        call attrib(idnc,idim,3,'wetfrac2',lname,'none',-6.5,6.5,0,
     &              itype)
        lname = 'Wetness fraction layer 3'
        call attrib(idnc,idim,3,'wetfrac3',lname,'none',-6.5,6.5,0,
     &              itype)
        lname = 'Wetness fraction layer 4'
        call attrib(idnc,idim,3,'wetfrac4',lname,'none',-6.5,6.5,0,
     &              itype)
        lname = 'Wetness fraction layer 5'
        call attrib(idnc,idim,3,'wetfrac5',lname,'none',-6.5,6.5,0,
     &              itype)
        lname = 'Wetness fraction layer 6'
        call attrib(idnc,idim,3,'wetfrac6',lname,'none',-6.5,6.5,0,
     &              itype)
     
        ! PH - Add wetfac to output for mbase=-19 option
        lname = 'Surface wetness fraction'
        call attrib(idnc,idim,3,'wetfac',lname,'none',-6.5,6.5,0,
     &              itype)

        lname = 'Sea ice depth'
        call attrib(idnc,idim,3,'siced',lname,'m',0.,65.,0,-1)
        lname = 'Sea ice fraction'
        call attrib(idnc,idim,3,'fracice',lname,'none',0.,6.5,0,itype)
        lname = '10m wind speed'
        call attrib(idnc,idim,3,'u10',lname,'m/s',0.,130.,0,itype)

        lname = 'Maximum precip rate in a timestep'
        call attrib(idnc,idim,3,'maxrnd',lname,'mm/day',0.,2600.,1,
     &              itype)
        lname = 'Maximum screen temperature'
        call attrib(idnc,idim,3,'tmaxscr',lname,'K',100.,425.,1,itype)
        lname = 'Minimum screen temperature'
        call attrib(idnc,idim,3,'tminscr',lname,'K',100.,425.,1,itype)
        lname = 'Maximum screen relative humidity'
        call attrib(idnc,idim,3,'rhmaxscr',lname,'%',0.,200.,1,itype)
        lname = 'Minimum screen relative humidity'
        call attrib(idnc,idim,3,'rhminscr',lname,'%',0.,200.,1,itype)
        lname = 'Maximum daily Cape'
        call attrib(idnc,idim,3,'capemax',lname,'J/kg',0.,20000.,1,
     &              itype) 
        lname = 'x-component max 10m wind'
        call attrib(idnc,idim,3,'u10max',lname,'m/s',-99.,99.,1,itype)
        lname = 'y-component max 10m wind'
        call attrib(idnc,idim,3,'v10max',lname,'m/s',-99.,99.,1,itype)
        lname = 'x-component max level_1 wind'
        call attrib(idnc,idim,3,'u1max',lname,'m/s',-99.,99.,1,itype)
        lname = 'y-component max level_1 wind'
        call attrib(idnc,idim,3,'v1max',lname,'m/s',-99.,99.,1,itype)
        lname = 'x-component max level_2 wind'
        call attrib(idnc,idim,3,'u2max',lname,'m/s',-99.,99.,1,itype)
        lname = 'y-component max level_2 wind'
        call attrib(idnc,idim,3,'v2max',lname,'m/s',-99.,99.,1,itype)
        lname = '3hr precipitation'
        call attrib(idnc,idim,3,'rnd03',lname,'mm',0.,1300.,1,itype)
        lname = '6hr precipitation'
        call attrib(idnc,idim,3,'rnd06',lname,'mm',0.,1300.,1,itype)
        lname = '9hr precipitation'
        call attrib(idnc,idim,3,'rnd09',lname,'mm',0.,1300.,1,itype)
        lname = '12hr precipitation'
        call attrib(idnc,idim,3,'rnd12',lname,'mm',0.,1300.,1,itype)
        lname = '15hr precipitation'
        call attrib(idnc,idim,3,'rnd15',lname,'mm',0.,1300.,1,itype)
        lname = '18hr precipitation'
        call attrib(idnc,idim,3,'rnd18',lname,'mm',0.,1300.,1,itype)
        lname = '21hr precipitation'
        call attrib(idnc,idim,3,'rnd21',lname,'mm',0.,1300.,1,itype)
        lname = '24hr precipitation'
        call attrib(idnc,idim,3,'rnd24',lname,'mm',0.,1300.,1,itype)
        if(nextout>=2) then  ! 6-hourly u10, v10, tscr, rh1
         mnam ='x-component 10m wind '
         nnam ='y-component 10m wind '
         call attrib(idnc,idim,3,'u10_06',mnam//'6hr','m/s',-99.,99.,1,
     &               itype)
         call attrib(idnc,idim,3,'v10_06',nnam//'6hr','m/s',-99.,99.,1,
     &               itype)
         call attrib(idnc,idim,3,'u10_12',mnam//'12hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'v10_12',nnam//'12hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'u10_18',mnam//'18hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'v10_18',nnam//'18hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'u10_24',mnam//'24hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'v10_24',nnam//'24hr','m/s',-99.,99.,
     &               1,itype)
         mnam ='tscrn 3-hrly'
         nnam ='rhum level_1 3-hrly'
         call attrib(idnc,idim,3,'tscr_06',mnam//'6hr', 'K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'tscr_12',mnam//'12hr','K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'tscr_18',mnam//'18hr','K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'tscr_24',mnam//'24hr','K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'rh1_06', nnam//'6hr', '%',-9.,200.,
     &               1,itype)
         call attrib(idnc,idim,3,'rh1_12', nnam//'12hr','%',-9.,200.,
     &               1,itype)
         call attrib(idnc,idim,3,'rh1_18', nnam//'18hr','%',-9.,200.,
     &               1,itype)
         call attrib(idnc,idim,3,'rh1_24', nnam//'24hr','%',-9.,200.,
     &               1,itype)
        endif     ! (nextout>=2)
        if(nextout>=3) then  ! also 3-hourly u10, v10, tscr, rh1
         call attrib(idnc,idim,3,'tscr_03',mnam//'3hr', 'K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'tscr_09',mnam//'9hr', 'K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'tscr_15',mnam//'15hr','K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'tscr_21',mnam//'21hr','K',100.,425.,
     &               1,itype)
         call attrib(idnc,idim,3,'rh1_03', nnam//'3hr', '%',-9.,200.,1,
     &               itype)
         call attrib(idnc,idim,3,'rh1_09', nnam//'9hr', '%',-9.,200.,1,
     &               itype)
         call attrib(idnc,idim,3,'rh1_15', nnam//'15hr','%',-9.,200.,1,
     &               itype)
         call attrib(idnc,idim,3,'rh1_21', nnam//'21hr','%',-9.,200.,1,
     &               itype)
         mnam ='x-component 10m wind '
         nnam ='y-component 10m wind '
         call attrib(idnc,idim,3,'u10_03',mnam//'3hr','m/s',-99.,99.,1,
     &               itype)
         call attrib(idnc,idim,3,'v10_03',nnam//'3hr','m/s',-99.,99.,1,
     &               itype)
         call attrib(idnc,idim,3,'u10_09',mnam//'9hr','m/s',-99.,99.,1,
     &               itype)
         call attrib(idnc,idim,3,'v10_09',nnam//'9hr','m/s',-99.,99.,1,
     &               itype)
         call attrib(idnc,idim,3,'u10_15',mnam//'15hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'v10_15',nnam//'15hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'u10_21',mnam//'21hr','m/s',-99.,99.,
     &               1,itype)
         call attrib(idnc,idim,3,'v10_21',nnam//'21hr','m/s',-99.,99.,
     &               1,itype)
        endif     ! (nextout>=3)

        lname = 'Average screen temperature'
        call attrib(idnc,idim,3,'tscr_ave',lname,'K',100.,425.,0,itype)
        lname = 'Avg cloud base'
        call attrib(idnc,idim,3,'cbas_ave',lname,'sigma',0.,1.1,0,itype)
        lname = 'Avg cloud top'
        call attrib(idnc,idim,3,'ctop_ave',lname,'sigma',0.,1.1,0,itype)
        lname = 'Avg dew flux'
        call attrib(idnc,idim,3,'dew_ave',lname,'W/m2',-100.,1000.,0,
     &              itype)
        lname = 'Avg evaporation'
        call attrib(idnc,idim,3,'evap',lname,'mm',-100.,100.,0,itype)
        lname = 'Avg potential "pan" evaporation'
        call attrib(idnc,idim,3,'epan_ave',lname,'W/m2',-1000.,10.e3,0,
     &              itype)
        lname = 'Avg potential evaporation'
        call attrib(idnc,idim,3,'epot_ave',lname,'W/m2',-1000.,10.e3,0,
     &              itype)
        lname = 'Avg latent heat flux'
        call attrib(idnc,idim,3,'eg_ave',lname,'W/m2',-1000.,3000.,0,
     &              itype)
        lname = 'Avg sensible heat flux'
        call attrib(idnc,idim,3,'fg_ave',lname,'W/m2',-3000.,3000.,0,
     &              itype)
        lname = 'Avg net radiation'
        call attrib(idnc,idim,3,'rnet_ave',lname,'none',-3000.,3000.,0,
     &              itype)
        lname = 'Avg flux into tgg1 layer'
        call attrib(idnc,idim,3,'ga_ave',lname,'W/m2',-1000.,1000.,0,
     &              itype)
        lname = 'Avg ice water path'
        call attrib(idnc,idim,3,'iwp_ave',lname,'kg/m2',0.,2.,0,itype)
        lname = 'Avg liquid water path'
        call attrib(idnc,idim,3,'lwp_ave',lname,'kg/m2',0.,2.,0,itype)
        lname = 'Low cloud ave'
        call attrib(idnc,idim,3,'cll',lname,'frac',0.,1.,0,itype)
        lname = 'Mid cloud ave'
        call attrib(idnc,idim,3,'clm',lname,'frac',0.,1.,0,itype)
        lname = 'Hi cloud ave'
        call attrib(idnc,idim,3,'clh',lname,'frac',0.,1.,0,itype)
        lname = 'Total cloud ave'
        call attrib(idnc,idim,3,'cld',lname,'frac',0.,1.,0,itype)
        lname = 'Avg soil moisture 1'
        call attrib(idnc,idim,3,'wb1_ave',lname,'m3/m3',0.,1.,0,itype)
        lname = 'Avg soil moisture 2'
        call attrib(idnc,idim,3,'wb2_ave',lname,'m3/m3',0.,1.,0,itype)
        lname = 'Avg soil moisture 3'
        call attrib(idnc,idim,3,'wb3_ave',lname,'m3/m3',0.,1.,0,itype)
        lname = 'Avg soil moisture 4'
        call attrib(idnc,idim,3,'wb4_ave',lname,'m3/m3',0.,1.,0,itype)
        lname = 'Avg soil moisture 5'
        call attrib(idnc,idim,3,'wb5_ave',lname,'m3/m3',0.,1.,0,itype)
        lname = 'Avg soil moisture 6'
        call attrib(idnc,idim,3,'wb6_ave',lname,'m3/m3',0.,1.,0,itype)
!          lname = 'Avg soil temperature 1'
!          call attrib(idnc,idim,3,'tgg1_ave',lname,'K',100.,425.,0)
!          lname = 'Avg soil temperature 2'
!          call attrib(idnc,idim,3,'tgg2_ave',lname,'K',100.,425.,0)
!          lname = 'Avg soil temperature 3'
!          call attrib(idnc,idim,3,'tgg3_ave',lname,'K',100.,425.,0)
!          lname = 'Avg soil temperature 4'
!          call attrib(idnc,idim,3,'tgg4_ave',lname,'K',100.,425.,0)
!          lname = 'Avg soil temperature 5'
!          call attrib(idnc,idim,3,'tgg5_ave',lname,'K',100.,425.,0)
!          lname = 'Avg soil temperature 6'
!          call attrib(idnc,idim,3,'tgg6_ave',lname,'K',100.,425.,0)
        lname = 'Avg surface temperature'
        call attrib(idnc,idim,3,'tsu_ave',lname,'K',100.,425.,0,itype)
        lname = 'Avg albedo'
        call attrib(idnc,idim,3,'alb_ave',lname,'none',0.,1.,0,itype)
        lname = 'Avg mean sea level pressure'
        call attrib(idnc,idim,3,'pmsl_ave',lname,'none',800.,1200.,0,
     &              itype)
        if (abs(nmlo)>0.and.abs(nmlo)<=9) then
          lname = 'Avg mixed layer depth'
          call attrib(idnc,idim,3,'mixd_ave',lname,'m',0.,1300.,0,
     &              itype)
        end if

        lname = 'Screen temperature'
        call attrib(idnc,idim,3,'tscrn',lname,'K',100.,425.,0,itype)
        lname = 'Screen mixing ratio'
        call attrib(idnc,idim,3,'qgscrn',lname,'kg/kg',0.,.06,0,itype)
        lname = 'Screen relative humidity'
        call attrib(idnc,idim,3,'rhscrn',lname,'%',0.,200.,0,itype)
        lname = 'Screen level wind speed'
        call attrib(idnc,idim,3,'uscrn',lname,'m/s',0.,65.,0,itype)
        lname = 'Net radiation'
        call attrib(idnc,idim,3,'rnet',lname,'W/m2',-3000.,3000.,0,
     &              itype)
        lname = 'Potential "pan" evaporation'
        call attrib(idnc,idim,3,'epan',lname,'W/m2',-1000.,10.e3,0,
     &              itype)
        lname = 'Latent heat flux'
        call attrib(idnc,idim,3,'eg',lname,'W/m2',-1000.,3000.,0,itype)
        lname = 'Sensible heat flux'
        call attrib(idnc,idim,3,'fg',lname,'W/m2',-3000.,3000.,0,itype)
        lname = 'x-component wind stress'
        call attrib(idnc,idim,3,'taux',lname,'N/m2',-50.,50.,0,itype)
        lname = 'y-component wind stress'
        call attrib(idnc,idim,3,'tauy',lname,'N/m2',-50.,50.,0,itype)
        if(nextout>=1) then
          if (myid==0) then
            write(6,*) 'nextout=',nextout
          end if
          lname = 'LW at TOA'
          call attrib(idnc,idim,3,'rtu_ave',lname,'W/m2',0.,800.,0,
     &                itype)
          lname = 'Clear sky LW at TOA'
          call attrib(idnc,idim,3,'rtc_ave',lname,'W/m2',0.,800.,0,
     &                itype)
          lname = 'LW downwelling at ground'
          call attrib(idnc,idim,3,'rgdn_ave',lname,'W/m2',-500.,1.e3,0,
     &                itype)
          lname = 'LW net at ground (+ve up)'
          call attrib(idnc,idim,3,'rgn_ave',lname,'W/m2',-500.,1000.,0,
     &                itype)
          lname = 'Clear sky LW at ground'
          call attrib(idnc,idim,3,'rgc_ave',lname,'W/m2',-500.,1000.,0,
     &                itype)
          lname = 'Solar in at TOA'
          call attrib(idnc,idim,3,'sint_ave',lname,'W/m2',0.,1600.,0,
     &                itype)
          lname = 'Solar out at TOA'
          call attrib(idnc,idim,3,'sot_ave',lname,'W/m2',0.,1000.,0,
     &                itype)
          lname = 'Clear sky SW out at TOA'
          call attrib(idnc,idim,3,'soc_ave',lname,'W/m2',0.,900.,0,
     &                itype)
          lname = 'Solar downwelling at ground'
          call attrib(idnc,idim,3,'sgdn_ave',lname,'W/m2',-500.,2.e3,0,
     &                itype)
          lname = 'Solar net at ground (+ve down)'
          call attrib(idnc,idim,3,'sgn_ave',lname,'W/m2',-500.,2000.,0,
     &                itype)
          lname = 'Sunshine hours'
          call attrib(idnc,idim,3,'sunhours',lname,'hrs',0.,64.5,0,
     &                itype)
          lname = 'Fraction of direct radiation'
          call attrib(idnc,idim,3,'fbeam_ave',lname,'none',-3.25,3.25,0,
     &                itype)
          lname = 'Surface pressure tendency'
          call attrib(idnc,idim,3,'dpsdt',lname,'hPa/day',-400.,400.,0,
     &                itype)
          lname = 'friction velocity'
          call attrib(idnc,idim,3,'ustar',lname,'m/s',0.,10.,0,itype)
        endif     ! (nextout>=1)
        if (nextout>=1.or.(nvmix==6.and.itype==-1)) then
          lname = 'PBL depth'
          call attrib(idnc,idim,3,'pblh',lname,'m',0.,13000.,0,itype)
        end if

        ! CABLE -----------------------------------------------------
        if (nsib==6.or.nsib==7) then
          if (nextout>=1.or.itype==-1) then
           if (ccycle==0) then
            lname = 'Carbon leaf pool'
            call attrib(idnc,idim,3,'cplant1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon wood pool'
            call attrib(idnc,idim,3,'cplant2',lname,'gC/m2',0.,65000.,
     &                  0,itype)
            lname = 'Carbon root pool'
            call attrib(idnc,idim,3,'cplant3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon soil fast pool'
            call attrib(idnc,idim,3,'csoil1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon soil slow pool'
            call attrib(idnc,idim,3,'csoil2',lname,'gC/m2',0.,6500.,
     &                  0,itype)
           else
            lname = 'Carbon leaf pool'
            call attrib(idnc,idim,3,'cplant1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen leaf pool'
            call attrib(idnc,idim,3,'nplant1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor leaf pool'
            call attrib(idnc,idim,3,'pplant1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon wood pool'
            call attrib(idnc,idim,3,'cplant2',lname,'gC/m2',0.,65000.,
     &                  0,itype)
            lname = 'Nitrogen wood pool'
            call attrib(idnc,idim,3,'nplant2',lname,'gC/m2',0.,65000.,
     &                  0,itype)
            lname = 'Phosphor wood pool'
            call attrib(idnc,idim,3,'pplant2',lname,'gC/m2',0.,65000.,
     &                  0,itype)
            lname = 'Carbon root pool'
            call attrib(idnc,idim,3,'cplant3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen root pool'
            call attrib(idnc,idim,3,'nplant3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor root pool'
            call attrib(idnc,idim,3,'pplant3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon met pool'
            call attrib(idnc,idim,3,'clitter1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen met pool'
            call attrib(idnc,idim,3,'nlitter1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor met pool'
            call attrib(idnc,idim,3,'plitter1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon str pool'
            call attrib(idnc,idim,3,'clitter2',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen str pool'
            call attrib(idnc,idim,3,'nlitter2',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor str pool'
            call attrib(idnc,idim,3,'plitter2',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon CWD pool'
            call attrib(idnc,idim,3,'clitter3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen CWD pool'
            call attrib(idnc,idim,3,'nlitter3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor CWD pool'
            call attrib(idnc,idim,3,'plitter3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon mic pool'
            call attrib(idnc,idim,3,'csoil1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen mic pool'
            call attrib(idnc,idim,3,'nsoil1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor mic pool'
            call attrib(idnc,idim,3,'psoil1',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon slow pool'
            call attrib(idnc,idim,3,'csoil2',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen slow pool'
            call attrib(idnc,idim,3,'nsoil2',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor slow pool'
            call attrib(idnc,idim,3,'psoil2',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Carbon pass pool'
            call attrib(idnc,idim,3,'csoil3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Nitrogen pass pool'
            call attrib(idnc,idim,3,'nsoil3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Phosphor pass pool'
            call attrib(idnc,idim,3,'psoil3',lname,'gC/m2',0.,6500.,
     &                  0,itype)
            lname = 'Prognostic LAI'
            call attrib(idnc,idim,3,'glai',lname,'none',0.,
     &                  13.,0,itype)
           end if
          end if
          if (nextout>=1.and.itype/=-1) then
            lname = 'Avg Net CO2 flux'
            call attrib(idnc,idim,3,'fnee_ave',lname,'gC/m2/s',-3.25E-3,
     &                  3.25E-3,0,itype)
            lname = 'Avg Photosynthesis CO2 flux'
            call attrib(idnc,idim,3,'fpn_ave',lname,'gC/m2/s',-3.25E-3,
     &                  3.25E-3,0,itype)
            lname = 'Avg Plant respiration CO2 flux'
            call attrib(idnc,idim,3,'frp_ave',lname,'gC/m2/s',-3.25E-3,
     &                  3.25E-3,0,itype)
            lname = 'Avg Soil respiration CO2 flux'
            call attrib(idnc,idim,3,'frs_ave',lname,'gC/m2/s',-3.25E-3,
     &                  3.25E-3,0,itype)
          end if
        end if

        ! URBAN -----------------------------------------------------
        if (nurban<=-1.or.(nurban>=1.and.itype==-1)) then
         lname = 'roof temperature lev 1'
         call attrib(idnc,idim,3,'rooftgg1',lname,'K',100.,425.,0,itype)
         lname = 'roof temperature lev 2'
         call attrib(idnc,idim,3,'rooftgg2',lname,'K',100.,425.,0,itype)
         lname = 'roof temperature lev 3'
         call attrib(idnc,idim,3,'rooftgg3',lname,'K',100.,425.,0,itype)
         lname = 'east wall temperature lev 1'
         call attrib(idnc,idim,3,'waletgg1',lname,'K',100.,425.,0,itype)
         lname = 'east wall temperature lev 2'
         call attrib(idnc,idim,3,'waletgg2',lname,'K',100.,425.,0,itype)
         lname = 'east wall temperature lev 3'
         call attrib(idnc,idim,3,'waletgg3',lname,'K',100.,425.,0,itype)
         lname = 'west wall temperature lev 1'
         call attrib(idnc,idim,3,'walwtgg1',lname,'K',100.,425.,0,itype)
         lname = 'west wall temperature lev 2'
         call attrib(idnc,idim,3,'walwtgg2',lname,'K',100.,425.,0,itype)
         lname = 'west wall temperature lev 3'
         call attrib(idnc,idim,3,'walwtgg3',lname,'K',100.,425.,0,itype)
         lname = 'road temperature lev 1'
         call attrib(idnc,idim,3,'roadtgg1',lname,'K',100.,425.,0,itype)
         lname = 'road temperature lev 2'
         call attrib(idnc,idim,3,'roadtgg2',lname,'K',100.,425.,0,itype)
         lname = 'road temperature lev 3'
         call attrib(idnc,idim,3,'roadtgg3',lname,'K',100.,425.,0,itype)
         lname = 'urban canyon soil moisture'
         call attrib(idnc,idim,3,'urbnsmc',lname,'m3/m3',0.,1.3,0,itype)
         lname = 'urban roof soil moisture'
         call attrib(idnc,idim,3,'urbnsmr',lname,'m3/m3',0.,1.3,0,itype)
         lname = 'urban roof water'
         call attrib(idnc,idim,3,'roofwtr',lname,'mm',0.,1.3,0,itype)
         lname = 'urban road water'
         call attrib(idnc,idim,3,'roadwtr',lname,'mm',0.,1.3,0,itype)
         lname = 'urban canyon leaf water'
         call attrib(idnc,idim,3,'urbwtrc',lname,'mm',0.,1.3,0,itype)
         lname = 'urban roof leaf water'
         call attrib(idnc,idim,3,'urbwtrr',lname,'mm',0.,1.3,0,itype)
         lname = 'urban roof snow'
         call attrib(idnc,idim,3,'roofsnd',lname,'mm',0.,1.3,0,itype)
         lname = 'urban road snow'
         call attrib(idnc,idim,3,'roadsnd',lname,'mm',0.,1.3,0,itype)
         lname = 'urban roof snow density'
         call attrib(idnc,idim,3,'roofden',lname,'kg/m3',0.,650.,0,
     &               itype)
         lname = 'urban road snow density'
         call attrib(idnc,idim,3,'roadden',lname,'kg/m3',0.,650.,0,
     &               itype)
         lname = 'urban roof snow albedo'
         call attrib(idnc,idim,3,'roofsna',lname,'none',0.,1.3,0,itype)
         lname = 'urban road snow albedo'
         call attrib(idnc,idim,3,'roadsna',lname,'none',0.,1.3,0,itype)
        end if
        
        ! STANDARD 3D VARIABLES -------------------------------------
        if (myid==0) then
          write(6,*) '3d variables'
        end if
        if(nextout>=4.and.nllp==3)then   ! N.B. use nscrn=1 for hourly output
          lname = 'Delta latitude'
          call attrib(idnc,dim,4,'del_lat',lname,'deg',-60.,60.,1,itype)
          lname = 'Delta longitude'
          call attrib(idnc,dim,4,'del_lon',lname,'deg',-180.,180.,1,
     &                itype)
          lname = 'Delta pressure'
          call attrib(idnc,dim,4,'del_p',lname,'hPa',-900.,900.,1,itype)
        endif  ! (nextout>=4.and.nllp==3)
        call attrib(idnc,dim,4,'temp','Air temperature','K',100.,350.,
     &              0,itype)
        lname= 'x-component wind'
        call attrib(idnc,dim,4,'u',lname,'m/s',-150.,150.,0,itype)
        lname= 'y-component wind'
        call attrib(idnc,dim,4,'v',lname,'m/s',-150.,150.,0,itype)
        lname= 'vertical velocity'
        call attrib(idnc,dim,4,'omega',lname,'Pa/s',-50.,50.,0,itype)
        lname= 'Water mixing ratio'
        call attrib(idnc,dim,4,'mixr',lname,'kg/kg',0.,.05,0,itype)
        
        ! MICROPHYSICS ----------------------------------------------
        if(ldr/=0)then
         call attrib(idnc,dim,4,'qfg','Frozen water','kg/kg',0.,.065,
     &               0,itype)
         call attrib(idnc,dim,4,'qlg','Liquid water','kg/kg',0.,.065,
     &               0,itype)
         call attrib(idnc,dim,4,'qrg','Rain','kg/kg',0.,.065,
     &               0,itype)
    !     call attrib(idnc,dim,4,'qsg','Snow','kg/kg',0.,.02,
    ! &               0,itype)
    !     call attrib(idnc,dim,4,'qgrau','Grauple','kg/kg',0.,.02,
    ! &               0,itype)
         call attrib(idnc,dim,4,'cfrac','Cloud fraction','none',0.,1.,
     &               0,itype)
         call attrib(idnc,dim,4,'cfrain','Rain fraction','none',0.,
     &               1.,0,itype)          
        endif
        
        ! TURBULENT MIXING ------------------------------------------
        if (nvmix==6.and.(nextout>=1.or.itype==-1))then
         call attrib(idnc,dim,4,'tke','Turbulent Kinetic Energy'
     &              ,'m2/s2',0.,65.,0,itype)
         call attrib(idnc,dim,4,'eps','Eddy dissipation rate'
     &              ,'m2/s3',0.,6.5,0,itype)
        end if

        ! TRACER ----------------------------------------------------
        if (ngas>0) then 
         do igas=1,ngas
           write(trnum,'(i3.3)') igas
!          rml 18/09/07 use tracmax from tracer.dat as previous formula
!                       wasn't always reliable
           trmax = tracmax(igas)
           trmin = tracmin(igas)
!          rml 19/09/07 use tracname as part of tracer long name
           lname = 'Tracer (inst.) '//trim(tracname(igas))
           call attrib(idnc,dim,4,'tr'//trnum,lname,'ppm',trmin,trmax,
     &                 0,-1) ! -1 = long
           lname = 'Tracer (average) '//trim(tracname(igas))
           call attrib(idnc,dim,4,'trav'//trnum,lname,'ppm',trmin,trmax
     &                 ,0,-1) ! -1 = long
!          rml 14/5/10 option to write out local time afternoon averages
           if (writetrpm)
     &     call attrib(idnc,dim,4,'trpm'//trnum,lname,'ppm',trmin,trmax
     &                 ,0,-1) ! -1 = long
         enddo ! igas loop
        endif  ! (ntrac.gt.0)

        ! AEROSOL ---------------------------------------------------
        if (iaero<=-2.or.(iaero>=2.and.itype==-1)) then  
          call attrib(idnc,dim,4,'dms','Dimethyl sulfide'
     &              ,'kg/kg',0.,6.5E-7,0,itype)
          call attrib(idnc,dim,4,'so2','Sulfur dioxide'
     &              ,'kg/kg',0.,6.5E-7,0,itype)
          call attrib(idnc,dim,4,'so4','Sulfate'
     &              ,'kg/kg',0.,6.5E-7,0,itype)
          call attrib(idnc,dim,4,'bco','Black carbon hydrophobic'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'bci','Black carbon hydrophilic'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'oco','Organic aerosol hydrophobic'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'oci','Organic aerosol hydrophilic'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'dust1','Dust 0.1-1 micrometers'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'dust2','Dust 1-2 micrometers'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'dust3','Dust 2-3 micrometers'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'dust4','Dust 3-6 micrometers'
     &              ,'kg/kg',0.,6.5E-6,0,itype)
          call attrib(idnc,dim,4,'seasalt1','Sea salt small'
     &              ,'1/m3',0.,6.5E9,0,itype)
          call attrib(idnc,dim,4,'seasalt2','Sea salt large'
     &              ,'1/m3',0.,6.5E7,0,itype)
        end if

        ! RESTART ---------------------------------------------------
        if(itype==-1)then   ! extra stuff just written for restart file
         lname= 'NHS adjustment to geopotential height'
         call attrib(idnc,dim,4,'zgnhs',lname,'m2/s2',-6.E5,6.E5,
     &               0,itype)
         lname= 'sdot: change in grid spacing per time step +.5'
         call attrib(idnc,dim,4,'sdot',lname,'1/ts',-3.,3.,0,itype) 
         lname= 'pslx: advective time rate of change of psl'
         call attrib(idnc,dim,4,'pslx',lname,'1/s',-1.E-3,1.E-3,0,
     &               itype)
         lname= 'savu'
         call attrib(idnc,dim,4,'savu',lname,'m/s',-1.E2,1.E2,0,
     &               itype)
         lname= 'savv'
         call attrib(idnc,dim,4,'savv',lname,'m/s',-1.E2,1.E2,0,
     &               itype)
         lname= 'savu1'
         call attrib(idnc,dim,4,'savu1',lname,'m/s',-1.E2,1.E2,0,
     &               itype)
         lname= 'savv1'
         call attrib(idnc,dim,4,'savv1',lname,'m/s',-1.E2,1.E2,0,
     &               itype)
         lname= 'savu2'
         call attrib(idnc,dim,4,'savu2',lname,'m/s',-1.E2,1.E2,0,
     &               itype)
         lname= 'savv2'
         call attrib(idnc,dim,4,'savv2',lname,'m/s',-1.E2,1.E2,0,
     &               itype)
         if (abs(nmlo)>=3.and.abs(nmlo)<=9) then
           do k=1,wlev
             write(lname,'("oldu1 ",I2)') k
             write(vname,'("oldu1",I2.2)') k
             call attrib(idnc,idim,3,vname,lname,'m/s',-100.,100.,0,
     &                   itype)
             write(lname,'("oldv1 ",I2)') k
             write(vname,'("oldv1",I2.2)') k
             call attrib(idnc,idim,3,vname,lname,'m/s',-100.,100.,0,
     &                   itype)
             write(lname,'("oldu2 ",I2)') k
             write(vname,'("oldu2",I2.2)') k
             call attrib(idnc,idim,3,vname,lname,'m/s',-100.,100.,0,
     &                   itype)
             write(lname,'("oldv2 ",I2)') k
             write(vname,'("oldv2",I2.2)') k
             call attrib(idnc,idim,3,vname,lname,'m/s',-100.,100.,0,
     &                   itype)
           end do
           lname= 'ipice'
           call attrib(idnc,idim,3,'ipice',lname,'Pa',0.,1.E6,0,
     &                 itype)
         end if
         lname = 'Soil ice lev 1'
         call attrib(idnc,idim,3,'wbice1',lname,'m3/m3',0.,1.,0,itype)
         lname = 'Soil ice lev 2'
         call attrib(idnc,idim,3,'wbice2',lname,'m3/m3',0.,1.,0,itype)
         lname = 'Soil ice lev 3'
         call attrib(idnc,idim,3,'wbice3',lname,'m3/m3',0.,1.,0,itype)
         lname = 'Soil ice lev 4'
         call attrib(idnc,idim,3,'wbice4',lname,'m3/m3',0.,1.,0,itype)
         lname = 'Soil ice lev 5'
         call attrib(idnc,idim,3,'wbice5',lname,'m3/m3',0.,1.,0,itype)
         lname = 'Soil ice lev 6'
         call attrib(idnc,idim,3,'wbice6',lname,'m3/m3',0.,1.,0,itype)
         if (nmlo==0) then ! otherwise already defined above
           lname = 'Snow temperature lev 1'
           call attrib(idnc,idim,3,'tggsn1',lname,'K',100.,425.,0,itype)
           lname = 'Snow temperature lev 2'
           call attrib(idnc,idim,3,'tggsn2',lname,'K',100.,425.,0,itype)
           lname = 'Snow temperature lev 3'
           call attrib(idnc,idim,3,'tggsn3',lname,'K',100.,425.,0,itype)
         end if
         lname = 'Snow mass lev 1'
         call attrib(idnc,idim,3,'smass1',lname,'K',0.,425.,0,itype)
         lname = 'Snow mass lev 2'
         call attrib(idnc,idim,3,'smass2',lname,'K',0.,425.,0,itype)
         lname = 'Snow mass lev 3'
         call attrib(idnc,idim,3,'smass3',lname,'K',0.,425.,0,itype)
         lname = 'Snow density lev 1'
         call attrib(idnc,idim,3,'ssdn1',lname,'K',0.,425.,0,itype)
         lname = 'Snow density lev 2'
         call attrib(idnc,idim,3,'ssdn2',lname,'K',0.,425.,0,itype)
         lname = 'Snow density lev 3'
         call attrib(idnc,idim,3,'ssdn3',lname,'K',0.,425.,0,itype)
         lname = 'Snow age'
         call attrib(idnc,idim,3,'snage',lname,'none',0.,20.,0,itype)   
         lname = 'Snow flag'
         call attrib(idnc,idim,3,'sflag',lname,'none',0.,4.,0,itype)
         lname = 'Solar net at ground (+ve down)'
         call attrib(idnc,idim,3,'sgsave',lname,'W/m2',-500.,2000.,0,
     &               itype)
         if (nsib==6.or.nsib==7) then
          call savetiledef(idnc,local,idim)
         end if
        endif  ! (itype==-1)
        
        if (myid==0) then
          write(6,*) 'finished defining attributes'
        end if
c       Leave define mode
        call ccnf_enddef(idnc)
        if (myid==0) then
          write(6,*) 'leave define mode'
        end if

        if(local)then
           ! Set these to global indices (relative to panel 0 in uniform decomp)
           do i=1,ipan
              xpnt(i) = float(i) + ioff(0)
           end do
           iduma(1)=1
           iduma(2)=il
           call ccnf_put_vara(idnc,ixp,iduma(1:1),iduma(2:2),
     &                             xpnt(1:il))
           i=1
           do n=1,npan
             do j=1,jpan
               ypnt(i) = float(j) + joff(n-noff) + (n-noff)*il_g
               i=i+1
             end do
           end do
           iduma(1)=1
           iduma(2)=jl
           call ccnf_put_vara(idnc,iyp,iduma(1:1),iduma(2:2),
     &                             ypnt(1:jl))
        else
           do i=1,il_g
              xpnt(i) = float(i)
           end do
           iduma(1)=1
           iduma(2)=il_g
           call ccnf_put_vara(idnc,ixp,iduma(1:1),iduma(2:2),xpnt)
           do j=1,jl_g
              ypnt(j) = float(j)
           end do
           iduma(1)=1
           iduma(2)=jl_g
           call ccnf_put_vara(idnc,iyp,iduma(1:1),iduma(2:2),ypnt)
        endif

        iduma(1)=1
        iduma(2)=kl
        call ccnf_put_vara(idnc,idlev,iduma(1:1),iduma(2:2),sig)

        call ccnf_inq_varid(idnc,'sigma',idv,tst)
        iduma(1)=1
        iduma(2)=kl
        call ccnf_put_vara(idnc,idv,iduma(1:1),iduma(2:2),sig)

        call ccnf_inq_varid(idnc,'lev',idv,tst)
        iduma(1)=1
        iduma(2)=kl
        call ccnf_put_vara(idnc,idv,iduma(1:1),iduma(2:2),sig)

        zsoil(1)=.5*zse(1)
        zsoil(2)=zse(1)+zse(2)*.5
        zsoil(3)=zse(1)+zse(2)+zse(3)*.5
        zsoil(4)=zse(1)+zse(2)+zse(3)+zse(4)*.5
        zsoil(5)=zse(1)+zse(2)+zse(3)+zse(4)+zse(5)*.5
        zsoil(6)=zse(1)+zse(2)+zse(3)+zse(4)+zse(5)+zse(6)*.5
        iduma(1)=1
        iduma(2)=ms
        call ccnf_put_vara(idnc,idms,iduma(1:1),iduma(2:2),zsoil)
        
        if (abs(nmlo)>0.and.abs(nmlo)<=9) then
          iduma(1)=1
          iduma(2)=wlev
          call ccnf_put_vara(idnc,idoc,iduma(1:1),iduma(2:2),
     &                            gosig)
        end if

        call ccnf_inq_varid(idnc,'ds',idv,tst)
        call ccnf_put_var1(idnc,idv,1,ds)
        call ccnf_inq_varid(idnc,'dt',idv,tst)
        call ccnf_put_var1(idnc,idv,1,dt)
       endif ! iarch==1
!      -----------------------------------------------------------      
       if (myid==0) then
         write(6,*) 'outcdf processing kdate,ktime,ktau,mtimer: ',
     &                                 kdate,ktime,ktau,mtimer
       end if
c      set time to number of minutes since start 
       call ccnf_inq_varid(idnc,'time',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,real(mtimer))
       call ccnf_inq_varid(idnc,'timer',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,timer)
       call ccnf_inq_varid(idnc,'mtimer',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,mtimer)
       call ccnf_inq_varid(idnc,'timeg',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,timeg)
       call ccnf_inq_varid(idnc,'ktau',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,ktau)
       call ccnf_inq_varid(idnc,'kdate',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,kdate)
       call ccnf_inq_varid(idnc,'ktime',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,ktime)
       call ccnf_inq_varid(idnc,'nstag',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,nstag)
       call ccnf_inq_varid(idnc,'nstagu',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,nstagu)
       idum=mod(ktau-nstagoff,max(abs(nstagin),1))
       idum=idum-max(abs(nstagin),1) ! should be -ve
       call ccnf_inq_varid(idnc,'nstagoff',idv,tst)
       call ccnf_put_var1(idnc,idv,iarch,idum)
       if ((nmlo<0.and.nmlo>=-9).or.
     &      (nmlo>0.and.nmlo<=9.and.itype==-1)) then
         idum=mod(ktau-nstagoffmlo,max(2*mstagf,1))
         idum=idum-max(2*mstagf,1) ! should be -ve
         call ccnf_inq_varid(idnc,'nstagoffmlo',idv,tst)
         call ccnf_put_var1(idnc,idv,iarch,idum)
       end if
       if (myid==0) then
         write(6,*) 'kdate,ktime,ktau=',kdate,ktime,ktau
         write(6,*) 'timer,timeg=',timer,timeg
         write(6,*) 'now write out variables'
       end if
       
      endif ! myid == 0 .or. local

      ! Export ocean data
      if (nmlo/=0.and.abs(nmlo)<=9) then
        mlodwn(:,:,1:2)=999.
        mlodwn(:,:,3:4)=0.
        micdwn=999.
        micdwn(:,9)=0.
        micdwn(:,10)=0.
        ocndep=0. ! ocean depth
        ocnheight=0. ! free surface height
        call mlosave(mlodwn,ocndep,ocnheight,micdwn,0)
      end if        

      !**************************************************************
      ! WRITE TIME-INVARIANT VARIABLES
      !**************************************************************

      if(ktau==0.or.itype==-1)then  ! also for restart file
        call histwrt3(zs,'zht',idnc,iarch,local,.true.)
        call histwrt3(he,'he',idnc,iarch,local,.true.)
        call histwrt3(em,'map',idnc,iarch,local,.true.)
        call histwrt3(f,'cor',idnc,iarch,local,.true.)
        call histwrt3(sigmu,'sigmu',idnc,iarch,local,.true.)
        aa(:)=real(isoilm(:))
        call histwrt3(aa,'soilt',idnc,iarch,local,.true.)
        aa(:)=real(ivegt(:))
        call histwrt3(aa,'vegt',idnc,iarch,local,.true.)
        if ((nmlo<0.and.nmlo>=-9).or.
     &      (nmlo>0.and.nmlo<=9.and.itype==-1)) then
          call histwrt3(ocndep,'ocndepth',idnc,iarch,local,.true.)
        end if
      endif ! (ktau==0.or.itype==-1) 

      !**************************************************************
      ! WRITE 3D VARIABLES (2D + Time)
      !**************************************************************

      ! BASIC -------------------------------------------------------
      lwrite=ktau>0
      if (nsib==6.or.nsib==7) then
        call histwrt3(rsmin,'rs',idnc,iarch,local,lwrite)
      else if (ktau==0.or.itype==-1) then
        call histwrt3(rsmin,'rsmin',idnc,iarch,local,.true.)
      end if
      call histwrt3(sigmf,'sigmf',idnc,iarch,local,.true.)
      call histwrt3(psl,'psf',idnc,iarch,local,.true.)
      tmpry=t(1:ifull,:)
      call mslp(aa,psl,zs,tmpry)
      aa=aa/100.
      call histwrt3(aa,'pmsl',idnc,iarch,local,.true.)      
      call histwrt3(zo,'zolnd',idnc,iarch,local,.true.)
      call histwrt3(vlai,'lai',idnc,iarch,local,.true.)
      call histwrt3(tss,'tsu',idnc,iarch,local,.true.)
      call histwrt3(tpan,'tpan',idnc,iarch,local,.true.)
!       scale up precip,precc,sno,runoff to mm/day (soon reset to 0 in globpe)
!       ktau in next line in case ntau (& thus ktau) < nwt 
      aa(:)=precip(1:ifull)*real(nperday)
     &     /real(min(nwt,max(ktau,1))) 
      call histwrt3(aa,'rnd',idnc,iarch,local,lwrite)
      aa(:)=precc(1:ifull)*real(nperday)
     &     /real(min(nwt,max(ktau,1)))
      call histwrt3(aa,'rnc',idnc,iarch,local,lwrite)
      aa(:)=sno(1:ifull)*real(nperday)
     &     /real(min(nwt,max(ktau,1)))
      call histwrt3(aa,'sno',idnc,iarch,local,lwrite)
      aa(:)=runoff(1:ifull)*real(nperday)
     &     /real(min(nwt,max(ktau,1)))
      call histwrt3(aa,'runoff',idnc,iarch,local,lwrite)
      aa(:)=swrsave*albvisnir(:,1)+(1.-swrsave)*albvisnir(:,2)
      call histwrt3(aa,'alb',idnc,iarch,local,.true.)
      call histwrt3(fwet,'fwet',idnc,iarch,local,lwrite)

      ! MLO ---------------------------------------------------------      
      if (nmlo/=0) then
        ocnheight=min(max(ocnheight,-130.),130.)
        do k=1,ms
          where (.not.land)
            tgg(:,k)=mlodwn(:,k,1)
          end where
        end do
        do k=1,3
          where (.not.land)
            tggsn(:,k)=micdwn(:,k)
          end where
        end do
        where (.not.land)
          fracice=micdwn(:,5)
          sicedep=micdwn(:,6)
          snowd=micdwn(:,7)*1000.
        end where
      end if

      call histwrt3(snowd,'snd',idnc,iarch,local,.true.)  ! long write
      call histwrt3(tgg(1,1),'tgg1',idnc,iarch,local,.true.)
      call histwrt3(tgg(1,2),'tgg2',idnc,iarch,local,.true.)
      call histwrt3(tgg(1,3),'tgg3',idnc,iarch,local,.true.)
      call histwrt3(tgg(1,4),'tgg4',idnc,iarch,local,.true.)
      call histwrt3(tgg(1,5),'tgg5',idnc,iarch,local,.true.)
      call histwrt3(tgg(1,6),'tgg6',idnc,iarch,local,.true.)
      
      if (abs(nmlo)<=9) then
       if (nmlo<0.or.(nmlo>0.and.itype==-1)) then
        do k=ms+1,wlev
          write(vname,'("tgg",I2.2)') k
          call histwrt3(mlodwn(:,k,1),vname,idnc,iarch,local,.true.)
        end do
        do k=1,wlev
          write(vname,'("sal",I2.2)') k
          call histwrt3(mlodwn(:,k,2),vname,idnc,iarch,local,.true.)
        end do
        do k=1,wlev
          write(vname,'("uoc",I2.2)') k
          call histwrt3(mlodwn(:,k,3),vname,idnc,iarch,local,.true.)
          write(vname,'("voc",I2.2)') k
          call histwrt3(mlodwn(:,k,4),vname,idnc,iarch,local,.true.)
        end do
        call histwrt3(ocnheight,'ocheight',idnc,iarch,local,.true.)
        call histwrt3(tggsn(:,1),'tggsn1',idnc,iarch,local,.true.)
        call histwrt3(tggsn(:,2),'tggsn2',idnc,iarch,local,.true.)
        call histwrt3(tggsn(:,3),'tggsn3',idnc,iarch,local,.true.)
        call histwrt3(micdwn(:,4),'tggsn4',idnc,iarch,local,.true.)
        call histwrt3(micdwn(:,8),'sto',idnc,iarch,local,.true.)
        call histwrt3(micdwn(:,9),'uic',idnc,iarch,local,.true.)
        call histwrt3(micdwn(:,10),'vic',idnc,iarch,local,.true.)
        call histwrt3(micdwn(:,11),'icesal',idnc,iarch,local,.true.)
        if (abs(nmlo)>=2) then
          call histwrt3(watbdy(1:ifull),'swater',idnc,iarch,local,
     &                  .true.)
          call histwrt3(salbdy(1:ifull),'ssalin',idnc,iarch,local,
     &                  .true.)
        end if
       end if
      end if

      ! SOIL --------------------------------------------------------
      !call histwrt3(wb(1,1),'wb1',idnc,iarch,local,.true.)
      !call histwrt3(wb(1,2),'wb2',idnc,iarch,local,.true.)
      !call histwrt3(wb(1,3),'wb3',idnc,iarch,local,.true.)
      !call histwrt3(wb(1,4),'wb4',idnc,iarch,local,.true.)
      !call histwrt3(wb(1,5),'wb5',idnc,iarch,local,.true.)
      !call histwrt3(wb(1,6),'wb6',idnc,iarch,local,.true.)
      aa(:)=(wb(:,1)-swilt(isoilm))/
     &      (sfc(isoilm)-swilt(isoilm))
      call histwrt3(aa,'wetfrac1',idnc,iarch,local,.true.)
      aa(:)=(wb(:,2)-swilt(isoilm))/
     &      (sfc(isoilm)-swilt(isoilm))
      call histwrt3(aa,'wetfrac2',idnc,iarch,local,.true.)
      aa(:)=(wb(:,3)-swilt(isoilm))/
     &      (sfc(isoilm)-swilt(isoilm))
      call histwrt3(aa,'wetfrac3',idnc,iarch,local,.true.)
      aa(:)=(wb(:,4)-swilt(isoilm))/
     &      (sfc(isoilm)-swilt(isoilm))
      call histwrt3(aa,'wetfrac4',idnc,iarch,local,.true.)
      aa(:)=(wb(:,5)-swilt(isoilm))/
     &      (sfc(isoilm)-swilt(isoilm))
      call histwrt3(aa,'wetfrac5',idnc,iarch,local,.true.)
      aa(:)=(wb(:,6)-swilt(isoilm))/
     &      (sfc(isoilm)-swilt(isoilm))
      call histwrt3(aa,'wetfrac6',idnc,iarch,local,.true.)
      
      ! PH - Add wetfac to output for mbase=-19 option
      call histwrt3(wetfac,'wetfac',idnc,iarch,local,.true.)
      
      ! SEAICE ------------------------------------------------------       
      call histwrt3(sicedep,'siced',idnc,iarch,local,.true.)
      call histwrt3(fracice,'fracice',idnc,iarch,local,.true.)
      
      ! DIAGNOSTICS -------------------------------------------------
      lwrite=ktau>0
      call histwrt3(u10,'u10',idnc,iarch,local,.true.)
      
      if(itype/=-1)then  ! these not written to restart file
       aa=rndmax(:)*86400./dt ! scale up to mm/day
       call histwrt3(aa,'maxrnd',idnc,iarch,local,lday)
       call histwrt3(tmaxscr,'tmaxscr',idnc,iarch,local,lday)
       call histwrt3(tminscr,'tminscr',idnc,iarch,local,lday)
       call histwrt3(rhmaxscr,'rhmaxscr',idnc,iarch,local,lday)
       call histwrt3(rhminscr,'rhminscr',idnc,iarch,local,lday)
       call histwrt3(capemax,'capemax',idnc,iarch,local,lday)
       call histwrt3(u10max,'u10max',idnc,iarch,local,lday)
       call histwrt3(v10max,'v10max',idnc,iarch,local,lday)
       call histwrt3(u1max,'u1max',idnc,iarch,local,lday)
       call histwrt3(v1max,'v1max',idnc,iarch,local,lday)
       call histwrt3(u2max,'u2max',idnc,iarch,local,lday)
       call histwrt3(v2max,'v2max',idnc,iarch,local,lday)
!      if writes done more than once per day, 
!      needed to augment accumulated 3-hourly rainfall in rnd06 to rnd21 
!      to allow for intermediate zeroing of precip()
!      but not needed from 17/9/03 with introduction of rnd24
       call histwrt3(rnd_3hr(1,1),'rnd03',idnc,iarch,local,lday)
       call histwrt3(rnd_3hr(1,2),'rnd06',idnc,iarch,local,lday)
       call histwrt3(rnd_3hr(1,3),'rnd09',idnc,iarch,local,lday)
       call histwrt3(rnd_3hr(1,4),'rnd12',idnc,iarch,local,lday)
       call histwrt3(rnd_3hr(1,5),'rnd15',idnc,iarch,local,lday)
       call histwrt3(rnd_3hr(1,6),'rnd18',idnc,iarch,local,lday)
       call histwrt3(rnd_3hr(1,7),'rnd21',idnc,iarch,local,lday)
       call histwrt3(rnd_3hr(1,8),'rnd24',idnc,iarch,local,lday)
       if(nextout>=2) then ! 6-hourly u10 & v10
         call histwrt3( u10_3hr(1,2), 'u10_06',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,2), 'v10_06',idnc,iarch,local,
     &                  lday)
         call histwrt3( u10_3hr(1,4), 'u10_12',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,4), 'v10_12',idnc,iarch,local,
     &                  lday)
         call histwrt3( u10_3hr(1,6), 'u10_18',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,6), 'v10_18',idnc,iarch,local,
     &                  lday)
         call histwrt3( u10_3hr(1,8), 'u10_24',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,8), 'v10_24',idnc,iarch,local,
     &                  lday)
         call histwrt3(tscr_3hr(1,2),'tscr_06',idnc,iarch,local,
     &                  lday)
         call histwrt3(tscr_3hr(1,4),'tscr_12',idnc,iarch,local,
     &                  lday)
         call histwrt3(tscr_3hr(1,6),'tscr_18',idnc,iarch,local,
     &                  lday)
         call histwrt3(tscr_3hr(1,8),'tscr_24',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,2), 'rh1_06',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,4), 'rh1_12',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,6), 'rh1_18',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,8), 'rh1_24',idnc,iarch,local,
     &                  lday)
       endif  ! (nextout>=2)
       if(nextout>=3) then  ! also 3-hourly u10 & v10
         call histwrt3(tscr_3hr(1,1),'tscr_03',idnc,iarch,local,
     &                  lday)
         call histwrt3(tscr_3hr(1,3),'tscr_09',idnc,iarch,local,
     &                  lday)
         call histwrt3(tscr_3hr(1,5),'tscr_15',idnc,iarch,local,
     &                  lday)
         call histwrt3(tscr_3hr(1,7),'tscr_21',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,1), 'rh1_03',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,3), 'rh1_09',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,5), 'rh1_15',idnc,iarch,local,
     &                  lday)
         call histwrt3( rh1_3hr(1,7), 'rh1_21',idnc,iarch,local,
     &                  lday)
         call histwrt3( u10_3hr(1,1), 'u10_03',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,1), 'v10_03',idnc,iarch,local,
     &                  lday)
         call histwrt3( u10_3hr(1,3), 'u10_09',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,3), 'v10_09',idnc,iarch,local,
     &                  lday)
         call histwrt3( u10_3hr(1,5), 'u10_15',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,5), 'v10_15',idnc,iarch,local,
     &                  lday)
         call histwrt3( u10_3hr(1,7), 'u10_21',idnc,iarch,local,
     &                  lday)
         call histwrt3( v10_3hr(1,7), 'v10_21',idnc,iarch,local,
     &                  lday)
       endif  ! nextout>=3
       if(nextout>=4.and.nllp==3) then  
        do k=1,klt
         do iq=1,ilt*jlt        
          tr(iq,k,ngas+1)=tr(iq,k,ngas+1)-rlatt(iq)*180./pi
          tr(iq,k,ngas+2)=tr(iq,k,ngas+2)-rlongg(iq)*180./pi
          if(tr(iq,k,ngas+2)>180.)
     &                       tr(iq,k,ngas+2)=tr(iq,k,ngas+2)-360.
          if(tr(iq,k,ngas+2)<-180.)
     &                       tr(iq,k,ngas+2)=tr(iq,k,ngas+2)+360.
          tr(iq,k,ngas+3)=tr(iq,k,ngas+3)-.01*ps(iq)*sig(k)  ! in hPa
         enddo
        enddo
!       N.B. does not yet properly handle across Grenwich Meridion
        tmpry=tr(1:ifull,:,ngas+1)	   
        call histwrt4(tmpry,'del_lat',idnc,iarch,local,
     &                .true.)
        tmpry=tr(1:ifull,:,ngas+2)
        call histwrt4(tmpry,'del_lon',idnc,iarch,local,
     &                .true.)
        tmpry=tr(1:ifull,:,ngas+3)
        call histwrt4(tmpry,'del_p',idnc,iarch,local,
     &                .true.)
       endif  ! (nextout>=4.and.nllp==3)
       ! only write these once per avg period
       call histwrt3(tscr_ave,'tscr_ave',idnc,iarch,local,lave)
       call histwrt3(cbas_ave,'cbas_ave',idnc,iarch,local,lave)
       call histwrt3(ctop_ave,'ctop_ave',idnc,iarch,local,lave)
       call histwrt3(dew_ave,'dew_ave',idnc,iarch,local,lave)
       call histwrt3(evap,'evap',idnc,iarch,local,lave)
       call histwrt3(epan_ave,'epan_ave',idnc,iarch,local,lave)
       call histwrt3(epot_ave,'epot_ave',idnc,iarch,local,lave)
       call histwrt3(eg_ave,'eg_ave',idnc,iarch,local,lave)
       call histwrt3(fg_ave,'fg_ave',idnc,iarch,local,lave)
       call histwrt3(rnet_ave,'rnet_ave',idnc,iarch,local,lave)
       call histwrt3(ga_ave,'ga_ave',idnc,iarch,local,lave)
       call histwrt3(riwp_ave,'iwp_ave',idnc,iarch,local,lave)
       call histwrt3(rlwp_ave,'lwp_ave',idnc,iarch,local,lave)
       call histwrt3(cll_ave,'cll',idnc,iarch,local,lrad)
       call histwrt3(clm_ave,'clm',idnc,iarch,local,lrad)
       call histwrt3(clh_ave,'clh',idnc,iarch,local,lrad)
       call histwrt3(cld_ave,'cld',idnc,iarch,local,lrad)
       call histwrt3(wb_ave(:,1),'wb1_ave',idnc,iarch,local,lave)
       call histwrt3(wb_ave(:,2),'wb2_ave',idnc,iarch,local,lave)
       call histwrt3(wb_ave(:,3),'wb3_ave',idnc,iarch,local,lave)
       call histwrt3(wb_ave(:,4),'wb4_ave',idnc,iarch,local,lave)
       call histwrt3(wb_ave(:,5),'wb5_ave',idnc,iarch,local,lave)
       call histwrt3(wb_ave(:,6),'wb6_ave',idnc,iarch,local,lave)
       !  call histwrt3(tgg_ave(:,1),'tgg1_ave',idnc,iarch,local,lave)
       !  call histwrt3(tgg_ave(:,2),'tgg2_ave',idnc,iarch,local,lave)
       !  call histwrt3(tgg_ave(:,3),'tgg3_ave',idnc,iarch,local,lave)
       !  call histwrt3(tgg_ave(:,4),'tgg4_ave',idnc,iarch,local,lave)
       !  call histwrt3(tgg_ave(:,5),'tgg5_ave',idnc,iarch,local,lave)
       !  call histwrt3(tgg_ave(:,6),'tgg6_ave',idnc,iarch,local,lave)
       call histwrt3(tsu_ave,'tsu_ave',idnc,iarch,local,lave)
       call histwrt3(alb_ave,'alb_ave',idnc,iarch,local,lrad)
       call histwrt3(psl_ave,'pmsl_ave',idnc,iarch,local,lave)
       if (nmlo/=0) then
         call histwrt3(mixdep_ave,'mixd_ave',idnc,iarch,local,lave)
       end if
       lwrite=ktau>0
       call histwrt3(tscrn,'tscrn',idnc,iarch,local,lwrite)
       call histwrt3(qgscrn,'qgscrn',idnc,iarch,local,lwrite)
       call histwrt3(rhscrn,'rhscrn',idnc,iarch,local,lwrite)
       call histwrt3(uscrn,'uscrn',idnc,iarch,local,lwrite)
       call histwrt3(rnet,'rnet',idnc,iarch,local,lwrite)
       call histwrt3(epan,'epan',idnc,iarch,local,lwrite)
       call histwrt3(eg,'eg',idnc,iarch,local,lwrite)
       call histwrt3(fg,'fg',idnc,iarch,local,lwrite)
       call histwrt3(taux,'taux',idnc,iarch,local,lwrite)
       call histwrt3(tauy,'tauy',idnc,iarch,local,lwrite)
c      "extra" outputs
       if(nextout>=1) then
         if(myid == 0 ) write(6,*) 'nextout, idnc: ',nextout,idnc
         call histwrt3(rtu_ave,'rtu_ave',idnc,iarch,local,lrad)
         call histwrt3(rtc_ave,'rtc_ave',idnc,iarch,local,lrad)
         call histwrt3(rgdn_ave,'rgdn_ave',idnc,iarch,local,lrad)
         call histwrt3(rgn_ave,'rgn_ave',idnc,iarch,local,lrad)
         call histwrt3(rgc_ave,'rgc_ave',idnc,iarch,local,lrad)
         call histwrt3(sint_ave,'sint_ave',idnc,iarch,local,lrad)
         call histwrt3(sot_ave,'sot_ave',idnc,iarch,local,lrad)
         call histwrt3(soc_ave,'soc_ave',idnc,iarch,local,lrad)
         call histwrt3(sgdn_ave,'sgdn_ave',idnc,iarch,local,lrad)
         call histwrt3(sgn_ave,'sgn_ave',idnc,iarch,local,lave)
         aa=sunhours/3600.
         call histwrt3(aa,'sunhours',idnc,iarch,local,lave)
         call histwrt3(fbeam_ave,'fbeam_ave',idnc,iarch,local,lrad)
         lwrite=ktau>0
         call histwrt3(dpsdt,'dpsdt',idnc,iarch,local,lwrite)
         call histwrt3(ustar,'ustar',idnc,iarch,local,lwrite)
       endif   ! nextout>=1
      endif    ! (ktau>0.and.itype/=-1)
      
      ! TURBULENT MIXING --------------------------------------------
      if (nextout>=1.or.(nvmix==6.and.itype==-1)) then
       call histwrt3(pblh,'pblh',idnc,iarch,local,.true.)
      end if

      ! CABLE -------------------------------------------------------
      if (nsib==6.or.nsib==7) then
        if (nextout>=1.or.itype==-1) then
         if (ccycle==0) then
          call histwrt3(cplant(:,1),'cplant1',idnc,iarch,local,.true.)
          call histwrt3(cplant(:,2),'cplant2',idnc,iarch,local,.true.)
          call histwrt3(cplant(:,3),'cplant3',idnc,iarch,local,.true.)
          call histwrt3(csoil(:,1),'csoil1',idnc,iarch,local,.true.)
          call histwrt3(csoil(:,2),'csoil2',idnc,iarch,local,.true.)
         else
          lwrite=mod(ktau,nperday)==0.or.ktau==ntau ! only write once per day
          do k=1,mplant
            write(vname,'("cplant",I1.1)') k
            call histwrt3(cplant(:,k),vname,idnc,iarch,local,lwrite)
            write(vname,'("nplant",I1.1)') k
            call histwrt3(niplant(:,k),vname,idnc,iarch,local,lwrite)
            write(vname,'("pplant",I1.1)') k
            call histwrt3(pplant(:,k),vname,idnc,iarch,local,lwrite)
          end do
          do k=1,mlitter
            write(vname,'("clitter",I1.1)') k
            call histwrt3(clitter(:,k),vname,idnc,iarch,local,lwrite)
            write(vname,'("nlitter",I1.1)') k
            call histwrt3(nilitter(:,k),vname,idnc,iarch,local,lwrite)
            write(vname,'("plitter",I1.1)') k
            call histwrt3(plitter(:,k),vname,idnc,iarch,local,lwrite)
          end do
          do k=1,msoil
            write(vname,'("csoil",I1.1)') k
            call histwrt3(csoil(:,k),vname,idnc,iarch,local,lwrite)
            write(vname,'("nsoil",I1.1)') k
            call histwrt3(nisoil(:,k),vname,idnc,iarch,local,lwrite)
            write(vname,'("psoil",I1.1)') k
            call histwrt3(psoil(:,k),vname,idnc,iarch,local,lwrite)
          end do
          call histwrt3(glai,'glai',idnc,iarch,local,lwrite)
         end if
        end if
        if (nextout>=1.and.itype/=-1) then
          aa=fpn_ave+frp_ave+frs_ave
          call histwrt3(aa,'fnee_ave',idnc,iarch,local,lave)
          call histwrt3(fpn_ave,'fpn_ave',idnc,iarch,local,lave)
          call histwrt3(frp_ave,'frp_ave',idnc,iarch,local,lave)
          call histwrt3(frs_ave,'frs_ave',idnc,iarch,local,lave)
        end if
      endif   

      ! URBAN -------------------------------------------------------
      if (nurban<=-1.or.(nurban>=1.and.itype==-1)) then
       atebdwn(:,:)=999. ! must be the same as spval in onthefly.f
       call atebsave(atebdwn,0)
       call histwrt3(atebdwn(:,1),'rooftgg1',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,2),'rooftgg2',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,3),'rooftgg3',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,4),'waletgg1',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,5),'waletgg2',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,6),'waletgg3',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,7),'walwtgg1',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,8),'walwtgg2',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,9),'walwtgg3',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,10),'roadtgg1',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,11),'roadtgg2',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,12),'roadtgg3',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,13),'urbnsmc',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,14),'urbnsmr',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,15),'roofwtr',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,16),'roadwtr',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,17),'urbwtrc',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,18),'urbwtrr',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,19),'roofsnd',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,20),'roadsnd',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,21),'roofden',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,22),'roadden',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,23),'roofsna',idnc,iarch,local,.true.)
       call histwrt3(atebdwn(:,24),'roadsna',idnc,iarch,local,.true.)
      end if

      ! **************************************************************
      ! WRITE 4D VARIABLES (3D + Time)
      ! **************************************************************

      ! ATMOSPHERE DYNAMICS ------------------------------------------
      if(myid == 0 ) write(6,*) 'netcdf save of 3d variables'
      lwrite=ktau>0
      tmpry=t(1:ifull,:)
      call histwrt4(tmpry,'temp',idnc,iarch,local,.true.)
      tmpry=u(1:ifull,:)
      call histwrt4(tmpry,'u',idnc,iarch,local,.true.)
      tmpry=v(1:ifull,:)
      call histwrt4(tmpry,'v',idnc,iarch,local,.true.)
      do k=1,kl
       do iq=1,ifull
        tmpry(iq,k)=ps(iq)*dpsldt(iq,k)
       enddo
      enddo
      call histwrt4(tmpry,'omega',idnc,iarch,local,lwrite)
      tmpry=qg(1:ifull,:)
      call histwrt4(tmpry,'mixr',idnc,iarch,local,.true.)
      
      ! MICROPHYSICS ------------------------------------------------
      if(ldr/=0)then
        tmpry=qfg(1:ifullw,:)
        call histwrt4(tmpry,'qfg',idnc,iarch,local,.true.)
        tmpry=qlg(1:ifullw,:)
        call histwrt4(tmpry,'qlg',idnc,iarch,local,.true.)
        tmpry=qrg(1:ifullw,:)
        call histwrt4(tmpry,'qrg',idnc,iarch,local,.true.)
        !call histwrt4(qsg(1:ifullw,:),'qsg',idnc,iarch,local,.true.)
        !call histwrt4(qgrau(1:ifullw,:),'qgrau',idnc,iarch,local,.true.)
        call histwrt4(cfrac,'cfrac',idnc,iarch,local,.true.)
        tmpry=cffall(1:ifullw,:)
        call histwrt4(tmpry,'cfrain',idnc,iarch,local,
     &                .true.)
      endif
      
      ! TURBULENT MIXING --------------------------------------------
      if (nvmix==6.and.(nextout>=1.or.itype==-1))then
        tmpry=tke(1:ifull,:)
        call histwrt4(tmpry,'tke',idnc,iarch,local,.true.)
        tmpry=eps(1:ifull,:)
        call histwrt4(tmpry,'eps',idnc,iarch,local,.true.)
      end if

      ! TRACERS -----------------------------------------------------
      if(ngas>0)then 
       do igas=1,ngas
        write(trnum,'(i3.3)') igas
        tmpry=tr(1:ilt*jlt,:,igas)+trback_g(igas)
        call histwrt4(tmpry,'tr'//trnum,
     &                idnc,iarch,local,.true.)
        tmpry=traver(:,:,igas)+trback_g(igas)
        call histwrt4(tmpry,'trav'//trnum,
     &                idnc,iarch,local,lave)
        ! rml 14/5/10 option to write out local time afternoon average
        if (writetrpm) then
          ! first divide by number of contributions to average
          do k=1,klt
            trpm(1:ifull,k,igas) = trpm(1:ifull,k,igas)/float(npm)
          enddo
          tmpry=trpm(:,:,igas)+trback_g(igas)
          call histwrt4(tmpry,'trpm'//trnum,
     &                  idnc,iarch,local,.true.)
        endif
       enddo ! igas loop
       ! reset arrays
       if (writetrpm) then
         trpm = 0.
         npm = 0.
       endif
      endif  ! (ngasc>0)

      ! AEROSOLS ----------------------------------------------------
      if (iaero<=-2.or.(iaero>=2.and.itype==-1)) then
        tmpry=xtg(1:ifull,:,1)
        call histwrt4(tmpry,'dms',idnc,iarch,local,.true.)
        tmpry=xtg(1:ifull,:,2)
        call histwrt4(tmpry,'so2',idnc,iarch,local,.true.)
        tmpry=xtg(1:ifull,:,3)
        call histwrt4(tmpry,'so4',idnc,iarch,local,.true.)
        tmpry=xtg(1:ifull,:,4)
        call histwrt4(tmpry,'bco',idnc,iarch,local,.true.)
        tmpry=xtg(1:ifull,:,5)
        call histwrt4(tmpry,'bci',idnc,iarch,local,.true.)
        tmpry=xtg(1:ifull,:,6)
        call histwrt4(tmpry,'oco',idnc,iarch,local,.true.)
        tmpry=xtg(1:ifull,:,7)
        call histwrt4(tmpry,'oci',idnc,iarch,local,.true.)
        tmpry=xtg(1:ifull,:,8)
        call histwrt4(tmpry,'dust1',idnc,iarch,local,
     &                .true.)
        tmpry=xtg(1:ifull,:,9)
        call histwrt4(tmpry,'dust2',idnc,iarch,local,
     &                .true.)
        tmpry=xtg(1:ifull,:,10)
        call histwrt4(tmpry,'dust3',idnc,iarch,local,
     &                .true.)
        tmpry=xtg(1:ifull,:,11)
        call histwrt4(tmpry,'dust4',idnc,iarch,local,
     &                .true.)
        tmpry=ssn(1:ifull,:,1)
        call histwrt4(tmpry,'seasalt1',idnc,iarch,local,
     &                .true.)
        tmpry=ssn(1:ifull,:,2)
        call histwrt4(tmpry,'seasalt2',idnc,iarch,local,
     &                .true.)
      end if

      !**************************************************************
      ! RESTART ONLY DATA
      !**************************************************************

      if(itype==-1)then
       call histwrt4(phi_nh,'zgnhs',idnc,iarch,local,.true.)
       call histwrt4(sdot(1,2),'sdot',idnc,iarch,local,.true.)
       tmpry=pslx(1:ifull,:)
       call histwrt4(tmpry,'pslx',idnc,iarch,local,.true.)
       call histwrt4(savu,'savu',idnc,iarch,local,.true.)
       call histwrt4(savv,'savv',idnc,iarch,local,.true.)
       call histwrt4(savu1,'savu1',idnc,iarch,local,.true.)
       call histwrt4(savv1,'savv1',idnc,iarch,local,.true.)
       call histwrt4(savu2,'savu2',idnc,iarch,local,.true.)
       call histwrt4(savv2,'savv2',idnc,iarch,local,.true.)
       if (abs(nmlo)>=3.and.abs(nmlo)<=9) then
         do k=1,wlev
           write(vname,'("oldu1",I2.2)') k
           call histwrt3(oldu1(:,k),vname,idnc,iarch,local,.true.)
           write(vname,'("oldv1",I2.2)') k
           call histwrt3(oldv1(:,k),vname,idnc,iarch,local,.true.)
           write(vname,'("oldu2",I2.2)') k
           call histwrt3(oldu2(:,k),vname,idnc,iarch,local,.true.)
           write(vname,'("oldv2",I2.2)') k
           call histwrt3(oldv2(:,k),vname,idnc,iarch,local,.true.)
         end do
         call histwrt3(ipice,'ipice',idnc,iarch,local,.true.)
       end if
       call histwrt3(wbice(1,1),'wbice1',idnc,iarch,local,.true.)
       call histwrt3(wbice(1,2),'wbice2',idnc,iarch,local,.true.)
       call histwrt3(wbice(1,3),'wbice3',idnc,iarch,local,.true.)
       call histwrt3(wbice(1,4),'wbice4',idnc,iarch,local,.true.)
       call histwrt3(wbice(1,5),'wbice5',idnc,iarch,local,.true.)
       call histwrt3(wbice(1,6),'wbice6',idnc,iarch,local,.true.)
       if (nmlo==0) then ! otherwise already written above
         call histwrt3(tggsn(1,1),'tggsn1',idnc,iarch,local,.true.)
         call histwrt3(tggsn(1,2),'tggsn2',idnc,iarch,local,.true.)
         call histwrt3(tggsn(1,3),'tggsn3',idnc,iarch,local,.true.)
       end if
       call histwrt3(smass(1,1),'smass1',idnc,iarch,local,.true.)
       call histwrt3(smass(1,2),'smass2',idnc,iarch,local,.true.)
       call histwrt3(smass(1,3),'smass3',idnc,iarch,local,.true.)
       call histwrt3(ssdn(1,1),'ssdn1',idnc,iarch,local,.true.)
       call histwrt3(ssdn(1,2),'ssdn2',idnc,iarch,local,.true.)
       call histwrt3(ssdn(1,3),'ssdn3',idnc,iarch,local,.true.)
       call histwrt3(snage,'snage',idnc,iarch,local,.true.)
       aa(:)=isflag(:)
       call histwrt3(aa,'sflag',idnc,iarch,local,.true.)
       call histwrt3(sgsave,'sgsave',idnc,iarch,local,.true.)       
       if (nsib==6.or.nsib==7) then
         call savetile(idnc,local,iarch)
       end if
      endif  ! (itype==-1)

      return
      end

      !--------------------------------------------------------------
      ! HIGH FREQUENCY OUTPUT FILES
      
      ! Here we buffer high frequency output.  However, it may be
      ! unnecessary if the netcdf library also buffers the output
      subroutine freqfile

      use arrays_m                          ! Atmosphere dyamics prognostic arrays
      use cc_mpi                            ! CC MPI routines
      use infile                            ! Input file routines
      use morepbl_m                         ! Additional boundary layer diagnostics
      use parmhdff_m                        ! Horizontal diffusion parameters
      use screen_m                          ! Screen level diagnostics
      use sigs_m                            ! Atmosphere sigma levels
      use tracers_m                         ! Tracer data
      
      implicit none

      include 'newmpar.h'                   ! Grid parameters
      include 'dates.h'                     ! Date data
      include 'filnames.h'                  ! Filenames
      include 'kuocom.h'                    ! Convection parameters
      include 'parm.h'                      ! Model configuration
      include 'parmdyn.h'                   ! Dynamics parameters
      include 'parmgeom.h'                  ! Coordinate data
      include 'parmhor.h'                   ! Horizontal advection parameters
      include 'parmvert.h'                  ! Vertical advection parameters

      integer leap
      common/leap_yr/leap                   ! Leap year (1 to allow leap years)
      
      integer, parameter :: freqvars = 5  ! number of variables to write
      integer, parameter :: nihead   = 54
      integer, parameter :: nrhead   = 14
      integer, dimension(nihead) :: nahead
      integer, dimension(nwt) :: datedat
      integer, dimension(4) :: adim
      integer, dimension(3) :: sdim
      integer, dimension(1) :: start,ncount
      integer, dimension(2) :: iduma
      integer ierr,ixp,iyp,izp,old_mode
      integer icy,icm,icd,ich,icmi,ics,ti
      integer i,j,n,fiarch
      integer, save :: fncid = -1
      integer, save :: idnt = 0
      integer, save :: idkdate = 0
      integer, save :: idktime = 0
      integer, save :: idmtimer = 0
      real, dimension(:,:,:), allocatable, save :: freqstore
      real, dimension(ifull) :: uas,vas,umag
      real, dimension(il_g) :: xpnt
      real, dimension(jl_g) :: ypnt
      real, dimension(nrhead) :: ahead
      double precision, dimension(nwt) :: tpnt
      logical, save :: first = .true.
      character(len=180) :: ffile
      character(len=40) :: lname
      character(len=33) :: grdtim
      character(len=20) :: timorg
      character(len=3), dimension(12) :: month

      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/      

      call start_log(outfile_begin)

      ! allocate arrays and open new file
      if (first) then
        if (myid==0) then
          write(6,*) "Initialise high frequency output"
        end if
        allocate(freqstore(ifull,nwt,freqvars))
        if (localhist) then
          write(ffile,"(a,'.',i6.6)") trim(surfile), myid
        else
          ffile=surfile
        end if
        if (myid==0.or.localhist) then
          call ccnf_create(ffile,fncid)
          ! Turn off the data filling
          call ccnf_nofill(fncid)
          ! Create dimensions
          if (localhist) then
             call ccnf_def_dim(fncid,'longitude',il,adim(1))
             call ccnf_def_dim(fncid,'latitude',jl,adim(2))
          else
             call ccnf_def_dim(fncid,'longitude',il_g,adim(1))
             call ccnf_def_dim(fncid,'latitude',jl_g,adim(2))
          endif
          call ccnf_def_dim(fncid,'lev',kl,adim(3))
          call ccnf_def_dimu(fncid,'time',adim(4))
          ! Define coords.
          call ccnf_def_var(fncid,'longitude','float',1,adim(1:1),ixp)
          call ccnf_put_att(fncid,ixp,'point_spacing',4,'even')
          call ccnf_put_att(fncid,ixp,'units',12,'degrees_east')
          call ccnf_def_var(fncid,'latitude','float',1,adim(2:2),iyp)
          call ccnf_put_att(fncid,iyp,'point_spacing',4,'even')
          call ccnf_put_att(fncid,iyp,'units',13,'degrees_north')
          call ccnf_def_var(fncid,'lev','float',1,adim(3:3),izp)
          call ccnf_put_att(fncid,izp,'positive',4,'down')
          call ccnf_put_att(fncid,izp,'point_spacing',6,'uneven')
          call ccnf_put_att(fncid,izp,'units',11,'sigma_level')
          call ccnf_def_var(fncid,'time','double',1,adim(4:4),idnt)
          call ccnf_put_att(fncid,idnt,'point_spacing',4,'even')
          icy=kdate/10000
          icm=max(1,min(12,(kdate-icy*10000)/100))
          icd=max(1,min(31,(kdate-icy*10000-icm*100)))
          if(icy<100)icy=icy+1900
          ich=ktime/100
          icmi=(ktime-ich*100)
          ics=0
          write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &                   icd,month(icm),icy,ich,icmi,ics
          call ccnf_put_att(fncid,idnt,'time_origin',20,timorg)
          write(grdtim,'("seconds since ",i4.4,"-",i2.2,"-",i2.2," ",
     &         2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
          call ccnf_put_att(fncid,idnt,'units',33,grdtim)
          if (leap==0) then
            call ccnf_put_att(fncid,idnt,'calendar',6,'noleap')
          end if
          call ccnf_def_var(fncid,'kdate','int',1,adim(4:4),idkdate)
          call ccnf_def_var(fncid,'ktime','int',1,adim(4:4),idktime)
          call ccnf_def_var(fncid,'mtimer','int',1,adim(4:4),idmtimer)
          ! header data
          ahead(1)=ds
          ahead(2)=0.  !difknbd
          ahead(3)=0.  ! was rhkuo for kuo scheme
          ahead(4)=0.  !du
          ahead(5)=rlong0     ! needed by cc2hist
          ahead(6)=rlat0      ! needed by cc2hist
          ahead(7)=schmidt    ! needed by cc2hist
          ahead(8)=0.  !stl2
          ahead(9)=0.  !relaxt
          ahead(10)=0.  !hourbd
          ahead(11)=tss_sh
          ahead(12)=vmodmin
          ahead(13)=av_vmod
          ahead(14)=epsp
          nahead(1)=il_g       ! needed by cc2hist
          nahead(2)=jl_g       ! needed by cc2hist
          nahead(3)=kl         ! needed by cc2hist
          nahead(4)=m
          nahead(5)=0          ! nsd not used now
          nahead(6)=io_in
          nahead(7)=nbd
          nahead(8)=0          ! not needed now  
          nahead(9)=mex
          nahead(10)=mup
          nahead(11)=nem
          nahead(12)=mtimer
          nahead(13)=0.
          nahead(14)=nint(dt)  ! needed by cc2hist
          nahead(15)=0         ! not needed now 
          nahead(16)=nhor
          nahead(17)=nkuo
          nahead(18)=khdif
          nahead(19)=kl        ! needed by cc2hist (was kwt)
          nahead(20)=0  !iaa
          nahead(21)=0  !jaa
          nahead(22)=nvad
          nahead(23)=0       ! not needed now      
          nahead(24)=0  !lbd
          nahead(25)=nrun
          nahead(26)=nrunx
          nahead(27)=khor
          nahead(28)=ksc
          nahead(29)=kountr
          nahead(30)=ndiur
          nahead(31)=0  ! spare
          nahead(32)=nhorps
          nahead(33)=nsoil
          nahead(34)=ms        ! needed by cc2hist
          nahead(35)=ntsur
          nahead(36)=nrad
          nahead(37)=kuocb
          nahead(38)=nvmix
          nahead(39)=ntsea
          nahead(40)=0  
          nahead(41)=nextout
          nahead(42)=ilt
          nahead(43)=ntrac     ! needed by cc2hist
          nahead(44)=nsib
          nahead(45)=nrungcm
          nahead(46)=ncvmix
          nahead(47)=ngwd
          nahead(48)=lgwd
          nahead(49)=mup
          nahead(50)=nritch_t
          nahead(51)=ldr
          nahead(52)=nevapls
          nahead(53)=nevapcc
          nahead(54)=nt_adv
          call ccnf_put_attg(fncid,'real_header',nrhead,ahead)
          call ccnf_put_attg(fncid,'int_header',nihead,nahead)
          if(localhist)then
            call ccnf_put_attg(fncid,'processor_num',myid)
            call ccnf_put_attg(fncid,'nproc',nproc)
#ifdef uniform_decomp
            call ccnf_put_attg(fncid,'decomp','uniform')
#else
            call ccnf_put_attg(fncid,'decomp','face')
#endif
          endif 
          ! define variables
          sdim(1:2)=adim(1:2)
          sdim(3)=adim(4)
          lname='x-component 10m wind'
          call attrib(fncid,sdim,3,'uas',lname,'m/s',-130.,130.,
     &                0,1)
          lname='y-component 10m wind'     
          call attrib(fncid,sdim,3,'vas',lname,'m/s',-130.,130.,
     &                0,1)
          lname='Screen temperature'     
          call attrib(fncid,sdim,3,'tscrn',lname,'K',100.,425.,
     &                0,1)
          lname='Precipitation'
          call attrib(fncid,sdim,3,'rnd',lname,'mm/day',0.,1300.,
     &                0,1)
          lname='Snowfall'
          call attrib(fncid,sdim,3,'sno',lname,'mm/day',0.,1300.,
     &                0,1)

          ! end definition mode
          call ccnf_enddef(fncid)
          if (localhist) then
           ! Set these to global indices (relative to panel 0 in uniform decomp)
           do i=1,ipan
              xpnt(i) = float(i) + ioff(0)
           end do
           iduma(1)=1
           iduma(2)=il
           call ccnf_put_vara(fncid,ixp,iduma(1:1),iduma(2:2),
     &                             xpnt(1:il))
           i=1
           do n=1,npan
             do j=1,jpan
               ypnt(i) = float(j) + joff(n-noff) + (n-noff)*il_g
               i=i+1
             end do
           end do
           iduma(1)=1
           iduma(2)=jl
           call ccnf_put_vara(fncid,iyp,iduma(1:1),iduma(2:2),
     &                             ypnt(1:jl))
          else
           do i=1,il_g
              xpnt(i) = float(i)
           end do
           iduma(1)=1
           iduma(2)=il_g
           call ccnf_put_vara(fncid,ixp,iduma(1:1),iduma(2:2),
     &                             xpnt)
           do j=1,jl_g
              ypnt(j) = float(j)
           end do
           iduma(1)=1
           iduma(2)=jl_g
           call ccnf_put_vara(fncid,iyp,iduma(1:1),iduma(2:2),
     &                             ypnt)
          end if
          iduma(1)=1
          iduma(2)=kl
          call ccnf_put_vara(fncid,izp,iduma(1:1),iduma(2:2),sig)
        end if
        first=.false.
        if (myid==0) then
          write(6,*) "Finished initialising high frequency output"
        end if
      end if
      
      ! store output
      ti=mod(ktau,nwt)
      if (ti==0) ti=nwt
      umag=sqrt(u(1:ifull,1)*u(1:ifull,1)+v(1:ifull,1)*v(1:ifull,1))
      uas=0.
      vas=0.
      where (umag>0.)
        uas=u10*u(1:ifull,1)/umag
        vas=u10*v(1:ifull,1)/umag
      end where
      freqstore(:,ti,1)=uas
      freqstore(:,ti,2)=vas
      freqstore(:,ti,3)=tscrn
      freqstore(:,ti,4)=condx*86400./dt
      freqstore(:,ti,5)=conds*86400./dt

      ! wtite data to file
      if (mod(ktau,nwt)==0) then
        if (myid==0.or.localhist) then
          if (myid==0) then
            write(6,*) "Write high frequency output"
          end if
          fiarch=ktau-nwt+1
          start(1)=fiarch
          ncount(1)=nwt
          do i=1,nwt
            tpnt(i)=real(ktau-nwt+i,8)*real(dt,8)
          end do
          call ccnf_put_vara(fncid,idnt,start,ncount,tpnt)
          do i=1,nwt
            datedat(i)=kdate
          end do
          call ccnf_put_vara(fncid,idkdate,start,ncount,datedat)
          do i=1,nwt
            datedat(i)=ktime
          end do
          call ccnf_put_vara(fncid,idktime,start,ncount,datedat)
          do i=1,nwt
            datedat(i)=mtimer+nint(real(i-nwt)*dt/60.)
          end do
          call ccnf_put_vara(fncid,idmtimer,start,ncount,datedat)
        end if
        call freqwrite(fncid,'uas',fiarch,nwt,localhist,
     &                 freqstore(:,:,1))
        call freqwrite(fncid,'vas',fiarch,nwt,localhist,
     &                 freqstore(:,:,2))
        call freqwrite(fncid,'tscrn',fiarch,nwt,localhist,
     &                 freqstore(:,:,3))
        call freqwrite(fncid,'rnd',fiarch,nwt,localhist,
     &                 freqstore(:,:,4))
        call freqwrite(fncid,'sno',fiarch,nwt,localhist,
     &                 freqstore(:,:,5))
      end if
     
      ! close file at end of run
      if (myid==0.or.localhist) then
        call ccnf_sync(fncid)
        if (ktau==ntau) then
          call ccnf_close(fncid)
        end if
      end if
      
      call end_log(outfile_end)
      
      return
      end subroutine freqfile
