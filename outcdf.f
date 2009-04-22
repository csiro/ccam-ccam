c=======================================================================
      subroutine outcdf(rundate,nmi,itype,ms_out)
!     itype=-1  for restart file
!            1  for outfile
!     N.B. subr outcdfs is down the bottom (for nscrn=1)
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg,mtimer
      include 'filnames.h'  ! list of files, read in once only
      include 'kuocom.h'
      include 'liqwpar.h'  ! ifullw
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmgeom.h' ! rlong0,rlat0,schmidt  
      include 'parmhor.h'  ! mhint, m_bs, nt_adv, ndept
      include 'parmvert.h'
      include 'tracers.h'  ! ngas, nllp, ntrac, tr
      character rundate*10
      integer nmi, itype, ms_out

      integer, parameter :: nihead=54
      integer nahead(nihead)

      integer, parameter :: nrhead=14
      real ahead(nrhead)

      include 'netcdf.inc'
      character cdffile*80

      integer ixp,iyp,idlev,idnt,idms
      common/cdfind/ixp,iyp,idlev,idnt,idms
 
      integer nhor,nhorps,khor,khdif,nhorjlm
      real hdiff,hdifmax
      common/parmhdff/nhor,nhorps,hdiff(kl),khor,khdif,hdifmax,nhorjlm

      integer dim(4),dims(4)
      integer xdim,ydim,zdim,tdim,msdim
      character timorg*20
      character grdtim*33
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/

      integer :: ndt, icy, icm, icd, ich, icmi, ics, idv, ier, imode
      integer, save :: idnc
      integer, save :: iarch=0, idnc0=0
      logical :: local ! Each processor writes its local region

      ndt=dt
      ! The localhist variable controls whether the local file option
      !  is used at all. In any case it's only used for the restfile.
      local = localhist .and. itype == 1 ! Only for outfile
c     print *,'entering outcdf for myid,local ',myid,local

      if(myid==0 .or. local)then !  #########################
!      File setup follows
       if(itype==1)then
c       itype=1 outfile
        iarch=iarch+1
        if(local)then
           write(cdffile,"(a,'.',i2.2)") trim(ofile), myid
        else
           cdffile=ofile
        endif
       else
c       itype=-1 restfile
        iarch=1
        if(local)then
           write(cdffile,"(a,'.',i2.2)") trim(restfile), myid
        else
           cdffile=restfile
        endif
        idnc=0
       endif ! ( itype==1)then

       write(6,'("outcdf itype,idnc,iarch,cdffile=",3i5," ",a80)')
     &                   itype,idnc,iarch,cdffile

       if(iarch==1)then
        print *,'nccre of ',cdffile
        idnc = nccre(cdffile, ncclob, ier)
        print *,'idnc,ier=',idnc,ier
c       Turn off the data filling
        imode = ncsfil(idnc,ncnofill,ier)
        print *,'imode=',imode
c       Create dimensions, lon, lat
        if(local)then
           xdim = ncddef(idnc, 'longitude', il, ier)
           ydim = ncddef(idnc, 'latitude', jl, ier)
        else
           xdim = ncddef(idnc, 'longitude', il_g, ier)
           ydim = ncddef(idnc, 'latitude', jl_g, ier)
        endif
        zdim= ncddef(idnc, 'lev', kl, ier)
        msdim= ncddef(idnc, 'zsoil', ms, ier)
        tdim= ncddef(idnc, 'time',ncunlim,ier)
        print *,"xdim,ydim,zdim,tdim"
        print *, xdim,ydim,zdim,tdim

c       define coords.
        ixp = ncvdef(idnc,'longitude',NCFLOAT,1,xdim,ier)
        call ncaptc(idnc,ixp,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnc,ixp,'units',NCCHAR,12,'degrees_east',ier)
        iyp = ncvdef(idnc,'latitude',NCFLOAT,1,ydim,ier)
        call ncaptc(idnc,iyp,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnc,iyp,'units',NCCHAR,13,'degrees_north',ier)
        print *,'ixp,iyp=',ixp,iyp

        idlev = ncvdef(idnc,'lev',NCFLOAT,1,zdim,ier)
        call ncaptc(idnc,idlev,'positive',NCCHAR,4,'down',ier)
        call ncaptc(idnc,idlev,'point_spacing',NCCHAR,6,'uneven',ier)
        call ncaptc(idnc,idlev,'units',NCCHAR,11,'sigma_level',ier)
        call ncaptc(idnc,idlev,'long_name',NCCHAR,11,'sigma_level',ier)
        print *,'idlev=',idlev

        idms = ncvdef(idnc,'zsoil',NCFLOAT,1,msdim,ier)
        call ncaptc(idnc,idms,'point_spacing',NCCHAR,6,'uneven',ier)
        call ncaptc(idnc,idms,'units',NCCHAR,1,'m',ier)
        print *,'idms=',idms

        print *,'tdim,idnc=',tdim,idnc
ccc     idnt = ncvdef(idnc,'time',NCLONG,1,tdim,ier)
        idnt = ncvdef(idnc,'time',NCFLOAT,1,tdim,ier)  ! Sept 2006
        print *,'idnt=',idnt
        call ncaptc(idnc,idnt,'point_spacing',NCCHAR,4,'even',ier)

        print *,'kdate,ktime,ktau=',kdate,ktime,ktau
        icy=kdate/10000
        icm=max(1,min(12,(kdate-icy*10000)/100))
        icd=max(1,min(31,(kdate-icy*10000-icm*100)))
        if(icy<100)icy=icy+1900
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(6,*) icy,icm,icd,ich,icmi,ics
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &               icd,month(icm),icy,ich,icmi,ics
        print *,'timorg=',timorg
        call ncaptc(idnc,idnt,'time_origin',NCCHAR,20,timorg,ier)

        write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",
     &       2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
        print *,'grdtim=',grdtim
        call ncaptc(idnc,idnt,'units',NCCHAR,33,grdtim,ier)

        dim(1) = xdim
        dim(2) = ydim
        dim(3) = zdim
        dim(4) = tdim

        dims(1) = xdim
        dims(2) = ydim
        dims(3) = msdim
        dims(4) = tdim

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
        nahead(14)=ndt       ! needed by cc2hist
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
        nahead(40)=ms_out     
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
        write(6,'("nahead=",(20i4))') nahead
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
        write(6,*) "ahead=",ahead
        call ncapt(idnc,ncglobal,'int_header',nclong,nihead,nahead,ier)
        if(ier.ne.0)write(6,*)"ncapt int idnc,ier=",idnc,ier
        call ncapt(idnc,ncglobal,'real_header',ncfloat,nrhead,ahead,ier)
        if(ier.ne.0)write(6,*)"ncapt real idnc,ier=",idnc,ier
        call ncaptc(idnc,ncglobal,'date_header',ncchar,10,rundate,ier)
        if(ier.ne.0)write(6,*)"ncaptc date idnc,ier=",idnc,ier

        idv=ncvdef(idnc,'ds',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef ds idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'du',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef du idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'rnml',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef rnml idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'tanl',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef tanl idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'stl1',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl1 idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'stl2',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl2 idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'dt',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef dt idnc,ier=",idnc,ier
       endif ! ( iarch=1)then

      endif ! (myid==0.or.local) #########################
      ! openhist writes some fields so needs to be called by all processes
      call openhist(iarch,itype,dim,local,idnc)

      if(myid==0.or.local)then
        call ncsnc(idnc,ier)
        if(ier.ne.0)write(6,*)"ncsnc idnc,ier=",idnc,ier
      endif    ! (myid==0.or.local)
      if(ktau.eq.ntau.and.(myid==0.or.localhist))then
        call ncclos(idnc,ier)
        write(6,*) "calling ncclos(idnc,ier) ",idnc,ier
      endif
      return   ! outcdf  
      end
c=======================================================================
      subroutine openhist(iarch,itype,dim,local,idnc)
      use ateb ! MJT urban
      use cc_mpi
      use cable_ccam, only : savetile ! MJT cable
      use define_dimensions, only : ncs, ncp ! MJT cable      
      implicit none

c     this routine creates attributes and writes output

      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'
      include 'carbpools.h' ! MJT cable
      include 'dates.h'    ! ktime,kdate,timer,timeg,xg,yg,mtimer
      include 'extraout.h' ! u10_3hr,v10_3hr
      include 'filnames.h' ! list of files, read in once only
      include 'histave.h'
      include 'kuocom.h'
      include 'liqwpar.h'  ! ifullw
      include 'map.h'
      include 'mapproj.h'
      include 'morepbl.h'
      include 'netcdf.inc'
      include 'nsibd.h' ! rsmin,ivegt,sigmf,tgg,tgf,ssdn,res,rmc,isoilm,ico2em
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmvert.h'
      include 'pbl.h'
      include 'prec.h'
      include 'raddiag.h'
      include 'scamdim.h'
      include 'screen.h'
      include 'sigs.h'
      include 'soil.h'
      include 'soilsnow.h'
      include 'soilv.h'   ! sfc,zse
      include 'tracers.h'
      include 'trcom2.h'
      include 'vegpar.h' ! MJT cable
      include 'version.h'
      include 'vvel.h'    ! sdot, dpsldt

      integer iarch, itype
      logical, intent(in) :: local
      character lname*40,mnam*21,nnam*21,expdesc*50
      integer dim(4)
      integer idim(3)
      real xpnt(il_g),ypnt(jl_g)

      integer ixp,iyp,idlev,idnt,idms
      common/cdfind/ixp,iyp,idlev,idnt,idms
      real dpsdt,dpsdtb,dpsdtbb
      common/dpsdt/dpsdt(ifull),dpsdtb(ifull),dpsdtbb(ifull) !globpe,adjust5,outcdf
      real pmsl,aa,bb,cc,dum2
      common/work2/pmsl(ifull),aa(ifull),bb(ifull),cc(ifull),
     &             dum2(ifull,14)
      real tmpry
      common/work3c/tmpry(ifull,kl)

      integer i, idkdate, idktau, idktime, idmtimer, idnteg, idnter,
     &     idv, ier, iq, isoil, j, k, igas, idnc
      real trmax, trmin
      character*3 mon(12),trnum
      real cfrac
      common/cfrac/cfrac(ifull,kl)     ! globpe,radriv90,vertmix,convjlm
      real zsoil(ms)
      real, dimension(ifull,1:12) :: urban ! MJT urban
      data mon/'JAN','FEB','MAR','APR','MAY','JUN'
     &        ,'JUL','AUG','SEP','OCT','NOV','DEC'/

      if(myid == 0 .or. local)then  !#########################
       print *,'openhist itype,iarch,idnc=',itype,iarch,idnc

c      if this is the first archive, set up some global attributes
       if(iarch==1) then
        print *,'dim=',dim
        idim(1)=dim(1)
        idim(2)=dim(2)
        idim(3)=dim(4)
        print *,'idim=',idim

c       Create global attributes
c       Model run number
        print *,'nrun=',nrun
        call ncapt(idnc,ncglobal,'nrun',nclong,1,nrun,ier)
        write(6,*)"nrun ier=",ier

c       Experiment description
        expdesc = 'CCAM model run'
        call ncaptc(idnc,ncglobal,'expdesc',ncchar,len_trim(expdesc),
     &              expdesc,ier)
        write(6,*)"expdesc ier=",ier

c       Model version
        call ncaptc(idnc,ncglobal,'version',ncchar,len_trim(version),
     &              version,ier)

        if(local)then
           ier = nf_put_att_int(idnc,nf_global,"processor_num",nf_int,
     &                          1,myid)
        endif           

c       Sigma levels
        print *,'sig=',sig
        call ncapt(idnc,ncglobal,'sigma',ncfloat,kl,sig,ier)

        lname = 'timer (hrs)'
        idnter = ncvdef(idnc,'timer',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnter,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'mtimer (mins)'
        idmtimer = ncvdef(idnc,'mtimer',nclong,1,dim(4),ier)
        call ncaptc(idnc,idmtimer,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'timeg (UTC)'
        idnteg = ncvdef(idnc,'timeg',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnteg,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'number of time steps from start'
        idktau = ncvdef(idnc,'ktau',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktau,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'year-month-day at start of run'
        idkdate = ncvdef(idnc,'kdate',nclong,1,dim(4),ier)
        call ncaptc(idnc,idkdate,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'hour-minute at start of run'
        idktime = ncvdef(idnc,'ktime',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktime,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        idv = ncvdef(idnc,'sigma', ncfloat, 1, dim(3),ier)
        call ncaptc(idnc,idv,'positive',ncchar
     &             ,len_trim('down'),'down',ier)

        print *,'define attributes of variables'

        lname ='Scaled Log Surface pressure'
        call attrib(idnc,idim,3,'psf',lname,'none',-1.3,0.2,0)

        lname ='Mean sea level pressure'
        call attrib(idnc,idim,3,'pmsl',lname,'hPa',800.,1200.,0)
        lname = 'Surface geopotential'
        call attrib(idnc,idim,2,'zht',lname,'m2/s2',-100.,90.e3,0)
c       call attrib(idnc,idim,2,'zht',lname,'m2/s2',-1.e6,90.e3,0) ! ocean too

c       For time invariant surface fields
        lname = 'Map factor'
        call attrib(idnc,idim,2,'map',lname,'none',.001,1500.,0)
        lname = 'Coriolis factor'
        call attrib(idnc,idim,2,'cor',lname,'1/sec',-1.5e-4,1.5e-4,0)
        lname = 'Rsmin'
        call attrib(idnc,idim,2,'rsmin',lname,'none',0.,200.,0)
        lname = 'Vegetation fraction'
        call attrib(idnc,idim,2,'sigmf',lname,'none',0.,1.,0)
        lname = 'Surface roughness'
        call attrib(idnc,idim,2,'zolnd',lname,'m',0.,40.,0)
        lname = 'Soil type'
        call attrib(idnc,idim,2,'soilt',lname,'none',0.,40.,0)
        lname = 'Vegetation type'
        call attrib(idnc,idim,2,'vegt',lname,'none',0.,44.,0)
        !lname = 'Initial wetness fraction layer 3' ! MJT delete
        !call attrib(idnc,idim,2,'wetfrac',lname,'none',-2.,5.,0)

c       For time varying surface fields
        lname = 'Surface temperature'
        call attrib(idnc,idim,3,'tsu',lname,'K',0.,350.,0)
        lname = 'Pan temperature'
        call attrib(idnc,idim,3,'tpan',lname,'K',0.,350.,0)
        lname = 'Precipitation'
        call attrib(idnc,idim,3,'rnd',lname,'mm/day',0.,1000.,0)
        lname = 'Convective precipitation'
        call attrib(idnc,idim,3,'rnc',lname,'mm/day',0.,1000.,0)
        call attrib(idnc,idim,3,'sno','snowfall','mm/day',0.,1000.,0)
        call attrib(idnc,idim,3,'runoff','Runoff','mm/day',0.,1000.,0)
        lname = '3hr precipitation'
        call attrib(idnc,idim,3,'rnd03',lname,'mm',0.,1000.,1)
        lname = '6hr precipitation'
        call attrib(idnc,idim,3,'rnd06',lname,'mm',0.,1000.,1)
        lname = '9hr precipitation'
        call attrib(idnc,idim,3,'rnd09',lname,'mm',0.,1000.,1)
        lname = '12hr precipitation'
        call attrib(idnc,idim,3,'rnd12',lname,'mm',0.,1000.,1)
        lname = '15hr precipitation'
        call attrib(idnc,idim,3,'rnd15',lname,'mm',0.,1000.,1)
        lname = '18hr precipitation'
        call attrib(idnc,idim,3,'rnd18',lname,'mm',0.,1000.,1)
        lname = '21hr precipitation'
        call attrib(idnc,idim,3,'rnd21',lname,'mm',0.,1000.,1)
        lname = '24hr precipitation'
        call attrib(idnc,idim,3,'rnd24',lname,'mm',0.,1000.,1)
        lname = 'Maximum precip rate in a timestep'
        call attrib(idnc,idim,3,'maxrnd',lname,'mm/day',0.,2000.,1)
        lname = 'Maximum screen temperature'
        call attrib(idnc,idim,3,'tmaxscr',lname,'K',100.,400.,1)
        lname = 'Minimum screen temperature'
        call attrib(idnc,idim,3,'tminscr',lname,'K',100.,400.,1)
        lname = 'Average screen temperature'
        call attrib(idnc,idim,3,'tscr_ave',lname,'K',100.,400.,0)
        lname = 'Screen temperature'
        call attrib(idnc,idim,3,'tscrn',lname,'K',100.,400.,0)
        lname = 'Screen mixing ratio'
        call attrib(idnc,idim,3,'qgscrn',lname,'kg/kg',0.,.06,0)
        lname = 'Maximum screen relative humidity'
        call attrib(idnc,idim,3,'rhmaxscr',lname,'%',0.,200.,1)
        lname = 'Minimum screen relative humidity'
        call attrib(idnc,idim,3,'rhminscr',lname,'%',0.,200.,1)
        lname = 'Maximum daily Cape'
        call attrib(idnc,idim,3,'capemax',lname,'J/kg',0.,20000.,0) 
        lname = 'x-component max 10m wind'
        call attrib(idnc,idim,3,'u10max',lname,'m/s',-99.,99.,1)
        lname = 'y-component max 10m wind'
        call attrib(idnc,idim,3,'v10max',lname,'m/s',-99.,99.,1)
        lname = 'x-component max level_1 wind'
        call attrib(idnc,idim,3,'u1max',lname,'m/s',-99.,99.,1)
        lname = 'y-component max level_1 wind'
        call attrib(idnc,idim,3,'v1max',lname,'m/s',-99.,99.,1)
        lname = 'x-component max level_2 wind'
        call attrib(idnc,idim,3,'u2max',lname,'m/s',-99.,99.,1)
        lname = 'y-component max level_2 wind'
        call attrib(idnc,idim,3,'v2max',lname,'m/s',-99.,99.,1)
        lname = '10m wind speed'
        call attrib(idnc,idim,3,'u10',lname,'m/s',0.,60.,0)
c       lname = '3m wind speed'
c       call attrib(idnc,idim,3,'u3',lname,'K',0.,60.,0)
        lname = 'Screen level wind speed'
        call attrib(idnc,idim,3,'uscrn',lname,'K',0.,40.,0)
        lname = 'Surface albedo'
        call attrib(idnc,idim,3,'alb',lname,'none',0.,1.,0)
        lname = 'Sea ice depth'
        call attrib(idnc,idim,3,'siced',lname,'m',0.,50.,0)
        lname = 'Sea ice fraction'
        call attrib(idnc,idim,3,'fracice',lname,'none',0.,1.,0)
        lname = 'Snow depth (liquid water)'
c       call attrib(idnc,idim,3,'snd',lname,'mm',0.,5000.,0)
        call attrl (idnc,idim,3,'snd',lname,'mm',0.,5000.,0)  ! long
        lname = 'Low cloud ave'
        call attrib(idnc,idim,3,'cll',lname,'frac',0.,1.,0)
        lname = 'Mid cloud ave'
        call attrib(idnc,idim,3,'clm',lname,'frac',0.,1.,0)
        lname = 'Hi cloud ave'
        call attrib(idnc,idim,3,'clh',lname,'frac',0.,1.,0)
        lname = 'Total cloud ave'
        call attrib(idnc,idim,3,'cld',lname,'frac',0.,1.,0)
        lname = 'x-component wind stress'
        call attrib(idnc,idim,3,'taux',lname,'N/m2',-50.,50.,0)
        lname = 'y-component wind stress'
        call attrib(idnc,idim,3,'tauy',lname,'N/m2',-50.,50.,0)
        !lname = 'Soil moisture as frac FC levels 1-2' ! MJT delete
        !call attrib(idnc,idim,3,'wbfshal',lname,'frac',0.,4.,0)
        !lname = 'Soil moisture as frac FC levels 3-4'
        !call attrib(idnc,idim,3,'wbfroot',lname,'frac',0.,4.,0)
        !lname = 'Soil moisture as frac FC levels 1-6'
        !call attrib(idnc,idim,3,'wbftot',lname,'frac',0.,4.,0)
        if ((nsib.eq.4).or.(nsib.eq.6)) then  ! MJT cable
          lname = 'Carbon leaf pool'
          call attrib(idnc,idim,3,'cplant1',lname,'none',0.,50000.,0)
          lname = 'Carbon wood pool'
          call attrib(idnc,idim,3,'cplant2',lname,'none',0.,50000.,0)
          lname = 'Carbon root pool'
          call attrib(idnc,idim,3,'cplant3',lname,'none',0.,50000.,0)
          lname = 'Carbon soil fast pool'
          call attrib(idnc,idim,3,'csoil1',lname,'none',0.,50000.,0)
          lname = 'Carbon soil slow pool'
          call attrib(idnc,idim,3,'csoil2',lname,'none',0.,50000.,0)
          lname = 'cansto'
          call attrib(idnc,idim,3,'cansto',lname,'none',0.,10.,0)
        endif
        if(nextout>=1) then
          print *,'nextout=',nextout
          lname = 'LW at TOA'
          call attrib(idnc,idim,3,'rtu_ave',lname,'W/m2',0.,800.,0)
          lname = 'Clear sky LW at TOA'
          call attrib(idnc,idim,3,'rtc_ave',lname,'W/m2',0.,800.,0)
          lname = 'LW downwelling at ground'
          call attrib(idnc,idim,3,'rgdn_ave',lname,'W/m2',-500.,1.e3,0)
          lname = 'LW net at ground (+ve up)'
          call attrib(idnc,idim,3,'rgn_ave',lname,'W/m2',-500.,1000.,0)
          lname = 'Clear sky LW at ground'
          call attrib(idnc,idim,3,'rgc_ave',lname,'W/m2',-500.,1000.,0)
          lname = 'Solar in at TOA'
          call attrib(idnc,idim,3,'sint_ave',lname,'W/m2',0.,1600.,0)
          lname = 'Solar out at TOA'
          call attrib(idnc,idim,3,'sot_ave',lname,'W/m2',0.,1000.,0)
          lname = 'Clear sky SW out at TOA'
          call attrib(idnc,idim,3,'soc_ave',lname,'W/m2',0.,900.,0)
          lname = 'Solar downwelling at ground'
          call attrib(idnc,idim,3,'sgdn_ave',lname,'W/m2',-500.,2.e3,0)
          lname = 'Solar net at ground (+ve down)'
          call attrib(idnc,idim,3,'sgn_ave',lname,'W/m2',-500.,2000.,0)
          lname = 'Surface pressure tendency'
          call attrib(idnc,idim,3,'dpsdt',lname,'hPa/day',-400.,400.,0)
          lname = 'PBL depth'
          call attrib(idnc,idim,3,'pblh',lname,'m',0.,6000.,0)
          lname = 'friction velocity'
          call attrib(idnc,idim,3,'ustar',lname,'m/s',0.,10.,0)
        endif     ! (nextout>=1)
        if(nextout>=2) then  ! 6-hourly u10, v10, tscr, rh1
         mnam ='x-component 10m wind '
         nnam ='y-component 10m wind '
         call attrib(idnc,idim,3,'u10_06',mnam//'6hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_06',nnam//'6hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'u10_12',mnam//'12hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_12',nnam//'12hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'u10_18',mnam//'18hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_18',nnam//'18hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'u10_24',mnam//'24hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_24',nnam//'24hr','m/s',-99.,99.,1)
         mnam ='tscrn 3-hrly'
         nnam ='rhum level_1 3-hrly'
         call attrib(idnc,idim,3,'tscr_06',mnam//'6hr', 'K',100.,400.,1)
         call attrib(idnc,idim,3,'tscr_12',mnam//'12hr','K',100.,400.,1)
         call attrib(idnc,idim,3,'tscr_18',mnam//'18hr','K',100.,400.,1)
         call attrib(idnc,idim,3,'tscr_24',mnam//'24hr','K',100.,400.,1)
         call attrib(idnc,idim,3,'rh1_06', nnam//'6hr', '%',-9.,200.,1)
         call attrib(idnc,idim,3,'rh1_12', nnam//'12hr','%',-9.,200.,1)
         call attrib(idnc,idim,3,'rh1_18', nnam//'18hr','%',-9.,200.,1)
         call attrib(idnc,idim,3,'rh1_24', nnam//'24hr','%',-9.,200.,1)
        endif     ! (nextout>=2)
        if(nextout>=3) then  ! also 3-hourly u10, v10, tscr, rh1
         call attrib(idnc,idim,3,'tscr_03',mnam//'3hr', 'K',100.,400.,1)
         call attrib(idnc,idim,3,'tscr_09',mnam//'9hr', 'K',100.,400.,1)
         call attrib(idnc,idim,3,'tscr_15',mnam//'15hr','K',100.,400.,1)
         call attrib(idnc,idim,3,'tscr_21',mnam//'21hr','K',100.,400.,1)
         call attrib(idnc,idim,3,'rh1_03', nnam//'3hr', '%',-9.,200.,1)
         call attrib(idnc,idim,3,'rh1_09', nnam//'9hr', '%',-9.,200.,1)
         call attrib(idnc,idim,3,'rh1_15', nnam//'15hr','%',-9.,200.,1)
         call attrib(idnc,idim,3,'rh1_21', nnam//'21hr','%',-9.,200.,1)
         mnam ='x-component 10m wind '
         nnam ='y-component 10m wind '
         call attrib(idnc,idim,3,'u10_03',mnam//'3hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_03',nnam//'3hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'u10_09',mnam//'9hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_09',nnam//'9hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'u10_15',mnam//'15hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_15',nnam//'15hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'u10_21',mnam//'21hr','m/s',-99.,99.,1)
         call attrib(idnc,idim,3,'v10_21',nnam//'21hr','m/s',-99.,99.,1)
        endif     ! (nextout>=3)

        lname = 'Soil temperature lev 1'
        call attrib(idnc,idim,3,'tgg1',lname,'K',100.,400.,0)
        lname = 'Soil temperature lev 2'
        call attrib(idnc,idim,3,'tgg2',lname,'K',100.,400.,0)
        lname = 'Soil temperature lev 3'
        call attrib(idnc,idim,3,'tgg3',lname,'K',100.,400.,0)
        lname = 'Soil temperature lev 4'
        call attrib(idnc,idim,3,'tgg4',lname,'K',100.,400.,0)
        lname = 'Soil temperature lev 5'
        call attrib(idnc,idim,3,'tgg5',lname,'K',100.,400.,0)
        lname = 'Soil temperature lev 6'
        call attrib(idnc,idim,3,'tgg6',lname,'K',100.,400.,0)
       ! lname = 'Soil moisture lev 1' ! MJT delete
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
        lname = 'Net radiation'
        call attrib(idnc,idim,3,'rnet',lname,'W/m2',-3000.,3000.,0)
        lname = 'Avg cloud base'
        call attrib(idnc,idim,3,'cbas_ave',lname,'sigma',0.,1.1,0)
        lname = 'Avg cloud top'
        call attrib(idnc,idim,3,'ctop_ave',lname,'sigma',0.,1.1,0)
        lname = 'Avg dew flux'
        call attrib(idnc,idim,3,'dew_ave',lname,'W/m2',-100.,1000.,0)
        lname = 'Avg evaporation'
        call attrib(idnc,idim,3,'evap',lname,'mm',-100.,100.,0)
        lname = 'Avg potential "pan" evaporation'
        call attrib(idnc,idim,3,'epan_ave',lname,'W/m2',-1000.,10.e3,0)
        lname = 'Potential "pan" evaporation'
        call attrib(idnc,idim,3,'epan',lname,'W/m2',-1000.,10.e3,0)
        lname = 'Avg potential evaporation'
        call attrib(idnc,idim,3,'epot_ave',lname,'W/m2',-1000.,10.e3,0)
        lname = 'Latent heat flux'
        call attrib(idnc,idim,3,'eg',lname,'W/m2',-1000.,3000.,0)
        lname = 'Avg latent heat flux'
        call attrib(idnc,idim,3,'eg_ave',lname,'W/m2',-1000.,3000.,0)
        lname = 'Sensible heat flux'
        call attrib(idnc,idim,3,'fg',lname,'W/m2',-3000.,3000.,0)
        lname = 'Avg sensible heat flux'
        call attrib(idnc,idim,3,'fg_ave',lname,'W/m2',-3000.,3000.,0)
        lname = 'Avg flux into tgg1 layer'
        call attrib(idnc,idim,3,'ga_ave',lname,'W/m2',-1000.,1000.,0)
        lname = 'Avg ice water path'
        call attrib(idnc,idim,3,'iwp_ave',lname,'kg/m2',0.,2.,0)
        lname = 'Avg liquid water path'
        call attrib(idnc,idim,3,'lwp_ave',lname,'kg/m2',0.,2.,0)

!       rml 16/02/06 set attributes for trNNN and travNNN
        if (ngas>0) then 
         do igas=1,ngas
           write(trnum,'(i3.3)') igas
           trmax=max(1.,10.*maxval(tr(:,:,igas))) !max to avoid trmax and trmin=0
           trmin=gasmin(igas) !gasmin needed in adjust5, set in tracers.h
           lname = 'Tracer (inst.) '//trnum
           call attrib(idnc,dim,4,'tr'//trnum,lname,'ppm',trmin,trmax,0)
           lname = 'Tracer (average) '//trnum
           call attrib(idnc,dim,4,'trav'//trnum,lname,'ppm',trmin,trmax
     &                 ,0)
         enddo ! igas loop
        endif  ! (ntrac.gt.0)

        print *,'3d variables'
        if(nextout>=4.and.nllp==3)then   ! N.B. use nscrn=1 for hourly output
          lname = 'Delta latitude'
          call attrib(idnc,dim,4,'del_lat',lname,'deg',-60.,60.,1)
          lname = 'Delta longitude'
          call attrib(idnc,dim,4,'del_lon',lname,'deg',-180.,180.,1)
          lname = 'Delta pressure'
          call attrib(idnc,dim,4,'del_p',lname,'hPa',-900.,900.,1)
        endif  ! (nextout>=4.and.nllp==3)
        call attrib(idnc,dim,4,'temp','Air temperature','K',100.,350.,0)
        lname= 'x-component wind'
        call attrib(idnc,dim,4,'u',lname,'m/s',-150.,150.,0)
        lname= 'y-component wind'
        call attrib(idnc,dim,4,'v',lname,'m/s',-150.,150.,0)
        lname= 'vertical velocity'
        call attrib(idnc,dim,4,'omega',lname,'Pa/s',-50.,50.,0)
        lname= 'Water mixing ratio'
        call attrib(idnc,dim,4,'mixr',lname,'kg/kg',0.,.05,0)
        if(ldr.ne.0)then
         call attrib(idnc,dim,4,'qfg','Frozen water','kg/kg',0.,.02,0)
         call attrib(idnc,dim,4,'qlg','Liquid water','kg/kg',0.,.02,0)
         call attrib(idnc,dim,4,'cfrac','Cloud fraction','none',0.,1.,0)
        endif

        if(itype==-1)then   ! extra stuff just written for restart file
         lname= 'sdot: change in grid spacing per time step +.5'
         call attrib(idnc,dim,4,'sdot',lname,'1/ts',-3.,3.,0) 
         lname = 'Soil ice lev 1'
         call attrib(idnc,idim,3,'wbice1',lname,'m3/m3',0.,1.,0)
         lname = 'Soil ice lev 2'
         call attrib(idnc,idim,3,'wbice2',lname,'m3/m3',0.,1.,0)
         lname = 'Soil ice lev 3'
         call attrib(idnc,idim,3,'wbice3',lname,'m3/m3',0.,1.,0)
         lname = 'Soil ice lev 4'
         call attrib(idnc,idim,3,'wbice4',lname,'m3/m3',0.,1.,0)
         lname = 'Soil ice lev 5'
         call attrib(idnc,idim,3,'wbice5',lname,'m3/m3',0.,1.,0)
         lname = 'Soil ice lev 6'
         call attrib(idnc,idim,3,'wbice6',lname,'m3/m3',0.,1.,0)
         lname = 'Snow temperature lev 1'
         call attrib(idnc,idim,3,'tggsn1',lname,'K',100.,400.,0)
         lname = 'Snow temperature lev 2'
         call attrib(idnc,idim,3,'tggsn2',lname,'K',100.,400.,0)
         lname = 'Snow temperature lev 3'
         call attrib(idnc,idim,3,'tggsn3',lname,'K',100.,400.,0)
         lname = 'Snow mass lev 1'
         call attrib(idnc,idim,3,'smass1',lname,'K',0.,400.,0)
         lname = 'Snow mass lev 2'
         call attrib(idnc,idim,3,'smass2',lname,'K',0.,400.,0)
         lname = 'Snow mass lev 3'
         call attrib(idnc,idim,3,'smass3',lname,'K',0.,400.,0)
         lname = 'Snow density lev 1'
         call attrib(idnc,idim,3,'ssdn1',lname,'K',0.,400.,0)
         lname = 'Snow density lev 2'
         call attrib(idnc,idim,3,'ssdn2',lname,'K',0.,400.,0)
         lname = 'Snow density lev 3'
         call attrib(idnc,idim,3,'ssdn3',lname,'K',0.,400.,0)
         lname = 'Snow age'
         call attrib(idnc,idim,3,'snage',lname,'none',0.,20.,0)   
         lname = 'Snow flag'
         call attrib(idnc,idim,3,'sflag',lname,'none',0.,4.,0)
         lname = 'Soil turbulent resistance' ! MJT cable
         call attrib(idnc,idim,3,'rtsoil',lname,'none',0.,9.e4,0) 
        endif  ! (itype==-1)

        !--------------------------------------------------------
        ! MJT urban
        if ((nurban.eq.-1).or.((nurban.eq.1).and.(itype==-1))) then
         lname = 'roof temperature lev 1'
         call attrib(idnc,idim,3,'rooftgg1',lname,'K',100.,400.,0)
         lname = 'roof temperature lev 2'
         call attrib(idnc,idim,3,'rooftgg2',lname,'K',100.,400.,0)
         lname = 'roof temperature lev 3'
         call attrib(idnc,idim,3,'rooftgg3',lname,'K',100.,400.,0)
         lname = 'east wall temperature lev 1'
         call attrib(idnc,idim,3,'waletgg1',lname,'K',100.,400.,0)
         lname = 'east wall temperature lev 2'
         call attrib(idnc,idim,3,'waletgg2',lname,'K',100.,400.,0)
         lname = 'east wall temperature lev 3'
         call attrib(idnc,idim,3,'waletgg3',lname,'K',100.,400.,0)
         lname = 'west wall temperature lev 1'
         call attrib(idnc,idim,3,'walwtgg1',lname,'K',100.,400.,0)
         lname = 'west wall temperature lev 2'
         call attrib(idnc,idim,3,'walwtgg2',lname,'K',100.,400.,0)
         lname = 'west wall temperature lev 3'
         call attrib(idnc,idim,3,'walwtgg3',lname,'K',100.,400.,0)
         lname = 'road temperature lev 1'
         call attrib(idnc,idim,3,'roadtgg1',lname,'K',100.,400.,0)
         lname = 'road temperature lev 2'
         call attrib(idnc,idim,3,'roadtgg2',lname,'K',100.,400.,0)
         lname = 'road temperature lev 3'
         call attrib(idnc,idim,3,'roadtgg3',lname,'K',100.,400.,0)
        end if
        !--------------------------------------------------------  

        !-------------------------------------------------------
        ! MJT CHANGE - add wetfrac1-6 and possibly delete wb1-6 above
        lname = 'Wetness fraction layer 1' ! 5. for frozen sand
        call attrib(idnc,idim,3,'wetfrac1',lname,'none',-2.,5.,0)
        lname = 'Wetness fraction layer 2'
        call attrib(idnc,idim,3,'wetfrac2',lname,'none',-2.,5.,0)
        lname = 'Wetness fraction layer 3'
        call attrib(idnc,idim,3,'wetfrac3',lname,'none',-2.,5.,0)
        lname = 'Wetness fraction layer 4'
        call attrib(idnc,idim,3,'wetfrac4',lname,'none',-2.,5.,0)
        lname = 'Wetness fraction layer 5'
        call attrib(idnc,idim,3,'wetfrac5',lname,'none',-2.,5.,0)
        lname = 'Wetness fraction layer 6'
        call attrib(idnc,idim,3,'wetfrac6',lname,'none',-2.,5.,0)
        !-------------------------------------------------------        
 
        print *,'finished defining attributes'
c       Leave define mode
        call ncendf(idnc,ier)
        print *,'leave define mode: ier=',ier

        if(local)then
           ! Set these to global indices
           do i=1,ipan
              xpnt(i) = float(i) + ioff
           end do
           call ncvpt(idnc,ixp,1,il,xpnt,ier)
           do j=1,jl
              ypnt(j) = float(j) + jl*((myid*il)/il_g)
           end do
           call ncvpt(idnc,iyp,1,jl,ypnt,ier)
        else
           do i=1,il_g
              xpnt(i) = float(i)
           end do
           call ncvpt(idnc,ixp,1,il_g,xpnt,ier)
           do j=1,jl_g
              ypnt(j) = float(j)
           end do
           call ncvpt(idnc,iyp,1,jl_g,ypnt,ier)
        endif

        call ncvpt(idnc,idlev,1,kl,sig,ier)

        idv = ncvid(idnc,'sigma',ier)
        call ncvpt(idnc,idv,1,kl,sig,ier)

        idv = ncvid(idnc,'lev',ier)
        call ncvpt(idnc,idv,1,kl,sig,ier)

        zsoil(1)=.5*zse(1)
        zsoil(2)=zse(1)+zse(2)*.5
        zsoil(3)=zse(1)+zse(2)+zse(3)*.5
        zsoil(4)=zse(1)+zse(2)+zse(3)+zse(4)*.5
        zsoil(5)=zse(1)+zse(2)+zse(3)+zse(4)+zse(5)*.5
        zsoil(6)=zse(1)+zse(2)+zse(3)+zse(4)+zse(5)+zse(6)*.5
        call ncvpt(idnc,idms,1,ms,zsoil,ier)

        idv = ncvid(idnc,'ds',ier)
        call ncvpt1(idnc,idv,1,ds,ier)
        idv = ncvid(idnc,'tanl',ier)
        call ncvpt1(idnc,idv,1,tanl,ier)
        idv = ncvid(idnc,'rnml',ier)
        call ncvpt1(idnc,idv,1,rnml,ier)
        idv = ncvid(idnc,'du',ier)
        call ncvpt1(idnc,idv,1,du,ier)
        idv = ncvid(idnc,'stl1',ier)
        call ncvpt1(idnc,idv,1,stl1,ier)
        idv = ncvid(idnc,'stl2',ier)
        call ncvpt1(idnc,idv,1,stl2,ier)
        idv = ncvid(idnc,'dt',ier)
        call ncvpt1(idnc,idv,1,dt,ier)
       endif ! iarch==1
!      -----------------------------------------------------------      
       print *,'outcdf processing kdate,ktime,ktau,mtimer: ',
     .                            kdate,ktime,ktau,mtimer
c      set time to number of minutes since start 
       idv = ncvid(idnc,'time',ier)
ccc    call ncvpt1(idnc,idv,iarch,mtimer,ier)
       call ncvpt1(idnc,idv,iarch,real(mtimer),ier)  ! Sept 2006

       idv = ncvid(idnc,'timer',ier)
       call ncvpt1(idnc,idv,iarch,timer,ier)
       idv = ncvid(idnc,'mtimer',ier)
       call ncvpt1(idnc,idv,iarch,mtimer,ier)
       idv = ncvid(idnc,'timeg',ier)
       call ncvpt1(idnc,idv,iarch,timeg,ier)
       idv = ncvid(idnc,'ktau',ier)
       call ncvpt1(idnc,idv,iarch,ktau,ier)
       idv = ncvid(idnc,'kdate',ier)
       call ncvpt1(idnc,idv,iarch,kdate,ier)
       idv = ncvid(idnc,'ktime',ier)
       call ncvpt1(idnc,idv,iarch,ktime,ier)
       print *,'kdate,ktime,ktau=',kdate,ktime,ktau
       print *,'timer,timeg=',timer,timeg
       print *,'now write out variables'
      endif ! myid == 0 .or. local

      if(ktau==0.or.itype==-1)then  ! also for restart file
!       write time-invariant fields      
        call histwrt3(em,'map',idnc,iarch,local)
        call histwrt3(f,'cor',idnc,iarch,local)
        call histwrt3(rsmin,'rsmin',idnc,iarch,local)
        call histwrt3(sigmf,'sigmf',idnc,iarch,local)
        !--------------------------------------------
        ! MJT urban
        aa=zolnd
        bb=zolnd/7.4 ! dummy
        call atebzo(ifull,aa(:),bb(:),0)
        call histwrt3(aa,'zolnd',idnc,iarch,local)
        !call histwrt3(zolnd,'zolnd',idnc,iarch,local)
        !-------------------------------------------- 
        do iq=1,ifull
         aa(iq)=isoilm(iq)
        enddo
        call histwrt3(aa,'soilt',idnc,iarch,local)
        do iq=1,ifull
!        N.B. subtract 31 to get sib values
         aa(iq)=ivegt(iq)
        enddo
        call histwrt3(aa,'vegt',idnc,iarch,local)
        !do iq=1,ifull ! MJT delete
        ! isoil=isoilm(iq)
        ! aa(iq)=(wb(iq,3)-swilt(isoil))/(sfc(isoil)-swilt(isoil))
        !enddo
        !call histwrt3(aa,'wetfrac',idnc,iarch,local)
      endif ! (ktau==0.or.itype==-1) 

      call histwrt3(zs,'zht',idnc,iarch,local)   ! always from 13/9/02
      call histwrt3(psl,'psf',idnc,iarch,local)
      do iq=1,ifull
        aa(iq)=pmsl(iq)/100.
      enddo
      call histwrt3(aa,'pmsl',idnc,iarch,local)
      call histwrt3(tss,'tsu',idnc,iarch,local)
      aa(:)=0.5*sum(albvisnir(:,:),2) ! MJT CHANGE albedo
      call atebalb1(1,ifull,aa(:),0) ! MJT urban
      call histwrt3(aa,'alb',idnc,iarch,local)
      call histwrt3(tgg(1,1),'tgg1',idnc,iarch,local)
      call histwrt3(tgg(1,2),'tgg2',idnc,iarch,local)
      call histwrt3(tgg(1,3),'tgg3',idnc,iarch,local)
      call histwrt3(tgg(1,4),'tgg4',idnc,iarch,local)
      call histwrt3(tgg(1,5),'tgg5',idnc,iarch,local)
      call histwrt3(tgg(1,6),'tgg6',idnc,iarch,local)
      !call histwrt3(wb(1,1),'wb1',idnc,iarch,local) ! MJT delete
      !call histwrt3(wb(1,2),'wb2',idnc,iarch,local)
      !call histwrt3(wb(1,3),'wb3',idnc,iarch,local)
      !call histwrt3(wb(1,4),'wb4',idnc,iarch,local)
      !call histwrt3(wb(1,5),'wb5',idnc,iarch,local)
      !call histwrt3(wb(1,6),'wb6',idnc,iarch,local)
      if ((nsib.eq.4).or.(nsib.eq.6)) then ! MJT cable
! rml: moved from section that isn't written to restart file
         !call histwrt3(sumpn,'sumpn',idnc,iarch,local)
         !call histwrt3(sumrp,'sumrp',idnc,iarch,local)
         !call histwrt3(sumrs,'sumrs',idnc,iarch,local)
         !call histwrt3(sumrd,'sumrd',idnc,iarch,local)
      call histwrt3(cplant(:,1),'cplant1',idnc,iarch,local)
      call histwrt3(cplant(:,2),'cplant2',idnc,iarch,local)
      call histwrt3(cplant(:,3),'cplant3',idnc,iarch,local)
      call histwrt3(csoil(:,1),'csoil1',idnc,iarch,local)
      call histwrt3(csoil(:,2),'csoil2',idnc,iarch,local)
      call histwrt3(cansto,'cansto',idnc,iarch,local)
      endif    
    !  do iq=1,ifull ! MJT delete
!   !   calculate wb/field_capacity;  up to 3.0 for sand (isoil=1)	   
    !   isoil=isoilm(iq)
    !   aa(iq)=(zse(1)*wb(iq,1)+zse(2)*wb(iq,2))/
    ! .	       ((zse(1)+zse(2))*sfc(isoil))
    !   bb(iq)=(zse(3)*wb(iq,3)+zse(4)*wb(iq,4))/
    ! .	       ((zse(3)+zse(4))*sfc(isoil))
    !   cc(iq)=(zse(1)*wb(iq,1)+zse(2)*wb(iq,2)+zse(3)*wb(iq,3)+
    ! .         zse(4)*wb(iq,4)+zse(5)*wb(iq,5)+zse(6)*wb(iq,6))/
    ! .	       ((zse(1)+zse(2)+zse(3)+zse(4)+zse(5)+zse(6))*sfc(isoil))
    !  enddo
    !  call histwrt3(aa,'wbfshal',idnc,iarch,local)
    !  call histwrt3(bb,'wbfroot',idnc,iarch,local)
    !  call histwrt3(cc,'wbftot',idnc,iarch,local)
      call histwrt3(sicedep,'siced',idnc,iarch,local)
      call histwrt3(fracice,'fracice',idnc,iarch,local)
c     call histwrt3(snowd,'snd',idnc,iarch,local)
      call histwrt3l(snowd,'snd',idnc,iarch,local)  ! long write
      
      if(ktau>0.and.nwt.ne.nperday.and.itype.ne.-1)then  ! reinstated July '05
!       scale up precip,precc,sno,runoff to mm/day (soon reset to 0 in globpe)
!       but, don't scale up for restart file as just done in previous write
!       ktau in next line in case ntau (& thus ktau) < nwt 
        precip=precip*real(nperday)/min(nwt,max(1,ktau))     
        precc =precc *real(nperday)/min(nwt,max(1,ktau))     
        sno   =sno   *real(nperday)/min(nwt,max(1,ktau))     
        runoff=runoff*real(nperday)/min(nwt,max(1,ktau))    
      endif   ! (ktau>0.and.nwt.ne.nperday.and.itype.ne.-1)
      call histwrt3(precip,'rnd',idnc,iarch,local)
      call histwrt3(precc,'rnc',idnc,iarch,local)
      call histwrt3(sno,'sno',idnc,iarch,local)
      call histwrt3(runoff,'runoff',idnc,iarch,local)
      call histwrt3(tpan,'tpan',idnc,iarch,local)
      
      if(ktau>0.and.itype.ne.-1)then  ! these not written to restart file
       if(mod(ktau,nperday)==0.or.ktau==ntau)then  ! only write once per day
         rndmax(:)=rndmax(:)*86400./dt ! scale up to mm/day
         call histwrt3(rndmax,'maxrnd',idnc,iarch,local)
         call histwrt3(tmaxscr,'tmaxscr',idnc,iarch,local)
         call histwrt3(tminscr,'tminscr',idnc,iarch,local)
         call histwrt3(rhmaxscr,'rhmaxscr',idnc,iarch,local)
         call histwrt3(rhminscr,'rhminscr',idnc,iarch,local)
         call histwrt3(capemax,'capemax',idnc,iarch,local)
         call histwrt3(u10max,'u10max',idnc,iarch,local)
         call histwrt3(v10max,'v10max',idnc,iarch,local)
         call histwrt3(u1max,'u1max',idnc,iarch,local)
         call histwrt3(v1max,'v1max',idnc,iarch,local)
         call histwrt3(u2max,'u2max',idnc,iarch,local)
         call histwrt3(v2max,'v2max',idnc,iarch,local)
!        if writes done more than once per day, 
!        needed to augment accumulated 3-hourly rainfall in rnd06 to rnd21 
!        to allow for intermediate zeroing of precip()
!        but not needed from 17/9/03 with introduction of rnd24
         call histwrt3(rnd_3hr(1,1),'rnd03',idnc,iarch,local)
         call histwrt3(rnd_3hr(1,2),'rnd06',idnc,iarch,local)
         call histwrt3(rnd_3hr(1,3),'rnd09',idnc,iarch,local)
         call histwrt3(rnd_3hr(1,4),'rnd12',idnc,iarch,local)
         call histwrt3(rnd_3hr(1,5),'rnd15',idnc,iarch,local)
         call histwrt3(rnd_3hr(1,6),'rnd18',idnc,iarch,local)
         call histwrt3(rnd_3hr(1,7),'rnd21',idnc,iarch,local)
         call histwrt3(rnd_3hr(1,8),'rnd24',idnc,iarch,local)
         if(nextout>=2) then ! 6-hourly u10 & v10
           call histwrt3( u10_3hr(1,2), 'u10_06',idnc,iarch,local)
           call histwrt3( v10_3hr(1,2), 'v10_06',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,2),'tscr_06',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,2), 'rh1_06',idnc,iarch,local)
           call histwrt3( u10_3hr(1,4), 'u10_12',idnc,iarch,local)
           call histwrt3( v10_3hr(1,4), 'v10_12',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,4),'tscr_12',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,4), 'rh1_12',idnc,iarch,local)
           call histwrt3( u10_3hr(1,6), 'u10_18',idnc,iarch,local)
           call histwrt3( v10_3hr(1,6), 'v10_18',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,6),'tscr_18',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,6), 'rh1_18',idnc,iarch,local)
           call histwrt3( u10_3hr(1,8), 'u10_24',idnc,iarch,local)
           call histwrt3( v10_3hr(1,8), 'v10_24',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,8),'tscr_24',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,8), 'rh1_24',idnc,iarch,local)
         endif  ! (nextout>=2)
         if(nextout>=3) then  ! also 3-hourly u10 & v10
           call histwrt3( u10_3hr(1,1), 'u10_03',idnc,iarch,local)
           call histwrt3( v10_3hr(1,1), 'v10_03',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,1),'tscr_03',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,1), 'rh1_03',idnc,iarch,local)
           call histwrt3( u10_3hr(1,3), 'u10_09',idnc,iarch,local)
           call histwrt3( v10_3hr(1,3), 'v10_09',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,3),'tscr_09',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,3), 'rh1_09',idnc,iarch,local)
           call histwrt3( u10_3hr(1,5), 'u10_15',idnc,iarch,local)
           call histwrt3( v10_3hr(1,5), 'v10_15',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,5),'tscr_15',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,5), 'rh1_15',idnc,iarch,local)
           call histwrt3( u10_3hr(1,7), 'u10_21',idnc,iarch,local)
           call histwrt3( v10_3hr(1,7), 'v10_21',idnc,iarch,local)
           call histwrt3(tscr_3hr(1,7),'tscr_21',idnc,iarch,local)
           call histwrt3( rh1_3hr(1,7), 'rh1_21',idnc,iarch,local)
         endif  ! nextout>=3
         if(nextout>=4.and.nllp==3) then  
c         print *,'before corrn ',(tr(idjd,nlv,ngas+k),k=1,3)
          do k=1,klt
           do iq=1,ilt*jlt        
            tr(iq,k,ngas+1)=tr(iq,k,ngas+1)-alat(iq)
            tr(iq,k,ngas+2)=tr(iq,k,ngas+2)-along(iq)
            if(tr(iq,k,ngas+2)>180.)
     &                         tr(iq,k,ngas+2)=tr(iq,k,ngas+2)-360.
            if(tr(iq,k,ngas+2)<-180.)
     &                         tr(iq,k,ngas+2)=tr(iq,k,ngas+2)+360.
            tr(iq,k,ngas+3)=tr(iq,k,ngas+3)-.01*ps(iq)*sig(k)  ! in hPa
           enddo
          enddo
c	   print *,'in outcdf ps, sig ',ps(idjd),sig(nlv)
c	   print *,'after corrn ',(tr(idjd,nlv,ngas+k),k=1,3)
!         N.B. does not yet properly handle across Grenwich Meridion	   
          call histwrt4(tr(1:ifull,:,ngas+1),'del_lat',idnc,iarch,local)
          call histwrt4(tr(1:ifull,:,ngas+2),'del_lon',idnc,iarch,local)
          call histwrt4(tr(1:ifull,:,ngas+3),'del_p',idnc,iarch,local)
         endif  ! (nextout>=4.and.nllp==3)
       endif    ! (mod(ktau,nperday)==0.or.ktau==ntau)
       if(mod(ktau,nperavg)==0.or.ktau==ntau)then 
!        only write these once per avg period
         call histwrt3(tscr_ave,'tscr_ave',idnc,iarch,local)
         call histwrt3(cbas_ave,'cbas_ave',idnc,iarch,local)
         call histwrt3(ctop_ave,'ctop_ave',idnc,iarch,local)
         call histwrt3(dew_ave,'dew_ave',idnc,iarch,local)
         call histwrt3(evap,'evap',idnc,iarch,local)
         call histwrt3(epan_ave,'epan_ave',idnc,iarch,local)
         call histwrt3(epot_ave,'epot_ave',idnc,iarch,local)
         call histwrt3(eg_ave,'eg_ave',idnc,iarch,local)
         call histwrt3(fg_ave,'fg_ave',idnc,iarch,local)
         call histwrt3(ga_ave,'ga_ave',idnc,iarch,local)
         call histwrt3(riwp_ave,'iwp_ave',idnc,iarch,local)
         call histwrt3(rlwp_ave,'lwp_ave',idnc,iarch,local)
         call histwrt3(cll_ave,'cll',idnc,iarch,local)
         call histwrt3(clm_ave,'clm',idnc,iarch,local)
         call histwrt3(clh_ave,'clh',idnc,iarch,local)
         call histwrt3(cld_ave,'cld',idnc,iarch,local)
       endif   ! (mod(ktau,nperavg)==0.or.ktau==ntau)
       call histwrt3(tscrn,'tscrn',idnc,iarch,local)
       call histwrt3(qgscrn,'qgscrn',idnc,iarch,local)
       call histwrt3(u10,'u10',idnc,iarch,local)
       call histwrt3(uscrn,'uscrn',idnc,iarch,local)
       call histwrt3(rnet,'rnet',idnc,iarch,local)
       call histwrt3(epan,'epan',idnc,iarch,local)
       call histwrt3(eg,'eg',idnc,iarch,local)
       call histwrt3(fg,'fg',idnc,iarch,local)
       call histwrt3(taux,'taux',idnc,iarch,local)
       call histwrt3(tauy,'tauy',idnc,iarch,local)
c      "extra" outputs
       if(nextout>=1) then
         if(myid == 0 ) print *,'nextout, idnc: ',nextout,idnc
         if(mod(ktau,nperavg)==0.or.ktau==ntau)then
           call histwrt3(rtu_ave,'rtu_ave',idnc,iarch,local)
           call histwrt3(rtc_ave,'rtc_ave',idnc,iarch,local)
           call histwrt3(rgdn_ave,'rgdn_ave',idnc,iarch,local)
           call histwrt3(rgn_ave,'rgn_ave',idnc,iarch,local)
           call histwrt3(rgc_ave,'rgc_ave',idnc,iarch,local)
           call histwrt3(sint_ave,'sint_ave',idnc,iarch,local)
           call histwrt3(sot_ave,'sot_ave',idnc,iarch,local)
           call histwrt3(soc_ave,'soc_ave',idnc,iarch,local)
           call histwrt3(sgdn_ave,'sgdn_ave',idnc,iarch,local)
           call histwrt3(sgn_ave,'sgn_ave',idnc,iarch,local)
         endif   ! (mod(ktau,nperavg)==0.or.ktau==ntau)
         call histwrt3(dpsdt,'dpsdt',idnc,iarch,local)
         call histwrt3(pblh,'pblh',idnc,iarch,local)
         call histwrt3(ustar,'ustar',idnc,iarch,local)
       endif   ! nextout>=1
      endif    ! (ktau>0.and.itype.ne.-1)

      if(myid == 0 ) print *,'netcdf save of 3d variables'
      call histwrt4(t(1:ifull,:),'temp',idnc,iarch,local)
      call histwrt4(u(1:ifull,:),'u',idnc,iarch,local)
      call histwrt4(v(1:ifull,:),'v',idnc,iarch,local)
      do k=1,kl
       do iq=1,ifull
        tmpry(iq,k)=ps(iq)*dpsldt(iq,k)
       enddo
      enddo
      call histwrt4(tmpry,'omega',idnc,iarch,local)  ! 3d variable
      call histwrt4(qg(1:ifull,:),'mixr',idnc,iarch,local)
      if(ldr.ne.0)then
        call histwrt4(qfg(1:ifullw,:),'qfg',idnc,iarch,local)
        call histwrt4(qlg(1:ifullw,:),'qlg',idnc,iarch,local)
        call histwrt4(cfrac,'cfrac',idnc,iarch,local)
      endif

!     rml 16/02/06 histwrt4 for trNNN and travNNN
      if(ngas>0)then 
       do igas=1,ngas
        write(trnum,'(i3.3)') igas
        call histwrt4(tr(1:ilt*jlt,:,igas),'tr'//trnum,idnc,iarch,local)
        call histwrt4(traver(:,:,igas),'trav'//trnum,idnc,iarch,local)
       enddo ! igas loop
      endif  ! (ngasc>0)

      if(itype==-1)then   ! extra stuff just needed for restart file
       call histwrt4(sdot(1,2),'sdot',idnc,iarch,local)
       call histwrt3(wbice(1,1),'wbice1',idnc,iarch,local)
       call histwrt3(wbice(1,2),'wbice2',idnc,iarch,local)
       call histwrt3(wbice(1,3),'wbice3',idnc,iarch,local)
       call histwrt3(wbice(1,4),'wbice4',idnc,iarch,local)
       call histwrt3(wbice(1,5),'wbice5',idnc,iarch,local)
       call histwrt3(wbice(1,6),'wbice6',idnc,iarch,local)
       call histwrt3(tggsn(1,1),'tggsn1',idnc,iarch,local)
       call histwrt3(tggsn(1,2),'tggsn2',idnc,iarch,local)
       call histwrt3(tggsn(1,3),'tggsn3',idnc,iarch,local)
       call histwrt3(smass(1,1),'smass1',idnc,iarch,local)
       call histwrt3(smass(1,2),'smass2',idnc,iarch,local)
       call histwrt3(smass(1,3),'smass3',idnc,iarch,local)
       call histwrt3(ssdn(1,1),'ssdn1',idnc,iarch,local)
       call histwrt3(ssdn(1,2),'ssdn2',idnc,iarch,local)
       call histwrt3(ssdn(1,3),'ssdn3',idnc,iarch,local)
       call histwrt3(snage,'snage',idnc,iarch,local)
       call histwrt3(rtsoil,'rtsoil',idnc,iarch,local) ! MJT cable       
       aa(:)=isflag(:)
       call histwrt3(aa,'sflag',idnc,iarch,local)
       if (nsib.eq.4.or.nsib.eq.6) call savetile(idnc,local,idim) ! MJT cable
      endif  ! (itype==-1)
      !---------------------------------------------------------
      ! MJT urban
      if ((nurban.eq.-1).or.((nurban.eq.1).and.(itype==-1))) then
       urban(:,:)=999. ! must be the same as spval in onthefly.f
       call atebsavem(ifull,urban,0)
       call histwrt3(urban(:,1),'rooftgg1',idnc,iarch,local)
       call histwrt3(urban(:,2),'rooftgg2',idnc,iarch,local)
       call histwrt3(urban(:,3),'rooftgg3',idnc,iarch,local)
       call histwrt3(urban(:,4),'waletgg1',idnc,iarch,local)
       call histwrt3(urban(:,5),'waletgg2',idnc,iarch,local)
       call histwrt3(urban(:,6),'waletgg3',idnc,iarch,local)
       call histwrt3(urban(:,7),'walwtgg1',idnc,iarch,local)
       call histwrt3(urban(:,8),'walwtgg2',idnc,iarch,local)
       call histwrt3(urban(:,9),'walwtgg3',idnc,iarch,local)
       call histwrt3(urban(:,10),'roadtgg1',idnc,iarch,local)
       call histwrt3(urban(:,11),'roadtgg2',idnc,iarch,local)
       call histwrt3(urban(:,12),'roadtgg3',idnc,iarch,local)
      end if
      !---------------------------------------------------------      
      
      !---------------------------------------------------------
      ! MJT CHANGE - Add wetfrac1-6 and possibly remove wb1-6 above
        aa(:)=(wb(:,1)-swilt(isoilm(:)))/
     &        (sfc(isoilm(:))-swilt(isoilm(:)))
        call histwrt3(aa,'wetfrac1',idnc,iarch,local)
        aa(:)=(wb(:,2)-swilt(isoilm(:)))/
     &        (sfc(isoilm(:))-swilt(isoilm(:)))
        call histwrt3(aa,'wetfrac2',idnc,iarch,local)
        aa(:)=(wb(:,3)-swilt(isoilm(:)))/
     &        (sfc(isoilm(:))-swilt(isoilm(:)))
        call histwrt3(aa,'wetfrac3',idnc,iarch,local)
        aa(:)=(wb(:,4)-swilt(isoilm(:)))/
     &        (sfc(isoilm(:))-swilt(isoilm(:)))
        call histwrt3(aa,'wetfrac4',idnc,iarch,local)
        aa(:)=(wb(:,5)-swilt(isoilm(:)))/
     &        (sfc(isoilm(:))-swilt(isoilm(:)))
        call histwrt3(aa,'wetfrac5',idnc,iarch,local)
        aa(:)=(wb(:,6)-swilt(isoilm(:)))/
     &        (sfc(isoilm(:))-swilt(isoilm(:)))
        call histwrt3(aa,'wetfrac6',idnc,iarch,local)       
      !---------------------------------------------------------

      return
      end
c=======================================================================
      subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax,daily)

      include 'netcdf.inc'

      integer*2 minv, maxv, missval   ! was integer*2
      parameter(minv = -32500, maxv = 32500, missval = -32501)
      integer cdfid, idv, dim(3), daily
      character name*(*), lname*(*), units*(*)
      real xmin, xmax
      integer, parameter :: vtype = ncshort

      idv = ncvdef(cdfid, name, vtype, ndim, dim, ier)
      if(ier.ne.0)then
        write(0,*) ier,' Error in variable declaration ', name
        stop
      endif

      call ncaptc(cdfid,idv,'long_name',ncchar,len_trim(lname),lname,
     &            ier)
      if(len_trim(units).ne.0)then
        call ncaptc(cdfid,idv,'units',ncchar,len_trim(units),units,ier)
      endif
      if(vtype == ncshort)then
        call ncapt(cdfid,idv,'valid_min'    ,ncshort,1,minv,ier)
        call ncapt(cdfid,idv,'valid_max'    ,ncshort,1,maxv,ier)
        call ncapt(cdfid,idv,'missing_value',ncshort,1,missval,ier)
!       scalef=(xmax-xmin)/float(maxv - minv)
        scalef=(xmax-xmin)/(real(maxv)-real(minv)) ! jlm fix for precision problems
        addoff=xmin-scalef*minv
        call ncapt(cdfid,idv,'add_offset',ncfloat,1,addoff,ier)
        call ncapt(cdfid,idv,'scale_factor',ncfloat,1,scalef,ier)
      endif
      call ncaptc(cdfid,idv,'FORTRAN_format',ncchar,5,'G11.4',ier)
      if(daily>0)then
        call ncaptc(cdfid,idv,'valid_time',ncchar,5,'daily',ier) 
      endif
      return
      end
c=======================================================================
      subroutine attrl(cdfid,dim,ndim,name,lname,units,xmin,xmax,daily)
!     this one for long variables, e.g. snd
      include 'netcdf.inc'
      integer cdfid, idv, dim(3), daily
      character name*(*), lname*(*), units*(*)

      idv = ncvdef(cdfid, name,NCFLOAT, ndim, dim, ier)
      if(ier.ne.0)then
        write(0,*) ier,' Error in variable declaration ', name
        stop
      endif

      call ncaptc(cdfid,idv,'long_name',ncchar,len_trim(lname),lname,
     &            ier)
      if(len_trim(units).ne.0)then
        call ncaptc(cdfid,idv,'units',ncchar,len_trim(units),units,ier)
      endif
      call ncaptc(cdfid,idv,'FORTRAN_format',ncchar,5,'G11.4',ier)
      if(daily>0)then
        call ncaptc(cdfid,idv,'valid_time',ncchar,5,'daily',ier) 
      endif
      return
      end
c=======================================================================
      subroutine histwrt3(var,sname,idnc,iarch,local)
c Write 2d+t fields from the savegrid array.

      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'netcdf.inc'

      integer idnc, iarch
      logical, intent(in) :: local
      integer mid, start(3), count(3)
      integer*2 ipack(ifull_g) ! was integer*2 
      character* (*) sname
c     character*8 sname
      integer*2 minv, maxv, missval ! was integer*2 
      parameter(minv = -32500, maxv = 32500, missval = -32501)

      real var(ifull)
      integer iq, ier, imn, imx, jmn, jmx, vtype
      real addoff, pvar, scale_f, varn, varx, xmax, xmin
      real, dimension(ifull_g) :: globvar

      if(local)then
         start = (/ 1, 1, iarch /)
         count = (/ il, jl, 1 /)
         mid = ncvid(idnc,sname,ier)

!        Check variable type
         ier = nf_inq_vartype(idnc, mid, vtype)
         if(vtype == ncshort)then
            call ncagt(idnc,mid,'add_offset',addoff,ier)
            call ncagt(idnc,mid,'scale_factor',scale_f,ier)

            xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
            xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems

            do iq=1,ifull
               pvar = max(xmin,min(xmax,var(iq)))
               ipack(iq)=nint((pvar-addoff)/scale_f)
               ipack(iq)=max(min(ipack(iq),maxv),minv)
            end do

            call ncvpt(idnc, mid, start, count, ipack, ier)
         else
            call ncvpt(idnc, mid, start, count, var, ier)
         endif
         if(ier.ne.0)then
           write(0,*) "in histwrt3 ier not zero",ier
           stop
         endif
      else

      if(myid == 0)then
         call ccmpi_gather(var, globvar)
         start(1) = 1
         start(2) = 1
         start(3) = iarch
         count(1) = il_g
         count(2) = jl_g
         count(3) = 1

c find variable index
         mid = ncvid(idnc,sname,ier)

!        Check variable type
         ier = nf_inq_vartype(idnc, mid, vtype)
         if(vtype == ncshort)then
            call ncagt(idnc,mid,'add_offset',addoff,ier)
            call ncagt(idnc,mid,'scale_factor',scale_f,ier)

            xmin=addoff+scale_f*minv
!           xmax=xmin+scale_f*float(maxv-minv)
            xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems

            do iq=1,ifull_g
               pvar = max(xmin,min(xmax,globvar(iq)))
               ipack(iq)=nint((pvar-addoff)/scale_f)
               ipack(iq)=max(min(ipack(iq),maxv),minv)
            end do

            call ncvpt(idnc, mid, start, count, ipack, ier)
         else
            call ncvpt(idnc, mid, start, count, globvar, ier)
         endif
         if(ier.ne.0)then
           write(0,*) "In histwrt3 ier not zero",ier
           stop
         endif
      else
         call ccmpi_gather(var)
      endif  ! (myid == 0) .. else ..

c     print *,'myid,ktau,nmaxpr,sname ',myid,ktau,nmaxpr,sname
      if(myid==0 .and. mod(ktau,nmaxpr)==0)then
         varn = minval(globvar)
         varx = maxval(globvar)
         ! This should work ???? but sum trick is more portable???
         ! iq = minloc(globvar,dim=1)
         iq = sum(minloc(globvar))
         ! Convert this 1D index to 2D
         imn = 1 + modulo(iq-1,il_g)
         jmn = 1 + (iq-1)/il_g
         iq = sum(maxloc(globvar))
         ! Convert this 1D index to 2D
         imx = 1 + modulo(iq-1,il_g)
         jmx = 1 + (iq-1)/il_g
         write(6,'("histwrt3 ",a7,i4,f12.4,2i4,f12.4,2i4,f12.4)')
     &             sname,iarch,varn,imn,jmn,varx,imx,jmx,
     &             globvar(id+(jd-1)*il_g)
      endif
      endif ! local  .. else ..

      return
      end ! histwrt3
c=======================================================================
      subroutine histwrt3l(var,sname,idnc,iarch,local)  ! long write
c Write 2d+t fields from the savegrid array.  long write (e.g. snd)

      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'netcdf.inc'

      integer idnc, iarch
      logical, intent(in) :: local
      integer mid, start(3), count(3)
      character* (*) sname
c     character*8 sname

      real var(ifull)
      integer iq, ier, imn, imx, jmn, jmx
      real  varn, varx
      real, dimension(ifull_g) :: globvar

      if(local)then
         start = (/ 1, 1, iarch /)
         count = (/ il, jl, 1 /)
         mid = ncvid(idnc,sname,ier)

!        Check variable type
         ier = nf_inq_vartype(idnc, mid, NCFLOAT)
         call ncvpt(idnc, mid, start, count, var, ier)
         if(ier.ne.0)then
           write(0,*) "in histwrt3l ier not zero"
           stop
         endif
      else

      if(myid == 0)then
         call ccmpi_gather(var, globvar)
         start(1) = 1
         start(2) = 1
         start(3) = iarch
         count(1) = il_g
         count(2) = jl_g
         count(3) = 1

c find variable index
         mid = ncvid(idnc,sname,ier)

!        Check variable type
c        ier = nf_inq_vartype(idnc, mid, NCFLOAT) ! removed 13/6/07
         call ncvpt(idnc, mid, start, count, globvar, ier)
         if(ier.ne.0)then
           write(0,*) "In histwrt3l ier not zero"
           stop
         endif
      else
         call ccmpi_gather(var)
      endif

      if(myid==0 .and. mod(ktau,nmaxpr)==0)then
         varn = minval(globvar)
         varx = maxval(globvar)
         ! This should work ???? but sum trick is more portable???
         ! iq = minloc(globvar,dim=1)
         iq = sum(minloc(globvar))
         ! Convert this 1D index to 2D
         imn = 1 + modulo(iq-1,il_g)
         jmn = 1 + (iq-1)/il_g
         iq = sum(maxloc(globvar))
         ! Convert this 1D index to 2D
         imx = 1 + modulo(iq-1,il_g)
         jmx = 1 + (iq-1)/il_g
         write(6,'("histwrt3l",a7,i4,f12.4,2i4,f12.4,2i4,f12.4)')
     &             sname,iarch,varn,imn,jmn,varx,imx,jmx,
     &             globvar(id+(jd-1)*il_g)
      endif
      endif ! local

      return
      end ! histwrt3l
c=======================================================================
      subroutine histwrt4(var,sname,idnc,iarch,local)
c Write 3d+t fields from the savegrid array.

      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'netcdf.inc'

      integer idnc, iarch
      logical, intent(in) :: local
      integer mid, start(4), count(4)
      integer*2 ipack(ifull_g,kl) ! was integer*2 
      character* (*) sname
c     character*8 sname
      integer*2 minv, maxv, missval ! was integer*2 
      parameter(minv = -32500, maxv = 32500, missval = -32501)

      real var(ifull,kl)
      integer ier, imx, jmx, k, kmx, iq, vtype
      real addoff, pvar, scale_f, varn, varx, xmax, xmin
      real, dimension(ifull_g,kl) :: globvar
      integer, dimension(2) :: max_result

      if(local)then
         start = (/ 1, 1, 1, iarch /)
         count = (/ il, jl, kl, 1 /)

c find variable index
         mid = ncvid(idnc,sname,ier)
!        Check variable type
         ier = nf_inq_vartype(idnc, mid, vtype)
         if(vtype == ncshort)then
            call ncagt(idnc,mid,'add_offset',addoff,ier)
            call ncagt(idnc,mid,'scale_factor',scale_f,ier)

            xmin=addoff+scale_f*minv
!           xmax=xmin+scale_f*float(maxv-minv)
            xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems
            do k=1,kl
               do iq=1,ifull_g
                  pvar = max(xmin,min(xmax,var(iq,k)))
                  ipack(iq,k)=nint((pvar-addoff)/scale_f)
                  ipack(iq,k)=max(min(ipack(iq,k),maxv),minv)
               end do
            end do
            call ncvpt(idnc, mid, start, count, ipack, ier)
         else
            call ncvpt(idnc, mid, start, count, var, ier)
         endif

      else ! not local
      if(myid == 0)then
         call ccmpi_gather(var, globvar)
         start(1) = 1
         start(2) = 1
         start(3) = 1
         start(4) = iarch
         count(1) = il_g
         count(2) = jl_g
         count(3) = kl
         count(4) = 1

c find variable index
         mid = ncvid(idnc,sname,ier)
!        Check variable type
         ier = nf_inq_vartype(idnc, mid, vtype)
         if(vtype == ncshort)then
            call ncagt(idnc,mid,'add_offset',addoff,ier)
            call ncagt(idnc,mid,'scale_factor',scale_f,ier)

            xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
            xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems
            do k=1,kl
               do iq=1,ifull_g
                  pvar = max(xmin,min(xmax,globvar(iq,k)))
                  ipack(iq,k)=nint((pvar-addoff)/scale_f)
                  ipack(iq,k)=max(min(ipack(iq,k),maxv),minv)
               end do
            end do
            call ncvpt(idnc, mid, start, count, ipack, ier)
         else
            call ncvpt(idnc, mid, start, count, globvar, ier)
         endif
      else
         call ccmpi_gather(var)
      endif

      if(myid==0 .and. mod(ktau,nmaxpr)==0)then
         varn = minval(globvar)
         varx = maxval(globvar)
         max_result = maxloc(globvar)
         kmx = max_result(2)
         iq = max_result(1)
         ! Convert this 1D index to 2D
         imx = 1 + modulo(iq-1,il_g)
         jmx = 1 + (iq-1)/il_g
         write(6,'("histwrt4 ",a7,i4,2f12.4,3i4,f12.4)')
     &     sname,iarch,varn,varx,imx,jmx,kmx,globvar(id+(jd-1)*il_g,nlv)
      endif
      endif ! local

      return
      end ! histwrt4
     
      subroutine mtimerget(mtimer,kdate1,ktime1,kdate2,ktime2) ! jlm
!     returns mtimer in minutes, corr. to (kdate2,ktime2) -  (kdate1,ktime1)    
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      common/leap_yr/leap  ! 1 to allow leap years
 
      if(leap.ne.0)then     
        write(0,*) 'leap years not catered for in mtimerget'
        stop
      endif
!     Set up number of minutes from beginning of year
!     For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
      jyear1=kdate1/10000
      jmonth=(kdate1-jyear1*10000)/100
      jday=kdate1-jyear1*10000-jmonth*100
      jhour=ktime1/100
      jmin=ktime1-jhour*100
      mstart1=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      jyear2=kdate2/10000
      jmonth=(kdate2-jyear2*10000)/100
      jday=kdate2-jyear2*10000-jmonth*100
      jhour=ktime2/100
      jmin=ktime2-jhour*100
      mstart2=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      mtimer=mstart2-mstart1+(jyear2-jyear1)*365*24*60
      return
      end
ce=======================================================================
      subroutine outcdfs(rundate)  ! for (hourly) scrnfile
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg,mtimer
      include 'filnames.h'  ! list of files, read in once only
      include 'kuocom.h'
      include 'liqwpar.h'  ! ifullw
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmgeom.h' ! rlong0,rlat0,schmidt  
      include 'parmhor.h'  ! mhint, m_bs, nt_adv, ndept
      include 'parmvert.h'
      include 'tracers.h'  ! ngas, nllp, ntrac, tr
      character rundate*10
 
      integer nhor,nhorps,khor,khdif,nhorjlm
      real hdiff,hdifmax
      common/parmhdff/nhor,nhorps,hdiff(kl),khor,khdif,hdifmax,nhorjlm
      integer, parameter :: nihead=54
      integer nahead(nihead)

      integer, parameter :: nrhead=14
      real ahead(nrhead)

      include 'netcdf.inc'
      character cdffile*80

      integer ixp,iyp,idlev,idnt,idms
      common/cdfind/ixp,iyp,idlev,idnt,idms

      integer dim(4),dims(4)
      integer xdim,ydim,zdim,tdim,msdim
      character timorg*20
      character grdtim*33
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/

      integer :: ndt, icy, icm, icd, ich, icmi, ics, idv, ier, imode
      integer, save :: idnc
      integer, save :: iarch=0, idnc0=0
      logical :: local ! Each processor writes its local region

      ndt=dt
      local = localhist 
      if(myid==0.or. local)then !  #########################
!      File setup follows
       iarch=iarch+1
       cdffile=scrnfile
       write(6,'("outcdfs idnc,iarch,cdffile=",2i5," ",a80)')
     &                    idnc,iarch,cdffile

       if(iarch==1)then
        print *,'nccre of ',cdffile
        idnc = nccre(cdffile, ncclob, ier)
        print *,'idnc,ier=',idnc,ier
c       Turn off the data filling
        imode = ncsfil(idnc,ncnofill,ier)
        print *,'imode=',imode
c       Create dimensions, lon, lat
        if(local)then
           xdim = ncddef(idnc, 'longitude', il, ier)
           ydim = ncddef(idnc, 'latitude', jl, ier)
        else
           xdim = ncddef(idnc, 'longitude', il_g, ier)
           ydim = ncddef(idnc, 'latitude', jl_g, ier)
        endif
        zdim= ncddef(idnc, 'lev', kl, ier)
        msdim= ncddef(idnc, 'zsoil', ms, ier)
        tdim= ncddef(idnc, 'time',ncunlim,ier)
        print *,"xdim,ydim,zdim,tdim"
        print *,xdim,ydim,zdim,tdim

c       define coords.
        ixp = ncvdef(idnc,'longitude',NCFLOAT,1,xdim,ier)
        call ncaptc(idnc,ixp,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnc,ixp,'units',NCCHAR,12,'degrees_east',ier)
        iyp = ncvdef(idnc,'latitude',NCFLOAT,1,ydim,ier)
        call ncaptc(idnc,iyp,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnc,iyp,'units',NCCHAR,13,'degrees_north',ier)
        print *,'ixp,iyp=',ixp,iyp

        idlev = ncvdef(idnc,'lev',NCFLOAT,1,zdim,ier)
        call ncaptc(idnc,idlev,'positive',NCCHAR,4,'down',ier)
        call ncaptc(idnc,idlev,'point_spacing',NCCHAR,6,'uneven',ier)
        call ncaptc(idnc,idlev,'units',NCCHAR,11,'sigma_level',ier)
        call ncaptc(idnc,idlev,'long_name',NCCHAR,11,'sigma_level',ier)
        print *,'idlev=',idlev

        idms = ncvdef(idnc,'zsoil',NCFLOAT,1,msdim,ier)
        call ncaptc(idnc,idms,'point_spacing',NCCHAR,6,'uneven',ier)
        call ncaptc(idnc,idms,'units',NCCHAR,1,'m',ier)
        print *,'idms=',idms

        print *,'tdim,idnc=',tdim,idnc
ccc     idnt = ncvdef(idnc,'time',NCLONG,1,tdim,ier)
        idnt = ncvdef(idnc,'time',NCFLOAT,1,tdim,ier)  ! Sept 2006
        print *,'idnt=',idnt
        call ncaptc(idnc,idnt,'point_spacing',NCCHAR,4,'even',ier)

        print *,'kdate,ktime,ktau=',kdate,ktime,ktau
        icy=kdate/10000
        icm=max(1,min(12,(kdate-icy*10000)/100))
        icd=max(1,min(31,(kdate-icy*10000-icm*100)))
        if(icy<100)icy=icy+1900
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(6,*) icy,icm,icd,ich,icmi,ics
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &               icd,month(icm),icy,ich,icmi,ics
        print *,'timorg=',timorg
        call ncaptc(idnc,idnt,'time_origin',NCCHAR,20,timorg,ier)

        write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",
     &       2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
        print *,'grdtim=',grdtim
        call ncaptc(idnc,idnt,'units',NCCHAR,33,grdtim,ier)

        dim(1) = xdim
        dim(2) = ydim
        dim(3) = zdim
        dim(4) = tdim

        dims(1) = xdim
        dims(2) = ydim
        dims(3) = msdim
        dims(4) = tdim

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
        nahead(13)=0
        nahead(14)=ndt       ! needed by cc2hist
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
        write(6,'("nahead=",(20i4))') nahead
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
        write(6,*) "ahead=",ahead
        call ncapt(idnc,ncglobal,'int_header',nclong,nihead,nahead,ier)
        if(ier.ne.0)write(6,*)"ncapt int idnc,ier=",idnc,ier
        call ncapt(idnc,ncglobal,'real_header',ncfloat,nrhead,ahead,ier)
        if(ier.ne.0)write(6,*)"ncapt real idnc,ier=",idnc,ier
        call ncaptc(idnc,ncglobal,'date_header',ncchar,10,rundate,ier)
        if(ier.ne.0)write(6,*)"ncaptc date idnc,ier=",idnc,ier

        idv=ncvdef(idnc,'ds',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef ds idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'du',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef du idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'rnml',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef rnml idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'tanl',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef tanl idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'stl1',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl1 idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'stl2',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl2 idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'dt',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef dt idnc,ier=",idnc,ier
       endif ! ( iarch=1)then

      endif ! (myid==0.or.local) #########################
      ! openhist writes some fields so needs to be called by all processes
      call openhists(iarch,dim,local,idnc)
c     print *,'after openhists for myid = ',myid

      if(myid==0.or.local)then
        call ncsnc(idnc,ier)
        if(ier.ne.0)write(6,*)"ncsnc idnc,ier=",idnc,ier
      endif    ! (myid==0.or.local)
      if(ktau.eq.ntau.and.(myid==0.or.localhist))then
        call ncclos(idnc,ier)
        write(6,*) "calling ncclos(idnc,ier) ",idnc,ier
      endif
      return   ! outcdfs
      end
c=======================================================================
      subroutine openhists(iarch,dim,local,idnc)
      use cc_mpi
      implicit none
c     this routine creates attributes and writes output
      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'
      include 'dates.h'    ! ktime,kdate,timer,timeg,xg,yg,mtimer
      include 'extraout.h' ! u10_3hr,v10_3hr
      include 'filnames.h' ! list of files, read in once only
      include 'histave.h'
      include 'kuocom.h'
      include 'liqwpar.h'  ! ifullw
      include 'map.h'
      include 'mapproj.h'
      include 'morepbl.h'
      include 'netcdf.inc'
      include 'nsibd.h' ! rsmin,ivegt,sigmf,tgg,tgf,ssdn,res,rmc,isoilm,ico2em
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmvert.h'
      include 'pbl.h'
      include 'prec.h'
      include 'raddiag.h'
      include 'scamdim.h'
      include 'screen.h'
      include 'sigs.h'
      include 'soil.h'
      include 'soilsnow.h'
      include 'soilv.h'   ! sfc,zse
      include 'tracers.h'
      include 'trcom2.h'
      include 'version.h'
      include 'vvel.h'    ! sdot, dpsldt

      integer iarch
      logical, intent(in) :: local
      character lname*40,expdesc*50
      integer dim(4)
      integer idim(3)
      real xpnt(il_g),ypnt(jl_g)

      integer ixp,iyp,idlev,idnt,idms
      common/cdfind/ixp,iyp,idlev,idnt,idms
      real pmsl,aa,bb,cc,dum2
      common/work2/pmsl(ifull),aa(ifull),bb(ifull),cc(ifull),
     &             dum2(ifull,14)
      real tmpry
      common/work3c/tmpry(ifull,kl)

      integer i, idkdate, idktau, idktime, idmtimer, idnteg, idnter,
     &     idv, ier, j, idnc
      character*3 mon(12)
      data mon/'JAN','FEB','MAR','APR','MAY','JUN'
     &        ,'JUL','AUG','SEP','OCT','NOV','DEC'/

      if(myid == 0 .or. local)then  !#########################
       print *,'openhists iarch,idnc=',iarch,idnc

c      if this is the first archive, set up some global attributes
       if(iarch==1) then
        print *,'dim=',dim
        idim(1)=dim(1)
        idim(2)=dim(2)
        idim(3)=dim(4)
        print *,'idim=',idim

c       Create global attributes
c       Model run number
        print *,'nrun=',nrun
        call ncapt(idnc,ncglobal,'nrun',nclong,1,nrun,ier)
        write(6,*)"nrun ier=",ier,idnc

c       Experiment description
        expdesc = 'CCAM model run'
        call ncaptc(idnc,ncglobal,'expdesc',ncchar,len_trim(expdesc),
     &              expdesc,ier)
        write(6,*)"expdesc ier=",ier,idnc

c       Model version
        call ncaptc(idnc,ncglobal,'version',ncchar,len_trim(version),
     &              version,ier)

        if(local)then
           ier = nf_put_att_int(idnc,nf_global,"processor_num",nf_int,
     &                          1,myid)
        endif     

c       Sigma levels
        print *,'sig=',sig
        call ncapt(idnc,ncglobal,'sigma',ncfloat,kl,sig,ier)

        lname = 'timer (hrs)'
        idnter = ncvdef(idnc,'timer',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnter,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'mtimer (mins)'
        idmtimer = ncvdef(idnc,'mtimer',nclong,1,dim(4),ier)
        call ncaptc(idnc,idmtimer,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'timeg (UTC)'
        idnteg = ncvdef(idnc,'timeg',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnteg,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'number of time steps from start'
        idktau = ncvdef(idnc,'ktau',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktau,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'year-month-day at start of run'
        idkdate = ncvdef(idnc,'kdate',nclong,1,dim(4),ier)
        call ncaptc(idnc,idkdate,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        lname = 'hour-minute at start of run'
        idktime = ncvdef(idnc,'ktime',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktime,'long_name',ncchar
     &             ,len_trim(lname),lname,ier)

        idv = ncvdef(idnc,'sigma', ncfloat, 1, dim(3),ier)
        call ncaptc(idnc,idv,'positive',ncchar
     &             ,len_trim('down'),'down',ier)

        print *,'define attributes of variables for scrnfile'
c       For time varying surface fields
        lname = 'Screen temperature'
        call attrib(idnc,idim,3,'tscrn',lname,'K',100.,400.,0)
        lname = 'Screen mixing ratio'
        call attrib(idnc,idim,3,'qgscrn',lname,'kg/kg',0.,.06,0)
        lname = 'Cape'
        call attrib(idnc,idim,3,'cape',lname,'J/kg',0.,20000.,0) ! 1 was daily
        lname = 'max 10m wind'
        call attrib(idnc,idim,3,'u10mx',lname,'m/s',0.,99.,0) ! 1 was daily attr
        lname = '10m wind speed'
        call attrib(idnc,idim,3,'u10',lname,'m/s',0.,60.,0)
        lname = 'u level 1'
        call attrib(idnc,idim,3,'u1',lname,'m/s',-60.,60.,0)
        lname = 'v level 1'
        call attrib(idnc,idim,3,'v1',lname,'m/s',-60.,60.,0)
        lname = 'u level 2'
        call attrib(idnc,idim,3,'u2',lname,'m/s',-60.,60.,0)
        lname = 'v level 2'
        call attrib(idnc,idim,3,'v2',lname,'m/s',-60.,60.,0)
        lname = 'v level 3'
        call attrib(idnc,idim,3,'u3',lname,'m/s',-60.,60.,0)
        lname = 'u level 3'
        call attrib(idnc,idim,3,'v3',lname,'m/s',-60.,60.,0)

        print *,'finished defining attributes'
c       Leave define mode
        call ncendf(idnc,ier)
        print *,'leave define mode: ier=',ier

        if(local)then
           ! Set these to global indices
           do i=1,ipan
              xpnt(i) = float(i) + ioff
           end do
           call ncvpt(idnc,ixp,1,il,xpnt,ier)
           do j=1,jl
              ypnt(j) = float(j) + jl*((myid*il)/il_g)
           end do
           call ncvpt(idnc,iyp,1,jl,ypnt,ier)
        else
           do i=1,il_g
              xpnt(i) = float(i)
           end do
           call ncvpt(idnc,ixp,1,il_g,xpnt,ier)
           do j=1,jl_g
              ypnt(j) = float(j)
           end do
           call ncvpt(idnc,iyp,1,jl_g,ypnt,ier)
        endif

        call ncvpt(idnc,idlev,1,kl,sig,ier)

        idv = ncvid(idnc,'sigma',ier)
        call ncvpt(idnc,idv,1,kl,sig,ier)

        idv = ncvid(idnc,'lev',ier)
        call ncvpt(idnc,idv,1,kl,sig,ier)

        idv = ncvid(idnc,'ds',ier)
        call ncvpt1(idnc,idv,1,ds,ier)
        call ncvpt1(idnc,idv,1,stl2,ier)
        idv = ncvid(idnc,'dt',ier)
        call ncvpt1(idnc,idv,1,dt,ier)
       endif ! iarch==1
!      -----------------------------------------------------------      
       print *,'outcdfs processing kdate,ktime,ktau,mtimer: ',
     .                            kdate,ktime,ktau,mtimer
c      set time to number of minutes since start 
       idv = ncvid(idnc,'time',ier)
ccc    call ncvpt1(idnc,idv,iarch,mtimer,ier)
       call ncvpt1(idnc,idv,iarch,real(mtimer),ier)  ! Sept 2006
       
       idv = ncvid(idnc,'timer',ier)
       call ncvpt1(idnc,idv,iarch,timer,ier)
       idv = ncvid(idnc,'mtimer',ier)
       call ncvpt1(idnc,idv,iarch,mtimer,ier)
       idv = ncvid(idnc,'timeg',ier)
       call ncvpt1(idnc,idv,iarch,timeg,ier)
       idv = ncvid(idnc,'ktau',ier)
       call ncvpt1(idnc,idv,iarch,ktau,ier)
       idv = ncvid(idnc,'kdate',ier)
       call ncvpt1(idnc,idv,iarch,kdate,ier)
       idv = ncvid(idnc,'ktime',ier)
       call ncvpt1(idnc,idv,iarch,ktime,ier)
       print *,'kdate,ktime,ktau=',kdate,ktime,ktau
       print *,'timer,timeg=',timer,timeg
       print *,'now write out variables'
      endif ! (myid == 0 .or. local)

      call histwrt3(tscrn,'tscrn',idnc,iarch,local)
      call histwrt3(qgscrn,'qgscrn',idnc,iarch,local)
      call histwrt3(cape,'cape',idnc,iarch,local)
      call histwrt3(u10,'u10',idnc,iarch,local)
      call histwrt3(u10mx,'u10mx',idnc,iarch,local)
      call histwrt3(u(1:ifull,1),'u1',idnc,iarch,local)
      call histwrt3(v(1:ifull,1),'v1',idnc,iarch,local)
      call histwrt3(u(1:ifull,2),'u2',idnc,iarch,local)
      call histwrt3(v(1:ifull,2),'v2',idnc,iarch,local)
      call histwrt3(u(1:ifull,3),'u3',idnc,iarch,local)
      call histwrt3(v(1:ifull,3),'v3',idnc,iarch,local)

      return   ! openhists
      end
