      subroutine outfile(iout,ik,jk,kk,psa,psm,rundate,nmi,nsnowout,
     .                   nwrite)
      use cc_mpi
      include 'newmpar.h'
      parameter (mev1=il+1-il/2*2,mev2=3-mev1)  ! from helmsol
c     mev1 = 1 for il even (2 for il odd)
c     mev2 = 2 for il even (1 for il odd)
      parameter (nwrite0=1) ! 0 original, 1 writes initial nveg,nsoil etc
      include 'aalat.h'
      include 'arrays.h'
      include 'darcdf.h'   ! idnc,ncid,idifil  - stuff for netcdf
      include 'dates.h'    ! mtimer
      include 'dava.h' ! davt    for writing out in precc
      include 'extraout.h'    ! taux,tauy,rtsave,rtclsave
      include 'filnames.h' ! list of files, read in once only
      include 'histave.h'
      include 'kuocom.h'
      include 'map.h'     ! passes thru omgf/ps in dpsldt array
      include 'morepbl.h'  ! fg,eg,runoff
      include 'nlin.h'
      include 'nsibd.h'
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmvert.h'
      include 'pbl.h'
      include 'prec.h'
      include 'scamdim.h' ! npmax
      include 'screen.h'  
      include 'sigs.h'
      include 'soil.h'  ! land,sice,sicedep,alb
      include 'soilv.h'  ! swilt,sfc
      include 'soilsnow.h'  ! tgg,wb,snowd
      include 'tracers.h'
      include 'vvel.h'
      common/work2/pmsl(ifull),tssout(ifull),seaice(ifull)
     .      ,aa(ifull),bb(ifull),cc(ifull),dum2(ifull,12)
      real omgf(ifull,kl)
      equivalence (omgf,dpsldt)
      dimension psa(*),psm(*)
      character rundate*10,formout*9,fileout*10,qgout*20
      character co2out*80,radonout*80,surfout*80
      integer nface6(4,3)   ! Faces to use in each phase    (c-cub)
      integer ioff6(4,3)    ! Starting offset for each face (c-cub)
      data nface6 / 0, 1, 3, 4,   0, 2, 3, 5,   1, 2, 4, 5 /
      data ioff6 / 1,mev1,mev2,2,  2,1,mev1,mev2,  mev2,2,1,mev1 / ! jlm general
      data nspare/0/

!!!      if(ktau.eq.0.and.io_rest.ne.0)then ! switch for cctocc4
!!!        print *,'in outfile opening surfile: ',surfile
!!!c       write surface information file at beginning of run
!!!        open(unit=77,file=surfile,form='unformatted',status='unknown')
!!!        do iq=1,ifull
!!!         aa(iq)=ivegt(iq)+.01
!!!         bb(iq)=isoilm(iq)+.01
!!!        enddo   ! iq loop
!!!        write (77) aa  ! ivegt
!!!        write (77) bb  ! isoilm
!!!        write (77) alb ! just starting albedo 
!!!        write (77) rsmin
!!!        write (77) zolnd
!!!        close (77)
!!!        do iphase = 1,3         !  mrd code to run fast on NEC
!!!         do iface=1,4
!!!          if = nface6(iface,iphase)
!!!          istart = ioff6(iface,iphase) ! 1 or 2
!!!          do j=1,il
!!!           istart = 3 - istart ! 1 or 2 alternately
!!!           do i=istart,il,2
!!!            iq=ind(i,j,if)
!!!            cc(iq)=iphase+.5               ! 3-colouring of grid cells
!!!            if(land(iq))cc(iq)=cc(iq)+3. ! bigger for land points
!!!           end do
!!!          end do
!!!         end do
!!!        end do
!!!      endif  ! (ktau.eq.0.and.io_rest.ne.0)
!!!      if(ktau.eq.nsnowout)then
!!!       write (fileout,'(6hsnowd. ,i4)')nrun  !  i.e. snowd.xxxx
!!!       open(unit=77,file=fileout,status='unknown',form='formatted')
!!!       write(77,'(i3,i4,2f6.1,f5.2,f9.0,''  snowd from run'',i5)')
!!!     .                        ifull,rlong0,rlat0,schmidt,ds,nrun
!!!       write (formout,'(1h(,i3,4hf8.1,1h))')il  !  i.e. (<il>f8.1)
!!!       write(77,formout) snowd
!!!       close(77)
!!!       print *,'snowout written for ktau =',ktau
!!!      endif  ! (ktau.eq.nsnowout)

      call start_log(outfile_begin)
      ndt=dt
      io_outt=io_out
      if(iout.eq.19)io_outt=io_rest  !  choice for writing restart file
      if ( myid==0 ) then
         print *,'ofile written for iout,kdate,ktime,mtimer: ',
     &                              iout,kdate,ktime,mtimer
      end if

      call mslp(pmsl,psl,zs,t(1:ifull,:))
c     reincorporate land mask into surface temperature
      do iq=1,ifull
       tssout(iq)=tss(iq)
       if(.not.land(iq))tssout(iq)=-tss(iq)
      enddo   ! iq loop
      if ( mydiag ) print *,'tssout: ',tssout(idjd)

      if(io_outt.eq.3)then
         print*, "Error, binary output not supported"
         stop
      end if
c---------------------------------------------------------------------------
      if(io_outt.eq.1)then  ! for netcdf
         ms_out = 0             ! Not set anywhere else ?????
         if(iout.eq.19)then
            if ( myid==0 ) print *,"restart write of data to cdf"
            call outcdf(rundate,nmi,-1,ms_out)
            if ( myid==0 ) then
               call ncclos(idnc,ier)
               write(6,*) "call ncclos(idnc,ier) ",idnc,ier
            end if
            call end_log(outfile_end)
            return              ! done with restart data
         endif                  !  (iout.eq.19)
         if ( myid==0 ) print *,'calling outcdf from outfile'
         call outcdf(rundate,nmi,1,ms_out)  

         if(ktau.eq.ntau.and.(myid==0.or.localhist))then
            call ncclos(idnc,ier)
            write(6,*) "calling ncclos(idnc,ier) ",idnc,ier
         endif
      endif ! (io_outt.eq.1)

      call end_log(outfile_end)
      return
      end
