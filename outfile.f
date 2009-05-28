      subroutine outfile(iout,rundate,nmi,nwrite)
      use cc_mpi
      include 'newmpar.h'
      parameter (mev1=il+1-il/2*2,mev2=3-mev1)  ! from helmsol
c     mev1 = 1 for il even (2 for il odd)
c     mev2 = 2 for il even (1 for il odd)
      parameter (nwrite0=1) ! 0 original, 1 writes initial nveg,nsoil etc
      include 'arrays.h'
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
      character rundate*10,qgout*20
      character co2out*80,radonout*80,surfout*80
      integer nface6(4,3)   ! Faces to use in each phase    (c-cub)
      integer ioff6(4,3)    ! Starting offset for each face (c-cub)
      data nface6 / 0, 1, 3, 4,   0, 2, 3, 5,   1, 2, 4, 5 /
      data ioff6 / 1,mev1,mev2,2,  2,1,mev1,mev2,  mev2,2,1,mev1 / ! jlm general
      data nspare/0/

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

      if(nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)then
         if(ktau.eq.nwrite/2.or.ktau.eq.nwrite)then
!        usually after first 24 hours, save soil variables for next run
            if(ktau.eq.nwrite)then  ! 24 hour write
               if(ktime.eq.1200)then
                  co2out=co2_12     ! 'co2.1200'
                  radonout=radon_12 ! 'radon.1200'
                  surfout=surf_12   ! 'current.1200'
                  qgout='qg_12'
               else
                  co2out=co2_00     !  'co2.0000'
                  radonout=radon_00 ! 'radon.0000'
                  surfout=surf_00   ! 'current.0000'
                  qgout='qg_00'
               endif
            else                    ! 12 hour write
               if(ktime.eq.1200)then
                  co2out=co2_00     !  'co2.0000'
                  radonout=radon_00 ! 'radon.0000'
                  surfout=surf_00   ! 'current.0000'
                  qgout='qg_00'
               else
                  co2out=co2_12     ! 'co2.1200'
                  radonout=radon_12 ! 'radon.1200'
                  surfout=surf_12   ! 'current.1200'
                  qgout='qg_12'
               endif
            endif               ! (ktau.eq.nwrite)
            if ( myid == 0 ) then
               print *,'writing current soil & snow variables to ',
     &                  surfout
               open(unit=77,file=surfout,form='formatted',
     &              status='unknown')
               write (77,*) kdate,ktime,' ktau = ',ktau
            end if
            call writeglobvar(77, wb, fmt='(14f6.3)')
            call writeglobvar(77, tgg, fmt='(12f7.2)')
            call writeglobvar(77, tss, fmt='(12f7.2)')
            call writeglobvar(77, snowd, fmt='(12f7.1)')
            call writeglobvar(77, sicedep, fmt='(12f7.1)')
            if ( myid == 0 ) close (77)
            if(ico2.gt.0)then
               ico2x=max(1,ico2)
               if ( myid == 0 ) then
                  open(unit=77,file=co2out,form='formatted',
     &                 status='unknown')
                  write (77,*) kdate,ktime,' ktau = ',ktau
               end if
               call writeglobvar(77, tr(:,:,ico2x), fmt='(12f7.2)')
               if ( myid==0 ) close (77)
            endif               ! (ico2.gt.0)
            if(iradon.gt.0)then
               iradonx=max(1,iradon)
               if ( myid == 0 ) then
                  open(unit=77,file=radonout,form='formatted',
     &                 status='unknown')
                  write (77,*) kdate,ktime,' ktau = ',ktau
               end if
               call writeglobvar (77, tr(:,:,iradonx), fmt='(12f7.1)')
               if ( myid == 0 ) close (77)
            endif               ! (ico2.gt.0)
            if(nrungcm.eq.-2.or.nrungcm.eq.-5)then
               if ( myid == 0 ) then
                  print *,'writing special qgout file: ',qgout
                  open(unit=77,file=qgout,form='unformatted',
     &                 status='unknown')
               end if
               call writeglobvar(77, qg)
               if ( myid == 0 ) close (77)
         endif  ! (nrungcm.eq.-2.or.nrungcm.eq.-5)
        endif  ! (ktau.eq.nwrite/2.or.ktau.eq.nwrite)
      endif    ! (nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)


c---------------------------------------------------------------------------
      if(io_outt.eq.1)then  ! for netcdf
         ms_out = 0             ! Not set anywhere else ?????
         if(iout.eq.19)then
            if ( myid==0 ) print *,"restart write of data to cdf"
            call outcdf(rundate,nmi,-1,ms_out)
            call end_log(outfile_end)
            return              ! done with restart data
         endif                  !  (iout.eq.19)
         if ( myid==0 ) print *,'calling outcdf from outfile'
         call outcdf(rundate,nmi,1,ms_out)  

      endif ! (io_outt.eq.1)

      call end_log(outfile_end)
      return
      end
