      subroutine outfile(iout,ik,jk,kk,psa,psm,rundate,nmi,nsnowout,
     .                   nwrite)
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
      common/histave/eg_ave(ifull),fg_ave(ifull),ga_ave(ifull),
     .    epot_ave(ifull),cbas_ave(ifull),ctop_ave(ifull),
     .    qscrn_ave(ifull),tmaxscr(ifull),tminscr(ifull),tscr_ave(ifull)
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
      ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,npanels

      if(ktau.eq.0.and.io_rest.ne.0)then ! switch for cctocc4
        print *,'in outfile opening surfile: ',surfile
c       write surface information file at beginning of run
        open(unit=77,file=surfile,form='unformatted',status='unknown')
        do iq=1,ifull
         aa(iq)=ivegt(iq)+.01
         bb(iq)=isoilm(iq)+.01
        enddo   ! iq loop
        write (77) aa  ! ivegt
        write (77) bb  ! isoilm
        write (77) alb ! just starting albedo 
        write (77) rsmin
        write (77) zolnd
        close (77)
        do iphase = 1,3         !  mrd code to run fast on NEC
         do iface=1,4
          if = nface6(iface,iphase)
          istart = ioff6(iface,iphase) ! 1 or 2
          do j=1,il
           istart = 3 - istart ! 1 or 2 alternately
           do i=istart,il,2
            iq=ind(i,j,if)
            cc(iq)=iphase+.5               ! 3-colouring of grid cells
            if(land(iq))cc(iq)=cc(iq)+3. ! bigger for land points
           end do
          end do
         end do
        end do
      endif  ! (ktau.eq.0.and.io_rest.ne.0)
      if(ktau.eq.nsnowout)then
       write (fileout,'(6hsnowd. ,i4)')nrun  !  i.e. snowd.xxxx
       open(unit=77,file=fileout,status='unknown',form='formatted')
       write(77,'(i3,i4,2f6.1,f5.2,f9.0,''  snowd from run'',i5)')
     .                        ifull,rlong0,rlat0,schmidt,ds,nrun
       write (formout,'(1h(,i3,4hf8.1,1h))')il  !  i.e. (<il>f8.1)
       write(77,formout) snowd
       close(77)
       print *,'snowout written for ktau =',ktau
      endif  ! (ktau.eq.nsnowout)

      ndt=dt
      io_outt=io_out
      if(iout.eq.19)io_outt=io_rest  !  choice for writing restart file
      print *,'ofile written for iout,kdate,ktime,mtimer: ',
     .                           iout,kdate,ktime,mtimer

      call mslp(pmsl,psl,zs,t)
c     reincorporate land mask into surface temperature
      do iq=1,ifull
       tssout(iq)=tss(iq)
       if(.not.land(iq))tssout(iq)=-tss(iq)
      enddo   ! iq loop
      print *,'tssout: ',tssout(idjd)

      if(io_outt.eq.3)then
      write(iout) kdate,ktime,ktau,ik,jk,kk,m,nsd,io_out,nx1,
     2 nps,mex,mup,nem,mtimer,nmi,ndt,npsav,rundate,nhor,nkuo,
     3 khdif,kwt,nx3,nx4,timer,timeg,ds,nvad,nqg
     4 ,nunst,nrun,nrunx,khor,ksc,kountr,ndiur,nspare,nhorps,ndum1
     5 ,ms,ntsur,nrad,kuocb,nvmix,ntsea,nonl,nextout,ilt,ntrac
     . ,dumknbd,rhkuo,dumdu,dumtanl,rlong0,rlat0,schmidt
      write(iout) sig
      write(iout) psl
c     following is actually pmsl
      write(iout) pmsl
      write(iout) zs
      write(iout) em
      write(iout) f
      write(iout) tssout
      write(iout)precip
      if(nqg.ge.8)then
        write(iout) ((tgg(iq,k),iq=1,ifull),k=1,ms)
        write(iout) (( wb(iq,k),iq=1,ifull),k=1,ms)
      else      ! old write, meant for nsib=0 or 1
c       write subsoil temps in this case, lower one(s) first
        write(iout) (tgg(iq,ms),iq=1,ifull)
        write(iout) (tgg(iq,2),iq=1,ifull)
        write(iout)(wb(iq,1),iq=1,ifull)
        write(iout)(wb(iq,ms),iq=1,ifull)
      endif    ! (nqg.ge.8)  .. else ..
      write(iout)alb    ! called albsav till 14/5/01
      if(ktau.eq.0.and.nwrite0.eq.1)then
        write(iout)davt  ! davt written instead of precc for ktau=0
        print *,'ktau,davt(idjd)',ktau,davt(idjd)
      else
        write(iout)precc
      endif
      if(nqg.ge.6)then
        write(iout)snowd
        if(ktau.eq.0.and.nwrite0.eq.1)then
          print *,' grid-lengths written instead of cloudlo for ktau=0'
          write(iout) (ds/em(iq),iq=1,ifull)
          print *,' wet-frac(3) written instead of cloudmi for ktau=0'
          isoil = isoilm(idjd)
          print *,'isoil,swilt,sfc,ssat,land ',
     .             isoil,swilt(isoil),sfc(isoil),ssat(isoil),land(idjd)
          do iq=1,ifull
c           if(land(iq))then
             isoil=isoilm(iq)
             aa(iq)=(wb(iq,3)-swilt(isoil))/(sfc(isoil)-swilt(isoil))
c           else
c             aa(iq)=0.
c           endif
          enddo   ! iq loop
          print *,'wb,wetfrac_out ',wb(idjd,3),aa(idjd)
          write(iout) aa  ! wet-frac(3) instead of cloudmi
          print *,' sigmf written instead of cloudhi for ktau=0'
          write(iout) (sigmf(iq),iq=1,ifull)
          write(iout)sicedep   ! keep this one
          print *,' isoilm written instead of tscrn for ktau=0'
          write(iout) (isoilm(iq)+.01,iq=1,ifull) ! instead of tscrn
          print *,' map colouring written instead of qgscrn for ktau=0'
          write(iout)cc
          print *,'cc# ',((cc(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,' zolnd written instead of u10 for ktau=0'
          write(iout) zolnd
        else
          write(iout)cloudlo
          write(iout)cloudmi
          write(iout)cloudhi
          write(iout)sicedep
          write(iout)tscrn
          write(iout)qgscrn
          write(iout)u10
         endif
      endif
      if(nqg.ge.9)then
        if(ktau.eq.0.and.nwrite0.eq.1)then
          print *,' ivegt written instead of eg for ktau=0'
!         N.B. subtract 31 to get sib values
          write(iout) (ivegt(iq)+.01,iq=1,ifull)
        else
          write(iout) eg
        endif
        write(iout) fg
        write(iout) taux             ! wind stress (x-direction)
        write(iout) tauy             ! wind stress (y-direction)
        write(iout) runoff           ! runoff
      endif
      if(nqg.ge.10)then
        write(iout) tmaxscr
        write(iout) tminscr
        write(iout) tscr_ave   ! ave now done in globpe
      endif
      if(nqg.ge.11)then
	write(iout) (( wbice(iq,k),iq=1,ifull),k=1,ms)
      endif
      if(nqg.ge.12)then
	write(iout) tggsn
	write(iout) smass
	write(iout) ssdn
	write(iout) ssdnn
	write(iout) osnowd
	write(iout) snage
	write(iout) isflag
      endif ! nqg = 12
      if(nqg.ge.13)then
        print *,'writing rtsave and was tscrn3hr'
        write(iout) rtsave
        write(iout) rtsave
      endif
      print *,'q, psl, pmsl, zs written for ktau=',ktau
      write(iout) ((t(iq,k),iq=1,ifull),k=1,kwt)
      write(iout) ((u(iq,k),iq=1,ifull),k=1,kwt)
      write(iout) ((v(iq,k),iq=1,ifull),k=1,kwt)
      write(iout) ((qg(iq,k),iq=1,ifull),k=1,kwt)
c     before 2/12/96 were multiplying qg by 1000. for unformatted output
      print *,'t, u, v written for ktau = ',ktau
      if(npsav.ne.0) write(iout) (psa(n),psm(n),n=1,npsav)
      if(nsd.eq.1)then    
        write(iout) ((sdot(iq,k),iq=1,ifull),k=2,kwt+1) ! 2 from 31/5/00
      endif
        if(ilt.gt.1)then
          do ntr=1,ntrac
           write(iout) ((tr(iq,k,ntr),iq=1,ifull),k=1,kwt)
          enddo
        endif
      endif    ! (io_outt.eq.3)
	 
	 
        if(ktau.gt.1.and.ktau.eq.abs(newsoilm))then
!         writes current wb,tgg (not read in during indata if -ve)
!         put current data (usually after 1st 12 hours) into 'smoist.dat'
          open(unit=77,file='smoist.dat',status='unknown'
     .                ,form='formatted')
          write(77,'(16f5.2)') wb     ! write current wb
          close(77)
          open(unit=77,file='tsoil.dat',status='unknown'
     .                ,form='formatted')
          write(77,'(11f7.2)') tgg   ! write current tgg
          close(77)
        endif                    ! (ktau.eq.abs(newsoilm))

      if(nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)then
        if(ktau.eq.nwrite/2.or.ktau.eq.nwrite)then
!        usually after first 24 hours, save soil variables for next run
         if(ktau.eq.nwrite)then     ! 24 hour write
          if(ktime.eq.1200)then
           co2out=co2_12      ! 'co2.1200'
           radonout=radon_12  ! 'radon.1200'
           surfout=surf_12    ! 'current.1200'
           qgout='qg_12'
          else
           co2out=co2_00      !  'co2.0000'
           radonout=radon_00  ! 'radon.0000'
           surfout=surf_00    ! 'current.0000'
           qgout='qg_00'
          endif
         else                       ! 12 hour write
          if(ktime.eq.1200)then
           co2out=co2_00      !  'co2.0000'
           radonout=radon_00  ! 'radon.0000'
           surfout=surf_00    ! 'current.0000'
           qgout='qg_00'
          else
           co2out=co2_12      ! 'co2.1200'
           radonout=radon_12  ! 'radon.1200'
           surfout=surf_12    ! 'current.1200'
           qgout='qg_12'
          endif
         endif  ! (ktau.eq.nwrite)
         print *,'writing current soil & snow variables to ',surfout
         open(unit=77,file=surfout,form='formatted',status='unknown')
         write (77,*) kdate,ktime,' ktau = ',ktau
         write (77,'(14f6.3)') wb
         write (77,'(12f7.2)') tgg
         write (77,'(12f7.2)') tss
         write (77,'(12f7.1)') snowd
         write (77,'(12f7.1)') sicedep
         close (77)
         if(ico2.gt.0)then
           open(unit=77,file=co2out,form='formatted',status='unknown')
           write (77,*) kdate,ktime,' ktau = ',ktau
           write (77,'(12f7.2)')((tr(iq,k,ico2),iq=1,ilt*jlt),k=1,klt)
           close (77)
         endif  ! (ico2.gt.0)
         if(iradon.gt.0)then
           open(unit=77,file=radonout,form='formatted',status='unknown')
           write (77,*) kdate,ktime,' ktau = ',ktau
           write (77,'(12f7.1)')((tr(iq,k,iradon),iq=1,ilt*jlt),k=1,klt)
           close (77)
         endif  ! (ico2.gt.0)
         if(nrungcm.eq.-2.or.nrungcm.eq.-5)then
	    print *,'writing special qgout file: ',qgout
           open(unit=77,file=qgout,form='unformatted',status='unknown')
           write (77) qg
           close (77)
         endif  ! (nrungcm.eq.-2.or.nrungcm.eq.-5)
        endif  ! (ktau.eq.nwrite/2.or.ktau.eq.nwrite)
      endif    ! (nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)

c---------------------------------------------------------------------------
      if(io_outt.eq.1)then  ! for netcdf
        if(iout.eq.19)then
         print *,"restart write of data to cdf"
         call outcdf(rundate,nmi,-1,ms_out)
         call ncclos(idnc,ier)
         write(6,*) "call ncclos(idnc,ier) ",idnc,ier
         return ! done with restart data
        endif    !  (iout.eq.19)
        print *,'calling outcdf from outfile'
        call outcdf(rundate,nmi,1,ms_out)  

        if(ktau.eq.ntau)then
           write(6,*) "calling ncclos(idnc,ier) ",idnc,ier
           call ncclos(idnc,ier)
        endif
      endif ! (io_outt.eq.1)

      return
      end
