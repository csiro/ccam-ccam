!     rml 21/02/06 deleted variables that were no longer used

      integer iradon,ico2,ngas,nllp,ntrac,ntracmax,                       &
     &        npwr,ilt,jlt,klt,ngasmax
      real tr,traver,radonem
!     parameter(iradon=0,ico2=1,ngas=9)  ! rlw setting
      parameter(iradon=0,ico2=0,ngas=0)
      parameter(nllp=0)              ! set this or next lines
!     parameter(nllp=3)

      parameter(ntrac=ngas+nllp)
      parameter(ntracmax=max(ntrac,1))      ! ntracmax >= 1
      parameter(ngasmax=max(ngas,1))  ! ngasmax >= 1

      parameter(npwr=min(ntrac,1))        ! i.e. 0 or 1 for no_gas/gas
      parameter(ilt=il**npwr,jlt=jl**npwr,klt=kl**npwr) ! gives one of next two   

      common/tracer/tr(ilt*jlt+npwr*iextra,klt,ntracmax),                 &
     &              traver(ilt*jlt,klt,ntrac)
!      suspect alpha and beta can go but need to confirm
!      real alpha,beta
!      common/c_tracer/ alpha,beta  

      common/emission/radonem(ilt*jlt)

      real gasmin(ngasmax)           ! used by mass fixer in adjust5
      data gasmin/ngasmax*-1000./    ! modify this line for specific gases        
!     data gasmin/-1000.,-1000.,-1000.,-1000.,-1000.,0.,0.,0.,-1000./ ! rlw
