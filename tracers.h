      integer iradon,ico2,iso2,iso4,ich4,io2,ngas,nllp,ntrac,ntracmax,
     .        npwr,ilt,jlt,klt,indso2,ilso2,jlso2,ngasmax,
     .        nsumtr,k_t,ico2em,idrydep,iwetdep,irain,nso2sour,nso2slev,
     .        nso2lev,iso2em,jso2em,iso2emindex,iso2lev
      real tr,sumtr,cem,alpha,beta,radonem,rain,
     .     so2background,co2fact,condrag,so2em,so2fact
      parameter(iradon=0,ico2=0,iso2=0,iso4=0,ich4=0,io2=0,ngas=0)
!     parameter(iradon=1,ico2=2,iso2=0,iso4=0,ich4=0,io2=0,ngas=2)
      parameter(nllp=0)              ! set this or next lines
!     parameter(nllp=3)

      parameter(ntrac=ngas+nllp)

      parameter(ntracmax=max(ntrac,1))      ! ntracmax >= 1
!     parameter(ntracmax=ntrac)             ! set this or next line
!     parameter(ntracmax=1)                 ! no gases & nllp=0

      parameter(npwr=min(ntrac,1))        ! i.e. 0 or 1
      parameter(ilt=il**npwr,jlt=jl**npwr,klt=kl**npwr) ! gives one of next two   
!     parameter(ilt= il,jlt= jl,klt=kl)      ! set this or next line
!     parameter(ilt=  1,jlt=  1,klt= 1)      ! no gases & nllp=0

      parameter(indso2=0,ilso2=1,jlso2=1)

      common/tracer/tr(ilt*jlt,klt,ntracmax),sumtr(ilso2*jlso2,-2:klt)
      common/c_tracer/ nsumtr,cem,k_t,alpha,beta  
      common/emission/radonem(ilt*jlt),ico2em(ilt*jlt)

      common/raindata/rain(ilso2*jlso2,kl)            
      parameter ( idrydep = -2,iwetdep= -1,irain= 0)

      parameter(nso2sour=   1,nso2slev=3)
      common/cnsib/nso2lev,so2background,co2fact,     
     -  condrag,iso2em(2,nso2sour),       
     -  jso2em(nso2sour),so2em(nso2sour),           
     -  iso2emindex(nso2sour),so2fact(nso2slev),iso2lev(0:nso2slev)   
     
      parameter(ngasmax=max(ngas,1))  ! ngasmax >= 1
!     parameter(ngasmax=ngas)         ! ngasmax >= 1
!     parameter(ngasmax=1)            ! ngasmax >= 1
      real gasmin(ngasmax)            ! used by mass fixer in adjust5
      data gasmin/ngasmax*0./         ! modify this line for specific gases          
