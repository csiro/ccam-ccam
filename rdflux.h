c     common block rdflux contains flux quantities computed by the
c     radiation code, used for diagnostic purposes: 

      real flx1e1(imax)     ! Flux at top for 0-160,1200-2200 cm-1 range
      real gxcts(imax)      ! Flux at top for 160-1200 cm-1 range 
      real fctsg(imax,nbly) ! CTS flux at ground. used to obtain gxcts
                            !    by bands.
      real flx1e1clr(imax)
      real gxctsclr(imax)

      common / rdflux / flx1e1, gxcts, fctsg
      common / clrflx / flx1e1clr, gxctsclr
c 
