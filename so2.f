c this code needs major changes to get i,j out of the common block
c     so2.f bundles together so2accumulate, so2dry, so2inisrc,
c                            so2sflux and so2vmix
      subroutine so2accumulate

* clobber the so2 tracer & accumulate it

      include 'const_phys.h'
      include 'newmpar.h'
      include 'tracers.h'

      do i=1,ijk
        if( tr(i,1,max(1,iso2)).lt.0. ) then
           tr(i,1,max(1,iso2)) = 0.
        else
           sumtr(i,1) = sumtr(i,1) + tr(i,1,max(1,iso2))
        end if
      end do
      nsumtr = nsumtr + 1

c      write(*,*)' so2accumulate: clobbered tr(so2), nsumtr =', nsumtr

      return
      end
      subroutine so2dry( kappa_d )

* ---
c ---  this routine computes vertical tracer (so2) transport with rain
c --- coupled with the modified routine conjob
* ---
* ............ updates dtr
*
      include 'const_phys.h'
      include 'newmpar.h'
      include 'parm.h'
      include 'sigs.h'          ! dsig
      include 'tracers.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)

      real    l, beta_0, f_0, f_k_max, kappa_w, kappa_d
      data    l/5.0e-4/
      data    f_k_max/0./, f_0/0.8/

      save f_k, f_k_max, f_0, l

      common/so2block/ i, j, k, beta_0, psong, drain, rainb, dtr
       iq=i+(j-1)*il

* note:
*	kappa_w --- deposition rate per time step

      idry = 1
      go to 10

      entry so2drywet( kappa_d, rk_wet, frac, dtrdry, dtrwt )

      idry = 0
      
10    if( tr(iq,k,max(1,iso2)).le.0. ) then
        tr(iq,k,max(1,iso2)) = 0.
        dtr = (-trsrc(iq,k)/kappa_d)*(exp(-kappa_d)-1.)
      else
        dtr = (tr(iq,k,max(1,iso2))-trsrc(iq,k)/kappa_d)*
     .                       (exp(-kappa_d)-1.)
      end if

      if(idry.eq.1) return           ! only dry deposition required

      dtrdry = dtr
      go to 20
    
*.................................................................................

      entry so2wet( kappa_d )
      idry = -1

20    f_k = f_k_max

      if( drain.gt.0. ) then         ! drain > 0   indicates cloud water

        rainonl = drain/l
        f_k     = f_0/(1.+beta_0/rainonl)
        f_k_max = max(f_k,f_k_max)

      else				! below cloud (evaporation)

        rainonl =  cem*rainb

      end if

      rainonl = rainonl*psong
      temp = rainonl/f_k
      kappa_w = kappa_d + temp

      if ( tr(iq,k,max(1,iso2)).le.0. ) tr(iq,k,max(1,iso2)) = 0.

      if( idry.eq.0 ) then
        frac = f_k
        rk_wet = temp
        dtrwt = f_k*(tr(iq,k,max(1,iso2))-trsrc(iq,k)/kappa_w)*
     .                       (exp(-kappa_w)-1.)
        dtrdry = (1.-f_k)*dtr
        dtr = dtrwt + dtrdry
      else
        dtr = f_k*(tr(iq,k,max(1,iso2))-trsrc(iq,k)/kappa_w)*
     .                       (exp(-kappa_w)-1.) +
     .      (1.-f_k)*dtr
      end if

      return
      end
      subroutine so2inisrc

* ----
c ----   the purpose of this routine is to initialise the tracer (so2) source block
c ----  from a list of point sources
* ----
      include 'const_phys.h'
      include 'newmpar.h'
      include 'arrays.h'
      include 'nsibd.h'
      include 'tracers.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)

      logical firstcall
      save    firstcall
      data    firstcall/.true./

*   ---- level 1 includes surface emissions as well as point sources

      if( firstcall ) then
        oneab = 1. - alpha - beta
        do k=1,nso2lev
          do ix=iso2lev(k-1)+1,iso2lev(k)
            so2em(iso2emindex(ix)) = so2em(iso2emindex(ix)) *
     .          oneab * (-so2fact(k))
          end do !      do ix=iso2lev(k-1)+1,iso2lev(k)

          oneab = 1. - beta	! at levels above one there is no immediate fallout at present

        end do !     do k=2,nso2lev
        firstcall = .not.firstcall
      end if

      do k=1,nso2lev
        do ix=iso2lev(k-1)+1,iso2lev(k)
          i = iso2em(1,iso2emindex(ix))
          j = iso2em(2,iso2emindex(ix))
          iq=i+(j-1)*il
          trsrc(iq,k) = trsrc(iq,k) + so2em(iso2emindex(ix))/ps(iq)
        end do !      do ix=iso2lev(k-1)+1,iso2lev(k)
      end do !     do k=2,nso2lev

      return
      end

      subroutine so2sflux
      include 'newmpar.h'
      include 'const_phys.h'
      include 'arrays.h'
      include 'map.h'
      include 'nsibd.h'
      include 'tracers.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)
      logical firstcall
      save    firstcall
c     data    firstcall/.true./,j/1/,k/1/
      data    firstcall/.true./

* so2background --- emissions per unit time step

      if( firstcall ) then
        so2background = so2background * (-so2fact(1))
        firstcall = .not.firstcall
      end if

      do iq=1,ifull

        if( land(iq) ) then
          trsrc(iq,1)   = so2background / ps(iq)
        else
          trsrc(iq,1)   = 0.
        endif

      end do ! iq=1,ifull

      do k=2,kl
      do iq=1,ifull
        trsrc(iq,k)   = 0.
      end do ! iq
      end do ! k


      return
      end
      subroutine so2vmix( updtr )

      include 'const_phys.h'
      include 'newmpar.h'
      include 'arrays.h'
      include 'map.h'
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'		! ktau
c      include 'rain.h'		! dsig
      include 'sigs.h'		! dsig
      include 'tracers.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)

* land use indices
      parameter ( ibrlfevrgrtr =  1,
     .            ibrlfdecidtr =  2,
     .            ibrlfneedltr =  3,
     .            ineedevrgrtr =  4,
     .            ineeddecidtr =  5,
     .            ibrlfgrdcvtr =  6,
     .            igrdcvonly   =  7,
     .            ibrlfgrdcvsh =  8,
     .            ibrlfbrsoish =  9,
     .            idwarfshrbtr = 10,
     .            ibaresoil    = 11,
     .            iwnwhbrlfdtr = 12,
     .            isea         = 13)


c     real vdrydepso2(13)
c     data vdrydepso2/6*0.7e-2,0.3e-2,0.5e-2,2*0.3e-2,0.2e-2,0.5e-2,
c    .                  1.0e-2/ ! so2 dry deposition velocities in m/s
* vegetation height map
      parameter ( ibaregrnd = 1,
     .            ilowveg   = 2,
     .            imedveg   = 3,
     .            ihighveg  = 4,
     .            iocean    = 5)

      real vdrydepso2(5)
      data vdrydepso2/0.2e-2, 0.3e-2, 0.5e-2, 0.7e-02, 1.0e-2/ ! so2 dry deposition velocities in m/s

      integer imap(44)
      data imap/7*ihighveg,5*imedveg,15*ilowveg,imedveg,3*ibaregrnd,
     .          6*ihighveg,ilowveg,imedveg,2*ilowveg,ibaregrnd,imedveg,
     .          iocean/
c      integer imap(13)
c      data imap/6*ihighveg,ilowveg,imedveg,2*ilowveg,ibaregrnd,imedveg,
c     .            iocean/

      integer i,j,k
      real    updtr(ifull,kl)
      real    kappa_d
      real    c, emis, r_t, r_tdt

      data    c/5.2/, emis/0.1/

      logical firstcall
      data    firstcall/.true./

      save r_t, r_tdt, firstcall

      common/so2block/ i, j, k, beta_0, psong, drain, rainb, dtr

* note:
* 	trsrc   --- source emission per time step
*	r_tdt   --- so2 <---> so4 conversion rate per time step
*	kappa_d --- dry deposition rate per time step

*  --- initialise conversion factors and coefficients
      if( firstcall ) then
        alpha = 0.
        beta  = 0.
        r_t   = 0.
        r_tdt = r_t*dt
        beta_0= dt*8.0e-5
        cem   = c*emis
        firstcall = .not.firstcall
      end if
* --- initialise trsrc's with so2 emission sources

      call so2inisrc

      do j=1,jl
        do i=1,il
          iq=i+(j-1)*il
          psong     = ps(iq)/grav

                                     
* ---- define dry deposition values
          if( land(iq) ) then
            vdryso2 =     tsigmf(iq) *vdrydepso2(imap(ivegt(iq))) +
     .                (1.-tsigmf(iq))*vdrydepso2(       ibaregrnd)
          else
            vdryso2 = vdrydepso2( imap(isea) )
          end if

* ---               rho dz = - ps dsigma / g
* ---               p = rho r t
* ---               condrag = - g dt / r dsig

          vdryso2 = - vdryso2 * condrag / t(iq,1)
         

* --- 	do so2 for all the layers above bottom one
          sumdtr  = 0.

          raint = 0.
          do k=kl,2,-1
*           rain(iq,k) = 0.  !!!!!!!!!!!!!!!!!!! temporary
            rainb = rain(iq,k)
            drain = rainb - raint
 
* --- increment the source by evaporated rain fraction
            if( drain.lt.0. ) then
              evfrac  = sumdtr*drain/raint
              trsrc(iq,k) = trsrc(iq,k) + evfrac/(-dsig(k))
            end if
            dtr = trsrc(iq,k)

            if( rainb.ne.0. ) then   ! drain > 0 indicates cloud water, < 0 rain evaporation
              call so2wet( r_tdt )
              sumdtr = sumdtr + (-dsig(k))*(dtr-trsrc(iq,k))
            else
              if(tr(iq,k,max(1,iso2)).lt.0.) tr(iq,k,max(1,iso2)) = 0.
              sumdtr = 0.
            end if

            updtr(iq,k) = tr(iq,k,max(1,iso2)) + dtr
            raint = rainb

          end do

* --- 	do so2 for the bottom layer
          k=1 
          kappa_d = r_tdt + vdryso2
*         rain(iq,k) = 0.  !!!!!!!!!!!!!!!!!!! temporary
          rainb   = rain(iq,k)
          drain   = rainb - raint

* --- increment the source by evaporated rain fraction
          if( drain.lt.0. ) then
            evfrac  = sumdtr*drain/raint
            trsrc(iq,k) = trsrc(iq,k) + evfrac/(-dsig(k))
          end if

          if( rainb.ne.0. ) then     ! drain > 0 indicates cloud water, < 0 rain evaporation
            call so2drywet( kappa_d, rk_wet, frac, dtrdry, dtrwt )

* --- accumulate the wet deposition into respective arrays
*     after converting to air mass equivalent 
* ---   rain in kg m^-2 (because mass mixing ratio is used)
            src = frac*trsrc(iq,k)
            tmp = psong*(src - dtrwt) * (-dsig(k))/(kappa_d + rk_wet)
            sumtr(iq,iwetdep) = sumtr(iq,iwetdep) +
     .                         (tmp*rk_wet - sumdtr*psong)
            sumtr(iq,irain)   = sumtr(iq,irain)   + rain(iq,k)*psong
            src = trsrc(iq,k) - src
          else
            call so2dry( kappa_d )
            src    = trsrc(iq,k)
            dtrdry = dtr
            tmp    = 0.
          end if

          tmp  = tmp + psong*(src - dtrdry)*(-dsig(k))/kappa_d

* --- accumulate the dry deposition into respective arrays
*     after converting to air mass equivalent 
          sumtr(iq,idrydep)=sumtr(iq,idrydep) + tmp*vdryso2
          updtr(iq,k) = tr(iq,k,max(1,iso2)) + dtr
        end do ! i=1,il
      end do   ! j=1,jl

* --- add to the dry deposition the tracer that drops out immediately after emission
        do ix=iso2lev(k-1)+1,iso2lev(k)
          i = iso2em(1,iso2emindex(ix))
          j = iso2em(2,iso2emindex(ix))
          sumtr(iq,idrydep) = sumtr(iq,idrydep) +
     .                                     alpha*so2em(iso2emindex(ix))
        end do !      do ix=iso2lev(k-1)+1,iso2lev(k)

      return
      end
