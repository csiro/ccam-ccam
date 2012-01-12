! from July 2006, use split scheme
!  version for globpex with tendencies
!  This is a simplified physics routine implementing the temperature
!  relaxation and wind drag for Held-Suarez dynamical core experiments. 
!  The drag coefficient depends on height and the temperature relaxation
!  time depends on latitude and height.

!  The equilibrium temperature depends only on latitude and height, but on 
!  the conformal-cubic grid this requires a full 3D array. To save this
!  it's recalculated for each point.

      subroutine hs_phys
      use arrays_m
      use latlong_m
      use nlin_m
      use sigs_m
      include 'newmpar.h'
      include 'parm.h'

      parameter(stefan=5.67e-8)

!     All coefficients are in units of inverse days
      real invday
      parameter (invday=1./86400.)
      real kf, ks, ka
      parameter(kf = 1. * invday)
!      parameter(kf = 2. * invday)
      parameter(ks = 0.25 * invday)
      parameter(ka = 0.025 * invday)

      real sig_b
      parameter(sig_b = 0.7)  ! Drag applied below this level
      real delty
      parameter(delty = 60.)    ! Pole to equator variation in equil temperature
      real deltheta
      parameter(deltheta = 10.) ! Vertical variation
      real kappa
      parameter(kappa = 2./7.)

      real kv, kt, teq

      do k=1,kl 
       do iq=1,ifull
           kt = ka +
     &          (ks-ka)*max(0., (sig(k)-sig_b)/(1.-sig_b)) *
     &          cos(rlatt(iq))**4
           teq = max ( 200.,
     &               (315. - delty*sin(rlatt(iq))**2 -
     &                deltheta*log(sig(k))*cos(rlatt(iq))**2)
     &                *sig(k)**kappa )
!          following for tendencies based on implicit a
!          tn(iq,k) =tn(iq,k)-kt*(t(iq,k)-teq)/(1.+dt*kt)  ! implicit form a
           t(iq,k) = (t(iq,k)*(1.-.5*dt*kt)+dt*kt*teq)/
     .                                  (1.+.5*dt*kt)  ! implicit form b
       end do
      end do

!     Winds have a height dependent drag
      do k=1,kl
       kv = kf * max(0., (sig(k)-sig_b)/(1.-sig_b))
       do iq=1,ifull
!          following for tendencies based on implicit a
c          un(iq,k) = un(iq,k) -kv*u(iq,k)/(1.+ dt * kv)
c          vn(iq,k) = vn(iq,k) -kv*v(iq,k)/(1.+ dt * kv)
           u(iq,k) = u(iq,k)*(1.-.5*dt*kv)/(1.+.5*dt*kv) ! im form b
           v(iq,k) = v(iq,k)*(1.-.5*dt*kv)/(1.+.5*dt*kv) ! im form b
       end do
      end do

      return
      end
