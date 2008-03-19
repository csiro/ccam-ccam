!     integer, parameter :: mydoubl2 = selected_real_kind (10,100)
!     real (kind=mydoubl2) x(ifull),y(ifull),z(ifull)
      real*8 x(ifull),y(ifull),z(ifull)
      real wts(ifull)
      common /xyz_g/ x,y,z,wts ! rlat and rlong gone
