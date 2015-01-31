subroutine ints(ntr,s,intsch,nface,xg,yg,nfield)

!     parameter (mhint):   0 for simple; 2 for Bessel (was called mbess)
!     jlm finds Bessel (mhint=2) gives bigger overshoots near discontinuities
!     this one includes Bermejo & Staniforth option
!     parameter (mh_bs): B&S on/off depending on value of nfield
!     can put sx into work array (only 2d)
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation
!     nfield: 1 (psl), 2 (u, v), 3 (T), 4 (gases)

use cc_mpi
use indices_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'    ! has mh_bs

integer, parameter :: ntest=0
integer, intent(in) :: ntr
integer, intent(in) :: intsch, nfield
integer idel, iq, jdel, nn
integer i, j, k, n, ip, jp
integer ii
integer, dimension(ifull,kl), intent(in) :: nface
real xxg, yyg
real, dimension(ifull,kl), intent(in) :: xg,yg ! now passed through call
real, dimension(ntr,-1:ipan+2,-1:jpan+2,1:npan,kl) :: sx
real, dimension(ifull+iextra,kl,ntr), intent(inout) :: s
real, dimension(ntr) :: c1, c2, c3, c4
real, dimension(ntr) :: a3, a4, sss
real, dimension(ntr) :: cmax, cmin
real, dimension(ntr,4) :: r

call START_LOG(ints_begin)

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

  do nn=1,ntr
    sx(nn,1:ipan,1:jpan,1:npan,1:kl) = reshape(s(1:ipan*jpan*npan,1:kl,nn), (/ipan,jpan,npan,kl/))
            
    ! this is intsb           EW interps done first
    ! first extend s arrays into sx - this one -1:il+2 & -1:il+2
    do k=1,kl
      do n=1,npan
        do j=1,jpan
          sx(nn,0,j,n,k)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),k,nn)
          sx(nn,-1,j,n,k)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),k,nn)
          sx(nn,ipan+1,j,n,k) = s(ie(j*ipan+(n-1)*ipan*jpan),k,nn)
          sx(nn,ipan+2,j,n,k) = s(iee(j*ipan+(n-1)*ipan*jpan),k,nn)
        enddo            ! j loop
        do i=1,ipan
          sx(nn,i,0,n,k)      = s(is(i+(n-1)*ipan*jpan),k,nn)
          sx(nn,i,-1,n,k)     = s(iss(i+(n-1)*ipan*jpan),k,nn)
          sx(nn,i,jpan+1,n,k) = s(in(i-ipan+n*ipan*jpan),k,nn)
          sx(nn,i,jpan+2,n,k) = s(inn(i-ipan+n*ipan*jpan),k,nn)
        enddo            ! i loop
        ! for ew interpolation, sometimes need (different from ns):
        ! (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
        ! (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

        sx(nn,-1,0,n,k)          = s(lwws(n),k,nn)
        sx(nn,0,0,n,k)           = s(iws(1+(n-1)*ipan*jpan),k,nn)
        sx(nn,0,-1,n,k)          = s(lwss(n),k,nn)
        sx(nn,ipan+1,0,n,k)      = s(ies(ipan+(n-1)*ipan*jpan),k,nn)
        sx(nn,ipan+2,0,n,k)      = s(lees(n),k,nn)
        sx(nn,ipan+1,-1,n,k)     = s(less(n),k,nn)
        sx(nn,-1,jpan+1,n,k)     = s(lwwn(n),k,nn)
        sx(nn,0,jpan+2,n,k)      = s(lwnn(n),k,nn)
        sx(nn,ipan+2,jpan+1,n,k) = s(leen(n),k,nn)
        sx(nn,ipan+1,jpan+2,n,k) = s(lenn(n),k,nn)
        sx(nn,0,jpan+1,n,k)      = s(iwn(1-ipan+n*ipan*jpan),k,nn)
        sx(nn,ipan+1,jpan+1,n,k) = s(ien(n*ipan*jpan),k,nn)
      enddo               ! n loop
    end do                ! k loop
  end do                  ! nn loop

! Loop over points that need to be calculated for other processes
  if(nfield<mh_bs)then
    if(mhint==2)then ! Bessel interp
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff
             
          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = r(nn,2) + 0.5*yyg*(r(nn,3)-r(nn,1)+yyg*(a3(nn)+yyg*a4(nn)))
          end do
        enddo         ! iq loop
      end do          ! ii loop
    else
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff
              
          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)
                    
          r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*r(nn,2)-yyg*r(nn,1)/3.)      &
                                         -yyg*(1.+yyg)*r(nn,4)/3.)+yyg*(1.+yyg)*(2.-yyg)*r(nn,3))/2.
          end do
        enddo         ! iq loop
      end do          ! ii loop
    endif             !  (mhint==2)
  else                ! (nfield<mh_bs)
    if(mhint==2)then ! Bessel interp
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff

          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          sss = r(:,2)+.5*yyg*(r(:,3)-r(:,1)+yyg*(a3+yyg*a4))
          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = min(max(cmin(nn),sss(nn)),cmax(nn)) ! Bermejo & Staniforth
          end do
        end do        ! iq loop
      end do          ! ii loop
    else
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff

          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
          r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          sss = ((1.-yyg)*((2.-yyg)*((1.+yyg)*r(:,2)-yyg*r(:,1)/3.)-yyg*(1.+yyg)*r(:,4)/3.)+yyg*(1.+yyg)*(2.-yyg)*r(:,3))/2.
          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = min(max(cmin(nn),sss(nn)),cmax(nn)) ! Bermejo & Staniforth
          end do
        end do        ! iq loop
      end do          ! ii loop
    endif             !  (mhint==2)
  endif               ! (nfield<mh_bs)  .. else ..

  call intssync_send(ntr)

  if(nfield<mh_bs)then
    if(mhint==2)then ! Bessel interp
      do k=1,kl
        do iq=1,ifull    ! non Berm-Stan option
          ! Convert face index from 0:npanels to array indices
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
          jdel=int(yg(iq,k))
          yyg=yg(iq,k)-jdel
          ! Now make them proper indices in this processor's region
          idel = idel - ioff
          jdel = jdel - joff
          n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          s(iq,k,:) = r(:,2)+.5*yyg*(r(:,3)-r(:,1)+yyg*(a3+yyg*a4))
        enddo         ! iq loop
      enddo           ! k loop
    else
      do k=1,kl
        do iq=1,ifull    ! non Berm-Stan option
          ! Convert face index from 0:npanels to array indices
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
          jdel=int(yg(iq,k))
          yyg=yg(iq,k)-jdel
          ! Now make them proper indices in this processor's region
          idel = idel - ioff
          jdel = jdel - joff
          n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)

          r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          s(iq,k,:) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*r(:,2)-yyg*r(:,1)/3.)-yyg*(1.+yyg)*r(:,4)/3.)+yyg*(1.+yyg)*(2.-yyg)*r(:,3))/2.
        enddo         ! iq loop
      enddo           ! k loop
    endif             !  (mhint==2)
  else                ! (nfield<mh_bs)
    if(mhint==2)then ! Bessel interp
      do k=1,kl
        do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
           jdel=int(yg(iq,k))
           yyg=yg(iq,k)-jdel
           ! Now make them proper indices in this processor's region
           idel = idel - ioff
           jdel = jdel - joff
           n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          sss = r(:,2)+.5*yyg*(r(:,3)-r(:,1)+yyg*(a3+yyg*a4))
          s(iq,k,:) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
        enddo           ! iq loop
      enddo             ! k loop
    else
      do k=1,kl
        do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
          jdel=int(yg(iq,k))
          yyg=yg(iq,k)-jdel
          ! Now make them proper indices in this processor's region
          idel = idel - ioff
          jdel = jdel - joff
          n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
          c2 = sx(:,idel  ,jdel,n,k)
          c3 = sx(:,idel+1,jdel,n,k)
          c4 = sx(:,idel+2,jdel,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
          c1 = sx(:,idel-1,jdel+1,n,k)
          c2 = sx(:,idel  ,jdel+1,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+2,jdel+1,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
          r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
!         r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!              -x*(1+x)*c4/3}
!              +x*(1+x)*(2-x)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel  ,jdel+nn-2,n,k)
            c3 = sx(:,idel+1,jdel+nn-2,n,k)
            r(:,nn) = (1.-xxg)*c2 +xxg*c3
          enddo         ! nn loop
          sss = ((1.-yyg)*((2.-yyg)*((1.+yyg)*r(:,2)-yyg*r(:,1)/3.)-yyg*(1.+yyg)*r(:,4)/3.)+yyg*(1.+yyg)*(2.-yyg)*r(:,3))/2.
          s(iq,k,:) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
        enddo         ! iq loop
      enddo           ! k loop
    endif             !  (mhint==2)
  endif               ! (nfield<mh_bs)  .. else ..
            
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
  do nn=1,ntr
    sx(nn,1:ipan,1:jpan,1:npan,1:kl) = reshape(s(1:ipan*jpan*npan,1:kl,nn), (/ipan,jpan,npan,kl/))
    do k=1,kl
      do n=1,npan
        do j=1,jpan
          sx(nn,0,j,n,k)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),k,nn)
          sx(nn,-1,j,n,k)     = s(iww(1+(j-1)*ipan+(n-1)*ipan*jpan),k,nn)
          sx(nn,ipan+1,j,n,k) = s(ie(j*ipan+(n-1)*ipan*jpan),k,nn)
          sx(nn,ipan+2,j,n,k) = s(iee(j*ipan+(n-1)*ipan*jpan),k,nn)
        enddo            ! j loop
        do i=1,ipan
          sx(nn,i,0,n,k)      = s(is(i+(n-1)*ipan*jpan),k,nn)
          sx(nn,i,-1,n,k)     = s(iss(i+(n-1)*ipan*jpan),k,nn)
          sx(nn,i,jpan+1,n,k) = s(in(i-ipan+n*ipan*jpan),k,nn)
          sx(nn,i,jpan+2,n,k) = s(inn(i-ipan+n*ipan*jpan),k,nn)
        enddo            ! i loop
!       for ns interpolation, sometimes need (different from ew):
!           (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

        sx(nn,-1,0,n,k)          = s(lsww(n),k,nn)
        sx(nn,0,0,n,k)           = s(isw(1+(n-1)*ipan*jpan),k,nn)
        sx(nn,0,-1,n,k)          = s(lssw(n),k,nn)
        sx(nn,ipan+2,0,n,k)      = s(lsee(n),k,nn)
        sx(nn,ipan+1,-1,n,k)     = s(lsse(n),k,nn)
        sx(nn,-1,jpan+1,n,k)     = s(lnww(n),k,nn)
        sx(nn,0,jpan+1,n,k)      = s(inw(1-ipan+n*ipan*jpan),k,nn)
        sx(nn,0,jpan+2,n,k)      = s(lnnw(n),k,nn)
        sx(nn,ipan+2,jpan+1,n,k) = s(lnee(n),k,nn)
        sx(nn,ipan+1,jpan+2,n,k) = s(lnne(n),k,nn)
        sx(nn,ipan+1,0,n,k)      = s(ise(ipan+(n-1)*ipan*jpan),k,nn)
        sx(nn,ipan+1,jpan+1,n,k) = s(ine(n*ipan*jpan),k,nn)
      enddo               ! n loop
    end do                ! k loop
  end do                  ! nn loop

  ! For other processes
  if(nfield<mh_bs)then
    if(mhint==2)then ! Bessel interp
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff
          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!           -y*(1+y)*c4/3}
!           +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop

          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = r(nn,2)+0.5*xxg*(r(nn,3)-r(nn,1) +xxg*(a3(nn)+xxg*a4(nn)))
          end do
        end do           ! iq loop
      end do             ! ii
    else
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff
          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)

          r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)
                
          r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!           -y*(1+y)*c4/3}
!           +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop

          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*r(nn,2)-xxg*r(nn,1)/3.)      & 
                                         -xxg*(1.+xxg)*r(nn,4)/3.)+xxg*(1.+xxg)*(2.-xxg)*r(nn,3))/2.
          end do
        end do           ! iq loop
      end do             ! ii
    endif                !  (mhint==2)
  else                   ! (nfield<mh_bs)
    if(mhint==2)then ! Bessel interp
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff
          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!            -y*(1+y)*c4/3}
!            +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop

          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          sss = r(:,2)+.5*xxg*(r(:,3)-r(:,1)+xxg*(a3+xxg*a4))

          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = min(max(cmin(nn),sss(nn)),cmax(nn)) ! Bermejo & Staniforth
          end do
        enddo         ! iq loop
      end do          ! ii

    else
      do ii=neighnum,1,-1
        do iq=1,drlen(ii)
          !  Convert face index from 0:npanels to array indices
          ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
          jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
          n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
          !  Need global face index in fproc call
          idel = int(dpoints(ii)%a(2,iq))
          xxg = dpoints(ii)%a(2,iq) - idel
          jdel = int(dpoints(ii)%a(3,iq))
          yyg = dpoints(ii)%a(3,iq) - jdel
          k = nint(dpoints(ii)%a(4,iq))
          idel = idel - ioff
          jdel = jdel - joff
          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth

          r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!            -y*(1+y)*c4/3}
!            +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop

          sss = ((1.-xxg)*((2.-xxg)*((1.+xxg)*r(:,2)-xxg*r(:,1)/3.)-xxg*(1.+xxg)*r(:,4)/3.)+xxg*(1.+xxg)*(2.-xxg)*r(:,3))/2.

          do nn=1,ntr
            sextra(ii)%a(nn+(iq-1)*ntr) = min(max(cmin(nn),sss(nn)),cmax(nn)) ! Bermejo & Staniforth
          end do
        enddo         ! iq loop
      end do          ! ii
    endif             !  (mhint==2)
  endif               ! (nfield<mh_bs)  .. else ..

  call intssync_send(ntr)

  if(nfield<mh_bs)then
    if(mhint==2)then ! Bessel interp
      do k=1,kl
        do iq=1,ifull    ! non Berm-Stan option
          ! Convert face index from 0:npanels to array indices
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
          jdel=int(yg(iq,k))
          yyg=yg(iq,k)-jdel
          ! Now make them proper indices in this processor's region
          idel = idel - ioff
          jdel = jdel - joff
          n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!              -y*(1+y)*c4/3}
!              +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop

          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          s(iq,k,:) = r(:,2)+.5*xxg*(r(:,3)-r(:,1)+xxg*(a3+xxg*a4))
        enddo         ! iq loop
      enddo           ! k loop
    else
      do k=1,kl
        do iq=1,ifull    ! non Berm-Stan option
          ! Convert face index from 0:npanels to array indices
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
          jdel=int(yg(iq,k))
          yyg=yg(iq,k)-jdel
          ! Now make them proper indices in this processor's region
          idel = idel - ioff
          jdel = jdel - joff
          n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)

          r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)

          r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!              -y*(1+y)*c4/3}
!              +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop

          s(iq,k,:) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*r(:,2)-xxg*r(:,1)/3.)-xxg*(1.+xxg)*r(:,4)/3.)+xxg*(1.+xxg)*(2.-xxg)*r(:,3))/2.
        enddo         ! iq loop
      enddo           ! k loop
    endif             !  (mhint==2)
  else                ! (nfield<mh_bs)
    if(mhint==2)then ! Bessel interp
      do k=1,kl
        do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
          jdel=int(yg(iq,k))
          yyg=yg(iq,k)-jdel
          ! Now make them proper indices in this processor's region
          idel = idel - ioff
          jdel = jdel - joff
          n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth

          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!            -y*(1+y)*c4/3}
!            +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop

          a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
          a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
          sss = r(:,2)+.5*xxg*(r(:,3)-r(:,1)+xxg*(a3+xxg*a4))
          s(iq,k,:) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
        enddo         ! iq loop
      enddo           ! k loop
    else
      do k=1,kl
        do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
          idel=int(xg(iq,k))
          xxg=xg(iq,k)-idel
          jdel=int(yg(iq,k))
          yyg=yg(iq,k)-jdel
          ! Now make them proper indices in this processor's region
          idel = idel - ioff
          jdel = jdel - joff
          n = nface(iq,k) + noff ! Make this a local index

          if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
            cycle      ! Will be calculated on another processor
          end if

          c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
          c2 = sx(:,idel,jdel  ,n,k)
          c3 = sx(:,idel,jdel+1,n,k)
          c4 = sx(:,idel,jdel+2,n,k)
          cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
          cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth

          r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

          c1 = sx(:,idel+1,jdel-1,n,k)
          c2 = sx(:,idel+1,jdel  ,n,k)
          c3 = sx(:,idel+1,jdel+1,n,k)
          c4 = sx(:,idel+1,jdel+2,n,k)
          cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
          cmax = max(cmax,c2,c3) ! Bermejo & Staniforth

          r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)-yyg*(1.+yyg)*c4/3.)+yyg*(1.+yyg)*(2.-yyg)*c3)/2.

!         r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
!            -y*(1+y)*c4/3}
!            +y*(1+y)*(2-y)*c3}/2
          do nn=1,4,3   ! N.B.
            c2 = sx(:,idel+nn-2,jdel  ,n,k)
            c3 = sx(:,idel+nn-2,jdel+1,n,k)
            r(:,nn) = (1.-yyg)*c2 +yyg*c3
          enddo         ! nn loop
                
          sss = ((1.-xxg)*((2.-xxg)*((1.+xxg)*r(:,2)-xxg*r(:,1)/3.)-xxg*(1.+xxg)*r(:,4)/3.)+xxg*(1.+xxg)*(2.-xxg)*r(:,3))/2.
          s(iq,k,:) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
        enddo         ! iq loop
      enddo           ! k loop
    endif             !  (mhint==2)
  endif               ! (nfield<mh_bs)  .. else ..

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)
      
call END_LOG(ints_end)
return
end subroutine ints

subroutine ints_bl(s,intsch,nface,xg,yg)  ! not usually called

use cc_mpi
use indices_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'    ! has mh_bs

integer, parameter :: ntest=0
integer idel, iq, jdel, nn
integer i, j, k, n, ip, jp
integer ii
integer, intent(in) :: intsch
integer, dimension(ifull,kl), intent(in) :: nface
real xxg, yyg
real, dimension(ifull,kl), intent(in) :: xg,yg      ! now passed through call
real, dimension(0:ipan+1,0:jpan+1,1:npan,kl) :: sx
real, dimension(ifull+iextra,kl,1), intent(inout) :: s
real, dimension(ifull+iextra,kl) :: duma

!     this one does bi-linear interpolation only
!     first extend s arrays into sx - this one -1:il+2 & -1:il+2
!                    but for bi-linear only need 0:il+1 &  0:il+1
call START_LOG(ints_begin)
call bounds(s,corner=.true.)
sx(1:ipan,1:jpan,1:npan,1:kl) = reshape(s(1:ipan*jpan*npan,1:kl,1), (/ipan,jpan,npan,kl/))
do k=1,kl
  do n=1,npan
    do j=1,jpan
      sx(0,j,n,k)      = s(iw(1+(j-1)*ipan+(n-1)*ipan*jpan),k,1)
      sx(ipan+1,j,n,k) = s(ie(j*ipan+(n-1)*ipan*jpan),k,1)
    enddo               ! j loop
    do i=1,ipan
      sx(i,0,n,k)      = s(is(i+(n-1)*ipan*jpan),k,1)
      sx(i,jpan+1,n,k) = s(in(i-ipan+n*ipan*jpan),k,1)
    enddo               ! i loop

    sx(0,0,n,k)           = s(iws(1+(n-1)*ipan*jpan),k,1)
    sx(ipan+1,0,n,k)      = s(ies(ipan+(n-1)*ipan*jpan),k,1)
    sx(0,jpan+1,n,k)      = s(iwn(1-ipan+n*ipan*jpan),k,1)
    sx(ipan+1,jpan+1,n,k) = s(ien(n*ipan*jpan),k,1)
  enddo                  ! n loop
enddo                     ! k loop

! Loop over points that need to be calculated for other processes
do ii=neighnum,1,-1
  do iq=1,drlen(ii)
    !  Convert face index from 0:npanels to array indices
    ip = min(il_g,max(1,nint(dpoints(ii)%a(2,iq))))
    jp = min(il_g,max(1,nint(dpoints(ii)%a(3,iq))))
    n = nint(dpoints(ii)%a(1,iq)) + noff ! Local index
    !  Need global face index in fproc call
    idel = int(dpoints(ii)%a(2,iq))
    xxg = dpoints(ii)%a(2,iq) - idel
    jdel = int(dpoints(ii)%a(3,iq))
    yyg = dpoints(ii)%a(3,iq) - jdel
    k = nint(dpoints(ii)%a(4,iq))
    idel = idel - ioff
    jdel = jdel - joff
    sextra(ii)%a(iq) = yyg*( xxg*sx(idel+1,jdel+1,n,k)+(1.-xxg)*sx(  idel,jdel+1,n,k))         &
                      +(1.-yyg)*(   xxg*sx(idel+1,  jdel,n,k)+(1.-xxg)*sx(  idel,  jdel,n,k))
  end do
end do

call intssync_send(1)

do k=1,kl
  do iq=1,ifull
    ! Convert face index from 0:npanels to array indices
    idel=int(xg(iq,k))
    xxg=xg(iq,k)-idel
    jdel=int(yg(iq,k))
    yyg=yg(iq,k)-jdel
    ! Now make them proper indices in this processor's region
    idel = idel - ioff
    jdel = jdel - joff
    n = nface(iq,k) + noff ! Make this a local index

    if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
      cycle            ! Will be calculated on another processor
    end if

    s(iq,k,1) =      yyg*(      xxg*sx(idel+1,jdel+1,n,k)+(1.-xxg)*sx(idel,jdel+1,n,k))   &
                    +(1.-yyg)*(      xxg*sx(idel+1,jdel,n,k)+(1.-xxg)*sx(idel,jdel,n,k))
  enddo                  ! iq loop
end do                    ! k

call intssync_recv(s)

call END_LOG(ints_end)
return
end subroutine ints_bl
