! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------ 
    
!     this routine modified sept '06 to accept smaller il_g via abs(ik)
!     essentially to temporarily provide xx4, yy4, ax...bz for onthefly  
!     note that ax6 etc not needed for onthefly      

module setxyz_m

private
public setxyz

contains
    
subroutine setxyz(ik,rlong0,rlat0,schmidtin,x,y,z,wts,ax,ay,az,bx,by,bz,xx4,yy4, &
                  id,jd,ktau,ds)
    
use cc_mpi, only : indx, ccmpi_abort
use const_phys
use indices_m
use jimcc_m
use latlong_m
use map_m
use newmpar_m
use utilities
use workglob_m

implicit none

!     schmidt included
!     sets up x, y, z on sphere and unit u,v vectors
!     note that x,y,z have been normalized by rearth, the radius of the earth
!     suffix 6 denotes hex (6)

integer, intent(in) :: ik  ! passed as argument. Actual i dimension.
                           ! if negative, suppress calc of rlat4, rlong4, indices,em_g       
integer, intent(in) :: id, jd, ktau
integer i,j,n,ikk,idjd_g,iq
integer m
integer iq11,iq12,iq13,iq22,iq32,iqcc,iqnn
integer imin,imax,jmin,jmax,numpts
integer iqm,iqp, n_n, n_e, n_w, n_s
integer iquadx
integer, save :: num = 0
real(kind=8), dimension(:), pointer :: x, y, z    ! avoid intent for pointers
real(kind=8), dimension(:,:), pointer :: xx4, yy4 ! avoid intent for pointers
real(kind=8) alf,den1,xx,yy,zz,x4_iq_m,y4_iq_m,z4_iq_m
real(kind=8) xin,yin,zin
real(kind=8), parameter :: one = 1._8
real, dimension(ik*ik*6), intent(inout) :: wts
real, dimension(ik*ik*6), intent(inout) :: ax,ay,az
real, dimension(ik*ik*6), intent(inout) :: bx,by,bz
real, dimension(:), allocatable :: axx,ayy,azz,bxx,byy,bzz
real, dimension(:,:), allocatable :: em4,ax4,ay4,az4
real, dimension(3,3) :: rotpole
real rlong0,rlat0,schmidt,schmidtin
real dsfact
real den, dot,eps,dx2,dy2,sumwts,ratmin,ratmax,rat
real rlatdeg,rlondeg
real, intent(inout) :: ds

if ( size(x)/=ik*ik*6 ) then
  write(6,*) "ERROR: x argument is invalid in setxyz"
  call ccmpi_abort(-1)
end if

if ( size(y)/=ik*ik*6 ) then
  write(6,*) "ERROR: y argument is invalid in setxyz"
  call ccmpi_abort(-1)
end if

if ( size(z)/=ik*ik*6 ) then
  write(6,*) "ERROR: z argument is invalid in setxyz"
  call ccmpi_abort(-1)
end if

if ( size(xx4,1)/=1+4*ik .or. size(xx4,2)/=1+4*ik ) then
  write(6,*) "ERROR: xx4 argument is invalid in setxyz"
  call ccmpi_abort(-1)
end if

if ( size(yy4,1)/=1+4*ik .or. size(yy4,2)/=1+4*ik ) then
  write(6,*) "ERROR: yy4 argument is invalid in setxyz"
  call ccmpi_abort(-1)
end if

allocate( axx(ik*ik*6), ayy(ik*ik*6), azz(ik*ik*6) )
allocate( bxx(ik*ik*6), byy(ik*ik*6), bzz(ik*ik*6) )
allocate( em4(1+4*ik,1+4*ik), ax4(1+4*ik,1+4*ik), ay4(1+4*ik,1+4*ik), az4(1+4*ik,1+4*ik) )

dsfact=0.

num=num+1
!     When using the ifull_g notation: in_g, ie_g, iw_g and is_g give the
!     indices for the n, e, w, s neighbours respectively
!     a, b denote unit vectors in the direction of x, y (e & n) respectively
idjd_g = id+il_g*(jd-1)  ! Global value
schmidt=abs(schmidtin)
ikk=abs(ik)
iquadx=1+ik*((8*npanels)/(npanels+4))
      
em4=0. ! for cray compiler
xx=0.  ! for cray compiler
yy=0.  ! for cray compiler
      
! MJT notes - indices are now defined in indices_m.f90 as
! functions to save memory with global arrays
      
!----------------------------------------------------------------------------
! calculate grid information using quadruple resolution grid
call jimcc(em4,ax4,ay4,az4,xx4,yy4,ikk)

if(ktau<=1)then
  write(6,*) 'xx4 first & last ',xx4(1,1),xx4(iquadx,iquadx)
  write(6,*) 'xx4 (5,5),(7,7),(9,9) ',xx4(5,5),xx4(7,7),xx4(9,9)
  write(6,*)'yy4 first & last ',yy4(1,1),yy4(iquadx,iquadx)
  write(6,*) 'yy4 (5,5),(7,7),(9,9) ',yy4(5,5),yy4(7,7),yy4(9,9)
  write(6,*) 'xx4, yy4 central',xx4(2*ikk+1,2*ikk+1),yy4(2*ikk+1,2*ikk+1)
endif  ! (ktau<=1)

! rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
! rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
! rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
rotpole = calc_rotpole(rlong0,rlat0)

if ( schmidtin>=0. ) then
  ! following just for rlong4, rlat4 (& x4,y4,z4)
  alf=(one-schmidt**2)/(one+schmidt**2)
  do m=1,4
    do j=1,ikk
      do i=1,ikk
        select case(m)
          case(1)
            xx=xx4(4*i-1-1,4*j-1-1)
            yy=yy4(4*i-1-1,4*j-1-1)
          case(2)
            xx=xx4(4*i-1-1,4*j-1+1)
            yy=yy4(4*i-1-1,4*j-1+1)
          case(3)
            xx=xx4(4*i-1+1,4*j-1-1)
            yy=yy4(4*i-1+1,4*j-1-1)
          case(4)
            xx=xx4(4*i-1+1,4*j-1+1)
            yy=yy4(4*i-1+1,4*j-1+1)
        end select
        ! set up x0, y0, z0 coords on cube -1 to 1
        ! avoids earlier equivalencing of x,x0  etc
        x(indx(i,j,0,ikk,ikk))= 1.
        y(indx(i,j,0,ikk,ikk))=xx
        z(indx(i,j,0,ikk,ikk))=yy
        x(indx(i,j,3,ikk,ikk))=-1.
        z(indx(i,j,3,ikk,ikk))=-xx
        y(indx(i,j,3,ikk,ikk))=-yy
        x(indx(i,j,1,ikk,ikk))=-yy
        y(indx(i,j,1,ikk,ikk))=xx
        z(indx(i,j,1,ikk,ikk))= 1.
        y(indx(i,j,4,ikk,ikk))=-yy
        x(indx(i,j,4,ikk,ikk))=xx
        z(indx(i,j,4,ikk,ikk))=-1.
        x(indx(i,j,2,ikk,ikk))=-yy
        y(indx(i,j,2,ikk,ikk))= 1.
        z(indx(i,j,2,ikk,ikk))=-xx
        z(indx(i,j,5,ikk,ikk))=yy
        y(indx(i,j,5,ikk,ikk))=-1.
        x(indx(i,j,5,ikk,ikk))=xx
      enddo  ! i loop
    enddo   ! j loop
    do iq=1,6*ik*ik
      call norm8(x(iq),y(iq),z(iq),den1) ! x, y, z are coords on sphere  -1 to 1
      x4_iq_m=x(iq)*schmidt*(1.+alf)/(1.+alf*z(iq))
      y4_iq_m=y(iq)*schmidt*(1.+alf)/(1.+alf*z(iq))
      z4_iq_m=(alf+z(iq))/(1.+alf*z(iq)) 
      ! here is calculation of rlong4, rlat4
      ! also provide latitudes and longitudes (-pi to pi)
      !if(rlong0==0..and.rlat0==90.)then
      !  xx=x4_iq_m
      !  yy=y4_iq_m
      !  zz=z4_iq_m
      !else
        ! x4(), y4(z), z4() are "local" coords with z4 out of central panel
        ! while xx, yy, zz are "true" Cartesian values
        ! xx is new x after rot by rlong0 then rlat0
        xx=rotpole(1,1)*x4_iq_m+rotpole(1,2)*y4_iq_m+rotpole(1,3)*z4_iq_m
        yy=rotpole(2,1)*x4_iq_m+rotpole(2,2)*y4_iq_m+rotpole(2,3)*z4_iq_m
        zz=rotpole(3,1)*x4_iq_m+rotpole(3,2)*y4_iq_m+rotpole(3,3)*z4_iq_m
      !endif
      rlat4(iq,m)=real(asin(zz))
      !if(yy/=0..or.xx/=0.)then
        rlong4(iq,m)=real(atan2(yy,xx))                       ! N.B. -pi to pi
        if(rlong4(iq,m)<0.)  then
          rlong4(iq,m)=rlong4(iq,m)+2.*pi ! 0 to 2*pi  09-25-1997
        end if
      !else
      !  rlong4(iq,m)=0.    ! a default value for NP/SP
      !endif
      ! convert long4 and lat4 (used by cctocc4) to degrees	
      rlat4(iq,m)=rlat4(iq,m)*180./pi
      rlong4(iq,m)=rlong4(iq,m)*180./pi
    enddo   ! iq loop
  enddo    ! m loop

  dsfact=4*ikk/(2.*pi)     ! con-cube
  ds=rearth/dsfact
  ! extend em4 to uppermost i and j rows
  do j=1,4*ikk
    em4(iquadx,j)=em4(1,j)
  enddo
  do i=1,4*ikk
    em4(i,iquadx)=em4(i,1)
  enddo
  do j=1,ikk
    do i=1,ikk
      do n=0,5
        ! average Purser em is pi/2
        em_g(indx(i,j,n,ikk,ikk))=pi/(2.*em4(4*i-1,4*j-1))
        emu_g(indx(i,j,n,ikk,ikk))=pi/(2.*em4(4*i+1,4*j-1))
        emv_g(indx(i,j,n,ikk,ikk))=pi/(2.*em4(4*i-1,4*j+1))
      enddo ! n loop
      ax(indx(i,j,0,ikk,ikk))=ax4(4*i-1,4*j-1)
      ay(indx(i,j,0,ikk,ikk))=ay4(4*i-1,4*j-1)
      az(indx(i,j,0,ikk,ikk))=az4(4*i-1,4*j-1)
    enddo  ! i loop
  enddo   ! j loop
  if(ktau<=1)then
    write(6,*) 'ax6 (1,1,0) & (2,2,0) ',ax(indx(1,1,0,ikk,ikk)),ax(indx(2,2,0,ikk,ikk))
    write(6,*) 'ay6 (1,1,0) & (2,2,0) ',ay(indx(1,1,0,ikk,ikk)),ay(indx(2,2,0,ikk,ikk))
    write(6,*) 'az6 (1,1,0) & (2,2,0) ',az(indx(1,1,0,ikk,ikk)),az(indx(2,2,0,ikk,ikk))
  endif ! (ktau<=1)
end if ! (schmidt>=0.)
      
do j=1,ikk
  do i=1,ikk
    xx=xx4(4*i-1,4*j-1)
    yy=yy4(4*i-1,4*j-1)

    ! set up x0, y0, z0 coords on cube -1 to 1
    ! avoids earlier equivalencing of x,x0  etc
    x(indx(i,j,0,ikk,ikk))= 1.
    y(indx(i,j,0,ikk,ikk))=xx
    z(indx(i,j,0,ikk,ikk))=yy
    x(indx(i,j,3,ikk,ikk))=-1.
    z(indx(i,j,3,ikk,ikk))=-xx
    y(indx(i,j,3,ikk,ikk))=-yy
    x(indx(i,j,1,ikk,ikk))=-yy
    y(indx(i,j,1,ikk,ikk))=xx
    z(indx(i,j,1,ikk,ikk))= 1.
    y(indx(i,j,4,ikk,ikk))=-yy
    x(indx(i,j,4,ikk,ikk))=xx
    z(indx(i,j,4,ikk,ikk))=-1.
    x(indx(i,j,2,ikk,ikk))=-yy
    y(indx(i,j,2,ikk,ikk))= 1.
    z(indx(i,j,2,ikk,ikk))=-xx
    z(indx(i,j,5,ikk,ikk))=yy
    y(indx(i,j,5,ikk,ikk))=-1.
    x(indx(i,j,5,ikk,ikk))=xx
  enddo  ! i loop
enddo   ! j loop
do iq=1,ik*ik*6
  call norm8(x(iq),y(iq),z(iq),den1) ! x, y, z are coords on sphere  -1 to 1
enddo   ! iq loop
if(ktau==0) write(6,*) 'basic grid length ds =',ds

!if (schmidt/=1.) then
  alf=(one-schmidt**2)/(one+schmidt**2)
  write(6,*) 'doing schmidt with schmidt,alf: ',schmidt,alf
  if (schmidtin>0.) then
    do iq=1,ik*ik*6
      xin=x(iq)
      yin=y(iq)
      zin=z(iq)
      x(iq)=xin*schmidt*(1.+alf)/(1.+alf*zin)
      y(iq)=yin*schmidt*(1.+alf)/(1.+alf*zin)
      z(iq)=(alf+zin)/(1.+alf*zin)
      em_g(iq)=real(real(em_g(iq),8)*real(schmidt,8)*(1.+alf*zin)/(1.-alf))
    enddo   ! iq loop
    do iq=1,ifull_g
      ! with schmidt, for ntang=1 or 2 must average em_g to get emu_g & emv_g
      emu_g(iq)=.5*(em_g(iq)+em_g(ie_g(iq)))
      emv_g(iq)=.5*(em_g(iq)+em_g(in_g(iq)))
    enddo   ! iq loop
  else
    do iq=1,ik*ik*6
      xin=x(iq)
      yin=y(iq)
      zin=z(iq)
      x(iq)=xin*schmidt*(1.+alf)/(1.+alf*zin)
      y(iq)=yin*schmidt*(1.+alf)/(1.+alf*zin)
      z(iq)=(alf+zin)/(1.+alf*zin)
    enddo   ! iq loop
  end if      
!endif      ! (schmidt.ne.1.)

write(6,*) 'ktau,ikk,schmidtin ',ktau,ikk,schmidtin
if(ktau==0.and.schmidtin>0.)then
  write(6,*) 'On each panel (ntang=0)_em_g for (1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
  do n=0,npanels
    iq11=indx(1,1,n,ikk,ikk)
    iq12=indx(1,2,n,ikk,ikk)
    iq13=indx(1,3,n,ikk,ikk)
    iq22=indx(2,2,n,ikk,ikk)
    iq32=indx(3,2,n,ikk,ikk)
    iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
    iqnn=indx(ikk,ikk,n,ikk,ikk)
    write(6,'(i3,7f8.3)') n,em_g(iq11),em_g(iq12),em_g(iq13),em_g(iq22),em_g(iq32),em_g(iqcc),em_g(iqnn)
  enddo
  write(6,*) 'On each panel (ntang=0)_emu_g for (1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
  do n=0,npanels
    iq11=indx(1,1,n,ikk,ikk)
    iq12=indx(1,2,n,ikk,ikk)
    iq13=indx(1,3,n,ikk,ikk)
    iq22=indx(2,2,n,ikk,ikk)
    iq32=indx(3,2,n,ikk,ikk)
    iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
    iqnn=indx(ikk,ikk,n,ikk,ikk)
    write(6,'(i3,7f8.3)') n,emu_g(iq11),emu_g(iq12),emu_g(iq13),emu_g(iq22),emu_g(iq32),emu_g(iqcc),emu_g(iqnn)
  enddo
  write(6,*) 'On each panel (ntang=0)_emv_g for (1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
  do n=0,npanels
    iq11=indx(1,1,n,ikk,ikk)
    iq12=indx(1,2,n,ikk,ikk)
    iq13=indx(1,3,n,ikk,ikk)
    iq22=indx(2,2,n,ikk,ikk)
    iq32=indx(3,2,n,ikk,ikk)
    iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
    iqnn=indx(ikk,ikk,n,ikk,ikk)
    write(6,'(i3,7f8.3)') n,emv_g(iq11),emv_g(iq12),emv_g(iq13),emv_g(iq22),emv_g(iq32),emv_g(iqcc),emv_g(iqnn)
  enddo
endif  ! (ktau.eq.0.and.schmidtin>0.)

! set up vectors in direction of u and v
if(schmidtin<0.)then
  ! same as below but avoids recalc in_g, ie_g, iw_g, is_g)  
  do n=0,5,2
    n_w=mod(n+5,6)
    n_e=mod(n+2,6)
    n_n=mod(n+1,6)
    n_s=mod(n+4,6)
    do j=1,ikk
      do i=1,ikk
        iq=indx(i,j,n,ikk,ikk)
        iqp=iq+1
        iqm=iq-1
        if(i==1)iqm=indx(ikk,j,n_w,ikk,ikk)
        if(i==ikk)iqp=indx(ikk+1-j,1,n_e,ikk,ikk)
        ax(iq)=real(x(iqp)-x(iqm))
        az(iq)=real(z(iqp)-z(iqm))
        ay(iq)=real(y(iqp)-y(iqm))
        iqp=iq+ikk
        iqm=iq-ikk
        if(j==1)iqm=indx(ikk,ikk+1-i,n_s,ikk,ikk)
        if(j==ikk)iqp=indx(i,1,n_n,ikk,ikk)
        bx(iq)=real(x(iqp)-x(iqm))
        by(iq)=real(y(iqp)-y(iqm))
        bz(iq)=real(z(iqp)-z(iqm))
      enddo
    enddo
  enddo
  do n=1,5,2
    n_w=mod(n+4,6)
    n_e=mod(n+1,6)
    n_n=mod(n+2,6)
    n_s=mod(n+5,6)
    do j=1,ikk
      do i=1,ikk
        iq=indx(i,j,n,ikk,ikk)
        iqp=iq+1
        iqm=iq-1
        if(i==1)iqm=indx(ikk+1-j,ikk,n_w,ikk,ikk)
        if(i==ikk)iqp=indx(1,j,n_e,ikk,ikk)
        ax(iq)=real(x(iqp)-x(iqm))
        ay(iq)=real(y(iqp)-y(iqm))
        az(iq)=real(z(iqp)-z(iqm))
        iqp=iq+ikk
        iqm=iq-ikk
        if(j==1)iqm=indx(i,ikk,n_s,ikk,ikk)
        if(j==ikk)iqp=indx(1,ikk+1-i,n_n,ikk,ikk)
        bx(iq)=real(x(iqp)-x(iqm))
        by(iq)=real(y(iqp)-y(iqm))
        bz(iq)=real(z(iqp)-z(iqm))
      enddo
    enddo
  enddo
else  !  usual with (schmidtin>0. (but equiv. to above)
  do iq=1,ifull_g
    ! first guess tang vectors by finite differences
    ax(iq)=real(x(ie_g(iq))-x(iw_g(iq)))
    ay(iq)=real(y(ie_g(iq))-y(iw_g(iq)))
    az(iq)=real(z(ie_g(iq))-z(iw_g(iq)))
    bx(iq)=real(x(in_g(iq))-x(is_g(iq)))
    by(iq)=real(y(in_g(iq))-y(is_g(iq)))
    bz(iq)=real(z(in_g(iq))-z(is_g(iq)))
  enddo   ! iq loop
endif   ! (schmidtin<0.  ... else)
! form axx and bxx tangential to the sphere
call cross3b(axx,ayy,azz, bx,by,bz, x,y,z, ik*ik*6)
call cross3(bxx,byy,bzz, x,y,z, ax,ay,az, ik*ik*6)
do iq=1,ikk*ikk*6
  call norm(axx(iq),ayy(iq),azz(iq),den)
  call norm(bxx(iq),byy(iq),bzz(iq),den)
  ! make sure they are perpendicular & normalize
  dot=axx(iq)*bxx(iq)+ayy(iq)*byy(iq)+azz(iq)*bzz(iq)
  eps=-dot/(1.+sqrt(1.-dot*dot))
  ax(iq)=axx(iq)+eps*bxx(iq)
  ay(iq)=ayy(iq)+eps*byy(iq)
  az(iq)=azz(iq)+eps*bzz(iq)
  bx(iq)=bxx(iq)+eps*axx(iq)
  by(iq)=byy(iq)+eps*ayy(iq)
  bz(iq)=bzz(iq)+eps*azz(iq)
  call norm(ax(iq),ay(iq),az(iq),den)
  call norm(bx(iq),by(iq),bz(iq),den)
enddo   ! iq loop
if (schmidtin>=0.) then  ! (schmidt<0. is to finish of stuff needed for onthefly
        
  do iq=1,ifull_g
    ! calculate inverse of emu_g & emv_g first
    dx2=real((x(ie_g(iq))-x(iq))**2+(y(ie_g(iq))-y(iq))**2+(z(ie_g(iq))-z(iq))**2)
    ! include arc-length corrn using 2*arcsin(theta/2)
    emu_g(iq)=sqrt(dx2)*(1.+dx2/24.) *dsfact
    dy2=real((x(in_g(iq))-x(iq))**2+(y(in_g(iq))-y(iq))**2+(z(in_g(iq))-z(iq))**2)
    emv_g(iq)=sqrt(dy2)*(1.+dy2/24.) *dsfact
  enddo   ! iq loop
  do iq=1,ifull_g   ! based on inverse values of emu_g & emv_g
    if (isv2_g(iq)<1.and.iwu2_g(iq)>ifull_g) then                                              ! MJT bug fix
      em_g(iq)=4./(emv_g(iwu2_g(iq)-ifull_g)+emu_g(iq)+emu_g(isv2_g(iq)+ifull_g)+emv_g(iq))    ! MJT bug fix
    else if (isv2_g(iq)<1) then                                                                ! MJT bug fix
      em_g(iq)=4./(emu_g(iwu2_g(iq))+emu_g(iq)+emu_g(isv2_g(iq)+ifull_g)+emv_g(iq))            ! MJT bug fix
    else if (iwu2_g(iq)>ifull_g) then                                                          ! MJT bug fix
      em_g(iq)=4./(emv_g(iwu2_g(iq)-ifull_g)+emu_g(iq)+emv_g(isv2_g(iq))+emv_g(iq))            ! MJT bug fix
    else                                                                                       ! MJT bug fix
      em_g(iq)=4./(emu_g(iwu2_g(iq))+emu_g(iq)+emv_g(isv2_g(iq))+emv_g(iq))                    ! MJT bug fix
    end if                                                                                     ! MJT bug fix
  enddo   ! iq loop
  do iq=1,ifull_g
    emu_g(iq)=1./emu_g(iq)
    emv_g(iq)=1./emv_g(iq)
  enddo   ! iq loop
      
  if(ktau==0)then
    do iq=ikk-2,ikk
      write(6,*) 'iq,em_g,emu_g,emv_g',iq,em_g(iq),emu_g(iq),emv_g(iq)
    enddo   ! iq loop
    if(id<=ikk.and.jd<=jl)then
      iq=id+ikk*(jd-1)
      write(6,*) 'values at idjd'
      write(6,*) 'iq,x,y,z',iq,x(iq),y(iq),z(iq)
      write(6,*) 'iq,ax,ay,az',iq,ax(iq),ay(iq),az(iq)
      write(6,*) 'iq,bx,by,bz',iq,bx(iq),by(iq),bz(iq)
      write(6,*) 'values at in_g(idjd)'
      write(6,*) 'iq,x,y,z',in_g(iq),x(in_g(iq)),y(in_g(iq)),z(in_g(iq))
      write(6,*) 'iq,ax,ay,az',in_g(iq),ax(in_g(iq)),ay(in_g(iq)),az(in_g(iq))
      write(6,*) 'iq,bx,by,bz',in_g(iq),bx(in_g(iq)),by(in_g(iq)),bz(in_g(iq))
      write(6,*) 'values at ie_g(idjd)'
      write(6,*) 'iq,x,y,z',ie_g(iq),x(ie_g(iq)),y(ie_g(iq)),z(ie_g(iq))
      write(6,*) 'iq,ax,ay,az',ie_g(iq),ax(ie_g(iq)),ay(ie_g(iq)),az(ie_g(iq))
      write(6,*) 'iq,bx,by,bz',ie_g(iq),bx(ie_g(iq)),by(ie_g(iq)),bz(ie_g(iq))
      write(6,*) 'values at iw_g(idjd)'
      write(6,*) 'iq,x,y,z',iw_g(iq),x(iw_g(iq)),y(iw_g(iq)),z(iw_g(iq))
      write(6,*) 'iq,ax,ay,az',iw_g(iq),ax(iw_g(iq)),ay(iw_g(iq)),az(iw_g(iq))
      write(6,*) 'iq,bx,by,bz',iw_g(iq),bx(iw_g(iq)),by(iw_g(iq)),bz(iw_g(iq))
      write(6,*) 'values at is_g(idjd)'
      write(6,*) 'iq,x,y,z',is_g(iq),x(is_g(iq)),y(is_g(iq)),z(is_g(iq))
      write(6,*) 'iq,ax,ay,az',is_g(iq),ax(is_g(iq)),ay(is_g(iq)),az(is_g(iq))
      write(6,*) 'iq,bx,by,bz',is_g(iq),bx(is_g(iq)),by(is_g(iq)),bz(is_g(iq))
    endif
  endif  ! (ktau==0)

  !     calculate approx areas around each grid point
  !     just used for error diagnostics
  !     now use 1/(em_g**2) to cope with schmidt, rotated and ocatagon coordinates
  sumwts=0.
  do iq=1,ifull_g
    wts(iq)=1./em_g(iq)**2
    sumwts=sumwts+wts(iq)
    ! cosa is dot product of unit vectors
    ! *** only useful as diagnostic for gnew
    ! cosa(iq)=ax(iq)*bx(iq)+ay(iq)*by(iq)+az(iq)*bz(iq)
  enddo   ! iq loop
  if(ktau==0)then
    write(6,*) 'sumwts/ifull_g ',sumwts/ifull_g  ! ideally equals 4*pi ??
    write(6,*) 'in setxyz rlong0,rlat0,schmidt ',rlong0,rlat0,schmidt
  endif  ! (ktau==0)

  do iq=1,ifull_g
    ! scale wts so sum over globe is 1.
    ! wts(iq)=wts(iq)/(6.*sumwts)  ! for old conf-cub defn
    wts(iq)=wts(iq)/sumwts
    ! also provide latitudes and longitudes (-pi to pi)
    !if(rlong0==0..and.rlat0==90.)then
    !  xx=x(iq)
    !  yy=y(iq)
    !  zz=z(iq)
    !else
      ! x(), y(z), z() are "local" coords with z out of central panel
      ! while xx, yy, zz are "true" Cartesian values
      ! xx is new x after rot by rlong0 then rlat0
      xx=rotpole(1,1)*x(iq)+rotpole(1,2)*y(iq)+rotpole(1,3)*z(iq)
      yy=rotpole(2,1)*x(iq)+rotpole(2,2)*y(iq)+rotpole(2,3)*z(iq)
      zz=rotpole(3,1)*x(iq)+rotpole(3,2)*y(iq)+rotpole(3,3)*z(iq)
    !endif
    ! f_g(iq)=2. *2.*pi *(z(iq)/rdiv) /86400.
    rlatt_g(iq)=real(asin(zz))
    f_g(iq)=real(2._8*2._8*real(pi,8)*zz/86400._8)  !  zz along "true" N-S axis
    !if(yy/=0..or.xx/=0.)then
      rlongg_g(iq)=real(atan2(yy,xx))               ! N.B. -pi to pi
      if(rlongg_g(iq)<0.) then
        rlongg_g(iq)=rlongg_g(iq)+2.*pi ! 0 to 2*pi  09-25-1997
      end if
    !else
    !  rlongg_g(iq)=0.    ! a default value for NP/SP
    !endif
  enddo   ! iq loop
  if(ktau==0)then
    write(6,*) 'At centre of the faces:'
    do n=0,npanels
      iq=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
      write(6,'(" n,iq,x,y,z,long,lat,f ",i2,i7,3f7.3,2f8.2,f9.5)') n,iq,x(iq),y(iq),z(iq),   &
          rlongg_g(iq)*180./pi,rlatt_g(iq)*180./pi,f_g(iq)
    enddo
    write(6,*) 'At mid-x along edges:'
    do n=0,npanels
      iq=indx((ikk+1)/2,1,n,ikk,ikk)
      write(6,'(" n,iq,x,y,z,long,lat,f_g ",i2,i7,3f7.3,2f8.2,f9.5)') n,iq,x(iq),y(iq),z(iq), &
         rlongg_g(iq)*180./pi,rlatt_g(iq)*180./pi,f_g(iq)
    enddo
    write(6,*) 'At mid-y along edges:'
    do n=0,npanels
      iq=indx(1,(ikk+1)/2,n,ikk,ikk)
      write(6,'(" n,iq,x,y,z,long,lat,f_g ",i2,i7,3f7.3,2f8.2,f9.5)') n,iq,x(iq),y(iq),z(iq), &
         rlongg_g(iq)*180./pi,rlatt_g(iq)*180./pi,f_g(iq)
    enddo
    write(6,*) 'On each panel final_em_g for (1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
    do n=0,npanels
      iq11=indx(1,1,n,ikk,ikk)
      iq12=indx(1,2,n,ikk,ikk)
      iq13=indx(1,3,n,ikk,ikk)
      iq22=indx(2,2,n,ikk,ikk)
      iq32=indx(3,2,n,ikk,ikk)
      iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
      iqnn=indx(ikk,ikk,n,ikk,ikk)
      write(6,'(i3,7f8.3)') n,em_g(iq11),em_g(iq12),em_g(iq13),em_g(iq22),em_g(iq32),em_g(iqcc),em_g(iqnn)
    enddo
    write(6,*) 'On each panel final_emu_g for (1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
    do n=0,npanels
      iq11=indx(1,1,n,ikk,ikk)
      iq12=indx(1,2,n,ikk,ikk)
      iq13=indx(1,3,n,ikk,ikk)
      iq22=indx(2,2,n,ikk,ikk)
      iq32=indx(3,2,n,ikk,ikk)
      iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
      iqnn=indx(ikk,ikk,n,ikk,ikk)
      write(6,'(i3,7f8.3)') n,emu_g(iq11),emu_g(iq12),emu_g(iq13),emu_g(iq22),emu_g(iq32),emu_g(iqcc),emu_g(iqnn)
    enddo
    write(6,*) 'On each panel final_emv_g for (1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
    do n=0,npanels
      iq11=indx(1,1,n,ikk,ikk)
      iq12=indx(1,2,n,ikk,ikk)
      iq13=indx(1,3,n,ikk,ikk)
      iq22=indx(2,2,n,ikk,ikk)
      iq32=indx(3,2,n,ikk,ikk)
      iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
      iqnn=indx(ikk,ikk,n,ikk,ikk)
      write(6,'(i3,7f8.3)') n,emv_g(iq11),emv_g(iq12),emv_g(iq13),emv_g(iq22),emv_g(iq32),emv_g(iqcc),emv_g(iqnn)
    enddo
  endif  ! (ktau.eq.0)
  do iq=1,ifull_g   ! set up Coriolis
    fu_g(iq)=(f_g(iq)+f_g(ie_g(iq)))*.5
    fv_g(iq)=(f_g(iq)+f_g(in_g(iq)))*.5
  enddo   ! iq loop
  do iq=1,ifull_g   ! average map factor derivs needed for nxmap=1
    dmdx_g(iq)=.5*(em_g(ie_g(iq))-em_g(iw_g(iq)))/ds  
    dmdy_g(iq)=.5*(em_g(in_g(iq))-em_g(is_g(iq)))/ds  
  enddo   ! iq loop

  ratmin=100.
  ratmax=0.
  do n=0,npanels
    do i=1,ikk
      iq=indx(i,ikk/2,n,ikk,ikk)
      rat=em_g(iq)/em_g(ie_g(iq))
      if(rat<ratmin)then
        ratmin=rat
        imin=i
      endif
      if(rat>ratmax)then
        ratmax=rat
        imax=i
      endif
    enddo
    if(num==1)then
      write(6,*) 'em_g ratio for j=ikk/2 on npanel ',n
      write (6,"(12f6.3)") (em_g(indx(i,ikk/2,n,ikk,ikk))/em_g(ie_g(indx(i,ikk/2,n,ikk,ikk))),i=1,ikk)
    endif
  enddo
  write(6,*) 'for j=ikk/2 & myid=0, ratmin,ratmax = ',ratmin,ratmax
  write(6,*) 'with imin,imax ',imin,imax
  ratmin=100.
  ratmax=0.
  do n=0,npanels
    do j=1,ikk
      do i=1,ikk
        iq=indx(i,j,n,ikk,ikk)
        rat=em_g(iq)/em_g(ie_g(iq))
        if(rat<ratmin)then
          ratmin=rat
          imin=i
          jmin=j
        endif
        if(rat>ratmax)then
          ratmax=rat
          imax=i
          jmax=j
        endif
      enddo
    enddo
  enddo
  write(6,*) 'for all j & myid=0, ratmin,ratmax = ',ratmin,ratmax
  write(6,*) 'with imin,jmin,imax,jmax ',imin,jmin,imax,jmax
  write (6,"('1st 10 ratios',10f6.3)") (em_g(iq)/em_g(ie_g(iq)),iq=1,10)
  numpts=0
  do iq=1,ifull_g
    rlatdeg=rlatt_g(iq)*180./pi
    rlondeg=rlongg_g(iq)*180./pi
    if(rlatdeg>20..and.rlatdeg<60..and.rlondeg>230..and.rlatdeg<300.)numpts=numpts+1
  enddo
  write(6,*) 'points in SGMIP region ',numpts

end if

deallocate( axx, ayy, azz )
deallocate( bxx, byy, bzz )
deallocate( em4, ax4, ay4, az4 )

return
end subroutine setxyz

subroutine norm(a,b,c,den)

implicit none

real, intent(inout) :: a,b,c,den

den=sqrt(a*a+b*b+c*c)
a=a/den
b=b/den
c=c/den

return
end subroutine norm

subroutine norm8(a,b,c,den)

implicit none

real(kind=8), intent(inout) :: a,b,c,den

den=sqrt(a*a+b*b+c*c)
a=a/den
b=b/den
c=c/den

return
end subroutine norm8
    
subroutine cross3(c1,c2,c3,a1,a2,a3,b1,b2,b3, ifull_g)
!     calculate vector components of c = a x b
!     where each RHS component represents 3 vector components
!     this one need not have contiguous memory in common

implicit none

integer, intent(in) :: ifull_g
real(kind=8), dimension(ifull_g) :: a1,a2,a3
real, dimension(ifull_g) :: b1,b2,b3
real, dimension(ifull_g) :: c1,c2,c3

c1=real(a2*real(b3,8)-real(b2,8)*a3)
c2=real(a3*real(b1,8)-real(b3,8)*a1)
c3=real(a1*real(b2,8)-real(b1,8)*a2)

return
end subroutine cross3

subroutine cross3b(c1,c2,c3,a1,a2,a3,b1,b2,b3, ifull_g)
!     calculate vector components of c = a x b
!     where each RHS component represents 3 vector components
!     this one need not have contiguous memory in common

implicit none

integer, intent(in) :: ifull_g
real(kind=8), dimension(ifull_g), intent(in) :: b1,b2,b3
real, dimension(ifull_g), intent(in) :: a1,a2,a3
real, dimension(ifull_g), intent(out) :: c1,c2,c3      

c1=real(real(a2,8)*b3-b2*real(a3,8))
c2=real(real(a3,8)*b1-b3*real(a1,8))
c3=real(real(a1,8)*b2-b1*real(a2,8))

return
end subroutine cross3b

end module setxyz_m