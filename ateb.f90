
! This code was originally based on the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)

! Usual pratice is:
!   call tebinit     ! to initalise state arrays, etc
!   call tebdefault  ! to set state arrays to default values
!   call tebload     ! to load previous state arrays (from tebsave)
!   call tebtype     ! to define urban type (or use tebfndef instead)
!   ...
!   call tebnewangle ! update solar zenith angle (use tebccangle for CCAM)
!   call tebalb      ! to return urban contrabution to albedo
!   ...
!   call tebcalc     ! to return fluxes, etc
!   call tebzom      ! includes urban contrabution to roughness length
!   ...
!   call tebsave     ! to save current state arrays (for use by tebload)
!   call tebend      ! to deallocate memory before quiting


! NOTES: 
!  Below are differences with TEB (Masson 2000) scheme and aTEB:
!
! - Currently snow is neglected for urban cover.
!
! - aTEB uses two walls instead of the TEB single wall.  Also, only up to 2nd order reflections are used for
!   longwave and short wave radation in aTEB.  In TEB, infinite reflections are used for shortwave, but not long wave.
!   The second wall allows the canyon to be orentated wrt the path of the sun (i.e., for planned urban model experiments),
!   as well as supporting the TEB approach where all orentiations are assumed.
!
! - aTEB allows model levels below the building height (i.e., inside the canyon).  This can probably be improved
!   over time.
!
! - traffic and industry sensible fluxes as well as building confort temperatures should be time varying.
!


module ateb

implicit none

private
public tebinit,tebcalc,tebend,tebzom,tebload,tebsave,tebdefault,tebtype,tebalb,tebfndef, &
       tebnewangle,tebccangle

integer ufull,maxtype
integer, dimension(:), allocatable :: ugrid,qgrid
real, dimension(:,:), allocatable :: rooftemp,walletemp,wallwtemp,roadtemp
real, dimension(:,:), allocatable :: roofadjt,walleadjt,wallwadjt,roadadjt
real, dimension(:), allocatable :: roofadjw,roadadjw
real, dimension(:), allocatable :: roofwater,roadwater
real, dimension(:), allocatable :: fnsigmabld,fnhwratio,fnindustryfg,fntrafficfg
real, dimension(:), allocatable :: fnzo,fnbldheight
real, dimension(:), allocatable :: vangle,hangle
real, dimension(3) :: roofdepth,walldepth,roaddepth
real, dimension(3) :: roofcp,wallcp,roadcp
real, dimension(3) :: rooflambda,walllambda,roadlambda
real, dimension(0:220) :: table
real, parameter :: aircp=1004.64   ! Specific heat of dry air
real, parameter :: bldtemp=291.16  ! Comfort temperature = 18deg C
real, parameter :: grav=9.80616    ! gravity
real, parameter :: lv=2.5104e6     ! Latent heat of vaporisation
real, parameter :: pi=3.1415927    ! pi
real, parameter :: rd=287.04       ! Gas constant for dry air
real, parameter :: sbconst=5.67e-8 ! Stefan-Boltzmann constant    
real, parameter :: zoroof=0.15     ! Roughness length for rooftops (see Masson 2000)
real, parameter :: roofemiss=0.90
real, parameter :: wallemiss=0.85
real, parameter :: roademiss=0.94
real, parameter :: roofalpha=0.15
real, parameter :: wallalpha=0.25
real, parameter :: roadalpha=0.08
real, parameter :: maxroofwater=1.
real, parameter :: maxroadwater=1.

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepainvres the arrays used by the TEB scheme
! This is a compulsory subroutine that must be called during
! model initalisation

subroutine tebinit(ifull,sigmau,rlon,rlat,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq,iq1,iqmark(1)
integer, dimension(ifull) :: utype
real, dimension(ifull), intent(in) :: sigmau,rlon,rlat
real, dimension(:), allocatable :: dis

ufull=count(sigmau.gt.0.)
if (ufull.eq.0) return

if (diag.ne.0) write(6,*) "Initialising aTEB"

allocate(ugrid(ufull),qgrid(ifull))
allocate(rooftemp(ufull,3),walletemp(ufull,3),wallwtemp(ufull,3),roadtemp(ufull,3))
allocate(roofadjt(ufull,3),walleadjt(ufull,3),wallwadjt(ufull,3),roadadjt(ufull,3))
allocate(roofadjw(ufull),roadadjw(ufull))
allocate(roofwater(ufull),roadwater(ufull))
allocate(fnhwratio(ufull),fnsigmabld(ufull),fnindustryfg(ufull))
allocate(fntrafficfg(ufull),fnbldheight(ufull),fnzo(ufull))
allocate(vangle(ufull),hangle(ufull))

! Other paramaters

roofcp(1)=2.11E6
roofcp(2)=0.28E6
roofcp(3)=0.29E6 ! insulation
wallcp(1)=1.55E6
wallcp(2)=1.55E6
wallcp(3)=0.29E6 ! insulation
roadcp(1)=1.94E6
roadcp(2)=1.28E6
roadcp(3)=1.28E6
roofdepth(1)=0.05
roofdepth(2)=0.4
roofdepth(3)=0.1
walldepth(1)=0.02
walldepth(2)=0.125
walldepth(3)=0.05
roaddepth(1)=0.05
roaddepth(2)=0.1
roaddepth(3)=1.
rooflambda(1)=1.51
rooflambda(2)=0.08
rooflambda(3)=0.05
walllambda(1)=0.9338
walllambda(2)=0.9338
walllambda(3)=0.05
roadlambda(1)=0.7454
roadlambda(2)=0.2513
roadlambda(3)=0.2513


iqu=0
do iq=1,ifull
  if (sigmau(iq).gt.0.) then
    iqu=iqu+1
    ugrid(iqu)=iq
  end if
end do

if (ufull.gt.0.) then
  allocate(dis(ufull))
  do iq=1,ifull
    do iqu=1,ufull ! since ufull << ifull
      iq1=ugrid(iqu)
      dis(iqu)=abs(rlon(iq)-rlon(iq1))
      if (dis(iqu).gt.pi) dis(iqu)=2.*pi-dis(iqu)
      dis(iqu)=dis(iqu)**2+(rlat(iq)-rlat(iq1))**2
    end do
    iqmark=minloc(dis)
    qgrid(iq)=iqmark(1)
  end do
  deallocate(dis)
else
  qgrid=1
end if

! Initialise state variables
rooftemp=bldtemp
walletemp=bldtemp
wallwtemp=bldtemp
roadtemp=bldtemp
roofwater=0.
roadwater=0.
vangle=0.
hangle=0.

utype=1
call tebtype(ifull,utype,0)

table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
table(219:220)=(/ 29845.0, 31169.0 /)  

return
end subroutine tebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine tebend(diag)

implicit none

integer, intent(in) :: diag

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Deallocating aTEB arrays"

deallocate(ugrid,qgrid)
deallocate(rooftemp,walletemp,wallwtemp,roadtemp)
deallocate(roofadjt,walleadjt,wallwadjt,roadadjt)
deallocate(roofadjw,roadadjw)
deallocate(roofwater,roadwater)
deallocate(fnhwratio,fnsigmabld,fnindustryfg)
deallocate(fntrafficfg,fnbldheight,fnzo)
deallocate(vangle,hangle)

return
end subroutine tebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine sets defaults for aTEB state arrays

subroutine tebdefault(ifull,ta,tb,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq
real, dimension(1:ifull), intent(in) :: ta,tb
real, dimension(1:ifull) :: tc

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Setting aTEB state array defaults"

tc=tb
where (tc.lt.250.)
  tc=ta
end where

do iqu=1,ufull
  iq=ugrid(iqu)
  rooftemp(iqu,1)=ta(iq)
  walletemp(iqu,1)=ta(iq)
  wallwtemp(iqu,1)=ta(iq)
  roadtemp(iqu,1)=ta(iq)
  rooftemp(iqu,2)=0.3*ta(iq)+0.7*bldtemp
  walletemp(iqu,2)=0.3*ta(iq)+0.7*bldtemp
  wallwtemp(iqu,2)=0.3*ta(iq)+0.7*bldtemp
  roadtemp(iqu,2)=0.3*ta(iq)+0.7*tc(iq)
  rooftemp(iqu,3)=bldtemp
  walletemp(iqu,3)=bldtemp
  wallwtemp(iqu,3)=bldtemp
  roadtemp(iqu,3)=tc(iq)
end do

return
end subroutine tebdefault

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine tebload(ifull,rooftempi,walletempi,wallwtempi,roadtempi,roofwateri,roadwateri,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq
real, dimension(ifull,3), intent(in) :: rooftempi,walletempi,wallwtempi,roadtempi
real, dimension(ifull), intent(in) :: roofwateri,roadwateri

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Load aTEB state arrays"

do iqu=1,ufull
  iq=ugrid(iqu)
  rooftemp(iqu,:)=rooftempi(iq,:)
  walletemp(iqu,:)=walletempi(iq,:)
  wallwtemp(iqu,:)=wallwtempi(iq,:)
  roadtemp(iqu,:)=roadtempi(iq,:)
  roofwater(iqu)=roofwateri(iq)
  roadwater(iqu)=roadwateri(iq)
end do

if (any(rooftemp(:,:).lt.250.)) then
  write(6,*) "ERROR: Missing input data for tebload"
end if

return
end subroutine tebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine tebtype(ifull,itype,diag)

implicit none

integer, intent(in) :: ifull,diag
integer, dimension(ifull), intent(in) :: itype
integer iqu,iq,utype

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Load aTEB type arrays"

do iqu=1,ufull
  iq=ugrid(iqu)
  utype=itype(iq)
  
  select case(utype)
    case DEFAULT
      ! default urban (ecosystems 007,152,153,154 with differences in sigmau only)
      fnhwratio(iqu)=0.21
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=5.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=10.
      fnzo(iqu)=1.
    case(2)
      ! Dense urban (ecosystems 151)
      fnhwratio(iqu)=0.82
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=10.
      fntrafficfg(iqu)=20.
      fnbldheight(iqu)=30.
      fnzo(iqu)=3.
    case(3)
      ! Industries (ecosystems 155)
      fnhwratio(iqu)=0.4
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=20.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=20.
      fnzo(iqu)=2.
    case(4)
      ! Road and rail (ecosystems 156)
      fnhwratio(iqu)=0.05
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=30.
      fnbldheight(iqu)=5.
      fnzo(iqu)=0.5
    case(5)
      ! Ports (ecosystems 157)
      fnhwratio(iqu)=0.82
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=20.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=20.
      fnzo(iqu)=2.
    case(6)
      ! Airport (ecosystems 158)
      fnhwratio(iqu)=0.02
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=10.
      fnzo(iqu)=0.01
    case(7)
      ! Construction (ecosystems 159)
      fnhwratio(iqu)=0.005
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=0.
      fnbldheight(iqu)=5.
      fnzo(iqu)=0.1
    case(8)
      ! Urban parks (ecosystems 160)
      fnhwratio(iqu)=0.005
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=0.
      fnbldheight(iqu)=5.
      fnzo(iqu)=0.5
    case(9)
      ! Sport facilities (ecosystems 161)
      fnhwratio(iqu)=0.11
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=0.
      fnbldheight(iqu)=10.
      fnzo(iqu)=1.
  end select
 
end do

return
end subroutine tebtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines the various urban parameters at each
! grid point, instead of by type
!

subroutine tebfndef(ifull,hwratioi,sigmabldi,industryfgi,trafficfgi,bldheighti,zoi,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq
real, dimension(ifull), intent(in) :: hwratioi,sigmabldi,industryfgi,trafficfgi,bldheighti,zoi

if (ufull.eq.0) return

do iqu=1,ufull
  iq=ugrid(iqu)
  fnhwratio(iqu)=hwratioi(iq)
  fnsigmabld(iqu)=sigmabldi(iq)
  fnindustryfg(iqu)=industryfgi(iq)
  fntrafficfg(iqu)=trafficfgi(iq)
  fnbldheight(iqu)=bldheighti(iq)
  fnzo(iqu)=zoi(iq)
end do

return
end subroutine tebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine tebsave(ifull,rooftempo,walletempo,wallwtempo,roadtempo,roofwatero,roadwatero,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq
real, dimension(ifull,3), intent(out) :: rooftempo,walletempo,wallwtempo,roadtempo
real, dimension(ifull), intent(out) :: roofwatero,roadwatero

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Save aTEB state arrays"

do iq=1,ifull
  iqu=qgrid(iq)
  rooftempo(iq,:)=rooftemp(iqu,:)
  walletempo(iq,:)=walletemp(iqu,:)
  wallwtempo(iq,:)=wallwtemp(iqu,:)
  roadtempo(iq,:)=roadtemp(iqu,:)
  roofwatero(iq)=roofwater(iqu)
  roadwatero(iq)=roadwater(iqu)
end do

return
end subroutine tebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to roughness
! length.

subroutine tebzom(ifull,zo,zmin,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq
real, intent(in) :: zmin
real, dimension(ifull), intent(inout) :: zo
real, dimension(ifull), intent(in) :: sigmau
real zotemp

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Blend urban roughness length"

do iqu=1,ufull
  iq=ugrid(iqu)
  zotemp=(1.-sigmau(iq))/log(zmin/zo(iq))**2+sigmau(iq)/log(zmin/fnzo(iqu))**2
  zo(iq)=zmin*exp(-1./sqrt(zotemp))
end do

return
end subroutine tebzom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.

subroutine tebalb(ifull,alb,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iq,iqu
real, dimension(ifull), intent(out) :: alb
real wallesg,wallwsg,roadsg,wallpsi,roadpsi

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Calculate urban albedo"

alb=0.

do iqu=1,ufull
  iq=ugrid(iqu)
  call getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,fnhwratio(iqu),vangle(iqu),hangle(iqu))
  alb(iq)=1.-(fnsigmabld(iqu)*(1.-roofalpha)+(1.-fnsigmabld(iqu))* &
          (fnhwratio(iqu)*(wallesg+wallwsg)*(1.-wallalpha)+roadsg*(1.-roadalpha)))
end do

return
end subroutine tebalb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
!

subroutine tebnewangle(is,ifull,cosin,azimuthin,diag)

implicit none

integer, intent(in) :: is,ifull,diag
integer iqu,iq,ip
real, dimension(ifull), intent(in) :: cosin
real, dimension(ifull), intent(in) :: azimuthin

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Update solar zenith angle and azimuth angle"

do iqu=1,ufull
  iq=ugrid(iqu)
  ip=iq-is+1
  if ((ip.ge.1).and.(ip.le.ifull)) then
    hangle(iqu)=0.5*pi-azimuthin(ip)
    vangle(iqu)=acos(cosin(ip))
  end if
end do

return
end subroutine tebnewangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tebccangle(is,ifull,cosin,rlon,rlat,fjd,slag,dhr,dlt,diag)

implicit none

integer, intent(in) :: is,ifull,diag
real, intent(in) :: fjd,slag,dhr,dlt
real, dimension(ifull), intent(in) :: cosin,rlon,rlat
real, dimension(ifull) :: hloc,azimuth,x,y

if (ufull.eq.0) return


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! from zenith.f
hloc(:) = 2.*pi*fjd+slag+pi+rlon(:)+dhr*pi/24.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! estimate azimuth angle
x(:)=sin(-hloc(:))*cos(dlt)
y(:)=-cos(-hloc(:))*cos(dlt)*sin(rlat(:))+cos(rlat(:))*sin(dlt)
azimuth(:)=atan2(x(:),y(:))

call tebnewangle(is,ifull,cosin,azimuth,diag)

return
end subroutine tebccangle


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for calculating urban flux contrabution

! ifull = number of horizontal grid points
! dt = model time step
! zmin = first model level height (m)
! sg = incomming short wave radiation
! rg = incomming long wave radiation
! rnd = incomming rainfall rate (kg/(m^2 s))
! rho = atmospheric density at first model level
! temp = atmospheric temperature at first model level
! ps = surface pinvressure
! pa = pressure at first model level
! umag = horizontal wind speed at first model level
! ofg = Input/Output sensible heat flux
! oeg = Input/Output latient heat flux
! ots = Input/Output surface temperature
! owf = Input/Output wetness fraction

subroutine tebcalc(ifull,ofg,oeg,ots,owf,ddt,zmin,sg,rg,rnd,rho,temp,mixr,ps,pa,umag,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq,j,k,firstcall
real, intent(in) :: ddt,zmin
real, dimension(ifull), intent(in) :: sg,rg,rnd,rho,temp,mixr,ps,pa,umag,sigmau
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf
real, dimension(3) :: roofga,wallega,wallwga,roadga
real, dimension(3) :: roofdumtemp,walledumtemp,wallwdumtemp,roaddumtemp,roofdumt,walledumt,wallwdumt,roaddumt
real, dimension(3) :: rooforgt,walleorgt,wallworgt,roadorgt
real rooforgw,roadorgw
real roofdumwat,roaddumwat,roofdumw,roaddumw
real roofsg,roofrg,rooffg,roofeg
real wallesg,wallwsg,wallerg,wallwrg,wallefg,wallwfg
real wallpsi,roadpsi
real roadsg,roadrg,roadfg,roadeg
real topfg,topeg
real roofinvres,rwinvres,topinvres
real roofdelta,roaddelta
real oldcanyontemp,newcanyontemp
real cd,roofqsat,roadqsat
real canyontemp,canyonmix,canyonu
real ctmax,ctmin,evctmax,evct,oldevct,evctdum
real efftrafffg
real qsats,qsata
real tempc,mixrc,sigr,tempr,mixrr,dzmin
data firstcall/0/
save firstcall

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Evaluating aTEB"

do iqu=1,ufull
  iq=ugrid(iqu)

  ! calculate shortwave radiation
  call getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,fnhwratio(iqu),vangle(iqu),hangle(iqu))
  roofsg=(1.-roofalpha)*sg(iq)
  wallesg=(1.-wallalpha)*wallesg*sg(iq)
  wallwsg=(1.-wallalpha)*wallwsg*sg(iq)
  roadsg=(1.-roadalpha)*roadsg*sg(iq)

  ! canyon floor
  tempc=temp(iq)*(ps(iq)/pa(iq))**(rd/aircp)
  call getqsat(qsats,tempc,ps(iq))
  call getqsat(qsata,temp(iq),pa(iq))
  mixrc=mixr(iq)*qsats/qsata

  ! canyon roof (MJT suggestion)
  sigr=exp(-grav*fnbldheight(iqu)/(rd*temp(iq))) 
  tempr=temp(iq)*(ps(iq)*sigr/pa(iq))**(rd/aircp)
  call getqsat(qsats,tempr,ps(iq)*sigr)
  mixrr=mixr(iq)*qsats/qsata

  ! estimate canyonu
  if (zmin.gt.fnbldheight(iqu)) then ! above canyon
    canyonu=(2./pi)*umag(iq)*exp(-0.25*fnhwratio(iqu)) &
            *log(fnbldheight(iqu)/(3.*fnzo(iqu)))/log((zmin-fnbldheight(iqu)*2./3.)/fnzo(iqu))
  else ! inside canyon (MJT suggestion)
    canyonu=(2./pi)*umag(iq)*exp(-0.5*fnhwratio(iqu)*(zmin/fnbldheight(iqu)-0.5))
  end if

  ! scale traffic sensible heat flux for canyon          
  efftrafffg=fntrafficfg(iqu)/(1.-fnsigmabld(iqu))

  ! prep predictor-corrector arrays
  roofdumtemp(:)=rooftemp(iqu,:)
  walledumtemp(:)=walletemp(iqu,:)
  wallwdumtemp(:)=wallwtemp(iqu,:)
  roaddumtemp(:)=roadtemp(iqu,:)
  roofdumwat=roofwater(iqu)
  roaddumwat=roadwater(iqu)

  do j=1,2 ! corrector-predictor loop
    
    ! estimate water converage
    roofdelta=(roofdumwat/maxroofwater)**(2./3.)
    roaddelta=(roaddumwat/maxroadwater)**(2./3.)
  
    ! calculate long wave radiation (note additional wall compared to TEB scheme)
    roofrg=roofemiss*(rg(iq)-sbconst*roofdumtemp(1)**4)
    wallerg=wallemiss*(rg(iq)*(wallpsi+(1.-roademiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                      +sbconst*walledumtemp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
		              +sbconst*wallwdumtemp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-roademiss)*wallpsi*(1.-roadpsi)) &
                      +sbconst*roaddumtemp(1)**4*(wallpsi*roademiss+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
    wallwrg=wallemiss*(rg(iq)*(wallpsi+(1.-roademiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                      +sbconst*wallwdumtemp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
		              +sbconst*walledumtemp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-roademiss)*wallpsi*(1.-roadpsi)) &
                      +sbconst*roaddumtemp(1)**4*(wallpsi*roademiss+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
    roadrg=roademiss*(rg(iq)*(roadpsi+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*roaddumtemp(1)**4*(-1.+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*0.5*(walledumtemp(1)**4+wallwdumtemp(1)**4) &
                      *(wallemiss*(1.-roadpsi)+wallemiss*(1.-wallemiss)*(1.-roadpsi)*(1.-2.*wallpsi)))

    ! calculate mixing ratios
    call getqsat(roofqsat,roofdumtemp(1),ps(iq)*sigr)
    call getqsat(roadqsat,roaddumtemp(1),ps(iq))

    ! calculate aerodynamic resistances 
    ! (should really use two model atmospheric levels here.  One in the canyon and one above roofs.
    ! However, for this to work the sensible heat flux must also be fed back into both atmospheric
    ! levels.  Since the host model only allows the sensible heat flux to be fed into the lowest model level,
    ! then we also estimate the interaction with the roof using the lowest model level)
    dzmin=max(abs(zmin-fnbldheight(iqu)),zoroof+1.)
    call getinvres(roofinvres,cd,zoroof,dzmin,roofdumtemp(1),tempr,umag(iq)) 
    rooffg=aircp*rho(iq)*(roofdumtemp(1)-tempr)*roofinvres

    ! diagnose canyon air temperature
    ctmax=max(tempc,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1))+1. ! max canyon temp
    ctmin=min(tempc,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1))    ! min canyon temp
    call solvecanyon(evctmax,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,fnzo(iqu) &
      ,zmin,ctmax,tempc,umag(iq),canyonu,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1) &
      ,rho(iq),efftrafffg,fnhwratio(iqu))
    canyontemp=0.5*(ctmax+ctmin)
    do k=1,2 ! bisect
      call solvecanyon(evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,fnzo(iqu) &
        ,zmin,canyontemp,tempc,umag(iq),canyonu,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1) &
        ,rho(iq),efftrafffg,fnhwratio(iqu))
      if ((evct*evctmax).lt.0.) then
        ctmin=canyontemp
      else
        ctmax=canyontemp
        evctmax=evct
      end if
      oldcanyontemp=canyontemp
      canyontemp=0.5*(ctmax+ctmin)
    end do

    do k=1,5 ! sectant
      oldevct=evct
      call solvecanyon(evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,fnzo(iqu) &
        ,zmin,canyontemp,tempc,umag(iq),canyonu,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1) &
        ,rho(iq),efftrafffg,fnhwratio(iqu))
      evctdum=evct-oldevct
      if (evctdum.eq.0.) exit    
      newcanyontemp=canyontemp-evct*(canyontemp-oldcanyontemp)/evctdum
      oldcanyontemp=canyontemp
      canyontemp=newcanyontemp
    end do
    canyontemp=min(max(canyontemp,ctmin),ctmax)    
      
    if (roofqsat.lt.mixrr) roofdelta=1. 
    if (roadqsat.lt.mixrc) roaddelta=1. ! equivilent to roadqsat.lt.canyonmix
    canyonmix=(roaddelta*roadqsat*rwinvres+mixrc*topinvres)/(roaddelta*rwinvres+topinvres)
 
    ! calculate latent heat flux
    roofeg=lv*min(rho(iq)*roofdelta*(roofqsat-mixrr)*roofinvres,roofdumwat+rnd(iq))
    roadeg=lv*min(rho(iq)*roaddelta*(roadqsat-canyonmix)*rwinvres,roaddumwat+rnd(iq))
    topeg=roadeg
  
    ! calculate condution heat flux (e.g., through walls)
    do k=1,2
      roofga(k)=2.*(roofdumtemp(k)-roofdumtemp(k+1))/(roofdepth(k)/rooflambda(k)+roofdepth(k+1)/rooflambda(k+1))
      wallega(k)=2.*(walledumtemp(k)-walledumtemp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
      wallwga(k)=2.*(wallwdumtemp(k)-wallwdumtemp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
      roadga(k)=2.*(roaddumtemp(k)-roaddumtemp(k+1))/(roaddepth(k)/roadlambda(k)+roaddepth(k+1)/roadlambda(k+1))
    end do
    roofga(3)=2.*rooflambda(3)*(roofdumtemp(3)-bldtemp)/(roofdepth(3))
    wallega(3)=2.*walllambda(3)*(walledumtemp(3)-bldtemp)/(walldepth(3))
    wallwga(3)=2.*walllambda(3)*(wallwdumtemp(3)-bldtemp)/(walldepth(3))
    roadga(3)=0.
  
    ! update urban temperatures
    roofdumt(1)=(roofsg+roofrg-rooffg-roofeg-roofga(1))/(roofcp(1)*roofdepth(1))
    walledumt(1)=(wallesg+wallerg-wallefg-wallega(1))/(wallcp(1)*walldepth(1))
    wallwdumt(1)=(wallwsg+wallwrg-wallwfg-wallwga(1))/(wallcp(1)*walldepth(1))
    roaddumt(1)=(roadsg+roadrg-roadfg-roadeg-roadga(1))/(roadcp(1)*roaddepth(1))
    do k=2,3
      roofdumt(k)=(roofga(k-1)-roofga(k))/(roofcp(k)*roofdepth(k))
      walledumt(k)=(wallega(k-1)-wallega(k))/(wallcp(k)*walldepth(k))
      wallwdumt(k)=(wallwga(k-1)-wallwga(k))/(wallcp(k)*walldepth(k))
      roaddumt(k)=(roadga(k-1)-roadga(k))/(roadcp(k)*roaddepth(k))
    end do
    roofdumw=(rnd(iq)-roofeg/lv)
    roaddumw=(rnd(iq)-roadeg/lv)

    ! predictor-corrector scheme
    if (j.eq.1) then ! predictor
      if (firstcall.eq.0) then
        roofadjt(iqu,:)=roofdumt(:)
        walleadjt(iqu,:)=walledumt(:)
        wallwadjt(iqu,:)=wallwdumt(:)
        roadadjt(iqu,:)=roaddumt(:)
        roofadjw(iqu)=roofdumw
        roadadjw(iqu)=roaddumw
      end if
      roofdumtemp(:)=rooftemp(iqu,:)+ddt*(1.5*roofdumt(:)-0.5*roofadjt(iqu,:))
      walledumtemp(:)=walletemp(iqu,:)+ddt*(1.5*walledumt(:)-0.5*walleadjt(iqu,:))
      wallwdumtemp(:)=wallwtemp(iqu,:)+ddt*(1.5*wallwdumt(:)-0.5*wallwadjt(iqu,:))
      roaddumtemp(:)=roadtemp(iqu,:)+ddt*(1.5*roaddumt(:)-0.5*roadadjt(iqu,:))
      rooforgt(:)=roofdumt(:)
      walleorgt(:)=walledumt(:)
      wallworgt(:)=wallwdumt(:)
      roadorgt(:)=roaddumt(:)
      roofdumwat=roofwater(iqu)+ddt*(1.5*roofdumw-0.5*roofadjw(iqu))
      roaddumwat=roadwater(iqu)+ddt*(1.5*roaddumw-0.5*roadadjw(iqu))
      rooforgw=roofdumw
      roadorgw=roadorgw
        
    else ! corrector
 
      roofdumtemp(:)=rooftemp(iqu,:)+(ddt/12.)*(5.*roofdumt(:)+8.*rooforgt(:)-roofadjt(iqu,:))
      walledumtemp(:)=walletemp(iqu,:)+(ddt/12.)*(5.*walledumt(:)+8.*walleorgt(:)-walleadjt(iqu,:))
      wallwdumtemp(:)=wallwtemp(iqu,:)+(ddt/12.)*(5.*wallwdumt(:)+8.*wallworgt(:)-wallwadjt(iqu,:))
      roaddumtemp(:)=roadtemp(iqu,:)+(ddt/12.)*(5.*roaddumt(:)+8.*roadorgt(:)-roadadjt(iqu,:))
      roofadjt(iqu,:)=rooforgt(:)
      walleadjt(iqu,:)=walleorgt(:)
      wallwadjt(iqu,:)=wallworgt(:)
      roadadjt(iqu,:)=roadorgt(:)
      roofdumwat=roofwater(iqu)+(ddt/12.)*(5.*roofdumw+8.*rooforgw-roofadjw(iqu))
      roaddumwat=roadwater(iqu)+(ddt/12.)*(5.*roaddumw+8.*roadorgw-roadadjw(iqu))
      roofadjw(iqu)=rooforgw
      roadadjw(iqu)=roadorgw
  
    end if
     
    roofdumtemp(:)=min(max(roofdumtemp(:),250.),350.)
    walledumtemp(:)=min(max(walledumtemp(:),250.),350.)
    wallwdumtemp(:)=min(max(wallwdumtemp(:),250.),350.)
    roaddumtemp(:)=min(max(roaddumtemp(:),250.),350.)
    roofdumwat=min(max(roofdumwat,0.),maxroofwater)
    roaddumwat=min(max(roaddumwat,0.),maxroadwater)
  end do
    
  rooftemp(iqu,:)=roofdumtemp(:)
  walletemp(iqu,:)=walledumtemp(:)
  wallwtemp(iqu,:)=wallwdumtemp(:)
  roadtemp(iqu,:)=roaddumtemp(:)
  roofwater(iqu)=roofdumwat
  roadwater(iqu)=roaddumwat

  ! calculate outputs
  ofg(iq)=(1.-sigmau(iq))*ofg(iq)+sigmau(iq)*(fnsigmabld(iqu)*rooffg+(1.-fnsigmabld(iqu))*topfg+fnindustryfg(iqu))
  oeg(iq)=(1.-sigmau(iq))*oeg(iq)+sigmau(iq)*(fnsigmabld(iqu)*roofeg+(1.-fnsigmabld(iqu))*topeg)
  ots(iq)=(1.-sigmau(iq))*ots(iq)+sigmau(iq)*(fnsigmabld(iqu)*rooftemp(iqu,1)+(1.-fnsigmabld(iqu))*canyontemp) !MJT - since this is what the atmosphere can 'see'
  owf(iq)=(1.-sigmau(iq))*owf(iq)+sigmau(iq)*(fnsigmabld(iqu)*roofdelta+(1.-fnsigmabld(iqu))*roaddelta)

end do
  
firstcall=1

return
end subroutine tebcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getqsat(qsat,temp,ps)

implicit none

real, intent(in) :: temp,ps
real, intent(out) :: qsat
real esatf,tdiff

tdiff=min(max( temp-123.16, 0.), 219.)
esatf=(1.-(tdiff-aint(tdiff)))*table(int(tdiff))+ (tdiff-aint(tdiff))*table(int(tdiff)+1)
qsat=.622*esatf/(ps-esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Based on Eva's sflux.f (i.e., vegetation)
! It is likely that these relations are modified for urban land-use
! (e.g., zot/zo is possibly an order of magnitude greater for urban than vegetation)

subroutine getinvres(invres,cd,zo,zmin,stemp,theta,umag)

implicit none

real, intent(in) :: zo,zmin,stemp,theta,umag
real, intent(out) :: invres,cd
real zolog,af,aft,xx,ri,fm,fh,factch,root,denma,denha,vmag
real rimax
real, parameter :: bprm = 5.
real, parameter :: chs=2.6
real, parameter :: cms=5.
real, parameter :: fmroot=0.57735
real, parameter :: vkar=0.4

rimax=(1./fmroot-1.)/bprm
factch=sqrt(7.4)

zolog=log(zmin/zo)
vmag=max(umag,0.2)

af=(vkar/zolog)**2
aft=vkar**2/(zolog*(2.+zolog))

xx=grav*zmin*(1.-stemp/theta)
ri=min(xx/vmag**2,rimax)

if (ri>0.) then
  fm=1./(1.+bprm*ri)**2
  fh=fm
else
  root=sqrt(-ri*zmin/zo)
  denma=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm *ri/denma
  denha=1.+chs*2.*bprm*factch*aft*root
  fh=1.-2.*bprm *ri/denha
endif

cd=af*fm
invres=aft*fh*umag

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate shortwave radiation coefficents (modified to include 2nd wall)

subroutine getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,hwratio,vangle,hangle)

implicit none

real, intent(in) :: hwratio,vangle,hangle
real, intent(out) :: wallesg,wallwsg,roadsg,wallpsi,roadpsi
real thetazero,walles,wallws,roads,ta,tc,tz,xa,xb,ya,yb

wallpsi=0.5*(hwratio+1.-sqrt(hwratio**2+1.))/hwratio
roadpsi=sqrt(hwratio**2+1.)-hwratio

! integrate through 180 instead of 360
if (vangle.ge.(0.5*pi)) then
  walles=0.
  wallws=1./hwratio
  roads=0.
else
  ta=tan(vangle)
  thetazero=asin(1./max(hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  tz=2.*thetazero
  xa=max(hangle-thetazero,0.)-max(hangle-pi+thetazero,0.)+pi-min(hangle+pi+thetazero,pi)
  xb=pi-tz-xa
  ya=2.-cos(min(thetazero,max(abs(hangle),0.)))-cos(min(thetazero,max(hangle-pi+thetazero,0.)))
  yb=tc-ya
  
  ! note that these terms now include the azimuth angle
  walles=((1./hwratio)*xa+ta*ya)/pi
  wallws=((1./hwratio)*xb+ta*yb)/pi
  roads=(tz-hwratio*ta*tc)/pi
end if

! note that these terms are truncated to 2nd order reflections, compared to TEB which uses infinte reflections.
wallesg=walles+roadalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*wallws+(wallalpha*(1.-2.*wallpsi))**2*walles &
        +roadalpha*wallalpha*wallpsi*(1.-roadpsi)*wallws+roadalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
wallwsg=wallws+roadalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*walles+(wallalpha*(1.-2.*wallpsi))**2*wallws &
        +roadalpha*wallalpha*wallpsi*(1.-roadpsi)*walles+roadalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
roadsg=roads+wallalpha*(1.-roadpsi)*0.5*(walles+wallws)+wallalpha*roadalpha*wallpsi*(1.-roadpsi)*roads &
       +wallalpha**2*(1.-roadpsi)*(1.-2.*wallpsi)*0.5*(walles+wallws)

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solvecanyon(evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,zo,zmin,ctemp &
  ,theta,umag,cu,walletemp,wallwtemp,roadtemp,rho,efftrafficfg,hwratio)

implicit none

real, intent(in) :: zo,zmin,ctemp,theta,umag,cu,walletemp,wallwtemp,roadtemp,rho,efftrafficfg,hwratio
real, intent(out) :: evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres
real cd,cw

call getinvres(topinvres,cd,zo,zmin,ctemp,theta,umag) 
cw=sqrt(cd)*umag ! diagnose canyonw
rwinvres=(11.8+4.2*sqrt(cu**2+cw**2))/(aircp*rho) ! From Rowley, et al (1930)
wallefg=aircp*rho*(walletemp-ctemp)*rwinvres
wallwfg=aircp*rho*(wallwtemp-ctemp)*rwinvres
roadfg=aircp*rho*(roadtemp-ctemp)*rwinvres
topfg=aircp*rho*(ctemp-theta)*topinvres
evct=topfg-(roadfg+hwratio*(wallefg+wallwfg)+efftrafficfg)

return
end subroutine solvecanyon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ateb
