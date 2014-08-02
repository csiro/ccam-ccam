module indices_m

implicit none

private
public iw_g,is_g,ise_g,ie_g,in_g,iwn_g,inw_g,isw_g,ies_g,iws_g
public ine_g,ien_g,inn_g,iss_g,iww_g,iee_g,iwu2_g,isv2_g
public lwws_g,lwss_g,lees_g,less_g,lwwn_g,lwnn_g,leen_g,lenn_g,lsww_g
public lssw_g,lsee_g,lsse_g,lnww_g,lnnw_g,lnee_g,lnne_g
public npann_g,npane_g,npanw_g,npans_g
public iw,is,ise,ie,ine,in,iwn,inw,isw,ies,iws,ien,inn,iss,iww,iee,iwu,isv
public ieu,inv,iwwu,issv,ieeu,innv
public iev,iwv,inu,isu
public lwws,lwss,lees,less,lwwn,lwnn,leen,lenn,lsww
public lssw,lsee,lsse,lnww,lnnw,lnee,lnne
public indices_init,indices_end

integer, dimension(:), allocatable, save :: iw,is,ise,ie,ine,in,iwn,inw,isw,ies,iws,ien,inn,iss,iww,iee,iwu,isv
integer, dimension(:), allocatable, save :: ieu,inv,iwwu,issv,ieeu,innv
integer, dimension(:), allocatable, save :: iev,iwv,inu,isu
integer, dimension(:), allocatable, save :: lwws,lwss,lees,less,lwwn,lwnn,leen,lenn,lsww
integer, dimension(:), allocatable, save :: lssw,lsee,lsse,lnww,lnnw,lnee,lnne

integer, parameter, dimension(0:5) :: npann_g = (/ 1, 103, 3, 105, 5, 101 /)
integer, parameter, dimension(0:5) :: npans_g = (/ 104, 0, 100, 2, 102, 4 /)
integer, parameter, dimension(0:5) :: npane_g = (/ 102, 2, 104, 4, 100, 0 /)
integer, parameter, dimension(0:5) :: npanw_g = (/ 5, 105, 1, 101, 3, 103 /)

contains

subroutine indices_init(ifull_g,ifull,iextra,npanels,npan)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra,npanels,npan

allocate(iw(ifull),is(ifull),ise(ifull))
allocate(ie(ifull),ine(ifull),in(ifull),iwn(ifull))
allocate(inw(ifull),isw(ifull),ies(ifull),iws(ifull))
allocate(ien(ifull),inn(ifull),iss(ifull),iww(ifull))
allocate(iee(ifull),iwu(ifull),isv(ifull),ieu(ifull))
allocate(inv(ifull),iwwu(ifull),issv(ifull),ieeu(ifull))
allocate(innv(ifull))
allocate(iev(ifull),iwv(ifull),inu(ifull),isu(ifull))
allocate(lwws(npan),lwss(npan))
allocate(lees(npan),less(npan),lwwn(npan),lwnn(npan))
allocate(leen(npan),lenn(npan),lsww(npan))
allocate(lssw(npan),lsee(npan),lsse(npan),lnww(npan))
allocate(lnnw(npan),lnee(npan),lnne(npan))

return
end subroutine indices_init

subroutine indices_end

implicit none

deallocate(iw,is,ise,ie,ine,in,iwn,inw,isw,ies,iws,ien,inn,iss,iww,iee,iwu,isv)
deallocate(ieu,inv,iwwu,issv,ieeu,innv)
deallocate(iev,iwv,inu,isu)
deallocate(lwws,lwss,lees,less,lwwn,lwnn,leen,lenn,lsww)
deallocate(lssw,lsee,lsse,lnww,lnnw,lnee,lnne)

return
end subroutine indices_end

function in_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (j==il_g) then
  if (npann_g(n)<100) then
    iqq=i+npann_g(n)*il_g*il_g
  else
    iqq=1+(il_g-i)*il_g+(npann_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq+il_g    
end if
end function in_g

function is_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (j==1) then
  if (npans_g(n)<100) then
    iqq=i+(il_g-1)*il_g+npans_g(n)*il_g*il_g
  else
    iqq=il_g+(il_g-i)*il_g+(npans_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq-il_g    
end if
end function is_g

function ie_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (i==il_g) then
  if (npane_g(n)<100) then
    iqq=1+(j-1)*il_g+npane_g(n)*il_g*il_g
  else
    iqq=il_g+1-j+(npane_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq+1
end if
end function ie_g

function iw_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (i==1) then
  if (npanw_g(n)<100) then
    iqq=il_g+(j-1)*il_g+npanw_g(n)*il_g*il_g
  else
    iqq=il_g+1-j+(il_g-1)*il_g+(npanw_g(n)-100)*il_g*il_g
  end if
else
  iqq=iq-1
end if
end function iw_g

function inn_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npann_g(n)>=100.and.j==il_g) then
  iqq=ie_g(in_g(iq))
else
  iqq=in_g(in_g(iq))
end if
end function inn_g

function iss_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=iw_g(is_g(iq))
else
  iqq=is_g(is_g(iq))
end if
end function iss_g

function iee_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npane_g(n)>=100.and.i==il_g) then
  iqq=in_g(ie_g(iq))
else
  iqq=ie_g(ie_g(iq))
end if
end function iee_g

function iww_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=is_g(iw_g(iq))
else
  iqq=iw_g(iw_g(iq))
end if
end function iww_g

function ine_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npane_g(n)>=100.and.i==il_g) then
  iqq=iw_g(ie_g(iq))
else
  iqq=in_g(ie_g(iq))
end if
end function ine_g

function ise_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npane_g(n)>=100.and.i==il_g) then
  iqq=ie_g(ie_g(iq))
else
  iqq=is_g(ie_g(iq))
end if
end function ise_g

function ien_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npann_g(n)>=100.and.j==il_g) then
  iqq=is_g(in_g(iq))
else
  iqq=ie_g(in_g(iq))
end if
end function ien_g

function iwn_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npann_g(n)>=100.and.j==il_g) then
  iqq=in_g(in_g(iq))
else
  iqq=iw_g(in_g(iq))
end if
end function iwn_g

function inw_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=iw_g(iw_g(iq))
else
  iqq=in_g(iw_g(iq))
end if
end function inw_g

function isw_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=ie_g(iw_g(iq))
else
  iqq=is_g(iw_g(iq))
end if
end function isw_g

function ies_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=is_g(is_g(iq))
else
  iqq=ie_g(is_g(iq))
end if
end function ies_g

function iws_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=in_g(is_g(iq))
else
  iqq=iw_g(is_g(iq))
end if
end function iws_g

function iwu2_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npanw_g(n)>=100.and.i==1) then
  iqq=iw_g(iq)+ifull_g
else
  iqq=iw_g(iq)
end if
end function iwu2_g

function isv2_g(iq) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: iq
integer iqq
integer n,i,j
n=(iq-1)/(il_g*il_g)
j=(iq-1-n*il_g*il_g)/il_g+1
i=iq-(j-1)*il_g-n*il_g*il_g
if (npans_g(n)>=100.and.j==1) then
  iqq=is_g(iq)-ifull_g
else
  iqq=is_g(iq)
end if
end function isv2_g

function leen_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=iss_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=iee_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function leen_g

function lenn_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=ise_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=ien_g(in_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lenn_g

function lwnn_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=ine_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=iwn_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lwnn_g

function lsee_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iwn_g(ie_g(il_g+n*il_g*il_g))
else
  iqq=ise_g(ie_g(il_g+n*il_g*il_g))
end if
end function lsee_g

function lnee_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iwn_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=ine_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnee_g

function lnne_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iww_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=inn_g(ie_g(il_g+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnne_g

function lsww_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=ies_g(iw_g(1+n*il_g*il_g))
else
  iqq=isw_g(iw_g(1+n*il_g*il_g))
end if
end function lsww_g

function lssw_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=iee_g(iw_g(1+n*il_g*il_g))
else
  iqq=iss_g(iw_g(1+n*il_g*il_g))
end if
end function lssw_g

function lnww_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=iws_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=inw_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnww_g

function lwws_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=inn_g(is_g(1+n*il_g*il_g))
else
  iqq=iww_g(is_g(1+n*il_g*il_g))
end if
end function lwws_g

function lwss_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=inw_g(is_g(1+n*il_g*il_g))
else
  iqq=iws_g(is_g(1+n*il_g*il_g))
end if
end function lwss_g

function less_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=isw_g(is_g(il_g+n*il_g*il_g))
else
  iqq=ies_g(is_g(il_g+n*il_g*il_g))
end if
end function less_g

function lwwn_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npann_g(n)>=100) then
  iqq=inn_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=iww_g(in_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lwwn_g

function lsse_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npane_g(n)>=100) then
  iqq=iee_g(ie_g(il_g+n*il_g*il_g))
else
  iqq=iss_g(ie_g(il_g+n*il_g*il_g))
end if
end function lsse_g

function lnnw_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npanw_g(n)>=100) then
  iqq=iww_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
else
  iqq=inn_g(iw_g(1+(il_g-1)*il_g+n*il_g*il_g))
end if
end function lnnw_g

function lees_g(n) result(iqq)
implicit none
include 'newmpar.h'
integer, intent(in) :: n
integer iqq
if (npans_g(n)>=100) then
  iqq=iss_g(is_g(il_g+n*il_g*il_g))
else
  iqq=iee_g(is_g(il_g+n*il_g*il_g))
end if
end function lees_g

end module indices_m
