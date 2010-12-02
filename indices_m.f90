module indices_m

implicit none

private
public iw_g,isw_g,is_g,ise_g,ie_g,ine_g,in_g,iwn_g,ien_g,inn_g,iss_g,iww_g,iee_g,iwu_g,isv_g
public iwu2_g,isv2_g,ieu2_g,inv2_g,iev2_g,inu2_g,ieu_g,inv_g,iwwu2_g,issv2_g,ieeu2_g,innv2_g
public lwws_g,lws_g,lwss_g,les_g,lees_g,less_g,lwwn_g,lwnn_g,leen_g,lenn_g,lsww_g
public lsw_g,lssw_g,lsee_g,lsse_g,lnww_g,lnw_g,lnnw_g,lnee_g,lnne_g
public npann_g,npane_g,npanw_g,npans_g
public iw,isw,is,ise,ie,ine,in,iwn,ien,inn,iss,iww,iee,iwu,isv
public ieu,inv,iwwu,issv,ieeu,innv
public lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen,lenn,lsww
public lsw,lssw,lsee,lsse,lnww,lnw,lnnw,lnee,lnne
public npann,npane,npanw,npans
public indices_init,indices_end

integer, dimension(:), allocatable, save :: iw_g,isw_g,is_g,ise_g,ie_g,ine_g,in_g,iwn_g,ien_g,inn_g,iss_g,iww_g,iee_g,iwu_g,isv_g
integer, dimension(:), allocatable, save :: iwu2_g,isv2_g,ieu2_g,inv2_g,iev2_g,inu2_g,ieu_g,inv_g,iwwu2_g,issv2_g,ieeu2_g,innv2_g
integer, dimension(:), allocatable, save :: lwws_g,lws_g,lwss_g,les_g,lees_g,less_g,lwwn_g,lwnn_g,leen_g,lenn_g,lsww_g
integer, dimension(:), allocatable, save :: lsw_g,lssw_g,lsee_g,lsse_g,lnww_g,lnw_g,lnnw_g,lnee_g,lnne_g
integer, dimension(:), allocatable, save :: npann_g,npane_g,npanw_g,npans_g
integer, dimension(:), allocatable, save :: iw,isw,is,ise,ie,ine,in,iwn,ien,inn,iss,iww,iee,iwu,isv
integer, dimension(:), allocatable, save :: ieu,inv,iwwu,issv,ieeu,innv
integer, dimension(:), allocatable, save :: lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen,lenn,lsww
integer, dimension(:), allocatable, save :: lsw,lssw,lsee,lsse,lnww,lnw,lnnw,lnee,lnne
integer, dimension(:), allocatable, save :: npann,npane,npanw,npans

contains

subroutine indices_init(ifull_g,ifull,iextra,npanels,npan)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra,npanels,npan

allocate(iw_g(ifull_g),isw_g(ifull_g),is_g(ifull_g),ise_g(ifull_g))
allocate(ie_g(ifull_g),ine_g(ifull_g),in_g(ifull_g),iwn_g(ifull_g))
allocate(ien_g(ifull_g),inn_g(ifull_g),iss_g(ifull_g),iww_g(ifull_g))
allocate(iee_g(ifull_g),iwu_g(ifull_g),isv_g(ifull_g),iwu2_g(ifull_g))
allocate(isv2_g(ifull_g),ieu2_g(ifull_g),inv2_g(ifull_g),iev2_g(ifull_g))
allocate(inu2_g(ifull_g),ieu_g(ifull_g),inv_g(ifull_g),iwwu2_g(ifull_g))
allocate(issv2_g(ifull_g),ieeu2_g(ifull_g),innv2_g(ifull_g))
allocate(lwws_g(0:npanels),lws_g(0:npanels),lwss_g(0:npanels),les_g(0:npanels))
allocate(lees_g(0:npanels),less_g(0:npanels),lwwn_g(0:npanels),lwnn_g(0:npanels))
allocate(leen_g(0:npanels),lenn_g(0:npanels),lsww_g(0:npanels),lsw_g(0:npanels))
allocate(lssw_g(0:npanels),lsee_g(0:npanels),lsse_g(0:npanels),lnww_g(0:npanels))
allocate(lnw_g(0:npanels),lnnw_g(0:npanels),lnee_g(0:npanels),lnne_g(0:npanels))
allocate(npann_g(0:13),npane_g(0:13),npanw_g(0:13),npans_g(0:13))
allocate(iw(ifull),isw(ifull),is(ifull),ise(ifull))
allocate(ie(ifull),ine(ifull),in(ifull),iwn(ifull))
allocate(ien(ifull),inn(ifull),iss(ifull),iww(ifull))
allocate(iee(ifull),iwu(ifull),isv(ifull),ieu(ifull))
allocate(inv(ifull),iwwu(ifull),issv(ifull),ieeu(ifull))
allocate(innv(ifull))
allocate(lwws(npan),lws(npan),lwss(npan),les(npan))
allocate(lees(npan),less(npan),lwwn(npan),lwnn(npan))
allocate(leen(npan),lenn(npan),lsww(npan),lsw(npan))
allocate(lssw(npan),lsee(npan),lsse(npan),lnww(npan))
allocate(lnw(npan),lnnw(npan),lnee(npan),lnne(npan))
allocate(npann(0:13),npane(0:13),npanw(0:13),npans(0:13))

npann_g=(/  1, 2,107,  4,106,  6,  7,109,  9,112, 11, 12,102,101/)
npane_g=(/103, 3,  4,105,  5,110,108,  8, 10, 11,100,113, 13,  0/)
npanw_g=(/13,113,112,  1,  2,  4,104,102,  7,107,  8,  9,109, 12/)
npans_g=(/110, 0,  1,100,  3,103,  5,  6,106,  8,105, 10, 11,111/)   

return
end subroutine indices_init

subroutine indices_end

implicit none

if (allocated(iw_g)) deallocate(iw_g)
if (allocated(isw_g)) deallocate(isw_g)
if (allocated(is_g)) deallocate(is_g)
if (allocated(ise_g)) deallocate(ise_g)
if (allocated(ie_g)) deallocate(ie_g)
if (allocated(ine_g)) deallocate(ine_g)
if (allocated(in_g)) deallocate(in_g)
if (allocated(iwn_g)) deallocate(iwn_g)
if (allocated(ien_g)) deallocate(ien_g)
if (allocated(inn_g)) deallocate(inn_g)
if (allocated(iss_g)) deallocate(iss_g)
if (allocated(iww_g)) deallocate(iww_g)
if (allocated(iee_g)) deallocate(iee_g)
if (allocated(iwu_g)) deallocate(iwu_g)
if (allocated(isv_g)) deallocate(isv_g)
if (allocated(iwu2_g)) deallocate(iwu2_g)
if (allocated(isv2_g)) deallocate(isv2_g)
if (allocated(ieu2_g)) deallocate(ieu2_g)
if (allocated(inv2_g)) deallocate(inv2_g)
if (allocated(iev2_g)) deallocate(iev2_g)
if (allocated(inu2_g)) deallocate(inu2_g)
if (allocated(ieu_g)) deallocate(ieu_g)
if (allocated(inv_g)) deallocate(inv_g)
if (allocated(iwwu2_g)) deallocate(iwwu2_g)
if (allocated(issv2_g)) deallocate(issv2_g)
if (allocated(ieeu2_g)) deallocate(ieeu2_g)
if (allocated(innv2_g)) deallocate(innv2_g)
if (allocated(lwws_g)) deallocate(lws_g)
if (allocated(lws_g)) deallocate(lws_g)
if (allocated(lwss_g)) deallocate(lwss_g)
if (allocated(les_g)) deallocate(les_g)
if (allocated(lees_g)) deallocate(lees_g)
if (allocated(less_g)) deallocate(less_g)
if (allocated(lwwn_g)) deallocate(lwwn_g)
if (allocated(lwnn_g)) deallocate(lwnn_g)
if (allocated(leen_g)) deallocate(leen_g)
if (allocated(lenn_g)) deallocate(lenn_g)
if (allocated(lsww_g)) deallocate(lsww_g)
if (allocated(lsw_g)) deallocate(lsw_g)
if (allocated(lssw_g)) deallocate(lssw_g)
if (allocated(lsee_g)) deallocate(lsee_g)
if (allocated(lsse_g)) deallocate(lsse_g)
if (allocated(lnww_g)) deallocate(lnww_g)
if (allocated(lnw_g)) deallocate(lnw_g)
if (allocated(lnnw_g)) deallocate(lnnw_g)
if (allocated(lnee_g)) deallocate(lnee_g)
if (allocated(lnne_g)) deallocate(lnne_g)
if (allocated(npann_g)) deallocate(npann_g)
if (allocated(npane_g)) deallocate(npane_g)
if (allocated(npanw_g)) deallocate(npanw_g)
if (allocated(npans_g)) deallocate(npans_g)
deallocate(iw,isw,is,ise,ie,ine,in,iwn,ien,inn,iss,iww,iee,iwu,isv)
deallocate(ieu,inv,iwwu,issv,ieeu,innv)
deallocate(lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen,lenn,lsww)
deallocate(lsw,lssw,lsee,lsse,lnww,lnw,lnnw,lnee,lnne)
deallocate(npann,npane,npanw,npans)

return
end subroutine indices_end

end module indices_m
