c     second line are used in interpolation routines
      integer iw,isw,is,ise,ie,ine,in,iwn,ien,inn,iss,iww,iee,iwu,isv,
     &        iwu2,isv2,ieu2,inv2,iev2,inu2,ieu,inv,
     &        lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen,lenn,lsww,  
     &        lsw ,lssw,lsee,lsse,lnww,lnw,lnnw,lnee,lnne,  
     &        npann,npane,npanw,npans  

      common/indices/iw(ifull),isw(ifull),is(ifull),ise(ifull),
     &               ie(ifull),ine(ifull),in(ifull),iwn(ifull),
     &               ien(ifull),
     &               inn(ifull),iss(ifull),iww(ifull),iee(ifull),
     &               iwu(ifull),isv(ifull),    ! these for sflux, vertmix
     &               iwu2(ifull),isv2(ifull),ieu2(ifull),inv2(ifull) ! div calcs
     &               ,iev2(ifull),inu2(ifull)  ! upglobal
     &               ,ieu(ifull),inv(ifull)    ! staguv3
     &               ,lwws(0:npanels),lws (0:npanels),lwss(0:npanels)  ! ints
     &               ,les (0:npanels),lees(0:npanels),less(0:npanels)  ! ints
     &               ,lwwn(0:npanels),lwnn(0:npanels),leen(0:npanels)  ! ints
     &               ,lenn(0:npanels)                ,lsww(0:npanels)  ! ints
     &               ,lsw (0:npanels),lssw(0:npanels),lsee(0:npanels)  ! ints
     &               ,lsse(0:npanels),lnww(0:npanels),lnw (0:npanels)  ! ints
     &               ,lnnw(0:npanels),lnee(0:npanels),lnne(0:npanels)  ! ints
     &               ,npann(0:13),npane(0:13),npanw(0:13),npans(0:13)  ! hordifg
