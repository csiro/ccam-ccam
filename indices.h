!     second line are used in interpolation routines
      integer iw,isw,is,ise,ie,ine,in,iwn,ien,inn,iss,iww,iee,iwu,isv,     &
     &        ieu,inv,iwwu,issv,ieeu,innv,                                 &
     &        lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen,lenn,lsww,        &
     &        lsw ,lssw,lsee,lsse,lnww,lnw,lnnw,lnee,lnne,                 &
     &        npann,npane,npanw,npans  

      common/indices/iw(ifull),isw(ifull),is(ifull),ise(ifull),            &
     &               ie(ifull),ine(ifull),in(ifull),iwn(ifull),            &
     &               ien(ifull),                                           &
     &               inn(ifull),iss(ifull),iww(ifull),iee(ifull),          &
     &               iwu(ifull),isv(ifull),                                & ! these for sflux, vertmix
     &               ieu(ifull),inv(ifull),                                & ! staguv3
     &               iwwu(ifull),issv(ifull),ieeu(ifull),innv(ifull),      & 
     &               lwws(npan),lws (npan),lwss(npan),                     & ! ints
     &               les (npan),lees(npan),less(npan),                     & ! ints
     &               lwwn(npan),lwnn(npan),leen(npan),                     & ! ints
     &               lenn(npan)                ,lsww(npan),                & ! ints
     &               lsw (npan),lssw(npan),lsee(npan),                     & ! ints
     &               lsse(npan),lnww(npan),lnw (npan),                     & ! ints
     &               lnnw(npan),lnee(npan),lnne(npan),                     & ! ints
     &               npann(0:13),npane(0:13),npanw(0:13),npans(0:13)      ! hordifg
