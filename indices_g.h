!     Global index arrays

!     second line are used in interpolation routines
      integer, dimension(ifull_g) :: iw_g, isw_g, is_g, ise_g, ie_g,     &
     &   ine_g, in_g, iwn_g, ien_g, inn_g, iss_g, iww_g, iee_g, iwu_g,   &
     &   isv_g, iwu2_g, isv2_g, ieu2_g, inv2_g, iev2_g, inu2_g, ieu_g,   &
     &   inv_g, iwwu2_g, issv2_g, ieeu2_g, innv2_g
      integer, dimension(0:npanels) :: lwws_g, lws_g, lwss_g, les_g,     &
     &   lees_g, less_g, lwwn_g, lwnn_g, leen_g, lenn_g, lsww_g,         &
     &   lsw_g, lssw_g, lsee_g, lsse_g, lnww_g, lnw_g, lnnw_g, lnee_g,   &
     &   lnne_g
      integer, dimension(0:13) :: npann_g, npane_g, npanw_g, npans_g

      common/indices_g/ iw_g,isw_g,is_g,ise_g,ie_g,ine_g,in_g,iwn_g,     &
     &     ien_g,inn_g,iss_g,iww_g,iee_g,iwu_g,isv_g,iwu2_g,isv2_g,      &
     &     ieu2_g,inv2_g,iwwu2_g,issv2_g,ieeu2_g,innv2_g,                &
     &     iev2_g,inu2_g,ieu_g, inv_g,lwws_g,lws_g,lwss_g,               &
     &     les_g,lees_g,less_g,lwwn_g,lwnn_g,leen_g,lenn_g,lsww_g,       &
     &     lsw_g,lssw_g,lsee_g,lsse_g,lnww_g,lnw_g,lnnw_g,lnee_g,lnne_g, &
     &     npann_g,npane_g,npanw_g,npans_g
