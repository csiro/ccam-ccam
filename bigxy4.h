c     real*8 xx4,yy4
      integer, parameter :: mydouble = selected_real_kind (10,100)
      real (kind=mydouble) xx4,yy4
      common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
