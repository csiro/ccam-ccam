      integer, parameter :: nstnmax=10, nstn2=3
      integer mstn,nstn,istn,jstn,iunp,nstn2,istn2,jstn2,iunp2,
     &        iaustw,iauste,iaustn,iausts
      real slat,slon,slat2,slon2
      common / trcom2 / mstn,nstn,slat(nstnmax),slon(nstnmax),
     &                  istn(nstnmax),jstn(nstnmax),iunp(nstnmax)
     & ,slat2(nstn2),slon2(nstn2),istn2(nstn2),jstn2(nstn2),iunp2(nstn2)
     & ,iaustw,iauste,iaustn,iausts
