! $Log: establ.h,v $
! Revision 1.2  2003/09/05 02:34:26  dix043
! Make sure variable names won't clash with anything else.
!
! Revision 1.1.1.1  2003/08/13 01:24:20  dix043
! Imported sources
!
! Revision 1.7  2000/11/14 03:11:35  rot032
! Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
!
! Revision 1.6  1998/12/10 00:55:29  ldr
! HBG changes to V5-1-21
!
! Revision 1.5  1996/11/28  00:59:35  mrd
! Add type declaration of pp.
!
! Revision 1.4  1996/03/21  03:18:28  ldr
! Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
!
! Revision 1.3  1995/06/30  02:44:39  ldr
! Changes to qcloud (run F88) and put qsat formula into ESTABL.f
!
! Revision 1.2  1994/12/07  23:53:36  mrd
! Explicitly declare all variables in include files.
!
! Revision 1.1  1991/03/14  09:49:41  ldr
! Initial revision
!
      real table, tablei
      common /es_table/ table(0:220)
      common /esitable/ tablei(0:220)

! Arithmetic statement functions to replace call to establ.
! T is temp in Kelvin, which should lie between 123.16 and 343.16;
! TDIFF is difference between T and 123.16, subject to 0 <= TDIFF <= 220

      real tdiff, establ, qsat, estabi, qsati
      ! Use local names that shouldn't conflict with anything else
      real t_, pp_
      tdiff(t_)=min(max( t_-123.16, 0.), 219.)

      ! MJT suggestion
      establ(t_) =                                                       &
     &       (1.-(tdiff(t_)-anint(tdiff(t_))))*table(nint(tdiff(t_)))    &
     &     + (tdiff(t_)-anint(tdiff(t_)))*table(nint(tdiff(t_))+1)
!      establ(t_) =                                                       &
!     &       (1.-(tdiff(t_)-aint(tdiff(t_))))*table(int(tdiff(t_)))      &
!     &     + (tdiff(t_)-aint(tdiff(t_)))*table(int(tdiff(t_))+1)
      qsat(pp_,t_) = epsil*establ(t_)/max(pp_-establ(t_),.1) !jlm strato
!     qsat(pp_,t_) = epsil*establ(t_)/(pp_-establ(t_)) !Usual formula
!     qsat(pp_,t) = epsil*establ(t_)/pp_ !Consistent with V4-5 to V4-7

! These give the ice values needed for the qcloud scheme
      estabi(t_) =                                                       &
     &      (1.-(tdiff(t_)-anint(tdiff(t_))))*tablei(nint(tdiff(t_)))    &
     &    + (tdiff(t_)-anint(tdiff(t_)))*tablei(nint(tdiff(t_))+1)
!      estabi(t_) =                                                       &
!     &      (1.-(tdiff(t_)-aint(tdiff(t_))))*tablei(int(tdiff(t_)))      &
!     &    + (tdiff(t_)-aint(tdiff(t_)))*tablei(int(tdiff(t_))+1)
      qsati(pp_,t_) = epsil*estabi(t_)/max(pp_-estabi(t_),.1) !jlm strato
!     qsati(pp_,t_) = epsil*estabi(t_)/(pp_-estabi(t_)) !Usual formula
!     qsati(pp_,t_) = epsil*estabi(t_)/pp_ !Consistent with V4-5 to V4-7
