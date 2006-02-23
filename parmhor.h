!     horizontal advection/staggering options (globpe, ints, staguv, upglobal)
      integer nt_adv, ndept
      common/parmhor/nt_adv,ndept
      integer, parameter :: mhint=0  !  0 for usual simple; 2 for Bessel  (in ints)
      
!     for RMIP1 m_bs was -2; during 2002 it was 2
      integer, parameter :: m_bs=-2  !  0 for B&S off     usually -2
!                               2 for B&S on (in ints)
!                              -2 on for gases only
!     m_bs is superseded on 23/7/03 by mh_bs
      integer, parameter :: mh_bs=3  !  5 for B&S off     usually 4
!                               4 for B&S on for qg, gases (in ints)
!                               3 for B&S on for T, qg, gases 
!                               2 for B&S on for u, v, T, qg, gases 
!                               1 for B&S on for psl, u, v, T, qg, gases 
      integer, parameter :: mstagpt=-3 !  2 for 2-point
!                               3 for 3-point
!                              -3 for other (new) 3-point
!                               4 for 4-point     original
!                               5 for 5-point
!                               7 for 7-term, 8 for 8-term

