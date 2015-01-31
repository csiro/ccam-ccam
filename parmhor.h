!     horizontal advection/staggering options (globpe, ints, staguv, upglobal)

      integer mh_bs, nt_adv
      common/parmhor/nt_adv,mh_bs  ! here from June '06

!     for RMIP1 m_bs was -2; during 2002 it was 2

!     integer, parameter :: m_bs=-2  !  0 for B&S off     usually -2

!                               2 for B&S on (in ints)

!                              -2 on for gases only

!     m_bs is superseded on 23/7/03 by mh_bs

!     integer, parameter :: mh_bs=4  !  5 for B&S off     usually 4

!                               4 for B&S on for qg, gases (in ints)

!                               3 for B&S on for T, qg, gases 

!                               2 for B&S on for u, v, T, qg, gases 

!                               1 for B&S on for psl, u, v, T, qg, gases
