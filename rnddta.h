c    common block bandta contains random band parameters for the LW 
c    calculations using 10 cm-1 wide bands. The 15 um CO2 complex
c    is 2 bands,560-670 and 670-800 cm-1. Ozone coefficients are
c    in 3 bands,670-800 (14.1 um),990-1070 and 1070-1200 (9.6 um).
c    The  (nblw) bands now include: 
c                56 bands, 10  cm-1 wide    0  -   560  cm-1
c                 2 bands, 15 um complex  560  -   670  cm-1
c                                         670  -   800  cm-1
c                 3 "continuum" bands     800  -   900  cm-1
c                                         900  -   990  cm-1
c                                        1070  -   1200 cm-1
c                 1 band for 9.6 um band  990  -   1070 cm-1
c               100 bands, 10 cm-1 wide  1200  -   2200 cm-1
c                 1 band for 4.3 um src  2270  -   2380 cm-1
c    Thus nblw presently equals    163
c    All bands are arranged in order of increasing wavenumber 
c 
      real arndm(nblw)  ! Random "a" parameter for (nblw) bands
      real brndm(nblw)  ! Random "b" parameter for (nblw) bands
      real betad(nblw)  ! Continuum coefficients for (nblw) bands
      real ap(nblw)     ! capphi coefficients for (nblw) bands 
      real bp(nblw)     
      real atp(nblw)    ! cappsi coefficients for (nblw) bands 
      real btp(nblw)
      real bandlo(nblw) ! Lowest frequency in each of (nblw) freq. bands 
      real bandhi(nblw) ! Highest frequency in each of (nblw) freq. bands
      real ao3rnd(3)    ! Random "a" parameter for ozone in (3) ozone bands
      real bo3rnd(3)    ! Random "b" parameter for ozone in (3) ozone bands
      real ab15(2)      ! The product arndm*brndm for the two bands
                        ! representing the 15 um band complex of CO2 

c     Data for arndm,brndm,ap,bp,atp,btp,ao3rnd,bo3rnd are obtained by
c     using the AFGL 1982 catalog. Continuum coefficients are from
c     Roberts (1976). 
      
      common / bandta / arndm, brndm, betad, ap, bp, atp, btp, bandlo, 
     &                  bandhi, ao3rnd, bo3rnd, ab15
c 
c    common block bdwide contains random band parameters for specific 
c    wide bands. At present,the information consists of 
c    1) random model parameters for the 15 um band,560-800 cm-1; 
c    2) the continuum coefficient for the 800-990,1070-1200 cm-1 band
c       specifically:  

      real awide        ! Random "a" parameter for  band 
      real bwide        ! Random "b" parameter for  band 
      real betawd       ! Continuum coefficients for band
      real apwd, bpwd   ! capphi coefficients for  band
      real atpwd, btpwd ! cappsi coefficients for band 
      real bdlowd       ! Lowest frequency in each  freq  band 
      real bdhiwd       ! Highest frequency in each freq  band 
      real ab15wd       ! The product arndm*brndm for the one band 
                        ! representing the 15 um band complex of CO2 
      real betinw       ! Cont. coefficient for a specified wide
                        ! freq. band (800-990 and 1070-1200 cm-1).
      real sko2d        ! 1./betinw, used in spa88 for cont. coeffs
      real skc1r        ! betawd/betinw, used for cont. coeff. for 
                        ! 15 um band in fst88
      real sko3r        ! Ratio of cont. coeff. for 9.9 um band to 
                        ! betinw, used for 9.6 um cont coeff in fst88

c     Data for awide,bwide,apwd,bpwd,atpwd,btpwd,ao3wd,bo3wd are
c     obtained by using the AFGL 1982 catalog. Continuum coefficients 
c     are from Roberts (1976).
      common / bdwide / awide,bwide,betawd, 
     &                  apwd,bpwd,atpwd,btpwd,
     &                  bdlowd,bdhiwd,betinw, 
     &                  ab15wd,sko2d,skc1r,sko3r
c 
c    common block bdcomb contains random band parameters for the lw 
c    calculations using combined wide frequency bands between 160 and 
c    1200 cm-1,as well as the 2270-2380 band for source calc. 
c        bands 1-8: combined wide frequency bands for 160-560 cm-1
c        bands 9-14: frequency bands,as in bandta (narrow bands)
c                    for 560-1200 cm-1
c        band  15:  frequency band 2270-2380 cm-1,used for source 
c                   calculation only
c        thus nbly presently equals   15
c 
c        bands are arranged in order of increasing wavenumber 

      real acomb(nbly)  ! Random "a" parameter for (nbly) bands
      real bcomb(nbly)  ! Random "b" parameter for (nbly) bands
      real betacm(nbly) ! Continuum coefficients for (nbly) bands
      real apcm(nbly)   ! capphi coefficients for (nbly) bands 
      real bpcm(nbly) 
      real atpcm(nbly)  ! cappsi coefficients for (nbly) bands 
      real btpcm(nbly)
      real bdlocm(nbly) ! Lowest frequency in each of (nbly) freq. bands 
      real bdhicm(nbly) ! Highest frequency in each of (nbly) freq. bands
      real ao3cm(3)     ! Random "a" parameter for ozone in (3) ozone bands
      real bo3cm(3)     ! Random "b" parameter for ozone in (3) ozone bands
      real ab15cm(2)    ! The product arndm*brndm for the two bands
                        ! representing the 15 um band complex of CO2 
      real betinc       ! Cont. coefficient for a specified wide
                        ! freq.band (800-990 and 1070-1200 cm-1).
      integer iband(40) ! Index no of the 40 wide bands used in
                        ! combined wide band calculations. In other
                        ! words, index telling which of the 40 wide 
                        ! bands between 160-560 cm-1 are included in 
                        ! each of the first 8 combined wide bands

c     Data for acomb,bcomb,apcm,bpcm,atpcm,btpcm,ao3cm,bo3cm are
c     obtained by using the AFGL 1982 catalog. Continuum coefficients 
c     are from Roberts (1976). iband index values are obtained by 
c     experimentation.
      common / bdcomb / iband, acomb, bcomb, betacm, apcm, bpcm, atpcm,  
     &                  btpcm, bdlocm, bdhicm, betinc, ao3cm, bo3cm,
     &                  ab15cm
