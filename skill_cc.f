      program skill_cc
c     computes MSL s1 skillscore and rms on BoM array of grid points
c     assumes data on 1 degree grid longs: 100 to 190, lats -55 to 0

      real anal(91,56),prog(91,56)

!     read in anal & prog for appropriate date/time (jjk to do)
      call score(anal,prog,bias,rms,score)
      print *,'bias, rms, score: ',bias,rms,score
      end

      subroutine score (anal,prog,bias,rms,score)
      real anal(91,56),prog(91,56)
      common /popf/ po(58), pf(58)
      dimension lts(58),lgs(58),sk(92),xlt(92),xln(92)
      dimension npt1(74),npt2(74),npt3(92),npt4(92)
         data npt1    /1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20,21
     &,23,24,25,26,27,28,29,31,32,33,34,35,36,38,39,40,41,42,43,45,46,47
     &,8,1,16,9,24,2,17,10,25,3,18,11,26,4,19,34,12,27,5,20,35,13,28,6,
     &21,36,14,29,7,22,37,15,30/
         data npt2    /2,3,4,5,6,7,9,10,11,12,13,14,15,17,18,19,20,21,22
     &,24,25,26,27,28,29,30,32,33,34,35,36,37,39,40,41,42,43,44,46,47,48
     &,23,16,31,24,38,17,32,25,39,18,33,26,40,19,34,45,27,41,20,35,46,28
     &,42,21,36,47,29,43,22,37,48,30,44/
         data npt3    /1,2,3,4,5,6,7,9,10,11,12,13,14,16,17,18,19,20,21
     &,22,24,25,26,27,28,29,31,32,33,34,35,36,37,39,40,41,42,43,44,46,
     &47,48,49,50,51,53,54,55,57,1,16,9,24,2,17,32,10,25,3,18,33,11,26,4
     &,19,34,12,27,42,5,20,35,13,28,43,6,21,36,50,14,29,44,7,22,37,51,
     &15,30,45,8,23,38/
         data npt4    /2,3,4,5,6,7,8,10,11,12,13,14,15,17,18,19,20,21,
     &22,23,25,26,27,28,29,30,32,33,34,35,36,37,38,40,41,42,43,44,45,47
     &,48,49,50,51,52,54,55,56,58,16,31,24,39,17,32,46,25,40,18,33,47,26
     &,41,19,34,48,27,42,53,20,35,49,28,43,54,21,36,50,57,29,44,55,22
     &,37,51,58,30,45,56,23,38,52/
      data lts / 8*15,7*20,8*25,7*30,8*35,7*40,7*45,4*50,2*55 / ! NB. SHem
      data lgs / 100,110,120,130,140,150,160,170,
     &           105,115,125,135,145,155,165,
     &           100,110,120,130,140,150,160,170,
     &           105,115,125,135,145,155,165,
     &           100,110,120,130,140,150,160,170,
     &           105,115,125,135,145,155,165,
     &                110,120,130,140,150,160,170,
     &                          135,145,155,165,
     &                                    150,160 /

      na=1
      mx=58
      ioffset=99
      joffset=56
      bias=0.
      rms=0.
      do m=1,mx
c      set up po and pf arrays
c      assumes data on 1 degree grid   longs: 100 to 190, lats -55 to 0
c      with      real anal(91,56),prog(91,56)
       po(m)=anal(lgs(m)-ioffset,joffset-lts(m))
       pf(m)=prog(lgs(m)-ioffset,joffset-lts(m))
       d = pf(m) - po(m)
       bias = bias + d
       rms = rms + d*d
      enddo
      bias= bias/mx
      rms = sqrt(rms/mx)
      efd=0.
      amxdif=0.
      do m=1,mx
c      if ( na .eq. 1 ) then
          i=npt3(m)
          j=npt4(m)
c       else         ! this part for 500 hPa height only
c         i=npt1(m)
c         j=npt2(m)
c       endif
        df=pf(i)-pf(j)
        do=po(i)-po(j)
        amxdif=amxdif+max(abs(do),abs(df))
        efd=efd+abs(df-do)
           if (na.eq.1) then
           xlt(m) = 0.5 * (lts(i)+lts(j))
           xln(m) = 0.5 * (lgs(i)+lgs(j))
           xdif = max(abs(do),abs(df))
           top = abs(df-do)
           sk(m) = 100.*top/xdif
           endif
      enddo
      score = (efd/amxdif)*100.
c          if (na.eq.1)then
c          do m=1,mx
c          print 6000,m,xlt(m),xln(m),sk(m)
c          enddo
c6000      format(1x,i3,3f6.1)
c          endif

      return
      end

