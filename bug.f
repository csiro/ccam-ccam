      parameter (il=48,jl=48*6)
      common/tst/istn(10),jstn(10),alb(il,jl),sigmf(il,jl)

      print *,'testing vopt bug; works OK with -Cvsafe'

c     istn(1)=48
c     jstn(1)=200
c     istn(2)=44
c     jstn(2)=203
c     istn(3)=47
c     jstn(3)=196
c     istn(4)=45
c     jstn(4)=202
c     istn(5)=38
c     jstn(5)=2
c     istn(6)=40
c     jstn(6)=203
c     istn(7)=37
c     jstn(7)=102
c     istn(8)=46
c     jstn(8)=114
c     istn(9)=48
c     jstn(9)=107
c     istn(10)=46
c     jstn(10)=201

      sigmf(37,102)=.33
      print *,'starting sigmf(37,102): ',sigmf(37,102)
c       the following do can be 1,10 or 7,7 but the only
c       faulty value (of sigmf) occurs for nn=7
c       N.B. the fault needs the call to the subroutine
c       and also needs the alb() lines to be set
c       Other values than .195 give the same problem

        do nn=6,8
         call testtest(nn,istn(nn),jstn(nn))
         ii=istn(nn)
         jj=jstn(nn)
         iq=istn(nn)+(jstn(nn)-1)*il
         if(nn.eq.2)alb(ii,jj)=.16   ! fix-up for Sydney
         if(nn.eq.3)alb(ii,jj)=.24   ! fix-up for Adelaide
         if(nn.eq.6)alb(ii,jj)=.18   ! fix-up for Brisbane
         if(nn.eq.7)alb(ii,jj)=.18   ! fix-up for Perth
         if(nn.eq.8)alb(ii,jj)=.12   ! fix-up for Darwin
         if(nn.eq.10)alb(ii,jj)=.16  ! fix-up for Albury
         if(nn.eq.2)sigmf(ii,jj)=.5    ! fix-up for Sydney
         if(nn.eq.3)sigmf(ii,jj)=.5    ! fix-up for Adelaide
         if(nn.eq.6)sigmf(ii,jj)=.5    ! fix-up for Brisbane
         if(nn.eq.7)sigmf(ii,jj)=.195  ! fix-up for Perth
         if(nn.eq.8)sigmf(ii,jj)=.5    ! fix-up for Darwin
         print *,nn,istn(nn),jstn(nn),iq,alb(ii,jj),sigmf(ii,jj)
        enddo  ! nn
        end

      subroutine testtest(nn,istn,jstn)
      if(nn.eq.1)istn=48
      if(nn.eq.1)jstn=200
      if(nn.eq.2)istn=44
      if(nn.eq.2)jstn=203
      if(nn.eq.3)istn=47
      if(nn.eq.3)jstn=196
      if(nn.eq.4)istn=45
      if(nn.eq.4)jstn=202
      if(nn.eq.5)istn=38
      if(nn.eq.5)jstn=2
      if(nn.eq.6)istn=40
      if(nn.eq.6)jstn=203
      if(nn.eq.7)istn=37
      if(nn.eq.7)jstn=102
      if(nn.eq.8)istn=46
      if(nn.eq.8)jstn=114
      if(nn.eq.9)istn=48
      if(nn.eq.9)jstn=107
      if(nn.eq.10)istn=46
      if(nn.eq.10)jstn=201
      return
      end

