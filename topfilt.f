      program topfilt   ! for conformal-cubic model
c   filt set via nfilt:  0  off                       e.g. save as topout30
c                        1  1111-4 scheme             e.g. save as topout30a
c                        2  1111-4 scheme repeated    e.g. save as topout30b
      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'parm.h'   ! to pass rlong0,rlat0,schmidt  to setxyz
      dimension rmsk(il,jl), tsd(il,jl), zss(il,jl),zssout(il,jl)
      character formout*9
      character header*47
      data id/10/,jd/10/
      open(10,file='topout0',form='formatted',status='old')
      read(10,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
     .          ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      print *,ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      if(ilx.ne.il.or.jlx.ne.jl)
     .            stop 'wrong topo file supplied'
      read (10,*) zss
      read (10,*) rmsk
      read (10,*) tsd
      call setxyz   ! for indices

      do nfilt=1,2
        if(nfilt.eq.1)
     .         open(20,file='topouta',form='formatted',status='unknown')
        if(nfilt.eq.2)
     .         open(20,file='topoutb',form='formatted',status='unknown')
        call filt (zss)  ! 1114 for conformal-cubic

c       reset ocean points via rmsk
        do j=1,jl
         do i=1,il
          zssout(i,j)=zss(i,j)
          if(rmsk(i,j).lt. .5)zssout(i,j)=-1.
         enddo
        enddo

c       write out g*zs(il,jl) to formatted file
        write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  filtered'',a37)')
     .            il,jl,rlong0,rlat0,schmidt,ds,header
        ilout=min(il,30)
        write (formout,'(1h(,i3,4hf7.0,1h))')ilout   !  i.e. (<il>f7.0)
        write(20,formout) zssout

c       write out land/sea mask to formatted file
        write (formout,'(1h(,i3,4hf4.1,1h))')ilout   !  i.e. (<il>f4.1)
        write(20,formout) rmsk

c       write out std.dev of top. to formatted file
        write (formout,'(1h(,i3,4hf6.0,1h))')ilout   !  i.e. (<il>f6.0)
        write(20,formout) tsd
        close (20)
      enddo
      end
      subroutine filt ( a)  ! just 1114  jlm
c     Applies a simple filter over data in 2 dimensions
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee

c input/output variables
      real a(il,jl)  ! input/output array
      real b(il,jl)  ! work array

      print *,'1114 filt called'

c     initialize b array to a
      do iq = 1, ifull
       b(iq,1) = a(iq,1)
      end do

c     1114 filter
      do iq = 1, ifull
       a(iq,1) = .125*(b(iw(iq),1)+b(ie(iq),1)+b(is(iq),1)+b(in(iq),1))
     .           +.5*b(iq,1)
      end do

      return
      end
      include 'setxyz.f'
      include 'jimcc.f'
