      program topmerge   ! for conformal-cubic model   jlm 14/7/01
      include 'newmpar.h'
      dimension rmsk(ifull), tsd(ifull), zss(ifull),zssout(ifull)
      dimension rmsk1(ifull), tsd1(ifull), zss1(ifull)
      character formout*9
      character header*47

      open(10,file='top_coarse',form='formatted',status='old')  ! coarse res
      open(15,file='top_fine',form='formatted',status='old')  ! fine res

      read(10,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
     .          ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      print *,ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      if(ilx.ne.il.or.jlx.ne.jl)
     .            stop 'wrong topo file supplied'
      read (10,*) zss1
      read (10,*) rmsk1
      read (10,*) tsd1

      read(15,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
     .          ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      print *,ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      if(ilx.ne.il.or.jlx.ne.jl)
     .            stop 'wrong topo file supplied'
      read (15,*) zss
      read (15,*) rmsk
      read (15,*) tsd

      open(20,file='topout0',form='formatted',status='unknown')
      do iq=1,ifull
        if(zss(iq).lt.-.1)then
          zss(iq)=zss1(iq)
          rmsk(iq)=rmsk1(iq)
          tsd(iq)=tsd1(iq)
        endif
      enddo

c       write out g*zs(ifull) to formatted file
        write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  merged'',a37)')
     .            il,jl,rlong0,rlat0,schmidt,dsx,header
        ilout=min(il,30)
        write (formout,'(1h(,i3,4hf7.0,1h))')ilout   !  i.e. (<il>f7.0)
        write(20,formout) zss

c       write out land/sea mask to formatted file
        write (formout,'(1h(,i3,4hf4.1,1h))')ilout   !  i.e. (<il>f4.1)
        write(20,formout) rmsk

c       write out std.dev of top. to formatted file
        write (formout,'(1h(,i3,4hf6.0,1h))')ilout   !  i.e. (<il>f6.0)
        write(20,formout) tsd
        close (20)
      end

