      subroutine clo89
c 
CDIR$ TASK COMMON CLDCOM
CDIR$ TASK COMMON RADISW
c CDIR$ TASK COMMON VTEMP

c     subroutine clo88 computes cloud transmission functions for the
c  longwave code,using code written by bert katz (301-763-8161).
c  and modified by dan schwarzkopf in december,1988.
c                inputs:          (common block)
c      camt,ktop,kbtm,nclds,emcld   radisw
c                output:  
c      cldfac                       cldcom
c 
c          called by:      radmn or model routine 
c          calls    : 
c 
      include 'newmpar.h'
      include 'rdparm.h'
      include 'radisw.h'
      include 'cldcom.h'
c 
      real tempc(lp1,lp1,imax),cldfip(lp1,lp1)
c 
      do 11 ip=1,imax
      if (nclds(ip).eq.0) then
        do 29 i=1,lp1*lp1
        tempc(i,1,ip)=1.
29      continue
      endif
      if (nclds(ip).ge.1) then
          xcld=1.-camt(ip,2)*emcld(ip,2)
           k1=ktop(ip,2)+1
           k2=kbtm(ip,2)
          do 31 i=1,lp1*lp1
              cldfip(i,1)=1.
31        continue
          do 41 k=k1,lp1
          do 41 kp=1,k2
               cldfip(kp,k)=xcld
41        continue
          do 43 k=1,k2
          do 43 kp=k1,lp1
              cldfip(kp,k)=xcld
43        continue
            do 61 i=1,lp1*lp1
          tempc(i,1,ip)=cldfip(i,1)
61        continue
      endif
      if (nclds(ip).ge.2) then
        do 21 nc=2,nclds(ip)
          xcld=1.-camt(ip,nc+1)*emcld(ip,nc+1)
           k1=ktop(ip,nc+1)+1
           k2=kbtm(ip,nc+1)
          do 32 i=1,lp1*lp1
              cldfip(i,1)=1.
32        continue
          do 42 k=k1,lp1
          do 42 kp=1,k2
               cldfip(kp,k)=xcld
42        continue
          do 44 k=1,k2
          do 44 kp=k1,lp1
              cldfip(kp,k)=xcld
44        continue
            do 62 i=1,lp1*lp1
          tempc(i,1,ip)=tempc(i,1,ip)*cldfip(i,1)
62        continue
21        continue
      endif
11    continue
      do 70 ip=1,imax
      do 70 i=1,lp1*lp1
         cldfac(ip,i,1)=tempc(i,1,ip)
70      continue
      return
      end 
