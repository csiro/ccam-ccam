      subroutine gettin(n)
      use arrays_m
c  gettina reads back preceding u,v,t,psl for n=2 (for nmi, assuming mex=1)
c  saves and re-reads initial arrays of t and psl
c  called only by darlam and vmodes
      include 'newmpar.h'
      include 'savuvt.h'
      if(n.eq.0)then
        savt(:,:)=t(1:ifull,:)
        savpsl(:)=psl(1:ifull)
      else    ! for n=1, n=2
        t(1:ifull,:)=savt(:,:)   ! for n=1, n=2
        psl(1:ifull)=savpsl(:)   ! for n=1, n=2
      endif

c     following called by vmodes during initialization
      if(n.eq.2)then
        u(1:ifull,:)=savu(:,:)  ! only for n=2 (VMI init)
        v(1:ifull,:)=savv(:,:)
      endif
      return
      end
