      subroutine gettin(n)
c  gettina reads back preceding u,v,t,psl for n=2 (for nmi, assuming mex=1)
c  saves and re-reads initial arrays of t and psl
c  called only by darlam and vmodes
      include 'newmpar.h'
      include 'arrays.h'
      include 'savuvt.h'
      if(n.eq.0)then
        savt(:,:)=t(:,:)
        savpsl(:)=psl(:)
      else    ! for n=1, n=2
        t(:,:)=savt(:,:)   ! for n=1, n=2
        psl(:)=savpsl(:)   ! for n=1, n=2
      endif

c     following called by vmodes during initialization
      if(n.eq.2)then
        u(:,:)=savu(:,:)  ! only for n=2 (VMI init)
        v(:,:)=savv(:,:)
      endif
      return
      end
