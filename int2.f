c This is version incorporating jlm and kjw versions as of August 1994.
c 
c called by startnst and original-style nestin
c Initial revision
c 
      subroutine int2(tt,ta,i,j,x,ii,y,jj,ktop)
c     Bilinear interpolation from array tt to ta
c     Only the bottom LH corner of tt is used; its dimension are set as
c     il,jl,kl nowadays when it is read via infile
      include 'newmpar.h'
      dimension tt(il,jl,kl),ta(il,jl,kl)
      do 2 k=1,ktop
2     ta(i,j,k)=(1.-x)*((1.-y)*tt(ii,jj,k)+y*tt(ii,jj+1,k))
     .             +x*((1.-y)*tt(ii+1,jj,k)+y*tt(ii+1,jj+1,k))
      return
      end
