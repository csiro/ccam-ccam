! eak 16/03/06
!     canopy parameters
      real, dimension(mxvt)   ::canst1,dleaf,ejmax,frac4,hc,rp20
      real, dimension(mxvt)   ::rpcoef,shelrb,vcmax,xfang
      real, dimension(mxvt)   ::tminvj,tmaxvj,vbeta
      real, dimension(mxvt)   :: wai,vegcf,extkn,rootbeta,xalbnir
      real, dimension(ifull)   ::cansto,vlai
      real, dimension(ifull)   ::sumpn,sumrp,sumrpw,sumrpr,sumrs,sumrd
      real, dimension(ifull)   ::dsumpn,dsumrp,dsumrs,dsumrd
      real, dimension(ifull)   ::rlai,rlaimax
      real, dimension(ifull,3)   ::rlai123
! rml 28/09/07 add c4 fraction for all gridpoints
      real, dimension(ifull) :: c4frac
      common/vegpar/canst1,dleaf,ejmax,frac4,hc,rp20,rpcoef,shelrb
      common/vegpar/xfang,vcmax,tminvj,tmaxvj,vbeta
      common/vegpar/wai,vegcf,extkn,rootbeta,xalbnir
      common/vegpar/cansto,vlai,c4frac
      common/vegpar/sumpn,sumrp,sumrpw,sumrpr,sumrs,sumrd
      common/vegpar/dsumpn,dsumrp,dsumrs,dsumrd
      common/vegpar/rlai,rlaimax,rlai123



