! eak 16/03/06
!     canopy parameters
      real, dimension(mxvt)   ::canst1,dleaf,ejmax,frac4,hc,rp20
      real, dimension(mxvt)   ::rpcoef,shelrb,vcmax,xfang
      real, dimension(mxvt)   ::tminvj,tmaxvj,vbeta
      real, dimension(ifull)   ::cansto,vlai,vlaimax
      real, dimension(ifull)   ::sumpn,sumrp,sumrpw,sumrpr,sumrs,sumrd
      real, dimension(ifull)   ::dsumpn,dsumrp,dsumrs,dsumrd
      real, dimension(ifull)   ::rlai,rlaimax
      real, dimension(ifull,3)   ::rlai123
      common/vegpar/canst1,dleaf,ejmax,frac4,hc,rp20,rpcoef,shelrb
      common/vegpar/xfang,vcmax,tminvj,tmaxvj,vbeta
      common/vegpar/cansto,vlai,vlaimax
      common/vegpar/sumpn,sumrp,sumrpw,sumrpr,sumrs,sumrd
      common/vegpar/dsumpn,dsumrp,dsumrs,dsumrd
      common/vegpar/rlai,rlaimax,rlai123



