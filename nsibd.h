c      data for biospheric model (elai added for nsib=5)

      real, dimension(ifull) :: rsmin, sigmf, tgf, res, rmc, tsigmf
      real, dimension(ifull) :: elai, sigmu,ualb  ! MJT CHANGE - add ualb delete tsigmu
      integer, dimension(ifull) :: ivegt, isoilm
      common/dnsib/ rsmin,ivegt,sigmf,tgf,res,rmc,isoilm,tsigmf,
     &              elai,sigmu,ualb    ! MJT CHANGE - add ualb delete tsigmu

