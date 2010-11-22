!      data for biospheric model

      real, dimension(ifull) :: rsmin, sigmf, tgf, res, tsigmf ! MJT aerosols - delete rmc
      real, dimension(ifull) :: sigmu  ! MJT cable - delete elai
      integer, dimension(ifull) :: ivegt, isoilm
      common/dnsib/ rsmin,ivegt,sigmf,tgf,res,isoilm,tsigmf,sigmu  ! MJT cable - delete elai ! MJT aerosols - delete rmc

