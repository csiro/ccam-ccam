c      data for biospheric model

      real, dimension(ifull) :: rsmin, sigmf, tgf, res, rmc, tsigmf
      integer, dimension(ifull) :: ivegt, isoilm
      common/dnsib/ rsmin,ivegt,sigmf,tgf,res,rmc,isoilm,tsigmf

ckf     -               tgg(ifull),ssdn(ifull)
ckf tgg now in soilsnow and has ms layers specidfied by scamdim.h
