      include 'jimcc.f'
      include 'latltoij.f'
      include 'setxyz.f'
      include 'newmpar.h'
      include 'parm.h'
      include 'xyzinfo.h'  ! x,y,z,rlat,rlong,wts
      include 'vecsuv.h'   ! vecsuv info
      data rlong0/0./,rlat0/90./,schmidt/1./

      call setxyz
      do ilat=-80,80,20
       do ilong=0,360,20
        print *
        print *,'ilong,ilat ',ilong,ilat
        grlong=ilong
        grlat=ilat
        call latltoij(grlong,grlat,xg,yg,nface)
       enddo
      enddo
      end
