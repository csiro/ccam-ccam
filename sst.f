c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:31  mrd
c Initial revision
c
c Revision 1.2  94/09/02  11:55:02  kjw
c This is version incorporating jlm and kjw versions as of August 1994.
c 
c Revision 1.1  91/07/19  11:41:48  ldr
c Month >12 fixed up Thu  11-12-1992
c Initial revision
c 
      subroutine sst(tss,land,ldays,month)
c A routine for nested model which will interpolate
c linearly between two sst data sets.
c LDAYS is a variable from CSIRO9 which is equal to the number of
c whole days of the month that have been completed.(Thus ldays = 0 at the
c start of a month, is updated to 1 at end of first day etc.) It is
c written to the file plotst which is read by the nested model,
c in place of the variable ktau for now.   (LDR 12/90)
c MONTH is the number of the month we are in (1<=month<=12)
c As coded, this assumes that all 12 months' sst's have been read from
c a file (probably lamb1.sst) at the start of the run into the big array SST
      include 'newmpar.h'
      logical land(il,jl)
      integer mdays(12),num
      real tss(il,jl)
      real bigsst(il,jl,12)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31/,num/0/
      save num,bigsst

      print *,'sst  called'
      if(num.eq.0)then
        num=1
        do 20 mon=1,12
          read(17,'(10f8.3)')((bigsst(i,j,mon),i=1,il),j=1,jl)
 20     continue
        print *,'bigsst',bigsst(25,15,1),bigsst(25,15,12)
      endif

c If end of month, swap sst's

      if(ldays.eq.mdays(month))then
        ldays=0
        month=month+1
        if(month.gt.12)month=1   ! Thu  11-12-1992
        print *,'End of month detected in routine sst'
      endif

c Ratio used for interpolation between previous and next months sst's
c Note that Hal's model uses ratlm=ldays/(mdays(month)-1.), but we suspect
c that this is probably wrong.

      ratlm=ldays/real(mdays(month))
      print *,'month=',month
      print *,'ldays=',ldays
      print *,'ratlm=',ratlm

c Interpolate, making sst's positive

      do 40 j=1,jl
        do 30 i=1,il
          if(.not.land(i,j))then
          if(month.eq.12) then
            tss(i,j)=bigsst(i,j,month)*(ratlm-1)-bigsst(i,j,1)*ratlm
          else
            tss(i,j)=bigsst(i,j,month)*(ratlm-1)
     1       - bigsst(i,j,month+1)*ratlm
          end if
          endif
 30     continue
 40   continue

      return
      end
