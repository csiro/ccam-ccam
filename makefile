FC = ifc
INC = -I /home/dix043/include_ifc  -I /home/dix043/mpich-1.2.2/include

# Need -fpp for a couple of routines
FFLAGS = -O -WB -w -fpp

LIBS = -L /home/dix043/lib -lnetcdf_ifc -lmrd90_ifc -Vaxlib -L/home/dix043/mpich-1.2.2/lib -lmpich_ifc
LDFLAGS = 

SRC =        adjust5.f     amipsst.f     co2.f         conjob.f    \
betts.f       bett_cuc.f    bettinit.f    bettrain.f    bettspli.f  \
convjlm.f     depts.f       esbda.f       gettin.f    \
globpe.f      gwdrag.f      hordifg.f     hs_phys.f     iabsdate.f  \
indata.f      infile.f      ints.f        helmsol.f     jimcc.f     \
mslp.f        nestin.f      nonlin.f         \
outcdf.f      outfile.f     pbldif.f      radriv90.f  \
retopo.f      scrnout.f     setxyz.f      sflux.f       so2.f       \
soilsnow.f    sst.f         staguv.f          \
trim.f        upglobal.f    updps.f       vadv30.f        \
vadvtvd.f     vertmix.f     \
esibda.f      icefall.f     leoncld.f     newcloud.f    newrain.f \
dummy.f90 \
cldblk.f      clddia.f      cldset.f      clo89.f     \
cloud.f  cloud2.f  co2_read.f    e1e288.f      e3v88.f       extras.f    \
fst88.f       hconst.f      lwr88.f       o3_read.f     o3set.f     \
resetd.f      spa88.f       swr99.f       table.f       zenith.f90 \
cc_mpi.f90  diag_m.f90

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $(FFLAGS) $(INC) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INC) $<

# Remove mod rule from Modula 2
%.o : %.mod

# Include the dependency-list created by makedepf90 below
include .depend

clean:
	rm *.o *.mod globpea .depend

depend .depend:
	makedepf90 -o globpea $(SRC) > .depend

