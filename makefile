#Makefile for GLOBPE on NEC SX4   f90 (from mrd)

# Specify any desired preprocessor or compiler flags here; -R2 for .L files.

FF = /usr/local/f90_230/bin/f90 -Yf/usr/local/f90_230/lib
FF = /usr/local/f90_253/bin/f90 -Yf/usr/local/f90_253/lib
FF =  /SX/usr/bin/sxf90
# Flags for all optimisation levels
#XFLAGS = -Wf"-O nodiv -pvctl noassume" -R2 -float0 -ew -ftrace
#XFLAGS = -g -Wf"-O nodiv -pvctl noassume loopcnt=24000" -R2 -float0 -ew -ftrace -L/usr/local/f90_240/lib
#XFLAGS = -Wf"-O nodiv -pvctl noassume loopcnt=24000" -R2 -float0 -ew -ftrace -L/usr/local/f90_230/lib
XFLAGS = -Wf"-O nodiv -pvctl noassume loopcnt=24000" -R2 -float0 -ew -ftrace 

# Standard optimisation level - vopt usually fine
OPT = -Chopt
OPT = -Cvsafe
OPT = -Csopt
OPT = -Cvopt

#LIBS = /usr/local/lib/libnetcdf64.a
LIBS = /usr/local/lib/libnetcdf3.5-ew.a
LIBS = /cs/u/csdar/csmrd/portal/lib_nec/libnetcdf3.5-ew.a
LIBS = /cs/u/csdar/csmrd/portal/lib_nec/libnetcdf.a

# The following lines specify extra compiler FLAGS. Uncomment one of them
# to use that option or add extra options as required. It is assumed that 
# these are for occasional use only.

# This next section specifies the object files. Only change this part if you
# add or remove new source code files. (no blanks at end of lines!)

OBJS =        adjust5.o     amipsst.o     co2.o         conjob.o    \
betts.o       bett_cuc.o    bettinit.o    bettrain.o    bettspli.o  \
convjlm.o     davies.o      depts.o       esbda.o       gettin.o    \
globpe.o      gwdrag.o      hordifg.o     hs_phys.o     iabsdate.o  \
indata.o      infile.o      ints.o        helmsol.o     jimcc.o     \
maxmin.o      mslp.o        nestin.o      nonlin.o      optmx.o     \
outcdf.o      outfile.o     pbldif.o      printa.o      radriv90.o  \
retopo.o      scrnout.o     setxyz.o      sflux.o       so2.o       \
soilsnow.o    sst.o         staguv.o          \
trim.o        upglobal.o    updps.o       vadv30.o        \
vadvtvd.o     vertmix.o     onthefly.o    latltoij.o  \
esibda.o      icefall.o     leoncld.o     newcloud.o    newrain.o 

NEWRADOBS =   cldblk.o      clddia.o      cldset.o      clo89.o     \
cloud.o  cloud2.o  co2_read.o    e1e288.o      e3v88.o       extras.o    \
fst88.o       hconst.o      lwr88.o       o3_read.o     o3set.o     \
resetd.o      spa88.o       swr99.o       table.o       zenith.o

# This section states that the executable 'model' depends on OBJS and gives
# the command which produces model. Note that the line containing the cf77
# command must begin with a TAB or make will get confused.
# first of the following will be done

hsmodel :$(OBJS) $(NEWRADOBS)
	$(FF) $(XFLAGS) $(OBJS) $(NEWRADOBS) $(LIBS) -o globpea.new9    

# This section gives the rules for building object modules.

.SUFFIXES:.f90
.f90.o:
	$(FF) -c $(OPT) $(XFLAGS) $<
.f.o:
	$(FF) -c $(OPT) $(XFLAGS) $<

# Files which require non-standard optimisation
# version 212 seems OK for infile (prev. vsafe)
# version 212 still needs special treatment for :adjust5(vopt), sflux,soilsnow(vsafe)
# and iabsdate( was faulty for hopt) for differing months
# version 212 needs sopt for indata (just nrungcm=1 section) - OK now
# 18/12/00 it now seems sflux & soilsnow OK, but use vopt to be cautious
# scrnout needs vsafe for C48 runs (but vopt fine with noassume)
#indata.o:indata.f
#	$(FF) -c -Csopt $(XFLAGS) $<
adjust5.o:adjust5.f
	$(FF) -c -Cvopt $(XFLAGS) $<
#fst88.o:fst88.f
#	$(FF) -c -Cvsafe $(XFLAGS) $<
#swr99.o:swr99.f
#	$(FF) -c -Cvsafe $(XFLAGS) $<
sflux.o:sflux.f
	$(FF) -c -Cvopt  $(XFLAGS) $<
soilsnow.o:soilsnow.f
	$(FF) -c -Cvopt  $(XFLAGS) $<
scrnout.o:scrnout.f
	$(FF) -c -Cvopt   $(XFLAGS) $<
staguv.o:staguv.f
	$(FF) -c -Cvopt   $(XFLAGS) $<
#iabsdate.o:iabsdate.f
#	$(FF) -c -Cvopt $(XFLAGS) $<

# Files which depend on include files....You should only need to change this
# section if you add or remove 'include' statements.
# N.B. no blanks at end of line after \

co2.o clddia.o globpe.o indata.o nestin.o outcdf.o outfile.o  \
 :aalat.h

adjust5.o amipsst.o clddia.o co2.o conjob.o convjlm.o davies.o gettin.o \
globpe.o gwdrag.o hordifg.o hs_phys.o indata.o mslp.o nestin.o nonlin.o \
outcdf.o outfile.o pbldif.o radriv90.o scrnout.o sflux.o so2.o \
soilsnow.o updps.o upglobal.o vadvtvd.o vertmix.o :arrays.h 

betts.o bett_cuc.o bettinit.o bettrain.o : betts1.h

depts.o jimcc.o jimco.o latltoij.o onthefly.o setxyz.o : bigxy4.h

clo89.o fst88.o radriv90.o spa88.o : cldcom.h

co2_read.o lwr88.o : co2dta.h

adjust5.o co2.o conjob.o convjlm.o depts.o globpe.o gwdrag.o indata.o \
mslp.o nonlin.o setxyz.o so2.o soilsnow.o upglobal.o \
vertmix.o : constant.h

icefall.f leoncld.o newcloud.o newrain.o : const_phys.h

icefall.f leoncld.o newcloud.o newrain.o : cparams.h

globpe.o infile.o outcdf.o outfile.o : darcdf.h

amipdata.o amipsst.o co2.o davies.o globpe.o indata.o nestin.o onthefly.o \
outcdf.o outfile.o radriv90.o sflux.o  vertmix.o : dates.h

conjob.o convjlm.o davies.o indata.o nestin.o outfile.o : dava.h

clddia.o davies.o nestin.o : davb.h

cloud.o co2.o globpe.o infile.o outcdf.o outfile.o pbldif.o \
radriv90.o sflux.o : extraout.h

amipsst.o globpe.o indata.o outcdf.o outfile.o : filnames.h

gwdrag.o indata.o sflux.o : gdrag.h

e1e288.o e3v88.o fst88.o hconst.o lwr88.o radriv90.o spa88.o swr99.o \
table.o : hcon.h

adjust5.o depts.o globpe.o helmsol.o hordifg.o indata.o ints.o nonlin.o \
onthefly.o setxyz.o staguv.o updps.o upglobal.o vadv30.o  \
vertmix.o : indices.h

e1e288.o fst88.o lwr88.o spa88.o : kdacom.h

clddia.o conjob.o convjlm.o globpe.o outcdf.o outfile.o pbldif.o \
radriv90.o vertmix.o : kuocom.h

amipdata.o globpe.o hs_phys.o indata.o nonlin.o onthefly.o radriv90.o \
setxyz.o sflux.o : latlong.h

convjlm.o globpe.o infile.o leoncld.o \
outcdf.o radriv90.o upglobal.o vadv30.o  vadvtvd.o vertmix.o : liqwpar.h

fst88.o radriv90.o spa88.o : lwout.h

adjust5.o amipsst.o clddia.o co2.o convjlm.o depts.o globpe.o hordifg.o \
indata.o mslp.o nestin.o nonlin.o onthefly.o outcdf.o outfile.o pbldif.o \
radriv90.o scrnout.o setxyz.o sflux.o so2.o staguv.o updps.o \
upglobal.o vadv30.o   vadvtvd.o vertmix.o : map.h

globpe.o infile.o outcdf.o :mapproj.h

adjust5.o betts.o clddia.o conjob.o convjlm.o globpe.o gwdrag.o indata.o \
infile.o leoncld.o nonlin.o outcdf.o outfile.o pbldif.o sflux.o soilsnow.o \
vertmix.o : morepbl.h

globpe.o infile.o outcdf.o outfile.o : netcdf.h

adjust5.o amipdata.o amipsst.o bett_cuc.o bettinit.o bettrain.o betts.o \
clddia.o clo89.o cloud.o co2.o co2_read.o conjob.o convjlm.o davies.o \
depts.o e1e288.o e3v88.o fst88.o gettin.o globpe.o gwdrag.o helmsol.o \
hordifg.o hs_phys.o indata.o infile.o int2.o ints.o jimcc.o jimco.o \
latltoij.o lwr88.o maxmin.o mslp.o nestin.o nonlin.o o3_read.o o3set.o \
onthefly.o optm.o outcdf.o outfile.o pbldif.o printa.o radriv90.o \
retopo.o scrnout.o setxyz.o sflux.o so2.o soilsnow.o spa88.o \
sst.o staguv.o swr99.o table.o trim.o updps.o upglobal.o \
vadv30.o   vadvtvd.o vertmix.o \
icefall.f leoncld.o newcloud.o newrain.o : newmpar.h

adjust5.o conjob.o convjlm.o globpe.o gwdrag.o hordifg.o hs_phys.o \
nonlin.o outfile.o soilsnow.o upglobal.o vadvtvd.o vertmix.o : nlin.h

amipsst.o co2.o globpe.o indata.o outcdf.o outfile.o radriv90.o scrnout.o \
sflux.o so2.o soilsnow.o vertmix.o : nsibd.h

icefall.f newcloud.o newrain.o : params.h

adjust5.o amipdata.o amipsst.o betts.o clddia.o cloud.o co2.o conjob.o \
convjlm.o davies.o depts.o globpe.o gwdrag.o helmsol.o hordifg.o \
hs_phys.o indata.o infile.o ints.o jimcc.o latltoij.o mslp.o nestin.o nonlin.o \
onthefly.o outcdf.o outfile.o pbldif.o radriv90.o retopo.o scrnout.o \
setxyz.o sflux.o so2.o soilsnow.o staguv.o updps.o \
upglobal.o vadv30.o   vadvtvd.o vertmix.o : parm.h

adjust5.o globpe.o helmsol.o indata.o latltoij.o nonlin.o outcdf.o \
outfile.o staguv.o upglobal.o vadv30.o  vadvtvd.o : parmdyn.h

globpe.o ints.o staguv.o upglobal.o : parmhor.h

adjust5.o globpe.o nonlin.o outcdf.o outfile.o vadv30.o  \
 vadvtvd.o : parmvert.h

globpe.o indata.o infile.o onthefly.o : parm_nqg.h

globpe.o upglobal.o : particle.h

adjust5.o amipsst.o clddia.o globpe.o gwdrag.o indata.o nestin.o outcdf.o \
outfile.o radriv90.o scrnout.o sflux.o vertmix.o : pbl.h

betts.o conjob.o convjlm.o globpe.o indata.o leoncld.o outcdf.o outfile.o \
sflux.o : prec.h

clo89.o cloud.o co2_read.o e1e288.o fst88.o lwr88.o radriv90.o spa88.o \
table.o : radisw.h

fst88.o radriv90.o spa88.o : rdflux.h

clo89.o cloud.o co2_read.o e1e288.o e3v88.o fst88.o lwr88.o radriv90.o \
spa88.o swr99.o table.o : rdparm.h

fst88.o lwr88.o spa88.o table.o : rnddta.h

co2.o globpe.o indata.o outcdf.o outfile.o radriv90.o scrnout.o sflux.o \
: scamdim.h

globpe.o infile.o outcdf.o outfile.o sflux.o : screen.h

adjust5.o betts.o clddia.o co2.o conjob.o convjlm.o davies.o globpe.o \
gwdrag.o hordifg.o hs_phys.o indata.o infile.o maxmin.o mslp.o nonlin.o \
onthefly.o outcdf.o outfile.o pbldif.o radriv90.o retopo.o scrnout.o \
sflux.o so2.o soilsnow.o updps.o upglobal.o vadv30.o  \
 vadvtvd.o vertmix.o newcloud.o : sigs.h

amipsst.o clddia.o conjob.o convjlm.o globpe.o gwdrag.o indata.o infile.o \
nestin.o outcdf.o outfile.o radriv90.o scrnout.o sflux.o soilsnow.o \
vertmix.o : soil.h

amipsst.o co2.o globpe.o indata.o infile.o outcdf.o outfile.o radriv90.o \
scrnout.o sflux.o soilsnow.o : soilsnow.h 

globpe.o indata.o outcdf.o outfile.o radriv90.o sflux.o soilsnow.o \
: soilv.h 

fst88.o radriv90.o spa88.o : srccom.h

globpe.o indata.o infile.o nestin.o onthefly.o : stime.h

radriv90.o : swocom.h

e1e288.o e3v88.o fst88.o table.o : tabcom.h

e1e288.o fst88.o lwr88.o radriv90.o spa88.o : tfcom.h

adjust5.o co2.o conjob.o convjlm.o globpe.o indata.o infile.o leoncld.o \
nonlin.o onthefly.o outcdf.o outfile.o sflux.o so2.o upglobal.o \
vadv30.o  vadvtvd.o vertmix.o : tracers.h

co2.o globpe.o indata.o outcdf.o sflux.o : trcom2.h

adjust5.o indata.o : vecs.h

adjust5.o depts.o globpe.o hordifg.o indata.o nonlin.o onthefly.o \
setxyz.o staguv.o upglobal.o : vecsuv.h

adjust5.o globpe.o staguv.o : vecsuva.h

adjust5.o co2.o convjlm.o globpe.o nonlin.o onthefly.o outcdf.o outfile.o \
sflux.o updps.o upglobal.o vadv30.o vadvtvd.o : vvel.h

adjust5.o conjob.o convjlm.o globpe.o hordifg.o nonlin.o updps.o \
upglobal.o vadv30.o vertmix.o : xarrs.h

depts.o globpe.o indata.o maxmin.o nonlin.o onthefly.o setxyz.o \
upglobal.o : xyzinfo.h

radriv90.o : swr99.o zenith.o
