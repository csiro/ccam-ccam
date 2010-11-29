FC = ifort

#FFLAGS = -O -fpp -I /tools/netcdf/3.6.0-p1/include -assume buffered_io -Dsimple_timer -Duniform_decomp
FFLAGS = -O -fpp -I /tools/netcdf/3.6.0-p1/include -assume buffered_io -Dsimple_timer
LIBS = -L /tools/netcdf/3.6.0-p1/lib -lnetcdf -lmpi

LDFLAGS = 

OBJS = adjust5.o amipsst.o conjob.o betts.o bett_cuc.o bettinit.o \
bettrain.o bettspli.o convjlm.o depts.o esbda.o gettin.o globpe.o gwdrag.o \
hordifg.o hs_phys.o iabsdate.o indata.o infile.o ints.o helmsol.o jimcc.o\
helmsor.o optmx.o\
mslp.o nestin.o nonlin.o outcdf.o outfile.o pbldif.o radriv90.o retopo.o \
scrnout.o setxyz.o sflux.o soilsnow.o staguv.o trim.o upglobal.o eig.o \
updps.o vadv30.o vadvtvd.o vertmix.o esibda.o icefall.o leoncld.o newcloud.o \
newrain.o latltoij.o cldblk.o clddia.o cldset.o clo89.o cloud.o \
cloud2.o co2_read.o e1e288.o e3v88.o extras.o fst88.o hconst.o lwr88.o \
o3_read.o resetd.o spa88.o swr99.o table.o zenith.o cc_mpi.o \
diag_m.o sumdd_m.o ilu_m.o davies.o utilities.o onthefly.o o3read_amip.o \
o3set_amip.o tracermodule.o timeseries.o trvmix.o  stacklimit.o \
xyzinfo_m.o vecsuv_m.o map_m.o latlong_m.o indices_m.o bigxy4_m.o \
arrays_m.o betts1_m.o carbpools_m.o cldcom_m.o co2dta_m.o dava_m.o \
davb_m.o extraout_m.o gdrag_m.o histave_m.o kdacom_m.o kuocomb_m.o \
liqwpar_m.o lwout_m.o morepbl_m.o nlin_m.o nsibd_m.o o3amip_m.o pbl_m.o \
permsurf_m.o prec_m.o raddiag_m.o radisw_m.o rdflux_m.o savuvt_m.o \
screen_m.o sigs_m.o soil_m.o \
cable_ccam2.o cable_air.o cable_albedo.o cable_canopy.o \
cable_carbon.o cable_define_dimensions.o cable_define_types.o \
cable_math_constants.o cable_other_constants.o cable_photosynthetic_constants.o \
cable_physical_constants.o cable_radiation.o cable_roughness.o cable_soilsnow.o \
ateb.o mlo.o tkeeps.o \
seaesfrad.o rad_utilities.o microphys_rad.o esfsw_driver.o esfsw_parameters.o \
longwave_params.o sealw99.o longwave_clouds.o longwave_fluxes.o longwave_tables.o \
optical_path.o gas_tf.o lw_gases_stdtf.o

globpea: $(OBJS)
	$(FC) -o globpea $(FFLAGS) $(LDFLAGS) $(OBJS) $(LIBS)

clean:
	rm *.o *.mod globpea

.SUFFIXES:.f90 .F90

esfsw_driver.o: esfsw_driver.f90
	$(FC)  -c -r8 -override-limits $(FFLAGS) $<
esfsw_parameters.o: esfsw_parameters.f90
	$(FC)  -c -r8 $(FFLAGS) $<
gas_tf.o: gas_tf.f90
	$(FC)  -c -r8 $(FFLAGS) $<
longwave_clouds.o: longwave_clouds.f90
	$(FC)  -c -r8 $(FFLAGS) $<
longwave_fluxes.o: longwave_fluxes.f90
	$(FC)  -c -r8 $(FFLAGS) $<
longwave_tables.o: longwave_tables.f90
	$(FC)  -c -r8 $(FFLAGS) $<
longwave_params.o: longwave_params.f90
	$(FC)  -c -r8 $(FFLAGS) $<
lw_gases_stdtf.o: lw_gases_stdtf.f90
	$(FC)  -c -r8 $(FFLAGS) $<
microphys_rad.o: microphys_rad.f90
	$(FC)  -c -r8 $(FFLAGS) $<
optical_path.o: optical_path.f90
	$(FC)  -c -r8 $(FFLAGS) $<
rad_utilities.o: rad_utilities.f90
	$(FC)  -c -r8 $(FFLAGS) $<
sealw99.o: sealw99.f90
	$(FC)  -c -r8 $(FFLAGS) $<
ateb.o: ateb.f90
	$(FC)  -c -override-limits $(FFLAGS) $<
cable_canopy.o: cable_canopy.F90
	$(FC)  -c -override-limits $(FFLAGS) $<
onthefly.o: onthefly.f
	$(FC)  -c -override-limits $(FFLAGS) $<	
stacklimit.o: stacklimit.c
	cc -c stacklimit.c


.f90.o:
	$(FC) -c $(FFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $<	
.f.o:
	$(FC) -c $(FFLAGS) $<

# Remove mod rule from Modula 2 so GNU make doesn't get confused
%.o : %.mod

# Dependencies
adjust5.o : adjust5.f xyzinfo_m.o xarrs.h vvel.h vecsuv_m.o vecs.h tracers.h sigs_m.o pbl_m.o parmvert.h parmdyn.h parm.h nlin_m.o morepbl_m.o map_m.o liqwpar_m.o kuocom.h indices_m.o const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o tracermodule.o tkeeps.o
amipsst.o : amipsst.f soilsnow.h soil_m.o pbl_m.o parm.h nsibd_m.o map_m.o filnames.h dates.h arrays_m.o newmpar.h cc_mpi.o mlo.o permsurf_m.o
bett_cuc.o : bett_cuc.f betts1_m.o newmpar.h 
bettinit.o : bettinit.f betts1_m.o newmpar.h 
bettrain.o : bettrain.f betts1_m.o newmpar.h 
betts.o : betts.f sigs_m.o prec_m.o parm.h morepbl_m.o betts1_m.o newmpar.h 
bettspli.o : bettspli.f 
cable_ccam2.o : zenith.o cable_define_dimensions.o cable_albedo.o cable_canopy.o cable_albedo.o cable_carbon.o cable_soilsnow.o tracers.h cc_mpi.o radisw_m.o
cable_canopy.o: cable_photosynthetic_constants.o cable_radiation.o cable_roughness.o cable_air.o cable_define_types.o cable_physical_constants.o
cable_photosynthetic_constants.o: cable_define_dimensions.o
cable_radiation.o: cable_math_constants.o cable_other_constants.o cable_define_types.o cable_physical_constants.o 
cable_albedo.o: cable_math_constants.o cable_other_constants.o cable_define_types.o cable_physical_constants.o 
cable_define_types.o: cable_define_dimensions.o
cc_mpi.o : cc_mpi.f90 sigs_m.o  parm.h indices_m.o indices_m.o latlong_m.o latlong_m.o vecsuv_m.o vecsuv_m.o map_m.o map_m.o xyzinfo_m.o xyzinfo_m.o newmpar.h sumdd_m.o 
cldblk.o : cldblk.f 
cldcom.o : cldcom.f 
clddia.o : clddia.f vvel.h soil_m.o sigs_m.o pbl_m.o parm.h morepbl_m.o map_m.o kuocom.h davb_m.o const_phys.h arrays_m.o newmpar.h cc_mpi.o 
cldset.o : cldset.f const_phys.h 
clo89.o : clo89.f cldcom_m.o radisw_m.o rdparm.h newmpar.h 
cloud2.o : cloud2.f radisw_m.o hcon.h rdparm.h params.h kuocom.h cparams.h const_phys.h newmpar.h 
cloud.o : cloud.f radisw_m.o rdparm.h parm.h extraout_m.o newmpar.h 
co2blk.o : co2blk.f 
co2dta.o : co2dta.f 
co2_read.o : co2_read.f radisw_m.o co2dta_m.o rdparm.h newmpar.h cc_mpi.o 
co2trn.o : co2trn.f 
conjob.o : conjob.f establ.h tracers.h soil_m.o sigs_m.o prec_m.o parm.h nlin_m.o morepbl_m.o kuocom.h dava_m.o const_phys.h arrays_m.o newmpar.h tkeeps.o kuocomb_m.o
convjlm.o : convjlm.f establ.h vvel.h tracers.h soil_m.o sigs_m.o prec_m.o parm.h nlin_m.o morepbl_m.o map_m.o liqwpar_m.o latlong_m.o kuocom.h dava_m.o const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o tkeeps.o
davies.o : davies.f dates.h sigs_m.o parm.h davb_m.o dava_m.o arrays_m.o newmpar.h cc_mpi.o 
diag_m.o : diag_m.f90 parm.h xyzinfo_m.o sigs_m.o newmpar.h cc_mpi.o 
depts.o : depts.f bigxy4_m.o xyzinfo_m.o vecsuv_m.o parm.h map_m.o indices_m.o const_phys.h newmpar.h cc_mpi.o 
drive.o : drive.f 
e1e288.o : e1e288.f tfcom.h kdacom_m.o tabcom.h radisw_m.o rdparm.h hcon.h newmpar.h 
e3v88.o : e3v88.f tabcom.h rdparm.h hcon.h newmpar.h 
esbda.o : esbda.f 
esibda.o : esibda.f 
esfsw_driver.o : esfsw_parameters.o rad_utilities.o
establ.o : establ.f 
extras.o : extras.f 
findij.o : findij.f newmpar.h 
findll.o : findll.f newmpar.h 
findnear.o : findnear.f 
fst88.o : fst88.f cldcom_m.o tfcom.h srccom.h kdacom_m.o lwout_m.o rdflux_m.o tabcom.h rnddta.h radisw_m.o rdparm.h hcon.h newmpar.h 
gas_tf.o : rad_utilities.o longwave_params.o
gettin.o : gettin.f savuvt_m.o arrays_m.o newmpar.h 
globpe.o : globpe.f mapproj.h establ.h xyzinfo_m.o xarrs.h vvel.h vecsuv_m.o trcom2.h tracers.h stime.h soilv.h soilsnow.h soil_m.o sigs_m.o screen_m.o scamdim.h savuvt_m.o raddiag_m.o prec_m.o pbl_m.o parmvert.h parm_nqg.h parmhor.h parmdyn.h parm.h nsibd_m.o nlin_m.o morepbl_m.o map_m.o liqwpar_m.o latlong_m.o kuocom.h indices_m.o histave_m.o filnames.h extraout_m.o dates.h darcdf.h const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o tracermodule.o timeseries.o seaesfrad.o betts1_m.o cldcom_m.o co2dta_m.o davb_m.o gdrag_m.o kdacom_m.o lwout_m.o o3amip_m.o rdflux_m.o
gwdrag.o : gwdrag.f soil_m.o sigs_m.o pbl_m.o parm.h morepbl_m.o nlin_m.o gdrag_m.o const_phys.h arrays_m.o newmpar.h 
hconst.o : hconst.f hcon.h 
helmsol.o : helmsol.f90 parmdyn.h parm.h indices_m.o newmpar.h sumdd_m.o ilu_m.o cc_mpi.o 
helmsor.o : helmsor.f parmdyn.h parm.h indices_m.o newmpar.h cc_mpi.o 
hordifg.o : hordifg.f vecsuv_m.o sigs_m.o parm.h nlin_m.o map_m.o indices_m.o const_phys.h arrays_m.o newmpar.h cc_mpi.o tkeeps.o
hs_phys.o : hs_phys.f sigs_m.o parm.h nlin_m.o latlong_m.o arrays_m.o newmpar.h 
iabsdate.o : iabsdate.f 
icefall.o : icefall.f params.h parm.h morepbl_m.o kuocom.h cparams.h const_phys.h newmpar.h cc_mpi.o
ilu_m.o : ilu_m.f90 indices_m.o newmpar.h cc_mpi.o 
indata.o : indata.f vecsuv_m.o xyzinfo_m.o vecs.h trcom2.h tracers.h stime.h soilv.h soilsnow.h soil_m.o sigs_m.o prec_m.o permsurf_m.o pbl_m.o parm_nqg.h parmdyn.h parm.h nsibd_m.o morepbl_m.o map_m.o liqwpar_m.o latlong_m.o indices_m.o gdrag_m.o filnames.h dava_m.o dates.h const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o tracermodule.o timeseries.o ateb.o cable_ccam2.o mlo.o tkeeps.o
infile.o : infile.f sigs_m.o tracers.h stime.h parm_nqg.h parm.h liqwpar_m.o kuocom.h darcdf.h newmpar.h diag_m.o cc_mpi.o mlo.o ateb.o tkeeps.o
int2.o : int2.f newmpar.h 
ints.o : ints.f indices_m.o parmhor.h parm.h newmpar.h cc_mpi.o 
jimcc.o : jimcc.f bigxy4_m.o parm.h  
latltoij.o : latltoij.f parmdyn.h parm.h const_phys.h bigxy4_m.o newmpar.h utilities.o 
leoncld.o : leoncld.f establ.h vvel.h tracers.h soil_m.o sigs_m.o prec_m.o parm.h nlin_m.o morepbl_m.o map_m.o latlong_m.o kuocom.h dava_m.o arrays_m.o cparams.h const_phys.h liqwpar_m.o newmpar.h cc_mpi.o diag_m.o 
longwave_clouds.o : rad_utilities.o
longwave_fluxes.o : rad_utilities.o
longwave_tables.o : rad_utilities.o longwave_params.o
lw_gases_stdtf.o : rad_utilities.o gas_tf.o
lwr88.o : lwr88.f tfcom.h rnddta.h kdacom_m.o co2dta_m.o radisw_m.o rdparm.h parm.h hcon.h newmpar.h 
microphys_rad.o : rad_utilities.o longwave_params.o esfsw_parameters.o
mpidummy.o : mpidummy.f90 
mslp.o : mslp.f sigs_m.o parm.h const_phys.h newmpar.h cc_mpi.o 
mtimerget.o : mtimerget.f 
nestin.o : nestin.f stime.h soilsnow.h soil_m.o sigs_m.o pbl_m.o parm.h map_m.o davb_m.o dava_m.o dates.h const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o cable_define_dimensions.o
newcloud.o : newcloud.f sigs_m.o parm.h params.h kuocom.h cparams.h const_phys.h newmpar.h 
newrain.o : newrain.f params.h morepbl_m.o kuocom.h cparams.h const_phys.h newmpar.h 
nonlin.o : nonlin.f xyzinfo_m.o xarrs.h vvel.h vecsuv_m.o tracers.h sigs_m.o savuvt_m.o parmvert.h parmdyn.h parm.h nlin_m.o morepbl_m.o map_m.o latlong_m.o liqwpar_m.o kuocom.h indices_m.o const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o tkeeps.o
o3_read.o : o3_read.f newmpar.h const_phys.h
onthefly.o : onthefly.f indices_m.o indices_m.o xyzinfo_m.o vvel.h vecsuv_m.o tracers.h stime.h sigs_m.o parm_nqg.h parm.h map_m.o latlong_m.o const_phys.h bigxy4_m.o newmpar.h utilities.o cc_mpi.o mlo.o tkeeps.o ateb.o
optical_path.o : rad_utilities.o longwave_params.o lw_gases_stdtf.o
outcdf.o : outcdf.f vvel.h version.h trcom2.h soilv.h soilsnow.h soil_m.o sigs_m.o screen_m.o scamdim.h raddiag_m.o prec_m.o pbl_m.o nsibd_m.o morepbl_m.o mapproj.h map_m.o histave_m.o extraout_m.o arrays_m.o tracers.h parmvert.h parmhor.h parmdyn.h parm.h liqwpar_m.o kuocom.h filnames.h dates.h darcdf.h newmpar.h cc_mpi.o ateb.o mlo.o tracermodule.o tkeeps.o
outfile.o : outfile.f vvel.h tracers.h soilsnow.h soilv.h soil_m.o sigs_m.o screen_m.o scamdim.h prec_m.o pbl_m.o parmvert.h parmdyn.h parm.h nsibd_m.o nlin_m.o morepbl_m.o map_m.o kuocom.h histave_m.o filnames.h extraout_m.o dava_m.o dates.h darcdf.h arrays_m.o newmpar.h cc_mpi.o 
pbldif.o : pbldif.f map_m.o sigs_m.o parm.h morepbl_m.o kuocom.h extraout_m.o const_phys.h arrays_m.o newmpar.h 
radriv90.o : radriv90.f establ.h tfcom.h swocom.h srccom.h rdflux_m.o raddiag_m.o radisw_m.o lwout_m.o hcon.h cldcom_m.o rdparm.h soilv.h soilsnow.h soil_m.o sigs_m.o scamdim.h pbl_m.o parm.h nsibd_m.o map_m.o liqwpar_m.o latlong_m.o kuocom.h extraout_m.o dates.h cparams.h const_phys.h arrays_m.o newmpar.h zenith.o swr99.o ateb.o mlo.o o3_read.o
rdparm.o : rdparm.f 
read_ht.o : read_ht.f 
resetd.o : resetd.f 
retopo.o : retopo.f sigs_m.o parm.h const_phys.h newmpar.h cc_mpi.o 
rnddta.o : rnddta.f 
scamrdn.o : scamrdn.f soilv.h soilsnow.h soil_m.o sigs_m.o scamdim.h pbl_m.o parm.h nsibd_m.o filnames.h const_phys.h arrays_m.o newmpar.h 
scrnout.o : scrnout.f establ.h soilsnow.h soil_m.o sigs_m.o scamdim.h prec_m.o pbl_m.o parm.h nsibd_m.o map_m.o liqwpar_m.o const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o morepbl_m.o
seaesfrad.o : rad_utilities.o microphys_rad.o esfsw_driver.o sealw99.o esfsw_parameters.o zenith.o ateb.o cable_ccam2.o mlo.o o3_read.o radisw_m.o
sealw99.o : rad_utilities.o longwave_params.o longwave_clouds.o longwave_fluxes.o longwave_tables.o optical_path.o gas_tf.o lw_gases_stdtf.o
setxyz.o : setxyz.f bigxy4_m.o indices_m.o vecsuv_m.o xyzinfo_m.o parm.h map_m.o latlong_m.o const_phys.h  utilities.o 
sflux.o : sflux.f latlong_m.o dates.h establ.h vvel.h trcom2.h tracers.h soilsnow.h soilv.h soil_m.o sigs_m.o screen_m.o scamdim.h savuvt_m.o prec_m.o permsurf_m.o pbl_m.o parm.h nsibd_m.o morepbl_m.o map_m.o liqwpar_m.o gdrag_m.o extraout_m.o const_phys.h arrays_m.o newmpar.h cc_mpi.o diag_m.o ateb.o cable_ccam2.o mlo.o
soilsnow.o : soilsnow.f nlin_m.o soil_m.o sigs_m.o arrays_m.o morepbl_m.o nsibd_m.o soilv.h soilsnow.h permsurf_m.o parm.h const_phys.h newmpar.h diag_m.o cc_mpi.o 
solargh.o : solargh.f 
spa88.o : spa88.f lwout_m.o cldcom_m.o tfcom.h kdacom_m.o srccom.h rdflux_m.o rnddta.h radisw_m.o rdparm.h hcon.h newmpar.h 
srccom.o : srccom.f 
sscam2.o : sscam2.f soilv.h soilsnow.h scamdim.h parm.h newmpar.h 
staguv.o : staguv.f vecsuv_m.o parmdyn.h parm.h map_m.o indices_m.o newmpar.h cc_mpi.o 
sumdd_m.o : sumdd_m.f90 
swr99.o : swr99.f rdparm.h hcon.h newmpar.h 
table.o : table.f tabcom.h radisw_m.o hcon.h rnddta.h rdparm.h newmpar.h 
timeseries.o: timeseries.f dates.h newmpar.h tracers.h extraout_m.o arrays_m.o soil_m.o prec_m.o vvel.h pbl_m.o morepbl_m.o soilsnow.h nsibd_m.o sigs_m.o tracermodule.o cable_define_dimensions.o carbpools_m.o
tracermodule.o: tracermodule.f newmpar.h tracers.h parm.h const_phys.h arrays_m.o sigs_m.o xyzinfo_m.o
trim.o : trim.f newmpar.h 
trvmix.o: trvmix.f newmpar.h const_phys.h parm.h sigs_m.o tracers.h arrays_m.o tracermodule.o
updps.o : updps.f xarrs.h vvel.h sigs_m.o parmdyn.h parm.h parmhor.h savuvt_m.o nlin_m.o map_m.o indices_m.o const_phys.h arrays_m.o newmpar.h cc_mpi.o 
upglobal.o : upglobal.f xyzinfo_m.o xarrs.h vvel.h vecsuv_m.o tracers.h sigs_m.o parmvert.h parmhor.h parmdyn.h parm.h nlin_m.o map_m.o liqwpar_m.o kuocom.h indices_m.o const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o tkeeps.o
utilities.o : utilities.f90 const_phys.h 
vadv30.o : vadv30.f xarrs.h vvel.h tracers.h sigs_m.o parmvert.h parmdyn.h parm.h map_m.o liqwpar_m.o kuocom.h indices_m.o arrays_m.o newmpar.h cc_mpi.o tkeeps.o
vadvtvd.o : vadvtvd.f xarrs.h vvel.h tracers.h sigs_m.o parmvert.h parmdyn.h parm.h map_m.o liqwpar_m.o kuocom.h arrays_m.o newmpar.h cc_mpi.o tkeeps.o
vertmix.o : vertmix.f nsibd_m.o establ.h tracers.h soil_m.o sigs_m.o savuvt_m.o permsurf_m.o screen_m.o pbl_m.o parm.h morepbl_m.o nlin_m.o map_m.o liqwpar_m.o kuocom.h indices_m.o dates.h const_phys.h arrays_m.o newmpar.h diag_m.o cc_mpi.o trvmix.o tkeeps.o
zenith.o : zenith.f90 
