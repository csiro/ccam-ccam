FC = mpif90
FCSCM = ifort
CC = cc

# Common compiler flags
ifneq ($(CUSTOM),yes)
NCFLAG = -I $(NETCDF_ROOT)/include
ifeq ($(NCCLIB),yes)
NCFLAG += -Dncclib
endif
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dusempi3
endif
FHOST = -O3 -xHost
FOVERRIDE =
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH =
MPISPECIAL =
ifeq ($(XEONPHI),yes)
FHOST = -O3 -xMIC-AVX512
endif
ifeq ($(BROADWELL),yes)
FHOST = -O3 -xCORE-AVX2
IPFLAG = -ip
IPOFLAG = -ipo
VTHRESH = -vec-threshold0
endif
ifeq ($(CASCADELAKE),yes)
FHOST = -O3 -xCASCADELAKE -fimf-use-svml
FOVERRIDE = -qoverride-limits
ZMM = -qopt-zmm-usage=high
IPFLAG = -ip
VTHRESH = -vec-threshold0
endif
# OpenMP compile flag
ifeq ($(OMP),yes)
OMPFLAG = -qopenmp -qno-openmp-simd
else
OMPFLAG =
endif
# Default intel compiler options
FFLAGS = $(FHOST) -ftz -fp-model precise -no-fma -traceback $(MPIFLAG) $(NCFLAG) $(OMPFLAG)
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf
ifneq ($(NCCLIB),yes)
LIBS += -lnetcdff
endif
PPFLAG90 = -fpp
PPFLAG77 = -fpp
PPFLAG90F = -fpp
REAL8FLAG = -r8
INT8FLAG = -i8
DEBUGFLAG = -check all -debug all -fpe0
endif

# Gfortran compiler options
ifeq ($(GFORTRAN),yes)
MPIFC = gfortran
MPIF77 = gfortran
FC = mpif90
FCSCM = gfortran
FHOST = -march=native
MPIFLAG =
#MPISPECIAL = -fallow-argument-mismatch
MPISPECIAL =
FFLAGS = -O3 -mtune=native -mveclibabi=svml $(FHOST) -fbacktrace $(MPIFLAG) $(NCFLAG)
FOVERRIDE =
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH =
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
PPFLAG90F =
REAL8FLAG = -fdefault-real-8
INT8FLAG = -fdefault-int-8
DEBUGFLAG = -g -Wall -Wextra -fbounds-check
endif

#PGFORTRAN
ifeq ($(PGI),yes)
MPIFC = pgfortran
MPIF77 = pgfortran
FC = pgfortran -I$(I_MPI_ROOT)/include64 -L$(I_MPI_ROOT)/lib64 -lmpi -lmpiif 
FCSCM = pgfortran
CC = pgcc
FHOST = -O3 -tp=haswell -fast
#FHOST = -g
MPIFLAG +=  
MPISPECIAL =
FFLAGS = $(FHOST) -Dpgi -DGPU -traceback $(MPIFLAG) $(NCFLAG)
#FFLAGS = $(FHOST) -Dpgi -D_GPU -traceback $(MPIFLAG) $(NCFLAG)
#FFLAGS += -Minfo=accel -acc -ta=host
#FFLAGS += -Minfo=accel -acc -ta=multicore
#FFLAGS += -Minfo=accel -acc -ta=nvidia:cc60
FFLAGS += -Minfo=accel -acc -ta=tesla:cc60
FOVERRIDE =
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH =
PPFLAG90 = -cpp
PPFLAG77 = -cpp
PPFLAG90F = -cpp
REAL8FLAG = -r8
INT8FLAG = -i8
DEBUGFLAG =
endif

# CRAY compiler options
ifeq ($(CRAY),yes)
FC = ftn
FCSCM = ftn
FFLAGS = -h noomp -Dusenc_mod
FOVERRIDE =
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH =
MPISPECIAL =
PPFLAG90 = -eZ
PPFLAG77 = -eZ
PPFLAG90F = -eZ
REAL8FLAG = -s real64
INT8FLAG = -s integer64
DEBUGFLAG =
endif

# IBM compiler options
#ifeq ($(IBM),yes)
#FC = xlf
#FCSCM = xlf
#CC = xlc
#FFLAGS = -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -qextname $(MPIFLAG)
#LIBS = -L $(NETCDF_ROOT)/lib -L /opt/ibmhpc/pecurrent/mpich/xlf/lib64 -lnetcdf -lnetcdff -lmpi -lmpigf
#MPIFLAG = -I /opt/ibmhpc/pecurrent/mpich/xlf/include64
#PPFLAG90 = -qsuffix=ccp=f90
#PPFLAG77 = -qsuffix=ccp=f
#PPFLAG90F = qsuffix=ccp=F90
#REAL8FLAG = -qrealsize=8
#INT8FLAG = -qintsize=8
#DEBUGFLAG = -q
#endif

# Options for building with VAMPIRTrace
ifeq ($(VT),yes)
FC = vtfort -vt:fc mpif90 -vt:inst manual
FFLAGS += -Dvampir -DVTRACE
else
FFLAGS += -Dsimple_timer
endif

# Testing - I/O and fpmodel
ifeq ($(TEST),yes)
FFLAGS += $(DEBUGFLAG)
endif

# Build with 64 ints/reals
ifeq ($(I8R8),yes)
FFLAGS += $(REAL8FLAG) $(INT8FLAG) -Di8r8
endif

# Use NetCDF F90 interface
ifeq ($(NCMOD),yes)
FFLAGS += -Dusenc_mod
endif

# Use Netcdf3
ifeq ($(NETCDF3),yes)
FFLAGS += -Dusenc3
endif

# CABLE, MLO, aTEB, etc
FFLAGS += -DCCAM

# Object files for dynamical model
OBJS = adjust5.o amipsst.o convjlm.o convjlm22.o depts.o estab.o gettin.o \
globpe.o gdrag_m.o hordifg.o hs_phys.o indata.o infile.o ints.o \
helmsolve.o jimcc.o nesting.o nonlin.o ensemble.o \
outcdf.o radriv90.o scrnout.o setxyz.o sflux.o \
soilsnow.o staguv.o upglobal.o eig.o updps.o vadvtvd.o \
vertmix.o leoncld.o cloudmod.o latltoij.o \
cldblk.o clddia.o clo89.o cloud.o cloud2.o co2_read.o e1e288.o \
e3v88.o fst88.o hconst.o lwr88.o ozoneread.o spa88.o \
swr99.o table.o zenith.o cc_acc.o cc_mpi.o cc_omp.o diag_m.o sumdd_m.o daviesnudge.o \
utilities.o onthefly.o tracermodule.o timeseries.o \
trvmix.o getopt_m.o usage_m.o const_phys.o \
betts.o bett_cuc.o bettinit.o bettrain.o bettspli.o \
xyzinfo_m.o vecsuv_m.o map_m.o latlong_m.o indices_m.o bigxy4_m.o \
arrays_m.o betts1_m.o carbpools_m.o cldcom_m.o co2dta_m.o cfrac_m.o \
dpsdt_m.o epst_m.o extraout_m.o histave_m.o kdacom_m.o \
kuocomb_m.o liqwpar_m.o lwout_m.o morepbl_m.o nharrs_m.o \
nlin_m.o nsibd_m.o parmhdff_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o \
radisw_m.o rdflux_m.o riverarrays_m.o savuvt_m.o savuv1_m.o sbar_m.o screen_m.o \
sigs_m.o soil_m.o soilsnow_m.o srccom_m.o swocom_m.o tabcom_m.o \
tbar2d_m.o tfcom_m.o tracers_m.o unn_m.o uvbar_m.o vecs_m.o vegpar_m.o vvel_m.o \
workglob_m.o work2_m.o work3_m.o work3b_m.o work3f_m.o work3lwr_m.o work3sav_m.o \
xarrs_m.o \
aerointerface.o aerosolldr.o \
cable_air.o cable_albedo.o cable_canopy.o cable_ccam2.o cable_ccam3.o cable_ccam4.o cable_common.o \
cable_data.o cable_define_types.o cable_radiation.o cable_roughness.o cable_soilsnow.o \
cable_optimiseJVratio.o cable_psm.o cable_pft_params.o cable_soil_params.o cable_constants.o \
cable_gw_hydro.o \
cable_sli_main.o cable_sli_numbers.o cable_sli_roots.o cable_sli_solve.o cable_sli_utils.o \
casa_cnp.o casa_variable.o POP.o casa_dimension.o casa_phenology.o casa_param.o \
ateb.o mlo.o river.o tkeeps.o \
seaesfrad.o rad_utilities.o microphys_rad.o esfsw_driver.o esfsw_parameters.o \
longwave_params.o sealw99.o longwave_clouds.o longwave_fluxes.o longwave_tables.o \
optical_path.o gas_tf.o lw_gases_stdtf.o \
mlodynamics.o mlodynamicsarrays_m.o mlodiffg.o mlostag.o mlodepts.o mloints.o \
darcdf_m.o dates_m.o filnames_m.o newmpar_m.o parm_m.o parmdyn_m.o parmgeom_m.o \
parmhor_m.o soilv_m.o stime_m.o \
netcdf_m.o mpif_m.o stacklimit.o

# Object files for single column mode
OBJSCM = aerointerface.o aerosolldr.o arrays_m.o ateb.o cable_air.o cable_albedo.o \
cable_canopy.o cable_ccam2.o cable_ccam3.o cable_ccam4.o cable_common.o cable_data.o \
cable_define_types.o cable_radiation.o cable_roughness.o cable_soilsnow.o carbpools_m.o \
cable_optimiseJVratio.o cable_psm.o cable_pft_params.o cable_soil_params.o cable_constants.o \
cable_gw_hydro.o \
cable_sli_main.o cable_sli_numbers.o cable_sli_roots.o cable_sli_solve.o cable_sli_utils.o \
casa_cnp.o casa_variable.o POP.o casa_dimension.o casa_phenology.o casa_param.o \
cc_mpi.o cc_omp.o cfrac_m.o cloudmod.o co2_read.o co2dta_m.o  \
const_phys.o convjlm.o convjlm22.o darcdf_m.o dates_m.o diag_m.o esfsw_driver.o esfsw_parameters.o estab.o \
extraout_m.o filnames_m.o gas_tf.o gdrag_m.o histave_m.o indices_m.o infile.o \
kuocomb_m.o latlong_m.o leoncld.o liqwpar_m.o longwave_clouds.o longwave_fluxes.o \
longwave_params.o longwave_tables.o lw_gases_stdtf.o map_m.o microphys_rad.o mlo.o \
mlodynamicsarrays_m.o morepbl_m.o netcdf_m.o newmpar_m.o nharrs_m.o nsibd_m.o \
optical_path.o ozoneread.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o \
pbl_m.o permsurf_m.o prec_m.o rad_utilities.o raddiag_m.o radisw_m.o \
riverarrays_m.o savuvt_m.o scm.o scmarrays_m.o screen_m.o scrnout.o \
seaesfrad.o sealw99.o sflux.o sigs_m.o soil_m.o soilsnow.o soilsnow_m.o soilv_m.o \
stime_m.o tkeeps.o tracers_m.o vecsuv_m.o vegpar_m.o vertmix.o vvel_m.o work2_m.o \
work3_m.o work3b_m.o work3f_m.o xyzinfo_m.o zenith.o \
getopt_m.o stacklimit.o

ifeq ($(SCM),yes)
FC = $(FCSCM)
FFLAGS += -Dscm
scm: $(OBJSCM)
	$(FC) -o scm $(FFLAGS) $(OBJSCM) $(LIBS)
else
globpea: $(OBJS)
	$(FC) -o globpea $(FFLAGS) $(OBJS) $(LIBS)
endif

clean:
	rm *.o *.i *.mod

.SUFFIXES:.f90 .F90

netcdf_m.o: netcdf_m.f90
	$(FC) -c $(PPFLAG90) $(NCFLAG) $<
mpif_m.o: mpif_m.f90
	$(FC) -c $(PPFLAG90) $(MPIFLAG) $(MPISPECIAL) $<
esfsw_driver.o: esfsw_driver.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(VTHRESH) $<
esfsw_parameters.o: esfsw_parameters.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
gas_tf.o: gas_tf.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
longwave_clouds.o: longwave_clouds.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
longwave_fluxes.o: longwave_fluxes.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
longwave_tables.o: longwave_tables.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
longwave_params.o: longwave_params.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
lw_gases_stdtf.o: lw_gases_stdtf.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
microphys_rad.o: microphys_rad.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
optical_path.o: optical_path.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
rad_utilities.o: rad_utilities.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
sealw99.o: sealw99.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_air.o: cable_air.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_albedo.o: cable_albedo.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_canopy.o: cable_canopy.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_ccam3.o: cable_ccam3.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_common.o: cable_common.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_constants.o: cable_constants.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_data.o: cable_data.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_define_types.o: cable_define_types.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_gw_hydro.o: cable_gw_hydro.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_optimiseJVratio.o: cable_optimiseJVratio.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_pft_params.o: cable_pft_params.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_psm.o: cable_psm.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_radiation.o: cable_radiation.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_roughness.o: cable_roughness.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_sli_main.o: cable_sli_main.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_sli_numbers.o: cable_sli_numbers.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_sli_roots.o: cable_sli_roots.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_sli_solve.o: cable_sli_solve.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_sli_utils.o: cable_sli_utils.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_soil_params.o: cable_soil_params.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
cable_soilsnow.o: cable_soilsnow.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
casa_cnp.o: casa_cnp.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
casa_dimension.o: casa_dimension.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
casa_param.o: casa_param.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
casa_phenology.o: casa_phenology.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
casa_variable.o: casa_variable.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
POP.o: POP.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $(IPFLAG) $<
estab.o: estab.f90
	$(FC) -c $(FFLAGS) $(IPOFLAG) $(PPFLAG90) $<
helmsolve.o: helmsolve.f90
	$(FC) -c $(PPFLAG90) $(FFLAGS) $(FOVERRIDE) $<
ints.o: ints.f90
	$(FC) -c $(FFLAGS) $(IPFLAG) $(ZMM) $(PPFLAG90) $<
leoncld.o: leoncld.f90
	$(FC) -c $(FFLAGS) $(IPOFLAG) $(PPFLAG90) $<
seaesfrad.o: seaesfrad.f90
	$(FC) -c $(FFLAGS) $(VTHRESH) $(PPFLAG90) $<
stacklimit.o: stacklimit.c
	$(CC) -c stacklimit.c
tkeeps.o: tkeeps.f90
	$(FC) -c $(FFLAGS) $(PPFLAG90) $<
vertmix.o: vertmix.f90
	$(FC) -c $(FFLAGS) $(IPOFLAG) $(PPFLAG90) $<
version.h: FORCE
	rm -f brokenver tmpver
	echo "      character(len=*), parameter :: version ='CCAM r'" > brokenver
	echo "      character(len=*), parameter :: version ='CCAM r`svnversion .`'" > tmpver
	grep exported tmpver || grep Unversioned tmpver || cmp tmpver brokenver || cmp tmpver version.h || mv tmpver version.h
FORCE:


.f90.o:
	$(FC) -c $(FFLAGS) $(IPFLAG) $(PPFLAG90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(IPFLAG) $(PPFLAG90F) $<	
.f.o:
	$(FC) -c $(FFLAGS) $(IPFLAG) $(PPFLAG77) $<

# Remove mod rule from Modula 2 so GNU make doesn't get confused
%.o : %.mod

# Dependencies
adjust5.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o dpsdt_m.o epst_m.o helmsolve.o indices_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o pbl_m.o sigs_m.o staguv.o tbar2d_m.o tracers_m.o vadvtvd.o vecsuv_m.o vecs_m.o vvel_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o kuocom.h
aerointerface.o : aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o cfrac_m.o cloudmod.o const_phys.o extraout_m.o infile.o kuocomb_m.o latlong_m.o liqwpar_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o ozoneread.o parm_m.o parmgeom_m.o pbl_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o tkeeps.o vegpar_m.o work2_m.o zenith.o kuocom.h
amipsst.o : arrays_m.o cc_mpi.o const_phys.o dates_m.o filnames_m.o infile.o latlong_m.o latltoij.o mlo.o newmpar_m.o nharrs_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o setxyz.o soil_m.o soilsnow_m.o workglob_m.o
ateb.o : cc_omp.o
bett_cuc.o : betts1_m.o newmpar_m.o
bettinit.o : betts1_m.o newmpar_m.o
bettrain.o : betts1_m.o newmpar_m.o
betts.o : betts1_m.o morepbl_m.o newmpar_m.o parm_m.o prec_m.o sigs_m.o
cable_air.o : cable_common.o cable_data.o cable_define_types.o 
cable_albedo.o : cable_common.o cable_data.o cable_define_types.o
cable_canopy.o : cable_air.o cable_common.o cable_data.o cable_define_types.o cable_gw_hydro.o cable_psm.o cable_radiation.o cable_roughness.o cable_sli_main.o cable_sli_utils.o
cable_common.o : cable_define_types.o cable_pft_params.o cable_soil_params.o
cable_ccam2.o : arrays_m.o cable_air.o cable_albedo.o cable_canopy.o cable_ccam3.o cable_ccam4.o cable_common.o cable_define_types.o cable_optimiseJVratio.o cable_radiation.o cable_roughness.o cable_sli_main.o cable_soilsnow.o carbpools_m.o casa_cnp.o casa_variable.o casa_phenology.o cc_mpi.o cc_omp.o const_phys.o darcdf_m.o dates_m.o estab.o extraout_m.o infile.o latlong_m.o liqwpar_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o POP.o prec_m.o raddiag_m.o radisw_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o tracers_m.o vegpar_m.o work2_m.o work3_m.o zenith.o
cable_ccam3.o : cable_data.o cable_define_types.o casa_phenology.o casa_variable.o POP.o
cable_ccam4.o : cable_ccam3.o cable_data.o cable_define_types.o casa_variable.o cc_omp.o newmpar_m.o POP.o
cable_constants.o : cable_define_types.o
cable_data.o : cable_constants.o
cable_gw_hydro.o : cable_define_types.o cable_common.o cable_soilsnow.o cable_data.o
cable_optimiseJVratio.o : cable_canopy.o cable_data.o cable_define_types.o POP.o
cable_pft_parms.o : cable_define_types.o
cable_psm.o : cable_common.o cable_define_types.o
cable_radiation.o : cable_common.o cable_data.o cable_define_types.o
cable_roughness.o : cable_common.o cable_data.o cable_define_types.o
cable_sli_main.o : cable_common.o cable_define_types.o cable_sli_numbers.o cable_sli_roots.o cable_sli_solve.o cable_sli_utils.o
cable_sli_numbers.o : cable_define_types.o
cable_sli_roots.o : cable_define_types.o cable_sli_numbers.o
cable_sli_solve.o : cable_define_types.o cable_sli_numbers.o cable_sli_utils.o
cable_sli_utils.o : cable_define_types.o cable_sli_numbers.o
cable_soil_parms.o : cable_define_types.o
cable_soilsnow.o : cable_common.o cable_data.o cable_define_types.o parm_m.o
carbpools_m.o : cable_define_types.o casa_variable.o parm_m.o
casa_cnp.o : cable_define_types.o casa_variable.o
casa_dimension.o : cable_define_types.o
casa_param.o : casa_dimension.o
cable_pft_params.o : cable_define_types.o
casa_phenology.o : casa_dimension.o
cable_soil_params.o : cable_define_types.o
casa_variable.o : cable_define_types.o casa_dimension.o casa_param.o
cc_mpi.o : cc_omp.o const_phys.o indices_m.o latlong_m.o map_m.o newmpar_m.o sumdd_m.o vecsuv_m.o workglob_m.o xyzinfo_m.o
ifneq ($(SCM),yes)
cc_mpi.o : mpif_m.o
endif
cc_omp.o : newmpar_m.o
clddia.o : arrays_m.o cc_mpi.o const_phys.o map_m.o morepbl_m.o newmpar_m.o parm_m.o pbl_m.o sigs_m.o soil_m.o vvel_m.o kuocom.h
clo89.o : cldcom_m.o newmpar_m.o parm_m.o radisw_m.o rdparm.h
cloud2.o : diag_m.o cc_mpi.o const_phys.o leoncld.o newmpar_m.o parm_m.o radisw_m.o sigs_m.o hcon.h kuocom.h rdparm.h
cloud.o : extraout_m.o newmpar_m.o parm_m.o radisw_m.o rdparm.h
cloudmod.o : const_phys.o newmpar_m.o parm_m.o sigs_m.o kuocom.h
co2_read.o : cc_mpi.o co2dta_m.o filnames_m.o infile.o newmpar_m.o parm_m.o radisw_m.o stime_m.o rdparm.h
convjlm.o : aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocomb_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o prec_m.o sigs_m.o soil_m.o tkeeps.o tracers_m.o vvel_m.o work2_m.o kuocom.h
convjlm22.o : aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocomb_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o prec_m.o sigs_m.o soil_m.o tkeeps.o tracers_m.o vvel_m.o work2_m.o kuocom.h
daviesnudge.o : aerosolldr.o arrays_m.o cc_mpi.o newmpar_m.o parm_m.o sigs_m.o
depts.o : bigxy4_m.o cc_mpi.o const_phys.o indices_m.o map_m.o newmpar_m.o parm_m.o parmhor_m.o parmgeom_m.o uvbar_m.o vecsuv_m.o work3f_m.o xyzinfo_m.o 
diag_m.o : cc_mpi.o newmpar_m.o parm_m.o sigs_m.o sumdd_m.o xyzinfo_m.o
e1e288.o : kdacom_m.o newmpar_m.o radisw_m.o tabcom_m.o tfcom_m.o hcon.h rdparm.h
e3v88.o : newmpar_m.o tabcom_m.o hcon.h rdparm.h
eig.o : cc_mpi.o const_phys.o newmpar_m.o vecs_m.o
ensemble.o : aerosolldr.o arrays_m.o cc_mpi.o dates_m.o filnames_m.o mlo.o newmpar_m.o onthefly.o outcdf.o parm_m.o pbl_m.o soil_m.o soilsnow_m.o xyzinfo_m.o
esfsw_driver.o : esfsw_parameters.o rad_utilities.o
esfsw_parameters.o : rad_utilities.o
estab.o : const_phys.o
fst88.o : cc_mpi.o cldcom_m.o diag_m.o kdacom_m.o lwout_m.o newmpar_m.o parm_m.o radisw_m.o rdflux_m.o srccom_m.o tabcom_m.o tfcom_m.o hcon.h rdparm.h rnddta.h
gas_tf.o : longwave_params.o rad_utilities.o
gdrag_m.o : arrays_m.o cc_mpi.o cc_omp.o const_phys.o liqwpar_m.o newmpar_m.o nharrs_m.o parm_m.o pbl_m.o sigs_m.o
gettin.o : arrays_m.o newmpar_m.o savuvt_m.o
globpe.o : aerointerface.o aerosolldr.o amipsst.o arrays_m.o ateb.o bigxy4_m.o cable_ccam2.o carbpools_m.o cc_acc.o cc_mpi.o cc_omp.o cfrac_m.o cloudmod.o const_phys.o convjlm.o convjlm22.o darcdf_m.o dates_m.o daviesnudge.o diag_m.o dpsdt_m.o ensemble.o epst_m.o estab.o extraout_m.o filnames_m.o gdrag_m.o getopt_m.o histave_m.o hs_phys.o indata.o indices_m.o infile.o kuocomb_m.o latlong_m.o leoncld.o liqwpar_m.o map_m.o mlo.o mlodiffg.o mlodynamics.o morepbl_m.o nesting.o newmpar_m.o nharrs_m.o nlin_m.o nsibd_m.o outcdf.o ozoneread.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o river.o savuvt_m.o savuv1_m.o sbar_m.o screen_m.o seaesfrad.o setxyz.o sflux.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o stime_m.o tbar2d_m.o timeseries.o tkeeps.o tracermodule.o tracers_m.o unn_m.o usage_m.o uvbar_m.o vecs_m.o vecsuv_m.o vegpar_m.o vertmix.o vvel_m.o workglob_m.o work2_m.o work3_m.o work3f_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o kuocom.h version.h
hconst.o : hcon.h
helmsolve.o : cc_mpi.o diag_m.o indices_m.o newmpar_m.o parm_m.o parmdyn_m.o parmgeom_m.o sumdd_m.o vecs_m.o
histave_m.o: parm_m.o
hordifg.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o dpsdt_m.o indices_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o parmhdff_m.o savuvt_m.o sigs_m.o tkeeps.o vecsuv_m.o vvel_m.o kuocom.h
hs_phys.o : arrays_m.o cc_omp.o latlong_m.o newmpar_m.o nlin_m.o parm_m.o sigs_m.o
indata.o : aerointerface.o aerosolldr.o amipsst.o arrays_m.o ateb.o bigxy4_m.o cable_ccam2.o cc_mpi.o const_phys.o convjlm.o convjlm22.o darcdf_m.o dates_m.o daviesnudge.o diag_m.o ensemble.o epst_m.o extraout_m.o filnames_m.o gdrag_m.o indices_m.o infile.o latlong_m.o latltoij.o liqwpar_m.o map_m.o mlo.o mlodynamics.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o onthefly.o parm_m.o parmdyn_m.o parmgeom_m.o pbl_m.o permsurf_m.o river.o seaesfrad.o sflux.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o stime_m.o timeseries.o tracermodule.o tracers_m.o vecs_m.o vecsuv_m.o vegpar_m.o vertmix.o xyzinfo_m.o kuocom.h
indices_m.o : newmpar_m.o
infile.o : cc_mpi.o dates_m.o netcdf_m.o newmpar_m.o parm_m.o parmgeom_m.o sigs_m.o
ints.o : cc_mpi.o indices_m.o newmpar_m.o parm_m.o parmhor_m.o
latltoij.o : const_phys.o newmpar_m.o parm_m.o parmdyn_m.o utilities.o 
leoncld.o : aerointerface.o aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o cfrac_m.o cloudmod.o const_phys.o kuocomb_m.o latlong_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o prec_m.o sigs_m.o soil_m.o work3f_m.o kuocom.h
longwave_clouds.o : rad_utilities.o
longwave_fluxes.o : rad_utilities.o
longwave_tables.o : longwave_params.o rad_utilities.o
lw_gases_stdtf.o : cc_mpi.o filnames_m.o infile.o gas_tf.o newmpar_m.o rad_utilities.o
lwr88.o : co2dta_m.o kdacom_m.o newmpar_m.o parm_m.o radisw_m.o tfcom_m.o work3lwr_m.o hcon.h rdparm.h rnddta.h
microphys_rad.o : esfsw_parameters.o longwave_params.o rad_utilities.o
mlo.o : cc_omp.o
mlodepts.o : bigxy4_m.o cc_mpi.o const_phys.o indices_m.o mlo.o newmpar_m.o parm_m.o parmgeom_m.o parmhor_m.o vecsuv_m.o xyzinfo_m.o
mlodiffg.o : cc_mpi.o const_phys.o indices_m.o map_m.o mlo.o mlodynamicsarrays_m.o newmpar_m.o parm_m.o soil_m.o vecsuv_m.o
mlodynamics.o : arrays_m.o bigxy4_m.o cc_mpi.o cc_omp.o const_phys.o helmsolve.o indices_m.o infile.o latlong_m.o map_m.o mlo.o mlodepts.o mlodiffg.o mlodynamicsarrays_m.o mloints.o mlostag.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o parmgeom_m.o parmhor_m.o soil_m.o soilsnow_m.o soilv_m.o vecsuv_m.o xyzinfo_m.o 
mloints.o : cc_mpi.o indices_m.o mlo.o newmpar_m.o parm_m.o parmhor_m.o
mlostag.o : cc_mpi.o indices_m.o mlo.o mlodynamicsarrays_m.o newmpar_m.o parm_m.o
nesting.o : aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o const_phys.o dates_m.o daviesnudge.o darcdf_m.o diag_m.o indices_m.o latlong_m.o liqwpar_m.o map_m.o mlo.o mlodynamicsarrays_m.o newmpar_m.o nharrs_m.o onthefly.o parm_m.o parmdyn_m.o parmgeom_m.o pbl_m.o savuvt_m.o savuv1_m.o sigs_m.o soil_m.o soilsnow_m.o stime_m.o vecsuv_m.o work3sav_m.o xyzinfo_m.o kuocom.h
nonlin.o : aerosolldr.o arrays_m.o cc_mpi.o const_phys.o diag_m.o epst_m.o hs_phys.o indices_m.o latlong_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o savuvt_m.o sigs_m.o staguv.o tbar2d_m.o tkeeps.o tracers_m.o unn_m.o vadvtvd.o vecsuv_m.o vvel_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o kuocom.h
onthefly.o : aerointerface.o aerosolldr.o ateb.o cable_ccam2.o casa_variable.o carbpools_m.o cc_mpi.o const_phys.o darcdf_m.o diag_m.o extraout_m.o histave_m.o infile.o latlong_m.o latltoij.o mlo.o mlodynamics.o mlodynamicsarrays_m.o mlostag.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmdyn_m.o parmgeom_m.o prec_m.o raddiag_m.o riverarrays_m.o savuvt_m.o savuv1_m.o screen_m.o setxyz.o sigs_m.o soil_m.o soilv_m.o stime_m.o tkeeps.o tracers_m.o utilities.o vecsuv_m.o vvel_m.o workglob_m.o work2_m.o xarrs_m.o kuocom.h
optical_path.o : longwave_params.o lw_gases_stdtf.o rad_utilities.o
outcdf.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o cable_ccam2.o cable_define_types.o casa_variable.o carbpools_m.o cc_mpi.o cfrac_m.o const_phys.o dates_m.o daviesnudge.o dpsdt_m.o extraout_m.o filnames_m.o gdrag_m.o histave_m.o infile.o latlong_m.o liqwpar_m.o map_m.o mlo.o mlodiffg.o mlodynamics.o mlodynamicsarrays_m.o mlostag.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o ozoneread.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o pbl_m.o prec_m.o raddiag_m.o river.o riverarrays_m.o savuvt_m.o savuv1_m.o screen_m.o seaesfrad.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o tkeeps.o tracermodule.o tracers_m.o vegpar_m.o vvel_m.o work2_m.o xarrs_m.o kuocom.h version.h
ozoneread.o : cc_mpi.o const_phys.o dates_m.o filnames_m.o infile.o latlong_m.o newmpar_m.o parm_m.o 
radriv90.o : aerointerface.o arrays_m.o ateb.o cc_mpi.o cfrac_m.o cldcom_m.o co2_read.o co2dta_m.o const_phys.o diag_m.o estab.o extraout_m.o infile.o kdacom_m.o kuocomb_m.o latlong_m.o liqwpar_m.o lwout_m.o mlo.o newmpar_m.o nharrs_m.o nsibd_m.o ozoneread.o parm_m.o pbl_m.o raddiag_m.o radisw_m.o rdflux_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o srccom_m.o swocom_m.o swr99.o tabcom_m.o tfcom_m.o work3f_m.o work3lwr_m.o zenith.o kuocom.h rdparm.h hcon.h
river.o : arrays_m.o cable_ccam2.o cc_mpi.o const_phys.o indices_m.o map_m.o newmpar_m.o nsibd_m.o parm_m.o riverarrays_m.o soil_m.o soilsnow_m.o soilv_m.o xyzinfo_m.o 
riverarrays_m.o : newmpar_m.o
scm.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o cable_ccam2.o carbpools_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o convjlm.o convjlm22.o dates_m.o estab.o extraout_m.o filnames_m.o gdrag_m.o getopt_m.o histave_m.o hs_phys.o infile.o kuocomb_m.o latlong_m.o leoncld.o liqwpar_m.o map_m.o mlo.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o pbl_m.o prec_m.o raddiag_m.o riverarrays_m.o radisw_m.o savuvt_m.o scmarrays_m.o screen_m.o seaesfrad.o sflux.o sigs_m.o soil_m.o soilv_m.o soilsnow_m.o stime_m.o tkeeps.o vegpar_m.o vertmix.o vvel_m.o work2_m.o work3_m.o work3f_m.o
scrnout.o : arrays_m.o cc_mpi.o cc_omp.o const_phys.o diag_m.o estab.o extraout_m.o liqwpar_m.o mlo.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o pbl_m.o permsurf_m.o prec_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o work2_m.o
seaesfrad.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o cc_mpi.o cfrac_m.o co2_read.o const_phys.o esfsw_driver.o esfsw_parameters.o estab.o extraout_m.o filnames_m.o infile.o latlong_m.o longwave_params.o microphys_rad.o mlo.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o ozoneread.o pbl_m.o raddiag_m.o radisw_m.o rad_utilities.o sealw99.o sigs_m.o soil_m.o soilsnow_m.o work3f_m.o zenith.o kuocom.h
sealw99.o : gas_tf.o longwave_clouds.o longwave_fluxes.o longwave_params.o longwave_tables.o lw_gases_stdtf.o optical_path.o rad_utilities.o
setxyz.o : cc_mpi.o const_phys.o indices_m.o jimcc.o latlong_m.o map_m.o newmpar_m.o utilities.o workglob_m.o 
sflux.o : arrays_m.o ateb.o cable_ccam2.o cc_mpi.o cc_omp.o const_phys.o dates_m.o diag_m.o estab.o extraout_m.o gdrag_m.o latlong_m.o liqwpar_m.o map_m.o mlo.o mlodynamicsarrays_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o riverarrays_m.o savuvt_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o vecsuv_m.o vegpar_m.o work2_m.o work3_m.o xyzinfo_m.o 
soilsnow.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o morepbl_m.o newmpar_m.o nsibd_m.o parm_m.o permsurf_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o work2_m.o work3_m.o work3b_m.o 
soilv_m.o : newmpar_m.o
spa88.o :  cldcom_m.o kdacom_m.o lwout_m.o newmpar_m.o radisw_m.o rdflux_m.o srccom_m.o tfcom_m.o hcon.h rdparm.h rnddta.h
staguv.o : cc_mpi.o indices_m.o map_m.o newmpar_m.o parm_m.o parmdyn_m.o vecsuv_m.o
swr99.o : newmpar_m.o parm_m.o hcon.h rdparm.h
table.o : newmpar_m.o radisw_m.o tabcom_m.o hcon.h rdparm.h rnddta.h
timeseries.o : arrays_m.o cable_define_types.o carbpools_m.o cc_mpi.o const_phys.o dates_m.o infile.o extraout_m.o morepbl_m.o newmpar_m.o nharrs_m.o parmgeom_m.o pbl_m.o prec_m.o sigs_m.o soil_m.o soilsnow_m.o tracermodule.o tracers_m.o vecsuv_m.o vegpar_m.o vvel_m.o xyzinfo_m.o 
tracermodule.o : arrays_m.o cc_mpi.o const_phys.o dates_m.o infile.o latlong_m.o newmpar_m.o parm_m.o sigs_m.o sumdd_m.o tracers_m.o xyzinfo_m.o 
trvmix.o : arrays_m.o cc_mpi.o cable_ccam2.o cable_define_types.o carbpools_m.o cc_mpi.o const_phys.o dates_m.o diag_m.o newmpar_m.o nsibd_m.o parm_m.o pbl_m.o sigs_m.o tracermodule.o tracers_m.o xyzinfo_m.o 
updps.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o indices_m.o map_m.o newmpar_m.o nlin_m.o parm_m.o parmdyn_m.o parmhor_m.o savuvt_m.o savuv1_m.o sigs_m.o staguv.o vecsuv_m.o vvel_m.o xarrs_m.o xyzinfo_m.o 
upglobal.o : aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o cfrac_m.o const_phys.o diag_m.o epst_m.o indices_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o parmhor_m.o sbar_m.o sigs_m.o staguv.o tkeeps.o tracers_m.o unn_m.o vadvtvd.o vecsuv_m.o vvel_m.o work3f_m.o xarrs_m.o xyzinfo_m.o kuocom.h
usage_m.o: cc_mpi.o
utilities.o : const_phys.o
vadvtvd.o : aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o cfrac_m.o diag_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o sigs_m.o tkeeps.o tracers_m.o vvel_m.o xarrs_m.o kuocom.h
vertmix.o : aerosolldr.o arrays_m.o cc_mpi.o cc_omp.o carbpools_m.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocomb_m.o liqwpar_m.o map_m.o mlo.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o pbl_m.o savuvt_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o tkeeps.o tracermodule.o tracers_m.o trvmix.o work2_m.o kuocom.h
