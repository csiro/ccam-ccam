FC = mpif90

VPATH = main/general:main/grid:main/parallel:main/file
VPATH += :dynamic/atmosphere:dynamic/ocean:dynamic/river
VPATH += :physics/gwdrag:physics/cloud:physics/convection
VPATH += :physics/radiation:physics/radiation/seaesf:physics/radiation/lhsf
VPATH += :physics/turbmix	
VPATH += :surface:surface/CABLE:surface/uclem
VPATH += :chemistry/aerosol:chemistry/tracers

INC = -I .


# Common compiler flags
ifneq ($(CUSTOM),yes)
NCFLAG = -I $(NETCDF_ROOT)/include
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dusempi3 -Dshare_ifullg
endif
FOPT = -O3
FHOST = -xHost
FOVERRIDE =
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH =
ifeq ($(XEONPHI),yes)
FHOST = -xMIC-AVX512
endif
ifeq ($(BROADWELL),yes)
FHOST = -xCORE-AVX2 -align array32byte
FOVERRIDE = -qoverride-limits
ZMM = -qopt-zmm-usage=high
VTHRESH = -vec-threshold0
IPFLAG =
endif
ifeq ($(SKYLAKE),yes)
FHOST = -xSKYLAKE-AVX512 -align array64byte
FOVERRIDE = -qoverride-limits
ZMM = -qopt-zmm-usage=high
VTHRESH = -vec-threshold0
IPFLAG =
endif
ifeq ($(CASCADELAKE),yes)
FHOST = -xCASCADELAKE -align array64byte
FOVERRIDE = -qoverride-limits
ZMM = -qopt-zmm-usage=high
VTHRESH = -vec-threshold0
IPFLAG =
endif
ifeq ($(SAPPHIRERAPIDS),yes)
FHOST = -xSAPPHIRERAPIDS -align array64byte
FOVERRIDE = -qoverride-limits
ZMM = -qopt-zmm-usage=high
VTHRESH = -vec-threshold0
IPFLAG =
endif
ifeq ($(ZEN3),yes)
FHOST = -axCORE-AVX2 -align array32byte
FOVERRIDE = -qoverride-limits
ZMM = -qopt-zmm-usage=high
VTHRESH = -vec-threshold0
IPFLAG =
endif
ifeq ($(MAGNUS),yes)
FC = ftn
FHOST = -xHost
FOVERRIDE = -qoverride-limits
ZMM = -qopt-zmm-usage=high
VTHRESH = -vec-threshold0
endif
# Default intel compiler options
FFLAGS = $(FHOST) -assume byterecl -ftz -fp-model precise -no-fma -traceback $(MPIFLAG)
ifeq ($(OMP),yes)
FFLAGS += -qopenmp -qno-openmp-simd
endif
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf
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
FOPT = -O3
FHOST = -fallow-argument-mismatch -march=native
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dusempi3 -Dshare_ifullg
endif
FFLAGS = -ftree-vectorize -fstack-arrays -lmvec $(FHOST) -fbacktrace $(MPIFLAG) -Wl,--as-needed -Wl,--disable-new-dtags  -Wl,--rpath -Wl,${LD_RUN_PATH}
FOVERRIDE =
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH =
ifeq ($(GPU),yes)
FFLAGS += -DGPU -foffload=nvptx-none
endif
ifeq ($(OMP),yes)
FFLAGS += -fopenmp
endif
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
PPFLAG90F =
REAL8FLAG = -fdefault-real-8
INT8FLAG = -fdefault-int-8
DEBUGFLAG = -g -Wall -Wextra -fbounds-check
endif

#SETONIX
ifeq ($(SETONIX),yes)
MPIFC = ftn
MPIF77 = ftn
FC = ftn
FOPT = -O3
FHOST = -march=native
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dusempi3
endif
NCFLAG =
FFLAGS = -mtune=native $(FHOST) -fbacktrace $(MPIFLAG) $(NCFLAG) -fallow-argument-mismatch -I /opt/cray/pe/mpich/8.1.27/ofi/gnu/9.1/include
LIB = -lnetcdf
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

#NVFORTRAN
ifeq ($(NVFORTRAN),yes)
MPIFC = nvfortran
MPIF77 = nvfortran
FC = mpifort
NCFLAG = -I $(NETCDF_ROOT)/include/GNU
LIBS = -L $(NETCDF_ROOT)/lib/GNU -lnetcdf
FOPT = -O4
FHOST = -fast -tp=host
ifeq ($(SKYLAKE),yes)
FHOST = -fast -tp=skylake-avx512
endif
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dusempi3 -Dshare_ifullg
endif
FFLAGS = $(FHOST) -traceback $(MPIFLAG) $(NCFLAG)
ifeq ($(GPU),yes)
#FFLAGS += -Minfo=accel -acc -gpu=cc60,cc70,cc80,fastmath,flushz -DGPU
FFLAGS += -Minfo=accel -acc -gpu=cuda12.8,fastmath,flushz -DGPU
endif
ifeq ($(OMP),yes)
FFLAGS += -mp
endif
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
DEBUGFLAG = -g -Mbounds
endif

# CRAY compiler options
ifeq ($(CRAY),yes)
FC = ftn
FFLAGS = -h noomp -h noacc
FOVERRIDE =
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH =
PPFLAG90 = -eZ
PPFLAG77 = -eZ
PPFLAG90F = -eZ
REAL8FLAG = -s real64
INT8FLAG = -s integer64
DEBUGFLAG = -g -R b -K trap=fp
endif

# MAUI compiler options
ifeq ($(MAUI),yes)
MPIFC = ftn
MPIF77 = ftn
FC = ftn
NCFLAG = -I $(NETCDF_ROOT)/include
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dusempi3 -Dshare_ifullg
endif
FOPT = -O3
FHOST = -xSKYLAKE-AVX512
FFLAGS = $(FHOST) -assume byterecl -ftz -fp-model precise -no-fma -traceback $(MPIFLAG) $(NCFLAG)
ifeq ($(OMP),yes)
FFLAGS += -qopenmp -qno-openmp-simd
endif
FOVERRIDE = -qoverride-limits
ZMM =
IPFLAG =
IPOFLAG =
VTHRESH = -vec-threshold0
PPFLAG90 = -fpp
PPFLAG77 = -fpp
PPFLAG90F = -fpp
REAL8FLAG = -r8
INT8FLAG = -i8
DEBUGFLAG = -check all -debug all -fpe0
endif

# llvm compiler options
ifeq ($(FLANG),yes)
MPIFC = mpifort
MPIF77 = mpif77
FC = mpifort
FOPT = -O3
FHOST = -march=native -fopenmp
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dusempi3 -Dshare_ifullg
endif
FFLAGS = -mtune=native $(FHOST) -fbacktrace $(MPIFLAG) $(NCFLAG)
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
REAL8FLAG = -r8
INT8FLAG = -i8
endif


# IBM compiler options
#ifeq ($(IBM),yes)
#FC = xlf
#MPIFLAG = -Dusempi3 -Dshare_ifullg
#FFLAGS = -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -qextname $(MPIFLAG)
#LIBS = -L $(NETCDF_ROOT)/lib -L /opt/ibmhpc/pecurrent/mpich/xlf/lib64 -lnetcdf -lnetcdff -lmpi -lmpigf
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
endif

# Testing - I/O and fpmodel
ifeq ($(TEST),yes)
FFLAGS += $(DEBUGFLAG) -Ddebug
endif

# Build with 64 ints/reals
ifeq ($(I8R8),yes)
FFLAGS += $(REAL8FLAG) $(INT8FLAG) -Di8r8
endif

# Use Netcdf3
ifeq ($(NETCDF3),yes)
FFLAGS += -Dusenc3
endif

#MPI 
ifeq ($(MPIMOD),yes)
FFLAGS += -Dusempimod
# The following fails on Cray due to a possible bug in the f08 module
#FFLAGS += -Dusempimod_f08 
endif

# COSP
ifeq ($(COSP),yes)
VPATH = ../../program/COSPv2.0/build
FFLAGS += -DCOSP -I $(VPATH)
LIBS += -lcosp -lsubcol -L $(VPATH)
endif

# CABLE, MLO, aTEB, etc
FFLAGS += -DCCAM $(INC)


# Object files for dynamical model
OBJS = adjust5.o amipsst.o convjlm.o convjlm22.o depts.o estab.o gettin.o \
globpe.o gdrag_m.o hordifg.o hs_phys.o indata.o infile.o ints.o \
helmsolve.o jimcc.o nesting.o nonlin.o ensemble.o \
outcdf.o radriv90.o scrnout.o setxyz.o sflux.o \
soilsnow.o staguv.o upglobal.o eig.o updps.o vadvtvd.o \
leoncld.o cloudmod.o latltoij.o module_aux_rad.o module_ctrl_microphysics.o \
module_mp_sbu_ylin.o \
cldblk.o clddia.o clo89.o cloud.o cloud2.o co2_read.o e1e288.o \
e3v88.o fst88.o hconst.o lwr88.o ozoneread.o spa88.o \
swr99.o table.o zenith.o cc_acc.o cc_mpi.o diag_m.o sumdd_m.o daviesnudge.o \
utilities.o onthefly.o tracermodule.o timeseries.o \
trvmix.o getopt_m.o usage_m.o const_phys.o \
betts.o bett_cuc.o bettinit.o bettrain.o bettspli.o \
xyzinfo_m.o vecsuv_m.o map_m.o latlong_m.o indices_m.o bigxy4_m.o \
arrays_m.o betts1_m.o carbpools_m.o cldcom_m.o co2dta_m.o cfrac_m.o \
dpsdt_m.o epst_m.o extraout_m.o histave_m.o kdacom_m.o \
kuocom_m.o liqwpar_m.o lwout_m.o morepbl_m.o nharrs_m.o \
nlin_m.o nsibd_m.o parmhdff_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o \
radisw_m.o rdflux_m.o riverarrays_m.o savuvt_m.o savuv1_m.o sbar_m.o screen_m.o \
sigs_m.o soil_m.o soilsnow_m.o srccom_m.o swocom_m.o tabcom_m.o \
tbar2d_m.o tfcom_m.o tracers_m.o unn_m.o uvbar_m.o vecs_m.o vegpar_m.o vvel_m.o \
workglob_m.o work2_m.o work3_m.o work3b_m.o work3f_m.o work3lwr_m.o work3sav_m.o \
xarrs_m.o \
aerointerface.o aerosolldr.o aerosol_arrays.o \
cable_air.o cable_albedo.o cable_canopy.o cable_ccam2.o cable_ccam3.o cable_ccam4.o cable_common.o \
cable_data.o cable_define_types.o cable_radiation.o cable_roughness.o cable_soilsnow.o \
cable_optimiseJVratio.o cable_psm.o cable_pft_params.o cable_soil_params.o cable_constants.o \
cable_gw_hydro.o \
cable_sli_main.o cable_sli_numbers.o cable_sli_roots.o cable_sli_solve.o cable_sli_utils.o \
casa_cnp.o casa_variable.o POP.o casa_dimension.o casa_phenology.o casa_param.o \
uclem.o uclem_ctrl.o uclem_parameters.o uclem_types.o mlo_ctrl.o mlo.o river.o \
module_ctrl_turbmix.o trimmix.o vertmix.o tkeeps.o \
seaesfrad.o rad_utilities.o microphys_rad.o esfsw_driver.o esfsw_parameters.o \
longwave_params.o sealw99.o longwave_clouds.o longwave_fluxes.o longwave_tables.o \
optical_path.o gas_tf.o lw_gases_stdtf.o \
mlodynamics.o mlodynamicsarrays_m.o mlodiffg.o mlostag.o mlodepts.o mloints.o mlovadvtvd.o \
darcdf_m.o dates_m.o filnames_m.o newmpar_m.o parm_m.o parmdyn_m.o parmgeom_m.o \
parmhor_m.o soilv_m.o stime_m.o \
netcdf_m.o parmvert_m.o module_aux_cosp.o module_ctrl_convection.o \
cu_gf_deep.o


globpea: $(OBJS)
	$(FC) -o globpea $(FOPT) $(FFLAGS) $(OBJS) $(LIBS)

clean:
	rm *.o *.i *.mod

.SUFFIXES:.f90 .F90


netcdf_m.o: netcdf_m.f90
	$(FC) -c $(PPFLAG90) $(NCFLAG) $<
esfsw_driver.o: esfsw_driver.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FOPT) $(FFLAGS) $(VTHRESH) $(IPFLAG) $<
esfsw_parameters.o gas_tf.o longwave_clouds.o longwave_fluxes.o longwave_tables.o longwave_params.o lw_gases_stdtf.o microphys_rad.o optical_path.o rad_utilities.o sealw99.o: %.o: %.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FOPT) $(FFLAGS) $(IPFLAG) $<
cable_air.o cable_albedo.o cable_canopy.o cable_common.o cable_constants.o cable_data.o cable_define_types.o cable_gw_hydro.o cable_optimiseJVratio.o cable_pft_params.o cable_psm.o cable_radiation.o cable_roughness.o cable_sli_main.o cable_sli_numbers.o cable_sli_roots.o cable_sli_solve.o cable_sli_utils.o cable_soil_params.o cable_soilsnow.o casa_cnp.o casa_dimension.o casa_param.o casa_phenology.o casa_variable.o POP.o: %.o: %.F90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FOPT) $(FFLAGS) $(IPFLAG) $<
module_mp_sbu_ylin.o: module_mp_sbu_ylin.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FOPT) $(FFLAGS) $(IPFLAG) $<
estab.o: estab.f90
	$(FC) -c $(FOPT) $(FFLAGS) $(IPOFLAG) $(PPFLAG90) $(IPFLAG) $<
helmsolve.o: helmsolve.f90
	$(FC) -c $(PPFLAG90) $(FOPT) $(FFLAGS) $(FOVERRIDE) $(IPFLAG) $<
ints.o: ints.f90
	$(FC) -c $(FOPT) $(FFLAGS) $(ZMM) $(PPFLAG90) $(IPFLAG) $<
leoncld.o: leoncld.f90
	$(FC) -c $(FOPT) $(FFLAGS) $(IPOFLAG) $(PPFLAG90) $(IPFLAG) $<
seaesfrad.o: seaesfrad.f90
	$(FC) -c $(FOPT) $(FFLAGS) $(VTHRESH) $(PPFLAG90) $(IPFLAG) $<
tkeeps.o: tkeeps.f90
	$(FC) -c $(FOPT) $(FFLAGS) $(PPFLAG90) $(IPFLAG) $<
vertmix.o: vertmix.f90
	$(FC) -c $(FOPT) $(FFLAGS) $(IPOFLAG) $(PPFLAG90) $(IPFLAG) $<
version.h: FORCE
	rm -f tmpver
	echo "character(len=*), parameter :: version= &" > tmpver
	echo "'CCAM `git log | head -3 | tail -1`" "`git log | head -1`' " >> tmpver
	cmp tmpver version.h || mv tmpver version.h
FORCE:


.f90.o:
	$(FC) -c $(FOPT) $(FFLAGS) $(PPFLAG90) $(IPFLAG) $<
.F90.o:
	$(FC) -c $(FOPT) $(FFLAGS) $(PPFLAG90F) $(IPFLAG) $<
.f.o:
	$(FC) -c $(FOPT) $(FFLAGS) $(PPFLAG77) $(IPFLAG) $< 

# Remove mod rule from Modula 2 so GNU make doesn't get confused
%.o : %.mod

# Dependencies
adjust5.o : aerosol_arrays.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o dpsdt_m.o epst_m.o helmsolve.o indices_m.o kuocom_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o pbl_m.o sigs_m.o staguv.o tbar2d_m.o tracers_m.o vecsuv_m.o vecs_m.o vvel_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o
aerointerface.o : aerosol_arrays.o aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o extraout_m.o infile.o kuocom_m.o latlong_m.o liqwpar_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o ozoneread.o parm_m.o parmgeom_m.o pbl_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o trimmix.o vegpar_m.o work2_m.o zenith.o
aerosolldr.o : aerosol_arrays.o newmpar_m.o
amipsst.o : arrays_m.o cc_mpi.o const_phys.o dates_m.o filnames_m.o infile.o latlong_m.o latltoij.o mlo_ctrl.o nesting.o newmpar_m.o nharrs_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o setxyz.o soil_m.o soilsnow_m.o workglob_m.o
bett_cuc.o : betts1_m.o newmpar_m.o
bettinit.o : betts1_m.o newmpar_m.o
bettrain.o : betts1_m.o newmpar_m.o
betts.o : betts1_m.o morepbl_m.o newmpar_m.o parm_m.o prec_m.o sigs_m.o
cable_air.o : cable_common.o cable_data.o cable_define_types.o 
cable_albedo.o : cable_common.o cable_data.o cable_define_types.o
cable_canopy.o : cable_air.o cable_common.o cable_data.o cable_define_types.o cable_gw_hydro.o cable_psm.o cable_radiation.o cable_roughness.o cable_sli_main.o cable_sli_utils.o
cable_common.o : cable_define_types.o cable_pft_params.o cable_soil_params.o
cable_ccam2.o : arrays_m.o cable_air.o cable_albedo.o cable_canopy.o cable_ccam3.o cable_ccam4.o cable_common.o cable_define_types.o cable_optimiseJVratio.o cable_radiation.o cable_roughness.o cable_sli_main.o cable_soilsnow.o carbpools_m.o casa_cnp.o casa_variable.o casa_phenology.o cc_mpi.o const_phys.o darcdf_m.o dates_m.o estab.o extraout_m.o infile.o latlong_m.o liqwpar_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o POP.o prec_m.o raddiag_m.o radisw_m.o riverarrays_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o tracers_m.o vegpar_m.o work2_m.o work3_m.o zenith.o
cable_ccam3.o : cable_data.o cable_define_types.o casa_phenology.o casa_variable.o POP.o
cable_ccam4.o : cable_ccam3.o cable_data.o cable_define_types.o casa_variable.o newmpar_m.o POP.o
cable_constants.o : cable_define_types.o
cable_data.o : cable_constants.o
cable_gw_hydro.o : cable_define_types.o cable_common.o cable_soilsnow.o cable_data.o
cable_optimiseJVratio.o : cable_canopy.o cable_data.o cable_define_types.o POP.o
cable_pft_params.o : cable_define_types.o
cable_psm.o : cable_common.o cable_define_types.o
cable_radiation.o : cable_common.o cable_data.o cable_define_types.o
cable_roughness.o : cable_common.o cable_data.o cable_define_types.o
cable_sli_main.o : cable_common.o cable_define_types.o cable_sli_numbers.o cable_sli_roots.o cable_sli_solve.o cable_sli_utils.o
cable_sli_numbers.o : cable_define_types.o
cable_sli_roots.o : cable_define_types.o cable_sli_numbers.o
cable_sli_solve.o : cable_define_types.o cable_sli_numbers.o cable_sli_utils.o
cable_sli_utils.o : cable_define_types.o cable_sli_numbers.o
cable_soil_params.o : cable_define_types.o
cable_soilsnow.o : cable_common.o cable_data.o cable_define_types.o parm_m.o
carbpools_m.o : cable_define_types.o casa_variable.o parm_m.o
casa_cnp.o : cable_common.o cable_define_types.o casa_variable.o
casa_dimension.o : cable_define_types.o
casa_param.o : casa_dimension.o
casa_phenology.o : casa_dimension.o
casa_variable.o : cable_define_types.o casa_dimension.o casa_param.o
cc_mpi.o : const_phys.o indices_m.o latlong_m.o map_m.o newmpar_m.o sumdd_m.o vecsuv_m.o workglob_m.o xyzinfo_m.o
clddia.o : arrays_m.o cc_mpi.o const_phys.o kuocom_m.o map_m.o morepbl_m.o newmpar_m.o parm_m.o pbl_m.o sigs_m.o soil_m.o vvel_m.o
clo89.o : cldcom_m.o parm_m.o radisw_m.o rdparm.h
cloud2.o : diag_m.o cc_mpi.o const_phys.o kuocom_m.o leoncld.o newmpar_m.o parm_m.o radisw_m.o sigs_m.o hcon.h rdparm.h
cloud.o : extraout_m.o newmpar_m.o parm_m.o radisw_m.o rdparm.h
cloudmod.o : const_phys.o estab.o parm_m.o sigs_m.o
co2_read.o : cc_mpi.o co2dta_m.o dates_m.o filnames_m.o infile.o newmpar_m.o parm_m.o radisw_m.o stime_m.o rdparm.h
convjlm.o : aerointerface.o aerosol_arrays.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocom_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o prec_m.o sigs_m.o soil_m.o tracers_m.o vvel_m.o work2_m.o
convjlm22.o : aerointerface.o aerosol_arrays.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocom_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o prec_m.o sigs_m.o soil_m.o tracers_m.o vvel_m.o work2_m.o
daviesnudge.o : aerosol_arrays.o arrays_m.o cc_mpi.o newmpar_m.o parm_m.o sigs_m.o
depts.o : bigxy4_m.o cc_acc.o cc_mpi.o const_phys.o indices_m.o map_m.o newmpar_m.o parm_m.o parmhor_m.o parmgeom_m.o uvbar_m.o vecsuv_m.o work3f_m.o xyzinfo_m.o 
diag_m.o : cc_mpi.o newmpar_m.o parm_m.o sigs_m.o sumdd_m.o xyzinfo_m.o
e1e288.o : kdacom_m.o radisw_m.o tabcom_m.o tfcom_m.o hcon.h rdparm.h
e3v88.o : tabcom_m.o hcon.h rdparm.h
eig.o : cc_mpi.o const_phys.o newmpar_m.o parm_m.o vecs_m.o
ensemble.o : aerosol_arrays.o arrays_m.o cc_mpi.o dates_m.o filnames_m.o mlo_ctrl.o newmpar_m.o onthefly.o outcdf.o parm_m.o pbl_m.o soil_m.o soilsnow_m.o xyzinfo_m.o
esfsw_driver.o : esfsw_parameters.o rad_utilities.o
esfsw_parameters.o : rad_utilities.o
estab.o : const_phys.o
fst88.o : cc_mpi.o cldcom_m.o diag_m.o kdacom_m.o lwout_m.o newmpar_m.o parm_m.o radisw_m.o rdflux_m.o srccom_m.o tabcom_m.o tfcom_m.o hcon.h rdparm.h rnddta.h
gas_tf.o : longwave_params.o rad_utilities.o
gdrag_m.o : cc_mpi.o arrays_m.o const_phys.o newmpar_m.o nharrs_m.o parm_m.o pbl_m.o sigs_m.o
gettin.o : arrays_m.o newmpar_m.o savuvt_m.o
globpe.o : aerointerface.o aerosol_arrays.o amipsst.o arrays_m.o bigxy4_m.o cc_acc.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o darcdf_m.o dates_m.o daviesnudge.o diag_m.o dpsdt_m.o ensemble.o epst_m.o estab.o extraout_m.o filnames_m.o gdrag_m.o getopt_m.o histave_m.o hordifg.o hs_phys.o indata.o indices_m.o infile.o kuocom_m.o latlong_m.o liqwpar_m.o map_m.o mlodynamics.o module_aux_rad.o module_ctrl_convection.o module_ctrl_microphysics.o module_ctrl_turbmix.o morepbl_m.o nesting.o newmpar_m.o nharrs_m.o nlin_m.o nsibd_m.o outcdf.o ozoneread.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o parmvert_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o river.o savuvt_m.o savuv1_m.o sbar_m.o screen_m.o seaesfrad.o setxyz.o sflux.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o staguv.o stime_m.o tbar2d_m.o timeseries.o tkeeps.o tracermodule.o tracers_m.o trvmix.o unn_m.o usage_m.o uvbar_m.o vecs_m.o vecsuv_m.o vegpar_m.o vvel_m.o workglob_m.o work2_m.o work3_m.o work3f_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o version.h
hconst.o : hcon.h
helmsolve.o : cc_mpi.o diag_m.o indices_m.o newmpar_m.o parm_m.o parmdyn_m.o parmgeom_m.o sumdd_m.o vecs_m.o
histave_m.o: parm_m.o
hordifg.o : aerosol_arrays.o arrays_m.o cc_acc.o cc_mpi.o cfrac_m.o const_phys.o dpsdt_m.o indices_m.o kuocom_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o parmhdff_m.o savuvt_m.o sigs_m.o tkeeps.o vecsuv_m.o vvel_m.o
hs_phys.o : arrays_m.o latlong_m.o newmpar_m.o nlin_m.o parm_m.o sigs_m.o
indata.o : aerosol_arrays.o amipsst.o arrays_m.o bigxy4_m.o cc_mpi.o const_phys.o darcdf_m.o dates_m.o daviesnudge.o diag_m.o eig.o ensemble.o epst_m.o extraout_m.o filnames_m.o gdrag_m.o indices_m.o infile.o kuocom_m.o latlong_m.o latltoij.o liqwpar_m.o map_m.o mlodynamics.o module_ctrl_convection.o module_ctrl_turbmix.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o onthefly.o parm_m.o parmdyn_m.o parmgeom_m.o pbl_m.o permsurf_m.o river.o seaesfrad.o sflux.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o stime_m.o timeseries.o tracermodule.o tracers_m.o vecs_m.o vecsuv_m.o vegpar_m.o xyzinfo_m.o
indices_m.o : newmpar_m.o
infile.o : cc_mpi.o dates_m.o netcdf_m.o newmpar_m.o parm_m.o parmgeom_m.o sigs_m.o
ints.o : cc_acc.o cc_mpi.o indices_m.o newmpar_m.o parm_m.o parmhor_m.o
latltoij.o : const_phys.o newmpar_m.o parm_m.o parmdyn_m.o utilities.o 
leoncld.o : const_phys.o estab.o parm_m.o prec_m.o sigs_m.o
longwave_clouds.o : rad_utilities.o
longwave_fluxes.o : rad_utilities.o
longwave_tables.o : cc_mpi.o filnames_m.o longwave_params.o rad_utilities.o
lw_gases_stdtf.o : cc_mpi.o filnames_m.o infile.o gas_tf.o newmpar_m.o rad_utilities.o
lwr88.o : co2dta_m.o kdacom_m.o parm_m.o radisw_m.o tfcom_m.o work3lwr_m.o hcon.h rdparm.h rnddta.h
microphys_rad.o : esfsw_parameters.o longwave_params.o rad_utilities.o
mlo.o : newmpar_m.o
mlo_ctrl.o : mlo.o newmpar_m.o
mlodepts.o : bigxy4_m.o cc_acc.o cc_mpi.o const_phys.o indices_m.o mlo_ctrl.o newmpar_m.o parm_m.o parmgeom_m.o parmhor_m.o vecsuv_m.o xyzinfo_m.o
mlodiffg.o : cc_acc.o cc_mpi.o const_phys.o indices_m.o map_m.o mlo_ctrl.o mlodynamicsarrays_m.o newmpar_m.o nharrs_m.o parm_m.o soil_m.o vecsuv_m.o
mlodynamics.o : arrays_m.o bigxy4_m.o cc_acc.o cc_mpi.o const_phys.o helmsolve.o indices_m.o infile.o latlong_m.o map_m.o mlo_ctrl.o mlodepts.o mlodiffg.o mlodynamicsarrays_m.o mloints.o mlostag.o mlovadvtvd.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o parmgeom_m.o parmhor_m.o soil_m.o soilsnow_m.o soilv_m.o sumdd_m.o vecsuv_m.o xyzinfo_m.o 
mloints.o : cc_acc.o cc_mpi.o indices_m.o mlo_ctrl.o newmpar_m.o parm_m.o parmhor_m.o
mlostag.o : cc_mpi.o indices_m.o mlo_ctrl.o mlodynamicsarrays_m.o newmpar_m.o parm_m.o
mlovadvtvd.o : cc_acc.o cc_mpi.o mlo_ctrl.o newmpar_m.o
module_aux_cosp.o : arrays_m.o cc_mpi.o cfrac_m.o const_phys.o estab.o kuocom_m.o latlong_m.o liqwpar_m.o map_m.o module_aux_rad.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o pbl_m.o prec_m.o raddiag_m.o sigs_m.o soil_m.o soilsnow_m.o work3f_m.o vvel_m.o
module_aux_rad.o : const_phys.o parm_m.o
module_ctrl_convection.o : aerointerface.o aerosol_arrays.o arrays_m.o cc_mpi.o convjlm.o convjlm22.o const_phys.o cu_gf_deep.o kuocom_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nlin_m.o parm_m.o prec_m.o sigs_m.o soil_m.o vvel_m.o
module_ctrl_microphysics.o : aerointerface.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o estab.o filnames_m.o kuocom_m.o latlong_m.o leoncld.o liqwpar_m.o map_m.o module_aux_cosp.o module_aux_rad.o module_mp_sbu_ylin.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o pbl_m.o prec_m.o raddiag_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o work3f_m.o vvel_m.o
module_ctrl_turbmix.o : arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o extraout_m.o kuocom_m.o liqwpar_m.o map_m.o mlo_ctrl.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o pbl_m.o savuvt_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o tkeeps.o trimmix.o work2_m.o vertmix.o
module_mp_sbu_ylin.o : cc_mpi.o
nesting.o : aerointerface.o aerosol_arrays.o arrays_m.o cc_acc.o cc_mpi.o const_phys.o dates_m.o daviesnudge.o darcdf_m.o diag_m.o indices_m.o kuocom_m.o latlong_m.o liqwpar_m.o map_m.o mlo_ctrl.o mlodynamics.o newmpar_m.o nharrs_m.o onthefly.o parm_m.o parmdyn_m.o parmgeom_m.o pbl_m.o savuvt_m.o savuv1_m.o sigs_m.o soil_m.o soilsnow_m.o stime_m.o vecsuv_m.o work3sav_m.o xyzinfo_m.o
nonlin.o : aerosol_arrays.o arrays_m.o cc_mpi.o const_phys.o diag_m.o epst_m.o hs_phys.o indices_m.o kuocom_m.o latlong_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o savuvt_m.o sigs_m.o staguv.o tbar2d_m.o tkeeps.o tracers_m.o unn_m.o vadvtvd.o vecsuv_m.o vvel_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o
onthefly.o : aerointerface.o cc_acc.o cc_mpi.o const_phys.o darcdf_m.o diag_m.o extraout_m.o histave_m.o infile.o kuocom_m.o latlong_m.o latltoij.o mlodynamics.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmdyn_m.o parmgeom_m.o prec_m.o raddiag_m.o riverarrays_m.o savuvt_m.o savuv1_m.o screen_m.o setxyz.o sflux.o sigs_m.o soil_m.o soilv_m.o stime_m.o tkeeps.o tracers_m.o utilities.o vecsuv_m.o vvel_m.o workglob_m.o work2_m.o xarrs_m.o
optical_path.o : cc_mpi.o filnames_m.o longwave_params.o lw_gases_stdtf.o rad_utilities.o
outcdf.o : aerointerface.o aerosol_arrays.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o dates_m.o daviesnudge.o dpsdt_m.o extraout_m.o filnames_m.o gdrag_m.o histave_m.o infile.o kuocom_m.o latlong_m.o leoncld.o liqwpar_m.o map_m.o mlodynamics.o module_aux_rad.o module_ctrl_microphysics.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o ozoneread.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o parmvert_m.o pbl_m.o prec_m.o raddiag_m.o river.o riverarrays_m.o savuvt_m.o savuv1_m.o screen_m.o sflux.o seaesfrad.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o staguv.o tkeeps.o tracermodule.o tracers_m.o vegpar_m.o vvel_m.o work2_m.o xarrs_m.o version.h
ozoneread.o : cc_mpi.o const_phys.o dates_m.o filnames_m.o infile.o latlong_m.o newmpar_m.o parm_m.o 
radriv90.o : aerointerface.o arrays_m.o cc_mpi.o cfrac_m.o cldcom_m.o co2_read.o co2dta_m.o const_phys.o diag_m.o estab.o extraout_m.o infile.o kdacom_m.o kuocom_m.o latlong_m.o liqwpar_m.o lwout_m.o mlo_ctrl.o newmpar_m.o nharrs_m.o nsibd_m.o ozoneread.o parm_m.o pbl_m.o raddiag_m.o radisw_m.o rdflux_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o srccom_m.o swocom_m.o swr99.o tabcom_m.o tfcom_m.o uclem_ctrl.o work3f_m.o work3lwr_m.o zenith.o rdparm.h hcon.h
river.o : arrays_m.o cc_mpi.o const_phys.o indices_m.o map_m.o newmpar_m.o nsibd_m.o parm_m.o riverarrays_m.o sflux.o soil_m.o soilsnow_m.o soilv_m.o xyzinfo_m.o 
riverarrays_m.o : newmpar_m.o
scrnout.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o estab.o extraout_m.o liqwpar_m.o mlo_ctrl.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o pbl_m.o permsurf_m.o prec_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o tkeeps.o work2_m.o
seaesfrad.o : aerointerface.o aerosol_arrays.o arrays_m.o cc_mpi.o cfrac_m.o co2_read.o const_phys.o esfsw_driver.o esfsw_parameters.o estab.o extraout_m.o filnames_m.o infile.o kuocom_m.o latlong_m.o longwave_params.o microphys_rad.o mlo_ctrl.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o ozoneread.o pbl_m.o raddiag_m.o radisw_m.o rad_utilities.o sealw99.o sigs_m.o soil_m.o soilsnow_m.o uclem_ctrl.o work3f_m.o zenith.o
sealw99.o : cc_mpi.o filnames_m.o gas_tf.o longwave_clouds.o longwave_fluxes.o longwave_params.o longwave_tables.o lw_gases_stdtf.o optical_path.o rad_utilities.o
setxyz.o : cc_mpi.o const_phys.o indices_m.o jimcc.o latlong_m.o map_m.o newmpar_m.o utilities.o workglob_m.o 
sflux.o : arrays_m.o cable_ccam2.o cc_mpi.o const_phys.o dates_m.o diag_m.o estab.o extraout_m.o gdrag_m.o latlong_m.o liqwpar_m.o map_m.o mlo_ctrl.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o riverarrays_m.o savuvt_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o uclem_ctrl.o vecsuv_m.o vegpar_m.o work2_m.o work3_m.o xyzinfo_m.o 
soilsnow.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o morepbl_m.o newmpar_m.o nsibd_m.o parm_m.o permsurf_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o work2_m.o work3_m.o work3b_m.o 
soilv_m.o : newmpar_m.o
spa88.o :  cldcom_m.o kdacom_m.o lwout_m.o radisw_m.o rdflux_m.o srccom_m.o tfcom_m.o hcon.h rdparm.h rnddta.h
staguv.o : cc_mpi.o indices_m.o map_m.o newmpar_m.o parm_m.o parmdyn_m.o vecsuv_m.o
swr99.o : parm_m.o hcon.h rdparm.h
table.o : radisw_m.o tabcom_m.o hcon.h rdparm.h rnddta.h
timeseries.o : arrays_m.o cable_define_types.o carbpools_m.o cc_mpi.o const_phys.o dates_m.o infile.o extraout_m.o morepbl_m.o newmpar_m.o nharrs_m.o parmgeom_m.o pbl_m.o prec_m.o sigs_m.o soil_m.o soilsnow_m.o tracermodule.o tracers_m.o vecsuv_m.o vegpar_m.o vvel_m.o xyzinfo_m.o 
tkeeps.o : mlo_ctrl.o
tracermodule.o : arrays_m.o cc_mpi.o const_phys.o dates_m.o infile.o latlong_m.o newmpar_m.o parm_m.o sigs_m.o sumdd_m.o tracers_m.o xyzinfo_m.o 
trvmix.o : arrays_m.o cable_ccam2.o cable_define_types.o carbpools_m.o const_phys.o dates_m.o diag_m.o morepbl_m.o newmpar_m.o nsibd_m.o parm_m.o pbl_m.o sigs_m.o tracermodule.o tracers_m.o xyzinfo_m.o 
uclem.o : uclem_parameters.o uclem_types.o
uclem_ctrl.o : newmpar_m.o uclem.o uclem_parameters.o uclem_types.o
updps.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o indices_m.o map_m.o newmpar_m.o nlin_m.o parm_m.o parmdyn_m.o parmhor_m.o savuvt_m.o savuv1_m.o sigs_m.o staguv.o vecsuv_m.o vvel_m.o xarrs_m.o xyzinfo_m.o 
upglobal.o : aerointerface.o aerosol_arrays.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o epst_m.o indices_m.o kuocom_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o parmhor_m.o sbar_m.o sigs_m.o staguv.o tkeeps.o tracers_m.o unn_m.o vadvtvd.o vecsuv_m.o vvel_m.o work3f_m.o xarrs_m.o xyzinfo_m.o
usage_m.o: cc_mpi.o
utilities.o : const_phys.o
vadvtvd.o : aerosol_arrays.o arrays_m.o cc_acc.o cc_mpi.o cfrac_m.o diag_m.o kuocom_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o parmvert_m.o sigs_m.o tkeeps.o tracers_m.o vvel_m.o xarrs_m.o
vertmix.o : arrays_m.o cc_mpi.o carbpools_m.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocom_m.o liqwpar_m.o map_m.o mlo_ctrl.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o pbl_m.o savuvt_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o trimmix.o work2_m.o
zenith.o : parm_m.o
