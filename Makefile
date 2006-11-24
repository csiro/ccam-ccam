# Makefile for offline CABLE LSM:
# either "make" (equivalent to "make netcdf") or "make text".
# Gab Abramowitz - gabsun@gmail.com
PROG = ./cable

FC = ifort #g95
FFLAGS =  
NCDIR = /usr/local/netcdf-3.6.1_ifort/

netcdf: $(PROG)
	$(PROG)

text: cable_txt # non-netcdf version of offline CABLE 
	./cable_txt	

$(PROG): cable_driver.o 
	$(FC) $(FFLAGS) -o $(PROG) cable_driver.o cable_cbm.o cable_input.o cable_output.o cable_parameters.o cable_checks.o cable_variables.o cable_soilsnow.o cable_carbon.o -L$(NCDIR)lib -lnetcdf

cable_driver.o: cable_driver.f90 cable_output.o cable_parameters.o
	$(FC) $(FFLAGS) -c cable_driver.f90

cable_txt: cable_drivertxt.o 
	$(FC) $(FFLAGS) -o cable_txt cable_drivertxt.o cable_cbm.o cable_checks.o cable_output_text.o cable_parameters.o cable_variables.o cable_soilsnow.o cable_carbon.o

cable_drivertxt.o: cable_drivertxt.f90 cable_output_text.o cable_parameters.o 
	$(FC) $(FFLAGS) -c cable_drivertxt.f90

cable_variables.o: cable_variables.f90
	$(FC) $(FFLAGS) -c cable_variables.f90

cable_soilsnow.o: cable_soilsnow.f90 cable_variables.o
	$(FC) $(FFLAGS) -c cable_soilsnow.f90

cable_carbon.o: cable_carbon.f90 cable_variables.o
	$(FC) $(FFLAGS) -c cable_carbon.f90

cable_parameters.o: cable_parameters.f90 cable_variables.o
	$(FC) $(FFLAGS) -c cable_parameters.f90

cable_cbm.o: cable_cbm.f90 cable_carbon.o cable_soilsnow.o cable_parameters.o
	$(FC) $(FFLAGS) -c cable_cbm.f90

cable_checks.o: cable_checks.f90 cable_cbm.o
	$(FC) $(FFLAGS) -c cable_checks.f90

cable_input.o: cable_input.f90 cable_checks.o
	$(FC) $(FFLAGS) -c cable_input.f90

cable_output.o: cable_output.f90 cable_input.o
	$(FC) $(FFLAGS) -c cable_output.f90

cable_output_text.o: cable_output_text.f90 cable_checks.o
	$(FC) $(FFLAGS) -c cable_output_text.f90
clean:
	rm -f *.o $(PROG) cable_txt *.mod
	ln -s $(NCDIR)src/f90/netcdf.mod ./

