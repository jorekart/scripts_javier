1) Calculate the initial heat fluxes with field line tracing
	- set i_start, i_end and the jump in job_fluid_loads_3D
 
2) create a new folder and copy the vtk files there

3) Convert the vtk files to .dat files for faster data processing
ifort -qopenmp -o  get_data_files get_data_files.f90 && ./get_data_files

4) Calculate the temperature with those files or do something else

