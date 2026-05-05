# Intel Fortran compiler
FC = ifort

# Base flags
RELEASEFLAGS = -O2 -fpp
DEBUGFLAGS   = -O0 -g -check all -warn all -traceback -fpp

# OpenMP flag
OMPFLAG = -qopenmp

# Choose build type (default = release)
BUILD ?= release
ifeq ($(BUILD),debug)
  FFLAGS = $(DEBUGFLAGS)
else
  FFLAGS = $(RELEASEFLAGS)
endif

ifeq ($(WALL),W_alloy)
  DEFINED = -DW_alloy
endif

ifeq ($(WALL),Be)
  DEFINED = -DBE_WALL
endif

# Module source
MODSRC = mod_utils.f90
OBJMOD = mod_utils.o

# Executables and their sources
EXE1 = get_data_files
EXE2 = analyze_wall_T
EXE3 = calculate_T

OBJ1 = $(OBJMOD) get_data_files.o
OBJ2 = $(OBJMOD) analyze_wall_T.o
OBJ3 = $(OBJMOD) calculate_T.o

# Default target
all: $(EXE1) $(EXE2) $(EXE3)

# Dependencies: programs need the module
get_data_files.o: $(OBJMOD)
analyze_wall_T.o: $(OBJMOD)
calculate_T.o: $(OBJMOD)

# Module rule
$(OBJMOD): $(MODSRC)
	$(FC) $(FFLAGS) $(DEFINED) -c $< -o $@

# Generic rule for objects
%.o: %.f90
	$(FC) $(FFLAGS) $(DEFINED) -c $< -o $@

# Link executables
$(EXE1): $(OBJ1)
	$(FC) $(FFLAGS) -o $@ $(OBJ1)

$(EXE2): $(OBJ2)
	$(FC) $(FFLAGS) $(DEFINED) -o $@ $(OBJ2)

$(EXE3): $(OBJ3)
	$(FC) $(FFLAGS) $(DEFINED) $(OMPFLAG) -o $@ $(OBJ3)

# Clean
clean:
	rm -f *.o *.mod $(EXE1) $(EXE2) $(EXE3)
