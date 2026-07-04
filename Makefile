# Configuration
FC      := gfortran
FITSDIR := /opt/local/lib
FLAGS   := -fopenmp -O3 -fdefault-real-8 -w -L$(FITSDIR) -llapack -lblas -lcfitsio -Wl,-rpath,$(FITSDIR)
MODULES := multifield_globals.f90 multifield_grf.f90 multifield_utils.f90

# Programs
PROG1 := ps_mode
PROG2 := ps_full
PROG3 := bkgd_onetraj
PROG4 := potential_cmap

# Default target
all: $(PROG1) $(PROG2) $(PROG3) $(PROG4)

# Build rules
$(PROG1): $(MODULES) $(PROG1).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG2): $(MODULES) $(PROG2).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG3): $(MODULES) $(PROG3).f90
	$(FC) $(FLAGS) $^ -o $@

$(PROG4): $(MODULES) $(PROG4).f90
	$(FC) $(FLAGS) $^ -o $@

# Clean targets
clean:
	@echo "Removing .o and .mod files..."
	@rm -f *.o *.mod

cleanout:
	@echo "Removing executables..."
	@rm -f $(PROG1)

cleanall: clean cleanout

# Phony targets
.PHONY: all clean cleanout cleanall

