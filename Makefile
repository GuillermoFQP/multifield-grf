# Configuration
FC      := gfortran
FITSDIR := /opt/local/lib
FLAGS   := -fopenmp -O3 -fdefault-real-8 -w -L$(FITSDIR) -llapack -lblas -lcfitsio -Wl,-rpath,$(FITSDIR)
SRCDIR  := src
MODULES := $(addprefix $(SRCDIR)/,multifield_globals.f90 multifield_grf.f90 multifield_utils.f90)

# Choose programs to build: use filenames from src/ without .f90.
# Add or remove names here; unlisted program sources are ignored by make all.
# You can also select programs for one invocation: make PROGRAMS="ps_mode ps_full"
PROGRAMS := ps_mode ps_full bkgd_onetraj potential_cmap efold_cmap-test

# Each recipe writes the same module files, so builds must run sequentially.
.NOTPARALLEL:
.DEFAULT_GOAL := all

# Default target
all: $(PROGRAMS)

# Shared rule: compile modules in dependency order, then the program.
$(PROGRAMS): %: $(MODULES) $(SRCDIR)/%.f90
	$(FC) $(FLAGS) $^ -o $@

# Clean targets
clean:
	@echo "Removing .o and .mod files..."
	@rm -f *.o *.mod

cleanout:
	@echo "Removing executables..."
	@rm -f $(PROGRAMS)

cleanall: clean cleanout

.PHONY: all clean cleanout cleanall
