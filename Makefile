FC      = gfortran
FFLAGS  = -O3 -march=native -funroll-loops -ffast-math -flto
SRCDIR  = src
TARGET  = ths_mc

# compile order: respects module dependencies
SRCS = \
	$(SRCDIR)/mod_params.f90 \
	$(SRCDIR)/mod_rng.f90 \
	$(SRCDIR)/mod_system.f90 \
	$(SRCDIR)/mod_overlap.f90 \
	$(SRCDIR)/mod_setup.f90 \
	$(SRCDIR)/mod_input.f90 \
	$(SRCDIR)/mod_pivot.f90 \
	$(SRCDIR)/mod_crank.f90 \
	$(SRCDIR)/mod_observables.f90 \
	$(SRCDIR)/mod_io.f90 \
	$(SRCDIR)/main.f90

OBJS = $(patsubst $(SRCDIR)/%.f90, %.o, $(SRCS))

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

debug: FFLAGS = -O0 -g -fcheck=all -fbacktrace -Wall
debug: $(TARGET)

clean:
	rm -f *.o *.mod $(TARGET)

.PHONY: clean debug
