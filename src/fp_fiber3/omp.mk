
# from PETESC_DIR/share/Makefile.user===============================================

# The following variable must either be a path to PETSc.pc or just "PETSc" if PETSc.pc
# has been installed to a system location or can be found in PKG_CONFIG_PATH.

PETSC_DIR=${PETSC_DIR_OMP}

PETSc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/PETSc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.


PACKAGES := $(PETSc.pc)

CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
# CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
# FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
# CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
# FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
INCLUDE_FLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lgomp


# User-----------------------------------------------

CFLAGS += -O3 -std=c99 -funroll-loops

OBJECTS = main.o parameter.o CalculateFunction.o CalculatePressure.o petsc_solver.o sub_function.o CalculateSolid.o
OBJECTS2 = main_i.o parameter.o CalculateFunction.o CalculatePressure.o petsc_solver.o velocity_implicit.o sub_function.o CalculateSolid.o
OBJECTS3 = main_t.o parameter.o CalculateFunction.o CalculatePressure.o petsc_solver.o sub_function.o CalculateSolid.o
OBJECTS4 = init.o parameter.o CalculateFunction.o CalculatePressure.o petsc_solver.o sub_function.o CalculateSolid.o
OBJECT_TEST = test_celllist.o
TARGETS = mps_omp

all : $(TARGETS)

.SUFFIXES:
.SUFFIXES: .c .o
.c.o:
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE_FLAGS) $<

mps_omp : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

mps_i_omp : $(OBJECTS2)
	$(CC) -o $@ $(OBJECTS2) $(LDFLAGS) $(LDLIBS)

mps_t_omp : $(OBJECTS3)
	$(CC) -o $@ $(OBJECTS3) $(LDFLAGS) $(LDLIBS)

init_omp : $(OBJECTS4)
	$(CC) -o $@ $(OBJECTS4) $(LDFLAGS) $(LDLIBS)

mps_gprof_omp : $(OBJECTS)
	$(CC) -pg -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

# Don't forget -pg in .c.o command

clean:
	rm -f core *~ *.o $(TARGETS) init mps_i mps_t
	# rm mps_i
	# rm mps_t
print:
	@echo CC=$(CC)
	@echo CFLAGS=$(CFLAGS)
	@echo INCLUDE_FLAGS=$(INCLUDE_FLAGS)
	@echo LDFLAGS=$(LDFLAGS)
	@echo LDLIBS=$(LDLIBS)


# =====[local]=========
# CC=gcc
# CFLAGS=-Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fno-stack-check -Qunused-arguments -fvisibility=hidden -g3
# INCLUDE_FLAGS=-I/opt/petsc/petsc-3.13.4/arch-darwin-c-debug/include -I/opt/petsc/petsc-3.13.4/include
# LDFLAGS=-L/opt/petsc/petsc-3.13.4/arch-darwin-c-debug/lib -Wl,-rpath,/opt/petsc/petsc-3.13.4/arch-darwin-c-debug/lib
# LDLIBS=-lpetsc -lm
# CFLAGS += -O3

# =====[cluster]=====
# CC=gcc
# CFLAGS=-fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden
# INCLUDE_FLAGS=-I/nfs/maxhome00/keno/lib/petsc-3.13.4/arch-linux2-c-debug/include -I/nfs/maxhome00/keno/lib/petsc-3.13.4/include
# LDFLAGS=-L/nfs/maxhome00/keno/lib/petsc-3.13.4/arch-linux2-c-debug/lib -Wl,-rpath,/nfs/maxhome00/keno/lib/petsc-3.13.4/arch-linux2-c-debug/lib
# LDLIBS=-lpetsc -lm
# CFLAGS += -O3 -std=c99


# -------[make print]--------
# INCLUDE_FLAGS -> (CPPFLAGS)

# @echo CXX=$(CXX)
# @echo FC=$(FC)
# @echo CXXFLAGS=$(CXXFLAGS)
# @echo FFLAGS=$(FFLAGS)

