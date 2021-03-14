
# from PETESC_DIR/share/Makefile.user===============================================

# The following variable must either be a path to PETSc.pc or just "PETSc" if PETSc.pc
# has been installed to a system location or can be found in PKG_CONFIG_PATH.
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
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm


# User-----------------------------------------------

CFLAGS += -O3 -std=c99 -funroll-loops

OBJECTS = main.o CalculateFunction.o
OBJECTS2 = read_data.o
TARGETS = 2d

all : $(TARGETS)

.SUFFIXES:
.SUFFIXES: .c .o
.c.o:
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE_FLAGS) $<

2d : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

read : $(OBJECTS2)
	$(CC) -o $@ $(OBJECTS2) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f core *~ *.o $(TARGETS)
	rm read
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

