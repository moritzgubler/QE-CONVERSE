# Makefile for QE-CONVERSE unit tests (test-drive framework)
#
# Usage:
#   make          – build test_runner
#   make run      – build and run all unit tests
#   make clean    – remove build artefacts
#   make distclean – alias for clean
#
# testdrive.F90 is vendored in this directory (Apache-2.0 OR MIT).
# No internet access is required.
# Requirements: the main QE-CONVERSE build must have been configured
# (./configure --with-qe-source=...) so that src/Makefile exists and
# contains QE_SOURCE.

# ---------------------------------------------------------------------------
# Derive QE_SOURCE from the parent src/Makefile
# ---------------------------------------------------------------------------
SRC_MAKEFILE := ../../src/Makefile
QE_SOURCE    := $(shell awk -F= '/^QE_SOURCE/{gsub(/ /,"",$$2); print $$2}' \
                    $(SRC_MAKEFILE))

ifeq ($(QE_SOURCE),)
$(error Cannot determine QE_SOURCE from $(SRC_MAKEFILE). \
        Run ./configure --with-qe-source=<path> first.)
endif

include $(QE_SOURCE)/make.inc

# ---------------------------------------------------------------------------
# Module search paths (same as src/Makefile)
# ---------------------------------------------------------------------------
MODFLAGS = \
    $(MOD_FLAG)$(QE_SOURCE)/upflib          \
    $(MOD_FLAG)$(QE_SOURCE)/Modules         \
    $(MOD_FLAG)$(QE_SOURCE)/MBD             \
    $(MOD_FLAG)$(QE_SOURCE)/FFTXlib/src     \
    $(MOD_FLAG)$(QE_SOURCE)/LAXlib          \
    $(MOD_FLAG)$(QE_SOURCE)/UtilXlib        \
    $(MOD_FLAG)$(QE_SOURCE)/XClib           \
    $(MOD_FLAG)$(QE_SOURCE)/external/devxlib/src \
    $(MOD_FLAG)$(QE_SOURCE)/KS_Solvers      \
    $(MOD_FLAG)$(QE_SOURCE)/PW/src          \
    $(MOD_FLAG)$(QE_SOURCE)/dft-d3          \
    $(MOD_FLAG)$(QE_SOURCE)/LR_Modules      \
    $(MOD_FLAG)../../src

# ---------------------------------------------------------------------------
# Libraries (same order as src/Makefile)
# ---------------------------------------------------------------------------
QEMODS = \
    $(QE_SOURCE)/Modules/libqemod.a              \
    $(QE_SOURCE)/MBD/libmbd.a                    \
    $(QE_SOURCE)/KS_Solvers/Davidson/libdavid.a  \
    $(QE_SOURCE)/KS_Solvers/CG/libcg.a           \
    $(QE_SOURCE)/KS_Solvers/libks_solvers.a      \
    $(QE_SOURCE)/FFTXlib/src/libqefft.a          \
    $(QE_SOURCE)/LAXlib/libqela.a                \
    $(QE_SOURCE)/UtilXlib/libutil.a              \
    $(QE_SOURCE)/upflib/libupf.a                 \
    $(QE_SOURCE)/dft-d3/libdftd3qe.a             \
    $(QE_SOURCE)/LR_Modules/liblrmod.a           \
    $(QE_SOURCE)/XClib/xc_lib.a

PWOBJS  = $(QE_SOURCE)/PW/src/libpw.a
LIBOBJS = $(QE_SOURCE)/UtilXlib/libutil.a

# ---------------------------------------------------------------------------
# Build rules
# ---------------------------------------------------------------------------
.PHONY: all run clean distclean

all: test_runner

run: test_runner
	./test_runner

testdrive.o: testdrive.F90
	$(F90) $(F90FLAGS) $(MODFLAGS) -c -o $@ $<

# Recompile util.f90 locally so we don't depend on the src/ build state
util.o: ../../src/util.f90
	$(F90) $(F90FLAGS) $(MODFLAGS) -c -o $@ $<

test_util.o: test_util.f90 testdrive.o util.o
	$(F90) $(F90FLAGS) $(MODFLAGS) -c -o $@ $<

test_runner.o: test_runner.f90 testdrive.o test_util.o
	$(F90) $(F90FLAGS) $(MODFLAGS) -c -o $@ $<

test_runner: test_runner.o test_util.o util.o testdrive.o \
             $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	    test_runner.o test_util.o util.o testdrive.o \
	    $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

clean:
	-rm -f *.o *.mod test_runner

distclean: clean
