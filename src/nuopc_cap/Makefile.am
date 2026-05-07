AUTOMAKE_OPTIONS = subdir-objects

AM_FCFLAGS = @ESMFFCFLAGS@

lib_LTLIBRARIES = libISSMNUOPCCap.la

libISSMNUOPCCap_la_SOURCES = issm_nuopc_cap.F90
libISSMNUOPCCap_la_LIBADD = ../c/libISSMCore.la ../c/libISSMModules.la $(ESMFLIB)

if !MSYS2
libISSMNUOPCCap_la_LIBADD += ../c/libISSMOverload.la
endif

if FORTRAN
libISSMNUOPCCap_la_LIBADD += $(FLIBS) $(FORTRANLIB)
endif

EXTRA_DIST = README.md
