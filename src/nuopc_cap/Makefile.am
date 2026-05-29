AUTOMAKE_OPTIONS = subdir-objects

AM_FCFLAGS = @ESMFFCFLAGS@

lib_LTLIBRARIES = libISSMNUOPCCap.la

libISSMNUOPCCap_la_SOURCES = issm_nuopc_cap.F90
libISSMNUOPCCap_la_LIBADD = ../c/libISSMCore.la $(ESMFLIB)

if !MSYS2
libISSMNUOPCCap_la_LIBADD += ../c/libISSMOverload.la
endif

if FORTRAN
libISSMNUOPCCap_la_LIBADD += $(FLIBS) $(FORTRANLIB)
endif

EXTRA_DIST = README.md

install-data-hook:
	$(MKDIR_P) $(DESTDIR)$(includedir)
	$(INSTALL_DATA) $(srcdir)/../c/main/issm_nuopc_capi.h $(DESTDIR)$(includedir)/issm_nuopc_capi.h
	@set -e; \
	for mod in $(builddir)/*.mod $(builddir)/*.MOD; do \
		if test -f "$$mod"; then \
			$(INSTALL_DATA) "$$mod" $(DESTDIR)$(includedir)/`basename "$$mod"`; \
		fi; \
	done

uninstall-hook:
	rm -f $(DESTDIR)$(includedir)/issm_nuopc_capi.h
	@set -e; \
	for mod in $(builddir)/*.mod $(builddir)/*.MOD; do \
		if test -f "$$mod"; then \
			rm -f $(DESTDIR)$(includedir)/`basename "$$mod"`; \
		fi; \
	done
