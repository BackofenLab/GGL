
################################################################
# pkg-config
################################################################

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ggl.pc

################################################################

include $(top_srcdir)/aminclude.am

doc: doxygen-doc
if DX_COND_html
	cp -f $(top_srcdir)/misc/logo/GGL-logo-100x100.png $(top_builddir)/$(DX_DOCDIR)/html/.
endif

# directory that contains all tests
TSTDIR = tests
			
tests: all
	@$(PRINTF) "-------------------------------------------------------------------------------\n"
	@$(PRINTF) " Running the tests :\n"
	@$(PRINTF) "-------------------------------------------------------------------------------\n"
	@make tests -s -C $(TSTDIR)
	@$(PRINTF) "-------------------------------------------------------------------------------\n"

test: tests


SUBDIRS = src misc data perl $(TSTDIR)
