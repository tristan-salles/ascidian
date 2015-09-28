CFGFILE = $(TOP)/config/Makefile.cfg.osx
include $(CFGFILE)

DIRPLUG = plugin

EXEC=bin/$(EXE)

BUILDDIR=$(TOP)/build
OBJDIR=$(BUILDDIR)/objects
MODDIR=$(BUILDDIR)/modules
LIBDIR=$(BUILDDIR)/lib

# C compiler to use
MPIC= $(CC) 
# C++ compiler to use
MPICPP= $(CPP) 
# FC compiler to use
F90=$(FC) 

FFLAGS= $(FCFLAGS) -I$(MODDIR)
	
LDFLAGS= -L$(LIBDIR)/

FOXFLAGS= $$($(FOX)/FoX-config --fcflags)
LDFOXFLAGS=$$($(FOX)/FoX-config --libs)

