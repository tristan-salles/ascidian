TOP=$(shell pwd)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

DIRMODS= OceanCurrents OceanWaves

SOURCES = SimOcean.f90
OBJS=$(SOURCES:.f90=.o)

.PHONY : all dist plugin dust clobber

all: dist

dist:
	@echo
	@echo "*************************************************"
	@echo "Ocean Plugin Makefile"
	@echo "*************************************************"
	@echo
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(OBJDIR)
	@mkdir -p $(MODDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p bin
	for i in $(DIRMODS) ; do   \
    	  ( cd $$i ; make dist) ;       \
	done
	@echo "*************************************************"
	@echo
	@echo "Build Ocean binary."
	@echo
	@echo "*************************************************"
	@$(if $(wildcard SimOcean.o),rm -f SimOcean.o,)
	make $(EXEC)

plugin :
	cd $(DIRPLUG); make plugin;
	@echo "*************************************************"
	@echo
	@echo "Ocean shared library created."
	@echo
	@echo "*************************************************"

$(EXEC) :	$(OBJS)
	$(F90) $(FFLAGS)  $(FOXFLAGS) $(H5FLAGS) -o $@ $^ $(LDFLAGS) -lOcean  $(H5LDFLAGS) $(H5LIBS) $(LDFOXFLAGS)
	#rm *.o
	@echo "*************************************************"
	@echo
	@echo "Ocean updated in ./bin/."
	@echo
	@echo "*************************************************"

%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) -c $< -o $@
	$(AR) $(LIBDIR)/libOcean.a $(OBJDIR)/*.o

dust :
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make dust) ;       \
	done
	$(foreach module,$(MODULES),cd $(module); cd - ; )
	rm -fv *~ *.bak *.o *.mod *.original

clobber : dust
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make clobber) ;   \
	done
	rm -rfv $(BUILDDIR)
