UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TOPDIR = $(shell pwd)

TARGET = OceanWave
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

SRCDIRs = $(TOPDIR)/SwanSrc
SOURCES = ShearVelocity.f90 SimWaves.f90  OceanOut.f90

OBJS=./*.o

FOBS = $(SOURCES:.f90=.o)

.PHONY: dist module dust clobber

dist : module

module :
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for Ocean Waves"
	@echo
	@echo "*************************************************"
	for i in $(SRCDIRs) ; do   \
    	( cd $$i ; make dist) ;       \
	done
	make $(TARGETLIB)

$(TARGETLIB): $(FOBS)
	$(AR) $(TARGETLIB) $(OBJDIR)/*.o
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) -c $< -o $@

dust :
	#rm -f $(SRCDIRs)/*.o $(SRCDIRs)/*.mod $(SRCDIRs)/*.original
	#rm -f $(SRCDIRs)/*~ $(SRCDIRs)/*_genmod*
	rm -fv *.bak *~ *.o *.mod *.original *_genmod*

clobber : dust
