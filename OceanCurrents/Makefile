UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = OceanCurrent
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES =  Interpolation.f90 CurrentClass.f90 InputReader.f90 \
					ParallelCirc.f90 CurrentIO.f90 CurrentCirc.f90

FOBJS=$(FSOURCES:.f90=.o)

OBJS=./*.o

.PHONY: dist module dust clobber

dist : module

module :
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for Ocean Currents"
	@echo
	@echo "*************************************************"
	make $(TARGETLIB)

$(TARGETLIB) : $(FOBJS)
	$(AR) $(TARGETLIB) $(OBJS)
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) -c $< -o $@

dust :
	rm -f *.o *.mod *.original
	rm -f *~ *_genmod*

clobber : dust
