# --------------------------------------------------------------
# GNUmakefile for MiniBall simluation example Ex04
# Gabriele Cosmo, 06/04/98. Modifications by Hans Boie 
# --------------------------------------------------------------
#CPPVERBOSE := true

name := TRexGeant4
G4TARGET := $(name)
G4EXLIB := true

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
#ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := -I$(shell root-config --incdir)

COMMON_DIR = $(HOME)/CommandLineInterface
LIB_DIR = $(HOME)/lib

MBHOME := $(HOME)/Miniball

MBINC := $(MBHOME)/include
ifndef MBLIB
MBLIB := $(MBHOME)/lib
endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all cleantmp cleandir tar
all: mblib lib bin

include $(G4INSTALL)/config/architecture.gmk

EXTRA_INCLUDES = -I$(MBINC) -I$(COMMON_DIR)
CPPFLAGS += -O0 -Wall -Wno-write-strings -Wno-deprecated -pedantic -g $(INCFLAGS) $(ROOTCFLAGS) $(EXTRA_INCLUDES) -DSIMULATION_PATH=\"$(PWD)\"
CPPFLAGS += -pthread
CPPFLAGS += -DG4MULTITHREADED
LDFLAGS  += $(ROOTLIBS) -Wl,--no-as-needed $(G4TMP)/$(G4SYSTEM)/$(G4TARGET)/Dictionary.o -lG4OpenGL
EXTRALIBS = -L$(MBLIB) -L$(LIB_DIR) -lMiniBall -lCommandLineInterface


EXTRA_LINK_DEPENDENCIES = $(G4TMP)/$(G4SYSTEM)/$(G4TARGET)/Dictionary.o

include $(G4INSTALL)/config/binmake.gmk

#------------- ROOT -------------#

DEPENDENCIES = include/TRexSettings.hh \
	include/Barrel.hh \
	include/Annular.hh \
	include/Germanium.hh \
	include/ParticleMC.hh \
	RootLinkDef.h


$(G4TMP)/$(G4SYSTEM)/$(G4TARGET)/Dictionary.o: Dictionary.cc
	@echo "... generating Dictionary.o"
	#$(CC) -fPIC $(CFLAGS) $(ROOTINC) -o $@ -c $< original
	$(CC) -fPIC $(CPPFLAGS) $(ROOTINC) -o $@ -c $<

Dictionary.cc: $(DEPENDENCIES)
	@echo "... updating Dictionary"
	rm -f Dictionary.*; rootcint -f $@ -c $(EXTRA_INCLUDES) -I$(CLHEP_INCLUDE_DIR) $(DEPENDENCIES)

#-------------  MB  _____________#

mblib:
	@echo "...updating libMiniBall"
	@cd $(MBHOME); make mblib

cleantmp:
	@echo "cleaning tmp files from " $(G4TMP)/$(G4SYSTEM)
	@cd $(G4TMP)/$(G4SYSTEM)/; rm $(G4TARGET)/*.d; rm $(G4TARGET)/*.o; rm MiniBall/*.d; rm MiniBall/*.o;


#---------- other stuff (e.g. tar-file) ----------

cleandir:
	@echo "cleaning unneccesary files from " $(PWD)
	rm AngularDistribution.root; rm OpticalPotential.root; rm test.root;

tar:
	@echo "creating zipped tar-ball ... "
	tar -chvzf Transfer.tar.gz src/*.cc include/*.hh Transfer.cc GNUmakefile RootLinkDef.h cmass.dat LevelFiles/DefaultLevels.dat Settings/Mg32Settings.dat AngularDistributions/79Zn-test.dat OGL_vis.mac OGLS_vis.mac ReadMe
