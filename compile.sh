#!/bin/bash

# script to get Dictionary.o into the shared library
# since the dependencies in Geant4 make it such that
# the Dictionary.o is compiled after the shared 
# library has been created, we delete the shared library,
# and the make it again (now that Dictionary.o has been
# created)

# this scripts passes all commandline flags onto the make command
# (mainly intended for -jN flags)
make $@ && make clean_libs && make $@
