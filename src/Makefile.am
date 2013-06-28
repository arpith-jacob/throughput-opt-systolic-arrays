# Copyright (C) 2009 Arpith C. Jacob <jarpith@cse.wustl.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

INCLUDES = -I$(top_srcdir)/include
LDADD = $(top_srcdir)/lib.linux64/libbarvinok-core.a $(top_srcdir)/lib.linux64/libbarvinok.a $(top_srcdir)/lib.linux64/libpiplibMP.a $(top_srcdir)/lib.linux64/libpolylibgmp.a $(top_srcdir)/lib.linux64/libntl.a -lgmp -lboost_program_options -lboost_filesystem

bin_PROGRAMS = enumerate
enumerate_SOURCES = enumerate-projections.cpp commandline-options.hpp  \
                    polyhedron-options.hpp   projection-solver.hpp  \
                    solutions.hpp   index-enumerator.hpp     \
                    projection-solution.hpp  schedule-ilp.hpp     \
                    throughput-ilp.hpp



