############################################################################
#
#   This file is a part of
#
#   MetPhoMod v2.0, the comprehensive three dimensional, prognostic
#		   mesoscale atmospheric summer smog model.
#
#   Copyright (C) 1996-1999, Silvan Perego
#
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, refer to the WEB page
#
#   http://www.gnu.org/copyleft/gpl.html
#
############################################################################
#
#	MetPhoMod Makefile
#
#	use this file by either typing:
#
#	make opt        -> to build the optimized, sequential version
#	make optp	-> to build the optimized, parallel version
#	make debug	-> to build the sequential, debugging version
#	make debugp	-> to build the parallel, debugging version
#
# Before using this file, change the user variables to values which
# reflect your system.
# On some system this Makefile does not work together with "make", try
# GNU-make "gmake"
#
# ***************************************************************************
#
# USER Variables: change settings according to your local installation
#
# ***************************************************************************
# Compiler: The C-compiler to be used for compilation. The compiler must be
# ANSI-C compliant.
CC=cc
#
# for GNU-C
#CC=gcc
#
# ***************************************************************************
# OPTIMIZE    compiler optimization flags. Use maximum optimization. 
#             IEEE floating point compliance
#             is often time consuming but not necessary for MetPhoMod.
#
# Solaris options
#OPTIMIZE=-fast -xO5
#
# Linux/GNU-C options
#OPTIMIZE=-O2 -ffast-math
#
# DEC/UNIX options
#OPTIMIZE=-O2 -std1
#
# AIX options
OPTIMIZE=-O3 -qtune=pwrx -qflttrap
#
# SGI/MIPS-4 options
#OPTIMIZE=-O3 -64
#
# ***************************************************************************
# DEBUG       compiler debugging flags
DEBUG=-g
#
# ***************************************************************************
# C compiler flags to be used in all versions
#GENFLAGS=
#
# For Solaris:
GENFLAGS=-ftrap=invalid
# ***************************************************************************
# LEXER program. A lex compatible program
LEX=lex
#
# or
#LEX=flex
#
# ***************************************************************************
# Syntax parser program. Must be yacc compatible
YACC=yacc
#
# or
#YACC=bison -y
#
# ***************************************************************************
# Libraries to be linked with the sequential version
SEQLIBS=-lnetcdf -lm
#
# for SUN Solaris -lnsl is required, as well.
#SEQLIBS=-lnetcdf -lm -lnsl
#
# ***************************************************************************
# Libraries to be linked with the parallel version
PARLIBS=-lnetcdf -lm -lpvm3 -lgpvm3
#
# For solaris, -lnsl and -lsocket are required as well.
#PARLIBS=-lnetcdf -lm -lpvm3 -lgpvm3 -lnsl -lsocket
#
# On some parallel machines, -lgpvm3 does not exist
#PARLIBS=-lnetcdf -lm -lpvm3
#
# ***************************************************************************
# NETCDF 3 home directory. This makefile assumes that includes are under 
# $(NETCDFHOME)/include and libs are under $(NETCDFHOME)/lib
NETCDFHOME=/usr/local
#
# ***************************************************************************
# PVM 3.3 home directory-> set the PVM_ROOT environment variable within your
# shell!
# (These variables have no effect, when building the sequential version.)
#PVM_ROOT=$(HOME)/pvm3
#
# PVM 3.3 system architecture. Refer to the PVM3 documentation
#
ARCH=SUN4SOL2
#
############################################################################
#
# NO CHANGES BELOW THIS LINE SHOULD BE NECESSARY!
#
############################################################################


opt :
	- rm -f .debug
	$(MAKE) meteochem CFLAGS="$(OPTIMIZE) -DSYSTEM=\\\""`uname`"\\\" -D"`uname`" -DDATE=\\\""`date +%D`"\\\" "$(DEFS) \
	   OBJDIR=optm STRIP=strip

optp :
	- rm -f .debug
	$(MAKE) parallel CFLAGS="$(OPTIMIZE) -DPARALLEL -DSYSTEM=\\\""`uname`"\\\" -D"`uname`" -DDATE=\\\""`date +%D`"\\\" "$(DEFS) \
	   OBJDIR=optmp STRIP=strip

debug :
	- rm -f .nodebug
	$(MAKE) meteochem CFLAGS="$(DEBUG) -DSYSTEM=\\\""`uname`"\\\" -D"`uname`" -DDATE=\\\""`date +%D`"\\\" "$(DEFS) OBJDIR=dbg STRIP=ls

debugp :
	- rm -f .debug
	$(MAKE) parallel CFLAGS="$(DEBUG) -DPARALLEL -DSYSTEM=\\\""`uname`"\\\" -D"`uname`" -DDATE=\\\""`date +%D`"\\\" "$(DEFS) \
	   OBJDIR=dbgp STRIP=ls

CFLAGS2=$(GENFLAGS) -I$(NETCDFHOME)/include -I$(PVM_ROOT)/include

OBJS = $(OBJDIR)/mcmain.o $(OBJDIR)/mcglobal.o $(OBJDIR)/mcdynamics.o \
       $(OBJDIR)/mcparse.o $(OBJDIR)/mcpress.o $(OBJDIR)/mcdep.o \
       $(OBJDIR)/sunpos.o $(OBJDIR)/mcground.o $(OBJDIR)/mcfloat.o \
       $(OBJDIR)/mcttt.o $(OBJDIR)/mcfilter.o $(OBJDIR)/mcemiss.o \
       $(OBJDIR)/mcsyntax.o $(OBJDIR)/mcinp.o \
       $(OBJDIR)/mccdfin.o $(OBJDIR)/mcprint.o $(OBJDIR)/mchelp.o \
       $(OBJDIR)/mcadvect.o $(OBJDIR)/mpdata.o $(OBJDIR)/mchemsyntax.o \
       $(OBJDIR)/mcheminp.o $(OBJDIR)/mchem.o $(OBJDIR)/mchemparse.o \
       $(OBJDIR)/matrix.o $(OBJDIR)/mcclouds.o $(OBJDIR)/mcgroundclasses.o \
       $(OBJDIR)/mcspecial.o $(OBJDIR)/mckeps.o $(OBJDIR)/mcnest.o \
       $(OBJDIR)/mcppm.o

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(CFLAGS2) $*.c -o $@

parallel : $(OBJS) $(OBJDIR)/mcparallel.o 
	$(CC)  $(OBJS) $(OBJDIR)/mcparallel.o -L$(NETCDFHOME)/lib -L$(PVM_ROOT)/lib/$(ARCH) \
	    $(PARLIBS) -o meteochem
	$(STRIP) meteochem

meteochem : $(OBJS)
	$(CC)  $(CFLAGS) $(OBJS) -L$(NETCDFHOME)/lib $(SEQLIBS) -o meteochem
	$(STRIP) meteochem

$(OBJDIR)/mcmain.o : mcmain.c mcglobal.h mcdynamics.h mcparse.h mcemiss.h mcdep.h \
           mcfloat.h mcttt.h mcprint.h mchelp.h mcadvect.h mccdfin.h \
           sunpos.h mcground.h mchemparse.h mchem.h mcparallel.h mcclouds.h \
           mckeps.h mcnest.h mcppm.h

$(OBJDIR)/mcglobal.o : mcglobal.c mcglobal.h mcground.h mcparse.h sunpos.h mcparse.h \
		       mcprint.h

$(OBJDIR)/mcdynamics.o : mcdynamics.c mcdynamics.h mcglobal.h mcground.h sunpos.h mcparallel.h

$(OBJDIR)/mcparse.o : mcparse.c mcparse.h mcglobal.h mcground.h sunpos.h mcdep.h \
            mchemparse.h mchem.h mcclouds.h mcgroundclasses.h mpdata.h mcparallel.h \
            mcnest.h mcppm.h mcemiss.h

$(OBJDIR)/mcpress.o : mcpress.c mcpress.h mcglobal.h sunpos.h mcground.h mcparallel.h

$(OBJDIR)/sunpos.o : sunpos.c sunpos.h

$(OBJDIR)/mcground.o : mcground.c mcground.h mcglobal.h sunpos.h

$(OBJDIR)/mcfloat.o : mcfloat.c mcfloat.h mcglobal.h mcprint.h mcparallel.h

$(OBJDIR)/mcttt.o : mcttt.c mcttt.h mcglobal.h mcpress.h mcground.h sunpos.h

$(OBJDIR)/mckeps.o : mckeps.c mckeps.h mcglobal.h sunpos.h mcground.h mchemparse.h

$(OBJDIR)/mcfilter.o : mcfilter.c mcfilter.h mcglobal.h

$(OBJDIR)/mcinp.o : mcinp.c mcglobal.h sunpos.h mcground.h mcsyntax.h mcparse.h mcprint.h

$(OBJDIR)/mcsyntax.o : mcsyntax.c mcsyntax.h mcglobal.h mcground.h sunpos.h mcparse.h \
             mccdfin.h mcprint.h mchemparse.h mcnest.h

$(OBJDIR)/mccdfin.o : mccdfin.c mcglobal.h sunpos.h mcground.h mcparse.h mccdfin.h \
		      mcparallel.h mcemiss.h

mcinp.c : mcinp.l
	$(LEX) -t mcinp.l | sed 's/yy/mci/g' > mcinp.c

mcsyntax.c mcsyntax.h : mcsyntax.y
	$(YACC) -d mcsyntax.y
	sed 's/yy/mci/g' y.tab.c > mcsyntax.c
	rm y.tab.c
	sed 's/yy/mci/g' y.tab.h > mcsyntax.h
	rm y.tab.h

$(OBJDIR)/mcprint.o : mcprint.c mcglobal.h sunpos.h mcground.h mcparse.h mcprint.h mcparallel.h \
              mchemparse.h mchem.h

$(OBJDIR)/mchelp.o : mchelp.c mcglobal.h sunpos.h mcground.h mcparse.h mchelp.h \
           mcsyntax.incl

mcsyntax.incl : mcsyntax.y
	sed '1,/%%/d;/{/,/}$$/d;/%%/,$$d;/^$$/d;s/\;$$//;s/^/\"/;s/$$/\\n\"/' mcsyntax.y > mcsyntax.incl

$(OBJDIR)/mpdata.o : mpdata.c mpdata.h mcglobal.h mcparallel.h

$(OBJDIR)/mcadvect.o : mcadvect.c mcadvect.h mpdata.h mcglobal.h sunpos.h mcground.h mchemparse.h

$(OBJDIR)/mchemparse.o : mchemparse.c mchemparse.h matrix.h mcglobal.h mcparse.h mchem.h sunpos.h mcground.h

$(OBJDIR)/mchem.o : mchem.c matrix.h mcglobal.h mcparse.h mchemparse.h mchem.h \
		    sunpos.h mcground.h mcparallel.h

$(OBJDIR)/mcclouds.o : mcclouds.c mcclouds.h mcglobal.h mcground.h mcpress.h

$(OBJDIR)/mchemsyntax.o : mchemsyntax.c mchemsyntax.h mchemparse.h mcglobal.h mcgroundclasses.h

$(OBJDIR)/mcheminp.o : mcheminp.c mchemsyntax.h mchemparse.h mcglobal.h

$(OBJDIR)/mcemiss.o : mcemiss.c mcemiss.h mcglobal.h sunpos.h mcground.h mcparse.h mcparallel.h

$(OBJDIR)/mcdep.o : mcdep.c mcdep.h mcglobal.h mcparse.h sunpos.h mcground.h

$(OBJDIR)/mcparallel.o : mcparallel.c mcparallel.h mcglobal.h sunpos.h mcground.h \
                         mccdfin.h mchem.h mcemiss.h mcdep.h mchemparse.h mcpress.h

$(OBJDIR)/mcgroundclasses.o : mcgroundclasses.c mcglobal.h sunpos.h mcground.h \
			      mcparse.h mcgroundclasses.h

mchemsyntax.c mchemsyntax.h : mchemsyntax.y
	$(YACC) -d mchemsyntax.y
	sed 's/yy/chi/g' y.tab.c > mchemsyntax.c
	rm y.tab.c
	sed 's/yy/chi/g' y.tab.h > mchemsyntax.h
	rm y.tab.h

mcheminp.c : mcheminp.l
	$(LEX) -t mcheminp.l | sed 's/yy/chi/g' > mcheminp.c

$(OBJDIR)/matrix.o : matrix.c matrix.h

$(OBJDIR)/mcspecial.o : mcspecial.c

$(OBJDIR)/mcnest.o : mcnest.c mcnest.h mcglobal.h mcparse.h mcground.h sunpos.h \
                     mcparallel.h mcprint.h

$(OBJDIR)/mcppm.o : mcppm.c mcppm.h mcglobal.h sunpos.h mcground.h mchemparse.h \
		    mcparallel.h
