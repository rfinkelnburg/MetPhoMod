############################################################################
#
#   This file is a part of
#
#   MetPhoMod v2.1, the comprehensive three dimensional, prognostic
#		    mesoscale atmospheric summer smog model.
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
#	$(MAKE) opt     -> to build the optimized, sequential version
#	$(MAKE) optp	-> to build the optimized, parallel version
#	$(MAKE) debug	-> to build the sequential, debugging version
#	$(MAKE) debugp	-> to build the parallel, debugging version
#
# Before using this file, change the user variables to values which
# reflect your system.
# On some system this Makefile does not work together with "$(MAKE)", try
# GNU-$(MAKE) "g$(MAKE)"
#
# ***************************************************************************
#
# USER Variables: change settings according to your local installation
#
# ***************************************************************************
# Compiler: The C++-compiler to be used for compilation. The compiler must be
# ANSI-C++ compliant.
CC=xlC -+
#
# for GNU-C
#CC=gcc
#
# ***************************************************************************
# CCOPT    C++ compiler optimization flags. Use maximum optimization. 
#          IEEE floating point compliance
#          is often time consuming but not necessary for MetPhoMod.
#
# Solaris options
#CCOPT=-fast -xO5
#
# Linux/GNU-C options
#CCOPT=-O2 -ffast-math
#
# DEC/UNIX options
#CCOPT=-O2 -std1
#
# AIX options
CCOPT=-O3 -qstrict -qarch=604 -qtune=604 -qflttrap=overflow:zerodivide:invalid
#
# SGI/MIPS-4 options
#CCOPT=-O3 -64
#
# ***************************************************************************
# FOPT compiler optimization flags
FOPT=-O3 -qtune=604
#
# ***************************************************************************
# CCDEBUG       compiler debugging flags
CCDEBUG=-g -qfullpath
#
# ***************************************************************************
# FDEBUG       compiler debugging flags
FDEBUG=-g -qfullpath -qflttrap=overflow:zerodivide:invalid
#
# ***************************************************************************
# CFLAGS2  general C flags, which should always be called.
CFLAGS2=-DSYSTEM=\"`uname`\" -DDATE=\"`date +%D`\" -DAIX -qfloat=nans -qflttrap=overflow:zerodivide:invalid:enable
#
# ***************************************************************************
# compiler flags to be used in all versions
GENFLAGS=
#
# For Solaris:
#GENFLAGS=-ftrap=invalid
# ***************************************************************************
# LEXER program. A lex compatible program
LEX=flex
#
# or
#LEX=flex
#
# ***************************************************************************
# Syntax parser program. Must be yacc compatible
YACC=bison -y
#
# or
#YACC=bison -y
#
# ***************************************************************************
# Libraries to be linked with the sequential version
SEQLIBS=-lnetcdf -lm -lxlf90
#
# for SUN Solaris -lnsl is required, as well.
#SEQLIBS=-lnetcdf -lm -lnsl
#
# ***************************************************************************
# Libraries to be linked with the parallel version
PARLIBS=-lnetcdf -lm -lxlf90 -lpvm3 -lgpvm3
#
# For solaris, -lnsl and -lsocket are required as well.
#PARLIBS=-lnetcdf -lm -lpvm3 -lgpvm3 -lnsl -lsocket
#
# On some parallel machines, -lgpvm3 does not exist
#PARLIBS=-lnetcdf -lm -lpvm3
#
# ***************************************************************************
# NETCDF 3 home directory. This $(MAKE)file assumes that includes are under 
# $(NETCDF)/include and libs are under $(NETCDF)/lib
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
ARCH=AIX
#
############################################################################
#
# NO CHANGES BELOW THIS LINE SHOULD BE NECESSARY!
#
############################################################################

opt :
	(cd mctwostream; \
	$(MAKE) FFLAGS="$(FOPT)" CFLAGS="$(CCOPT) $(CFLAGS2)" GENFLAGS=$(GENFLAGS) OBJDIR=optm)
	$(MAKE) meteochem FFLAGS="$(FOPT)" CFLAGS="$(CCOPT)" \
	   OBJDIR=optm STRIP=strip

optp :
	(cd mctwostream; \
	$(MAKE) FFLAGS="$(FOPT)" CFLAGS="$(CCOPT) $(INCLUDE) -DPARALLEL" \
	   GENFLAGS=$(GENFLAGS) OBJDIR=optmp)
	$(MAKE) parallel FFLAGS="$(FOPT)" CFLAGS="$(CCOPT) -DPARALLEL" \
	   OBJDIR=optmp STRIP=strip ARCH=RS6K

debug :
	(cd mctwostream; \
	$(MAKE) FFLAGS="$(FDEBUG)" CFLAGS="$(CCDEBUG) $(CFLAGS2)" GENFLAGS=$(GENFLAGS) OBJDIR=dbg)
	$(MAKE) meteochem FFLAGS="$(FDEBUG)" CFLAGS="$(CCDEBUG)" OBJDIR=dbg STRIP="ls -l"

debugp :
	(cd mctwostream; \
	$(MAKE) FFLAGS="$(FDEBUG)" CFLAGS="$(CCDEBUG) $(INCLUDE) -DPARALLEL" \
	   GENFLAGS=$(GENFLAGS) OBJDIR=dbgp)
	$(MAKE) parallel FFLAGS="$(FDEBUG)" CFLAGS="$(CCDEBUG) -DPARALLEL" \
	   OBJDIR=dbgp STRIP="ls -l" ARCH=RS6K

# opt:=CFLAGS=-O2 -I/opt/include
# debug:=CFLAGS=-g -I/opt/include

INCLUDE=-I. -I$(NETCDFHOME)/include -I$(PVM_ROOT)/include

#vvvvvvvvvvvvvvvvvvvvvvv replace these lines vvvvvvvvvvvvvvvvvvvvvvv

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
       $(OBJDIR)/mcppm.o \
       $(OBJDIR)/mc_module.o $(OBJDIR)/mc_variable.o $(OBJDIR)/mc_section.o \
       $(OBJDIR)/mc_commands.o $(OBJDIR)/puffemit.o

LIBS = mctwostream/$(OBJDIR)/libmctwostream.a

SRC = mcmain.c mcglobal.c mcdynamics.c mcparse.c mcpress.c mcdep.c \
      sunpos.c mcground.c mcfloat.c mcttt.c mcfilter.c mcemiss.c \
      mcsyntax.c mcinp.c mccdfin.c mcprint.c mchelp.c \
      mcadvect.c mpdata.c mchemsyntax.c mcheminp.c mchem.c mchemparse.c \
      matrix.c mcclouds.c mcgroundclasses.c mcspecial.c mckeps.c mcnest.c \
      mcppm.c \
      mc_module.cc mc_variable.cc mc_section.cc mc_commands.cc \
      puffemit.cc


#vvvvvvvvvvvvvvvvvvvvvvv keep these lines vvvvvvvvvvvvvvvvvvvvvvv

depend :
	makedepend $(INCLUDE) -p'$$(OBJDIR)/' $(SRC)
	(cd mctwostream; $(MAKE) depend)

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(GENFLAGS) $(INCLUDE) -DFERRET_FORMAT $*.c -o $@

$(OBJDIR)/%.o : %.cc
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(GENFLAGS) $(INCLUDE) -DFERRET_FORMAT $*.cc -o $@

$(OBJDIR)/%.o : %.f
	f77 -qextname -c $(FFLAGS) $(GENFLAGS) $*.f -o $@

parallel : $(OBJS) $(OBJDIR)/mcparallel.o 
	$(CC) $(CFLAGS) $(CFLAGS2) $(OBJS) $(OBJDIR)/mcparallel.o -L$(NETCDFHOME)/lib -L$(PVM_ROOT)/lib/$(ARCH) \
	    $(LIBS) $(PARLIBS) -o meteochem
	$(STRIP) meteochem

meteochem : $(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(CFLAGS2) -bloadmap:loadmap $(OBJS) -L$(NETCDFHOME)/lib $(LIBS) $(SEQLIBS) -o meteochem
	$(STRIP) meteochem

mcinp.c : mcinp.l
	$(LEX) -t mcinp.l | sed 's/yy/mci/g' > mcinp.c

mcsyntax.c mcsyntax.h : mcsyntax.y
	$(YACC) -d mcsyntax.y
	sed 's/yy/mci/g' y.tab.c > mcsyntax.c
	rm y.tab.c
	sed 's/yy/mci/g' y.tab.h > mcsyntax.h
	rm y.tab.h

$(OBJDIR)/mchelp.o : mcsyntax.incl

mcsyntax.incl : mcsyntax.y
	sed '1,/%%/d;/{/,/}$$/d;/%%/,$$d;/^$$/d;s/\;$$//;s/^/\"/;s/$$/\\n\"/' mcsyntax.y > mcsyntax.incl

mchemsyntax.c mchemsyntax.h : mchemsyntax.y
	$(YACC) -d mchemsyntax.y
	sed 's/yy/chi/g' y.tab.c > mchemsyntax.c
	rm y.tab.c
	sed 's/yy/chi/g' y.tab.h > mchemsyntax.h
	rm y.tab.h

mcheminp.c : mcheminp.l
	$(LEX) -t mcheminp.l | sed 's/yy/chi/g' > mcheminp.c

# DO NOT DELETE THIS LINE -- $(MAKE) depend depends on it.

$(OBJDIR)/mcmain.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcmain.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcmain.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcmain.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcmain.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcmain.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcmain.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcmain.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcmain.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcmain.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcmain.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcmain.o: /usr/include/math.h mcglobal.h mcpress.h mcdynamics.h
$(OBJDIR)/mcmain.o: mcprint.h mcclouds.h mcttt.h mckeps.h mcfloat.h
$(OBJDIR)/mcmain.o: mcfilter.h mchelp.h mcemiss.h mcadvect.h sunpos.h
$(OBJDIR)/mcmain.o: mcground.h mcparse.h mc_group.hh mccdfin.h mcdep.h
$(OBJDIR)/mcmain.o: mchemparse.h mchem.h mcnest.h mcppm.h mc_module.hh
$(OBJDIR)/mcmain.o: mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mcglobal.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcglobal.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/m_types.h /usr/include/sys/signal.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/localedef31.h /usr/include/stdio.h
$(OBJDIR)/mcglobal.o: /usr/include/va_list.h /usr/include/string.h
$(OBJDIR)/mcglobal.o: /usr/include/malloc.h /usr/include/math.h mcglobal.h
$(OBJDIR)/mcglobal.o: sunpos.h mcground.h mcparse.h mc_group.hh mcprint.h
$(OBJDIR)/mcdynamics.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcdynamics.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcdynamics.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcdynamics.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcdynamics.o: /usr/include/float.h /usr/include/string.h
$(OBJDIR)/mcdynamics.o: /usr/include/math.h /usr/include/malloc.h mcglobal.h
$(OBJDIR)/mcdynamics.o: sunpos.h mcground.h mcpress.h mcdynamics.h mcfilter.h
$(OBJDIR)/mcparse.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcparse.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcparse.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcparse.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcparse.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcparse.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcparse.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcparse.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcparse.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcparse.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcparse.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcparse.o: /usr/include/stdarg.h /usr/include/math.h mcglobal.h
$(OBJDIR)/mcparse.o: sunpos.h mcground.h mcgroundclasses.h mcpress.h
$(OBJDIR)/mcparse.o: mcparse.h mc_group.hh mcemiss.h mcfilter.h mcdep.h
$(OBJDIR)/mcparse.o: mchemparse.h mchem.h mcclouds.h mcnest.h mpdata.h
$(OBJDIR)/mcparse.o: mcppm.h mcrely.h mc_module.hh mc_variable.hh
$(OBJDIR)/mcparse.o: mc_section.hh mc_commands.hh
$(OBJDIR)/mcpress.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcpress.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcpress.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcpress.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcpress.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcpress.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcpress.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcpress.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcpress.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcpress.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcpress.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcpress.o: /usr/include/math.h sunpos.h mcglobal.h mcground.h
$(OBJDIR)/mcpress.o: mcpress.h
$(OBJDIR)/mcdep.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcdep.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcdep.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcdep.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcdep.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcdep.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcdep.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcdep.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcdep.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcdep.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcdep.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcdep.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcdep.o: mcparse.h mc_group.hh mcdep.h
$(OBJDIR)/sunpos.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/sunpos.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/sunpos.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/sunpos.o: /usr/include/float.h /usr/include/math.h sunpos.h
$(OBJDIR)/mcground.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcground.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcground.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcground.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcground.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcground.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcground.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcground.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcground.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcground.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcground.o: /usr/include/sys/localedef31.h /usr/include/math.h
$(OBJDIR)/mcground.o: mcglobal.h mcpress.h sunpos.h mcground.h mcfloat.h
$(OBJDIR)/mcground.o: mctwostream/mctwostream.h mc_module.hh mc_variable.hh
$(OBJDIR)/mcground.o: mc_section.hh mc_commands.hh mctwostream/mcphotodiss.h
$(OBJDIR)/mcground.o: mctwostream/mcapd.h mc_group.hh
$(OBJDIR)/mcfloat.o: /usr/include/signal.h /usr/include/sys/signal.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/mstsave.h /usr/include/sys/m_types.h
$(OBJDIR)/mcfloat.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mcfloat.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/localedef.h /usr/include/sys/lc_core.h
$(OBJDIR)/mcfloat.o: /usr/include/locale.h /usr/include/sys/localedef31.h
$(OBJDIR)/mcfloat.o: /usr/include/string.h /usr/include/unistd.h
$(OBJDIR)/mcfloat.o: /usr/include/standards.h /usr/include/sys/access.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/lockf.h /usr/include/sys/stat.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/mode.h /usr/include/fptrap.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/fp_cpusync.h mcglobal.h mcfloat.h
$(OBJDIR)/mcfloat.o: mcprint.h
$(OBJDIR)/mcttt.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcttt.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcttt.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcttt.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcttt.o: /usr/include/sys/m_types.h /usr/include/sys/signal.h
$(OBJDIR)/mcttt.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcttt.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcttt.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mcttt.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcttt.o: /usr/include/sys/localedef31.h /usr/include/stdio.h
$(OBJDIR)/mcttt.o: /usr/include/va_list.h /usr/include/math.h
$(OBJDIR)/mcttt.o: /usr/include/malloc.h mcglobal.h mcpress.h sunpos.h
$(OBJDIR)/mcttt.o: mcground.h mcttt.h
$(OBJDIR)/mcfilter.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcfilter.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcfilter.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcfilter.o: /usr/include/math.h mcglobal.h mcfilter.h
$(OBJDIR)/mcemiss.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcemiss.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcemiss.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcemiss.o: mcglobal.h sunpos.h mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mcemiss.o: mcemiss.h
$(OBJDIR)/mcsyntax.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcsyntax.o: /usr/include/float.h /usr/include/string.h
$(OBJDIR)/mcsyntax.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/signal.h /usr/include/sys/context.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/m_param.h /usr/include/sys/mstsave.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/localedef.h /usr/include/sys/lc_core.h
$(OBJDIR)/mcsyntax.o: /usr/include/locale.h /usr/include/sys/localedef31.h
$(OBJDIR)/mcsyntax.o: mcglobal.h sunpos.h mcground.h mcgroundclasses.h
$(OBJDIR)/mcsyntax.o: mcparse.h mc_group.hh mccdfin.h mcprint.h mchemparse.h
$(OBJDIR)/mcsyntax.o: mcnest.h mc_commands.hh
$(OBJDIR)/mcinp.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcinp.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcinp.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcinp.o: /usr/include/float.h /usr/include/string.h
$(OBJDIR)/mcinp.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcinp.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcinp.o: /usr/include/sys/signal.h /usr/include/sys/context.h
$(OBJDIR)/mcinp.o: /usr/include/sys/m_param.h /usr/include/sys/mstsave.h
$(OBJDIR)/mcinp.o: /usr/include/sys/localedef.h /usr/include/sys/lc_core.h
$(OBJDIR)/mcinp.o: /usr/include/locale.h /usr/include/sys/localedef31.h
$(OBJDIR)/mcinp.o: mcglobal.h sunpos.h mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mcinp.o: mcprint.h mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mcinp.o: mc_commands.hh mcsyntax.h
$(OBJDIR)/mccdfin.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mccdfin.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/m_types.h /usr/include/sys/signal.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/localedef31.h /usr/include/stdio.h
$(OBJDIR)/mccdfin.o: /usr/include/va_list.h /usr/include/string.h
$(OBJDIR)/mccdfin.o: /usr/local/include/netcdf.h /usr/include/errno.h
$(OBJDIR)/mccdfin.o: mcglobal.h sunpos.h mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mccdfin.o: mccdfin.h mcemiss.h
$(OBJDIR)/mcprint.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcprint.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcprint.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcprint.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcprint.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcprint.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcprint.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcprint.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcprint.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcprint.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcprint.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcprint.o: /usr/include/math.h /usr/local/include/netcdf.h
$(OBJDIR)/mcprint.o: /usr/include/errno.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcprint.o: mcparse.h mc_group.hh mcprint.h mchemparse.h mchem.h
$(OBJDIR)/mchelp.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mchelp.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mchelp.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mchelp.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mchelp.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mchelp.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mchelp.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mchelp.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mchelp.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mchelp.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mchelp.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mchelp.o: mcglobal.h sunpos.h mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mchelp.o: mchelp.h mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mchelp.o: mc_commands.hh mctwostream/mctwostream.h
$(OBJDIR)/mchelp.o: mctwostream/mcphotodiss.h mctwostream/mcapd.h
$(OBJDIR)/mchelp.o: mcsyntax.incl
$(OBJDIR)/mcadvect.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcadvect.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcadvect.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcadvect.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcadvect.o: mchemparse.h mpdata.h
$(OBJDIR)/mpdata.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mpdata.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mpdata.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mpdata.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mpdata.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mpdata.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mpdata.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mpdata.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mpdata.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mpdata.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mpdata.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mpdata.o: /usr/include/stdarg.h /usr/include/math.h mcglobal.h
$(OBJDIR)/mpdata.o: mpdata.h
$(OBJDIR)/mchemsyntax.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mchemsyntax.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mchemsyntax.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/context.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/m_param.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/localedef.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mchemsyntax.o: mcglobal.h mchemparse.h mctwostream/mctwostream.h
$(OBJDIR)/mchemsyntax.o: mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mchemsyntax.o: mc_commands.hh mctwostream/mcphotodiss.h
$(OBJDIR)/mchemsyntax.o: mctwostream/mcapd.h mc_group.hh
$(OBJDIR)/mcheminp.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcheminp.o: /usr/include/float.h /usr/include/string.h
$(OBJDIR)/mcheminp.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/signal.h /usr/include/sys/context.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/m_param.h /usr/include/sys/mstsave.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/localedef.h /usr/include/sys/lc_core.h
$(OBJDIR)/mcheminp.o: /usr/include/locale.h /usr/include/sys/localedef31.h
$(OBJDIR)/mcheminp.o: mcglobal.h mchemsyntax.h mchemparse.h
$(OBJDIR)/mchem.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mchem.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mchem.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mchem.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mchem.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mchem.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mchem.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mchem.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mchem.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mchem.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mchem.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mchem.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mchem.o: mchemparse.h mchem.h matrix.h mcparallel.h
$(OBJDIR)/mchemparse.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mchemparse.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mchemparse.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/localedef.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mchemparse.o: /usr/include/math.h matrix.h mcglobal.h sunpos.h
$(OBJDIR)/mchemparse.o: mcground.h mcparse.h mc_group.hh mchemparse.h
$(OBJDIR)/matrix.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/matrix.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/matrix.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/matrix.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/matrix.o: /usr/include/sys/m_types.h /usr/include/sys/signal.h
$(OBJDIR)/matrix.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/matrix.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/matrix.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/matrix.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/matrix.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/matrix.o: /usr/include/math.h matrix.h
$(OBJDIR)/mcclouds.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcclouds.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcclouds.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcclouds.o: /usr/include/math.h sunpos.h mcglobal.h mcground.h
$(OBJDIR)/mcclouds.o: mcpress.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/types.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/inttypes.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/m_types.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/resource.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/time.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/signal.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/context.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/m_param.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/localedef.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/localedef31.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/string.h mcglobal.h sunpos.h
$(OBJDIR)/mcgroundclasses.o: mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mcgroundclasses.o: mcgroundclasses.h
$(OBJDIR)/mcspecial.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcspecial.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/m_types.h /usr/include/sys/signal.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/localedef.h /usr/include/sys/limits.h
$(OBJDIR)/mcspecial.o: /usr/include/float.h /usr/include/sys/lc_core.h
$(OBJDIR)/mcspecial.o: /usr/include/locale.h /usr/include/sys/localedef31.h
$(OBJDIR)/mcspecial.o: /usr/include/math.h
$(OBJDIR)/mckeps.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mckeps.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mckeps.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mckeps.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mckeps.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mckeps.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mckeps.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mckeps.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mckeps.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mckeps.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mckeps.o: /usr/include/sys/localedef31.h /usr/include/math.h
$(OBJDIR)/mckeps.o: mcglobal.h sunpos.h mcground.h mchemparse.h
$(OBJDIR)/mcnest.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcnest.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcnest.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcnest.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcnest.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcnest.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcnest.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcnest.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcnest.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcnest.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcnest.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcnest.o: /usr/include/math.h /usr/local/include/netcdf.h
$(OBJDIR)/mcnest.o: /usr/include/errno.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcnest.o: mcparse.h mc_group.hh mcnest.h mcprint.h
$(OBJDIR)/mcppm.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcppm.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcppm.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcppm.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcppm.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcppm.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcppm.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcppm.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcppm.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcppm.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcppm.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcppm.o: /usr/include/math.h mcglobal.h mcppm.h sunpos.h mcground.h
$(OBJDIR)/mcppm.o: mchemparse.h
$(OBJDIR)/mc_module.o: mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mc_module.o: mc_commands.hh /usr/include/stddef.h
$(OBJDIR)/mc_module.o: /usr/include/standards.h /usr/include/stdio.h
$(OBJDIR)/mc_module.o: /usr/include/va_list.h /usr/include/sys/types.h
$(OBJDIR)/mc_module.o: /usr/include/sys/inttypes.h /usr/include/sys/m_types.h
$(OBJDIR)/mc_module.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mc_module.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mc_module.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mc_module.o: /usr/include/sys/signal.h /usr/include/sys/context.h
$(OBJDIR)/mc_module.o: /usr/include/sys/m_param.h /usr/include/sys/mstsave.h
$(OBJDIR)/mc_module.o: /usr/include/sys/localedef.h
$(OBJDIR)/mc_module.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mc_module.o: /usr/include/sys/localedef31.h mcparse.h mc_group.hh
$(OBJDIR)/mc_variable.o: mc_variable.hh /usr/include/stddef.h
$(OBJDIR)/mc_variable.o: /usr/include/standards.h /usr/include/stdio.h
$(OBJDIR)/mc_variable.o: /usr/include/va_list.h /usr/include/sys/types.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/inttypes.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mc_variable.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/context.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/m_param.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/localedef.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mc_variable.o: mcground.h mc_group.hh mcparse.h mc_group.t
$(OBJDIR)/mc_section.o: mc_section.hh /usr/include/stddef.h
$(OBJDIR)/mc_section.o: /usr/include/standards.h /usr/include/stdio.h
$(OBJDIR)/mc_section.o: /usr/include/va_list.h /usr/include/sys/types.h
$(OBJDIR)/mc_section.o: /usr/include/sys/inttypes.h
$(OBJDIR)/mc_section.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mc_section.o: /usr/include/float.h /usr/include/string.h mcglobal.h
$(OBJDIR)/mc_section.o: mcrely.h mcnest.h
$(OBJDIR)/mc_commands.o: mc_commands.hh /usr/include/stddef.h
$(OBJDIR)/mc_commands.o: /usr/include/standards.h /usr/include/stdio.h
$(OBJDIR)/mc_commands.o: /usr/include/va_list.h /usr/include/sys/types.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/inttypes.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mc_commands.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/context.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/m_param.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/localedef.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mc_commands.o: mc_group.t mcparse.h mc_group.hh
$(OBJDIR)/puffemit.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/puffemit.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/puffemit.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/puffemit.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/puffemit.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/puffemit.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/puffemit.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/puffemit.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/puffemit.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/puffemit.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/puffemit.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/puffemit.o: mcglobal.h mcground.h mcparse.h mc_group.hh
$(OBJDIR)/puffemit.o: mc_section.hh mc_group.t mc_commands.hh mc_variable.hh
$(OBJDIR)/puffemit.o: mc_module.hh