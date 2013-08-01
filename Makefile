#vvvvvvvvvvvvvvvvvvvvvvv keep these lines vvvvvvvvvvvvvvvvvvvvvvv

all : optp

normal :
	if [ -f .nodebug ] ; \
	then make opt ; \
	else make debug ; \
	fi ;

opt :
	- rm -f .debug
	touch .nodebug
	(cd mctwostream; make opt $(DEFS))
	make meteochem FFLAGS=-O3 CFLAGS="-O3 -ffast-math "$(DEFS) \
	   STRIP=ls

optp :
	- rm -f .debug
	touch .nodebug
	(cd mctwostream && make libmctwostream.a)
	make parallel FFLAGS=-O3 CFLAGS="-O2 -ffast-math -I. -I/usr/local/include -I/home/perego/pvm3/include -DLINUX -DPARALLEL -DSYSTEM=\\\""`uname`"\\\" -DDATE=\\\""`date +%D`"\\\" "$(DEFS) \
	   STRIP=strip ARCH=LINUX

debug :
	- rm -f .nodebug
	touch .debug
	(cd mctwostream; make debug)
	make meteochem FFLAGS=-g CFLAGS="-g "$(DEFS) STRIP=ls

debugp :
	- rm -f .debug
	touch .nodebug
	make parallel FFLAGS=-g CFLAGS="-g -I. -I/usr/local/include -I/home/perego/pvm3/include -DLINUX -DPARALLEL -DSYSTEM=\\\""`uname`"\\\" -DDATE=\\\""`date +%D`"\\\" "$(DEFS) \
	   STRIP=ls ARCH=LINUX

clean :
	rm -f $(OBJS)
	(cd mctwostream; make clean)

# opt:=CFLAGS=-O2 -I/opt/include
# debug:=CFLAGS=-g -I/opt/include

CC=g++
INCLUDE=-I. -I/usr/local/include -I/home/perego/pvm3/include
CFLAGS2=-DSYSTEM=\"`uname`\" -DDATE=\"`date +%D`\" -DLINUX -fno-implicit-templates -fpermissive

OBJDIR?= .

$(OBJDIR)/mchelp.o: mcsyntax.incl

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

LIBS = mctwostream/libmctwostream.a

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
	(cd mctwostream; make depend)

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(INCLUDE) -DFERRET_FORMAT $*.c -o $@

$(OBJDIR)/%.o : %.cc
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(INCLUDE) -DFERRET_FORMAT $*.cc -o $@

$(OBJDIR)/%.o : %.f
	f77 -c $(FFLAGS) $*.f -o $@

parallel : $(OBJS) $(OBJDIR)/mcparallel.o 
	$(CC)  $(OBJS) $(OBJDIR)/mcparallel.o -L/usr/local/lib -L/home/perego/pvm3/lib/$(ARCH) \
	    -lm -lnetcdf -lpvm3 $(LIBS) -lf2c -o meteochem
	$(STRIP) meteochem

meteochem : $(OBJS) $(LIBS)
	$(CC) $(OBJS) -L/usr/local/lib -lm -lnetcdf $(LIBS) -lf2c -o meteochem
	$(STRIP) meteochem

mcinp.c : mcinp.l
	flex -P mci -+ -o $@ $<

mcsyntax.c mcsyntax.h : mcsyntax.y
	bison -y -d mcsyntax.y
	sed 's/yy/mci/g' y.tab.c > mcsyntax.c
	rm y.tab.c
	sed 's/yy/mci/g' y.tab.h > mcsyntax.h
	rm y.tab.h

mcsyntax.incl : mcsyntax.y
	sed '1,/%%/d;/{/,/}$$/d;/%%/,$$d;/^$$/d;s/\;$$//;s/^/\"/;s/$$/\\n\"/' mcsyntax.y > mcsyntax.incl

mchemsyntax.c mchemsyntax.h : mchemsyntax.y
	bison -y -d mchemsyntax.y
	sed 's/yy/chi/g' y.tab.c > mchemsyntax.c
	rm y.tab.c
	sed 's/yy/chi/g' y.tab.h > mchemsyntax.h
	rm y.tab.h

mcheminp.c : mcheminp.l
	flex -P chi -+ -o $@ $<

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(OBJDIR)/mcmain.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcmain.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcmain.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcmain.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcmain.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcmain.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcmain.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcmain.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcmain.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcmain.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcmain.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcmain.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcmain.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcmain.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcmain.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcmain.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcmain.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcmain.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcmain.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcmain.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcmain.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcmain.o: /usr/include/bits/mathcalls.h mcglobal.h mcpress.h
$(OBJDIR)/mcmain.o: mcdynamics.h mcprint.h mcclouds.h mcttt.h mckeps.h
$(OBJDIR)/mcmain.o: mcfloat.h mcfilter.h mchelp.h mcemiss.h mcadvect.h
$(OBJDIR)/mcmain.o: sunpos.h mcground.h mcparse.h mc_group.hh mccdfin.h
$(OBJDIR)/mcmain.o: mcdep.h mchemparse.h mchem.h mcnest.h mcppm.h
$(OBJDIR)/mcmain.o: mc_module.hh mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mcglobal.o: /usr/include/stddef.h /usr/include/stdlib.h
$(OBJDIR)/mcglobal.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcglobal.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/types.h /usr/include/bits/types.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/typesizes.h /usr/include/time.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcglobal.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcglobal.o: /usr/include/stdio.h /usr/include/libio.h
$(OBJDIR)/mcglobal.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
$(OBJDIR)/mcglobal.o: /usr/include/xlocale.h /usr/include/malloc.h
$(OBJDIR)/mcglobal.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcglobal.o: /usr/include/bits/mathcalls.h mcglobal.h sunpos.h
$(OBJDIR)/mcglobal.o: mcground.h mcparse.h mc_group.hh mcprint.h
$(OBJDIR)/mcdynamics.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcdynamics.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcdynamics.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcdynamics.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/types.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mcdynamics.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
$(OBJDIR)/mcdynamics.o: /usr/include/xlocale.h /usr/include/math.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/huge_val.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcdynamics.o: /usr/include/bits/mathcalls.h /usr/include/malloc.h
$(OBJDIR)/mcdynamics.o: mcglobal.h sunpos.h mcground.h mcpress.h mcdynamics.h
$(OBJDIR)/mcdynamics.o: mcfilter.h
$(OBJDIR)/mcparse.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcparse.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcparse.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcparse.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcparse.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcparse.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcparse.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcparse.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcparse.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcparse.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcparse.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcparse.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcparse.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcparse.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcparse.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcparse.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcparse.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcparse.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcparse.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcparse.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcparse.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcparse.o: /usr/include/bits/mathcalls.h mcglobal.h sunpos.h
$(OBJDIR)/mcparse.o: mcground.h mcgroundclasses.h mcpress.h mcparse.h
$(OBJDIR)/mcparse.o: mc_group.hh mcemiss.h mcfilter.h mcdep.h mchemparse.h
$(OBJDIR)/mcparse.o: mchem.h mcclouds.h mcnest.h mpdata.h mcppm.h mcrely.h
$(OBJDIR)/mcparse.o: mc_module.hh mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mcparse.o: mc_group.t
$(OBJDIR)/mcpress.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcpress.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcpress.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcpress.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcpress.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcpress.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcpress.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcpress.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcpress.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcpress.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcpress.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcpress.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcpress.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcpress.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcpress.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcpress.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcpress.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcpress.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcpress.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcpress.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcpress.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcpress.o: /usr/include/bits/mathcalls.h sunpos.h mcglobal.h
$(OBJDIR)/mcpress.o: mcground.h mcpress.h
$(OBJDIR)/mcdep.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcdep.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcdep.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcdep.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcdep.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcdep.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcdep.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcdep.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcdep.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcdep.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcdep.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcdep.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcdep.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcdep.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcdep.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcdep.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcdep.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcdep.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcdep.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcdep.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcdep.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcdep.o: /usr/include/bits/mathcalls.h mcglobal.h sunpos.h
$(OBJDIR)/mcdep.o: mcground.h mcparse.h mc_group.hh mcdep.h
$(OBJDIR)/sunpos.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/sunpos.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/sunpos.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/sunpos.o: /usr/include/gnu/stubs-64.h /usr/include/stddef.h
$(OBJDIR)/sunpos.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/sunpos.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/sunpos.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/sunpos.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
$(OBJDIR)/sunpos.o: /usr/include/bits/huge_val.h
$(OBJDIR)/sunpos.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/sunpos.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/sunpos.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/sunpos.o: /usr/include/bits/mathcalls.h sunpos.h
$(OBJDIR)/mcground.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcground.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcground.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcground.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcground.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcground.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcground.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcground.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcground.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcground.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcground.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcground.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcground.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcground.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcground.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcground.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcground.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcground.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcground.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcground.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcground.o: /usr/include/bits/mathcalls.h mcglobal.h mcpress.h
$(OBJDIR)/mcground.o: sunpos.h mcground.h mcfloat.h mctwostream/mctwostream.h
$(OBJDIR)/mcground.o: mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mcground.o: mc_commands.hh mctwostream/mcphotodiss.h
$(OBJDIR)/mcground.o: mctwostream/mcapd.h mc_group.hh
$(OBJDIR)/mcfloat.o: /usr/include/signal.h /usr/include/features.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/mcfloat.o: /usr/include/gnu/stubs-64.h /usr/include/bits/sigset.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/signum.h /usr/include/time.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/siginfo.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/sigaction.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/sigcontext.h /usr/include/stddef.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/sigstack.h /usr/include/sys/ucontext.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/sigthread.h /usr/include/stdio.h
$(OBJDIR)/mcfloat.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcfloat.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/types.h /usr/include/sys/select.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/select.h /usr/include/bits/time.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/sysmacros.h /usr/include/alloca.h
$(OBJDIR)/mcfloat.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcfloat.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/environments.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/confname.h /usr/include/getopt.h
$(OBJDIR)/mcfloat.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcfloat.o: /usr/include/bits/resource.h mcglobal.h mcfloat.h
$(OBJDIR)/mcfloat.o: mcprint.h
$(OBJDIR)/mcttt.o: /usr/include/stddef.h /usr/include/stdlib.h
$(OBJDIR)/mcttt.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcttt.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcttt.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcttt.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcttt.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcttt.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcttt.o: /usr/include/sys/types.h /usr/include/bits/types.h
$(OBJDIR)/mcttt.o: /usr/include/bits/typesizes.h /usr/include/time.h
$(OBJDIR)/mcttt.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcttt.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcttt.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcttt.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcttt.o: /usr/include/stdio.h /usr/include/libio.h
$(OBJDIR)/mcttt.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcttt.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcttt.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
$(OBJDIR)/mcttt.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
$(OBJDIR)/mcttt.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcttt.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcttt.o: /usr/include/bits/mathcalls.h /usr/include/malloc.h
$(OBJDIR)/mcttt.o: mcglobal.h mcpress.h sunpos.h mcground.h mcttt.h
$(OBJDIR)/mcfilter.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcfilter.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcfilter.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcfilter.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcfilter.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcfilter.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcfilter.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcfilter.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcfilter.o: /usr/include/bits/mathcalls.h mcglobal.h mcfilter.h
$(OBJDIR)/mcemiss.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcemiss.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcemiss.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcemiss.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcemiss.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcemiss.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcemiss.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcemiss.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcemiss.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcemiss.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcemiss.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcemiss.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcemiss.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcemiss.o: /usr/include/string.h /usr/include/xlocale.h mcglobal.h
$(OBJDIR)/mcemiss.o: sunpos.h mcground.h mcparse.h mc_group.hh mcemiss.h
$(OBJDIR)/mcsyntax.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/mcsyntax.o: /usr/include/gnu/stubs-64.h /usr/include/stddef.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcsyntax.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcsyntax.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
$(OBJDIR)/mcsyntax.o: /usr/include/xlocale.h /usr/include/stdlib.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcsyntax.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcsyntax.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcsyntax.o: mcglobal.h sunpos.h mcground.h mcgroundclasses.h
$(OBJDIR)/mcsyntax.o: mcparse.h mc_group.hh mccdfin.h mcprint.h mchemparse.h
$(OBJDIR)/mcsyntax.o: mcnest.h mc_commands.hh mc_group.t
$(OBJDIR)/mcinp.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcinp.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mcinp.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/mcinp.o: /usr/include/gnu/stubs-64.h /usr/include/stddef.h
$(OBJDIR)/mcinp.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcinp.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcinp.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcinp.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
$(OBJDIR)/mcinp.o: /usr/include/xlocale.h /usr/include/stdlib.h
$(OBJDIR)/mcinp.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcinp.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcinp.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcinp.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcinp.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcinp.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcinp.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcinp.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcinp.o: mcglobal.h sunpos.h mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mcinp.o: mcprint.h mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mcinp.o: mc_commands.hh mcsyntax.h mc_group.t
$(OBJDIR)/mccdfin.o: /usr/include/stddef.h /usr/include/stdlib.h
$(OBJDIR)/mccdfin.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mccdfin.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/types.h /usr/include/bits/types.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/typesizes.h /usr/include/time.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mccdfin.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mccdfin.o: /usr/include/stdio.h /usr/include/libio.h
$(OBJDIR)/mccdfin.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mccdfin.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
$(OBJDIR)/mccdfin.o: /usr/include/xlocale.h /usr/include/netcdf.h
$(OBJDIR)/mccdfin.o: /usr/include/errno.h /usr/include/bits/errno.h
$(OBJDIR)/mccdfin.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
$(OBJDIR)/mccdfin.o: /usr/include/asm-generic/errno.h
$(OBJDIR)/mccdfin.o: /usr/include/asm-generic/errno-base.h mcglobal.h
$(OBJDIR)/mccdfin.o: sunpos.h mcground.h mcparse.h mc_group.hh mccdfin.h
$(OBJDIR)/mccdfin.o: mcemiss.h
$(OBJDIR)/mcprint.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcprint.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcprint.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcprint.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcprint.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcprint.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcprint.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcprint.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcprint.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcprint.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcprint.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcprint.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcprint.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcprint.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcprint.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcprint.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcprint.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcprint.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcprint.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcprint.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcprint.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcprint.o: /usr/include/bits/mathcalls.h /usr/include/netcdf.h
$(OBJDIR)/mcprint.o: /usr/include/errno.h /usr/include/bits/errno.h
$(OBJDIR)/mcprint.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
$(OBJDIR)/mcprint.o: /usr/include/asm-generic/errno.h
$(OBJDIR)/mcprint.o: /usr/include/asm-generic/errno-base.h mcglobal.h
$(OBJDIR)/mcprint.o: sunpos.h mcground.h mcparse.h mc_group.hh mcprint.h
$(OBJDIR)/mcprint.o: mchemparse.h mchem.h
$(OBJDIR)/mchelp.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mchelp.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mchelp.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mchelp.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mchelp.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mchelp.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mchelp.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mchelp.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mchelp.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mchelp.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mchelp.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mchelp.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mchelp.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mchelp.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mchelp.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mchelp.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mchelp.o: /usr/include/string.h /usr/include/xlocale.h mcglobal.h
$(OBJDIR)/mchelp.o: sunpos.h mcground.h mcparse.h mc_group.hh mchelp.h
$(OBJDIR)/mchelp.o: mc_module.hh mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mchelp.o: mctwostream/mctwostream.h mctwostream/mcphotodiss.h
$(OBJDIR)/mchelp.o: mctwostream/mcapd.h mcsyntax.incl
$(OBJDIR)/mcadvect.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcadvect.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcadvect.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcadvect.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcadvect.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcadvect.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcadvect.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcadvect.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcadvect.o: /usr/include/bits/mathcalls.h mcglobal.h sunpos.h
$(OBJDIR)/mcadvect.o: mcground.h mchemparse.h mpdata.h
$(OBJDIR)/mpdata.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mpdata.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mpdata.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mpdata.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mpdata.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mpdata.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mpdata.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mpdata.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mpdata.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mpdata.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mpdata.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mpdata.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mpdata.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mpdata.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mpdata.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mpdata.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mpdata.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mpdata.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mpdata.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mpdata.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mpdata.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mpdata.o: /usr/include/bits/mathcalls.h mcglobal.h mpdata.h
$(OBJDIR)/mchemsyntax.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mchemsyntax.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/cdefs.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/wordsize.h
$(OBJDIR)/mchemsyntax.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/types.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mchemsyntax.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/sys_errlist.h
$(OBJDIR)/mchemsyntax.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/endian.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/byteswap.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mchemsyntax.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mchemsyntax.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mchemsyntax.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mchemsyntax.o: /usr/include/xlocale.h mcglobal.h mchemparse.h
$(OBJDIR)/mchemsyntax.o: mctwostream/mctwostream.h mc_module.hh
$(OBJDIR)/mchemsyntax.o: mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mchemsyntax.o: mctwostream/mcphotodiss.h mctwostream/mcapd.h
$(OBJDIR)/mchemsyntax.o: mc_group.hh
$(OBJDIR)/mcheminp.o: /usr/include/errno.h /usr/include/features.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/mcheminp.o: /usr/include/gnu/stubs-64.h /usr/include/bits/errno.h
$(OBJDIR)/mcheminp.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
$(OBJDIR)/mcheminp.o: /usr/include/asm-generic/errno.h
$(OBJDIR)/mcheminp.o: /usr/include/asm-generic/errno-base.h
$(OBJDIR)/mcheminp.o: /usr/include/FlexLexer.h /usr/include/stdio.h
$(OBJDIR)/mcheminp.o: /usr/include/stddef.h /usr/include/bits/types.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mcheminp.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
$(OBJDIR)/mcheminp.o: /usr/include/xlocale.h /usr/include/stdlib.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcheminp.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcheminp.o: mcglobal.h mchemsyntax.h mchemparse.h
$(OBJDIR)/mcheminp.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/environments.h
$(OBJDIR)/mcheminp.o: /usr/include/bits/confname.h /usr/include/getopt.h
$(OBJDIR)/mchem.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mchem.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mchem.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mchem.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mchem.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mchem.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mchem.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mchem.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mchem.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mchem.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mchem.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mchem.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mchem.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mchem.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mchem.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mchem.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mchem.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mchem.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mchem.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mchem.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mchem.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mchem.o: /usr/include/bits/mathcalls.h mcglobal.h sunpos.h
$(OBJDIR)/mchem.o: mcground.h mchemparse.h mchem.h matrix.h mcparallel.h
$(OBJDIR)/mchemparse.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mchemparse.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mchemparse.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mchemparse.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/types.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mchemparse.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/endian.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
$(OBJDIR)/mchemparse.o: /usr/include/time.h /usr/include/sys/select.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mchemparse.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mchemparse.o: /usr/include/xlocale.h /usr/include/math.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/huge_val.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mchemparse.o: /usr/include/bits/mathcalls.h matrix.h mcglobal.h
$(OBJDIR)/mchemparse.o: sunpos.h mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mchemparse.o: mchemparse.h
$(OBJDIR)/matrix.o: /usr/include/stddef.h /usr/include/stdlib.h
$(OBJDIR)/matrix.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/matrix.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/matrix.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/matrix.o: /usr/include/bits/waitflags.h
$(OBJDIR)/matrix.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/matrix.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/matrix.o: /usr/include/sys/types.h /usr/include/bits/types.h
$(OBJDIR)/matrix.o: /usr/include/bits/typesizes.h /usr/include/time.h
$(OBJDIR)/matrix.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/matrix.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/matrix.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/matrix.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/matrix.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/matrix.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/matrix.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/matrix.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/matrix.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/matrix.o: /usr/include/bits/mathcalls.h matrix.h
$(OBJDIR)/mcclouds.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcclouds.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcclouds.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcclouds.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcclouds.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcclouds.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcclouds.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcclouds.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcclouds.o: /usr/include/bits/mathcalls.h sunpos.h mcglobal.h
$(OBJDIR)/mcclouds.o: mcground.h mcpress.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/features.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/predefs.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/cdefs.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/wordsize.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/gnu/stubs.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/types.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/typesizes.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/wchar.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/sys_errlist.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stdlib.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/waitstatus.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/endian.h /usr/include/bits/endian.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/byteswap.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/select.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/select.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/sigset.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/time.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/xlocale.h mcglobal.h sunpos.h
$(OBJDIR)/mcgroundclasses.o: mcground.h mcparse.h mc_group.hh
$(OBJDIR)/mcgroundclasses.o: mcgroundclasses.h
$(OBJDIR)/mcspecial.o: /usr/include/stddef.h /usr/include/stdlib.h
$(OBJDIR)/mcspecial.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcspecial.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/endian.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/types.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/typesizes.h /usr/include/time.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcspecial.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcspecial.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcspecial.o: /usr/include/bits/mathcalls.h
$(OBJDIR)/mckeps.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mckeps.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mckeps.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mckeps.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mckeps.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mckeps.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mckeps.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mckeps.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mckeps.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mckeps.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mckeps.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mckeps.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mckeps.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mckeps.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mckeps.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mckeps.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mckeps.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mckeps.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mckeps.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mckeps.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mckeps.o: /usr/include/bits/mathcalls.h mcglobal.h sunpos.h
$(OBJDIR)/mckeps.o: mcground.h mchemparse.h
$(OBJDIR)/mcnest.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcnest.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcnest.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcnest.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcnest.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcnest.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcnest.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcnest.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcnest.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcnest.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcnest.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcnest.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcnest.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcnest.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcnest.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcnest.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcnest.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcnest.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcnest.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcnest.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcnest.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcnest.o: /usr/include/bits/mathcalls.h /usr/include/netcdf.h
$(OBJDIR)/mcnest.o: /usr/include/errno.h /usr/include/bits/errno.h
$(OBJDIR)/mcnest.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
$(OBJDIR)/mcnest.o: /usr/include/asm-generic/errno.h
$(OBJDIR)/mcnest.o: /usr/include/asm-generic/errno-base.h mcglobal.h sunpos.h
$(OBJDIR)/mcnest.o: mcground.h mcparse.h mc_group.hh mcnest.h mcprint.h
$(OBJDIR)/mcppm.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcppm.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcppm.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcppm.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcppm.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcppm.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcppm.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcppm.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcppm.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcppm.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcppm.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcppm.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcppm.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcppm.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcppm.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcppm.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcppm.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcppm.o: /usr/include/math.h /usr/include/bits/huge_val.h
$(OBJDIR)/mcppm.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcppm.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
$(OBJDIR)/mcppm.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
$(OBJDIR)/mcppm.o: /usr/include/bits/mathcalls.h mcglobal.h mcppm.h sunpos.h
$(OBJDIR)/mcppm.o: mcground.h mchemparse.h
$(OBJDIR)/mc_module.o: mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mc_module.o: mc_commands.hh /usr/include/stddef.h
$(OBJDIR)/mc_module.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mc_module.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mc_module.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/mc_module.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
$(OBJDIR)/mc_module.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mc_module.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mc_module.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mc_module.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mc_module.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mc_module.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mc_module.o: /usr/include/bits/endian.h
$(OBJDIR)/mc_module.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
$(OBJDIR)/mc_module.o: /usr/include/time.h /usr/include/sys/select.h
$(OBJDIR)/mc_module.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
$(OBJDIR)/mc_module.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
$(OBJDIR)/mc_module.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mc_module.o: mcparse.h mc_group.hh mc_group.t
$(OBJDIR)/mc_variable.o: mc_variable.hh /usr/include/stddef.h
$(OBJDIR)/mc_variable.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/wordsize.h
$(OBJDIR)/mc_variable.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/types.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mc_variable.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/sys_errlist.h
$(OBJDIR)/mc_variable.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/endian.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/byteswap.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mc_variable.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mc_variable.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mc_variable.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mc_variable.o: /usr/include/xlocale.h mcground.h mc_group.hh
$(OBJDIR)/mc_variable.o: mcparse.h mc_group.t
$(OBJDIR)/mc_section.o: mc_section.hh /usr/include/stddef.h
$(OBJDIR)/mc_section.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mc_section.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mc_section.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/mc_section.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
$(OBJDIR)/mc_section.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mc_section.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mc_section.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mc_section.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
$(OBJDIR)/mc_section.o: /usr/include/xlocale.h mcglobal.h mcrely.h mcnest.h
$(OBJDIR)/mc_commands.o: mc_commands.hh /usr/include/stddef.h
$(OBJDIR)/mc_commands.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/wordsize.h
$(OBJDIR)/mc_commands.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/types.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mc_commands.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/sys_errlist.h
$(OBJDIR)/mc_commands.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/endian.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/byteswap.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mc_commands.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mc_commands.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mc_commands.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mc_commands.o: /usr/include/xlocale.h mc_group.t mc_group.hh
$(OBJDIR)/mc_commands.o: mcparse.h
$(OBJDIR)/puffemit.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/puffemit.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/puffemit.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/puffemit.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/puffemit.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/puffemit.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/puffemit.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/puffemit.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/puffemit.o: /usr/include/bits/waitflags.h
$(OBJDIR)/puffemit.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/puffemit.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/puffemit.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/puffemit.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/puffemit.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/puffemit.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/puffemit.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/puffemit.o: /usr/include/string.h /usr/include/xlocale.h mcglobal.h
$(OBJDIR)/puffemit.o: mcground.h mcparse.h mc_group.hh mc_section.hh
$(OBJDIR)/puffemit.o: mc_group.t mc_commands.hh mc_variable.hh mc_module.hh
