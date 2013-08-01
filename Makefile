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
CFLAGS2=-DSYSTEM=\"`uname`\" -DDATE=\"`date +%D`\" -DLINUX -fpermissive

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
	bison -p chi -d mchemsyntax.y
	mv mchemsyntax.tab.c mchemsyntax.c
	mv mchemsyntax.tab.h mchemsyntax.h

mcheminp.c mcheminp.h : mcheminp.l
	flex -P chi --header-file=$*.h -o $*.c $<

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(OBJDIR)/mcmain.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcmain.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcmain.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcmain.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcmain.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcmain.o: /usr/include/math.h mcglobal.h mcpress.h mcdynamics.h
$(OBJDIR)/mcmain.o: mcprint.h mcclouds.h mcttt.h mckeps.h mcfloat.h
$(OBJDIR)/mcmain.o: mcfilter.h mchelp.h mcemiss.h mcadvect.h sunpos.h
$(OBJDIR)/mcmain.o: mcground.h mcparse.h mc_group.hh mccdfin.h mcdep.h
$(OBJDIR)/mcmain.o: mchemparse.h mchem.h mcnest.h mcppm.h mc_module.hh
$(OBJDIR)/mcmain.o: mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mcglobal.o: /usr/include/stdlib.h /usr/include/features.h
$(OBJDIR)/mcglobal.o: /usr/include/stdc-predef.h /usr/include/alloca.h
$(OBJDIR)/mcglobal.o: /usr/include/stdio.h /usr/include/libio.h
$(OBJDIR)/mcglobal.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcglobal.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcglobal.o: /usr/include/malloc.h /usr/include/math.h mcglobal.h
$(OBJDIR)/mcglobal.o: sunpos.h mcground.h mcparse.h mc_group.hh mcprint.h
$(OBJDIR)/mcdynamics.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcdynamics.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcdynamics.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcdynamics.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcdynamics.o: /usr/include/math.h /usr/include/malloc.h mcglobal.h
$(OBJDIR)/mcdynamics.o: sunpos.h mcground.h mcpress.h mcdynamics.h mcfilter.h
$(OBJDIR)/mcparse.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcparse.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcparse.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcparse.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcparse.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcparse.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcparse.o: mcgroundclasses.h mcpress.h mcparse.h mc_group.hh
$(OBJDIR)/mcparse.o: mcemiss.h mcfilter.h mcdep.h mchemparse.h mchem.h
$(OBJDIR)/mcparse.o: mcclouds.h mcnest.h mpdata.h mcppm.h mcrely.h
$(OBJDIR)/mcparse.o: mc_module.hh mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mcparse.o: mc_group.t
$(OBJDIR)/mcpress.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcpress.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcpress.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcpress.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcpress.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcpress.o: /usr/include/math.h sunpos.h mcglobal.h mcground.h
$(OBJDIR)/mcpress.o: mcpress.h
$(OBJDIR)/mcdep.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcdep.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcdep.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcdep.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcdep.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcdep.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcdep.o: mcparse.h mc_group.hh mcdep.h
$(OBJDIR)/sunpos.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/sunpos.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/sunpos.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/sunpos.o: /usr/include/math.h sunpos.h
$(OBJDIR)/mcground.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcground.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcground.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcground.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcground.o: /usr/include/math.h mcglobal.h mcpress.h sunpos.h
$(OBJDIR)/mcground.o: mcground.h mcfloat.h mctwostream/mctwostream.h
$(OBJDIR)/mcground.o: mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mcground.o: mc_commands.hh mctwostream/mcphotodiss.h
$(OBJDIR)/mcground.o: mctwostream/mcapd.h mc_group.hh
$(OBJDIR)/mcfloat.o: /usr/include/signal.h /usr/include/features.h
$(OBJDIR)/mcfloat.o: /usr/include/stdc-predef.h /usr/include/time.h
$(OBJDIR)/mcfloat.o: /usr/include/stdio.h /usr/include/libio.h
$(OBJDIR)/mcfloat.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcfloat.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcfloat.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcfloat.o: /usr/include/unistd.h /usr/include/getopt.h mcglobal.h
$(OBJDIR)/mcfloat.o: mcfloat.h mcprint.h
$(OBJDIR)/mcttt.o: /usr/include/stdlib.h /usr/include/features.h
$(OBJDIR)/mcttt.o: /usr/include/stdc-predef.h /usr/include/alloca.h
$(OBJDIR)/mcttt.o: /usr/include/stdio.h /usr/include/libio.h
$(OBJDIR)/mcttt.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcttt.o: /usr/include/math.h /usr/include/malloc.h mcglobal.h
$(OBJDIR)/mcttt.o: mcpress.h sunpos.h mcground.h mcttt.h
$(OBJDIR)/mcfilter.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcfilter.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcfilter.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcfilter.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcfilter.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcfilter.o: /usr/include/math.h mcglobal.h mcfilter.h
$(OBJDIR)/mcemiss.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcemiss.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcemiss.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcemiss.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcemiss.o: /usr/include/string.h /usr/include/xlocale.h mcglobal.h
$(OBJDIR)/mcemiss.o: sunpos.h mcground.h mcparse.h mc_group.hh mcemiss.h
$(OBJDIR)/mcinp.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcinp.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcinp.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcinp.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcinp.o: /usr/include/stdlib.h /usr/include/alloca.h mcglobal.h
$(OBJDIR)/mcinp.o: sunpos.h mcground.h mcparse.h mc_group.hh mcprint.h
$(OBJDIR)/mcinp.o: mc_module.hh mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mcinp.o: mc_group.t
$(OBJDIR)/mccdfin.o: /usr/include/stdlib.h /usr/include/features.h
$(OBJDIR)/mccdfin.o: /usr/include/stdc-predef.h /usr/include/alloca.h
$(OBJDIR)/mccdfin.o: /usr/include/stdio.h /usr/include/libio.h
$(OBJDIR)/mccdfin.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mccdfin.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mccdfin.o: /usr/include/netcdf.h /usr/include/errno.h mcglobal.h
$(OBJDIR)/mccdfin.o: sunpos.h mcground.h mcparse.h mc_group.hh mccdfin.h
$(OBJDIR)/mccdfin.o: mcemiss.h
$(OBJDIR)/mcprint.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcprint.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcprint.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcprint.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcprint.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcprint.o: /usr/include/math.h /usr/include/netcdf.h
$(OBJDIR)/mcprint.o: /usr/include/errno.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcprint.o: mcparse.h mc_group.hh mcprint.h mchemparse.h mchem.h
$(OBJDIR)/mchelp.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mchelp.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mchelp.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mchelp.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mchelp.o: /usr/include/string.h /usr/include/xlocale.h mcglobal.h
$(OBJDIR)/mchelp.o: sunpos.h mcground.h mcparse.h mc_group.hh mchelp.h
$(OBJDIR)/mchelp.o: mc_module.hh mc_variable.hh mc_section.hh mc_commands.hh
$(OBJDIR)/mchelp.o: mctwostream/mctwostream.h mctwostream/mcphotodiss.h
$(OBJDIR)/mchelp.o: mctwostream/mcapd.h
$(OBJDIR)/mcadvect.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcadvect.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcadvect.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcadvect.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcadvect.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcadvect.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcadvect.o: mchemparse.h mpdata.h
$(OBJDIR)/mpdata.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mpdata.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mpdata.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mpdata.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mpdata.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mpdata.o: /usr/include/math.h mcglobal.h mpdata.h
$(OBJDIR)/mchem.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mchem.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mchem.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mchem.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mchem.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mchem.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mchem.o: mchemparse.h mchem.h matrix.h mcparallel.h
$(OBJDIR)/mchemparse.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mchemparse.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mchemparse.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mchemparse.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mchemparse.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mchemparse.o: /usr/include/math.h matrix.h mcglobal.h sunpos.h
$(OBJDIR)/mchemparse.o: mcground.h mcparse.h mc_group.hh mchemparse.h
$(OBJDIR)/matrix.o: /usr/include/stdlib.h /usr/include/features.h
$(OBJDIR)/matrix.o: /usr/include/stdc-predef.h /usr/include/alloca.h
$(OBJDIR)/matrix.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/matrix.o: /usr/include/math.h matrix.h
$(OBJDIR)/mcclouds.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcclouds.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcclouds.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcclouds.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcclouds.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcclouds.o: /usr/include/math.h sunpos.h mcglobal.h mcground.h
$(OBJDIR)/mcclouds.o: mcpress.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcgroundclasses.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcgroundclasses.o: mcglobal.h sunpos.h mcground.h mcparse.h
$(OBJDIR)/mcgroundclasses.o: mc_group.hh mcgroundclasses.h
$(OBJDIR)/mcspecial.o: /usr/include/stdlib.h /usr/include/features.h
$(OBJDIR)/mcspecial.o: /usr/include/stdc-predef.h /usr/include/alloca.h
$(OBJDIR)/mcspecial.o: /usr/include/math.h
$(OBJDIR)/mckeps.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mckeps.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mckeps.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mckeps.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mckeps.o: /usr/include/math.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mckeps.o: mchemparse.h
$(OBJDIR)/mcnest.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcnest.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcnest.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcnest.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcnest.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcnest.o: /usr/include/math.h /usr/include/netcdf.h
$(OBJDIR)/mcnest.o: /usr/include/errno.h mcglobal.h sunpos.h mcground.h
$(OBJDIR)/mcnest.o: mcparse.h mc_group.hh mcnest.h mcprint.h
$(OBJDIR)/mcppm.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcppm.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcppm.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcppm.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcppm.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcppm.o: /usr/include/math.h mcglobal.h mcppm.h sunpos.h mcground.h
$(OBJDIR)/mcppm.o: mchemparse.h
$(OBJDIR)/mc_module.o: mc_module.hh mc_variable.hh mc_section.hh
$(OBJDIR)/mc_module.o: mc_commands.hh /usr/include/stdio.h
$(OBJDIR)/mc_module.o: /usr/include/features.h /usr/include/stdc-predef.h
$(OBJDIR)/mc_module.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mc_module.o: /usr/include/wchar.h /usr/include/stdlib.h
$(OBJDIR)/mc_module.o: /usr/include/alloca.h mcparse.h mc_group.hh mc_group.t
$(OBJDIR)/mc_variable.o: mc_variable.hh /usr/include/stdio.h
$(OBJDIR)/mc_variable.o: /usr/include/features.h /usr/include/stdc-predef.h
$(OBJDIR)/mc_variable.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mc_variable.o: /usr/include/wchar.h /usr/include/stdlib.h
$(OBJDIR)/mc_variable.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mc_variable.o: /usr/include/xlocale.h mcground.h mc_group.hh
$(OBJDIR)/mc_variable.o: mcparse.h mc_group.t
$(OBJDIR)/mc_section.o: mc_section.hh /usr/include/stdio.h
$(OBJDIR)/mc_section.o: /usr/include/features.h /usr/include/stdc-predef.h
$(OBJDIR)/mc_section.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mc_section.o: /usr/include/wchar.h /usr/include/string.h
$(OBJDIR)/mc_section.o: /usr/include/xlocale.h mcglobal.h mcrely.h mcnest.h
$(OBJDIR)/mc_commands.o: mc_commands.hh /usr/include/stdio.h
$(OBJDIR)/mc_commands.o: /usr/include/features.h /usr/include/stdc-predef.h
$(OBJDIR)/mc_commands.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mc_commands.o: /usr/include/wchar.h /usr/include/stdlib.h
$(OBJDIR)/mc_commands.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mc_commands.o: /usr/include/xlocale.h mc_group.t mcparse.h
$(OBJDIR)/mc_commands.o: mc_group.hh
$(OBJDIR)/puffemit.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/puffemit.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/puffemit.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/puffemit.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/puffemit.o: /usr/include/string.h /usr/include/xlocale.h mcglobal.h
$(OBJDIR)/puffemit.o: mcground.h mcparse.h mc_group.hh mc_section.hh
$(OBJDIR)/puffemit.o: mc_group.t mc_commands.hh mc_variable.hh mc_module.hh
