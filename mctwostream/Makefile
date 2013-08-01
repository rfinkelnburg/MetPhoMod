OBJDIR?= .
OBJS= $(OBJDIR)/mctwostream.o $(OBJDIR)/mcphotodiss.o \
      $(OBJDIR)/mcabsorb.o $(OBJDIR)/mcextrad.o\
      $(OBJDIR)/mc_inter3.o $(OBJDIR)/mc_rtlink.o $(OBJDIR)/mc_schu.o \
      $(OBJDIR)/mc_seto2.o $(OBJDIR)/mc_sphers.o $(OBJDIR)/mc_ps2str.o \
      $(OBJDIR)/mc_fchap.o $(OBJDIR)/mcapd.o

SRC= mcphotodiss.c mctwostream.c mcabsorb.c mcextrad.c mcapd.c

TARGET= libmctwostream.a

CC=g++
F77=f77

CFLAGS2=-I..

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(GENFLAGS) $(INCLUDE) -DFERRET_FORMAT $*.c -o $@

$(OBJDIR)/%.o : %.cc
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(GENFLAGS) $(INCLUDE) -DFERRET_FORMAT $*.cc -o $@

$(OBJDIR)/%.o : %.f
	$(F77) -c $(FFLAGS) $(GENFLAGS) $*.f -o $@

$(TARGET) : $(OBJS)
	ar cr $(TARGET) $(OBJS)

clean :
	rm -f $(OBJS) $(TARGET)

depend :
	makedepend -I.. -p'$$(OBJDIR)/' $(SRC)

$(OBJDIR)/mc_inter3.o : mc_inter3.f params

$(OBJDIR)/mc_rtlink.o : mc_rtlink.f params

$(OBJDIR)/mc_schu.o : mc_schu.f params

$(OBJDIR)/mc_seto2.o : mc_schu.f params

$(OBJDIR)/mc_sphers.o : mc_schu.f params

$(OBJDIR)/mc_fchap.o : mc_fchap.f params

$(OBJDIR)/mc_ps2str.o : mc_ps2str.f params

# DO NOT DELETE

$(OBJDIR)/mcphotodiss.o: mcphotodiss.h mcapd.h /usr/include/stddef.h
$(OBJDIR)/mcphotodiss.o: ../mc_group.hh /usr/include/stdio.h
$(OBJDIR)/mcphotodiss.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/cdefs.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/wordsize.h
$(OBJDIR)/mcphotodiss.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/types.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mcphotodiss.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/sys_errlist.h
$(OBJDIR)/mcphotodiss.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/endian.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/byteswap.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mcphotodiss.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mcphotodiss.o: /usr/include/xlocale.h /usr/include/math.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/huge_val.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/huge_vall.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/mathdef.h
$(OBJDIR)/mcphotodiss.o: /usr/include/bits/mathcalls.h ../mc_group.t
$(OBJDIR)/mcphotodiss.o: mcextrad.h ../mcglobal.h ph_o3data.c ph2data.c
$(OBJDIR)/mcphotodiss.o: ph3data.c ph4data.c ph6data.c ph7data.c ph11data.c
$(OBJDIR)/mcphotodiss.o: ph12data.c ph13data.c ph14data.c
$(OBJDIR)/mctwostream.o: mctwostream.h ../mcground.h ../mcglobal.h
$(OBJDIR)/mctwostream.o: ../sunpos.h ../mc_module.hh ../mc_variable.hh
$(OBJDIR)/mctwostream.o: ../mcparse.h ../mc_group.hh ../mc_section.hh
$(OBJDIR)/mctwostream.o: ../mc_commands.hh mcphotodiss.h mcapd.h
$(OBJDIR)/mctwostream.o: /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mctwostream.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/cdefs.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/wordsize.h
$(OBJDIR)/mctwostream.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/types.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mctwostream.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/sys_errlist.h
$(OBJDIR)/mctwostream.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/endian.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/byteswap.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/pthreadtypes.h
$(OBJDIR)/mctwostream.o: /usr/include/alloca.h /usr/include/string.h
$(OBJDIR)/mctwostream.o: /usr/include/xlocale.h /usr/include/math.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/huge_val.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/huge_valf.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/huge_vall.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/mathdef.h
$(OBJDIR)/mctwostream.o: /usr/include/bits/mathcalls.h ../mchemparse.h
$(OBJDIR)/mctwostream.o: ../mchem.h mcabsorb.h mcextrad.h
$(OBJDIR)/mcabsorb.o: mcabsorb.h /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcabsorb.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcabsorb.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcabsorb.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcabsorb.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcabsorb.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcabsorb.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcabsorb.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcabsorb.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcabsorb.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcabsorb.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcabsorb.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcabsorb.o: ../mcparse.h ../mc_group.hh mcextrad.h ph_o3data.c
$(OBJDIR)/mcextrad.o: mcextrad.h /usr/include/stddef.h /usr/include/stdio.h
$(OBJDIR)/mcextrad.o: /usr/include/features.h /usr/include/bits/predefs.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
$(OBJDIR)/mcextrad.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
$(OBJDIR)/mcextrad.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
$(OBJDIR)/mcextrad.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcextrad.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcextrad.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcextrad.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcextrad.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcextrad.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcextrad.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcextrad.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcextrad.o: ../mcparse.h ../mc_group.hh
$(OBJDIR)/mcapd.o: mcapd.h /usr/include/stddef.h ../mc_group.hh
$(OBJDIR)/mcapd.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcapd.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
$(OBJDIR)/mcapd.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
$(OBJDIR)/mcapd.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
$(OBJDIR)/mcapd.o: /usr/include/bits/typesizes.h /usr/include/libio.h
$(OBJDIR)/mcapd.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcapd.o: /usr/include/bits/stdio_lim.h
$(OBJDIR)/mcapd.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
$(OBJDIR)/mcapd.o: /usr/include/bits/waitflags.h
$(OBJDIR)/mcapd.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
$(OBJDIR)/mcapd.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
$(OBJDIR)/mcapd.o: /usr/include/sys/types.h /usr/include/time.h
$(OBJDIR)/mcapd.o: /usr/include/sys/select.h /usr/include/bits/select.h
$(OBJDIR)/mcapd.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
$(OBJDIR)/mcapd.o: /usr/include/sys/sysmacros.h
$(OBJDIR)/mcapd.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
$(OBJDIR)/mcapd.o: ../mc_group.t
