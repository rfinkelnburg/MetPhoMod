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

$(OBJDIR)/mcphotodiss.o: mcphotodiss.h mcapd.h ../mc_group.hh
$(OBJDIR)/mcphotodiss.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcphotodiss.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcphotodiss.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcphotodiss.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcphotodiss.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mcphotodiss.o: /usr/include/math.h ../mc_group.t mcextrad.h
$(OBJDIR)/mcphotodiss.o: ../mcglobal.h ph_o3data.c ph2data.c ph3data.c
$(OBJDIR)/mcphotodiss.o: ph4data.c ph6data.c ph7data.c ph11data.c ph12data.c
$(OBJDIR)/mcphotodiss.o: ph13data.c ph14data.c
$(OBJDIR)/mctwostream.o: mctwostream.h ../mcground.h ../mcglobal.h
$(OBJDIR)/mctwostream.o: ../sunpos.h ../mc_module.hh ../mc_variable.hh
$(OBJDIR)/mctwostream.o: ../mcparse.h ../mc_group.hh ../mc_section.hh
$(OBJDIR)/mctwostream.o: ../mc_commands.hh mcphotodiss.h mcapd.h
$(OBJDIR)/mctwostream.o: /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mctwostream.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mctwostream.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mctwostream.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mctwostream.o: /usr/include/string.h /usr/include/xlocale.h
$(OBJDIR)/mctwostream.o: /usr/include/math.h ../mchemparse.h ../mchem.h
$(OBJDIR)/mctwostream.o: mcabsorb.h mcextrad.h
$(OBJDIR)/mcabsorb.o: mcabsorb.h /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcabsorb.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcabsorb.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcabsorb.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcabsorb.o: ../mcparse.h ../mc_group.hh mcextrad.h ph_o3data.c
$(OBJDIR)/mcextrad.o: mcextrad.h /usr/include/stdio.h /usr/include/features.h
$(OBJDIR)/mcextrad.o: /usr/include/stdc-predef.h /usr/include/libio.h
$(OBJDIR)/mcextrad.o: /usr/include/_G_config.h /usr/include/wchar.h
$(OBJDIR)/mcextrad.o: /usr/include/stdlib.h /usr/include/alloca.h
$(OBJDIR)/mcextrad.o: ../mcparse.h ../mc_group.hh
$(OBJDIR)/mcapd.o: mcapd.h ../mc_group.hh /usr/include/stdio.h
$(OBJDIR)/mcapd.o: /usr/include/features.h /usr/include/stdc-predef.h
$(OBJDIR)/mcapd.o: /usr/include/libio.h /usr/include/_G_config.h
$(OBJDIR)/mcapd.o: /usr/include/wchar.h /usr/include/stdlib.h
$(OBJDIR)/mcapd.o: /usr/include/alloca.h ../mc_group.t
