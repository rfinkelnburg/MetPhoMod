OBJS= $(OBJDIR)/mctwostream.o $(OBJDIR)/mcphotodiss.o \
      $(OBJDIR)/mcabsorb.o $(OBJDIR)/mcextrad.o\
      $(OBJDIR)/mc_inter3.o $(OBJDIR)/mc_rtlink.o $(OBJDIR)/mc_schu.o \
      $(OBJDIR)/mc_seto2.o $(OBJDIR)/mc_sphers.o $(OBJDIR)/mc_ps2str.o \
      $(OBJDIR)/mc_fchap.o $(OBJDIR)/mcapd.o

SRC= mcphotodiss.c mctwostream.c mcabsorb.c mcextrad.c mcapd.c

CC=xlC -+

CFLAGS2=-I..

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(GENFLAGS) $(INCLUDE) -DFERRET_FORMAT $*.c -o $@

$(OBJDIR)/%.o : %.cc
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(GENFLAGS) $(INCLUDE) -DFERRET_FORMAT $*.cc -o $@

$(OBJDIR)/%.o : %.f
	f77 -qextname -c $(FFLAGS) $(GENFLAGS) $*.f -o $@

libtwostream.a : $(OBJS)
	ar -r $(OBJDIR)/libmctwostream.a $(OBJS)

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
$(OBJDIR)/mcphotodiss.o: /usr/include/standards.h ../mc_group.hh
$(OBJDIR)/mcphotodiss.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcphotodiss.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/context.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/m_param.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/localedef.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcphotodiss.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mcphotodiss.o: /usr/include/math.h ../mc_group.t mcextrad.h
$(OBJDIR)/mcphotodiss.o: ../mcglobal.h ph_o3data.c ph2data.c ph3data.c
$(OBJDIR)/mcphotodiss.o: ph4data.c ph6data.c ph7data.c ph11data.c ph12data.c
$(OBJDIR)/mcphotodiss.o: ph13data.c ph14data.c
$(OBJDIR)/mctwostream.o: mctwostream.h ../mcground.h ../mcglobal.h
$(OBJDIR)/mctwostream.o: ../sunpos.h ../mc_module.hh ../mc_variable.hh
$(OBJDIR)/mctwostream.o: ../mcparse.h ../mc_group.hh ../mc_section.hh
$(OBJDIR)/mctwostream.o: ../mc_commands.hh mcphotodiss.h mcapd.h
$(OBJDIR)/mctwostream.o: /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mctwostream.o: /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mctwostream.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/context.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/m_param.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/mstsave.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/localedef.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mctwostream.o: /usr/include/sys/localedef31.h /usr/include/string.h
$(OBJDIR)/mctwostream.o: /usr/include/math.h ../mchemparse.h ../mchem.h
$(OBJDIR)/mctwostream.o: mcabsorb.h mcextrad.h
$(OBJDIR)/mcabsorb.o: mcabsorb.h /usr/include/stddef.h
$(OBJDIR)/mcabsorb.o: /usr/include/standards.h /usr/include/stdio.h
$(OBJDIR)/mcabsorb.o: /usr/include/va_list.h /usr/include/sys/types.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/inttypes.h /usr/include/sys/m_types.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mcabsorb.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/signal.h /usr/include/sys/context.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/m_param.h /usr/include/sys/mstsave.h
$(OBJDIR)/mcabsorb.o: /usr/include/sys/localedef.h /usr/include/sys/lc_core.h
$(OBJDIR)/mcabsorb.o: /usr/include/locale.h /usr/include/sys/localedef31.h
$(OBJDIR)/mcabsorb.o: ../mcparse.h ../mc_group.hh mcextrad.h ph_o3data.c
$(OBJDIR)/mcextrad.o: mcextrad.h /usr/include/stddef.h
$(OBJDIR)/mcextrad.o: /usr/include/standards.h /usr/include/stdio.h
$(OBJDIR)/mcextrad.o: /usr/include/va_list.h /usr/include/sys/types.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/inttypes.h /usr/include/sys/m_types.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/limits.h /usr/include/float.h
$(OBJDIR)/mcextrad.o: /usr/include/stdlib.h /usr/include/sys/wait.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/resource.h /usr/include/sys/time.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/signal.h /usr/include/sys/context.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/m_param.h /usr/include/sys/mstsave.h
$(OBJDIR)/mcextrad.o: /usr/include/sys/localedef.h /usr/include/sys/lc_core.h
$(OBJDIR)/mcextrad.o: /usr/include/locale.h /usr/include/sys/localedef31.h
$(OBJDIR)/mcextrad.o: ../mcparse.h ../mc_group.hh
$(OBJDIR)/mcapd.o: mcapd.h /usr/include/stddef.h /usr/include/standards.h
$(OBJDIR)/mcapd.o: ../mc_group.hh /usr/include/stdio.h /usr/include/va_list.h
$(OBJDIR)/mcapd.o: /usr/include/sys/types.h /usr/include/sys/inttypes.h
$(OBJDIR)/mcapd.o: /usr/include/sys/m_types.h /usr/include/sys/limits.h
$(OBJDIR)/mcapd.o: /usr/include/float.h /usr/include/stdlib.h
$(OBJDIR)/mcapd.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
$(OBJDIR)/mcapd.o: /usr/include/sys/time.h /usr/include/sys/signal.h
$(OBJDIR)/mcapd.o: /usr/include/sys/context.h /usr/include/sys/m_param.h
$(OBJDIR)/mcapd.o: /usr/include/sys/mstsave.h /usr/include/sys/localedef.h
$(OBJDIR)/mcapd.o: /usr/include/sys/lc_core.h /usr/include/locale.h
$(OBJDIR)/mcapd.o: /usr/include/sys/localedef31.h ../mc_group.t
