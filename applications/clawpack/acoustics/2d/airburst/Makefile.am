## Makefile.am in applications/clawpack/acoustics_2d_airburst

bin_PROGRAMS += applications/clawpack/acoustics/2d/airburst/airburst

applications_clawpack_acoustics_2d_airburst_airburst_SOURCES = \
	applications/clawpack/acoustics/2d/airburst/airburst_user.cpp \
	applications/clawpack/acoustics/2d/airburst/airburst_user.h \
	applications/clawpack/acoustics/2d/airburst/airburst.cpp \
	applications/clawpack/acoustics/2d/airburst/fclaw2d_map_pillowdisk5.c \
	applications/clawpack/acoustics/2d/airburst/setprob.f \
	applications/clawpack/acoustics/2d/airburst/computedOverPressure_250MT.f90 \
	applications/clawpack/acoustics/2d/airburst/user_4.6/colaux.f90 \
	applications/clawpack/acoustics/2d/airburst/user_4.6/qinit.f \
	applications/clawpack/acoustics/2d/airburst/user_4.6/setaux.f \
	applications/clawpack/acoustics/2d/airburst/user_4.6/bc2.f \
	applications/clawpack/acoustics/2d/airburst/user_4.6/src2.f \
	applications/clawpack/acoustics/2d/airburst/user_5.0/tag4refinement.f90 \
	applications/clawpack/acoustics/2d/airburst/user_5.0/tag4coarsening.f90 \
	applications/clawpack/acoustics/2d/airburst/user_5.0/setaux.f90 \
	applications/clawpack/acoustics/2d/rp/rpn2ac.f \
	applications/clawpack/acoustics/2d/rp/acoustics_user_fort.h \
	applications/clawpack/acoustics/2d/rp/rpt2ac.f \
	applications/clawpack/acoustics/2d/rp/rpn2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/rpt2_acoustics.f90 \
	applications/clawpack/acoustics/2d/rp/rpn2_acoustics_mapped.f90 \
	applications/clawpack/acoustics/2d/rp/rpt2_acoustics_mapped.f90


applications_clawpack_acoustics_2d_airburst_airburst_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_CLAWPACK46_CPPFLAGS) $(FCLAW_CLAWPACK5_CPPFLAGS)

applications_clawpack_acoustics_2d_airburst_airburst_LDADD = \
        $(FCLAW_CLAWPACK46_LDADD) \
        $(FCLAW_CLAWPACK5_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(LDADD)
