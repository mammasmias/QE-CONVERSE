# Makefile for qe-gipaw

QE_SOURCE=/home/sfioccola/Desktop/q-e-qe-7.2_clean
include $(QE_SOURCE)/make.inc
DFLAGS += 
INCLUDE_DIRS = -I$(QE_SOURCE)/KS_Solvers

# location of needed modules
BASEMOD_FLAGS= $(MOD_FLAG)$(QE_SOURCE)/upflib \
               $(MOD_FLAG)$(QE_SOURCE)/Modules \
               $(MOD_FLAG)$(QE_SOURCE)/MBD \
               $(MOD_FLAG)$(QE_SOURCE)/FFTXlib/src \
               $(MOD_FLAG)$(QE_SOURCE)/LAXlib \
               $(MOD_FLAG)$(QE_SOURCE)/UtilXlib \
               $(MOD_FLAG)$(QE_SOURCE)/XClib \
	       $(MOD_FLAG)$(QE_SOURCE)/external/devxlib/src \
	       $(MOD_FLAG)$(QE_SOURCE)/KS_Solvers

MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)$(QE_SOURCE)/PW/src \
           $(MOD_FLAG)$(QE_SOURCE)/dft-d3 \
           $(MOD_FLAG)$(QE_SOURCE)/LR_Modules \
           $(MOD_FLAG).

CONVERSE_OBJS = electrons_gipaw.o \
		set_dvrs.o \
                c_bands_gipaw.o \
                h_psi_gipaw.o \
                stop_code.o \
                gipaw_module.o \
                nmr_routine.o \
                nmr_module.o \
                gipaw_setup.o \
                init_us_2_no_phase.o \
                compute_u_kq.o \
                compute_dudk_new.o \
                wfcinit_gipaw.o \
                rotate_wfc_gipaw.o \
                rotate_wfc_k_gipaw.o \
                init_gipaw_1.o \
                init_gipaw_2.o \
		orbital_magnetization.o \
                calc_orbital_magnetization.o \
		newscf.o \
                paw_gipaw.o \
		qe-converse.o \
                util.o


LIBOBJS = $(QE_SOURCE)/UtilXlib/libutil.a

QEMODS = $(QE_SOURCE)/Modules/libqemod.a \
         $(QE_SOURCE)/MBD/libmbd.a \
         $(QE_SOURCE)/KS_Solvers/Davidson/libdavid.a \
         $(QE_SOURCE)/KS_Solvers/CG/libcg.a \
         $(QE_SOURCE)/KS_Solvers/PPCG/libppcg.a \
         $(QE_SOURCE)/KS_Solvers/libks_solvers.a \
         $(QE_SOURCE)/FFTXlib/src/libqefft.a \
         $(QE_SOURCE)/LAXlib/libqela.a \
         $(QE_SOURCE)/UtilXlib/libutil.a \
         $(QE_SOURCE)/upflib/libupf.a \
         $(QE_SOURCE)/dft-d3/libdftd3qe.a \
         $(QE_SOURCE)/LR_Modules/liblrmod.a \
         $(QE_SOURCE)/XClib/xc_lib.a


PWOBJS = $(QE_SOURCE)/PW/src/libpw.a

all: qe-converse.x

make_lr_modules:
	make -C $(QE_SOURCE)/LR_Modules

qe-converse.x: $(CONVERSE_OBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ $(CONVERSE_OBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	(cd ../bin; ln -sf ../src/qe-converse.x .)

clean:
	-/bin/rm -f qe-converse.x *.o *.F90 *__genmod.f90 *.d *.mod *.i *.L

distclean: clean
	-/bin/rm -f Makefile make.depend

xml_routines.f90: xmlmacros.m4 xml_routines.f90.in
	m4 xmlmacros.m4 xml_routines.f90.in >xml_routines.f90

include make.depend

# DO NOT DELETE

