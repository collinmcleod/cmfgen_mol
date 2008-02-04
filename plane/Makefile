#
# Include local system definitions

include ../Makefile_definitions

OPTION=$(FD)

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a
#
# Entries for library
#

OBJ = $(LIB)(comp_add_edd_facs.o)\
      $(LIB)(drhsdchi_pp.o)\
      $(LIB)(expn.o)\
      $(LIB)(gauleg.o)\
      $(LIB)(get_ibound.o)\
      $(LIB)(edd_j_hub_var_v1.o)\
      $(LIB)(edd_jrel_var_v1.o)\
      $(LIB)(fgrey_norel_v1.o)\
      $(LIB)(jgrey_hub_ddt_v2.o)\
      $(LIB)(out_jh.o)\
      $(LIB)(get_jh_at_prev_time_step.o)\
      $(LIB)(mom_j_ddt_v1.o)\
      $(LIB)(mom_jrel_v3.o)\
      $(LIB)(mom_jrel_grey_v1.o)\
      $(LIB)(var_mom_j_ddt_v1.o)\
      $(LIB)(pp_form_cmf_v2.o)\
      $(LIB)(pp_var_mom_cmf_v1.o)\
      $(LIB)(pp_edd_var_cmf_v1.o)\
      $(LIB)(pp_mom_cmf_v1.o)\
      $(LIB)(fcomp_pp_v2.o)\
      $(LIB)(mom_j_pp_v1.o)\
      $(LIB)(var_jrel_v1.o)\
      $(LIB)(var_mom_pp_v1.o)


all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(OPTION) $<
	ar ruv $(LIB) $*.o

%.o : %.f 
	$(F90) -c $(OPTION) $<
	ar ruv $(LIB) $*.o

$(LIB)(drhsdchi_pp.o) : drhsdchi_pp.f $(LIB_MOD)(mod_trap_derivatives.o)
	$(F90) -c $(OPTION) drhsdchi_pp.f
	ar ruv $(LIB) drhsdchi_pp.o

$(LIB)(var_jrel_v1.o) : var_jrel_v1.f $(LIB)(mod_var_jrel.o)
	$(F90) -c $(OPTION) var_jrel_v1.f
	ar ruv $(LIB) var_jrel_v1.o

$(LIB)(edd_jrel_var_v1.o) : edd_jrel_var_v1.f $(LIB)(mod_var_jrel.o)
	$(F90) -c $(OPTION) edd_jrel_var_v1.f
	ar ruv $(LIB) edd_jrel_var_v1.o

$(LIB_MOD)(mod_trap_derivatives.o) :
	(cd ../subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm -f *.a
