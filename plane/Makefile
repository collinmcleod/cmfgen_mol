#
# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a
#
# Entries for library
#

OBJ = $(LIB)(drhsdchi_pp.o)\
      $(LIB)(expn.o)\
      $(LIB)(gauleg.o)\
      $(LIB)(get_ibound.o)\
      $(LIB)(pp_form_cmf_v2.o)\
      $(LIB)(pp_var_mom_cmf_v1.o)\
      $(LIB)(pp_edd_var_cmf_v1.o)\
      $(LIB)(pp_mom_cmf_v1.o)\
      $(LIB)(fcomp_pp_v2.o)\
      $(LIB)(mom_j_pp_v1.o)\
      $(LIB)(var_mom_pp_v1.o)


all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG) $<
	ar ruv $(LIB) $*.o

$(LIB)(drhsdchi_pp.o) : drhsdchi_pp.f $(LIB_MOD)(mod_trap_derivatives.o)
	$(F90) -c $(FG) drhsdchi_pp.f
	ar ruv $(LIB) drhsdchi_pp.o

$(LIB_MOD)(mod_trap_derivatives.o) :
	(cd ../subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm -f *.a
