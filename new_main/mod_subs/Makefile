# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
MOD = $(LIB_MOD)(mod_cmfgen.o)\
      $(LIB_MOD)(ang_qw_mod.o)\
      $(LIB_MOD)(opac_mod.o)\
      $(LIB_MOD)(control_variable_mod.o)\
      $(LIB_MOD)(steq_data_mod.o)\
      $(LIB_MOD)(two_phot_mod.o)\
      $(LIB_MOD)(mod_lev_dis_blk.o)\
      $(LIB_MOD)(line_vec_mod.o)\
      $(LIB_MOD)(chg_exch_mod_v3.o)

OBJ = $(LIB)(adjust_r_grid_v2.o)\
      $(LIB)(adjust_den_vecs.o)\
      $(LIB)(comp_opac.o)\
      $(LIB)(comp_var_opac.o)\
      $(LIB)(init_opac_mod.o)\
      $(LIB)(rd_control_variables.o)\
      $(LIB)(set_frequency_grid.o)\
      $(LIB)(set_ang_qw.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(MOD) $(OBJ)	

# General rules to make obects in a library and object files

$(LIB_MOD)(%.o) : %.f 
	$(F90) -c $(FG) $< 
	ar ruv $(LIB_MOD) $*.o

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG) $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(comp_opac.o) : comp_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) comp_opac.f
	ar ruv $(LIB) comp_opac.o

$(LIB)(comp_var_opac.o) : comp_var_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) comp_var_opac.f
	ar ruv $(LIB) comp_var_opac.o

$(LIB)(rd_control_variables.o) : rd_control_variables.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) rd_control_variables.f
	ar ruv $(LIB) rd_control_variables.o

clean :
	rm -f *.o 
