# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
OBJ = $(LIB)(add_par_to_full_v2.o)\
      $(LIB)(ba_update_v6.o)\
      $(LIB)(ba_update_v7.o)\
      $(LIB)(cmf_blkband_v2.o)\
      $(LIB)(cmf_blkband_v3.o)\
      $(LIB)(create_iv_links_v2.o)\
      $(LIB)(determine_nse.o)\
      $(LIB)(do_ng_accel_band_v2.o)\
      $(LIB)(eval_adiabatic_v3.o)\
      $(LIB)(evalse_qwvj_v6.o)\
      $(LIB)(evalse_x_qwvj_v4.o)\
      $(LIB)(eval_temp_ddt_v1.o)\
      $(LIB)(fixpop_in_ba_v2.o)\
      $(LIB)(fixpop_t.o)\
      $(LIB)(fixpop_in_ba_v3.o)\
      $(LIB)(generate_full_matrix_v3.o)\
      $(LIB)(get_pops_at_prev_time_step_v2.o)\
      $(LIB)(prevent_low_t.o)\
      $(LIB)(read_ba_data_v2.o)\
      $(LIB)(read_bcd_mat.o)\
      $(LIB)(read_time_model_v1.o)\
      $(LIB)(set_ba_storage.o)\
      $(LIB)(set_imp_vec.o)\
      $(LIB)(solveba_v6.o)\
      $(LIB)(solveba_v7.o)\
      $(LIB)(solveba_v8.o)\
      $(LIB)(steq_multi_v7.o)\
      $(LIB)(steqne_v4.o)\
      $(LIB)(steq_advec_v4.o)\
      $(LIB)(steq_co_mov_deriv_v1.o)\
      $(LIB)(store_ba_data_v2.o)\
      $(LIB)(var_op_v8.o)\
      $(LIB)(var_x_opa_eta_v3.o)\
      $(LIB)(var_x_opa_eta_v4.o)\
      $(LIB)(vsebyj_multi_v6.o)\
      $(LIB)(vsebyj_x_v6.o)\
      $(LIB)(wr_asci_steq.o)

all : $(LIB) AUTO CHG TWO
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(OBJ)	

STEQ=$(LIB_MOD)(steq_data_mod.o)
CHG=$(LIB_MOD)(chg_exch_mod_v3.o)
MOD=$(LIB_MOD)(mod_cmfgen.o)
DIS=$(LIB_MOD)(mod_lev_dis_blk.o)

$(LIB)(%.o) : %.f  $(STEQ)
	$(F90) -c $(FG) $<
	ar ruv $(LIB) $*.o

# General rules to make obects in a library and object files
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(determine_nse.o) : determine_nse.f $(STEQ) $(CHG) $(MOD)
	$(F90) -c $(FG) determine_nse.f
	ar ruv $(LIB) determine_nse.o

$(LIB)(create_iv_links_v2.o) : create_iv_links_v2.f $(STEQ) $(MOD)
	$(F90) -c $(FG)  create_iv_links_v2.f
	ar ruv $(LIB) create_iv_links_v2.o

$(LIB)(fixpop_in_ba_v3.o) : fixpop_in_ba_v3.f $(STEQ) $(MOD)
	$(F90) -c $(FG) fixpop_in_ba_v3.f
	ar ruv $(LIB) fixpop_in_ba_v3.o

$(LIB)(generate_full_matrix_v3.o) : generate_full_matrix_v3.f $(STEQ) $(MOD)
	$(F90) -c $(FG)  generate_full_matrix_v3.f
	ar ruv $(LIB) generate_full_matrix_v3.o

$(LIB)(set_ba_storage.o) : set_ba_storage.f $(STEQ) $(MOD)
	$(F90) -c $(FG)  set_ba_storage.f
	ar ruv $(LIB) set_ba_storage.o

$(LIB)(eval_temp_ddt_v1.o) : eval_temp_ddt_v1.f $(STEQ) $(MOD)
	$(F90) -c $(FG)  eval_temp_ddt_v1.f
	ar ruv $(LIB) eval_temp_ddt_v1.o

$(LIB)(steq_co_mov_deriv_v1.o) : steq_co_mov_deriv_v1.f $(STEQ) $(MOD)
	$(F90) -c $(FG)  steq_co_mov_deriv_v1.f
	ar ruv $(LIB) steq_co_mov_deriv_v1.o

$(LIB)(var_op_v8.o) : var_op_v8.f $(DIS)
	$(F90) -c $(FG)  var_op_v8.f
	ar ruv $(LIB) var_op_v8.o

$(STEQ):
	(cd ../mod_subs ; make)

$(MOD):
	(cd ../mod_subs ; make)

$(DIS):
	(cd ../mod_subs ; make)

$(CHG):
	(cd chg ; make)

AUTO:
	(cd auto; make)
CHG:
	(cd chg; make)
TWO:
	(cd two; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
