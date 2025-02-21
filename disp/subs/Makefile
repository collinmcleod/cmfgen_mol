# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)/libdisp.a


# Entries for library

OBJ = $(LIB)(betanew.o)\
      $(LIB)(arnaud_cross_v3_disp.o)\
      $(LIB)(branch_rat.o)\
      $(LIB)(comp_jbar.o)\
      $(LIB)(compion.o)\
      $(LIB)(determ_clump_pos.o)\
      $(LIB)(dlogvec.o)\
      $(LIB)(edgewr.o)\
      $(LIB)(fg_ham.o)\
      $(LIB)(fqcomp_v2.o)\
      $(LIB)(gam_abs.o)\
      $(LIB)(get_col_summary_v1.o)\
      $(LIB)(get_dtaulam.o)\
      $(LIB)(get_filename.o)\
      $(LIB)(get_flux_deficit.o)\
      $(LIB)(get_line_indices.o)\
      $(LIB)(hamtuvgh.o)\
      $(LIB)(ham_formsol.o)\
      $(LIB)(ham_varlam.o)\
      $(LIB)(jfeaunew.o)\
      $(LIB)(log_mon_interp.o)\
      $(LIB)(modify_density.o)\
      $(LIB)(map_plot.o)\
      $(LIB)(momham.o)\
      $(LIB)(new_writedc_v6.o)\
      $(LIB)(op_spec_file_v2.o)\
      $(LIB)(pallet.o)\
      $(LIB)(plt_arn.o)\
      $(LIB)(plt_profs.o)\
      $(LIB)(rd_ion_pop_v3.o)\
      $(LIB)(rd_non_therm_spec.o)\
      $(LIB)(rd_pop_dim.o)\
      $(LIB)(read_arnaud_ion_data_disp.o)\
      $(LIB)(recom_chk_v2.o)\
      $(LIB)(set_dc_or_pop_or_tx.o)\
      $(LIB)(set_dc_or_pop_or_tx_v2.o)\
      $(LIB)(set_two_phot_disp_v3.o)\
      $(LIB)(setdc_or_pop.o)\
      $(LIB)(set_fine_ray_grid.o)\
      $(LIB)(setforng.o)\
      $(LIB)(setrec.o)\
      $(LIB)(setvec.o)\
      $(LIB)(smooth_v2.o)\
      $(LIB)(sobjbar.o)\
      $(LIB)(starpcyg_tst.o)\
      $(LIB)(terminal_lu.o)\
      $(LIB)(tst_clump_law.o)\
      $(LIB)(two_phot_opac_disp_v3.o)\
      $(LIB)(update_pops.o)\
      $(LIB)(wrbeta.o)\
      $(LIB)(write_line.o)\
      $(LIB)(writscrat.o)\
      $(LIB)(wr_atom_sum.o)\
      $(LIB)(wr_col.o)\
      $(LIB)(wr_col_lines.o)\
      $(LIB)(wr_col_rates.o)\
      $(LIB)(wr_col_sl.o)\
      $(LIB)(wr_spec_sum_v2.o)\
      $(LIB)(wr_crit.o)\
      $(LIB)(write_line_12may98.o)\
      $(LIB)(write_line_old.o)\
      $(LIB)(write_line_multi.o)\
      $(LIB)(xvecfeau.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	$(NT)

NT =  $(LIB_DIR)libmod.a(mod_non_therm.o)

GEN =  $(LIB_DIR)libsubs.a(gen_in.o)\
	$(LIB_DIR)libsubs.a(usr_option.o)\
	$(LIB_DIR)libsubs.a(usr_hidden.o)\
	$(LIB_DIR)libsubs.a(mod_color_pen_def.o)\
	$(LIB_DIR)libmod.a(mod_lev_dis_blk.o)\
	$(LIB_DIR)libsubs.a(wr_string.o)

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $<
	ar $(AR_OPTION) $(LIB) $*.o

$(LIB)(%.o) : %.f90 ../mod_disp.o $(NT)
	$(F90) -c $(FFREE)  -I../ -I../../new_main/subs/non_therm $<
	ar $(AR_OPTION) $(LIB) $*.o

$(LIB)(branch_rat.o) : branch_rat.f ../mod_disp.o  $(GEN)
	$(F90) -c  $(FG) -I../ branch_rat.f
	ar $(AR_OPTION) $(LIB) branch_rat.o

$(LIB)(comp_jbar.o) : comp_jbar.f ../mod_disp.o  $(GEN)
	$(F90) -c  $(FG) -I../ comp_jbar.f
	ar $(AR_OPTION) $(LIB) comp_jbar.o

$(LIB)(disp_bethe_approx_v5.o) : disp_bethe_approx_v5.f ../mod_disp.o  $(GEN)
	$(F90) -c  $(FG) -I../ disp_bethe_approx_v5.f
	ar $(AR_OPTION) $(LIB) disp_bethe_approx_v5.o

$(LIB)(map_plot.o) : map_plot.f $(GEN)
	$(F90) -c  $(FD) -O1 map_plot.f 
	ar $(AR_OPTION) $(LIB) map_plot.o

$(LIB)(plt_profs.o) : plt_profs.f ../mod_disp.o $(GEN)
	$(F90) -c  $(FG) -I../ -I../../mod  plt_profs.f 
	ar $(AR_OPTION) $(LIB) plt_profs.o

$(LIB)(wr_atom_sum.o) : wr_atom_sum.f ../mod_disp.o $(GEN)
	$(F90) -c  $(FG) -I../ -I../../mod  wr_atom_sum.f 
	ar $(AR_OPTION) $(LIB) wr_atom_sum.o

$(LIB)(wr_spec_sum_v2.o) : wr_spec_sum_v2.f ../mod_disp.o $(GEN)
	$(F90) -c  $(FG) -I../ wr_spec_sum_v2.f
	ar $(AR_OPTION) $(LIB) wr_spec_sum_v2.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

TRAP = ../../lib/libsubs.a(mod_trap_derivatives.o)

$(LIB)(get_dtaulam.o)  : get_dtaulam.f  $(TRAP)

$(TRAP):
	(cd ../../subs; make )

$(GEN):
	(cd ../../tools; make )

$(NT):
	(cd ../../new_main/subs/non_therm; make )

clean :
	rm -f *.o
	rm -f $(LIB_DIR)/libdisp.a

clean_lib :
	rm -f $(LIB)
