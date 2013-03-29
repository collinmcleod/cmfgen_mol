# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)/libdisp.a

#
# Entries for library
#
OBJ = $(LIB)(betanew.o)\
      $(LIB)(branch_rat.o)\
      $(LIB)(chk_string.o)\
      $(LIB)(comp_jbar.o)\
      $(LIB)(compion.o)\
      $(LIB)(dlogvec.o)\
      $(LIB)(edgewr.o)\
      $(LIB)(fg_ham.o)\
      $(LIB)(fqcomp_v2.o)\
      $(LIB)(get_dtaulam.o)\
      $(LIB)(get_filename.o)\
      $(LIB)(get_line_indices.o)\
      $(LIB)(hamtuvgh.o)\
      $(LIB)(ham_formsol.o)\
      $(LIB)(ham_varlam.o)\
      $(LIB)(jfeaunew.o)\
      $(LIB)(map_plot.o)\
      $(LIB)(momham.o)\
      $(LIB)(new_writedc_v5.o)\
      $(LIB)(op_spec_file_v2.o)\
      $(LIB)(pallet.o)\
      $(LIB)(rd_ion_pop_v3.o)\
      $(LIB)(rd_non_therm_spec.o)\
      $(LIB)(rd_pop_dim.o)\
      $(LIB)(rd_rvtj.o)\
      $(LIB)(rd_rvtj_v2.o)\
      $(LIB)(recom_chk_v2.o)\
      $(LIB)(set_dc_or_pop_or_tx.o)\
      $(LIB)(set_dc_or_pop_or_tx_v2.o)\
      $(LIB)(setdc_or_pop.o)\
      $(LIB)(setforng.o)\
      $(LIB)(setrec.o)\
      $(LIB)(setvec.o)\
      $(LIB)(smooth_v2.o)\
      $(LIB)(sobjbar.o)\
      $(LIB)(starpcyg_tst.o)\
      $(LIB)(terminal_lu.o)\
      $(LIB)(update_pops.o)\
      $(LIB)(wrbeta.o)\
      $(LIB)(write_line.o)\
      $(LIB)(writscrat.o)\
      $(LIB)(wr_col.o)\
      $(LIB)(wr_col_rates.o)\
      $(LIB)(wr_spec_sum_v2.o)\
      $(LIB)(wr_crit.o)\
      $(LIB)(xvecfeau.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $<
	ar ruv $(LIB) $*.o

$(LIB)(wr_spec_sum_v2.o) : wr_spec_sum_v2.f ../mod_disp.o 
	$(F90) -c  $(FG) -I../ wr_spec_sum_v2.f
	ar ruv $(LIB) wr_spec_sum_v2.o

$(LIB)(branch_rat.o) : branch_rat.f ../mod_disp.o 
	$(F90) -c  $(FG) -I../ branch_rat.f
	ar ruv $(LIB) branch_rat.o

# 
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)
#

GEN = ../../lib/libsubs.a(mod_trap_derivatives.o)

$(LIB)(get_dtaulam.o)  : get_dtaulam.f  $(GEN)

$(GEN):
	(cd ../../subs; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
