#Include file with system dependent definitions

include ../Makefile_definitions

#
# These names are for the fortran linker
#
LN1=tools
LN2=subs
LN3=pgplt
LN4=newsubs

#
# These names refer to the full file name
#
FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a

#
# We only make dispgen.exe if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.
#
 
OBJ1  =  rd_f_to_s_ids.o rd_f_to_s_ids_v2.o rd_energy.o write_osc_v2.o
 
EXE   =  $(EXE_DIR)wr_f_to_s.exe\
         $(EXE_DIR)count_phot_data.exe\
         $(EXE_DIR)chk_phot.exe\
         $(EXE_DIR)pack_osc.exe\
         $(EXE_DIR)do_ng_v2.exe\
         $(EXE_DIR)land_col_merge.exe\
         $(EXE_DIR)land_multi_merge.exe\
         $(EXE_DIR)n_col_merge.exe\
         $(EXE_DIR)n_multi_merge.exe\
         $(EXE_DIR)rewrite_scr.exe\
         $(EXE_DIR)rewrite_dc.exe\
         $(EXE_DIR)rev_rvsig.exe\
         $(EXE_DIR)append_dc.exe\
         $(EXE_DIR)append_r_dc.exe\
         $(EXE_DIR)rev_hydro_file.exe\
         $(EXE_DIR)rev_hydro_turb.exe\
         $(EXE_DIR)mass_sc.exe\
         $(EXE_DIR)plt_gf.exe\
         $(EXE_DIR)plt_cool.exe\
         $(EXE_DIR)plt_j_comp.exe\
         $(EXE_DIR)plt_net.exe\
         $(EXE_DIR)solve_ba_mat.exe\
         $(EXE_DIR)tlusty_vel.exe
 
all: OBJ ALLEXE 

ALLEXE : $(EXE)
OBJ    :  $(OBJ1) $(FN1) $(FN2) $(FN3)

$(EXE_DIR)wr_f_to_s.exe : wr_f_to_s.o $(OBJ1)
	$(F90) wr_f_to_s.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN3) -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)wr_f_to_s.exe
#
$(EXE_DIR)chk_phot.exe : chk_phot.o recom_chk_v2.o get_exc_freq.o $(OBJ1)
	$(F90) chk_phot.o recom_chk_v2.o get_exc_freq.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN4) -l$(LN3) -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)chk_phot.exe
#
$(EXE_DIR)pack_osc.exe : pack_osc.o $(OBJ1)
	$(F90) pack_osc.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN4) -l$(LN3) -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)pack_osc.exe
#
# Give the dependencies of ?.exe and how it is created.
#
$(EXE_DIR)%.exe : %.o $(OBJ1) 
	$(F90) $*.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN3) -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)$*.exe

$(EXE_DIR)solve_ba_mat.exe : solve_ba_mat.o
	$(F90) solve_ba_mat.o -L$(LIB_DIR) $(LAPACK) $(BLAS) -l$(LN2) -l$(LN1) \
              -o $(EXE_DIR)solve_ba_mat.exe

# General rule to compile a fortran ".f" file.

%.o : %.f
	$(F90) -c $(FG)  $<

USR = $(LIB_DIR)libtools.a(usr_option.o)
GEN = $(LIB_DIR)libtools.a(gen_in.o)

wr_f_to_s.o     : wr_f_to_s.f  $(USR)

$(USR):
	(cd ../tools; make )
$(GEN):
	(cd ../tools; make )

 
$(FN1) : lib1
$(FN2) : lib2
$(FN3) : lib3
$(FN3) : lib4

# We now MAKE the required libraries.

lib1:
	(cd ../tools; make)

lib2:
	(cd ../subs; make )
	(cd ../unix; make )
	(cd ../stark; make )

lib3:
	(cd ../pgplt; make )

lib4:
	(cd ../newsubs; make )

clean :
	rm -f *.o
