#Include file with system dependent definitions

include ../Makefile_definitions

#
# These names are for the fortran linker
#
LN1=tools
LN2=disp
LN3=subs
LN4=newsubs
LN5=pgplt
LN6=mod
LN7=dev
#
# These names refer to the full file name
#
FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a
FN5=$(LIB_DIR)lib$(LN5).a
FN6=$(LIB_DIR)lib$(LN6).a
FN7=$(LIB_DIR)lib$(LN7).a
#
# We only make dispgen.exe if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.
#
 
OBJ1  =  cmf_flux_v5.o cmf_flux_sub_v5.o mod_cmf_obs.o\
        rd_model_file.o  ins_line_v5.o obs_frame_sub_v5.o\
        mon_interp_fast.o ins_line_obs_v4.o rd_cont_j_obs.o 
 
OBJ2  =  obs_frame.o obs_frame_sub_v2.o mon_interp_fast.o
 
all: $(EXE_DIR)cmf_flux.exe $(EXE_DIR)obs_frame.exe

#
# Give the dependencies of dispgen.exe and how it is created.
#

$(EXE_DIR)cmf_flux.exe :  $(OBJ1) $(FN1) $(FN2) $(FN3) $(FN4) $(FN5) $(FN6) $(FN7)
	$(F90) $(OBJ1) \
              -L$(LIB_DIR) -l$(LN7) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN5) -l$(LN1) -l$(LN6) \
               $(LOCLIB) \
              -o $(EXE_DIR)cmf_flux.exe

$(EXE_DIR)obs_frame.exe :  $(OBJ2) $(FN1) $(FN2) $(FN3) $(FN4)
	$(F90) $(OBJ2) \
              -L$(LIB_DIR) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) \
               $(LOCLIB) \
              -o $(EXE_DIR)obs_frame.exe

# General rule to compile a fortran ".f" file.

%.o : %.f
	$(F90) -c $(FG)  $<

USR = ../lib/libtools.a(usr_option.o)
GEN = ../lib/libmod.a(mod_lev_dis_blk.o)

cmf_flux_v5.o     : cmf_flux_v5.f      mod_cmf_obs.o   $(USR)
cmf_flux_sub_v5.o : cmf_flux_sub_v5.f  mod_cmf_obs.o   $(GEN) \
   COMP_JCONT_V4.INC  EVAL_LTE_INC_V4.INC\
   OPACITIES_V4.INC   XOPAC_V4.INC

$(GEN):
	(cd ../subs; make )
$(USR):
	(cd ../tools; make )

 
$(FN1) : lib1
$(FN2) : lib2
$(FN3) : lib3
$(FN4) : lib4
$(FN5) : lib5
$(FN6) : lib6
$(FN7) : lib7

# We now MAKE the required libraries.

lib1:
	(cd ../tools; make)
lib2:
	(cd ../disp/subs; make)
lib3:
	(cd ../subs; make )
	(cd ../unix; make )
	(cd ../stark; make )
lib4:
	(cd ../newsubs; make)
lib5:
	(cd ../pgplt; make)
lib6:
	(cd ../new_main/mod_subs; make)
lib7:
	(cd ../new_main/subs; make)

clean :
	rm -f *.o
