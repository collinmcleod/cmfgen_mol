# Include local system definitions

include ../Makefile_definitions

#
# These names are for the fortran linker
#
LN1=mod
LN2=subs
LN3=pgplt
LN4=disp
#
# These names refer to the full file name
#
FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a

#
# We only make dispgen.exe if there has been a changen to the main *.f files,
# or if the libraries have chnaged since the last link.
# 

all : $(EXE_DIR)dispgen.exe

#
# Give the depeedencies of dispgen.exe and how it is created.
#

$(EXE_DIR)dispgen.exe : mod_disp.o dispgen.o maingen.o pltphot_sub.o $(FN1) $(FN2) $(FN3) $(FN4)
	$(F90) dispgen.o mod_disp.o maingen.o pltphot_sub.o \
             -L$(LIB_DIR) -l$(LN4) -l$(LN3) -l$(LN2) -l$(LN1)  \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)/dispgen.exe

mod_disp.o : mod_disp.f
	$(F90) -c $(FG) mod_disp.f 

#Module dependencies for DISPGEN, MAINGEN

GEN =  $(LIB_DIR)libsubs.a(gen_in.o)\
               $(LIB_DIR)libsubs.a(usr_option.o)\
               $(LIB_DIR)libsubs.a(usr_hidden.o)\
               $(LIB_DIR)libsubs.a(mod_color_pen_def.o)\
               $(LIB_DIR)libmod.a(mod_lev_dis_blk.o)\
               $(LIB_DIR)libsubs.a(wr_string.o)

dispgen.o : dispgen.f mod_disp.o $(GEN)
	$(F90) -c $(FG) dispgen.f 


maingen.o : maingen.f mod_disp.o $(GEN) \
               EVAL_LTE_FULL.INC  OPACITIES.INC \
               PAR_OPACITIES.INC  XOPAC.INC PAR_XOPAC.INC
	$(F90) -c  $(FD) maingen.f 

pltphot_sub.o : pltphot_sub.f mod_disp.o $(GEN)
	$(F90) -c  $(FD) pltphot_sub.f 

$(GEN):
	(cd ../tools; make )
	(cd ../new_main/mod_subs; make )


# We now MAKE the required libraries. Note that this double
# dependence is required to ensure that the libraries are updated.
# Otherwise MAKE simple checks the date of the library.
 
$(FN1): lib1
$(FN2): lib2
$(FN3): lib3
$(FN4): lib4

lib1:
	(cd ../new_main/mod_subs; make )
lib2:
	(cd ../tools; make )
	(cd ../subs; make )
	(cd ../newsubs; make )
lib3:
	(cd ../pgplt; make )
lib4:
	(cd subs; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
