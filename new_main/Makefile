# TEST MAKEFILE

#Include file with system dependent definitions

include ../Makefile_definitions

#
# These names are for the fortran linker
#
LN1=tools
LN2=pgplt
LN3=subs
LN4=newsubs
LN5=dev
#
# These names refer to the full file name
#
FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a
FN5=$(LIB_DIR)lib$(LN5).a
#
# We only make dispgen.exe if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.
#
  
all : $(EXE_DIR)cmfgen_dev.exe

mod: mod_cmfgen.o

#
# Give the dependencies of dispgen.exe and how it is created.
#

$(EXE_DIR)cmfgen_dev.exe : cmfgen.o cmfgen_sub.o mod_cmfgen.o\
                                  $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) cmfgen.o cmfgen_sub.o mod_cmfgen.o \
             -L$(LIB_DIR) -l$(LN5) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) \
               $(LOCLIB) $(LAPACK) $(BLAS) \
              -o $(EXE_DIR)cmfgen_dev.exe

mod_cmfgen.o : mod_cmfgen.f 
	$(F90) -c $(FG) mod_cmfgen.f 

cmfgen.o : cmfgen.f mod_cmfgen.o
	$(F90) -c $(FG) cmfgen.f 

GEN = $(LIB_DIR)libsubs.a(mod_lev_dis_blk.o)
STEQ = $(LIB_DIR)libdev.a(steq_data_mod.o)

cmfgen_sub.o : cmfgen_sub.f mod_cmfgen.f $(GEN) $(STEQ) \
    COMP_JCONT_V4.INC    LINEGEN.INC         VARCONT.INC        XOPAC_V4.INC \
    EVAL_LTE_INC_V4.INC  OPACITIES_V4.INC    VAROPAC_V4.INC \
    GET_J_CHI.INC        SUP_TO_FULL_V4.INC  VAR_XOPACH_V4.INC
	$(F90) -c $(FG) -Isubs cmfgen_sub.f 

$(GEN):
	(cd ../subs; make )
 
$(STEQ):
	(cd subs; make mod)
 
$(FN1) : lib1
$(FN2) : lib2
$(FN3) : lib3
$(FN4) : lib4
$(FN5) : lib5

# We now MAKE the required libraries.

lib1:
	(cd ../tools; make )
lib2:
	(cd ../pgplt; make )
lib3:
	(cd ../subs; make )
lib4:
	(cd ../newsubs; make )
lib5:
	(cd subs; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
