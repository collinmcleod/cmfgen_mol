# Include local system definitions

include ../Makefile_definitions

#
# These names are for the fortran linker
#
LN1=tools
LN2=pgplt
LN3=subs
LN4=newsubs
LN5=disp
#
# These names refer to the full file name
#
FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a
FN5=$(LIB_DIR)lib$(LN5).a

#
# We only make dispgen.exe if there has been a changen to the main *.f files,
# or if the libraries have chnaged since the last link.
# 

all : $(EXE_DIR)dispgen.exe

#
# Give the depeedencies of dispgen.exe and how it is created.
#

$(EXE_DIR)dispgen.exe : mod_disp.o dispgen.o maingen.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5) 
	$(F90) dispgen.o mod_disp.o maingen.o \
             -L$(LIB_DIR) -l$(LN5) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)/dispgen.exe

mod_disp.o : mod_disp.f
	$(F90) -c $(FG) mod_disp.f 

# Module dependencies for DISPGEN

DGEN =  $(LIB_DIR)libtools.a(gen_in.o)\
             $(LIB_DIR)libtools.a(usr_option.o)

dispgen.o : dispgen.f mod_disp.o $(DGEN)
	$(F90) -c $(FG) dispgen.f 

#Module dependencies for MAINGEN

MGEN =  $(LIB_DIR)libtools.a(gen_in.o) $(LIB_DIR)libtools.a(usr_option.o)\
               $(LIB_DIR)libtools.a(usr_hidden.o)\
               $(LIB_DIR)libsubs.a(mod_lev_dis_blk.o)

maingen.o : maingen.f mod_disp.o $(MGEN)
	$(F90) -c $(FG) maingen.f 

$(DGEN):
	(cd ../tools; make )

$(MGEN):
	(cd ../tools; make )
	(cd ../subs; make )

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
