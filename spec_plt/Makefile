#Include file with system dependent definitions

include ../Makefile_definitions

#
# These names are for the fortran linker
#
LN1=tools
LN2=pgplt
LN3=subs
LN4=disp
LN5=spec
#
# These names refer to the full file name
#
FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a
FN5=$(LIB_DIR)lib$(LN5).a
#
EXE  = $(EXE_DIR)plt_jh.exe\
       $(EXE_DIR)plt_jh_cur.exe\
       $(EXE_DIR)plt_ip.exe\
       $(EXE_DIR)plt_phot_raw.exe\
       $(EXE_DIR)plt_rj.exe\
       $(EXE_DIR)plt_scr.exe\
       $(EXE_DIR)new_plt_rj.exe\
       $(EXE_DIR)plt_spec.exe

all : OBJ ALLEXE


# Give the dependencies of dispgen.exe and how it is created.
#

OBJ       : $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
ALLEXE    : $(EXE)

#
# A MACRO deinition could probably be used here, but it is difficult to get it to work on all
# UNIX systems, hence we take the easy but tediuos method.
#

$(EXE_DIR)plt_spec.exe : plt_spec.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90)  plt_spec.o \
             -L$(LIB_DIR) -l$(LN5) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)plt_spec.exe

$(EXE_DIR)plt_jh.exe : plt_jh.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) plt_jh.o \
             -L$(LIB_DIR) -lspec  -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)plt_jh.exe

$(EXE_DIR)plt_jh_cur.exe : plt_jh_cur.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) plt_jh_cur.o \
             -L$(LIB_DIR) -lspec  -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)plt_jh_cur.exe

$(EXE_DIR)plt_ip.exe : plt_ip.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) plt_ip.o \
             -L$(LIB_DIR) -l$(LN5) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)plt_ip.exe

$(EXE_DIR)plt_phot_raw.exe : plt_phot_raw.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) plt_phot_raw.o \
             -L$(LIB_DIR) -l$(LN5)  -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)plt_phot_raw.exe

$(EXE_DIR)new_plt_rj.exe : new_plt_rj.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) new_plt_rj.o \
             -L$(LIB_DIR) -l$(LN5)  -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)new_plt_rj.exe

$(EXE_DIR)plt_scr.exe : plt_scr.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) plt_scr.o \
             -L$(LIB_DIR) -l$(LN5)  -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)plt_scr.exe

$(EXE_DIR)plt_rj.exe : plt_rj.o $(FN1) $(FN2) $(FN3) $(FN4) $(FN5)
	$(F90) plt_rj.o \
             -L$(LIB_DIR) -l$(LN5)  -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) -l$(LN3) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)plt_rj.exe

%.o : %.f
	$(F90) -c $(FG) $<

SUBS:
	(cd subs; make)

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
	(cd ../disp/subs; make ) 
lib5:
	(cd subs; make ) 

#clean :
#	rm -f *.o

clean_lib :
	rm -f $(LIB)
