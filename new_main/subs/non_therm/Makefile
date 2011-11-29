# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored
# FX =  -Mextend -Mbackslash    -Ktrap=fp

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries in a library
#
OBJ = $(LIB)(read_arnaud_ion_data.o)\
      $(LIB)(arnaud_cross_v3.o)\
      $(LIB)(bethe_approx_v3.o)\
      $(LIB)(intsigc.o)\
      $(LIB)(rd_nt_crosec_sclfac_v2.o)\
      $(LIB)(se_ba_non_therm_v2.o)\
      $(LIB)(set_xkt_array.o)\
      $(LIB)(get_lelec.o)\
      $(LIB)(total_bethe_rate_v2.o)\
      $(LIB)(electron_non_therm_spec.o)

OBJ2 = $(LIB_MOD)(mod_non_therm.o)

all : $(LIB)
	rm -f *.o

$(LIB) : $(OBJ2) $(OBJ)

NT =  $(LIB_DIR)libmod.a(mod_non_therm.o)
ST =  $(LIB_DIR)libmod.a(steq_data_mod.o)
CMF = $(LIB_DIR)libmod.a(mod_cmfgen.o)
CNT = $(LIB_DIR)libmod.a(control_variable_mod.o)

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f90 $(NT) $(ST) $(CMF) $(CNT)
	$(F90) -c $(FFREE) -I..  $<
	ar ruv $(LIB) $*.o

$(LIB)(mod_non_therm.o) : mod_non_therm.f90
	$(F90) -c $(FFREE) -I..  mod_non_therm.f90
	ar ruv $(LIB_MOD) mod_non_therm.o

$(NT):
	$(F90) -c $(FFREE) -I..  mod_non_therm.f90
	ar ruv $(LIB_DIR)libmod.a mod_non_therm.o

$(ST):
	(cd ../../mod_subs ; make)

$(CMF):
	(cd ../../mod_subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
