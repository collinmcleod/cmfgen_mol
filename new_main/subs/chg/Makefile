# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored
# FX =  -Mextend -Mbackslash    -Ktrap=fp

LIB=$(LIB_DIR)libdev.a

#
# Entries in a library
#
OBJ = $(LIB)(chg_exch_mod_v3.o)\
      $(LIB)(eval_chg_rates_v3.o)\
      $(LIB)(set_chg_exch_v3.o)\
      $(LIB)(set_chg_lev_id_v3.o)\
      $(LIB)(verify_chg_exch_v3.o)\
      $(LIB)(steq_ba_chg_exch_v3.o)\
      $(LIB)(wr_chg_cool_v3.o)

all : $(LIB)
	rm -f *.o

$(LIB) : $(OBJ)	

UP = $(LIB_DIR)libdev.a(steq_data_mod.o)

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(LIB)(chg_exch_mod_v3.o) 
	$(F90) -c $(FG) $< 
	ar ruv $(LIB) $*.o

# Specific rules

mod: $(LIB)(steq_ba_chg_exch_v3.o)

$(LIB)(steq_ba_chg_exch_v3.o) : steq_ba_chg_exch_v3.f $(LIB)(chg_exch_mod_v3.o) $(UP) 
	$(F90) -c $(FG) -I.. steq_ba_chg_exch_v3.f
	ar ruv $(LIB) steq_ba_chg_exch_v3.o

$(LIB)(chg_exch_mod_v3.o) : chg_exch_mod_v3.f 
	$(F90) -c $(FG) chg_exch_mod_v3.f
	ar ruv $(LIB) chg_exch_mod_v3.o

$(UP):
	(cd .. ; make mod)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
