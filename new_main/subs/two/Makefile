# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a

#
# Entries for library
#
OBJ = $(LIB)(two_phot_mod.o)\
      $(LIB)(steq_ba_two_phot_rate_v3.o)\
      $(LIB)(two_phot_rate.o)\
      $(LIB)(set_two_phot_v2.o)\
      $(LIB)(two_phot_opac.o)\
      $(LIB)(two_phot_var_opac.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(LIB)(two_phot_mod.o) 
	$(F90) -c  $(FG) $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(two_phot_mod.o) : two_phot_mod.f 
	$(F90) -c $(FG) two_phot_mod.f 
	ar ruv $(LIB) two_phot_mod.o

UP = $(LIB)(steq_data_mod.o)

$(LIB)(steq_ba_two_phot_rate_v3.o) : steq_ba_two_phot_rate_v3.f $(LIB)(two_phot_mod.o) $(UP) 
	$(F90) -c  $(FG) -I.. steq_ba_two_phot_rate_v3.f
	ar ruv $(LIB) steq_ba_two_phot_rate_v3.o

$(UP):
	(cd .. ; make mod)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
