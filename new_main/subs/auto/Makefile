# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a

#
# Entries for library
#
OBJ = $(LIB)(rd_auto_v1.o)\
      $(LIB)(steq_auto_v2.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

UP = $(LIB)(steq_data_mod.o)

$(LIB)(steq_auto_v2.o) : steq_auto_v2.f  $(UP) 
	$(F90) -c $(FG) -I.. steq_auto_v2.f
	ar ruv $(LIB) steq_auto_v2.o

$(UP):
	(cd .. ; make mod)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
