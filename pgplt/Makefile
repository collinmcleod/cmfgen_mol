# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libpgplt.a

#
# Entries in a library
#
OBJ = $(LIB)(mod_curve_data.o)\
      $(LIB)(amoeba.o)\
      $(LIB)(amotry.o)\
      $(LIB)(change_color.o)\
      $(LIB)(change_pen.o)\
      $(LIB)(curve.o)\
      $(LIB)(do_vec_op.o)\
      $(LIB)(draw_right_y_axis.o)\
      $(LIB)(fntick.o)\
      $(LIB)(get_gramon_min_max.o)\
      $(LIB)(gramon_pgplot.o)\
      $(LIB)(hist_adj.o)\
      $(LIB)(justify_convert_v2.o)\
      $(LIB)(mon_num.o)\
      $(LIB)(monbord_v3.o)\
      $(LIB)(moninside.o)\
      $(LIB)(spacing.o)\
      $(LIB)(vectorpen.o)\
      $(LIB)(gaus_fit_data.o)\
      $(LIB)(do_gaus_fit.o)\
      $(LIB)(gaus_fit_func.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $<
	ar ruv $(LIB) $*.o

# 
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)
#
# NB: change_folor.f, change_pen.f, gramon_pgplot.f, vectorpen.f
#     depend on ../tools/gen_in

GEN = ../lib/libtools.a(gen_in.o)

$(LIB)(change_color.o)    : change_color.f  $(GEN)
$(LIB)(change_pen.o)      : change_pen.f    $(GEN)

$(LIB)(gramon_pgplot.o)   : gramon_pgplot.f $(LIB)(mod_curve_data.o) $(GEN) 
$(LIB)(curve.o)           : curve.f         $(LIB)(mod_curve_data.o) $(GEN)
$(LIB)(vectorpen.o)       : vectorpen.f     $(LIB)(mod_curve_data.o) $(GEN)

$(LIB)(gaus_fit_func.o)   : gaus_fit_func.f $(LIB)(gaus_fit_data.o)
$(GEN):
	(cd ../tools; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
