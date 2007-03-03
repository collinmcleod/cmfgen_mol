# TEST MAKEFILE

#Include local system dependent definitions

include Makefile_definitions

#
# We will access the makefile in each local directory tp create the
# libraries and executables.
#
all : d_blas d_lpack d_tools d_unix d_subs d_plane d_newsubs d_pgplt d_disp\
         d_spec_plt d_new_main d_obs d_misc 

# We now MAKE the required libraries and executables.

d_blas:
	(cd blas; make)
d_lpack:
	(cd lpack; make)
d_plane:
	(cd plane; make )
d_tools:
	(cd tools; make )
d_subs:
	(cd subs; make )
d_unix:
	(cd unix; make )
d_newsubs:
	(cd newsubs; make ) 
d_pgplt:
	(cd pgplt; make ) 

# The following will create the executables

d_disp:
	(cd disp; make ) 
d_spec_plt:
	(cd spec_plt; make )
d_new_main:
	(cd new_main; make )
d_obs:
	(cd obs; make )
d_misc:
	(cd misc; make)
d_lte:
	(cd lte_hydro; make)
#
# To use the following command enter"
#         make -i clean
# The -i is necessary in case some files don't exist.

clean:
	rm -f lib/*.a
	rm -f */*.o
	rm -f */*/*.o
	rm -f */*.mod
	rm -f */*/*.mod
	rm -f */*/*/*.mod
	rm -f exe/*.exe
