#
include ../Makefile_definitions

LIB = $(LIB_DIR)libmy_blas.a

OBJ =  chpr2.o   dcabs1.o  dsymm.o	izamax.o  sscal.o   xerbla.o  zhpmv.o \
caxpy.o   crotg.o   dcopy.o   dsymv.o	lsame.o   sspmv.o   zaxpy.o   zhpr.o \
ccopy.o   cscal.o   ddot.o    dsyr.o	sasum.o   sspr.o    zcopy.o   zhpr2.o \
cdotc.o   dgbmv.o   dsyr2.o	saxpy.o   sspr2.o   zdotc.o   zrotg.o \
cdotu.o   csscal.o  dgemm.o   dsyr2k.o	scasum.o  sswap.o   zdotu.o   zscal.o \
cgbmv.o   cswap.o   dgemv.o   dsyrk.o	scnrm2.o  ssymm.o   zswap.o \
cgemm.o   csymm.o   dger.o    dtbmv.o	scopy.o   ssymv.o   zdscal.o  zsymm.o \
cgemv.o   csyr2k.o  dnrm2.o   dtbsv.o	sdot.o	  ssyr.o    zgbmv.o   zsyr2k.o \
cgerc.o   csyrk.o   drot.o    dtpmv.o   ssyr2.o   zgemm.o   zsyrk.o \
cgeru.o   ctbmv.o   drotg.o   dtpsv.o	sgbmv.o   ssyr2k.o  zgemv.o   ztbmv.o \
chbmv.o   ctbsv.o   dtrmm.o	sgemm.o   ssyrk.o   zgerc.o   ztbsv.o \
chemm.o   ctpmv.o   dtrmv.o	sgemv.o   stbmv.o   zgeru.o   ztpmv.o \
chemv.o   ctpsv.o   dsbmv.o   dtrsm.o	sger.o	  stbsv.o   zhbmv.o   ztpsv.o \
cher.o	  ctrmm.o   dscal.o   dtrsv.o	snrm2.o   stpmv.o   zhemm.o   ztrmm.o \
cher2.o   ctrmv.o   dzasum.o  srot.o	 stpsv.o   zhemv.o   ztrmv.o \
cher2k.o  ctrsm.o   dspmv.o   dznrm2.o	srotg.o   strmm.o   zher.o    ztrsm.o \
cherk.o   ctrsv.o   dspr.o    icamax.o  strmv.o   zher2.o   ztrsv.o \
chpmv.o   dasum.o   dspr2.o   idamax.o  strsm.o   zher2k.o \
chpr.o	  daxpy.o   dswap.o   isamax.o	ssbmv.o   strsv.o   zherk.o 


all : $(LIB)
	rm -f *.o

$(LIB) : $(LIB)($(OBJ))

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f
	$(F90) -c $(FG)  $<
	ar $(AR_OPTION) $(LIB) $*.o

#
# Specific dependency of source files which uses external
# fortran moudles. (The order of compilation matters.)


clean :
	rm -f *.o

clean_lib :
	rm -f $(NIB)
