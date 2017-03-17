SOURCEF = caxpy.f ccopy.f cdotc.f cdotu.f cgbmv.f cgemm.f cgemv.f	\
cgerc.f cgeru.f chbmv.f chemm.f chemv.f cher2.f cher2k.f cher.f		\
cherk.f chpmv.f chpr2.f chpr.f crotg.f cscal.f csrot.f csscal.f		\
cswap.f csymm.f csyr2k.f csyrk.f ctbmv.f ctbsv.f ctpmv.f ctpsv.f	\
ctrmm.f ctrmv.f ctrsm.f ctrsv.f dasum.f daxpy.f dcabs1.f dcopy.f	\
ddot.f dgbmv.f dgemm.f dgemv.f dger.f dnrm2.f drot.f drotg.f drotm.f	\
drotmg.f dsbmv.f dscal.f dsdot.f dspmv.f dspr2.f dspr.f dswap.f		\
dsymm.f dsymv.f dsyr2.f dsyr2k.f dsyr.f dsyrk.f dtbmv.f dtbsv.f		\
dtpmv.f dtpsv.f dtrmm.f dtrmv.f dtrsm.f dtrsv.f dzasum.f dznrm2.f	\
icamax.f idamax.f isamax.f izamax.f lsame.f sasum.f saxpy.f scabs1.f	\
scasum.f scnrm2.f scopy.f sdot.f sdsdot.f sgbmv.f sgemm.f sgemv.f	\
sger.f snrm2.f srot.f srotg.f srotm.f srotmg.f ssbmv.f sscal.f sspmv.f	\
sspr2.f sspr.f sswap.f ssymm.f ssymv.f ssyr2.f ssyr2k.f ssyr.f ssyrk.f	\
stbmv.f stbsv.f stpmv.f stpsv.f strmm.f strmv.f strsm.f strsv.f		\
xerbla.f zaxpy.f zcopy.f zdotc.f zdotu.f zdrot.f zdscal.f zgbmv.f	\
zgemm.f zgemv.f zgerc.f zgeru.f zhbmv.f zhemm.f zhemv.f zher2.f		\
zher2k.f zher.f zherk.f zhpmv.f zhpr2.f zhpr.f zrotg.f zscal.f zswap.f	\
zsymm.f zsyr2k.f zsyrk.f ztbmv.f ztbsv.f ztpmv.f ztpsv.f ztrmm.f	\
ztrmv.f ztrsm.f ztrsv.f

OBJSF = caxpy.o ccopy.o cdotc.o cdotu.o cgbmv.o cgemm.o cgemv.o		\
cgerc.o cgeru.o chbmv.o chemm.o chemv.o cher2.o cher2k.o cher.o		\
cherk.o chpmv.o chpr2.o chpr.o crotg.o cscal.o csrot.o csscal.o		\
cswap.o csymm.o csyr2k.o csyrk.o ctbmv.o ctbsv.o ctpmv.o ctpsv.o	\
ctrmm.o ctrmv.o ctrsm.o ctrsv.o dasum.o daxpy.o dcabs1.o dcopy.o	\
ddot.o dgbmv.o dgemm.o dgemv.o dger.o dnrm2.o drot.o drotg.o drotm.o	\
drotmg.o dsbmv.o dscal.o dsdot.o dspmv.o dspr2.o dspr.o dswap.o		\
dsymm.o dsymv.o dsyr2.o dsyr2k.o dsyr.o dsyrk.o dtbmv.o dtbsv.o		\
dtpmv.o dtpsv.o dtrmm.o dtrmv.o dtrsm.o dtrsv.o dzasum.o dznrm2.o	\
icamax.o idamax.o isamax.o izamax.o lsame.o sasum.o saxpy.o scabs1.o	\
scasum.o scnrm2.o scopy.o sdot.o sdsdot.o sgbmv.o sgemm.o sgemv.o	\
sger.o snrm2.o srot.o srotg.o srotm.o srotmg.o ssbmv.o sscal.o sspmv.o	\
sspr2.o sspr.o sswap.o ssymm.o ssymv.o ssyr2.o ssyr2k.o ssyr.o ssyrk.o	\
stbmv.o stbsv.o stpmv.o stpsv.o strmm.o strmv.o strsm.o strsv.o		\
xerbla.o zaxpy.o zcopy.o zdotc.o zdotu.o zdrot.o zdscal.o zgbmv.o	\
zgemm.o zgemv.o zgerc.o zgeru.o zhbmv.o zhemm.o zhemv.o zher2.o		\
zher2k.o zher.o zherk.o zhpmv.o zhpr2.o zhpr.o zrotg.o zscal.o zswap.o	\
zsymm.o zsyr2k.o zsyrk.o ztbmv.o ztbsv.o ztpmv.o ztpsv.o ztrmm.o	\
ztrmv.o ztrsm.o ztrsv.o

lib:
	$(FC) -c $(FOPTFLAGS) $(SOURCEF)
	$(AR) $(AR_FLAGS) ../$(LIBNAME) $(OBJSF)
	$(RM) $(OBJSF)
