# --- COMPILER ----------------------------------------
CC = mpiicc 
#intel mpiicc
#gnu mpicc
# --- CFLAGS -----------------------------------------
CFLAGS_gnu = -std=gnu99 -Wall -pedantic -O3 -ffast-math -msse4.2  -fopenmp -g3
CFLAGS_intel = -std=gnu99 -Wall -pedantic -O3 -xHost -qopenmp -g3 -fp-model fast=2 # -vec-report
CFLAGS = $(CFLAGS_intel)

# --- DO NOT CHANGE -----------------------------------
CPP = cpp
MAKEDEP = $(CPP) -MM
SRCDIR = src
BUILDDIR = build
BINDIR=bin
LIBDIR=lib
INCDIR=include
TSTDIR=tests
DOCDIR=doc
GSRCDIR = $(BUILDDIR)/gsrc
SRC = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.c,$(wildcard $(SRCDIR)/*.c)))
TSTS = $(patsubst %.c,%,$(wildcard $(TSTDIR)/*.c))
LIB =  $(LIBDIR)/libDDalphaAMG.a $(LIBDIR)/libDDalphaAMG_devel.a $(INCDIR)/DDalphaAMG.h
SRCGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.c))
GSRCFLT = $(patsubst %_generic.c,$(GSRCDIR)/%_float.c,$(SRCGEN))
GSRCDBL = $(patsubst %_generic.c,$(GSRCDIR)/%_double.c,$(SRCGEN))
GSRC = $(patsubst %,$(GSRCDIR)/%,$(SRC)) $(GSRCFLT) $(GSRCDBL)
HEA = $(patsubst $(SRCDIR)/%,%,$(filter-out %_generic.h,$(wildcard $(SRCDIR)/*.h)))
HEAGEN = $(patsubst $(SRCDIR)/%,%,$(wildcard $(SRCDIR)/*_generic.h))
GHEAFLT = $(patsubst %_generic.h,$(GSRCDIR)/%_float.h,$(HEAGEN))
GHEADBL = $(patsubst %_generic.h,$(GSRCDIR)/%_double.h,$(HEAGEN))
GHEA = $(patsubst %,$(GSRCDIR)/%,$(HEA)) $(GHEAFLT) $(GHEADBL)
OBJ = $(patsubst $(GSRCDIR)/%.c,$(BUILDDIR)/%.o,$(GSRC))
OBJDB = $(patsubst %.o,%_devel.o,$(OBJ))
DEP = $(patsubst %.c,%.dep,$(GSRC))

# --- FLAGS FOR HDF5 ---------------------------------
# H5FLAGS=-DHAVE_HDF5 /usr/include
# H5LIB=-lhdf5 -lz

# --- FLAGS FOR LIME ---------------------------------
LIMEDIR=${HOME}/src/c-lime/install
LIMEFLAGS=-DHAVE_LIME -I${LIMEDIR}/include
LIMELIB= -L${LIMEDIR}/lib -llime

# --- FLAGS FOR EFENCE --------------------------------- 
EFENCEDIR=${HOME}/src/electric-fence
EFENCEFLAGS=-I${EFENCEDIR}
EFENCELIB= # -L${EFENCEDIR} -lefence

# --- FLAGS FOR FABULOUS --------------------------------- 
FABULOUSDIR=${HOME}/src/fabulous
FABULOUSFLAGS=-DHAVE_FABULOUS -I${FABULOUSDIR}/src/api/include
FABULOUSLIB=-L${FABULOUSDIR}/build/src/api -lfabulous

# --- FLAGS FOR LAPACK ---------------------------------
LAPACKDIR=/onyx/buildsets/noe190301/software/imkl/2018.1.163-iimpi-2018a/mkl
#LAPACKDIR=/onyx/buildsets/noe190301/software/ScaLAPACK/2.0.2-gompi-2017b-OpenBLAS-0.2.20
LAPACKFLAGS=-I${LAPACKDIR}/include
#LAPACKLIB=-L${LAPACKDIR}/lib -l-lscalapack
LAPACKLIB=-L${LAPACKDIR}/lib/intel64 -L${LAPACKDIR}/mkl/lib/intel64 -lmkl

# Available flags:
# -DPARAMOUTPUT -DTRACK_RES -DFGMRES_RESTEST -DPROFILING
# -DSINGLE_ALLREDUCE_ARNOLDI
# -DCOARSE_RES -DSCHWARZ_RES -DTESTVECTOR_ANALYSIS -DDEBUG
# -DOPTIMIZE -DSSE -DAVX -DAVX2 -DAVX512
OPT_VERSION_FLAGS = $(CFLAGS) $(LIMEFLAGS) ${FABULOUSFLAGS} $(H5FLAGS) ${LAPACKFLAGS} ${EFENCEFLAGS} -DPARAMOUTPUT -DTRACK_RES -DOPENMP -DPROFILING #-DDEBUG
DEVEL_VERSION_FLAGS = $(CFLAGS) $(LIMEFLAGS) ${FABULOUSFLAGS} ${LAPACKFLAGS} ${EFENCEFLAGS} -DDEBUG -DPARAMOUTPUT -DTRACK_RES -DFGMRES_RESTEST -DPROFILING -DCOARSE_RES -DSCHWARZ_RES -DTESTVECTOR_ANALYSIS -DOPENMP

all: execs library exec-tests
execs: $(BINDIR)/DDalphaAMG $(BINDIR)/DDalphaAMG_devel
library: $(LIB)
exec-tests: $(TSTS)
documentation: $(DOCDIR)/user_doc.pdf
install: copy

.PHONY: all wilson library
.SUFFIXES:
.SECONDARY:

$(BINDIR)/DDalphaAMG : $(OBJ)
	$(CC) $(OPT_VERSION_FLAGS) -o $@ $(OBJ) $(H5LIB) $(LIMELIB) $(LAPACKLIB) $(EFENCELIB) $(FABULOUSLIB) -lm

DDalphaAMG : $(BINDIR)/DDalphaAMG
	ln -sf $(BINDIR)/$@ $@

DDalphaAMG_devel: $(BINDIR)/DDalphaAMG_devel
	ln -sf $(BINDIR)/$@ $@

$(BINDIR)/DDalphaAMG_devel : $(OBJDB)
	$(CC) -g $(DEVEL_VERSION_FLAGS) -o $@ $(OBJDB) $(H5LIB) $(LIMELIB) $(LAPACKLIB) $(EFENCELIB) $(FABULOUSLIB) -lm

$(LIBDIR)/libDDalphaAMG.a: $(OBJ)
	ar rc $@ $(OBJ)
	ar d $@ main.o
	ranlib $@

$(LIBDIR)/libDDalphaAMG_devel.a: $(OBJDB)
	ar rc $@ $(OBJDB)
	ar d $@ main.o
	ranlib $@

$(TSTDIR)/%: $(LIB) $(TSTDIR)/%.c
	$(CC) $(CFLAGS) -o $@ $@.c -I$(INCDIR) -L$(LIBDIR) -lDDalphaAMG $(LIMELIB) $(LAPACKLIB) $(EFENCELIB) $(FABULOUSLIB) -lm 

$(DOCDIR)/user_doc.pdf: $(DOCDIR)/user_doc.tex $(DOCDIR)/user_doc.bib
	( cd $(DOCDIR); pdflatex user_doc; bibtex user_doc; pdflatex user_doc; pdflatex user_doc; )

$(INCDIR)/%: $(SRCDIR)/%
	cp $(SRCDIR)/`basename $@` $@

$(BUILDDIR)/%.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) $(OPT_VERSION_FLAGS) -c $< -o $@

$(BUILDDIR)/%_devel.o: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) -g $(DEVEL_VERSION_FLAGS) -c $< -o $@

$(BUILDDIR)/%.optrpt: $(GSRCDIR)/%.c $(SRCDIR)/*.h
	$(CC) $(OPT_VERSION_FLAGS) -qopt-report=3 -c $< 

$(GSRCDIR)/%.h: $(SRCDIR)/%.h $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.h: $(SRCDIR)/%_generic.h $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

$(GSRCDIR)/%.c: $(SRCDIR)/%.c $(firstword $(MAKEFILE_LIST))
	cp $< $@

$(GSRCDIR)/%_float.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f float.sed $< > $@

$(GSRCDIR)/%_double.c: $(SRCDIR)/%_generic.c $(firstword $(MAKEFILE_LIST))
	sed -f double.sed $< > $@

%.dep: %.c $(GHEA)
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1.o $@ : ,g' > $@
	$(MAKEDEP) $< | sed 's,\(.*\)\.o[ :]*,$(BUILDDIR)/\1_devel.o $@ : ,g' >> $@

copy: $(BINDIR) $(LIBDIR) $(INCDIR) 
	cp -r $(BINDIR)/ $(PREFIX)
	cp -r $(LIBDIR)/ $(PREFIX)
	cp -r $(INCDIR)/ $(PREFIX)

clean:
	rm -f $(BUILDDIR)/*.o
	rm -f $(GSRCDIR)/*

cleanall: clean
	rm -f $(BINDIR)/*
	rm -f $(LIBDIR)/*
	rm -f $(INCDIR)/*
	rm -f $(SRCDIR)/*~
	rm -f $(TSTDIR)/*~
	rm -f $(TSTS)
	rm -f *~

-include $(DEP)
