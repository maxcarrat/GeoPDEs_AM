OCTFILES=basisfun.oct bspeval.oct nrb_srf_basisfun__.oct surfderivcpts.oct basisfunder.oct \
curvederivcpts.oct nrb_srf_basisfun_der__.oct surfderiveval.oct bspderiv.oct \
nrbsurfderiveval.oct tbasisfun.oct

MKOCTFILE ?= mkoctfile

all: $(OCTFILES)

low_level_functions.o: low_level_functions.cc
	$(MKOCTFILE) -c $<

%.oct:  %.cc low_level_functions.o
	$(MKOCTFILE) $< low_level_functions.o

clean:
	-rm -f *.o core octave-core *.oct *~

