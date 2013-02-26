FILES = README COPYING pcubature.c hcubature.c cubature.h clencurt.h vwrapper.h converged.h test.c clencurt_gen.c ChangeLog

# CFLAGS = -pg -O3 -fno-inline-small-functions -Wall -ansi -pedantic
# CFLAGS = -g -Wall -ansi -pedantic
CFLAGS = -O3 -Wall -ansi -pedantic

all: htest ptest

htest: test.c hcubature.c cubature.h converged.h vwrapper.h
	cc $(CFLAGS) -o $@ test.c hcubature.c -lm

ptest: test.c pcubature.c cubature.h clencurt.h converged.h vwrapper.h
	cc $(CFLAGS) -DPCUBATURE -o $@ test.c pcubature.c -lm

clencurt.h: clencurt_gen.c # only depend on .c file so end-users don't re-gen
	make clencurt_gen
	./clencurt_gen 12 > $@

clencurt_gen: clencurt_gen.c
	cc $(CFLAGS) -o $@ clencurt_gen.c -lfftw3l -lm

ChangeLog:
	darcs changes --summary > $@

dist:
	rm -f ChangeLog
	make ChangeLog
	(d=cubature-`date +%Y%m%d`; rm -rf $$d $$d.tgz; mkdir $$d; cp $(FILES) $$d; tar czf $$d.tgz $$d; rm -rf $$d)

clean:
	rm -f htest ptest clencurt_gen

maintainer-clean:
	make clean
	rm -f clencurt.h
