FILES = README COPYING cubature.c cubature.h test.c ChangeLog

# CFLAGS = -pg -O3 -fno-inline-small-functions -Wall -ansi -pedantic
# CFLAGS = -g -Wall -ansi -pedantic
CFLAGS = -O3 -Wall -ansi -pedantic

all: hcubature pcubature

hcubature: test.c hcubature.c cubature.h converged.c
	cc $(CFLAGS) -o $@ test.c hcubature.c -lm

pcubature: test.c pcubature.c cubature.h clencurt.h converged.c
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
	rm -f hcubature scubature pcubature clencurt_gen

maintainer-clean:
	make clean
	rm -f clencurt.h
