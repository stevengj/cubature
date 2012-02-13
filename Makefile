FILES = README COPYING cubature.c cubature.h test.c ChangeLog

CFLAGS = -g -Wall -ansi -pedantic

all: cubature scubature pcubature

cubature: test.c cubature.c cubature.h
	cc $(CFLAGS) -o $@ test.c cubature.c -lm

pcubature: test.c pcubature.c cubature.h clencurt.h
	cc $(CFLAGS) -DPCUBATURE -o $@ test.c pcubature.c

scubature: test.c scubature.c cubature.h redblack.h
	cc $(CFLAGS) -DSCUBATURE -o $@ test.c scubature.c -lfftw3 -lm

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
	rm -f cubature scubature pcubature clencurt_gen

maintainer-clean:
	make clean
	rm -f clencurt.h
