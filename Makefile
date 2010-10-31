FILES = README COPYING cubature.c cubature.h test.c ChangeLog

CFLAGS = -g -Wall -ansi -pedantic

all: cubature scubature

cubature: test.c cubature.c cubature.h
	cc $(CFLAGS) -o $@ test.c cubature.c -lm

scubature: test.c scubature.c cubature.h redblack.h
	cc $(CFLAGS) -DSCUBATURE -o $@ test.c scubature.c -lfftw3 -lm

ChangeLog:
	darcs changes --summary > $@

dist:
	rm -f ChangeLog
	make ChangeLog
	(d=cubature-`date +%Y%m%d`; rm -rf $$d $$d.tgz; mkdir $$d; cp $(FILES) $$d; tar czf $$d.tgz $$d; rm -rf $$d)

clean:
	rm -f cubature scubature
