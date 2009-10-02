FILES = README COPYING cubature.c cubature.h

test: cubature.c cubature.h
	cc -DTEST_INTEGRATOR -o $@ cubature.c -lm

dist:
	(d=cubature-`date +%Y%m%d`; rm -rf $$d $$d.tgz; mkdir $$d; cp $(FILES) $$d; tar czf $$d.tgz $$d; rm -rf $$d)
