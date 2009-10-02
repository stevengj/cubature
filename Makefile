FILES = README COPYING cubature.c cubature.h

test: cubature.c cubature.h
	cc -DTEST_INTEGRATOR -o $@ cubature.c -lm

dist:
	(d=cubature-`date +%Y%m%d`; rm -rf $$d $$d.tar.gz; mkdir $$d; cp $(FILES) $$d; tar czf $$d.tar.gz $$d; rm -rf $$d)
