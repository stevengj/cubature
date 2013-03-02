FILES = README COPYING pcubature.c hcubature.c cubature.h clencurt.h vwrapper.h converged.h test.c clencurt_gen.c ChangeLog NEWS

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
	./clencurt_gen 19 > $@

clencurt_gen: clencurt_gen.c
	cc $(CFLAGS) -o $@ clencurt_gen.c -lfftw3l -lm

ChangeLog:
	darcs changes --summary > $@

dist:
	rm -f ChangeLog
	make ChangeLog
	(d=cubature-`head -n 1 NEWS | cut -d' ' -f2`; rm -rf $$d $$d.tgz; mkdir $$d; cp $(FILES) $$d; tar czf $$d.tgz $$d; rm -rf $$d)

clean:
	rm -f htest ptest clencurt_gen *.o

dll32:
	make clean
	i586-mingw32msvc-gcc -c -O3 hcubature.c
	i586-mingw32msvc-gcc -c -O3 pcubature.c
	i586-mingw32msvc-gcc -shared -o libcubature32-`head -n 1 NEWS | cut -d' ' -f2`.dll hcubature.o pcubature.o
	make clean

dll64:
	make clean
	amd64-mingw32msvc-gcc -c -O3 hcubature.c
	amd64-mingw32msvc-gcc -c -O3 pcubature.c
	amd64-mingw32msvc-gcc -shared -o libcubature64-`head -n 1 NEWS | cut -d' ' -f2`.dll hcubature.o pcubature.o
	make clean

dylib64:
	make clean
	gcc -fPIC -c -O3 hcubature.c
	gcc -fPIC -c -O3 pcubature.c
	gcc -dynamiclib hcubature.o pcubature.o -o libcubature64-`head -n 1 NEWS | cut -d' ' -f2`.dylib
	make clean

dylib32:
	make clean
	gcc -m32 -fPIC -c -O3 hcubature.c
	gcc -m32 -fPIC -c -O3 pcubature.c
	gcc -m32 -dynamiclib hcubature.o pcubature.o -o libcubature32-`head -n 1 NEWS | cut -d' ' -f2-.dylib
	make clean

maintainer-clean:
	make clean
	rm -f clencurt.h
