# Top level makefile for qe-converse

all:    build

build:
	@echo "Building qe-converse..."
	$(MAKE) -C src

clean:
	@echo "Cleaning qe-converse..."
	if test -s src/Makefile ; then ( $(MAKE) -C src clean ); fi
	-/bin/rm -f bin/qe-converse.x

distclean:
	$(MAKE) -C src distclean
	-/bin/rm -f config.log config.status makedeps.sh
	-/bin/rm -Rf autom4te.cache

