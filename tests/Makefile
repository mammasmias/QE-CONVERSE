# Top-level test Makefile for QE-CONVERSE
#
# Targets:
#   make unit                 – build and run Fortran unit tests
#   make integration          – run integration tests
#   make all                  – unit + integration
#   make clean                – clean unit test build artefacts
#
# Variables:
#   QE_SRC_DIR   Path to the Quantum ESPRESSO source tree (contains bin/pw.x).
#                Passed in automatically by the top-level 'make test' target.

.PHONY: all unit integration clean

QE_SRC_DIR ?=

export QECONVERSE := $(abspath $(CURDIR)/../bin/qe-converse.x)
ifneq ($(QE_SRC_DIR),)
  export PW := $(abspath $(QE_SRC_DIR)/bin/pw.x)
endif

INTEGRATION_TESTS := integration/CO+ integration/quartz integration/nacl_nmr integration/CO+_kpools integration/nacl_nmr_kpool integration/licoo2_U

all: unit integration

unit:
	$(MAKE) -C unit run

integration:
	@for d in $(INTEGRATION_TESTS); do \
	    echo "=== $$d ==="; \
	    (cd $$d && bash run.sh) || exit 1; \
	done

clean:
	$(MAKE) -C unit clean
