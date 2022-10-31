VERSIONS = std gam ncl gpu gpu_ncl

.PHONY: all veryclean versions $(VERSIONS)

all: std gam ncl

versions: $(VERSIONS)

$(VERSIONS):
	if [ ! -d build/$@ ] ; then mkdir build/$@  ; fi; \
	cp src/makefile src/.objects makefile.include build/$@ ; \
	if [ -f src/makedeps.awk ] ; then cp src/makedeps.awk build/$@ ; fi; \
	$(MAKE) -C build/$@ VERSION=$@ dependencies -j1 ; \
	$(MAKE) -C build/$@ VERSION=$@ all

veryclean:
	rm -rf build/*

test:
	if [ -f ./testsuite/tests/cleanall ] ; then \
	  cd ./testsuite/tests ; ./cleanall ; fi
	$(MAKE) -C testsuite test

test_all:
	if [ -f ./testsuite/tests/cleanall ] ; then \
	  cd ./testsuite/tests ; ./cleanall ; fi
	$(MAKE) -C testsuite test_all
