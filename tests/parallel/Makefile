#Define variables

DIRS   = $(shell ls -d */ | grep -v CVS)
DIRS_A = $(shell ls -d */ | grep -v CVS | grep    "\-a")
DIRS_B = $(shell ls -d */ | grep -v CVS | grep -v "\-a")

#Define targets

.PHONY: test-a test-b rebuild-a rebuild-b setup

setup: ;
	-@for subdir in $(DIRS); do \
		rm $$subdir/adios_config.xml 2> /dev/null ; \
		ln -s ../../../bin/adios_config.xml $$subdir 2> /dev/null ; \
		ln -s -f ../make.common $$subdir/makefile 2> /dev/null ; \
		$(MAKE) --no-print-directory -e -C $$subdir $@; done

test-a rebuild-a: ;
	-@for subdir in $(DIRS_A) ; do \
		$(MAKE) --no-print-directory -e -C $$subdir $@;  done

test-b rebuild-b: ;
	-@for subdir in $(DIRS_B) ; do \
		$(MAKE) --no-print-directory -e -C $$subdir $@;  done

clean: ;
	-@for subdir in $(DIRS); do \
		$(MAKE) --no-print-directory -e -C $$subdir $@; done

