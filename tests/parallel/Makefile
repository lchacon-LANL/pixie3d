#Define variables

DIRS = $(shell ls | grep -v ake | grep -v CVS)

#Define targets

.PHONY: test rebuild setup

setup: ;
	-@for subdir in $(DIRS); do \
		ln -s ../make.common $$subdir/makefile; \
		$(MAKE) --no-print-directory -C $$subdir $@; done

test rebuild: ;
	-@for subdir in $(DIRS) ; do \
		$(MAKE) --no-print-directory -C $$subdir $@;  done
