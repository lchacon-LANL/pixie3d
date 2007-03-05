#Define variables

DIRS = $(shell ls | grep -v ake)

#Define targets

.PHONY: test rebuild setup

setup: ;
	-@for subdir in $(DIRS); do \
		ln -s ../make.common $$subdir/makefile;  done

rebuild: ;
	-@for subdir in $(DIRS) ; do \
		$(MAKE) --no-print-directory -C $$subdir rebuild;  done

test: ;
	-@for subdir in $(DIRS) ; do \
		$(MAKE) --no-print-directory -C $$subdir test;  done
