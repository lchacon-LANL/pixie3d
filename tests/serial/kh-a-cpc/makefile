#Define variables

-include make.inc

#Define targets

.PHONY: test rebuild setup

setup: ;
	-@ln -s ../make.common makefile
	-@ln -s ../../../bin/$(EXEC) .
	-@ln -s ../../../bin/$(PLOT) .

rebuild: ;
	@echo "Rebuilding test in $(CURDIR)"
	-@$(EXEC) > pixie3d.out

test: ;
	@echo
	@echo "Running test in $(CURDIR)"
	-@mv record.bin record.bak
	-@$(EXEC) > pixie3d.tst
ifneq (,$(strip $(shell diff pixie3d.tst pixie3d.out 2> /dev/null)))
	@echo "Differences found in $(CURDIR)"
#	@echo $(shell diff pixie3d.tst pixie3d.out 2> /dev/null)
else
	@echo "No differences found in $(CURDIR)"
	-@rm *.bak pixie3d.tst
endif
