#Define variables

-include make.inc

LOCALDIR := $(shell basename $(CURDIR))

SEREXEC = ./run.test

#Define targets

.PHONY: test test-a test-b rebuild rebuild-a rebuild-b setup

setup: ;
	-@ln -s -f ../../bin/pixeq_xfer.x

rebuild-a rebuild-b rebuild: ;
	@echo "Rebuilding test in $(LOCALDIR)"
	-@rm pixie3d.out pixie3d.tst record.* 2> /dev/null || true
	@$(SEREXEC) > pixie3d.out

test-a test-b test: ;
	@echo
	@echo "Running test in $(LOCALDIR)"
	-@mv record.bin_proc0 record.bak 2> /dev/null  || true
	@$(SEREXEC) > pixeq_xfer.tst
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] || echo "==> Differences found in $(LOCALDIR) <==")  || true
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] || diff -bB pixeq_xfer.out pixeq_xfer.tst) || true
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] && echo "No differences found in $(LOCALDIR)") || true
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] && rm record.bak pixeq_xfer.tst 2> /dev/null) || true

run: ;
	@$(SEREXEC)

clean: ;
	-@rm -r -f draw*.in *.bin* pixeq_xfer.tst* pixeq_xfer.err* *.vlog *.bp *.ps *.orig *diag* draw*txt grid*txt *log* *dSYM 2>/dev/null
