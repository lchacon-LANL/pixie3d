#Define variables

LOCALDIR = $(shell basename $(CURDIR))

PAREXEC=./run.test

#Define targets

.PHONY: test test-a test-b rebuild rebuild-a rebuild-b setup clean

setup: ;
	-@ln -s -f ../../../bin/adios_config.xml
	-@ln -s -f ../../../bin/pixeq_xfer.petsc.x

rebuild rebuild-a rebuild-b: ;
	@echo "Rebuilding parallel test in $(LOCALDIR)"
	-@rm -f pixeq_xfer.out* pixeq_xfer.tst* pixeq_xfer.err* record.* 2> /dev/null || true
	@$(PAREXEC) > pixeq_xfer.out

test test-a test-b: ;
	@echo
	@echo "Running parallel test in $(LOCALDIR)"
	@$(PAREXEC) > pixeq_xfer.tst 2> pixeq_xfer.err || true
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] || echo "==> Differences found in $(LOCALDIR) <==") || true
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] || diff -bB  pixeq_xfer.out pixeq_xfer.tst) || true
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] && echo "No differences found in $(LOCALDIR)") || true
	-@([ -z "`diff -bwB pixeq_xfer.out pixeq_xfer.tst`" ] && rm pixeq_xfer.tst 2> /dev/null) || true

run: ;
	@$(PAREXEC)

clean: ;
	-@rm -r -f draw*.in *.bin* pixeq_xfer.tst* pixeq_xfer.err* *.vlog *.bp *.ps *.orig *diag* draw*txt grid*txt *log* *dSYM 2>/dev/null
