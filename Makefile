MAKEFILE=Makefile

.PHONY: help

#------------------------------------------------------------------
# Targets
#------------------------------------------------------------------


help:
	@echo ""
	@echo "- Available targets:"
	@echo ""
	@perl -ne 'if(	/^(\w+):/){print "\t",$$1,"\n"}' $(MAKEFILE)
	@echo ""
	@echo ""

check:
	@rm -f src/*.o src/*.so; rm -f dbfmcl.Rcheck/dbfmcl/libs/dbfmcl.so; rm -rf ./dbfmcl.Rcheck
	@rm -rf /tmp/dbfmcl; mkdir -p /tmp/dbfmcl; cp -r ./* /tmp/dbfmcl
	@R CMD check /tmp/dbfmcl

install:
	@rm -f src/*.o src/*.so
	@echo "library(roxygen2); roxygen2::roxygenise()" | R --slave
	@R CMD Install .

test:
	@echo "devtools::test()" | R --slave
