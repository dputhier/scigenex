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

clean:
	@rm -f src/*.o src/*.so; rm -f dbfmcl.Rcheck/dbfmcl/libs/dbfmcl.so; rm -rf ./dbfmcl.Rcheck
	@rm -rf /tmp/dbfmcl; 

check: clean
	@mkdir -p /tmp/dbfmcl; cp -r ./* /tmp/dbfmcl
	@R CMD check /tmp/dbfmcl

doc:
	@echo "library(roxygen2); roxygen2::roxygenise()" | R --slave

install:
	@rm -f src/*.o src/*.so
	@R CMD Install .

test:
	@echo "devtools::test()" | R --slave

all: doc install check test
