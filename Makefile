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
	@rm -f src/*.o src/*.so; rm -f scigenex.Rcheck/dbfmcl/libs/dbfmcl.so
	@rm -rf /tmp/dbfmcl; rm -rf *dbf_out.txt; rm -rf *mcl_out.txt 

check: clean
	@rm -rf /tmp/scigenex; mkdir -p /tmp/scigenex; cp -r ./* /tmp/scigenex; cd /tmp/scigenex; \
	find ..Rcheck/ -delete; find ./dbfmcl.Rcheck -delete \
	rm -f src/*.o src/*.so; rm -f scigenex.Rcheck/dbfmcl/libs/dbfmcl.so; \
	R CMD check .

doc:
	@echo ">>> Creating a package documentation"
	@echo "library(roxygen2); roxygen2::roxygenise()" | R --slave

install:
	@echo ">>> Installing..."
	@rm -f src/*.o src/*.so
	@R CMD INSTALL .

test:
	@echo ">>> Testing package"
	@rm -rf `ls tests/testthat/| grep -v \R$$`
	@echo "devtools::test()" | R --slave

coverage:
	@echo "Checking coverage"
	@echo "usethis::use_github_action('test-coverage'); cov <- covr::package_coverage(); print(as.data.frame(cov))" | R --slave

all: doc install check test
