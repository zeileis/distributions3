

# Extracting current package version
VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

.PHONY: document
document:
	Rscript -e "devtools::document()"

.PHONY: clean
clean:
	-rm src/*.so
	-rm src/*.o
	-rm vignettes/*.html
	-rm -r vignettes/*_files/

.PHONY: install build
build: clean document
	@echo Building current version: $(VERSION)
	(cd ../ && R CMD build distributions3)

install: build
	@echo Installing current version: $(VERSION)
	(cd ../ && R CMD INSTALL distributions3_$(VERSION).tar.gz)

check: build
	@echo "Checking current version (w/o --as-cran): $(VERSION)"
	(cd ../ && R CMD check distributions3_$(VERSION).tar.gz)
	##@echo Checking current version: $(VERSION)
	##(cd ../ && R CMD check --as-cran distributions3_$(VERSION).tar.gz)

.PHONY: test
test: clean install
	Rscript -e "testthat::test_local()"

.PHONY: coverage
coverage: install
	Rscript -e "covr::report(covr::package_coverage(line_exclusions = list('src/init.c'), function_exclusions = list('message\\\\s*\\\\(', 'plot.distributions3')), file = \"_coverage.html\")"

