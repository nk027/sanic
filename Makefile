PKG       := sanic
VERSION   := $(shell sed -n 's/^Version: *//p' DESCRIPTION)
TARBALL   := $(PKG)_$(VERSION).tar.gz
R         := Rscript
RUN_TESTS := RUN_TESTS=true
CHECK_ENV := _R_CHECK_FORCE_SUGGESTS_=false
# Vignettes render to HTML here, so no LaTeX is required.
VIGNETTES ?=

# Load local tokens (e.g. GOCOV_TOKEN) from .env when present.
ifneq (,$(wildcard .env))
include .env
export GOCOV_TOKEN
endif

.PHONY: all document install build check check-cran test coverage upload clean

all: check coverage

document:
	$(R) -e 'roxygen2::roxygenise()'

install: document
	R CMD INSTALL --no-multiarch --install-tests .

build: document
	R CMD build $(VIGNETTES) .

check: build
	$(RUN_TESTS) $(CHECK_ENV) R CMD check --no-manual $(VIGNETTES) $(TARBALL)

# Full submission check; requires LaTeX for the PDF reference manual.
check-cran: build
	$(RUN_TESTS) R CMD check --as-cran --run-donttest $(TARBALL)

test: install
	$(RUN_TESTS) $(R) -e 'tinytest::test_package("$(PKG)", testdir = "tests/tinytest", at_home = TRUE)'

coverage:
	@mkdir -p scratch
	$(RUN_TESTS) $(R) -e 'cov <- covr::package_coverage(); print(cov); covr::to_cobertura(cov, filename = "scratch/coverage.xml")'
upload:
	go run github.com/gocov/gocov/cmd/gocov@latest upload -format cobertura scratch/coverage.xml
clean:
	rm -rf $(PKG)_*.tar.gz $(PKG).Rcheck ..Rcheck $(PKG)_$(VERSION) scratch/coverage.xml
