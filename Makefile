# Default Snakemake options. Change here or override by setting
# an environment variable as needed.
OPTIONS         ?= -c1


# Local vars
FORCE_OPTIONS   := -c1 --use-conda
TEMPORARIES     := .showyourwork src/ms.pdf src/__latexindent*.tex
CONDA           := $(shell conda -V 2&> /dev/null && echo 1 || echo 0)
JINJA2          := $(shell conda list jinja2 | grep "jinja2 " && echo 1 || echo 0)
SNAKEMAKE       := $(shell snakemake -v 2&> /dev/null && echo 1 || echo 0)
WORKDIR         := ..
CLEAN_SYW       := rm -rf $(TEMPORARIES)
CLEAN_SM        := snakemake $(OPTIONS) $(FORCE_OPTIONS) ms.pdf --delete-all-output

.PHONY:  ms.pdf clean report reserve lint fast version update snakemake_setup conda_setup jinja2_setup Makefile


# Default target: generate the article
ms.pdf: snakemake_setup
	@cd $(WORKDIR);\
	snakemake $(FORCE_OPTIONS) $(OPTIONS) ms.pdf


# Ensure conda is setup
conda_setup:
	@if [ "$(CONDA)" = "0" ]; then \
		echo "Conda package manager not found. Please install it from anaconda.com/products/individual.";\
		false;\
	fi


# Ensure jinja2 is installed
jinja2_setup: conda_setup
	@if [ "$(JINJA2)" = "0" ]; then \
		echo "Jinja2 not found. Installing it using conda...";\
		conda install -c conda-forge jinja2==2.11.3; \
	fi


# Ensure Snakemake is setup
snakemake_setup: conda_setup jinja2_setup
	@if [ "$(SNAKEMAKE)" = "0" ]; then \
		echo "Snakemake not found. Installing it using conda...";\
		conda install -c defaults -c conda-forge -c bioconda mamba==0.17.0 snakemake-minimal==6.12.3; \
	fi
	

# Remove all intermediates, outputs, and temporaries
clean: snakemake_setup
	@cd $(WORKDIR);\
	{ $(CLEAN_SM) && $(CLEAN_SYW); } || { $(CLEAN_SYW) && $(CLEAN_SM); }
	

# Generate a workflow report
report: snakemake_setup
	@cd $(WORKDIR);\
	snakemake $(FORCE_OPTIONS) $(OPTIONS) ms.pdf --report


# Update to the latest version of showyourwork
update: snakemake_setup
	@git fetch --all --tags
	git checkout $(shell git describe --tags `git rev-list --tags --max-count=1`)


# Pre-reserve a Zenodo DOI
reserve: snakemake_setup
	@cd workflow; \
	python -c "import helpers; helpers.zenodo.reserve()"


# Lint the user's repository
lint: snakemake_setup
	@cd workflow; \
	python -c "import helpers; helpers.linter.lint()"


# Fast build (never re-generate dependencies available on Zenodo)
fast: snakemake_setup
	@cd $(WORKDIR);\
	snakemake $(FORCE_OPTIONS) $(OPTIONS) --config="download_only=true" ms.pdf


# Print the showyourwork version
version: snakemake_setup
	@python -c "import re; print(re.search('tag/(.*?)\"', open('README.md', 'r').read()).groups()[0])"


# Catch-all target: route all unknown targets to Snakemake
%: Makefile snakemake_setup
	@cd $(WORKDIR);\
	snakemake $(FORCE_OPTIONS) $(OPTIONS) $@