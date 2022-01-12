.PHONY: pdf workflow_setup snakemake_setup conda_setup Makefile


# PATHS
HERE            := $(realpath $(dir $(realpath $(lastword $(MAKEFILE_LIST)))))
USER 			:= $(realpath $(dir $(HERE)))
SETUPFILE       := $(realpath $(HERE)/workflow/setup.smk)
BUILDFILE       := $(realpath $(HERE)/workflow/build.smk)

# Default Snakemake options (user can override)
OPTIONS         ?=

# Always enforce these Snakemake options
FORCE_OPTIONS   := -c1 --use-conda -d $(USER)

# Boolean flag: is conda installed?
CONDA           := $(shell conda -V 2&> /dev/null && echo 1 || echo 0)

# Boolean flag: is snakemake installed?
SNAKEMAKE       := $(shell snakemake -v 2&> /dev/null && echo 1 || echo 0)


# Default target: generate the article
pdf: workflow_setup
	@snakemake $(FORCE_OPTIONS) $(OPTIONS) -s $(BUILDFILE) pdf


# Ensure conda is setup
conda_setup:
	@if [ "$(CONDA)" = "0" ]; then \
		echo "Conda package manager not found. Please install it from anaconda.com/products/individual.";\
		false;\
	fi


# Ensure Snakemake is setup
snakemake_setup: conda_setup
	@if [ "$(SNAKEMAKE)" = "0" ]; then \
		echo "Snakemake not found. Installing it using conda...";\
		conda install -c defaults -c conda-forge -c bioconda mamba==0.17.0 snakemake-minimal==6.12.3; \
	fi


# Setup stage of the workflow
workflow_setup: snakemake_setup
	@snakemake $(FORCE_OPTIONS) $(OPTIONS) -s $(SETUPFILE) setup


# Catch-all target: route all unknown targets to Snakemake
%: Makefile workflow_setup
	@snakemake $(FORCE_OPTIONS) $(OPTIONS) -s $(BUILDFILE) $@