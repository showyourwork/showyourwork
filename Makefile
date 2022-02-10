.PHONY: pdf reserve clean preprocess snakemake_setup conda_setup _update_cache _restore_cache Makefile

# PATHS
HERE            	:= $(realpath $(dir $(realpath $(lastword $(MAKEFILE_LIST)))))
PREPROCESS      	:= $(realpath $(HERE)/workflow/preprocess.smk)
USER 				:= $(realpath $(dir $(HERE)))
CACHE               := $(abspath $(USER)/.showyourwork/cache)

# Default Snakemake options (user can override)
OPTIONS         	?= 

# Pinned package versions
MAMBA_VERSION   	:= 0.17.0
SNAKEMAKE_VERSION 	:= 6.12.3

# Always enforce these Snakemake options
FORCE_OPTIONS   	:= -c1 --use-conda --cache -d $(USER)

# Boolean flag: is conda installed?
CONDA_INSTALLED     := $(shell conda -V 2&> /dev/null && echo 1 || echo 0)

# Boolean flag: is snakemake installed?
SNAKEMAKE_INSTALLED := $(shell snakemake -v 2&> /dev/null && echo 1 || echo 0)

# Error handlers
ERROR_HANDLER        = python workflow/utils/scripts/error_handler.py $$?

# Snakemake command
SNAKEMAKE            = SNAKEMAKE_OUTPUT_CACHE=$(CACHE) snakemake $(FORCE_OPTIONS) $(OPTIONS)


# Default target: generate the article
pdf: preprocess
	@$(SNAKEMAKE); $(ERROR_HANDLER)


debug:
	echo $(SNAKEMAKE_INSTALLED)



# Ensure conda is setup
conda_setup:
	@if [ "$(CONDA_INSTALLED)" = "0" ]; then \
		echo "Conda package manager not found. Please install it from anaconda.com/products/individual."; \
		false; \
	fi


# Ensure Snakemake is setup
snakemake_setup: conda_setup
	@if [ "$(SNAKEMAKE_INSTALLED)" = "0" ]; then \
		echo "Snakemake not found. Installing it using conda..."; \
		conda install -c defaults -c conda-forge -c bioconda mamba==$(MAMBA_VERSION) snakemake-minimal==$(SNAKEMAKE_VERSION); \
	fi


# Pre-processing step
preprocess: snakemake_setup
	@$(SNAKEMAKE) -s $(PREPROCESS); $(ERROR_HANDLER)


# Clean
clean: snakemake_setup
	@$(SNAKEMAKE) --config debug=true --delete-all-output
	@$(SNAKEMAKE) --config debug=true -s $(PREPROCESS) --delete-all-output
	@rm -rf $(USER)/.showyourwork


# Update the cache on CI (internal command)
_update_cache:
	@python workflow/utils/scripts/cache.py --update


# Restore the cache on CI (internal command)
_restore_cache:
	@python workflow/utils/scripts/cache.py --restore


# Catch-all target: route all unknown targets to Snakemake
%: Makefile preprocess
	@$(SNAKEMAKE) $@; $(ERROR_HANDLER)