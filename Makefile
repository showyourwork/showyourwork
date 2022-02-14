.PHONY: pdf reserve clean arxiv preprocess install_deps conda_setup _update_cache _restore_cache Makefile

# PATHS
HERE            	:= $(realpath $(dir $(realpath $(lastword $(MAKEFILE_LIST)))))
PREPROCESS      	:= $(realpath $(HERE)/workflow/preprocess.smk)
USER 				:= $(realpath $(dir $(HERE)))
CACHE               := $(abspath $(USER)/.showyourwork/cache)

# Default Snakemake options (user can override)
OPTIONS         	?= 

# Pinned package versions
MAMBA_VERSION   	:= 0.17.0
SNAKEMAKE_VERSION 	:= 6.15.5
JINJA2_VERSION      := 3.0.2

# Always enforce these Snakemake options
FORCE_OPTIONS   	:= -c1 --use-conda --reason --cache -d $(USER)

# Skip user prompt in conda install when running on CI
CONDA_YES           := $(shell [ "${CI}" == "true" ] && echo -y)

# Test if things are installed
CONDA_TEST     		:= $(shell command -v conda 2> /dev/null)
SNAKEMAKE_TEST 		:= $(shell command -v snakemake 2> /dev/null)
JINJA2_TEST         := $(shell conda list jinja2 | grep "jinja2 " 2> /dev/null)

# Get installed versions
SNAKEMAKE_INSTALLED_VERSION = $(shell snakemake -v 2> /dev/null)
JINJA2_INSTALLED_VERSION    = $(shell python -c "import jinja2; print(jinja2.__version__)" 2> /dev/null)

# Error handlers
ERROR_HANDLER        = python workflow/utils/scripts/error_handler.py $$?

# Snakemake command
SNAKEMAKE            = SNAKEMAKE_OUTPUT_CACHE=$(CACHE) snakemake $(FORCE_OPTIONS) $(OPTIONS)


# Default target: generate the article
pdf: preprocess
	@$(SNAKEMAKE); $(ERROR_HANDLER)


# Ensure conda is setup
conda_setup:
	@[ "${CONDA_TEST}" ] || ( \
		echo "Conda package manager not found. Please install it from anaconda.com/products/individual."; \
		false \
	)


# Ensure Snakemake & other dependencies are installed
install_deps: conda_setup
	@[ "${SNAKEMAKE_TEST}" ] || ( \
		echo "Snakemake not found. Installing it using conda..."; \
		conda install $(CONDA_YES) -c defaults -c conda-forge -c bioconda mamba==$(MAMBA_VERSION) snakemake-minimal==$(SNAKEMAKE_VERSION) \
	);
	@[ "${SNAKEMAKE_INSTALLED_VERSION}" == "${SNAKEMAKE_VERSION}" ] || ( \
		echo "Snakemake version ${SNAKEMAKE_INSTALLED_VERSION} found, but showyourwork requires version ${SNAKEMAKE_VERSION}. Installing it using conda..."; \
		conda install $(CONDA_YES) -c defaults -c conda-forge -c bioconda mamba==$(MAMBA_VERSION) snakemake-minimal==$(SNAKEMAKE_VERSION) \
	);
	@[ "${JINJA2_TEST}" ] || ( \
		echo "Jinja2 not found. Installing it using conda..."; \
		conda install $(CONDA_YES) -c conda-forge jinja2==$(JINJA2_VERSION) \
	);
	@[ "${JINJA2_INSTALLED_VERSION}" == "${JINJA2_VERSION}" ] || ( \
		echo "Jinja2 version ${JINJA2_INSTALLED_VERSION} found, but showyourwork requires version ${JINJA2_VERSION}. Installing it using conda..."; \
		conda install $(CONDA_YES) -c conda-forge jinja2==$(JINJA2_VERSION) \
	)


# Pre-processing step
preprocess: install_deps
	@$(SNAKEMAKE) -s $(PREPROCESS); $(ERROR_HANDLER)


# Generate the arxiv tarball
arxiv: preprocess
	@$(SNAKEMAKE) syw__arxiv_entrypoint; $(ERROR_HANDLER)


# Clean
clean: install_deps
	@$(SNAKEMAKE) --config debug=true --delete-all-output
	@$(SNAKEMAKE) --config debug=true -s $(PREPROCESS) --delete-all-output
	@rm -rf $(USER)/arxiv.tar.gz
	@rm -rf $(USER)/.showyourwork


# Update the cache on CI (internal command)
_update_cache:
	@cd $(USER) && python showyourwork/workflow/utils/scripts/cache.py --update


# Restore the cache on CI (internal command)
_restore_cache:
	@cd $(USER) && python showyourwork/workflow/utils/scripts/cache.py --restore


# Catch-all target: route all unknown targets to Snakemake
%: Makefile preprocess
	@$(SNAKEMAKE) $@; $(ERROR_HANDLER)