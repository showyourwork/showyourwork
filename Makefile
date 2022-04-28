# Makefile use is deprecated as of v1.0.0. If the user tries to build an
# article usign `make`, this displays a warning telling them to migrate the
# repository to the new showyourwork standard.

.PHONY: default

default:
	@printf "\e[31mWARNING: showyourwork no longer uses Makefiles to build articles.\e[0m\n"
	@printf "\e[31mIf you recently updated your installation from v0.2.3 or earlier, please visit\e[0m\n"
	@printf "\n"
	@printf "      showyourwork.readthedocs.io/en/latest/migrating\n"
	@printf "\n"
	@printf "\e[31mfor instructions on how to migrate your article to the new showyourwork.\e[0m\n"

%:
	@printf "\e[31mWARNING: showyourwork no longer uses Makefiles to build articles.\e[0m\n"
	@printf "\e[31mIf you recently updated your installation from v0.2.3 or earlier, please visit\e[0m\n"
	@printf "\n"
	@printf "      showyourwork.readthedocs.io/en/latest/migrating\n"
	@printf "\n"
	@printf "\e[31mfor instructions on how to migrate your article to the new showyourwork.\e[0m\n"