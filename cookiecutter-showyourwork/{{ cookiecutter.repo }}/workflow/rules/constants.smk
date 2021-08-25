from pathlib import Path, PurePosixPath


# Convert a path object to a string containing the POSIX path
POSIX = lambda path: str(PurePosixPath(path))


# TODO: Better way of getting abs paths?
WORKFLOW = "https://github.com/rodluger/showyourwork/raw/restructure/cookiecutter-showyourwork/%7B%7B%20cookiecutter.repo%20%7D%7D/workflow/" # DEBUG Path(workflow.modules["showyourwork"].snakefile).absolute().parents[0]
USER = Path(workflow.basedir)


# Relative paths (to top level of user's repo)
TEX = Path("tex")
FIGURES = Path("figures")
TEMP = Path(".showyourwork")
GITHUB = Path(".github")


# Dummy file dependency for figures w/ unknown parent scripts
UNKNOWN_SCRIPT = "unknown-script"


# Temporary tex files
TMPTEXFILE = ".showyourwork-xml-ms"
SYWTEXFILE = ".showyourwork-ms"



# DEBUG
AUXFILES = ['tex/xstring.sty', 'tex/sywxml.sty', 'tex/xstring.tex', 'tex/showyourwork.pdf']
CLASSFILES = {'aastex631': ['orcid-ID.png', 'aasjournal.bst', 'revtex4-1.cls', 'aastex631.cls']}


'''
DEBUG

# Auxiliary files we copy over to the user's `tex/` directory
AUXFILES = [POSIX(TEX / file.name) for file in (WORKFLOW / "resources" / "tex").glob("*.*")]

# Class-specific auxiliary files
CLASSFILES = {}
for folder in (WORKFLOW / "resources" / "classes").glob("*"):
    CLASSFILES[folder.name] = [file.name for file in (WORKFLOW / "resources" / "classes" / folder).glob("*.*")]
'''
