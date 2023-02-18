SYWPLUG__TEX_RESOURCE = partial(
    package_data, "showyourwork.plugins.tex", "workflow"
)

def sywplug__tex_local_or_provided_style(document):
    """Get the path to the showyourwork.sty file. We prefer to use the one
    provided by the user if it exists, but will provide our own if not.
    """
    path = (SYW__REPO_PATHS.root / document).parent / "showyourwork.sty"
    if path.is_file():
        return path
    else:
        return SYWPLUG__TEX_RESOURCE("resources", "showyourwork.sty")
