rule tectonic:
    """
    Installs the latest version of tectonic. This rule isn't
    used unless the user explicitly asks for ``tectonic_latest``
    in the config file. We originally implemented this rule when
    issues with ``archive.org`` caused the released versions of
    ``tectonic`` to fail, but an unreleased patch was functioning.
    This rule is probably of very limited use in general.
    
    """
    output:
        posix(relpaths.temp / "tectonic")
    params:
        TEMP=relpaths.temp,
        OS=config["tectonic_os"]
    script:
        "../scripts/tectonic.py"