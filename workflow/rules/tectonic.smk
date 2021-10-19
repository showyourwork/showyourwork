rule tectonic:
    """
    Install the latest version of tectonic.
    
    """
    output:
        posix(relpaths.temp / "tectonic")
    params:
        TEMP=relpaths.temp,
        OS=config["tectonic_os"]
    script:
        "../scripts/tectonic.py"