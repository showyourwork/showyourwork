rule tectonic:
    output:
        POSIX(TEMP / "tectonic")
    params:
        TEMP=TEMP,
        OS=tectonic_os
    script:
        "../scripts/tectonic.py"