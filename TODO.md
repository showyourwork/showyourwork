To-do
-----

A list of things we should fix/implement before the first major release.

- [ ] Regarding Zenodo deposits with tarballs containing lots of files:
      consider supporting directories (instead of specific files) in the ``contents`` mapping of
      the config file.
- [ ] Don't show errors on make clean
- [ ] Histeresis when adding a new figure to the tex, building, and _then_ adding the script
- [ ] Histeresis when building article, adding a free-floating figure (with no Snakemake rule),
      and building the script -- no error is thrown, and pdf isn't recompiled. Similar behavior
      when running `clean` -- the repo is not cleaned, presumably because of a silent error
      in the DAG generation?
- [ ] Add option to create new Zenodo deposit for a branch
- [ ] Test uploads of very large files
- [ ] Implement showyourwork publish (for Zenodo)
- [ ] Support .zip files on Zenodo