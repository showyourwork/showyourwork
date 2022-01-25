- [ ] If the user changes, say, the ID of a Zenodo deposit that provides a certain
      file, and that file already exists locally, the workflow currently does NOT
      re-download it. Think about how to make this happen automatically.

- [ ] Use the following `yaml` to test our Zenodo interface:

      ```yaml
      dependencies:
          src/figure-scripts/sample_figure.py: 
              - src/data/A.dat
              - src/data/baz/B.dat
              - src/data/foo/C.dat
              - src/data/foo/bar/H.dat
              - src/data/G.dat
              - src/data/C.dat
              - src/data/bar/H.dat
              - src/data/G_copy.dat

      zenodo_sandbox:
          976786:
              contents:
                  A.dat:                                    # auto mapping to src/data/A.dat
                  B.dat: src/data/baz/B.dat                 # explicit mapping to local subfolder (create if needed)
                  foo.tar.gz:                               # remote tarfiles behave like folders
                      foo:                                  # files are nested inside `foo` in this tarball
                          C.dat:                            # auto mapping to src/data/foo/C.dat
                          bar:                              # subfolder
                              H.dat:                        # auto mapping to src/data/foo/bar/H.dat
                              G.dat: src/data/G.dat         # change location
                  foo2.tar.gz:                              # same as `foo.tar.gz`, but no top-level `foo` folder
                      C.dat:                                # auto mapping to src/data/C.dat
                      bar:                                  # subfolder
                          H.dat:                            # auto mapping to src/data/bar/H.dat
                          G.dat: src/data/G_copy.dat        # rename
      ```