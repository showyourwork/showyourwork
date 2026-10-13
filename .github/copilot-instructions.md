# showyourwork contributor instructions

## Setup, checks, and documentation

- Use Python 3.11+ and install the editable development environment with:
  ```bash
  python -m pip install -e ".[dev]"
  ```
- Run the unit suite with `python -m pytest tests/unit`. Run one unit test with
  `python -m pytest tests/unit/test_config.py::test_edit_yaml_preserves_data`.
- Run the local integration suite with `python -m pytest tests/integration`.
  Run one integration test with
  `python -m pytest tests/integration/test_default.py::TestDefault::test_local`.
  These tests create projects below `tests/integration/sandbox`; set `DEBUG=1`
  and pass `-s` to expose the generated project's command output while debugging.
- Remote-marked integration tests are skipped unless `--remote` is supplied.
  They create and push temporary GitHub repositories. To run them against a
  personal account, set `GH_API_KEY` and use:
  ```bash
  python -m pytest tests/integration --remote --no-org \
    --action-spec git+https://github.com/showyourwork/showyourwork.git
  ```
  Zenodo tests are skipped unless both `SANDBOX_TOKEN` and `ZENODO_TOKEN` are
  available; use `--require-zenodo` only when those credentials are configured.
- Format and lint Python with `ruff format .` and `ruff check .`; `pre-commit
  run --all-files` also applies Ruff and formats Snakefiles with `snakefmt`.
- Build the Sphinx documentation with `cd docs && make html` after installing
  the `docs` extra (or the `dev` extra).

## Architecture

- `src/showyourwork/cli/` is the Click CLI. Public commands delegate to
  `cli/commands/`, which sets the `SNAKEMAKE_RUN_TYPE` context and invokes the
  packaged Snakemake workflows.
- A build always has two stages. `workflow/prep.smk` renders the project's
  Jinja-capable `showyourwork.yml`, parses the manuscript into a workflow graph,
  and writes `.showyourwork/config.json`. `workflow/build.smk` consumes that
  generated configuration, includes the project's root `Snakefile`, executes
  the graph, and produces `<ms_name>.pdf`.
- Built-in Snakemake rules live in `workflow/rules/`; scripts and TeX resources
  supporting those rules are packaged alongside them. Keep a workflow rule,
  its script/resource dependency, and the corresponding configuration handling
  consistent when changing behavior.
- `cookiecutter-showyourwork/` is the source for repositories created by
  `showyourwork setup`. Changes to its workflows or `showyourwork.yml` should
  be covered by `tests/unit/test_workflow_templates.py` and, where relevant,
  an integration test that customizes a generated project.
- Integration tests subclass `TemporaryShowyourworkRepository`. They create a
  fresh project through the public `showyourwork setup` command, customize it,
  build it locally, and optionally push it to GitHub Actions. Prefer extending
  that harness over constructing ad hoc fixture repositories.

## Workflow conventions

- The generated project's root `Snakefile` is the supported extension point.
  User-defined rules are deliberately ordered ahead of built-in `syw__*` rules.
  Do not use Snakemake `run:` directives in user rules: the package rejects
  them to preserve isolated, reproducible execution. Use `script:` or `shell:`
  and declare the rule's environment.
- Keep paths relative to an article repository root. The standard project
  layout is `src/scripts` for generators, `src/data` for data,
  `src/tex/figures` for generated figures, and `src/tex/output` for generated
  manuscript inputs; `paths.user()` defines these locations centrally.
- The preprocessing stage discovers figure provenance from TeX. A figure
  produced by a script should use one `\script{...}` declaration and a label in
  its figure environment. Generated text included by the manuscript should use
  `\variable{output/...}` and have a Snakefile rule that produces the file.
  Static, non-generated figures belong in `src/static`.
- Configuration is merged and normalized in `config.py`; preserve the
  preprocessing/build-stage distinction when adding configuration because
  `.showyourwork/config.json` is the handoff between stages.
- A user-facing change requires a Towncrier fragment in `docs/changes/` named
  `<pull_request_number>.<type>.rst`, where type is `bugfix`, `feature`,
  `maintenance`, `api`, `optimization`, or `documentation`.
