from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
from showyourwork.test_util import run_snakemake


@pytest.mark.parametrize("prefix", ["", "syw__", "sywplug_"])
def test_rule_order(prefix: str) -> None:
    with TemporaryDirectory() as d:
        root = Path(d)
        with open(root / "Snakefile", "w") as f:
            f.write(
                """
from showyourwork.rule_order import fix_rule_order

rule a:
    output: "a.txt"
    shell: "echo a > {output}"

rule %s__a:
    output: "a.txt"
    shell: "echo syw__a > {output}"

fix_rule_order(config, workflow)
"""
                % prefix
            )
        if prefix == "":
            # If no prefix is used then snakemake should fail since the rule
            # order is ambiguous.
            with pytest.raises(RuntimeError):
                run_snakemake(root / "Snakefile", ["a.txt"], cwd=root)
        else:
            run_snakemake(root / "Snakefile", ["a.txt"], cwd=root)
            assert (root / "a.txt").read_text().strip() == "a"
