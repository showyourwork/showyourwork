rule docstrings:
    """
    Generate simple documentation for the rules in this workflow.

    """
    run:
        with open(abspaths.workflow.parents[0] / "docs" / "rules.rst", "w") as f:
            for rule in workflow.rules:
                header = "rule {}".format(rule.name)
                print(header, file=f)
                print("^" * len(header), file=f)
                if rule.docstring:
                    print(rule.docstring.strip(), file=f)
                print("\n", file=f)