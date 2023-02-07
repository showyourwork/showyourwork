from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import snakemake


def fix_rule_order(workflow: "snakemake.workflow.Workflow") -> None:
    """Update the rule order for all rules defined by the workflow

    The logic here is that we want all user-defined rules to be higher priority
    than showyourwork (and plugin) rules, and we want all plugin rules to be
    higher priority than the base showyourwork rules. Whether or not a rule is a
    showyourwork or plugin rule is determined by the prefix of the rule name,
    ``syw__`` and ``sywplug__`` respectively.
    """
    syw_rules = []
    sywplug_rules = []
    user_rules = []
    for r in workflow.rules:
        if r.name.startswith("syw__"):
            syw_rules.append(r)
        elif r.name.startswith("sywplug__"):
            sywplug_rules.append(r)
        else:
            user_rules.append(r)

    # All user-defined rules should be higher priority than any showyourwork or
    # plugin rules.
    for ur in user_rules:
        for sr in syw_rules:
            workflow.ruleorder(ur.name, sr.name)
        for sr in sywplug_rules:
            workflow.ruleorder(ur.name, sr.name)

        if not ur.message:
            ur.message = f"Running user rule {ur.name}..."

    # All plugin rules should be higher priority than any showyourwork rules.
    for pr in sywplug_rules:
        for sr in syw_rules:
            workflow.ruleorder(pr.name, sr.name)
