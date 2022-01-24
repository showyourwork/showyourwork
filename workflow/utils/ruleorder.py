import snakemake

__all__ = ["process_user_rules"]


def process_user_rules(user_snakefile):
    """
    Give all user rules higher priority than showyourwork rules
    and add a message if they don't have one

    """
    syw_rules = []
    usr_rules = []
    for r in snakemake.workflow.workflow.rules:
        if r.snakefile == user_snakefile:
            usr_rules.append(r)
        else:
            syw_rules.append(r)
    for ur in usr_rules:
        for sr in syw_rules:
            snakemake.workflow.workflow.ruleorder(ur.name, sr.name)
        if not ur.message:
            ur.message = f"Running user rule {ur.name}..."