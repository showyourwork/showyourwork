import inspect
from showyourwork.dependencies import simplify_dependency_tree

def get_document_dependencies(doc):
    checkpoint_name = utils.rule_name(
        "check", "dependencies", document=doc
    )

    def impl(*_):
        getattr(checkpoints, checkpoint_name).get()

        # Include the dependencies extracted from the document.
        with open(SYW__WORK_PATHS.dependencies_for(doc), "r") as f:
            dependencies = json.load(f)
        files = list(dependencies.get("unlabeled", [])) + list(dependencies.get("files", []))
        for figure in dependencies.get("figures", {}).values():
            files.extend(figure)

        # And also any files that are explicitly listed in the document's
        # dependencies in the cofing file.
        files.extend(SYW__DOCUMENTS[doc])

        # Save the document dependencies to the "config" object for downstream
        # usage.
        if "_document_dependencies" not in config:
            config["_document_dependencies"] = {}
        config["_document_dependencies"][doc] = files

        return files
    return impl

def ensure_all_document_dependencies(*_):
    from collections import defaultdict

    # This checkpoint call serves two purposes: (1) it makes sure that we have
    # extracted the list of all dependencies from the document, and (2) it
    # ensures that the DAG of jobs has been constructed.
    for doc in SYW__DOCUMENTS:
        checkpoint_name = utils.rule_name(
            "check", "dependencies", document=doc
        )
        getattr(checkpoints, checkpoint_name).get()

    # Walk up the call stack to find an object called "dag"... yeah, this is a
    # hack, but we haven't found a better approach yet!
    dag = None
    for level in inspect.stack():
        dag = level.frame.f_locals.get("dag", None)
        if dag is not None:
            break

    # If "dag" is still None, then we couldn't find it. We shouldn't ever hit
    # this (until snakemake renames the variable...), but we have a check just
    # to be sure.
    if dag is None:
        raise RuntimeError(
            "Could not find DAG object in call stack. This error shouldn't "
            "ever be hit, but you found it! Please report the issue on the "
            "showyourwork GitHub page."
        )

    # Map the full tree of data dependencies. Here we're collecting the
    # "parents" for every file that has a rule defined.
    parents = defaultdict(set)
    for job in dag.jobs:
        for output in job.output:
            parents[output] |= set(str(f) for f in job.input)
    parents = {k: list(sorted(v)) for k, v in parents.items()}

    # Store the dependency tree in the global "config" object. This is also a
    # bit of a hack, but this gives us a nice way to provide downstream rules
    # with access to the computed dependency tree.
    config["_dependency_tree"] = parents
    config["_dependency_tree_simple"] = simplify_dependency_tree(
        parents, SYW__REPO_PATHS.root, SYW__WORK_PATHS.root
    )

    return []

for doc in SYW__DOCUMENTS:
    name = paths.path_to_rule_name(doc)
    checkpoint:
        name:
            utils.rule_name("check", "dependencies", document=doc)
        input:
            SYW__WORK_PATHS.dependencies_for(doc)
        output:
            touch(SYW__WORK_PATHS.flag(f"dependencies_{name}"))

rule:
    name:
        utils.rule_name("dag")
    input:
        [get_document_dependencies(doc) for doc in SYW__DOCUMENTS]
    output:
        touch(SYW__WORK_PATHS.flag("dag"))

rule:
    name:
        utils.rule_name("save", "dependencies")
    message:
        "Saving document dependency tree as JSON"
    input:
        rules.syw__dag.output,
        ensure_all_document_dependencies
    output:
        SYW__WORK_PATHS.root / "dependency_tree.json",
        SYW__WORK_PATHS.root / "dependency_tree_simple.json"
    run:
        import json
        with open(output[0], "w") as f:
            json.dump(config["_dependency_tree"], f, indent=2)
        with open(output[1], "w") as f:
            json.dump(config["_dependency_tree_simple"], f, indent=2)
