"""
This Snakefile contains default rules for copying the documents from the project
root to the work directory. These can be overloaded by users or plugins to
generate a working copy of the documents differently.
"""

for doc in SYW__DOCUMENTS:
    rule:
        """
        Copy a document from the project root to the work directory.
        """
        name:
            utils.rule_name("copy", "doc", document=doc)
        input:
            SYW__REPO_PATHS.root / doc
        output:
            SYW__WORK_PATHS.root / doc
        run:
            utils.copy_file_or_directory(input[0], output[0])
