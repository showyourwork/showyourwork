from pathlib import Path
from typing import Dict, List, Optional, Set


def relative_or_skip(
    input_path: str, repo_path: Path, skip_path: Optional[Path] = None
) -> Optional[str]:
    path = Path(input_path).resolve()
    if skip_path is not None and path.is_relative_to(skip_path):
        return None
    try:
        return str(path.relative_to(repo_path))
    except ValueError:
        return None


def get_including_directory(tree: Dict[str, List[str]], query: str) -> List[str]:
    if query in tree:
        return tree[query]

    # If the query is not in the tree, then it might be a file in a directory
    path = Path(query)
    for key in tree:
        # TODO(dfm): Would we ever expect there to be multiple rules that would
        # satisfy this condition? I think no because the rule order should have
        # already been accounted for and otherwise ambiguous rules would have
        # already errored out.
        if path.is_relative_to(key):
            return tree[key]

    # If we get here we don't have any rules for this file - assume it is static
    return []


def simplify_dependency_tree(
    full_tree: Dict[str, List[str]], repo_path: Path, skip_path: Optional[Path] = None
) -> Dict[str, List[str]]:
    tree: Dict[str, Set[str]] = {}
    for file, parents in full_tree.items():
        path = relative_or_skip(file, repo_path, skip_path)
        if path is None:
            continue
        tree[str(path)] = set()
        for parent in parents:
            todo = [parent]
            while len(todo):
                query = todo.pop()
                qpath = relative_or_skip(query, repo_path, skip_path)
                if qpath is None:
                    todo += get_including_directory(full_tree, query)
                else:
                    tree[str(path)] |= {str(qpath)}
    return {k: list(sorted(v)) for k, v in tree.items()}
