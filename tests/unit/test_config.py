import yaml

from showyourwork.config import edit_yaml, edit_yaml_roundtrip, expand_dependency_directories

START_YAML = """# Top-level comment
title: "Example"
enabled: true
items:
  - one
  - two
nested:
  # Nested comment
  value: 42
"""


def test_edit_yaml_preserves_data(tmp_path):
    config_file = tmp_path / "config.yml"
    config_file.write_text(START_YAML)

    with edit_yaml(config_file) as config:
        config["nested"]["value"] = 43

    before = yaml.safe_load(START_YAML)
    after = yaml.safe_load(config_file.read_text())

    assert after["title"] == before["title"]
    assert after["enabled"] == before["enabled"]
    assert after["items"] == before["items"]
    assert after["nested"]["value"] == 43


def test_edit_yaml_roundtrip_preserves_entire_file(tmp_path):
    config_file = tmp_path / "config.yml"
    config_file.write_text(START_YAML)

    with edit_yaml_roundtrip(config_file) as config:
        # Dummy change
        config["nested"]["value"] = 43
        config["nested"]["value"] = 42

    assert config_file.read_text() == START_YAML


def test_expand_dependency_directories_skips_external_paths(tmp_path):
    """expand_dependency_directories must not raise ValueError for files
    outside the repo root (e.g. symlinks to system paths like /boot/System.map-*)."""
    # Build a fake repo root with a dependency directory
    repo_root = tmp_path / "repo"
    dep_dir = repo_root / "src" / "data" / "inputs"
    dep_dir.mkdir(parents=True)

    # An in-repo file that should be included
    in_repo_file = dep_dir / "data.csv"
    in_repo_file.write_text("x,y\n1,2\n")

    # Simulate an external path by creating a symlink inside the dep dir that
    # points to a file outside the repo root (mimics /boot/System.map-* on CI)
    external_dir = tmp_path / "external"
    external_dir.mkdir()
    external_file = external_dir / "System.map-6.17.0-1020-azure"
    external_file.write_text("kernel map")
    external_link = dep_dir / "external_link"
    external_link.symlink_to(external_dir)

    result = expand_dependency_directories("src/data/inputs", repo_root=repo_root)

    # Only the in-repo file should be returned; external paths must be skipped
    assert result == ["src/data/inputs/data.csv"]


def test_expand_dependency_directories_in_repo_files_unchanged(tmp_path):
    """expand_dependency_directories returns all in-repo files unchanged."""
    repo_root = tmp_path / "repo"
    dep_dir = repo_root / "src" / "data"
    dep_dir.mkdir(parents=True)

    (dep_dir / "a.txt").write_text("a")
    (dep_dir / "b.txt").write_text("b")
    sub = dep_dir / "sub"
    sub.mkdir()
    (sub / "c.txt").write_text("c")

    result = expand_dependency_directories("src/data", repo_root=repo_root)

    assert sorted(result) == ["src/data/a.txt", "src/data/b.txt", "src/data/sub/c.txt"]
