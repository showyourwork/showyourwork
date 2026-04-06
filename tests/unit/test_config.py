import yaml

from showyourwork.config import edit_yaml, edit_yaml_roundtrip

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
