from helpers import TemporaryShowyourworkRepository


class TestDagGeneration(TemporaryShowyourworkRepository):
    """Ensure a local build produces `dag.pdf`."""

    # Keep this local only; remote covered elsewhere
    local_build_only = True

    def customize(self):
        """Enable DAG rendering in the source config."""
        config_file = self.cwd / "showyourwork.yml"
        with open(config_file) as f:
            lines = f.readlines()

        # Find and replace the render: false line in the dag section
        new_lines = []
        for line in lines:
            if "render: false" in line and new_lines and "dag" in new_lines[-1]:
                new_lines.append(line.replace("render: false", "render: true"))
            else:
                new_lines.append(line)

        with open(config_file, "w") as f:
            f.writelines(new_lines)

    def check_build(self):
        dag = self.cwd / "dag.pdf"
        assert dag.exists(), "dag.pdf was not generated"
        assert dag.stat().st_size > 0, "dag.pdf is empty"
