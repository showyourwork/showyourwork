import os
import shutil
import tempfile
from unittest.mock import patch

import yaml
from cookiecutter.main import cookiecutter

from showyourwork.cli.commands.setup import setup


def render_and_validate_workflow(
    template_dir,
    workflow_name,
    context,
    expected_action_version,
    action_path="showyourwork/showyourwork-action",
):
    temp_dir = tempfile.mkdtemp()
    try:
        cookiecutter(
            template_dir,
            no_input=True,
            extra_context=context,
            output_dir=temp_dir,
        )
        repo_name = context["repo"]
        workflow_path = os.path.join(
            temp_dir,
            repo_name,
            ".github",
            "workflows",
            workflow_name,
        )
        assert os.path.exists(workflow_path), f"Workflow {workflow_name} not found!"
        with open(workflow_path) as f:
            workflow_yaml = yaml.safe_load(f)
        # Check that the correct action version is injected
        found = False
        for job in workflow_yaml.get("jobs", {}).values():
            for step in job.get("steps", []):
                uses = step.get("uses", "")
                if uses.startswith(action_path):
                    assert uses.endswith(f"@{expected_action_version}"), (
                        f"Expected action version @{expected_action_version} "
                        f"in uses: {uses}"
                    )
                    found = True
        assert found, (
            f"No step found using {action_path}@{expected_action_version} "
            f"in {workflow_name}"
        )
    finally:
        shutil.rmtree(temp_dir)


class BaseSetupTest:
    """Base class for setup command tests with common functionality."""

    def setup_method(self):
        """Set up a temporary directory for each test."""
        self.temp_dir = tempfile.mkdtemp()
        try:
            self.original_cwd = os.getcwd()
        except FileNotFoundError:
            # If current directory doesn't exist, use the project root
            self.original_cwd = os.path.dirname(
                os.path.dirname(os.path.dirname(__file__))
            )

    def teardown_method(self):
        """Clean up after each test."""
        os.chdir(self.original_cwd)
        if hasattr(self, "temp_dir"):
            shutil.rmtree(self.temp_dir)

    def run_setup_with_git_mocks(self, git_responses, version="v0.4.3.dev1+test"):
        """
        Run setup command with mocked git responses.

        Args:
            git_responses: Dict mapping git commands to their outputs
            version: Version string to mock

        Returns:
            Tuple of (mock_logger_instance, repo_path, workflow_content)
        """
        with patch("showyourwork.cli.commands.setup.requests"):
            with patch("showyourwork.cli.commands.setup.__version__", version):
                with patch(
                    "showyourwork.cli.commands.setup.get_stdout"
                ) as mock_get_stdout:
                    with patch(
                        "showyourwork.cli.commands.setup.get_logger"
                    ) as mock_logger:

                        def mock_git_commands(cmd, **kwargs):
                            cmd = cmd.strip()
                            return git_responses.get(cmd, "")

                        mock_get_stdout.side_effect = mock_git_commands
                        mock_logger_instance = mock_logger.return_value

                        os.chdir(self.temp_dir)

                        # Run setup command
                        setup(
                            slug="testuser/testrepo",
                            cache=False,
                            overleaf_id=None,
                            ssh=False,
                            action_spec=None,  # Auto-detect
                            action_version="v1",
                        )

                        # Read generated workflow
                        repo_path = os.path.join(self.temp_dir, "testrepo")
                        workflow_path = os.path.join(
                            repo_path, ".github", "workflows", "build.yml"
                        )
                        with open(workflow_path) as f:
                            workflow_content = f.read()

                        return mock_logger_instance, repo_path, workflow_content

    def assert_warning_logged(self, mock_logger_instance, expected_text):
        """Assert that a warning containing the expected text was logged."""
        mock_logger_instance.warning.assert_called()
        warning_calls = [
            call.args[0] for call in mock_logger_instance.warning.call_args_list
        ]
        assert any(
            expected_text in msg for msg in warning_calls
        ), f"Expected warning containing '{expected_text}', got: {warning_calls}"

    def assert_action_spec_in_workflow(self, workflow_content, expected_spec):
        """Assert that the workflow contains the expected action spec."""
        assert (
            expected_spec in workflow_content
        ), f"Expected '{expected_spec}' in workflow content"


def test_build_workflow_valid():
    context = {
        "user": "testuser",
        "repo": "testrepo",
        "name": "Test User",
        "showyourwork_version": "v0.4.3",
        "cache_sandbox_doi": "",
        "overleaf_id": None,
        "year": 2025,
        "action_spec": None,
        "action_version": "v1",
    }
    template_dir = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "../../src/showyourwork/cookiecutter-showyourwork",
        )
    )
    render_and_validate_workflow(template_dir, "build.yml", context, "v1")


def test_build_pull_request_workflow_valid():
    context = {
        "user": "testuser",
        "repo": "testrepo",
        "name": "Test User",
        "showyourwork_version": "v0.4.3",
        "cache_sandbox_doi": "",
        "overleaf_id": None,
        "year": 2025,
        "action_spec": None,
        "action_version": "v1",
    }
    template_dir = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "../../src/showyourwork/cookiecutter-showyourwork",
        )
    )
    render_and_validate_workflow(template_dir, "build-pull-request.yml", context, "v1")


def test_process_pull_request_workflow_valid():
    context = {
        "user": "testuser",
        "repo": "testrepo",
        "name": "Test User",
        "showyourwork_version": "v0.4.3",
        "cache_sandbox_doi": "",
        "overleaf_id": None,
        "year": 2025,
        "action_spec": None,
        "action_version": "v1",
    }
    template_dir = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "../../src/showyourwork/cookiecutter-showyourwork",
        )
    )
    render_and_validate_workflow(
        template_dir,
        "process-pull-request.yml",
        context,
        "v1",
        action_path="showyourwork/showyourwork-action/process-pull-request",
    )


class TestSetupCommandGitScenarios(BaseSetupTest):
    """Test setup command behavior in complex git scenarios using the base class."""

    def test_clean_pypi_version(self):
        """Test that clean release versions use PyPI spec."""
        git_responses = {
            "git status --porcelain": "",  # Clean working tree
        }

        mock_logger, repo_path, workflow_content = self.run_setup_with_git_mocks(
            git_responses,
            "v0.4.3",  # Clean release version
        )

        # Should use PyPI spec for clean release + clean working tree
        self.assert_action_spec_in_workflow(workflow_content, "showyourwork==v0.4.3")

    def test_development_version_fallback(self):
        """Test development version with git command failures falls back to main."""
        # Empty git responses will cause all git commands to return empty strings
        git_responses = {}

        mock_logger, repo_path, workflow_content = self.run_setup_with_git_mocks(
            git_responses, "v0.4.3.dev1+local"
        )

        # Should fallback to default showyourwork git URL
        expected = "git+https://github.com/showyourwork/showyourwork"
        self.assert_action_spec_in_workflow(workflow_content, expected)

    def test_fork_with_upstream_remote(self):
        """Test proper handling of fork with upstream remote."""
        git_responses = {
            "git status --porcelain": "",
            "git branch --show-current": "my-branch",
            "git rev-parse --abbrev-ref --symbolic-full-name @{u}": (
                "upstream/my-branch"
            ),
            "git rev-parse HEAD": "abc123local",
            "git rev-parse upstream/my-branch": "abc123local",  # Same commits
            "git config --get remote.upstream.url": "https://github.com/upstream/showyourwork.git",
        }

        mock_logger, repo_path, workflow_content = self.run_setup_with_git_mocks(
            git_responses, "v0.4.3.dev1+fork"
        )

        # Should use upstream remote URL with the branch
        expected = "git+https://github.com/upstream/showyourwork@my-branch"
        self.assert_action_spec_in_workflow(workflow_content, expected)

    def test_explicit_git_action_spec(self):
        """Test that setup command handles explicit git action_spec correctly."""
        # For explicit action_spec, we bypass the git logic entirely
        with patch("showyourwork.cli.commands.setup.requests"):
            with patch("showyourwork.cli.commands.setup.get_stdout") as mock_get_stdout:
                # Return empty string for all git commands
                mock_get_stdout.return_value = ""

                os.chdir(self.temp_dir)

                # Run setup with explicit git action_spec (bypasses auto-detection)
                setup(
                    slug="testuser/testrepo",
                    cache=False,
                    overleaf_id=None,
                    ssh=False,
                    action_spec="git+https://github.com/user/showyourwork@feature",
                    action_version="v1",
                )

                # Check that the repo was created
                repo_path = os.path.join(self.temp_dir, "testrepo")
                assert os.path.exists(repo_path), "Repository not created!"

                # Check the generated workflow contains our explicit git spec
                workflow_path = os.path.join(
                    repo_path, ".github", "workflows", "build.yml"
                )
                assert os.path.exists(workflow_path), "build.yml not found!"

                with open(workflow_path) as f:
                    workflow_content = f.read()

                # Should contain our explicit git action_spec
                expected_spec = "git+https://github.com/user/showyourwork@feature"
                self.assert_action_spec_in_workflow(workflow_content, expected_spec)

    def test_uncommitted_changes(self):
        """Test setup command behavior with uncommitted changes."""
        git_responses = {
            "git status --porcelain": "M  some_file.py\n?? another.txt",
            "git branch --show-current": "feature-branch",
            "git rev-parse --abbrev-ref --symbolic-full-name @{u}": (
                "upstream/feature-branch"
            ),
            "git rev-parse HEAD": "abc123local",
            "git rev-parse upstream/feature-branch": "def456remote",
            "git config --get remote.upstream.url": (
                "https://github.com/myuser/showyourwork.git"
            ),
        }

        mock_logger, repo_path, workflow_content = self.run_setup_with_git_mocks(
            git_responses, "v0.4.3.dev1+uncommitted"
        )

        # Verify warning was logged about uncommitted changes
        self.assert_warning_logged(mock_logger, "uncommitted changes")

        # Should use git spec pointing to remote branch
        expected = "git+https://github.com/myuser/showyourwork@feature-branch"
        self.assert_action_spec_in_workflow(workflow_content, expected)

    def test_unpushed_commits(self):
        """Test setup command behavior with unpushed commits."""
        git_responses = {
            "git status --porcelain": "",  # Clean working tree
            "git branch --show-current": "my-feature",
            "git rev-parse --abbrev-ref --symbolic-full-name @{u}": (
                "origin/my-feature"
            ),
            "git rev-parse HEAD": "abc123newer",
            "git rev-parse origin/my-feature": "def456older",
            "git config --get remote.origin.url": (
                "git@github.com:forkuser/showyourwork.git"
            ),
        }

        mock_logger, repo_path, workflow_content = self.run_setup_with_git_mocks(
            git_responses, "v0.4.3.dev1+unpushed"
        )

        # Verify warning about unpushed commits
        self.assert_warning_logged(mock_logger, "unpushed commits")

        # Check that git spec points to the remote branch
        expected = "git+https://github.com/forkuser/showyourwork@my-feature"
        self.assert_action_spec_in_workflow(workflow_content, expected)

    def test_no_tracking_branch(self):
        """Test setup when local branch has no remote tracking branch."""
        git_responses = {
            "git status --porcelain": "",
            "git branch --show-current": "local-only-branch",
            "git rev-parse --abbrev-ref --symbolic-full-name @{u}": "",  # No upstream
        }

        mock_logger, repo_path, workflow_content = self.run_setup_with_git_mocks(
            git_responses, "v0.4.3.dev1+local"
        )

        # Verify warning about no tracking branch
        self.assert_warning_logged(mock_logger, "no remote tracking branch")

        # Should fallback to default git URL
        expected = "git+https://github.com/showyourwork/showyourwork"
        self.assert_action_spec_in_workflow(workflow_content, expected)

    def test_both_uncommitted_and_unpushed(self):
        """Test setup with both uncommitted changes AND unpushed commits."""
        git_responses = {
            "git status --porcelain": "M  modified.py\nA  added.py",
            "git branch --show-current": "messy-branch",
            "git rev-parse --abbrev-ref --symbolic-full-name @{u}": (
                "upstream/messy-branch"
            ),
            "git rev-parse HEAD": "local123ahead",
            "git rev-parse upstream/messy-branch": "remote456behind",
            "git config --get remote.upstream.url": (
                "https://github.com/upstream/showyourwork.git"
            ),
        }

        mock_logger, repo_path, workflow_content = self.run_setup_with_git_mocks(
            git_responses, "v0.4.3.dev1+messy"
        )

        # Verify the combined warning message
        self.assert_warning_logged(
            mock_logger, "uncommitted changes AND unpushed commits"
        )

        # Should still use remote git spec
        expected = "git+https://github.com/upstream/showyourwork@messy-branch"
        self.assert_action_spec_in_workflow(workflow_content, expected)
