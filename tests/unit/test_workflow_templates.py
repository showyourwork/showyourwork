import os
import shutil
import tempfile

import yaml
from cookiecutter.main import cookiecutter


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


def test_build_workflow_valid():
    context = {
        "user": "testuser",
        "repo": "testrepo",
        "name": "Test User",
        "showyourwork_version": "1.0.0",
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
        "showyourwork_version": "1.0.0",
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
        "showyourwork_version": "1.0.0",
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
