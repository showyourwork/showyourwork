from tests.integration.helpers.temp_repo import TemporaryShowyourworkRepository


class ExampleRepository(TemporaryShowyourworkRepository):
    __test__ = False


def test_repo_name_without_run_id(monkeypatch):
    monkeypatch.delenv("GITHUB_RUN_ID", raising=False)

    repo = ExampleRepository()

    assert repo.repo == "example-repository"


def test_repo_name_with_run_id(monkeypatch):
    monkeypatch.setenv("GITHUB_RUN_ID", "123456")

    repo = ExampleRepository()

    assert repo.repo == "example-repository-123456"


def test_repo_name_changes_between_runs(monkeypatch):
    repo = ExampleRepository()

    monkeypatch.setenv("GITHUB_RUN_ID", "111111")
    first_name = repo.repo

    monkeypatch.setenv("GITHUB_RUN_ID", "222222")
    second_name = repo.repo

    assert first_name != second_name
