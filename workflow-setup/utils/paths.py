from pathlib import Path


workflow = Path(__file__).absolute().parents[1]

rules = workflow / "rules"
checkpoints = workflow / "checkpoints"
