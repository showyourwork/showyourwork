from typing import Any, Dict, Optional

_CONFIG: Dict[str, Any] = {}


def configure(_config: Optional[Dict[str, Any]] = None, **config: Any) -> None:
    if _config is not None:
        _CONFIG.update(_config)
    _CONFIG.update(config)
