import importlib

for plugin in config.get("plugins", []):
    mod = importlib.import_module(plugin)
    if not hasattr(mod, "snakefiles"):
        # TODO(dfm): log warning
        continue
    for snakefile in mod.snakefiles():
        include: snakefile
