import importlib

for plugin in config.get("plugins", []):
    mod = importlib.import_module(plugin)
    if not hasattr(mod, "snakefile"):
        # TODO(dfm): log warning
        continue
    for snakefile in mod.snakefiles():
        include: snakefile
