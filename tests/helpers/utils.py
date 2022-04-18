import yaml
from contextlib import contextmanager


@contextmanager
def edit_yaml(file):
    with open(file, "r") as f:
        contents = yaml.load(f, Loader=yaml.CLoader)
    try:
        yield contents
    finally:
        with open(file, "w") as f:
            print(yaml.dump(contents, Dumper=yaml.CDumper), file=f)