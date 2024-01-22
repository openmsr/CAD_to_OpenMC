import os
import contextlib
from pathlib import Path
from shutil import rmtree

@contextlib.contextmanager
def mesher_datadir(path,delete_intermediate=True, movein=False):
    """Changes working directory and returns to previous on exit."""
    prev_cwd = Path.cwd()
    Path(path).mkdir(parents=True, exist_ok=False)
    if movein:
        os.chdir(path)
    try:
        yield
    finally:
        g=list(Path().glob('*.h5m')) + list(Path().glob('*.vtk'))
        for f in g:
            os.replace(f,prev_cwd / f)
        os.chdir(prev_cwd)
        if (delete_intermediate):
            rmtree(path)
