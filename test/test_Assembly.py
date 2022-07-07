import pytest
import h5massembly.assembly as ab
import pathlib as pl

class TestRun:
  def __init__(self):
    infile = 'examples/7pin.step'
    self.a = ab.Assembly(verbose=2)
    self.a.stp_files = [infile]
    self.a.import_stp_files()
    self.h5m = 'out.h5m'

  def test_gmsh(self):
    self.a.solids_to_h5m(backend='gmsh')
    h5p = pl.Path('dagmc.h5m')
    assert h5p.exists()
    h5p.unlink()
    self.cleanup

  def test_stl(self):
    self.a.solids_to_h5m(backend='stl')
    h5p = pl.Path('dagmc.h5m')
    assert h5p.exists()
    self._cleanup()

  def _cleanup(self):
    pwd=pl.Path('.')
    for v in pwd.glob('volume_*.stl'):
      v.unlink()
    (pwd / 'dagmc.h5,').unlink()
    (pwd / 'dagmc.vtk,').unlink()


t=TestRun()
t.test_stl()
t.test_gmsh()
