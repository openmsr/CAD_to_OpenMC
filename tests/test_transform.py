import pytest
import CAD_to_OpenMC.assembly as ab
from tests.harnessRun import HarnessRun
import pathlib as pl

class HarnessTrans(HarnessRun):
  def __init__(self):
    super().__init__()

  def run(self,merge=False, **kwargs):
    print(f'stl: {self.a.stp_files})')
    print(self.h5m)

    if merge:
      self.merge()

    #self.a.solids_to_h5m(backend='stl')
    h5p = pl.Path('out_transformed.h5m')
    self.a.solids_to_h5m(backend='stl',h5m_filename=str(h5p), kwargs)
    assert h5p.exists()
    assert self.is_validh5m(h5p)
    self.cleanup()

  def cleanup(self):
    super().cleanup()
    pwd=pl.Path('.')
    for v in pwd.glob("vol*_face*"):
      v.unlink()

def test_scale():
  t = HarnessTrans()
  t.run(scale = 2.0)

def test_translate():
  t = HarnessTrans()
  t.run(translate = [[1,2],[0.0, 0.0, -30.0],[3],[0.0, 0.0, 30]])
