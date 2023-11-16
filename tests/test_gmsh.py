import pytest
import CAD_to_OpenMC.assembly as ab
from tests.harnessRun import HarnessRun
import pathlib as pl

class HarnessGmsh(HarnessRun):
  def __init__(self):
    super().__init__()

  def run(self):
    print(f'gmsh: {self.a.stp_files})')
    self.a.merge_all()
    ab.mesher_config['refine']=False
    h5p = pl.Path('out_gmsh.h5m')
    self.a.solids_to_h5m(backend='gmsh',h5m_filename=str(h5p))
    assert h5p.exists()
    assert self.is_validh5m(h5p)
    #self.cleanup()

def testgmsh():
  t = HarnessGmsh()
  t.run()
