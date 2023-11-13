import pathlib
import os
import openmc

import CAD_to_OpenMC.assembly as ab

from tests.OMC_DAGMC_run import DAGMC_template
from tests.harnessRun import HarnessRun

class OMC_DAGMC_harness(HarnessRun):
  def __init__(self, step):
    self.path=pathlib.Path(step)
    self.h5m=self.path.with_suffix('.h5m')
    self.nuclear_lib=pathlib.Path('tests/nuclear_data_testlib/cross_sections.xml').absolute()

  def run(self):
    aa=ab.Assembly([str(self.path)],verbose=2)
    aa.run(backend='stl2',merge=True,h5m_filename=self.h5m)
    assert self.h5m.exists()
    tt=DAGMC_template(self.h5m)
    tt.run()
    assert pathlib.Path('statepoint.5.h5').exists()
    tt.cleanup()
    self.cleanup()

def test_h5m_neutronics_p1():
  o=OMC_DAGMC_harness('examples/pincell1.step')
  openmc.config['cross_sections']=str(o.nuclear_lib)
  o.run()

def test_h5m_neutronics_p2():
  o=OMC_DAGMC_harness('examples/pincell2.step')
  openmc.config['cross_sections']=str(o.nuclear_lib)
  o.run()

def test_h5m_neutronics_tors():
  o=OMC_DAGMC_harness('examples/toroids.step')
  openmc.config['cross_sections']=str(o.nuclear_lib)
  o.run()

if __name__=='__main__':
  test_h5m_neutronics_p1()
  test_h5m_neutronics_p2()
  test_h5m_neutronics_tors()
