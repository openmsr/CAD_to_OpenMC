import pytest
import src.h5massembly.assembly as asmb
import pathlib as pl

class TestRun()

  def __init__(self):
    pass

  def test_brep_export(self):
    infile = 'examples/7pins.step'
    self.a = asmb.Assembly(verbose=2)
    self.a.stp_files = [infile]
    self.a.import_stp_file()
    self.ofile = 'out.brep'
    self.h5m = 'out.h5m'
    assert a.export_brep(ofile)

  def test_gmsh(self):
    brep_to_h5m(ofile, h5m_filename=h5m, samples=20, min_size_mesh=0.5, max_size_mesh=10)
    h5p = pl.Path(h5m)
    assert h5p.exists()


  def test_stl(self):
    pass
