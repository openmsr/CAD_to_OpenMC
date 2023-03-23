import pytest
import CAD_to_OpenMC.assembly as ab
import pathlib as pl

class TestRun:
  def __init__(self):
    infile = 'examples/7pin.step'
    self.a = ab.Assembly(verbose=2)
    self.a.stp_files = [infile]
    self.a.import_stp_files()
    self.h5m = 'out.h5m'

  def merge(self):
    self.a.merge_all()

  def cleanup(self, h5m=None):
    pwd=pl.Path('.')
    for v in pwd.glob('volume_*.stl'):
      v.unlink()
    h5p=pl.Path(h5m) if h5m else pl.Path(self.h5m)
    try:
      (pwd / h5p ).unlink()
      (pwd / h5p.with_suffix('vtk')).unlink()
    except:
      pass

  def is_validh5m(self,file):
    with open(file,"rb") as f:
      magic_bytes=f.read(8)
      if(magic_bytes!=b'\x89HDF\x0d\x0a\x1a\x0a'):
         print(f'ERROR: generated file {file} does not appear to be a hdf-file. Did you compile the moab libs with HDF enabled?')
         return False
    return True
