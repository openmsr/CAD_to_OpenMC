import pytest
import CAD_to_OpenMC.assembly as ab
from tests.harnessRun import HarnessRun
import pathlib as pl
import sys

class HarnessCqSTL2(HarnessRun):
    def __init__(self):
        super().__init__()

    def run(self,merge=False):
        if merge:
            self.merge()

        #self.a.solids_to_h5m(backend='stl')
        h5p = pl.Path('out_cqstl.h5m')
        self.a.solids_to_h5m(backend='stl2',h5m_filename=str(h5p))
        assert h5p.exists()
        assert self.is_validh5m(h5p)
        self.cleanup()

    def cleanup(self):
        super().cleanup()
        pwd=pl.Path('.')
        for v in pwd.glob("vol*_face*"):
            v.unlink()

def testcq():
    t = HarnessCqSTL2()
    t.run()

def testcq_wmerge():
    t = HarnessCqSTL2()
    t.run(merge=True)

if __name__=='__main__':
    testcq()
    testcq_wmerge()
