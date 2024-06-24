import pytest
import CAD_to_OpenMC.assembly as ab
from tests.harnessRun import HarnessRun
import pathlib as pl
import sys

class HarnessDB(HarnessRun):
    def __init__(self):
        super().__init__()

    def run(self,merge=False):
        if merge:
            self.merge()

        h5p = pl.Path('out_db.h5m')
        self.a.solids_to_h5m(backend='db',h5m_filename=str(h5p))
        assert h5p.exists()
        assert self.is_validh5m(h5p)
        self.cleanup()

    def cleanup(self):
        super().cleanup()
        pwd=pl.Path('.')
        for v in pwd.glob("vol*_face*"):
            v.unlink()

def testcq():
    t = HarnessDB()
    t.run()

def testcq_wmerge():
    t = HarnessDB()
    t.run(merge=True)

if __name__=='__main__':
    testdb()
    testdb_wmerge()
