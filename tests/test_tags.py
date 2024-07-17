import pytest
import CAD_to_OpenMC.assembly as ab
from tests.harnessRun import HarnessRun
import pathlib as pl
import sys

class HarnessDB(HarnessRun):
    def __init__(self):
        super().__init__()

    def run(self,merge=False, cleanup=True):
        if merge:
            self.merge()

        h5p = pl.Path('out_db.h5m')

        self.a.solids_to_h5m(backend='db',h5m_filename=str(h5p))
        assert h5p.exists()
        assert self.is_validh5m(h5p)

        if cleanup:
            self.cleanup()


    def check_tags(self,extra_tags=[]):
        if self.a.tags is not None:
            for tag in self.a.tags.values():
                p1=sp.run(['grep','-qa',tag,str(h5p)])
                assert p1.returncode == 0

        for tag in check_tags:
            p1=sp.run(['grep','-qa',tag,str(h5p)])
            assert p1.returncode == 0

    def cleanup(self):
        super().cleanup()
        pwd=pl.Path('.')
        for v in pwd.glob("vol*_face*"):
            v.unlink()

def testdb_wtags():
    t = HarnessDB()
    t.a.tags={'h2.*':'water','zirconium':'Zi','u[0-9]o','uranium_oxide'}
    t.run(merge=True, cleanup=False)
    t.check_tags()
    t.cleanup()

def testdb_wpartialtags():
    t = HarnessDB()
    t.a.tags={'h2.*':'water','u[0-9]o','uranium_oxide'}
    t.run(merge=True, cleanup=False)
    t.check_tags(['zirconium'])
    t.cleanup()

if __name__=='__main__':
    testdb()
    testdb_wmerge()
