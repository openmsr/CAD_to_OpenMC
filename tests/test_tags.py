import pytest
import CAD_to_OpenMC.assembly as ab
from tests.harnessRun import HarnessRun
import pathlib as pl
import subprocess as sp
import sys

class HarnessDB(HarnessRun):
    def __init__(self):
        super().__init__(infile='examples/step_files/pincell1.step')
        self.h5p = pl.Path('out_db.h5m')

    def run(self,merge=False, cleanup=True, tags=None, sequential_tags=None):
        self.a.import_stp_files(tags=tags,sequential_tags=sequential_tags)
        if merge:
            self.merge()

        self.a.solids_to_h5m(backend='db',h5m_filename=str(self.h5p))
        assert self.h5p.exists()
        assert self.is_validh5m(self.h5p)

        if cleanup:
            self.cleanup()


    def check_tags(self,tags_to_check=[]):
        for tag in tags_to_check:
            p1=sp.run(['grep','-qa',tag,str(self.h5p)])
            assert p1.returncode == 0

    def cleanup(self):
        super().cleanup()
        pwd=pl.Path('.')
        for v in pwd.glob("vol*_face*"):
            v.unlink()

def testdb_seqtags():
    stags=['mat0','mat1','mat2']
    t = HarnessDB()
    t.run(merge=True, cleanup=False, sequential_tags=stags)
    t.check_tags(['mat0','mat2'])
    t.cleanup()

def testdb_wtags():
    tags={'h2.*':'water','zirconium':'Zi','uo[0-9]':'uranium_oxide'}
    t = HarnessDB()
    t.run(merge=True, cleanup=False, tags=tags)
    t.check_tags(tags.values())
    t.cleanup()

def testdb_wpartialtags():
    tags={'h2.*':'water','uo[0-9]':'uranium_oxide'}
    t = HarnessDB()
    t.run(merge=True, cleanup=False, tags=tags)
    t.check_tags(['zirconium'])
    t.cleanup()

def testdb_switch_delims():
    t = HarnessDB()
    t.a.set_tag_delim('o')
    t.run(merge=True, cleanup=False)
    t.check_tags(['h2','zirc','u'])
    t.cleanup()

if __name__=='__main__':
    testdb()
    testdb_wmerge()
