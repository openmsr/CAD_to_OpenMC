import pytest
import assembly as asmb
import pathlib as pl

def test_test():
    assert 1==1

@pytest.mark.parametrize("infile",['examples/7pins.step'])
def test_brep_export(infile: str="in.step"):
    file=pl.Path(infile)
    a=asmb.Assembly(verbose=2)
    a.stp_files=[infile]
    a.import_stp_file()
    ofile='out.brep'
    h5m='out.h5m'
    assert a.export_brep(ofile)
    a.brep_to_h5m(ofile, h5m_filename=h5m, samples=20, min_size_mesh=0.5, max_size_mesh=10)
    h5p=pl.Path(h5m)
    assert h5p.exists()
