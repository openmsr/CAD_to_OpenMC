import pytest
import assembly as asmb
import pathlib as pl

def test_test():
    assert 1==1

@pytest.mark.parametrize("infile",['1pin_hollow_contact.step'])
def test_assembly(infile="in.step"):
    file=pl.Path(infile)
    a=asmb.Assembly()
    a.import_stp_file(file)
    string=a.export_stl(file.with_suffix('.stl'))
    print(string)
    assert string == file.with_suffix('.stl')
