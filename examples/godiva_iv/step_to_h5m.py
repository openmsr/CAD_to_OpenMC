import CAD_to_OpenMC.assembly as ab
import pathlib as pl
for fn in ['GIV_core_detailed.step','GIV_BR_detailed.step','GIV_CR_detailed.step']:
    h5m_fn=pl.Path(fn).with_suffix('.h5m')
    A = ab.Assembly([fn])
    A.verbose=2
    A.set_tag_delim('@\ ')
    A.run(backend='db',h5m_filename=h5m_fn)
