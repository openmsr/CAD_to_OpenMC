import openmc
import openmc.lib
import glob
import os

class DAGMC_template():
  def __init__(self,dagmc=None):
    self.dagmc=dagmc
    self.materials=None
    self.bld_sts()
    self.bld_mats()
    self.bld_geom()

  def bld_geom(self):
    if self.dagmc:
      du=openmc.DAGMCUniverse(self.dagmc,auto_geom_ids=True)
    else:
      if not self.materials:
        self.bld_mats()
      s1=openmc.Sphere(r=1.0)
      c1=openmc.Cell(region=-s1,fill=self.materials[0])
      c2=openmc.Cell(region=+s1,fill=self.materials[-1])
      du=openmc.Universe()
      du.add_cells([c1,c2])
    bdr=openmc.Sphere(r=200.0,boundary_type='vacuum')
    scn=openmc.Cell(region=-bdr, fill=du)
    ru=openmc.Universe()
    ru.add_cell(scn)
    self.geometry=openmc.Geometry(ru)

  def bld_sts(self):
    sets=openmc.Settings()
    sets.batches=5
    sets.particles=1000
    sets.incactive=2
    src=openmc.Source()
    src.space=openmc.stats.Point()
    src.energy=openmc.stats.Discrete([2.0e6],[1.0])
    src.angle=openmc.stats.Isotropic()
    sets.source=src
    self.settings=sets

  def bld_mats(self):
    if self.materials is not None:
      return
    uo2=openmc.Material(name='uo2')
    uo2.add_nuclide('U235',1.0/3.0*0.05,'ao')
    uo2.add_nuclide('U238',1.0/3.0*0.95,'ao')
    uo2.add_element('O',2.0/3.0,'ao')
    uo2.set_density('g/cc',10.97)

    zirconium=openmc.Material(name='zirconium')
    zirconium.add_element('Zr',1.0,'ao')
    zirconium.set_density('g/cc',6.52)

    h2o=openmc.Material(name='h2o')
    h2o.add_element('O',1.0/3.0,'ao')
    h2o.add_element('H',2.0/3.0,'ao')
    h2o.set_density('g/cc',1.0)
    self.materials=openmc.Materials([uo2,zirconium,h2o])

  def export_to_xml(self):
    self.materials.export_to_xml()
    self.settings.export_to_xml()
    self.geometry.export_to_xml()

  def run(self):
    self.export_to_xml()
    openmc.lib.init()
    openmc.lib.run()
    openmc.lib.finalize()

  def cleanup(self):
    g1=glob.glob('*.xml')
    for g in g1:
      os.unlink(g)
    g2=glob.glob('*.h5')
    for g in g2:
      os.unlink(g)
    g3=glob.glob('*.vtk')
    for g in g3:
      os.unlink(g)
