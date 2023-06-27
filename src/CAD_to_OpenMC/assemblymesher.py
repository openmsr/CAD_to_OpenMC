import abc
from .assemblymesher_base import assemblymesher
from .object_factory import ObjectFactory

nogmsh=False
try:
  from .assemblymesher_gmsh import MesherGMSHBuilder
except (ImportError,OSError) as e:
  nogmsh=e

nocq=False
try:
  from .assemblymesher_cq import MesherCQSTLBuilder
except (ImportError,OSError) as e:
  nocq=e

nocq2=False
try:
  from .assemblymesher_cq2 import MesherCQSTL2Builder
except (ImportError,OSError) as e:
  npcq2=e

if (nogmsh and nocq and nocq2):
  raise ImportError("Could not import any of the mesher backends")

class MesherFactory(ObjectFactory):
  def get(self,mesher_id,**kwargs):
    return self.create(mesher_id, **kwargs)

meshers=MesherFactory()
if(not nogmsh):
  meshers.register_builder('gmsh',MesherGMSHBuilder())
if(not nocq):
  meshers.register_builder('stl',MesherCQSTLBuilder())
if(not nocq2):
  meshers.register_builder('stl2',MesherCQSTL2Builder())
