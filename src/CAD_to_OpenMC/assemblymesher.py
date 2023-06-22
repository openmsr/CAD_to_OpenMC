import abc
from .assemblymesher_base import *
from .object_factory import ObjectFactory

nogmsh=False
try:
  from .assemblymesher_gmsh import *
except ImportError as e:
  nogmsh=e

nocq=False
try:
  from .assemblymesher_cq import *
except ImportError as e:
  nocq=e

nocq2=False
try:
  from .assemblymesher_cq2 import *
except ImportError as e:
  npcq2=e

if (nogmsh and nocq and nocq2):
  raise ImportError("Could not import any of the mesher backends")

class MesherFactory(ObjectFactory):
  def get(self,mesher_id,**kwargs):
    return self.create(mesher_id, **kwargs)

meshers=MesherFactory()
meshers.register_builder('gmsh',MesherGMSHBuilder())
meshers.register_builder('stl',MesherCQSTLBuilder())
meshers.register_builder('stl2',MesherCQSTL2Builder())
