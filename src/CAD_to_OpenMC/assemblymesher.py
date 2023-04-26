import abc
from .assemblymesher_base import *
from .object_factory import ObjectFactory
from .assemblymesher_gmsh import *
from .assemblymesher_cq import *
from .assemblymesher_cq2 import *

class MesherFactory(ObjectFactory):
  def get(self,mesher_id,**kwargs):
    return self.create(mesher_id, **kwargs)

meshers=MesherFactory()
meshers.register_builder('gmsh',MesherGMSHBuilder())
meshers.register_builder('stl',MesherCQSTLBuilder())
meshers.register_builder('stl2',MesherCQSTL2Builder())
