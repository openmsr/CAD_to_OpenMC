import object_factory
from assemblymesher_gmsh import *
from assemblymesher_cq import *

class MesherFactory(object_factory.ObjectFactory):
  def get(self,mesher_id,**kwargs):
    return self.create(mesher_id, **kwargs)

meshers=MesherFactory()
meshers.register_builder('gmsh',MesherGMSHBuilder())
meshers.register_builder('stl',MesherCQSTLBuilder())
