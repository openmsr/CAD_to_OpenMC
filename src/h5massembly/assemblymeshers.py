class MesherFactory():
  def __init__(self):
    self._creators={}

  def register_creator(self,tag,creator):
    self._creators[tag]=creator

  def get_mesher(self,tag):
    creator=self._creators.get(tag)
    if not tag:
      raise ValueError(tag)
    return creator

factory=MesherFactory()
factory.register_backend('gmsh',MesherBuilderGMSH)
factory.register_backend('stl',MesherBuilderCQSTL)

