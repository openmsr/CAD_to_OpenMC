import cadquery2 as cq

class MesherCQSTL:
  def __init__(self, tolerance, angular_tolerance, default, entities):
    self.tolerance=tolerance
    self.angular_tolerance=angular_tolerance
    self.entities=entities
    self.verbose=2

  def generate_stls(self):
    stls=[]
    #created a cq-compund from list of entities
    for i,e in enumerate(self.entities):
      j=i+1
      filename=f"volume_{j}.stl"
      status=cq.exporters.export(e.solid,filename,exportType="STL",tolerance=self.tolerance,angularTolerance=self.angular_tolerance)
      if(self.verbose>1):
        print(f"INFO: cq export to file {filename}:{status}")
      e.stl=filename
    return stls

class MesherCQSTLBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, tolerance, angular_tolerance, default, entities, **_ignored):
    if not self._instance:
      self._instance = MesherCQSTL(tolerance, angular_tolerance, default, entities)
    return self._instance
