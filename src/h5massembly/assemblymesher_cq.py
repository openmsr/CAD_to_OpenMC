import cadquery2 as cq
from assemblymesher import *

class MesherCQSTL:
  def __init__(self, tolerance, angular_tolerance, default, solids, verbose):
    self.tolerance=tolerance
    self.angular_tolerance=angular_tolerance
    self.solids=solids
    self.verbose=verbose

  def generate_stls(self):
    stls=[]
    #If not merged we should operate directly on the entities objects
    #For now this is a hack relying on the fact that merged is a compund object.

    for i,s in enumerate(self.solids):
      j=i+1
      filename=f"volume_{j}.stl"
      status=cq.exporters.export(s,filename,exportType="STL",tolerance=self.tolerance,angularTolerance=self.angular_tolerance)
      if(self.verbose>1):
        print(f"INFO: cq export to file {filename}:{status}")
      stls.append((j,filename))
    return stls


class MesherCQSTLBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, tolerance, angular_tolerance, default, solids, verbose, **_ignored):
    if not self._instance:
      self._instance = MesherCQSTL(tolerance, angular_tolerance, default, solids, verbose)
    return self._instance
