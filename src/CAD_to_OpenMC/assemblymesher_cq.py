import cadquery2 as cq
import subprocess as sp
import pathlib as pl
import hashlib as hl

class MesherCQSTL:
  def __init__(self, tolerance, angular_tolerance, default, refine, entities):
    self.tolerance=tolerance
    self.angular_tolerance=angular_tolerance
    self.entities=entities
    self.verbose=2
    self.refine=refine

  def generate_stls(self):
    #created a cq-compund from list of entities
    for i,e in enumerate(self.entities):
      j=i+1
      filename=f"volume_{j}.stl"
      cq.exporters.export(e.solid,filename,exportType="STL",tolerance=self.tolerance,angularTolerance=self.angular_tolerance)
      if(self.verbose>1):
        print(f"INFO: cq export to file {filename}")
      e.stl=filename
      if(self.refine):
        self._refine_stls(e.stl)

  def _refine_stls(self,stl):
    try:
      if(self.refine.startswith('mmg')):
        self._mmgs_refine_stls(stl)
      else:
        print(f'ERROR: Unknown mesh refinement algorithm: {self.refine}')
    except:
        self._mmgs_refine_stls(stl)

  def _mmgs_refine_stls(self,stl):
    #call mmgs to refine the mesh
    try:
      cp=sp.run(['mmgs_O3','--help'],capture_output=True)
    except:
      print(f'WARNING: Cannot find mmgs mesh refinement tool. Vol. {stl} will not be remeshed. Did you forget to include it on the PATH?')
      return
    import gmsh
    gmsh.initialize()
    #for now we use gmsh as a converter from stl to inria mesh format - this should be changed.
    stlp=pl.Path(stl)
    gmsh.open(str(stlp))
    gmsh.write( str(stlp.with_suffix('.mesh')) )
    cp=sp.run(['mmgs_O3','-hausd','0.1','-optim','-in',stlp.with_suffix('.mesh'),'-out',stlp.with_suffix('.o.mesh')], capture_output=True)
    print(cp.stdout)

    gmsh.open(str(stlp.with_suffix('.o.mesh')))
    gmsh.write( str(stlp) )
    gmsh.finalize()


  def _mesh_surfaces(self):
    #loop over all surfaces in all entities
    #and generate meshes (refined or otherwise)
    facehashtable={}
    for i,e in enumerate(self.entities):
      j=i+1
      volname= f"volume_{j}.stl"
      k=0
      for f in e.solid.Faces():

        m=hl.sha256()
        hh=m.update(f)
        if hh not in keys(facehashtable):
          facehastable.update(hh,facename)
        else:

        #export to stl
        #refine?
          #if so reimport






class MesherCQSTLBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, tolerance, angular_tolerance, default, refine, entities, **_ignored):
    if not self._instance:
      self._instance = MesherCQSTL(tolerance, angular_tolerance, default,refine, entities)
    return self._instance
