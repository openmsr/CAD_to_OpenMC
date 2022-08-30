import cadquery2 as cq
import subprocess as sp
import pathlib as pl
import hashlib as hl
from .stl_utils import *
from . import meshutils

class MesherCQSTL:
  def __init__(self, tolerance, angular_tolerance, default, refine, entities):
    self.tolerance=tolerance
    self.angular_tolerance=angular_tolerance
    self.entities=entities
    self.verbose=2
    self.refine=refine

  def generate_stls(self):
    self._mesh_surfaces()
    return
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

  def _refine_volumefaces(self,stls,output_stl='out.stl'):
    #refine a volume consistng of a set of stl-files which is expected to contains its parts
    totalverts=None
    totaltriangles=None
    for i,vf in enumerate(volumefaces):
      buf,n=read_stl(vf)
      verts=buffer2vertices(buf)
      triangles=buffer2triangles(buf,verts)
      if (all_verts is None):
        all_verts=verts
        all_triangles=triangles
        all_tlabels=(i+1)*np.ones((triangles.shape[0]),dtype='uint')
        all_vlabels=1*np.ones((verts.shape[0]),dtype='uint')
        vertex_count=verts.shape[0]
      else:
        all_verts=np.vstack((all_verts,verts))
        all_triangles=np.vstack((all_triangles,triangles)
        all_tlabels=np.vstack((all_tlabels,(i+1)*np.ones((triangles.shape[0],1),dtype='uint')))
        all_vlabels=np.vstack((all_vlabels,(i+1)*np.ones((verts.shape[0],1),dtype='uint')))

    meshutils.write_dotmesh(stlp.with_suffix('.mesh'))
    meshutils.write_dummy_dotsol(stlp.with_suffix('.sol'),allverts.shape[0])

    cp=sp.run(['mmgs_O3','-hausd','0.1','-optim',,'-sol',stlp.with_suffix('.sol'),'-in',stlp.with_suffix('.mesh'),'-out',stlp.with_suffix('.o.mesh'),'-keep-ref'], capture_output=True)
    print(cp.stdout)
    import gmsh
    gmsh.initialize()
    gmsh.open(str(stlp.with_suffix('.o.mesh')))
    gmsh.write( str(stlp) )
    gmsh.finalize()

  def _mesh_surfaces(self):
    #loop over all surfaces in all entities
    #and generate meshes (refined or otherwise)
    facehashtable={}
    for i,e in enumerate(self.entities):
      volname= f"volume_{i+1}.stl"
      k=0
      volumefaces=[]
      for j,f in enumerate(e.solid.Faces()):
        hh=hash(f)
        if hh in facehashtable.keys():
          #surface is in table - use that file for this volume
          volumefaces.append(facehashtable[hh])
          print(f'reusing {hh} {facehashtable[hh]}')
        else:
          facename=f'vol_{i+1}_face{j}.stl'
          facehashtable[hh]=facename
          print(len(facehashtable.keys()))
          cq.exporters.export(f,facename, exportType="STL", tolerance=self.tolerance, angularTolerance=self.angular_tolerance)
          if(True or self.verbose>1):
            print(f"INFO: cq export to file {facename}")
          volumefaces.append(facename)
      if (self.refine):
        self._refine_volumefaces(volumenfaces,stl_file=volname)
      else:
        #merge the stls to a single .stl
        merge_stl(volname, volumefaces,of='bin')
      e.stl=volname
      #we have now a full stl-description of volume with completely imprinted surfaces
          #surface is not in table - we need to mesh (and possibly remesh) it
          #and put it into the main facetable as well as the local table for this volume.
          #refine?
          #if so reimport

class MesherCQSTLBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, tolerance, angular_tolerance, default, refine, entities, **_ignored):
    if not self._instance:
      self._instance = MesherCQSTL(tolerance, angular_tolerance, default,refine, entities)
    return self._instance
