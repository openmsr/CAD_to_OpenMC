import cadquery2 as cq
import subprocess as sp
import pathlib as pl
import hashlib as hl
from .stl_utils import *
from . import meshutils

class MesherCQSTL:
  def __init__(self, tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default, refine, entities):
    self.tolerance=tolerance
    self.angular_tolerance=angular_tolerance
    self.entities=entities
    self.verbose=2
    self.refine=refine
    self.default=default
    self.min_mesh_size=min_mesh_size
    self.max_mesh_size=max_mesh_size

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
    stlp=pl.Path(stl)
    with open(stl) as fp:
      buf,n=read_stl_ascii(fp)
    verts=buffer2vertices(buf)
    tris=buffer2triangles(buf,verts)
    edges=np.array(meshutils.find_edges(tris))
    meshutils.write_dotmesh(str(stlp.with_suffix('.mesh')),verts,tris,edges=edges,required_edges='all')
    cp=sp.run(['mmgs_O3','-hmin',f'{self.min_mesh_size}','-hmax',f'{self.max_mesh_size}','-optim','-in',stlp.with_suffix('.mesh'),'-out',stlp.with_suffix('.o.mesh')], capture_output=True)
    print(cp.stdout)

    import gmsh
    gmsh.initialize()
    gmsh.open(str(stlp.with_suffix('.o.mesh')))
    gmsh.write( str(stlp) )
    gmsh.finalize()

  def _refine_volumefaces(self,stls,output_stl='out.stl'):
    #refine a volume consistng of a set of stl-files which is expected to contains its parts
    all_verts=None
    for i,vf in enumerate(stls):
      print(f'processing: {vf}')
      #if this is reuse of an already refined surface we should extract the information from there instead.
      #check if the file exists
      vfp=pl.Path(vf)
      if False and vfp.with_suffix('.o.mesh').exists():
        print(f'found refined file {str(vfp.with_suffix(".o.mesh"))} extracting information from that')
        import gmsh
        gmsh.initialize()
        gmsh.open(str(vfp.with_suffix('.o.mesh')))
        gmsh.write(str(vfp.with_suffix('.tmp.stl')))
        gmsh.finalize()

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
        #these new triangle labels need an offset of the current vertex count
        all_triangles=np.vstack((all_triangles,triangles+vertex_count))
        all_tlabels=np.append(all_tlabels,(i+1)*np.ones((triangles.shape[0]),dtype='uint'))
        all_vlabels=np.append(all_vlabels,1*np.ones((verts.shape[0]),dtype='uint'))
        vertex_count+=verts.shape[0]
    stlp=pl.Path(output_stl)
    meshutils.write_dotmesh(stlp.with_suffix('.mesh'),all_verts,all_triangles,vertex_labels=all_vlabels, triangle_labels=all_tlabels)
    meshutils.write_dummy_dotsol(stlp.with_suffix('.sol'),all_verts.shape[0])

    cp=sp.run(['mmgs_O3','-ls','-sol',stlp.with_suffix('.sol'),'-in',stlp.with_suffix('.mesh'),'-out',stlp.with_suffix('.o.mesh'),'-keep-ref'], capture_output=True)
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
            self._refine_stls(facename)
      #merge the stls to a single .stl
      merge_stl(volname, volumefaces,of='bin')
      e.stl=volname

class MesherCQSTLBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default, refine, entities, **_ignored):
    if not self._instance:
      self._instance = MesherCQSTL(tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default,refine, entities)
    return self._instance
