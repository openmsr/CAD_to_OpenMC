import cadquery as cq
import subprocess as sp
import pathlib as pl
import hashlib as hl
from CAD_to_OpenMC import cdtemp

single_thread_override=True
try:
  import multiprocessing.pool as mp
except:
  single_thread_override=True

from .stl_utils import *
from . import meshutils



#following pattern here: https://thelaziestprogrammer.com/python/multiprocessing-pool-a-global-solution
#to enable multiprocessing meshings

class MesherCQSTL:
  #these need to be class attributes to avoid pickling when spawning a multiprocessing-pool
  cq_mesher_entities=None
  cq_mesher_tolerance=None
  cq_mesher_ang_tolereance=None
  cq_mesher_min_mesh_size=None
  cq_mesher_max_mesh_size=None

  cq_mesher_faceHash={}

  def __init__(self, tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default, refine, threads, entities):
    self._set_meshpars(tolerance,angular_tolerance,min_mesh_size,max_mesh_size)
    self.verbose=2
    self.refine=refine
    self._set_entities(entities)
    self.default=default
    self.min_mesh_size=min_mesh_size
    self.max_mesh_size=max_mesh_size
    self.threads=threads

  @classmethod
  def _set_meshpars(cls,tol,ang_tol, min_sz, max_sz):
    cls.cq_mesher_tolerance=tol
    cls.cq_mesher_ang_tolerance=ang_tol
    cls.cq_mesher_min_mesh_size=min_sz
    cls.cq_mesher_max_mesh_size=max_sz

  @classmethod
  def _set_entities(cls,entities):
    cls.cq_mesher_entities=entities

  def generate_stls(self):
    self._mesh_surfaces()
    return
    #created a cq-compund from list of entities
    for i,e in enumerate(self.cq_mesher_entities):
      j=i+1
      filename=f"volume_{j}.stl"
      e.solid.exportStl(filename,ascii=True,tolerance=self.tolerance,angularTolerance=self.angular_tolerance)
      if(self.verbose>1):
        print(f"INFO: cq export to file {filename}")
      e.stl=filename
      if(self.refine):
        self._refine_stls(e.stl)

  @classmethod
  def _refine_stls(cls,stl,refine):
    try:
      if(refine.startswith('mmg')):
        cls._mmgs_refine_stl(stl)
      else:
        print(f'ERROR: Unknown mesh refinement algorithm: {self.refine}')
    except:
        cls._mmgs_refine_stl(stl)

  @classmethod
  def _mmgs_refine_stl(cls,stl):
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
    cp=sp.run(['mmgs_O3','-hmin',f'{cls.cq_mesher_min_mesh_size}','-hmax',f'{cls.cq_mesher_max_mesh_size}','-optim','-in',stlp.with_suffix('.mesh'),'-out',stlp.with_suffix('.o.mesh')], capture_output=True)
    print(cp.stdout)

    import gmsh
    gmsh.initialize()
    gmsh.open(str(stlp.with_suffix('.o.mesh')))
    gmsh.write( str(stlp) )
    gmsh.finalize()

  def _mesh_surfaces(self):
    #loop over all surfaces in all entities
    #and generate meshes (refined or otherwise)
    cwd=pl.Path.cwd()
    #create a workplace in tmp
    with cdtemp() as mngr:
      for i,e in enumerate(self.cq_mesher_entities):
        volname= f"volume_{i+1}.stl"
        k=0
        mpargs=[(j,i,self.refine) for j,f in enumerate(e.solid.Faces())]
        if (single_thread_override):
          output=[]
          for args in mpargs:
            output.append(self._mesh_single(*args))
        else:
          pool=mp.Pool(processes=self.threads)
          output=pool.starmap(self._mesh_single, mpargs)
        volumefaces=[]
        for o in output:
          self.cq_mesher_faceHash[o[0]]=o[1]
          volumefaces.append(o[1])
        #merge the stls to a single .stl in the working directory
        merge_stl(str(cwd / volname), volumefaces,of='bin')
        e.stl=volname

  @classmethod
  def _mesh_single(cls, fid, vid, refine):
    f=cls.cq_mesher_entities[vid].solid.Faces()[fid]
    hh=hash(f)
    if hh in cls.cq_mesher_faceHash.keys():
      #surface is in table - use that file for this volume
      print(f'reusing {hh} {cls.cq_mesher_faceHash[hh]}')
      return(hh,cls.cq_mesher_faceHash[hh])
    else:
      facename=f'vol_{vid+1}_face{fid}.stl'
      cls.cq_mesher_faceHash[hh]=facename
      f.exportStl(facename, tolerance=cls.cq_mesher_tolerance, angularTolerance=cls.cq_mesher_ang_tolerance, ascii=True)
      if(True or self.verbose>1):
        print(f"INFO: cq export to file {facename}")
      if (refine):
        cls._refine_stls(facename,refine)
      return(hh,facename)

class MesherCQSTLBuilder:
  def __init__(self):
    self._instance = None

  def __call__(self, tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default, refine, threads, entities, **_ignored):
    if not self._instance:
      self._instance = MesherCQSTL(tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default,refine, threads, entities)
    return self._instance
