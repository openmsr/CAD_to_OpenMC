import cadquery as cq
import subprocess as sp
import pathlib as pl
import hashlib as hl
from CAD_to_OpenMC import cdtemp
from .assemblymesher_base import assemblymesher

single_thread_override=True
try:
  import multiprocessing.pool as mp
except:
  single_thread_override=True

from .stl_utils import *
from .meshutils import *

#following pattern here: https://thelaziestprogrammer.com/python/multiprocessing-pool-a-global-solution
#to enable multiprocessing meshings

class MesherCQSTL(assemblymesher):
  #these need to be class attributes to avoid pickling when spawning a multiprocessing-pool
  cq_mesher_entities=None
  cq_mesher_tolerance=None
  cq_mesher_ang_tolereance=None
  cq_mesher_min_mesh_size=None
  cq_mesher_max_mesh_size=None

  cq_mesher_faceHash={}

  def __init__(self, tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default, refine, threads, entities):
    self._set_meshpars(tolerance,angular_tolerance,min_mesh_size,max_mesh_size)
    self._clear_face_hashtable()
    self.refine=refine
    self._set_entities(entities)
    self.default=default
    self.min_mesh_size=min_mesh_size
    self.max_mesh_size=max_mesh_size
    self.threads=threads

  @property
  def refine(self):
    return self._refine

  @refine.setter
  def refine(self,ref):
    if(ref==True or ref!=0):
      self._refine=True
    else:
      self._refine=False

  @classmethod
  def _set_meshpars(cls,tol,ang_tol, min_sz, max_sz):
    cls.cq_mesher_tolerance=tol
    cls.cq_mesher_ang_tolerance=ang_tol
    cls.cq_mesher_min_mesh_size=min_sz
    cls.cq_mesher_max_mesh_size=max_sz

  @classmethod
  def _clear_face_hashtable(cls):
    cls.cq_mesher_faceHash={}

  @classmethod
  def _set_entities(cls,entities):
    cls.cq_mesher_entities=entities

  def generate_stls(self):
    stls=self._mesh_surfaces()
    return stls

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
    edges=np.array(find_edges(tris))
    write_dotmesh(str(stlp.with_suffix('.mesh')),verts,tris,edges=edges,required_edges='all')
    cp=sp.run(['mmgs_O3','-hmin',f'{cls.cq_mesher_min_mesh_size}','-hmax',f'{cls.cq_mesher_max_mesh_size}','-optim','-in',stlp.with_suffix('.mesh'),'-out',stlp.with_suffix('.o.mesh')], capture_output=True)
    if(cls.verbosity_level and cls.verbosity_level>1):
      print(cp.stdout.decode())

    import gmsh
    gmsh.initialize()
    gmsh.open(str(stlp.with_suffix('.o.mesh')))
    gmsh.write( str(stlp) )
    gmsh.finalize()

  def _mesh_surfaces(self):
    #loop over all surfaces in all entities
    #and generate meshes (refined or otherwise)
    stls=[]
    cwd=pl.Path.cwd()
    #create a workplace in tmp
    with cdtemp() as mngr:
      for i,e in enumerate(self.cq_mesher_entities):
        volname= f"volume_{i+1}.stl"
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
        stls.append(volname)
    # clear the hash table
    self._clear_face_hashtable()
    return stls

  @classmethod
  def _mesh_single(cls, fid, vid, refine):
    f=cls.cq_mesher_entities[vid].solid.Faces()[fid]
    hh=hash(f)
    if hh in cls.cq_mesher_faceHash.keys():
      #surface is in table - use that file for this volume
      if (cls.verbosity_level):
        print(f'INFO: mesher reusing {hh} {cls.cq_mesher_faceHash[hh]}')
      return(hh,cls.cq_mesher_faceHash[hh])
    else:
      facename=f'vol_{vid+1}_face{fid}.stl'
      cls.cq_mesher_faceHash[hh]=facename
      f.exportStl(facename, tolerance=cls.cq_mesher_tolerance, angularTolerance=cls.cq_mesher_ang_tolerance, ascii=True)
      if(cls.verbosity_level>1):
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
    else:
      #We are reusing a mesher instance. Hence reset the parameters and clear the hashtable
      self._instance._set_entities(entities)
      self._instance._set_meshpars(tolerance, angular_tolerance, min_mesh_size, max_mesh_size)
      self._instance._clear_face_hashtable()
    return self._instance
