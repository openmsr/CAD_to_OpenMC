import cadquery as cq
import subprocess as sp
import pathlib as pl
from .assemblymesher_base import assemblymesher

single_thread_override=False
try:
  import multiprocessing as mp
  manager = mp.Manager()
  lock = manager.Lock()
except:
  single_thread_override=True

from .stl_utils import *
from . import meshutils

class MesherCQSTL2(assemblymesher):
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
    if(len(entities)<=threads):
      self.threads=len(entities)
    else:
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
    return self._mesh_surfaces()

  def _mesh_surfaces(self):
    #loop over all surfaces in all entities
    #and generate meshes (refined or otherwise)
    mpargs=[]
    if (single_thread_override or self.threads==1):
      face_hash_table={}
    else:
      #manager=mp.Manager()
      face_hash_table=manager.dict()
    k=0
    for i,e in enumerate(self.cq_mesher_entities):
      mpargs.extend([(k+j,j,i,self.refine,hash(f), face_hash_table) for j,f in enumerate(e.solid.Faces())])

    #we have a set of mesh jobs - scatter those
    if (single_thread_override or self.threads==1):
      output=[]
      for args in mpargs:
        output.append(self._mesh_single_nothread(*args))
    else:
      pool=mp.Pool(processes=self.threads)
      output=pool.starmap(self._mesh_single, mpargs)

    #process the list of meshed faces.
    stls=[]
    for i,e in enumerate(self.cq_mesher_entities):
      face_stls=[]
      for k,v in face_hash_table.items():
        vids=v[1] # the volumes that this face belongs to
        if i in vids:
          #this face is part of this volume
          face_stls.append(v)
      stls.append(face_stls)
    return stls

  @classmethod
  def _mesh_single_nothread(cls, global_fid, fid, vid, refine, hh, faceHash):
    f=cls.cq_mesher_entities[vid].solid.Faces()[fid]
    if hh in faceHash.keys():
      #surface is in table - simply add the vid to the hash-table
      faceHash[hh][1].append(vid)
      if (cls.verbosity_level):
        print(f'INFO: mesher reusing {hh} {faceHash[hh][1]}')
      return(hh,faceHash[hh])
    else:
      facefilename=f'vol_{vid+1}_face{global_fid:04}.stl'
      faceHash[hh]=[facefilename,manager.list([vid])]
      f.exportStl(facefilename, tolerance=cls.cq_mesher_tolerance, angularTolerance=cls.cq_mesher_ang_tolerance, ascii=True)
      if(cls.verbosity_level>1):
        print(f"INFO: cq export to file {facefilename}")
      if (refine):
        cls._refine_stls(facefilename,refine)
      return(hh,faceHash[hh])

  @classmethod
  def _mesh_single(cls, global_fid, fid, vid, refine, hh, faceHash):
    f=cls.cq_mesher_entities[vid].solid.Faces()[fid]
    if hh in faceHash.keys():
      #surface is in table - simply add the vid to the hash-table
      with lock:
        faceHash[hh][1].append(vid)
      if (cls.verbosity_level):
        print(f'INFO: mesher reusing {hh} {faceHash[hh][1]}')
      return(hh,faceHash[hh])
    else:
      facefilename=f'vol_{vid+1}_face{global_fid:04}.stl'
      with lock:
        faceHash[hh]=[facefilename,manager.list([vid])]
      f.exportStl(facefilename, tolerance=cls.cq_mesher_tolerance, angularTolerance=cls.cq_mesher_ang_tolerance, ascii=True)
      if(cls.verbosity_level>1):
        print(f"INFO: cq export to file {facefilename}")
      if (refine):
        cls._refine_stls(facefilename,refine)
      return(hh,faceHash[hh])

class MesherCQSTL2Builder:
  def __init__(self):
    self._instance = None

  def __call__(self, tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default, refine, threads, entities, **_ignored):
    if not self._instance:
      self._instance = MesherCQSTL2(tolerance, angular_tolerance, min_mesh_size, max_mesh_size, default,refine, threads, entities)
    else:
      #We are reusing a mesher instance. Hence reset the parameters and clear the hashtable
      self._instance._set_entities(entities)
      self._instance._set_meshpars(tolerance, angular_tolerance, min_mesh_size, max_mesh_size)
      self._instance._clear_face_hashtable()
    return self._instance
