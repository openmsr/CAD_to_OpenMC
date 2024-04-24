import cadquery as cq
import subprocess as sp
import pathlib as pl
import OCP
from .assemblymesher_base import assemblymesher
import os

from .stl_utils import *
from . import meshutils

class MesherCQSTL2(assemblymesher):
    # these need to be class attributes to avoid pickling when spawning a multiprocessing-pool
    cq_mesher_entities = None
    cq_mesher_tolerance = None
    cq_mesher_ang_tolereance = None
    cq_mesher_min_mesh_size = None
    cq_mesher_max_mesh_size = None

    cq_mesher_faceHash = {}

    #writer object from OCP
    wr=OCP.StlAPI.StlAPI_Writer()
    wr.ASCIIMode=True

    def __init__(
        self,
        tolerance,
        angular_tolerance,
        min_mesh_size,
        max_mesh_size,
        default,
        refine,
        threads,
        entities,
    ):
        self._set_meshpars(tolerance, angular_tolerance, min_mesh_size, max_mesh_size)
        self._clear_face_hashtable()
        self.refine = refine
        self._set_entities(entities)
        self.default = default
        self.min_mesh_size = min_mesh_size
        self.max_mesh_size = max_mesh_size
        if len(entities) <= threads:
            self.threads = len(entities)
        else:
            self.threads = threads

    @property
    def refine(self):
        return self._refine

    @refine.setter
    def refine(self, ref):
        if ref == True or ref != 0:
            self._refine = True
        else:
            self._refine = False

    @classmethod
    def _set_meshpars(cls, tol, ang_tol, min_sz, max_sz):
        cls.cq_mesher_tolerance = tol
        cls.cq_mesher_ang_tolerance = ang_tol
        cls.cq_mesher_min_mesh_size = min_sz
        cls.cq_mesher_max_mesh_size = max_sz

    @classmethod
    def _clear_face_hashtable(cls):
        cls.cq_mesher_faceHash = {}

    @classmethod
    def _set_entities(cls, entities):
        cls.cq_mesher_entities = entities

    def generate_stls(self):
        return self._mesh_surfaces()

    @classmethod
    def surface_hash(self,surface):
        part1=surface
        part2=surface.Center().toTuple()
        return hash( (part1,part2) )

    def _mesh_surfaces(self):
        # loop over all surfaces in all entities
        # and generate meshes (refined or otherwise)
        mpargs = []
        face_hash_table={}

        k = 0
        for i, e in enumerate(self.cq_mesher_entities):
            if self.verbosity_level:
                print(f"INFO: triangulating solid {i}")
            e.solid=self._triangulate_solid(e.solid,self.cq_mesher_tolerance,self.cq_mesher_ang_tolerance)
            mpargs.extend(
                [
                    [k + j, j, i, self.refine, self.surface_hash(f), face_hash_table]
                    for j, f in enumerate(e.solid.Faces())
                ]
            )

        # we have a set of mesh jobs - scatter those
        for args in mpargs:
            self._mesh_single(*args)

        # process the list of meshed faces.
        stls = []
        for i, e in enumerate(self.cq_mesher_entities):
            face_stls = []
            for k, v in face_hash_table.items():
                vids = v[1:]  # the volumes that this face belongs to
                if i in vids:
                    # this face is part of this volume
                    face_stls.append([v[0],v[1:]])
            stls.append(face_stls)
        return stls

    def _triangulate_solid(self, solid, tol: float = 1e-3, atol: float = 1e-1):
        """ create a mesh by means of the underlying OCCT IncrementalMesh
            on a single solid. This will later be split into surfaces.
            This has to be done since otherwise a single solid can get leaky
            when its surfaces do not connect properly
        """
        solid.mesh(tol,atol)
        return solid

    @classmethod
    def _mesh_single(cls, global_fid, fid, vid, refine, hh, faceHash):
        f = cls.cq_mesher_entities[vid].solid.Faces()[fid]
        if hh in faceHash.keys():
            # surface is in table - simply add the vid to the hash-table
            done=True
            ffn=faceHash[hh][0]
            previd=faceHash[hh][1]
            faceHash[hh]=[ffn,previd,vid]
        else:
            done=False

        if done:
            if cls.verbosity_level:
                print(f"INFO: mesher reusing {hh} ({faceHash[hh][0]},{faceHash[hh][1:]})")
            return
        else:
            facefilename = f"vol_{vid+1}_face{global_fid:04}.stl"
            faceHash[hh] = [facefilename, vid, -1]
            status=cls.wr.Write(f.wrapped,facefilename)
            if cls.verbosity_level > 1:
                print(f"INFO: cq export to file {facefilename}")
            if refine:
                cls._refine_stls(facefilename, refine)
            return

class MesherCQSTL2Builder:
    def __init__(self):
        self._instance = None

    def __call__(
        self,
        tolerance,
        angular_tolerance,
        min_mesh_size,
        max_mesh_size,
        default,
        refine,
        threads,
        entities,
        **_ignored,
    ):
        if not self._instance:
            self._instance = MesherCQSTL2(
                tolerance,
                angular_tolerance,
                min_mesh_size,
                max_mesh_size,
                default,
                refine,
                threads,
                entities,
            )
        else:
            # We are reusing a mesher instance. Hence reset the parameters and clear the hashtable
            self._instance._set_entities(entities)
            self._instance._set_meshpars(
                tolerance, angular_tolerance, min_mesh_size, max_mesh_size
            )
            self._instance._clear_face_hashtable()
        return self._instance
