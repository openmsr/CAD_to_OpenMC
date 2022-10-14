import gmsh
import cadquery as cq
import numpy as np
import OCP
from collections.abc import Iterable
import pathlib as pl
from typing import List, Optional, Tuple, Union

import meshio
import trimesh

import tempfile
import re
import os
from pymoab import core, types

from .assemblymesher import *

mesher_config={
#general opts
  'default':False,
  'vetoed':None,
#stl-backend opts
  'tolerance':0.1,
  'angular_tolerance':0.2,
#gmsh-backend opts
  'min_mesh_size':0.1,
  'max_mesh_size':10,
  'curve_samples':20,
  'mesh_algorithm':1,
  'threads':4,
  'radial_threshold':0,
  'refine':2
}

#these are dummies that we still need to define
def _replace(filename, string1, string2):
    pass

class Entity:
    """This class is a shallow container simply to allow iterating of the geometry model
    where each object is merely a solid and a material tag.
    This is to emulate some properties of the paramak shape class, but avoid depending on that.
    At some point it should be able to import also Shapes.
    For now it can just get a solid as input
    """
    def __init__(self,solid=None,idx=0, tag='vacuum'):
        if solid is not None:
            self.solid=solid
            self.idx=idx
            self.tag=tag
            #extract some parameters from the solid for later.
            self.bb=solid.BoundingBox()
            self.center = solid.Center()
            self.volume = solid.Volume()

    def similarity(self,center:tuple=(0,0,0),bb:tuple=(0,0,0),volume:float=1,
            tolerance=1e-2)->bool:
        """method returns a value for the similary between the entity and the 3 parameters
           cms, bb, and volume"""
        cms_rel_dist= np.linalg.norm([self.center.x-center[0], self.center.y-center[1], self.center.z-center[2]])/np.linalg.norm(center)
        bb_rel_dist=np.linalg.norm([self.bb.xlen-bb[0],self.bb.ylen-bb[1],self.bb.zlen-bb[2]])/np.linalg.norm(bb)
        vol_rel_dist=np.abs(self.volume-volume)/volume
        return cms_rel_dist + bb_rel_dist + vol_rel_dist

    def similar(self,center:tuple=(0,0,0),bb:tuple=(0,0,0),volume:float=1,
            tolerance=1e-2)->bool:
        """method checks whether the entity can be regard as similar with another entities parameters"""
        cms_close=np.linalg.norm([self.center.x-center[0], self.center.y-center[1], self.center.z-center[2]])/np.linalg.norm(center)<tolerance
        bb_close=np.linalg.norm([self.bb.xlen-bb[0],self.bb.ylen-bb[1],self.bb.zlen-bb[2]])/np.linalg.norm(bb)<tolerance
        vol_close=np.abs(self.volume-volume)/volume<tolerance
        return (cms_close and bb_close and vol_close)

    def export_stp(self):
        """export the entity to a step-file using its tag as filename through cadquery export"""
        pass


def idx_similar(entity_list,center,bounding_box,volume):
    """returns the index in the solid_list for which a solid is similar in terms of bounding box, cms, and volume
       If no similar object is found return -1.
    """
    idx_found=[]
    found=False
    for i,ent in enumerate(entity_list):
      if ent.similar([center.x,center.y,center.z],[bounding_box.xlen,bounding_box.ylen, bounding_box.zlen],volume, tolerance=1e1):
        found=True
        idx_found.append(i)
    if(len(idx_found)>1):
      #we have multiple matches - pick the best one
      dsmall=1e9
      for i,ent in enumerate([entity_list[idx] for idx in idx_found]):
        d=ent.similarity([center.x,center.y,center.z],[bounding_box.xlen,bounding_box.ylen, bounding_box.zlen],volume)
        if(d<dsmall):
          dsmall=d
          end_idx=idx_found[i]
    elif(len(idx_found)==1):
        end_idx=idx_found[0]
        print(f'INFO: Found index at {end_idx}')
    else:
        end_idx=-1
        print('INFO: No similar object found')
    return end_idx

class Assembly:
    """This class encapsulates a set of geometries defined by step-files
    addtionally it provides access to meshing-utilities, and export to a DAGMC-enabled
    h5m scene, which may be used for neutronics.
    This class is based on (and borrows heavily from) the paramak package.
    """
    def __init__(self, stp_files=[], stl_files=[], verbose:int = 1, default_tag='vacuum'):
        self.stp_files=stp_files
        self.stl_files=stl_files
        self.entities=[]
        self.verbose=verbose
        self.default_tag=default_tag

    def import_stp_files(self,tags:dict=None,default_tag:str='vacuum', scale=0.1,translate=[],rotate=[]):
        tags_set=0
        #clear list to avoid double-import
        self.entities=[]

        for stp in self.stp_files:
            solid = self.load_stp_file(stp,scale,translate,rotate)

            ents=[]
            #try if solid is iterable
            try:
                for s in solid:
                    e = Entity(solid=s)
                    ents.append(e)
            except:
                e = Entity(solid=solid)
                ents.append(e)

            if(tags is None):
                #also import using gmsh to extract the material tags from the labels in the step files
                gmsh.initialize()
                vols=gmsh.model.occ.importShapes(stp)
                gmsh.model.occ.synchronize()
                for (e,v) in zip(ents,vols):
                    vid=v[1]
                    try:
                        s=gmsh.model.getEntityName(3,vid)
                        part=s.split('/')[-1]
                        g=re.match(r'^([^\s_@]+)',part)
                        tag=g[0]
                        if(self.verbose>1):
                            print(f"INFO: Tagging volume #{vid} label:{s} with material {tag}")
                    except:
                        tag=default_tag
                    e.tag=tag
                    tags_set=tags_set+1
                gmsh.finalize()
            elif (tags):
                #tag objects according to the tags dictionary.
                gmsh.initialize()
                vols=gmsh.model.occ.importShapes(stp)
                gmsh.model.occ.synchronize()
                for (e,v) in zip(ents,vols):
                    vid=v[1]
                    try:
                        s=gmsh.model.getEntityName(3,vid)
                        part=s.split('/')[-1]
                        tag=None
                        for k in tags.keys():
                            g=re.match(k,part)
                            if (g):
                                tag=tags[k]
                                break
                        if tag is None:
                            tag=self.default_tag
                        else:
                            tags_set=tags_set+1
                        if(self.verbose>1):
                            print(f"INFO: Tagging volume #{vid} label:{s} with material {tag}")
                    except:
                        tag=default_tag
                    e.tag=tag
                gmsh.finalize()

            self.entities.extend(ents)
        if(tags_set!=len(self.entities)):
           print(f"WARNING: {len(self.entities)-tags_set} volumes were tagged with the default ({default_tag}) material.")

    def load_stp_file(self,filename: str, scale_factor: float = 0.1,translate: list = [],rotate: list = []):
        """Loads a stp file and makes the 3D solid and wires available for use.

        Args:
            filename: the filename of the file containing the step-definition.
            scale_factor: a scaling factor to apply to the geometry that can be
                used to increase the size or decrease the size of the geometry.
                Useful when converting the geometry to cm for use in neutronics
                simulations. The default (0.1) corresponds to the standard setting of OpenCASCADE
                which assumes mm to the the standard of OpenMC which is cm.

        Returns:
            [CadQuery.solid]
        """
        #import _all_ the shapes in the file - i.e. may return a list
        part = cq.importers.importStep(str(filename)).vals()

        # apply apply_transforms
        scaled_part = self.apply_transforms(part,filename,scale_factor,translate,rotate)

        solid=[]
        #serialize
        #Solids() returns a list even if the part is not a Compund object
        try:
            for p in scaled_part:
                solid.extend(p.Solids())
        except:
            solid.extend(scaled_part.Solids())

        return solid

    def apply_transforms(self,part,filename,scale_factor,translate,rotate):
        #scale the shapes even if the factor is 1.
        if(self.verbose!=0):
            print(f'INFO: {str(filename)} imported - scaling')
        try:
            transformed_part = [p.scale(scale_factor) for p in part]
        except:
            transformed_part = part.scale(scale_factor)

        # translation
        if translate:
            if(self.verbose!=0):
                print(f'INFO: {str(filename)} imported - applying translation(s)')
            try:
                vols = translate[::2]
                translations = translate[1::2]
                for v in enumerate(vols):
                    if(self.verbose>1):
                        print(f"INFO: Applying translation: {translations[v[0]]} to vol(s) {v[1]}")
                    for vol in v[1]:
                        transformed_part[vol-1] = transformed_part[vol-1].translate(translations[v[0]])
            except:
                transformed_part = transformed_part.translate(translate[1])

        # rotation
        if rotate:
            if(self.verbose!=0):
                print(f'INFO: {str(filename)} imported - applying rotation(s)')
            try:
                vols = rotate[::4]
                for v in enumerate(vols):
                    if(self.verbose>1):
                        print(f"INFO: Applying rotation: {rotate[v[0]+3]} degrees about ax {rotate[v[0]+1]},{rotate[v[0]+2]} to vol(s) {v[1]}\n")
                    for vol in v[1]:
                        transformed_part[vol-1] = transformed_part[vol-1].rotate(rotate[4*v[0]+1],rotate[4*v[0]+2],rotate[4*v[0]+3])
            except:
                transformed_part = transformed_part.rotate(rotate[1],rotate[2],rotate[3])

        return transformed_part

    #export entire assembly to stp
    def export_stp(
        self,
        filename: Union[List[str], str] = None,
        mode: Optional[str] = "solid",
        units: Optional[str] = "mm",
    ) -> Union[List[str], str]:
        """Exports the assembly as a stp file or files.

        Args:
            filename: Accepts a single filename as a string which exports the
                full reactor model to a single file. Alternativley filename can
                also accept a list of strings where each string is the filename
                of the the individual shapes that make it up. This will result
                in separate files for each shape in the reactor. Defaults to
                None which uses the Reactor.name with '.stp' appended to the end
                of each entry.
            mode: the object to export can be either 'solid' which exports 3D
                solid shapes or the 'wire' which exports the wire edges of the
                shape.
            units: the units of the stp file, options are 'cm' or 'mm'.
                Default is mm.
        Returns:
            The stp filename(s) created
        """

        if isinstance(filename, str):

            # exports a single file for the whole model
            assembly = cq.Assembly(name="reactor")
            for entry in self.entities:
                if entry.color is None:
                    assembly.add(entry.solid)
                else:
                    assembly.add(entry.solid, color=cq.Color(*entry.color))

            assembly.save(filename, exportType="STEP")

            if units == "cm":
                _replace(filename, "SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)")

            return [filename]

        if filename is None:
            if None in self.name:
                msg = (
                    "Shape.name is None and therefore it can't be used "
                    "to name a stp file. Try setting Shape.name for all "
                    "shapes in the reactor"
                )
                raise ValueError(msg)
            filename = [f"{name}.stp" for name in self.name]

        # exports the assembly solid as a separate stp files
        if len(filename) != len(self.entities):
            msg = (
                f"The Assembly contains {len(self.shapes_and_components)} "
                f"Shapes and {len(filename)} filenames have be provided. "
                f"The names of the shapes are {self.name}"
            )
            raise ValueError(msg)

        for stp_filename, entry in zip(filename, self.entities):

            entry.export_stp(
                filename=stp_filename,
                mode=mode,
                units=units,
                verbose=False,
            )

            if units == "cm":
                _replace(stp_filename, "SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)")

        return filename

    #See issue 4 - we should clean up the parameter-interface to gmsh (and friends)
    def solids_to_h5m(self,brep_filename: str = None, h5m_filename:str="dagmc.h5m", samples: int =100,
            min_mesh_size:float =0.1, max_mesh_size:float =1.0,delete_intermediate_stl_files:bool=False,
            backend:str="gmsh", stl_tol:float=0.1, stl_ang_tol:float=0.2, threads:int=1, heal:bool=True, gmsh_default_opts=False):
        """calls the lower level gmsh functions in order"""

        mesher_config['entities']=self.entities
        meshgen=meshers.get(backend,**mesher_config)
        stl_list=meshgen.generate_stls()
        if(heal):
          stl_list=self.heal_stls(stl_list)
        self.stl2h5m(stl_list,h5m_filename,True)

    def tag_geometry_with_mats(self,volumes,implicit_complement_material_tag,graveyard, default_tag='vacuum'):
        """Tag all volumes with materials coming from the step files
        """
        volume_mat_list = {}
        offset=0
        #do auto-tags for all but the last volumue , which might be a graveyard
        for entry in volumes[:-1]:
            tagid=entry[1]
            try:
                s=gmsh.model.getEntityName(3,tagid)
                part=s.split('/')[-1]
                g=re.match(r'^([^\s_@]+)',part)
                volume_mat_list[tagid]=g[0]
            except:
                volume_mat_list[tagid]=default_tag
        if graveyard:
            #we have added a graveyard at the end of the list
            #tag it accordingly
            volume_mat_list[volumes[-1][1]]='Graveyard'
        else:
            #assume that the graveyard (if needed) has been added beforehand
            #use the tag from the file for the last element which may not be
            #the graveyard
            try:
                tagid=volumes[-1][1]
                s=gmsh.model.getEntityName(3,tagid)
                mat=s.split('/')
                volume_mat_list[tagid]=mat[1]
            except:
                volume_mat_list[tagid]=default_tag
        return volume_mat_list

    def stl2h5m(self,stls:list,h5m_file:str='dagmc.h5m', vtk:bool=False) -> str:
        """function that export the list of stls that we have presumably generated somehow
        and merges them into a DAGMC h5m-file by means of the MOAB-framework.
        """

        if(self.verbose>0):
            print("INFO: reassembling stl-files into h5m structure")
        h5m_p=pl.Path(h5m_file)
        moab_core,moab_tags = self.init_moab()

        sid,vid = (1,1)
        for e in self.entities:
            if (self.verbose>1):
                print(f"INFO: add stl-file \"{e.stl}\" with tag \"{e.tag}\" to MOAB structure")
            moab_core = self.add_stl_to_moab_core(moab_core,sid,vid,e.tag, moab_tags, e.stl)
            vid += 1
            sid += 1

        all_sets = moab_core.get_entities_by_handle(0)

        file_set = moab_core.create_meshset()

        moab_core.add_entities(file_set, all_sets)
        if(self.verbose>0):
            print(f"INFO: writing geometry to h5m: \"{h5m_file}\".")
        moab_core.write_file(str(h5m_p))

        self.check_h5m_file(h5m_file)

        if(vtk):
            moab_core.write_file(str(h5m_p.with_suffix('.vtk')))

        return str(h5m_p)

    def check_h5m_file(self,h5m_file:str='dagmc.h5m'):
      with open(h5m_file,"rb") as f:
        magic_bytes=f.read(8)
        if(magic_bytes!=b'\x89HDF\x0d\x0a\x1a\x0a'):
          print(f'ERROR: generated file {h5mfile} does not appear to be a hdf-file. Did you compile the moab libs with HDF enabled?')
          exit(-1)

    def add_stl_to_moab_core(self, moab_core: core.Core, surface_id: int, volume_id: int, material_name: str, tags: dict, stl_filename: str,
) -> core.Core:
        """
        Appends a set of surfaces (comprising a volume) from an stl-file to a moab.Core object and returns the updated object

        Args:
            moab_core: A moab core object
            surface_id: the id number to apply to the surface
            volume_id: the id numbers to apply to the volumes
            material_name: the material tag name to add. the value provided will
                be prepended with "mat:" unless it is "reflective" which is
                a special case and therefore will remain as is.
            tags: A dictionary of the MOAB tags
            stl_filename: the filename of the stl file to load into the moab core

        Returns:
            An updated pymoab.core.Core() instance
        """

        surface_set = moab_core.create_meshset()
        volume_set = moab_core.create_meshset()

        # recent versions of MOAB handle this automatically
        # but best to go ahead and do it manually
        moab_core.tag_set_data(tags["global_id"], volume_set, volume_id)
        moab_core.tag_set_data(tags["global_id"], surface_set, surface_id)

        # set geom IDs
        moab_core.tag_set_data(tags["geom_dimension"], volume_set, 3)
        moab_core.tag_set_data(tags["geom_dimension"], surface_set, 2)

        # set category tag values
        moab_core.tag_set_data(tags["category"], volume_set, "Volume")
        moab_core.tag_set_data(tags["category"], surface_set, "Surface")

        # establish parent-child relationship
        moab_core.add_parent_child(volume_set, surface_set)

        # Set surface sense
        # This should be fixed - we should know which volume comes next, instead of just setting it to be 0
        sense_data = [volume_set, np.uint64(0)]
        moab_core.tag_set_data(tags["surf_sense"], surface_set, sense_data)

        # load the stl triangles/vertices into the surface set
        moab_core.load_file(stl_filename, surface_set)

        group_set = moab_core.create_meshset()
        moab_core.tag_set_data(tags["category"], group_set, "Group")

        # reflective is a special case that should not have mat: in front
        if not material_name == "reflective":
            dag_material_tag = "mat:{}".format(material_name)
        else:
            dag_material_tag = material_name

        moab_core.tag_set_data(tags["name"], group_set, dag_material_tag)
        moab_core.tag_set_data(tags["geom_dimension"], group_set, 4)

        # add the volume to this group set
        moab_core.add_entity(group_set, volume_set)

        return moab_core

    def init_moab(self):
        """Creates a MOAB Core instance which can be built up by adding sets of
        triangles to the instance
        Returns:
        (pymoab Core): A pymoab.core.Core() instance
        (pymoab tag_handle): A pymoab.core.tag_get_handle() instance
        """
        # create pymoab instance
        moab_core = core.Core()

        tags = dict()
        sense_tag_name = "GEOM_SENSE_2"
        sense_tag_size = 2
        tags["surf_sense"] = moab_core.tag_get_handle(sense_tag_name, sense_tag_size,
            types.MB_TYPE_HANDLE, types.MB_TAG_SPARSE, create_if_missing=True,
        )
        tags["category"] = moab_core.tag_get_handle(
            types.CATEGORY_TAG_NAME, types.CATEGORY_TAG_SIZE, types.MB_TYPE_OPAQUE,
            types.MB_TAG_SPARSE, create_if_missing=True,
        )
        tags["name"] = moab_core.tag_get_handle(
            types.NAME_TAG_NAME, types.NAME_TAG_SIZE,
            types.MB_TYPE_OPAQUE, types.MB_TAG_SPARSE, create_if_missing=True,
        )
        tags["geom_dimension"] = moab_core.tag_get_handle(
            types.GEOM_DIMENSION_TAG_NAME, 1,
            types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, create_if_missing=True,
        )
        # Global ID is a default tag, just need the name to retrieve
        tags["global_id"] = moab_core.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
        return moab_core, tags

    def merge_all(self):
        #merging a single object does not really make sense
        if len(self.entities)>1:
          #extract cq solids backend algorithm
          unmerged=[e.solid for e in self.entities]
          #do merge
          merged=self._merge_solids(unmerged, fuzzy_value=1e-2)
          #the merging process may result in extra volumes.
          #We need to make sure these are at the end of the list
          #If not this results in a loss of volumes in the end.
          print("INFO: reordering volumes")
          tmp_ents=[]
          for solid in merged.Solids():
            center=solid.Center()
            bb=solid.BoundingBox()
            vol=solid.Volume()
            idx=idx_similar(self.entities,center,bb,vol)
            ent=self.entities[idx]
            #have to use the newly created solid to get the merged entries instead.
            ent.solid=solid
            tmp_ents.append(ent)
          self.entities=tmp_ents

    def _merge_solids(self,solids,fuzzy_value):
        """merge a set of cq-solids
           returns as cq-compound object
        """
        bldr = OCP.BOPAlgo.BOPAlgo_Splitter()
        bldr.SetFuzzyValue(fuzzy_value)
        #loop trough all objects in geometry and split and merge them accordingly
        #shapes should be a compund cq object or a list thereof
        for shape in solids:
          # checks if solid is a compound as .val() is not needed for compunds
          if isinstance(shape, cq.occ_impl.shapes.Compound):
            bldr.AddArgument(shape.wrapped)
          else:
              try:
                  bldr.AddArgument(shape.val().wrapped)
              except:
                  bldr.AddArgument(shape.wrapped)
        bldr.SetParallelMode_s(True)
        bldr.SetNonDestructive(False)

        if(self.verbose>1):
            print("INFO: Commence perform step of merge")
        bldr.Perform()

        if(self.verbose>1):
            print("INFO: Commence image step of merge")
        bldr.Images()

        if(self.verbose>1):
            print("INFO: Generate compound shape")
        merged = cq.Compound(bldr.Shape())

        return merged

    def heal_stls(self,stls):
        if(self.verbose>0):
            print("INFO: checking surfaces and reparing normals")

        healed=[]
        for e in self.entities:
            stl=e.stl
            mesh = trimesh.load_mesh(stl)
            if (self.verbose>1):
                print("INFO: stl-file", stl, ": mesh is watertight", mesh.is_watertight)
            trimesh.repair.fix_normals(
                mesh
            )  # reqired as gmsh stl export from brep can get the inside outside mixed up
            new_filename = stl[:-4] + "_with_corrected_face_normals.stl"
            mesh.export(new_filename)
            e.stl=new_filename

    def tag_stls(self,stls):
        stl_tagged=[]
        for (stl,e) in zip(stl_list,self.entities):
            try:
                stl_tagged.append((stl[0],stl[1],e.tag))
            except:
                print("WARNING: list of material tags is exhausted. Tagging volume {stl[0]},{stl[1]} with \'vacuum\'")
                stl_tagged.append(stl[0],stl[1],'vacuum')

###############
    def export_brep(self, filename: str, merge: bool = True, step: bool =False):
        """Exports a brep file for the Assembly
        This requires serializing the assembly

        Args:
            filename: the filename of exported the brep file.
            merged: if the surfaces should be merged (True) or not (False).
            step: Use step fileformat instead of the default brep

        Returns:
            filename of the brep created. Also stored in self.brep_filename
        """

        path_filename = pl.Path(filename)
        self.brep_filename=filename

        if not step and path_filename.suffix != ".brep":
            msg = "When exporting a brep file the filename must end with .brep"
            raise ValueError(msg)
        elif step and path_filename.suffix not in (".step",".stp"):
            msg = "When exporting a step file the filename must end with .step or .stp"
            raise ValueError(msg)
        path_filename.parents[0].mkdir(parents=True, exist_ok=True)

        if not merge or len(self.entities)<=1:
            #merging a single volume does not make sense
            #If we don't merge, construct a compund solid for later.
            self.merged=cq.Compound.makeCompound([e.solid for e in self.entities])
            rval=self.merged.exportBrep(str(path_filename))
        else:
            #the merge surface returns a cq-compound object.
            self.merged = self.merge_surfaces()
            rval=self.merged.exportBrep(str(path_filename))
            #the merging process may result in extra volumes.
            #We need to make sure these are at the end of the list
            #If not this results in a loss ofvolumes in the end.
            print("INFO: reordering volumes")
            idxs=[]
            for solid in self.merged.Solids():
              center=solid.Center()
              bb=solid.BoundingBox()
              vol=solid.Volume()
              idx=idx_similar(self.entities,center,bb,vol)
              idxs.append(idx)
            #reorder
            ents=[self.entities[i] for i in idxs if i!=-1]
            self.entities=ents

        return rval

    def merge_two(self,solid1,solid2):
        """ Checks two surfaces if their BB overlap. If so merge the two - and return a list of the results
        """
