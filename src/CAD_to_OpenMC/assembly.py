import cadquery as cq
import numpy as np
import OCP
import pathlib as pl
from typing import List, Optional, Union

from itertools import zip_longest

import trimesh

import re
import os
import math
from datetime import datetime

from pymoab import core, types

import CAD_to_OpenMC.assemblymesher as am
from CAD_to_OpenMC.datadirectory import mesher_datadir

try:
    import gmsh

    nogmsh = False
except ImportError as e:
    print(f"Warning: import gmsh failed ({e})- material tag-list, must be supplied.")
    nogmsh = True

mesher_config = {
    # general opts
    "default": False,
    "vetoed": None,
    # stl-backend opts
    "tolerance": 0.1,
    "angular_tolerance": 0.2,
    # gmsh-backend opts
    "min_mesh_size": 0.1,
    "max_mesh_size": 10,
    "curve_samples": 20,
    "mesh_algorithm": 1,
    "threads": 4,
    "radial_threshold": 0,
    "refine": 0,
    "verbose": 0,
}


# these are dummies that we still need to have defined
def _replace(filename, string1, string2):
    pass


class Entity:
    """This class is a shallow container simply to allow iterating of the geometry model
    where each object is merely a solid and a material tag.
    This is to emulate some properties of the paramak shape class, but avoid depending on that.
    At some point it should be able to import also Shapes.
    For now it can just get a solid as input
    """

    def __init__(self, solid=None, idx=0, tag="vacuum"):
        if solid is not None:
            self.solid = solid
            self.idx = idx
            self.tag = tag
            # extract some parameters from the solid for later.
            self.bb = solid.BoundingBox()
            self.center = solid.Center()
            self.volume = solid.Volume()

    def similarity(
        self,
        center: tuple = (0, 0, 0),
        bb: tuple = (0, 0, 0),
        volume: float = 1,
        tolerance=1e-2,
    ) -> bool:
        """method returns a value for the similarity between the entity and the 3 parameters
        cms, bb, and volume"""
        cms_rel_dist = np.linalg.norm(
            [
                self.center.x - center[0],
                self.center.y - center[1],
                self.center.z - center[2],
            ]
        ) / np.linalg.norm(center)
        bb_rel_dist = np.linalg.norm(
            [self.bb.xlen - bb[0], self.bb.ylen - bb[1], self.bb.zlen - bb[2]]
        ) / np.linalg.norm(bb)
        vol_rel_dist = np.abs(self.volume - volume) / volume
        return cms_rel_dist + bb_rel_dist + vol_rel_dist

    def similar(
        self,
        center: tuple = (0, 0, 0),
        bb: tuple = (0, 0, 0),
        volume: float = 1,
        tolerance=1e-2,
    ) -> bool:
        """method checks whether the entity can be regard as similar with another entities parameters"""
        cms_close = (
            np.linalg.norm(
                [
                    self.center.x - center[0],
                    self.center.y - center[1],
                    self.center.z - center[2],
                ]
            )
            / np.linalg.norm(center)
            < tolerance
        )
        bb_close = (
            np.linalg.norm(
                [self.bb.xlen - bb[0], self.bb.ylen - bb[1], self.bb.zlen - bb[2]]
            )
            / np.linalg.norm(bb)
            < tolerance
        )
        vol_close = np.abs(self.volume - volume) / volume < tolerance
        return cms_close and bb_close and vol_close

    def export_stp(self):
        """export the entity to a step-file using its tag as filename through cadquery export"""
        pass


def idx_similar(entity_list, center, bounding_box, volume):
    """returns the index in the solid_list for which a solid is similar in terms of bounding box, cms, and volume
    If no similar object is found return -1.
    """
    idx_found = []
    for i, ent in enumerate(entity_list):
        if ent.similar(
            [center.x, center.y, center.z],
            [bounding_box.xlen, bounding_box.ylen, bounding_box.zlen],
            volume,
            tolerance=1e1,
        ):
            idx_found.append(i)
    if len(idx_found) > 1:
        # we have multiple matches - pick the best one
        dsmall = 1e9
        for i, ent in enumerate([entity_list[idx] for idx in idx_found]):
            d = ent.similarity(
                [center.x, center.y, center.z],
                [bounding_box.xlen, bounding_box.ylen, bounding_box.zlen],
                volume,
            )
            if d < dsmall:
                dsmall = d
                end_idx = idx_found[i]
    elif len(idx_found) == 1:
        end_idx = idx_found[0]
        print(f"INFO: Found index at {end_idx}")
    else:
        end_idx = -1
        print("INFO: No similar object found")
    return end_idx


def similar_solids(solid1_vol, solid1_bb, solid1_c, solid2_vol, solid2_bb, solid2_c):
    """This function compares two solids and reports their similarity constant.
    defined as the sum of:
    1. cubic root difference in volume
    2. difference of bounding box diagonal
    3. difference in location vector.
    """
    dV = math.pow(
        math.fabs(solid1_vol - solid2_vol), 0.3333333333333333333333333333333333
    )
    dBB = math.fabs(solid1_bb.DiagonalLength - solid2_bb.DiagonalLength)
    c1 = solid1_c
    c2 = solid2_c
    dCntr = math.sqrt(
        (c1.x - c2.x) * (c1.x - c2.x)
        + (c1.y - c2.y) * (c1.y - c2.y)
        + (c1.z - c2.z) * (c1.z - c2.z)
    )
    return dV + dBB + dCntr


class Assembly:
    """This class encapsulates a set of geometries defined by step-files
    addtionally it provides access to meshing-utilities, and export to a DAGMC-enabled
    h5m scene, which may be used for neutronics.
    This class is based on (and borrows heavily from) logic from the paramak package.
    """

    def __init__(
        self,
        stp_files=[],
        stl_files=[],
        verbose: int = 1,
        default_tag="vacuum",
        implicit_complement=None,
    ):
        self.stp_files = stp_files
        self.stl_files = stl_files
        self.entities = []
        self.verbose = verbose
        self.default_tag = default_tag
        self.delete_intermediate = False
        self.cleanup = False
        self.datadir="."
        self.tags = None
        self.sequential_tags = None
        self.implicit_complement = implicit_complement

    @classmethod
    def hdf5_in_moab(cls):
        a = Assembly()
        try:
            a.dummy_h5m()
        except RuntimeError as e:
            print(
                "Warning: Can't write a hdf-file. Did you compile moab with hdf5 support?"
                "The resulting meshed file will actually be a vtk-file."
                f"{e}"
            )
        return True

    def dummy_h5m(self):
        mbc, mbt = self.init_moab()
        mbc.write_file("dummy.h5m")
        try:
            self.check_h5m_file("dummy.h5m")
        except RuntimeError as e:
            raise e
        os.unlink("dummy.h5m")
        return True

    def run(
        self,
        backend: str = "stl",
        h5m_filename: str = "dagmc.h5m",
        merge: bool = True,
        imprint: bool = False,
        **kwargs: dict,
    ):
        """convenience function that assumes the stp_files field is set, etc and simply runs the mesher with the set options"""
        self.import_stp_files(**kwargs)
        if merge:
            self.merge_all()
        elif imprint:
            self.imprint_all()
        self.solids_to_h5m(backend=backend, h5m_filename=h5m_filename, **kwargs)

    def import_stp_files(
        self,
        tags: dict = None,
        sequential_tags: iter = None,
        match_anywhere: bool = False,
        default_tag: str = "vacuum",
        scale: float = 0.1,
        translate: iter = [],
        rotate: iter = [],
        vol_skip: iter = [],
        **kwargs: dict,
    ):
        """
        Import a list of step-files.

        Args:
          tags: dictionary containing pairs of reg.exp. patterns and material tags. If not None, entities with
              names matching the patterns will be assigned to corresponding tag. If no patterns match
              the default_tag will be applied
          match_anywhere: match patterns anywhere in the entitiy name
          default_tag: The material tag that will be applied if no patterns are matched.
          scale: overall scaling factor applied to all parts
          translate: Translation vector to apply to all parts in the step-file.
          rotate: Rotation angles to apply to the parts in the step-file.
        """

        tags_set = 0
        # clear list to avoid double-import
        self.entities = []

        message = "Need gmsh python module installed to extract material tags from step-file. please supply a 'sequential_tags'-list instead"
        assert (nogmsh and tags is None and sequential_tags is not None) or (
            not nogmsh
        ), message
        # if no gmsh module was imported we must rely on explicit sequential tags, so check they're there.
        i = 1
        for stp in self.stp_files:
            solid = self.load_stp_file(stp, scale, translate, rotate)

            ents = []
            # try if solid is iterable
            try:
                for j in range(len(solid)):
                    if j + i not in vol_skip:
                        s = solid[j]
                        e = Entity(solid=s)
                        ents.append(e)
                i = i + len(solid)
            except Exception as _e:
                e = Entity(solid=solid)
                ents.append(e)

            if tags is None and sequential_tags is None:
                # also import using gmsh to extract the material tags from the labels in the step files
                gmsh.initialize()
                vols = gmsh.model.occ.importShapes(stp)
                gmsh.model.occ.synchronize()
                for e, v in zip(ents, vols):
                    vid = v[1]
                    try:
                        s = gmsh.model.getEntityName(3, vid)
                        part = s.split("/")[-1]
                        g = re.match(r"^([^\s_@]+)", part)
                        tag = g[0]
                        if self.verbose > 1:
                            print(
                                f"INFO: Tagging volume #{vid} label:{s} with material {tag}"
                            )
                    except Exception as _e:
                        tag = default_tag
                    e.tag = tag
                    tags_set = tags_set + 1
                gmsh.finalize()
            elif tags:
                # tag objects according to the tags dictionary.
                gmsh.initialize()
                vols = gmsh.model.occ.importShapes(stp)
                gmsh.model.occ.synchronize()
                for e, v in zip(ents, vols):
                    vid = v[1]
                    try:
                        s = gmsh.model.getEntityName(3, vid)
                        part = s.split("/")[-1]
                        tag = None
                        for k in tags.keys():
                            if match_anywhere:
                                g = re.search(k, part)
                            else:
                                g = re.match(k, part)
                            if g:
                                tag = tags[k]
                                break
                        if tag is None:
                            tag = self.default_tag
                        else:
                            tags_set = tags_set + 1
                        if self.verbose > 1:
                            print(
                                f"INFO: Tagging volume #{vid} label:{s} with material {tag}"
                            )
                    except Exception as _e:
                        tag = default_tag
                    e.tag = tag
                gmsh.finalize()
            elif sequential_tags:
                for e, t in zip_longest(
                    ents, sequential_tags[tags_set:], fillvalue=self.default_tag
                ):
                    # Apply tags to the imported volumes in the sequence they get imported.
                    if e == default_tag:
                        # this means we have exhausted the ents list
                        break
                    e.tag = t
                tags_set += len(ents)

            self.entities.extend(ents)
        if tags_set != len(self.entities):
            print(
                f"WARNING: {len(self.entities)-tags_set} volumes were tagged with the default ({default_tag}) material."
            )

    def load_stp_file(
        self,
        filename: str,
        scale_factor: float = 0.1,
        translate: list = [],
        rotate: list = [],
    ):
        """Loads a stp file and makes the 3D solid and wires available for use.

        Args:
            filename: the filename of the file containing the step-definition.
            scale_factor: a scaling factor to apply to the geometry that can be
                used to increase the size or decrease the size of the geometry.
                Useful when converting the geometry to cm for use in neutronics
                simulations. The default (0.1) corresponds to the standard setting of OpenCASCADE
                which assumes mm to the the standard of OpenMC which is cm.
            translate: optional, iterable of 3 float iterable. Translation vector to
                apply to parts in the model before meshing occurs. If list of list then
                each translation in the list is applied to 1 part in sequence.

        Returns:
            [CadQuery.solid]
        """
        # import _all_ the shapes in the file - i.e. may return a list
        part = cq.importers.importStep(str(filename)).vals()

        # apply apply_transforms
        scaled_part = self.apply_transforms(
            part, filename, scale_factor, translate, rotate
        )

        solid = []
        # serialize
        # Solids() returns a list even if the part is not a Compund object
        try:
            for p in scaled_part:
                solid.extend(p.Solids())
        except Exception as _e:
            solid.extend(scaled_part.Solids())

        return solid

    def apply_transforms(self, part, filename, scale_factor, translate, rotate):
        # scale the shapes even if the factor is 1.
        if self.verbose != 0:
            print(f"INFO: {str(filename)} imported - scaling")
        try:
            transformed_part = [p.scale(scale_factor) for p in part]
        except Exception as _e:
            transformed_part = part.scale(scale_factor)

        # translation
        if translate:
            if self.verbose != 0:
                print(f"INFO: {str(filename)} imported - applying translation(s)")
            try:
                vols = translate[::2]
                translations = translate[1::2]
                for v in enumerate(vols):
                    if self.verbose > 1:
                        print(
                            f"INFO: Applying translation: {translations[v[0]]} to vol(s) {v[1]}"
                        )
                    for vol in v[1]:
                        transformed_part[vol - 1] = transformed_part[vol - 1].translate(
                            translations[v[0]]
                        )
            except Exception as _e:
                transformed_part = transformed_part.translate(translate[1])

        # rotation
        if rotate:
            if self.verbose != 0:
                print(f"INFO: {str(filename)} imported - applying rotation(s)")
            try:
                vols = rotate[::4]
                for v in enumerate(vols):
                    if self.verbose > 1:
                        print(
                            f"INFO: Applying rotation: {rotate[v[0]+3]} degrees about ax {rotate[v[0]+1]},{rotate[v[0]+2]} to vol(s) {v[1]}\n"
                        )
                    for vol in v[1]:
                        transformed_part[vol - 1] = transformed_part[vol - 1].rotate(
                            rotate[4 * v[0] + 1],
                            rotate[4 * v[0] + 2],
                            rotate[4 * v[0] + 3],
                        )
            except Exception as _e:
                transformed_part = transformed_part.rotate(
                    rotate[1], rotate[2], rotate[3]
                )

        return transformed_part

    # export entire assembly to stp
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
                _replace(
                    filename, "SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)"
                )

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
                _replace(
                    stp_filename, "SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)"
                )

        return filename

    def print_summary(self):
        # output a summary of the meshing results
        print(f'SUMMARY: {"solid_id":8} {"material_tag":16} {"stl-file":16}')
        for i, a in zip(range(len(self.entities)), self.entities):
            if not isinstance(a.stl, str):
                print(
                    f"SUMMARY: {i+1:8} {a.tag:16} "
                    + " ".join([f"{stl[0]:16}" for stl in a.stl])
                )
            else:
                print(f"SUMMARY: {i+1:8} {a.tag:16} {a.stl:16}")

    def _datadir_name(self,h5m_filename=""):
        h5mf=pl.Path(h5m_filename)
        if self.datadir==".":
            self.datadir=datetime.now().strftime(f"{h5mf.stem}_%Y%m%d_%H%M%S.%f")
        if (self.verbose):
            print(f"INFO: storing temporary data in directory: {self.datadir}")
        return self.datadir

    def solids_to_h5m(
        self,
        brep_filename: str = None,
        h5m_filename: str = "dagmc.h5m",
        samples: int = 100,
        backend: str = "gmsh",
        heal: bool = True,

        **kwargs: dict,
    ):
        with mesher_datadir(self._datadir_name(h5m_filename),self.cleanup) as datadir:
            mesher_config["entities"] = self.entities
            meshgen = am.meshers.get(backend, **mesher_config)
            meshgen.set_verbosity(self.verbose)
            stl_list = meshgen.generate_stls()

            for e, s in zip(self.entities, stl_list):
                e.stl = s
            if self.verbose:
                self.print_summary()
            if backend == "stl2":
                self.stl2h5m_byface(h5m_filename, True)
            else:
                if heal:
                    stl_list = self.heal_stls(stl_list)
                self.stl2h5m(stl_list, h5m_filename, True)
            #self.cleanup()

    def stl2h5m_byface(self, h5m_file: str = "dagmc.h5m", vtk: bool = False) -> str:
        """create a h5m-file with a moab structure and fills
        it with the dagmc structure using the pymoab framework. Optionally creates
        a vtk-file.
        """
        if self.verbose > 0:
            print("INFO: reassembling stl-files into h5m structure")
        h5m_p = pl.Path(h5m_file)
        mbcore, mbtags = self.init_moab()
        mbcore = self.add_entities_to_moab_core(mbcore, mbtags)

        all_sets = mbcore.get_entities_by_handle(0)
        file_set = mbcore.create_meshset()

        mbcore.add_entities(file_set, all_sets)

        if self.verbose > 0:
            print(f'INFO: writing geometry to h5m: "{h5m_file}".')
        mbcore.write_file(str(h5m_p))

        self.check_h5m_file(h5m_file)
        if vtk:
            if self.verbose > 0:
                print(f'INFO: writing geometry to vtk: ' + str(h5m_p.with_suffix(".vtk")) )
            mbcore.write_file(str(h5m_p.with_suffix(".vtk")))

        return str(h5m_p)

    def stl2h5m(
        self, stls: list, h5m_file: str = "dagmc.h5m", vtk: bool = False
    ) -> str:
        """function that exports the list of stls that we have presumably generated somehow
        and merges them into a DAGMC h5m-file by means of the MOAB-framework.
        """

        if self.verbose > 0:
            print("INFO: reassembling stl-files into h5m structure")
        h5m_p = pl.Path(h5m_file)
        moab_core, moab_tags = self.init_moab()

        sid, vid = (1, 1)
        is_last=False
        for e in self.entities:
            if self.verbose > 1:
                print(
                    f'INFO: add stl-file "{e.stl}" with tag "{e.tag}" to MOAB structure'
                )
            if (e == self.entities[-1]):
                is_last=True

            moab_core = self.add_stl_to_moab_core(
                moab_core, sid, vid, e.tag, moab_tags, e.stl,
                last=is_last
            )
            vid += 1
            sid += 1
            if self.delete_intermediate:
                p = pl.Path(e.stl)
                p.unlink()

        all_sets = moab_core.get_entities_by_handle(0)

        file_set = moab_core.create_meshset()

        moab_core.add_entities(file_set, all_sets)

        if self.verbose > 0:
            print(f'INFO: writing geometry to h5m: "{h5m_file}".')
        moab_core.write_file(str(h5m_p))

        self.check_h5m_file(h5m_file)

        if vtk:
            moab_core.write_file(str(h5m_p.with_suffix(".vtk")))

        return str(h5m_p)

    def check_h5m_file(self, h5m_file: str = "dagmc.h5m"):
        with open(h5m_file, "rb") as f:
            magic_bytes = f.read(8)
            if magic_bytes != b"\x89HDF\x0d\x0a\x1a\x0a":
                raise RuntimeError(
                    "Generated h5m-file does not appear to be a hdf-file"
                )

    def add_entities_to_moab_core(self, mbcore: core.Core, mbtags: dict):
        vsets = []
        glob_id = 0
        for i in range(len(self.entities)):
            vset = mbcore.create_meshset()
            vsets.append(vset)
            glob_id += 1
            mbcore.tag_set_data(mbtags["global_id"], vset, glob_id)
            mbcore.tag_set_data(mbtags["geom_dimension"], vset, 3)
            mbcore.tag_set_data(mbtags["category"], vset, "Volume")

        faces_added = {}
        sid = 0
        gid = 0
        for i, e in enumerate(self.entities):
            for j, T in enumerate(e.stl):
                f, sense = T
                if f not in faces_added:
                    fset = mbcore.create_meshset()
                    sid += 1
                    glob_id += 1
                    faces_added[f] = fset

                    mbcore.tag_set_data(mbtags["global_id"], fset, glob_id)
                    mbcore.tag_set_data(mbtags["geom_dimension"], fset, 2)
                    mbcore.tag_set_data(mbtags["category"], fset, "Surface")

                    mbcore.add_parent_child(vsets[i], fset)
                    if len(sense) == 2:
                        mbcore.tag_set_data(
                            mbtags["surf_sense"],
                            fset,
                            np.array(
                                [vsets[sense[0]], vsets[sense[1]]], dtype="uint64"
                            ),
                        )
                    else:
                        mbcore.tag_set_data(
                            mbtags["surf_sense"],
                            fset,
                            np.array([vsets[sense[0]], 0], dtype="uint64"),
                        )
                    mbcore.load_file(f, fset)
                else:
                    # this face has already been added so only add a parent child relation here
                    fset = faces_added[f]
                    mbcore.add_parent_child(vsets[i], fset)
            # make this a group, this could ideally be a set of volumes with the same material
            gset = mbcore.create_meshset()
            gid += 1
            glob_id += 1
            mbcore.tag_set_data(mbtags["category"], gset, "Group")
            # reflective is a special case that should not have mat: in front
            if not e.tag == "reflective":
                dagmc_material_tag = f"mat:{e.tag}"
            else:
                dagmc_material_tag = e.tag

            mbcore.tag_set_data(mbtags["name"], gset, dagmc_material_tag)
            mbcore.tag_set_data(mbtags["geom_dimension"], gset, 4)

            # add the volume to this group set
            mbcore.add_entity(gset, vsets[i])

        # if wanted add an implicit complement material
        if ( self.implicit_complement is not None ):
            gset = mbcore.create_meshset()
            mbcore.tag_set_data(mbtags["category"], gset, "Group")
            mbcore.tag_set_data(mbtags["name"], gset, f"mat:{self.implicit_complement}_comp")
            mbcore.tag_set_data(mbtags["geom_dimension"], gset, 4)
            mbcore.add_entity(gset, vsets[-1])

        # finally set the faceting tolerance tag
        # this needs to be attached to some entity so create a dummy vertex
        vh=mbcore.create_vertices(np.array([1,1,1],dtype='float64'))
        mbcore.tag_set_data(
            mbtags["faceting_tol"], vh, np.array((mesher_config["tolerance"],))
        )

        return mbcore

    def add_stl_to_moab_core(
        self,
        mbcore: core.Core,
        surface_id: int,
        volume_id: int,
        material_name: str,
        mbtags: dict,
        stl_filename: str,
        last: bool = False
    ) -> core.Core:
        """
        Appends a set of surfaces (comprising a volume) from an stl-file to a moab.Core object and returns the updated object

        Args:
            mbcore: A moab core object
            surface_id: the id number to apply to the surface
            volume_id: the id numbers to apply to the volumes
            material_name: the material tag name to add. the value provided will
                be prepended with "mat:" unless it is "reflective" which is
                a special case and therefore will remain as is.
            mbtags: A dictionary of the MOAB tags
            stl_filename: the filename of the stl file to load into the moab core

        Returns:
            An updated pymoab.core.Core() instance
        """

        surface_set = mbcore.create_meshset()
        volume_set = mbcore.create_meshset()

        # recent versions of MOAB handle this automatically
        # but best to go ahead and do it manually
        mbcore.tag_set_data(mbtags["global_id"], volume_set, volume_id)
        mbcore.tag_set_data(mbtags["global_id"], surface_set, surface_id)

        # set geom IDs
        mbcore.tag_set_data(mbtags["geom_dimension"], volume_set, 3)
        mbcore.tag_set_data(mbtags["geom_dimension"], surface_set, 2)
        # set category tag values
        mbcore.tag_set_data(mbtags["category"], volume_set, "Volume")
        mbcore.tag_set_data(mbtags["category"], surface_set, "Surface")

        # establish parent-child relationship
        mbcore.add_parent_child(volume_set, surface_set)

        # Set surface sense
        # This should be fixed - we should know which volume comes next, instead of just setting it to be 0
        sense_data = [volume_set, np.uint64(0)]
        mbcore.tag_set_data(mbtags["surf_sense"], surface_set, sense_data)

        # load the stl triangles/vertices into the surface set
        mbcore.load_file(stl_filename, surface_set)

        group_set = mbcore.create_meshset()
        mbcore.tag_set_data(mbtags["category"], group_set, "Group")

        # reflective is a special case that should not have mat: in front
        if not material_name == "reflective":
            dag_material_tag = "mat:{}".format(material_name)
        else:
            dag_material_tag = material_name

        mbcore.tag_set_data(mbtags["name"], group_set, dag_material_tag)

        mbcore.tag_set_data(mbtags["geom_dimension"], group_set, 4)

        # add the volume to this group set
        mbcore.add_entity(group_set, volume_set)

        # if this is the last element also add some extra information
        if ( last ):
            # if wanted add an implicit complement material
            if ( self.implicit_complement is not None ):
                gset = mbcore.create_meshset()
                mbcore.tag_set_data(mbtags["category"], gset, "Group")
                mbcore.tag_set_data(mbtags["name"], gset, f"mat:{self.implicit_complement}_comp")
                mbcore.tag_set_data(mbtags["geom_dimension"], gset, 4)
                mbcore.add_entity(gset, volume_set)

            # finally set the faceting tolerance tag
            # this needs to be attached to some entity so create a dummy vertex
            vh=mbcore.create_vertices(np.array([1,1,1],dtype='float64'))
            mbcore.tag_set_data(
                mbtags["faceting_tol"], vh, np.array((mesher_config["tolerance"],))
            )

        return mbcore

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
        tags["surf_sense"] = moab_core.tag_get_handle(
            sense_tag_name,
            sense_tag_size,
            types.MB_TYPE_HANDLE,
            types.MB_TAG_SPARSE,
            create_if_missing=True,
        )
        tags["category"] = moab_core.tag_get_handle(
            types.CATEGORY_TAG_NAME,
            types.CATEGORY_TAG_SIZE,
            types.MB_TYPE_OPAQUE,
            types.MB_TAG_SPARSE,
            create_if_missing=True,
        )
        tags["name"] = moab_core.tag_get_handle(
            types.NAME_TAG_NAME,
            types.NAME_TAG_SIZE,
            types.MB_TYPE_OPAQUE,
            types.MB_TAG_SPARSE,
            create_if_missing=True,
        )
        tags["geom_dimension"] = moab_core.tag_get_handle(
            types.GEOM_DIMENSION_TAG_NAME,
            1,
            types.MB_TYPE_INTEGER,
            types.MB_TAG_DENSE,
            create_if_missing=True,
        )
        tags["faceting_tol"] = moab_core.tag_get_handle(
            "FACETING_TOL",
            1,
            types.MB_TYPE_DOUBLE,
            types.MB_TAG_DENSE,
            create_if_missing=True,
        )
        # Global ID is a default tag, just need the name to retrieve
        tags["global_id"] = moab_core.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
        return moab_core, tags

    def merge_all(self):
        # merging a single object does not really make sense
        if len(self.entities) > 0:
            # extract cq solids backend algorithm
            unmerged = [e.solid for e in self.entities]

            # Pre-calculate and cache the volume, bounding box, and center of each
            unmerged_vols = [x.Volume() for x in unmerged]
            unmerged_bbs = [x.BoundingBox() for x in unmerged]
            unmerged_centers = [x.Center() for x in unmerged]

            # do merge
            # merged = self._merge_solids(unmerged, fuzzy_value=1e-6)
            merged = self._merge_solids_full(unmerged, fuzzy_value=1e-6)
            # the merging process may result in extra volumes.
            # We need to make sure these are at the end of the list
            # If not this results in a loss of volumes in the end.
            print("INFO: reordering volumes after merge")
            tmp_ents = []

            # Figure of which of the merged solids best corresponds to
            # each of the unmerged volumes.
            merged_solids = merged.Solids()

            # Pre-calculate and cache the volume, bounding box, and center of each
            merged_vols = [x.Volume() for x in merged_solids]
            merged_bbs = [x.BoundingBox() for x in merged_solids]
            merged_centers = [x.Center() for x in merged_solids]

            for j, orig in enumerate(unmerged):
                d_small = 1e9
                i_small = -1
                if self.verbose > 1:
                    print(
                        f"INFO: {len(merged_solids)} merged solids left in list of originally {len(merged.Solids())}"
                    )
                for i, ms in enumerate(merged_solids):
                    d = similar_solids(
                        unmerged_vols[j],
                        unmerged_bbs[j],
                        unmerged_centers[j],
                        merged_vols[i],
                        merged_bbs[i],
                        merged_centers[i],
                    )
                    if d < d_small:
                        i_small, d_small = i, d
                if i_small == -1:
                    print(
                        f"WARNING: Could not find a matching merged volume for volume {j+1}.",
                        end=" ",
                    )
                    print(
                        "This volume/entity will be skipped. Please examine the output volume carefully."
                    )
                else:
                    # Transfer the picked solid to the list of merged solids, removing (pop) it from the list
                    # This to avoid going through the whole list more than once.
                    ent = self.entities[j]
                    ent.solid = merged_solids.pop(i_small)
                    merged_vols.pop(i_small)
                    merged_bbs.pop(i_small)
                    merged_centers.pop(i_small)
                    tmp_ents.append(ent)
            self.entities = tmp_ents

    def imprint_all(self):
        if len(self.entities) > 0:
            unmerged = [e.solid for e in self.entities]
            merged = self.imprint_solids(unmerged)
            for e, m in zip(self.entities, merged):
                e.solid = m

    def _merge_solids_full(self, solids, fuzzy_value=1e-6):
        """
        adds all solids on the solids-vector to the algorithm and
        merges them as arguments. This is stabler, but more resource
        intesive
        """
        bldr = OCP.BOPAlgo.BOPAlgo_MakeConnected()
        for shape in solids:
            if isinstance(shape, cq.occ_impl.shapes.Compound):
                bldr.AddArgument(shape.wrapped)
            else:
                try:
                    bldr.AddArgument(shape.val().wrapped)
                except Exception as _e:
                    bldr.AddArgument(shape.wrapped)

        bldr.SetRunParallel(False)
        bldr.Perform()

        if self.verbose > 1:
            print("INFO: Generate compound shape")
        merged = cq.Compound(bldr.Shape())
        return merged

    def _merge_solids(self, solids, fuzzy_value):
        """merge a set of cq-solids
        returns as cq-compound object
        """
        bldr = OCP.BOPAlgo.BOPAlgo_Splitter()
        bldr.SetFuzzyValue(fuzzy_value)
        # loop through all objects in geometry and split and merge them accordingly
        # shapes should be a compund cq object or a list thereof
        for i, shape in enumerate(solids):
            if self.verbose:
                print(f"splitting obj {i} of {len(solids)}")
            # checks if solid is a compound as .val() is not needed for compunds
            if isinstance(shape, cq.occ_impl.shapes.Compound):
                bldr.AddArgument(shape.wrapped)
            else:
                try:
                    bldr.AddArgument(shape.val().wrapped)
                except Exception as _e:
                    bldr.AddArgument(shape.wrapped)
            for j, shape2 in enumerate(solids):
                if i == j:
                    # don't split an object with itself
                    continue
                print(f"splitting object {i} with object {j}")
                if isinstance(shape2, cq.occ_impl.shapes.Compound):
                    bldr.AddArgument(shape2.wrapped)
                else:
                    try:
                        bldr.AddArgument(shape2.val().wrapped)
                    except Exception as _e:
                        bldr.AddArgument(shape2.wrapped)
                bldr.Perform()

        if self.verbose > 1:
            print("INFO: Commence perform step of merge")
        bldr.Perform()

        if self.verbose > 1:
            print("INFO: Commence image step of merge")
        bldr.Images()

        if self.verbose > 1:
            print("INFO: Generate compound shape")
        merged = cq.Compound(bldr.Shape())

        return merged

    def imprint_solids(self, solids):
        merged = solids
        for i in range(len(merged)):
            s0 = merged[i]
            if self.verbose > 0:
                print(f"INFO: Self_imprinting {i}")
            s0 = self.imprint_solid_on_self(s0)
            merged[i] = s0

        for i in range(len(merged)):
            if self.verbose > 0:
                print(f"INFO: Imprint on solid {i}")
            s0 = merged[i]
            for j in range(0, len(merged)):
                if j == i:
                    continue
                    # imprinting solid on itself - meaning we imprint its faces on each other.
                    s0 = self.imprint_solid_on_self(s0)
                    s0.mesh(mesher_config["tolerance"])
                    merged[i] = s0
                else:
                    s1 = merged[j]
                if self.verbose > 1:
                    print(f"INFO: Imprinting solid {j} on {i}")
                result = self.imprint_solid_on_solid(s0, s1)
                if result is None:
                    if self.verbose > 1:
                        print(
                            f"INFO: solid {j}'s bounding box does not touch/overlap that of solid {i}. Skipping imprinting"
                        )
                    merged[i] = s0
                    merged[j] = s1
                else:
                    if self.verbose > 1:
                        print(
                            f"INFO: solid {j} is possibly connected to solid {i} - perform imprint."
                        )
                    s0_i = result
                    for k, idx in enumerate([i, j]):
                        try:
                            merged[idx] = s0_i.Solids()[k]
                        except Exception as _e:
                            merged[idx] = merged[idx]
        return merged

    def imprint_solid_on_self(self, s0):
        bldr = OCP.BOPAlgo.BOPAlgo_MakeConnected()
        for fc in s0.Faces():
            bldr.AddArgument(fc.wrapped)
        bldr.SetUseOBB(True)
        bldr.Perform()
        s1 = cq.Shape(bldr.Shape())
        return s1

    def overlap_bounding_boxes(self, solid1, solid2):
        (bb1, bb2) = (solid1.BoundingBox(), solid2.BoundingBox())
        outside = (
            (bb2.xmin > bb1.xmax)
            or (bb2.xmax < bb1.xmin)
            or (bb2.ymin > bb1.ymax)
            or (bb2.ymax < bb1.ymin)
            or (bb2.zmin > bb1.zmax)
            or (bb2.zmax < bb1.zmin)
        )
        return not outside

    def imprint_solid_on_solid(self, solid0, solid1):
        # if the bounding boxes of the solids don't overlap
        # don't do anything
        if not self.overlap_bounding_boxes(solid0, solid1):
            return None

        bldr = OCP.BOPAlgo.BOPAlgo_MakeConnected()
        if isinstance(solid0, cq.occ_impl.shapes.Compound):
            bldr.AddArgument(solid0.wrapped)
        else:
            try:
                bldr.AddArgument(solid0.val().wrapped)
            except Exception as _e:
                bldr.AddArgument(solid0.wrapped)
        if isinstance(solid1, cq.occ_impl.shapes.Compound):
            bldr.AddArgument(solid1.wrapped)
        else:
            try:
                bldr.AddArgument(solid1.val().wrapped)
            except Exception as _e:
                bldr.AddArgument(solid1.wrapped)
        bldr.SetUseOBB(True)
        bldr.Perform()
        return cq.Compound(bldr.Shape())

    def heal_stls(self, stls):
        if self.verbose > 0:
            print("INFO: checking surfaces and reparing normals")

        for e in self.entities:
            stl = e.stl
            mesh = trimesh.load_mesh(stl)
            if self.verbose > 1:
                print("INFO: stl-file", stl, ": mesh is watertight", mesh.is_watertight)
            trimesh.repair.fix_normals(
                mesh
            )  # reqired as gmsh stl export from brep can get the inside outside mixed up
            new_filename = stl[:-4] + "_with_corrected_face_normals.stl"
            mesh.export(new_filename)
            if self.delete_intermediate:
                p = pl.Path(e.stl)
                p.unlink()
            e.stl = new_filename

    def get_all_tags(self):
        # extract_all_tags from the list of entities
        return [e.tag for e in self.entitites]

    def get_unique_tags(self):
        # extract a set of unique tags
        return {self.get_all_tags()}
