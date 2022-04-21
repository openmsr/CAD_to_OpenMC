import gmsh
import cadquery as cq
import OCP

from paramak import Shapes


class Assembly:
    """This class encapsulates a set of geometries defined by step-files
    addtionally it provides access to meshing-utilities, and export to a DAGMC-enabled
    h5m scene, which may be used for neutronics.
    This class is based on (and borrows heavily from) the paramak package. 
    """
    def __init__():
        pass

    def load_stp_file(filename: str, scale_factor: float = 1.0):
        """Loads a stp file and makes the 3D solid and wires available for use.

        Args:
            filename: the filename used to save the html graph.
            scale_factor: a scaling factor to apply to the geometry that can be
                used to increase the size or decrease the size of the geometry.
                Useful when converting the geometry to cm for use in neutronics
                simulations.

        Returns:
            CadQuery.solid, CadQuery.Wires: solid and wires belonging to the object
        """

        part = cq.importers.importStep(str(filename)).val()

        scaled_part = part.scale(scale_factor)
        solid = scaled_part
        wire = scaled_part.Wires()
        return solid, wire

    #export reactor to stp
    def export_stp(
        self,
        filename: Union[List[str], str] = None,
        mode: Optional[str] = "solid",
        units: Optional[str] = "mm",
    ) -> Union[List[str], str]:
        """Exports the 3D reactor model as a stp file or files.

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
            for entry in self.shapes_and_components:
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

        # exports the reactor solid as a separate stp files
        if len(filename) != len(self.shapes_and_components):
            msg = (
                f"The Assembly contains {len(self.shapes_and_components)} "
                f"Shapes and {len(filename)} filenames have be provided. "
                f"The names of the shapes are {self.name}"
            )
            raise ValueError(msg)

        for stp_filename, entry in zip(filename, self.shapes_and_components):

            entry.export_stp(
                filename=stp_filename,
                mode=mode,
                units=units,
                verbose=False,
            )

            if units == "cm":
                _replace(stp_filename, "SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)")

        return filename

    def export_stl(
        self,
        filename: Union[List[str], str] = None,
        tolerance: float = 0.001,
        angular_tolerance: float = 0.1,
    ) -> Union[str, List[str]]:
        """Writes stl files (CAD geometry) for each Shape object in the reactor

        Args:
            filename: Accepts a single filename as a string which exports the
                full reactor model to a single file. Alternativley filename can
                also accept a list of strings where each string is the filename
                of the the individual shapes that make it up. This will result
                in separate files for each shape in the reactor. Defaults to
                None which uses the Reactor.name with '.stl' appended to the end
                of each entry.
            tolerance (float):  the precision of the faceting
            include_graveyard: specify if the graveyard will be included or
                not. If True the the Reactor.make_graveyard will be called
                using Reactor.graveyard_size and Reactor.graveyard_offset
                attribute values.

        Returns:
            list: a list of stl filenames created
        """

        if isinstance(filename, str):

            path_filename = Path(filename)

            if path_filename.suffix != ".stl":
                path_filename = path_filename.with_suffix(".stl")

            path_filename.parents[0].mkdir(parents=True, exist_ok=True)

            # add an include_graveyard that add graveyard if requested
            exporters.export(
                self.solid,
                str(path_filename),
                exportType="STL",
                tolerance=tolerance,
                angularTolerance=angular_tolerance,
            )
            return str(path_filename)

        if filename is None:
            if None in self.name:
                msg = (
                    "Shape.name is None and therefore it can't be used "
                    "to name a stl file. Try setting Shape.name for all "
                    "shapes in the reactor"
                )
                raise ValueError()
            filename = [f"{name}.stl" for name in self.name]

        # exports the reactor solid as a separate stl files
        if len(filename) != len(self.shapes_and_components):
            msg = (
                f"The Reactor contains {len(self.shapes_and_components)} "
                f"Shapes and {len(filename)} filenames have be provided. "
                f"The names of the shapes are {self.name}"
            )
            raise ValueError(msg)

        for stl_filename, entry in zip(filename, self.shapes_and_components):

            entry.export_stl(
                filename=stl_filename,
                tolerance=tolerance,
                verbose=False,
            )

        return filename

    def merge(self):
        """calls the merge capabilities from OCP/OCCT to merge and imprint surfaces in a geometry
            
        """
        
    def export_brep(self,filename = None, merge=True) -> str:
        """Export the gcurrent geometry as a brep-file.
        
        """



#This should be separated into a few things
    def export_dagmc_h5m(
        self,
        filename: str = "dagmc.h5m",
        min_mesh_size: float = 5,
        max_mesh_size: float = 20,
        exclude: List[str] = None,
        verbose=False,
        volume_atol=0.000001,
        center_atol=0.000001,
        bounding_box_atol=0.000001,
    ) -> str:
        """Export a DAGMC compatible h5m file for use in neutronics simulations.
        This method makes use of Gmsh to create a surface mesh of the geometry.
        MOAB is used to convert the meshed geometry into a h5m with parts tagged by
        using the reactor.shape_and_components.name properties. You will need
        Gmsh installed and MOAB installed to use this function. Acceptable
        tolerances may need increasing to match reactor parts with the parts
        in the intermediate Brep file produced during the process

        Args:
            filename: the filename of the DAGMC h5m file to write
            min_mesh_size: the minimum mesh element size to use in Gmsh. Passed
                into gmsh.option.setNumber("Mesh.MeshSizeMin", min_mesh_size)
            max_mesh_size: the maximum mesh element size to use in Gmsh. Passed
                into gmsh.option.setNumber("Mesh.MeshSizeMax", max_mesh_size)
            exclude: A list of shape names to not include in the exported
                geometry. 'plasma' is often excluded as not many neutron
                interactions occur within a low density plasma.
            volume_atol: the absolute volume tolerance to allow when matching
                parts in the intermediate brep file with the cadquery parts
            center_atol: the absolute center coordinates tolerance to allow
                when matching parts in the intermediate brep file with the
                cadquery parts
            bounding_box_atol: the absolute volume tolerance to allow when
                matching parts in the intermediate brep file with the cadquery
                parts
        """

        # a local import is used here as these packages need CQ master to work
        #pick relevant functions from these packages
        
        #get volume tag
        #a material_tag could simply be appended to the cq-objects (no?)


        tmp_brep_filename = tempfile.mkstemp(suffix=".brep", prefix="paramak_")[1]

        # saves the reactor as a Brep file with merged surfaces
        self.export_brep(tmp_brep_filename)

        # brep file is imported
        brep_file_part_properties = bpf.get_brep_part_properties(tmp_brep_filename)

        if verbose:
            print("brep_file_part_properties", brep_file_part_properties)

        shape_properties = {}
        for shape_or_compound in self.shapes_and_components:
            sub_solid_descriptions = []

            # checks if the solid is a cq.Compound or not
            if isinstance(shape_or_compound.solid, cq.occ_impl.shapes.Compound):
                iterable_solids = shape_or_compound.solid.Solids()
            else:
                iterable_solids = shape_or_compound.solid.val().Solids()

            for sub_solid in iterable_solids:
                part_bb = sub_solid.BoundingBox()
                part_center = sub_solid.Center()
                sub_solid_description = {
                    "volume": sub_solid.Volume(),
                    "center": (part_center.x, part_center.y, part_center.z),
                    "bounding_box": (
                        (part_bb.xmin, part_bb.ymin, part_bb.zmin),
                        (part_bb.xmax, part_bb.ymax, part_bb.zmax),
                    ),
                }
                sub_solid_descriptions.append(sub_solid_description)
            shape_properties[shape_or_compound.name] = sub_solid_descriptions

        if verbose:
            print("shape_properties", shape_properties)

        # request to find part ids that are mixed up in the Brep file
        # using the volume, center, bounding box that we know about when creating the
        # CAD geometry in the first place
        key_and_part_id = bpf.get_dict_of_part_ids(
            brep_part_properties=brep_file_part_properties,
            shape_properties=shape_properties,
            volume_atol=volume_atol,
            center_atol=center_atol,
            bounding_box_atol=bounding_box_atol,
        )

        if verbose:
            print(f"key_and_part_id={key_and_part_id}")

        # allows components like the plasma to be removed
        if isinstance(exclude, Iterable):
            for name_to_remove in exclude:
                key_and_part_id = {key: val for key, val in key_and_part_id.items() if val != name_to_remove}

        brep_to_h5m(
            brep_filename=tmp_brep_filename,
            volumes_with_tags=key_and_part_id,
            h5m_filename=filename,
            min_mesh_size=min_mesh_size,
            max_mesh_size=max_mesh_size,
            delete_intermediate_stl_files=True,
        )

        # temporary brep is deleted
        os.remove(tmp_brep_filename)

        return filename


    def tag_geometry_with_mats(volumes,implicit_complement_material_tag,graveyard):
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
                g=re.match("^([^\s_@]+)",part)
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

    #gmsh interface function
    def gmsh_gy_meshsize(graveyard,divisions=5, field=1):
        """set the mesh size field close (and outside) to the graveyard to something coarse(r)"""
        gmsh.model.mesh.field.add("Box",field)
        gmsh.model.mesh.field.setNumber(field,"VIn",10)
        gmsh.model.mesh.field.setNumber(field,"VOut",graveyard/divisions)
        
        gmsh.model.mesh.field.setNumber(field,"XMin",-(graveyard*0.99)/2.0)
        gmsh.model.mesh.field.setNumber(field,"XMax", (graveyard*0.99)/2.0)
        gmsh.model.mesh.field.setNumber(field,"YMin",-(graveyard*0.99)/2.0)
        gmsh.model.mesh.field.setNumber(field,"YMax", (graveyard*0.99)/2.0)
        gmsh.model.mesh.field.setNumber(field,"ZMin",-(graveyard*0.99)/2.0)
        gmsh.model.mesh.field.setNumber(field,"ZMax", (graveyard*0.99)/2.0)
        return field

    def export_stls(self):
        """export all the optionally merged volumes as stl-files
        and returns the list of files. This list may subsequently be iterated upon
        to merge into a h5m-file. Expects that the geometry has been surface-mesh by gmsh
        so we have a list of volumes to operate on.
        We do this be greating gmsh physical groups and we may export only 1 group."""
        stls=[]
        for dim,vid in self.meshed_volumes:
           if (dim!=3):
               #appears not to be a volume - skip
               continue
           ents = gmsh.model.getAdjancencies(dim,vid)
           ps = gmsh.model.setPhysicalName(2,ents[1],f'surfaces_on_volume_{vid}')
           filename=f'volume_{vid}.stl'
           gmsh.write(filename)
           stls.append((vid,filename))
           gmsh.model.removePhysicalGroups([]) # remove group again
        return stls
    
    def heal_stls(self,stls):
        healed=[]
        for stl in stls:
            vid,fn=stls
            mesh = trimesh.load_mesh(filename)
            if (self.verbose):
                print("file", fn, ": mesh is watertight", mesh.is_watertight)
            trimesh.repair.fix_normals(
                mesh
            )  # reqired as gmsh stl export from brep can get the inside outside mixed up
            new_filename = fn[:-4] + "_with_corrected_face_normals.stl"
            mesh.export(new_filename)
            healed.append((vid,new_filename))
        return healed

    #this bit is picked from stl_to_h5m

    def stl2h5m(self,stls,h5m_file='dagmc.h5m') -> str:
        """function that export the lits of stls that we have presumably generated somehow
        and merges them into a DAGMC h5m-file by means of the MOAB-framework.
        """
        h5m_p=pl.Path(h5m_file)
        moab_core,moab_tags = init_moab()

        sid,vid = (1,1)
        for sfn,mtag in stls:
            moab_core = add_stl_to_moab_core(moab_core,sid,vid,mtag, moab_tags, sfn)
            vid += 1
            sid += 1
        
        all_sets = moab_core.get_entities_by_handle(0)

        file_set = moab_core.create_meshset()

        moab_core.add_entities(file_set, all_sets)

        moab_core.write_file(str(path_filename))

        return str(path_filename)

    def add_stl_to_moab_core( moab_core: core.Core, surface_id: int, volume_id: int, material_name: str, tags: dict, stl_filename: str,
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

        # set surface sense
        #This should be fixed I think - we should know which volume comes next
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

    def init_moab():
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

##########from paramak

    def export_brep(self, filename: str, merge: bool = True):
        """Exports a brep file for the Reactor.solid.

        Args:
            filename: the filename of exported the brep file.
            merged: if the surfaces should be merged (True) or not (False).

        Returns:
            filename of the brep created
        """

        path_filename = Path(filename)

        if path_filename.suffix != ".brep":
            msg = "When exporting a brep file the filename must end with .brep"
            raise ValueError(msg)

        path_filename.parents[0].mkdir(parents=True, exist_ok=True)

        if not merge:
            self.solid.exportBrep(str(path_filename))
        else:
            #the merge surface returns a cq-compound object.
            merged = self.merge_surfaces()
            merged.exportBrep(str(path_filename))

        return str(path_filename)

    def merge_surfaces(self):
        """Run through the assembly and merge concurrent surfaces.
        """
        bldr = OCP.BOPAlgo.BOPAlgo_Splitter()
        #loop trough all objects in geometry and split and merge them accordingly
        #shapes should be a compund cq object or a list thereof
        for shape in self.shapes:
          # checks if solid is a compound as .val() is not needed for compunds
          if isinstance(shape.solid, cq.occ_impl.shapes.Compound):
            bldr.AddArgument(shape.solid.wrapped)
          else:
            bldr.AddArgument(shape.solid.val().wrapped)

        bldr.SetNonDestructive(True)

        bldr.Perform()

        bldr.Images()

        merged = cq.Compound(bldr.Shape())

        return merged
    
    def export_h5m(self, merge_surfaces=False):
        pass

    def import_from_step(self,filename="in.step"):
        """Import geometry to the shape list through ocp/occt from the
           given filename"""
         
