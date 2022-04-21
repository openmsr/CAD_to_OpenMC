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
    
    @property
    def largest_dimension(self):
        """Calculates a bounding box for the Reactor and returns the largest
        absolute value of the largest dimension of the bounding box"""
        largest_dimension = 0

        if self.largest_shapes is None:
            shapes_to_bound = self.shapes_and_components
        else:
            shapes_to_bound = self.largest_shapes

        for component in shapes_to_bound:
            largest_dimension = max(largest_dimension, component.largest_dimension)
        # self._largest_dimension = largest_dimension
        return largest_dimension

    @largest_dimension.setter
    def largest_dimension(self, value):
        self._largest_dimension = value

    @property
    def largest_shapes(self):
        return self._largest_shapes

    @largest_shapes.setter
    def largest_shapes(self, value):
        if not isinstance(value, (list, tuple, type(None))):
            raise ValueError("paramak.Reactor.largest_shapes should be a " "list of paramak.Shapes")
        self._largest_shapes = value
    
    @property
    def name(self):
        """Returns a list of names of the individual Shapes that make up the
        reactor"""

        all_names = []
        for shape in self.shapes_and_components:
            all_names.append(shape.name)

        return all_names
    
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
                f"The Reactor contains {len(self.shapes_and_components)} "
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
        from brep_to_h5m import brep_to_h5m
        import brep_part_finder as bpf

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
            import OCP

            bldr = OCP.BOPAlgo.BOPAlgo_Splitter()

            for shape in self.shapes_and_components:
                # checks if solid is a compound as .val() is not needed for compunds
                if isinstance(shape.solid, cq.occ_impl.shapes.Compound):
                    bldr.AddArgument(shape.solid.wrapped)
                else:
                    bldr.AddArgument(shape.solid.val().wrapped)

            bldr.SetNonDestructive(True)

            bldr.Perform()

            bldr.Images()

            merged = cq.Compound(bldr.Shape())

            merged.exportBrep(str(path_filename))

        return str(path_filename)

    def merge_surfaces(self):
        """Run through the assembly and merge concurrent surfaces.
        """

    def export_h5m(self, merge_surfaces=False):
        pass


    def 
