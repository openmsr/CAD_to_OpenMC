#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(
    prog="c2omc_step2h5m.py",
    description="""
    The c2omc_step2h5m.py script is part of the CAD_to_OpenMC distribution. This script
    is intended as quick executable to access most of the features within CAD_to_OpenMC.
    """,
)
parser.add_argument("stepfile", nargs='*')
parser.add_argument(
    "--tol", "-t", default=1e-3, type=float, help="Meshing tolerance. [m]"
)
parser.add_argument(
    "--atol", "-a", default=1e-2, type=float, help="Angular meshing tolerance. [rad]"
)
parser.add_argument(
    "--nomerge",
    action="store_true",
    help="Do not perform a merge operation. speeds up meshing for geometries with no shared surfaces.",
)
parser.add_argument(
    "--imprint",
    action="store_true",
    help="Instead of merging, do a round-robin imprint operation (experimental).",
)
parser.add_argument(
    "--threads",
    "-n",
    type=int,
    default=1,
    help="The mesher will use THREADS simultaneous processes.",
)
parser.add_argument(
    "--verbose",
    "-V",
    type=int,
    default=0,
    choices=[0, 1, 2],
    help="Level of verbosity. 0-2 where 2 is a lot of output.",
)
parser.add_argument(
    "--backend",
    "-b",
    choices=["stl2", "stl", "gmsh", "db"],
    default="stl2",
    help="The meshing backend to use. Generally stl2 or db (experimental) performs the best.",
)
parser.add_argument(
    "--implc", "-i", default=None, help="Material tag for the implicit complement."
)
parser.add_argument(
    "--volskip",
    "-k",
    type=int,
    nargs="*",
    default=[],
    help="Skip meshing volume nbr VOLSKIP. May be repeated to skip multiple volumes.",
)
parser.add_argument(
    "--translate",
    "-r",
    help="Translation vector for step-geometry/ies. this will be applied during import. if a list of translation vectors is supplied, each vector is mapped to the volumes in the series.",
)
parser.add_argument(
    "--rotate",
    "-R",
    help="Rotation of step-geometry. this will be applied during import. if a list of rotations is supplied, each vector is mapped to the volumes in the series.",
)
parser.add_argument(
    "--cleanup",
    action='store_true',
    help="Delete intermediate files instead of leaving them in place."
)

args = parser.parse_args()

import pathlib as pl
import CAD_to_OpenMC.assembly as ab

ab.mesher_config["threads"] = args.threads
ab.mesher_config["tolerance"] = args.tol
ab.mesher_config["angular_tolerance"] = args.atol

for sf in args.stepfile:
    a = ab.Assembly(
        [sf],
    )

    a.verbose = args.verbose
    a.cleanup = args.cleanup
    if args.implc is not None:
        a.implicit_complement = args.implc

    a.run(
        merge=(not args.nomerge),
        imprint=args.imprint,
        h5m_filename=pl.Path(sf).with_suffix(".h5m").name,
        backend=args.backend,
        vol_skip=args.volskip,
    )
