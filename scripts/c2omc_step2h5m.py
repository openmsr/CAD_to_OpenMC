#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(
    prog="step2h5m.py",
    description="""
    The step2h5m.py script is part of the CAD_to_OpenMC distribution. This script
    is intended as quick executable to access most of the features within CAD_to_OpenMC.
    """,
)
parser.add_argument("stepfile")
parser.add_argument(
    "--tol", "-t", default=1e-3, type=float, help="meshing tolerance. [m]"
)
parser.add_argument(
    "--atol", "-a", default=1e-2, type=float, help="angular meshing tolerance. [rad]"
)
parser.add_argument(
    "--nomerge",
    action="store_true",
    help="do not perform a merge operation. speeds up meshing for geometries with no shared surfaces.",
)
parser.add_argument(
    "--imprint",
    action="store_true",
    help="instead of merging, do a round-robin imprint operation (experimental).",
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
    help="level of verbosity. 0-2 where 2 is a lot of output.",
)
parser.add_argument(
    "--backend",
    "-b",
    choices=["stl2", "stl", "gmsh"],
    default="stl2",
    help="the meshing backend to use. generally stl2 performs the best.",
)
parser.add_argument(
    "--implc", "-i", default=None, help="material tag for the implicit complement."
)
parser.add_argument(
    "--volskip",
    "-k",
    type=int,
    nargs="*",
    default=[],
    help="skip meshing volume nbr VOLSKIP. may be repeated to skip multiple volumes.",
)
parser.add_argument(
    "--translate",
    "-r",
    help="translation vector for step-geometry/ies. this will be applied during import. if a list of translation vectors is supplied, each vector is mapped to the volumes in the series.",
)
parser.add_argument(
    "--rotate",
    "-R",
    help="rotation of step-geometry. this will be applied during import. if a list of rotations is supplied, each vector is mapped to the volumes in the series.",
)
parser.add_argument(
    "--cleanup",
    action='store_true',
    help="leave intermediate files in place instead of deleting them."
)

args = parser.parse_args()

import pathlib as pl
import CAD_to_OpenMC.assembly as ab

a = ab.Assembly(
    [args.stepfile],
)

ab.mesher_config["threads"] = args.threads
ab.mesher_config["tolerance"] = args.tol
ab.mesher_config["angular_tolerance"] = args.atol

a.verbose = args.verbose
a.cleanup = args.cleanup
if args.implc is not None:
    a.implicit_complement = args.implc
a.run(
    merge=(not args.nomerge),
    imprint=args.imprint,
    h5m_filename=pl.Path(args.stepfile).with_suffix(".h5m").name,
    backend=args.backend,
    vol_skip=args.volskip,
)
