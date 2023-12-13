#!/usr/bin/env python
import pathlib as pl
import argparse

import CAD_to_OpenMC.assembly as ab

parser = argparse.ArgumentParser(prog='step2h5m.py', description="""
    The step2h5m.py script is part of the CAD_to_OpenMC distribution. This script
    is intended as quick executable to access most of the features within CAD_to_OpenMC.
    """)
parser.add_argument('stepfile')
parser.add_argument('--tol',default=1e-3,type=float, help='meshing tolerance. [m]')
parser.add_argument('--atol',default=1e-2,type=float, help='angular meshing tolerance. [rad]')
parser.add_argument('--nomerge',action='store_true', help='do not perform a merge operation. speeds up meshing for geometries with no shared surfaces.')
parser.add_argument('--imprint', action='store_true', help='instead of merging, do a round-robin imprint operation (experimental).')
parser.add_argument('--threads', type=int, default=1, help='The mesher will use THREADS simultaneous processes.')
parser.add_argument('--verbose','-V',type=int, default=0, choices=[0,1,2], help='level of verbosity. 0-2 where 2 is a lot of output.')
parser.add_argument('--backend',choices=['stl2','stl','gmsh'],default='stl2', help='the meshing backend to use. generally stl2 performs the best.')
parser.add_argument('--implc',default=None,help='material tag for the implicit complement.')
args=parser.parse_args()

a=ab.Assembly([args.stepfile])

ab.mesher_config['threads']=args.threads
ab.mesher_config['tolerance']=args.tol
ab.mesher_config['angular_tolerance']=args.atol
a.verbose=args.verbose
if args.implc is not None:
    a.implicit_complement=args.implc
a.run(merge=(not args.nomerge), imprint=args.imprint, h5m_filename=pl.Path(args.stepfile).with_suffix('.h5m'), backend=args.backend)
